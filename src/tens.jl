abstract type AbstractTens{order, dim, T <: Number} <: AbstractArray{T, order} end

@pure get_order(::AbstractTens{order, dim, T}) where {order, dim, T} = order
@pure get_dim(::AbstractTens{order, dim, T}) where {order, dim, T} = dim
@pure Base.eltype(::Type{AbstractTens{order, dim, T}}) where {order, dim, T} = T


"""
    Tens{order,dim,T,A<:AbstractArray}

Tensor type of any order defined by
- a multidata of components (of any type heriting from `AbstractArray`, e.g. `Tensor` or `SymmetricTensor`)
- a basis of `AbstractBasis` type
- a tuple of variances (covariant `:cov` or contravariant `:cont`) of length equal to the `order` of the tensor

# Examples
```julia
julia> ℬ = Basis(Sym[1 0 0; 0 1 0; 0 1 1]) ;

julia> T = Tens(metric(ℬ,:cov),ℬ,(:cov,:cov))
Tens{2, 3, Sym, SymmetricTensor{2, 3, Sym, 6}}
# data: 3×3 SymmetricTensor{2, 3, Sym, 6}:
 1  0  0
 0  2  1
 0  1  1
# basis: 3×3 Tensor{2, 3, Sym, 9}:
 1  0  0
 0  1  0
 0  1  1
# var: (:cov, :cov)

julia> components(T,(:cont,:cov),b)
3×3 Matrix{Sym}:
 1  0  0
 0  1  0
 0  0  1
```
"""
struct Tens{order, dim, T, A <: AbstractArray} <: AbstractTens{order, dim, T}
    data::A
    basis::Basis
    var::NTuple{order, Symbol}
    function Tens(
            data::AbstractArray{T, order},
            basis::Basis{dim},
            var::NTuple{order, Symbol} = ntuple(_ -> :cont, Val(order)),
        ) where {order, dim, T}
        newdata = tensor_or_array(data)
        return new{order, dim, T, typeof(newdata)}(newdata, basis, var)
    end
    Tens(data::AbstractArray, basis::RotatedBasis, args...) = TensRotated(data, basis)
    Tens(data::AbstractArray, basis::OrthogonalBasis, args...) =
        TensOrthogonal(data, basis, args...)
    Tens(
        data::AbstractArray{T},
        basis::OrthonormalBasis = CanonicalBasis{size(data)[1], T}(),
        args...,
    ) where {T} = TensOrthonormal(data, basis)
    Tens(data::AbstractArray, var::NTuple, basis::AbstractBasis = CanonicalBasis{size(data, 1), eltype(data)}()) = Tens(data, basis, var)
    Tens(data::AbstractArray{T, 0}, args...) where {T} = T(data[1])
    Tens(data::T, args...) where {T} = data
end

proj_tens(sym::Symbol, A::AbstractArray) = proj_tens(Val(sym), A)

# Internal helper: loop over symmetry projections and return the first match.
# `proj_fn(sym)` must return `(projected, d, drel)`.
function _best_sym_loop(newt::AbstractTens{order, dim, T}, proj, ε, proj_fn) where {order, dim, T}
    for sym in proj
        (projt, d, drel) = proj_fn(sym)
        if iszero(d) || drel < ε
            return projt, d, drel, sym
        end
    end
    return newt, zero(T), zero(T), :ANISO
end

# ── Structured-tensor reference extraction for the cheap path ────────────────
# If `t` is a structured container that already knows its axis/frame, reuse it
# directly; otherwise derive a candidate from the Kelvin-Mandel eigenstructure
# (or fall back to `e₃` / the canonical frame for non-numeric element types).

_default_TI_axis(t) = _candidate_TI_axis(Array(get_array(t)))
_default_ORTHO_frame(t) = _candidate_ORTHO_frame(Array(get_array(t)))

for order in (2, 4)
    @eval begin
        """
            best_sym_tens(t; proj=(:ISO, :TI, :ORTHO), ε=1e-6, optimize_angles=false)
        
        Find the best (most restrictive) symmetry of tensor `t` by trying each
        symmetry class in `proj` (from most to least symmetric) and accepting the
        first whose relative projection error is below `ε`.
        
        - `optimize_angles=false` (**default**, cheap path, no optimisation): the
          `:ISO` projection is closed-form; for `:TI` the symmetry axis is taken
          from `t` itself (if it is a structured TI container) or derived from the
          Kelvin-Mandel eigenstructure (otherwise); for `:ORTHO` the material frame
          is taken from `t` or derived likewise.  No external optimiser needed.
        - `optimize_angles=true`: the `:TI` axis and `:ORTHO` frame are found by
          nonlinear optimisation (multistart L-BFGS) — requires `using NLopt`.
        
        Returns `(projected, d, drel, sym)` where `sym ∈ {:ISO, :TI, :ORTHO, :ANISO}`.
        
        **Behaviour change (vs. pre-2026 versions):** the default no-argument call
        no longer throws when NLopt is absent; set `optimize_angles=true` to restore
        the previous angle-optimised behaviour.
        
        # Examples
        ```julia
        julia> n = [0., 0., 1.];
        
        julia> C = tens_TI(10., 3., 2.5, 12., 2., n);
        
        julia> _, _, _, sym = best_sym_tens(C);
        
        julia> sym === :TI
        true
        ```
        
        See also [`best_sym_tens(t, n_or_frame)`](@ref) for fixed-axis/frame use,
        [`proj_tens`](@ref).
        """
        function best_sym_tens(
                t::AbstractTens{$order, dim, T};
                proj = (:ISO, :TI, :ORTHO),
                ε = 1.0e-6,
                optimize_angles::Bool = false,
            ) where {dim, T}
            basis = relevant_OrthonormalBasis(get_basis(t))
            newt = change_tens(t, basis)
            A = Array(get_array(newt))
            if optimize_angles
                return _best_sym_loop(newt, proj, ε, sym -> proj_tens(sym, A))
            else
                # Cheap path: reuse structured references when available, otherwise
                # derive candidates from the Kelvin-Mandel eigendecomposition.
                n_default = hasmethod(axis, Tuple{typeof(t)}) ? axis(t) : _default_TI_axis(newt)
                frame_default = hasmethod(frame, Tuple{typeof(t)}) ? frame(t) : _default_ORTHO_frame(newt)
                return _best_sym_loop(
                    newt, proj, ε,
                    function (sym)
                        if sym === :ISO
                            return proj_tens(sym, A)
                        elseif sym === :TI
                            return proj_tens(sym, A, n_default)
                        else    # :ORTHO
                            return proj_tens(sym, A, frame_default)
                        end
                    end,
                )
            end
        end

        """
            best_sym_tens(t, n_or_frame; proj=(:ISO, :TI, :ORTHO), ε=1e-6)
        
        Find the best symmetry of tensor `t` with a **fixed** symmetry axis `n`
        (for TI) or material frame `frame` (for ORTHO).  No rotation optimisation
        is performed.
        
        - `n_or_frame`: a vector (axis for TI) or `OrthonormalBasis{3}` (frame for
          ORTHO). For ISO projection the extra argument is ignored.
        
        Returns `(projected, d, drel, sym)`.
        
        # Examples
        ```julia
        julia> n = [0., 0., 1.];
        
        julia> C = tens_TI(10., 3., 2.5, 12., 2., n);
        
        julia> _, _, drel, sym = best_sym_tens(C, n);
        
        julia> sym == :TI && drel < 1e-12
        true
        ```
        """
        function best_sym_tens(
                t::AbstractTens{$order, dim, T},
                n_or_frame;
                proj = (:ISO, :TI, :ORTHO),
                ε = 1.0e-6,
            ) where {dim, T}
            basis = relevant_OrthonormalBasis(get_basis(t))
            newt = change_tens(t, basis)
            A = Array(get_array(newt))
            return _best_sym_loop(
                newt, proj, ε,
                sym -> sym == :ISO ? proj_tens(sym, A) : proj_tens(sym, A, n_or_frame)
            )
        end
    end
end

struct TensRotated{order, dim, T, A <: AbstractArray} <: AbstractTens{order, dim, T}
    data::A
    basis::RotatedBasis
    function TensRotated(
            data::AbstractArray{T, order},
            basis::RotatedBasis{dim},
        ) where {order, dim, T}
        newdata = tensor_or_array(data)
        return new{order, dim, T, typeof(newdata)}(newdata, basis)
    end
end

struct TensCanonical{order, dim, T, A <: AbstractArray} <: AbstractTens{order, dim, T}
    data::A
    function TensCanonical(data::AbstractArray{T, order}) where {order, T}
        newdata = tensor_or_array(data)
        return new{order, size(data, 1), T, typeof(newdata)}(newdata)
    end
end

struct TensOrthogonal{order, dim, T, A <: AbstractArray} <: AbstractTens{order, dim, T}
    data::A
    basis::OrthogonalBasis
    var::NTuple{order, Symbol}
    function TensOrthogonal(
            data::AbstractArray{T, order},
            basis::OrthogonalBasis{dim},
            var::NTuple{order, Symbol} = ntuple(_ -> :cont, Val(order)),
        ) where {order, dim, T}
        newdata = tensor_or_array(data)
        return new{order, dim, T, typeof(newdata)}(newdata, basis, var)
    end
end

const TensOrthonormal{order, dim, T, A} =
    Union{TensRotated{order, dim, T, A}, TensCanonical{order, dim, T, A}}
const AllTensOrthogonal{order, dim, T, A} = Union{TensOrthonormal{order, dim, T, A}, TensOrthogonal{order, dim, T, A}}
const TensVar{order, dim, T, A} = Union{Tens{order, dim, T, A}, TensOrthogonal{order, dim, T, A}}
const TensBasis{order, dim, T, A} =
    Union{Tens{order, dim, T, A}, TensRotated{order, dim, T, A}, TensOrthogonal{order, dim, T, A}}
const TensArray{order, dim, T, A} = Union{
    Tens{order, dim, T, A},
    TensRotated{order, dim, T, A},
    TensCanonical{order, dim, T, A},
    TensOrthogonal{order, dim, T, A},
}

TensOrthonormal(data::AbstractArray, basis::RotatedBasis) = TensRotated(data, basis)
TensOrthonormal(data::AbstractArray, ::CanonicalBasis) = TensCanonical(data)

TensOrthonormal(t::TensOrthonormal) = t
function TensOrthonormal(t::TensOrthogonal{order, dim, T}) where {order, dim, T}
    m = Array(get_array(t))
    ℬ = get_basis(t)
    onℬ = relevant_OrthonormalBasis(ℬ)
    Λ = Dict(:cov => inv.(ℬ.λ), :cont => ℬ.λ)
    for ind in CartesianIndices(m)
        m[ind] *= prod([Λ[get_var(t, i)][ind[i]] for i in 1:order])
    end
    return Tens(m, onℬ)
end

@inline Base.size(t::TensArray) = size(get_array(t))
@inline Base.getindex(t::TensArray, ind...) = getindex(get_array(t), ind...)
@pure datatype(::TensArray{order, dim, T, A}) where {order, dim, T, A} = A

@inline Base.zero(t::AbstractTens) = Tens(zero.(get_array(t)), get_basis(t), get_var(t))

# This function aims at storing the table of components in the `Tensor` type whenever possible
# Convert a raw array into the best matching Tensors.jl type (Vec, Tensor, or
# SymmetricTensor).  Falls through to the input unchanged for types that are
# already Tensors.AllTensors or for orders not in {1,2,4}.
tensor_or_array(tab::AbstractArray{T, 1}) where {T} = Vec{size(tab, 1)}(tab)
for order in (2, 4)
    @eval function tensor_or_array(tab::AbstractArray{T, $order}) where {T}
        dim = size(tab, 1)
        newtab = Tensor{$order, dim}(tab)
        if Tensors.issymmetric(newtab)
            newtab = convert(SymmetricTensor{$order, dim}, newtab)
        end
        return newtab
    end
    @eval tensor_or_array(tab::Tensor{$order, dim}) where {dim} = Tensors.issymmetric(tab) ? convert(SymmetricTensor{$order, dim}, tab) : tab
end
tensor_or_array(tab::Tensors.AllTensors) = tab
tensor_or_array(tab::AbstractArray) = tab

##############################
# Utility/Accessor Functions #
##############################

get_array(t::TensArray) = t.data
get_basis(t::TensBasis) = t.basis
get_basis(::TensCanonical{order, dim, T}) where {order, dim, T} = CanonicalBasis{dim, T}()
get_var(::TensOrthonormal{order}) where {order} = ntuple(_ -> :cov, Val(order))
get_var(::TensOrthonormal, i::Integer) = :cov
get_var(t::TensVar) = t.var
get_var(t::TensVar, i::Integer) = t.var[i]


#####################
# Display Functions #
#####################
for OP in (:show, :print, :display)
    @eval begin
        Base.$OP(U::FourthOrderTensor) = $OP(tomandel(U))

        function Base.$OP(t::AbstractTens)
            $OP(typeof(t))
            if ndims(t) == 4
                print("→ KM: ")
                $OP(KM(t))
            else
                print("→ array: ")
                $OP(get_array(t))
            end
            print("→ basis: ")
            $OP(vecbasis(get_basis(t)))
            print("→ var: ")
            return $OP(get_var(t))
        end
        function Base.$OP(t::TensOrthonormal)
            $OP(typeof(t))
            if ndims(t) == 4
                print("→ KM: ")
                $OP(KM(t))
            else
                print("→ array: ")
                $OP(get_array(t))
            end
            print("→ basis: ")
            return $OP(vecbasis(get_basis(t)))
        end

        # Base.$OP(t::AbstractTens{order,dim,T}; vec = '𝐞', coords = ntuple(i -> i, dim)) where {order,dim,T} = intrinsic(t; vec= vec, coords = coords)
    end
end

intrinsic(t::T) where {T} = println(t)

function intrinsic(t::AbstractTens{order, dim, T}; vec = '𝐞', coords = ntuple(i -> i, dim)) where {order, dim, T}
    ind = CartesianIndices(t)
    ℬ = get_basis(t)
    firstprint = true
    s = ""
    for i in ind
        x = t[i]
        if !iszero(x)
            if !firstprint
                s *= " + "
            end
            if !isone(x)
                s *= "(" * string(x) * ")"
            end
            j = Tuple(i)
            for k in 1:order
                s *= strvecbasis(ℬ, coords[j[k]], invvar(get_var(t, k)); vec = vec)
                if k < order
                    s *= "⊗"
                end
            end
            firstprint = false
        end
    end
    return if length(s) > 0
        println(s)
    else
        println(0)
    end
end


########################
# Component extraction #
########################

"""
    components(t::AbstractTens{order,dim,T},ℬ::AbstractBasis{dim},var::NTuple{order,Symbol})
    components(t::AbstractTens{order,dim,T},ℬ::AbstractBasis{dim})
    components(t::AbstractTens{order,dim,T},var::NTuple{order,Symbol})

Extract the components of a tensor for new variances and/or in a new basis

# Examples
```julia
julia> ℬ = Basis(Sym[0 1 1; 1 0 1; 1 1 0]) ;

julia> TV = Tens(Tensor{1,3}(i->symbols("v\$i",real=true)))
TensND.TensCanonical{1, 3, Sym, Vec{3, Sym}}
# data: 3-element Vec{3, Sym}:
 v₁
 v₂
 v₃
# basis: 3×3 TensND.LazyIdentity{3, Sym}:
 1  0  0
 0  1  0
 0  0  1
# var: (:cont,)

julia> factor.(components(TV, ℬ, (:cont,)))
3-element Vector{Sym}:
 -(v1 - v2 - v3)/2
  (v1 - v2 + v3)/2
  (v1 + v2 - v3)/2

julia> components(TV, ℬ, (:cov,))
3-element Vector{Sym}:
 v₂ + v₃
 v₁ + v₃
 v₁ + v₂

julia> simplify.(components(TV, normalize(ℬ), (:cov,)))
3-element Vector{Sym}:
 sqrt(2)*(v2 + v3)/2
 sqrt(2)*(v1 + v3)/2
 sqrt(2)*(v1 + v2)/2

julia> TT = Tens(Tensor{2,3}((i,j)->symbols("t\$i\$j",real=true)))
TensND.TensCanonical{2, 3, Sym, Tensor{2, 3, Sym, 9}}
# data: 3×3 Tensor{2, 3, Sym, 9}:
 t₁₁  t₁₂  t₁₃
 t₂₁  t₂₂  t₂₃
 t₃₁  t₃₂  t₃₃
# basis: 3×3 TensND.LazyIdentity{3, Sym}:
 1  0  0
 0  1  0
 0  0  1
# var: (:cont, :cont)

julia> components(TT, ℬ, (:cov,:cov))
3×3 Matrix{Sym}:
 t₂₂ + t₂₃ + t₃₂ + t₃₃  t₂₁ + t₂₃ + t₃₁ + t₃₃  t₂₁ + t₂₂ + t₃₁ + t₃₂
 t₁₂ + t₁₃ + t₃₂ + t₃₃  t₁₁ + t₁₃ + t₃₁ + t₃₃  t₁₁ + t₁₂ + t₃₁ + t₃₂
 t₁₂ + t₁₃ + t₂₂ + t₂₃  t₁₁ + t₁₃ + t₂₁ + t₂₃  t₁₁ + t₁₂ + t₂₁ + t₂₂

julia> factor.(components(TT, ℬ, (:cont,:cov)))
3×3 Matrix{Sym}:
 -(t12 + t13 - t22 - t23 - t32 - t33)/2  …  -(t11 + t12 - t21 - t22 - t31 - t32)/2
  (t12 + t13 - t22 - t23 + t32 + t33)/2      (t11 + t12 - t21 - t22 + t31 + t32)/2
  (t12 + t13 + t22 + t23 - t32 - t33)/2      (t11 + t12 + t21 + t22 - t31 - t32)/2
```
"""
components(t::AbstractTens) = get_array(t)

components(t::TensOrthonormal, ::NTuple) = get_array(t)

function components(
        t::TensOrthogonal{order, dim, T},
        var::NTuple{order, Symbol},
    ) where {order, dim, T}
    if isequal(var, get_var(t))
        return get_array(t)
    else
        m = Array(get_array(t))
        ℬ = get_basis(t)
        g_or_G = ntuple(i -> isequal(get_var(t, i), var[i]) ? I : metric(ℬ, var[i]), order)
        for ind in CartesianIndices(m)
            m[ind] *= prod([g_or_G[i][ind[i], ind[i]] for i in 1:order])
        end
        return m
    end
end

for B in (AbstractBasis, OrthogonalBasis, OrthonormalBasis)
    @eval function components(
            t::TensOrthogonal{order, dim, T},
            ℬ::$B{dim},
            var::NTuple{order, Symbol},
        ) where {order, dim, T}
        if isequal(ℬ, get_basis(t))
            return components(t, var)
        else
            return components(TensOrthonormal(t), ℬ, var)
        end
    end
end

function components(t::Tens{order, dim, T}, var::NTuple{order, Symbol}) where {order, dim, T}
    if isequal(var, get_var(t))
        return get_array(t)
    else
        m = Array(get_array(t))
        ec1 = ntuple(i -> i, Val(order))
        newcp = order + 1
        for i in 1:order
            if !isequal(get_var(t, i), var[i])
                g_or_G = metric(get_basis(t), var[i])
                ec2 = (i, newcp)
                ec3 = ntuple(j -> j ≠ i ? j : newcp, Val(order))
                m = T.(einsum(EinCode((ec1, ec2), ec3), (m, g_or_G)))
            end
        end
        return m
    end
end

function components(
        t::AbstractTens{order, dim, T},
        ℬ::AbstractBasis{dim},
        var::NTuple{order, Symbol},
    ) where {order, dim, T}
    if isequal(ℬ, get_basis(t))
        return components(t, var)
    else
        bb = Dict{Tuple{Symbol, Symbol}, AbstractMatrix}()
        for v1 in (:cov, :cont), v2 in (:cov, :cont)
            if v1 ∈ get_var(t) && v2 ∈ var
                bb[v1, v2] = vecbasis(get_basis(t), invvar(v1))' * vecbasis(ℬ, v2)
            end
        end
        m = Array(get_array(t))
        ec1 = ntuple(i -> i, Val(order))
        newcp = order + 1
        for i in 1:order
            c = bb[get_var(t, i), var[i]]
            if c ≠ 1I
                ec2 = (i, newcp)
                ec3 = ntuple(j -> j ≠ i ? j : newcp, Val(order))
                m = T.(einsum(EinCode((ec1, ec2), ec3), (m, c)))
            end
        end
        return m
    end
end

function components(
        t::AbstractTens{order, dim, T},
        ℬ::OrthogonalBasis{dim},
        var::NTuple{order, Symbol},
    ) where {order, dim, T}
    if isequal(ℬ, get_basis(t))
        return components(t, var)
    else
        m = Array(components(t, relevant_OrthonormalBasis(ℬ)))
        Λ = Dict(:cont => inv.(ℬ.λ), :cov => ℬ.λ)
        for ind in CartesianIndices(m)
            m[ind] *= prod([Λ[var[i]][ind[i]] for i in 1:order])
        end
        return m
    end
end

components(t::AbstractTens{order, dim, T}, ℬ::AbstractBasis{dim}) where {order, dim, T} =
    components(t, ℬ, get_var(t))

components(t::AbstractTens{order, dim, T}, ℬ::OrthonormalBasis{dim}) where {order, dim, T} =
    components(t, ℬ, ntuple(_ -> :cont, Val(order)))

function components(
        t::TensOrthonormal{order, dim, T},
        ℬ::OrthonormalBasis{dim},
    ) where {order, dim, T}
    if isequal(ℬ, get_basis(t))
        return get_array(t)
    else
        bb = vecbasis(get_basis(t))' * vecbasis(ℬ)
        m = Array(get_array(t))
        ec1 = ntuple(i -> i, Val(order))
        newcp = order + 1
        for i in 1:order
            if bb ≠ 1I
                ec2 = (i, newcp)
                ec3 = ntuple(j -> j ≠ i ? j : newcp, Val(order))
                m = T.(einsum(EinCode((ec1, ec2), ec3), (m, bb)))
            end
        end
        return m
    end
end

components(
    t::TensOrthonormal{order, dim, T},
    basis::OrthonormalBasis{dim},
    ::NTuple{order, Symbol},
) where {order, dim, T} = components(t, basis)

"""
    components_canon(t::AbstractTens)

Extract the components of a tensor in the canonical basis
"""
components_canon(t::AbstractTens) =
    components(t, CanonicalBasis{get_dim(t), eltype(t)}(), get_var(t))

components_canon(t::TensOrthonormal) = components(t, CanonicalBasis{get_dim(t), eltype(t)}())

"""
    change_tens(t::AbstractTens{order,dim,T},ℬ::AbstractBasis{dim},var::NTuple{order,Symbol})
    change_tens(t::AbstractTens{order,dim,T},ℬ::AbstractBasis{dim})
    change_tens(t::AbstractTens{order,dim,T},var::NTuple{order,Symbol})

Rewrite the same tensor with components corresponding to new variances and/or to a new basis

```julia
julia> ℬ = Basis(Sym[0 1 1; 1 0 1; 1 1 0]) ;

julia> TV = Tens(Tensor{1,3}(i->symbols("v\$i",real=true)))
TensND.TensCanonical{1, 3, Sym, Vec{3, Sym}}
# data: 3-element Vec{3, Sym}:
 v₁
 v₂
 v₃
# basis: 3×3 TensND.LazyIdentity{3, Sym}:
 1  0  0
 0  1  0
 0  0  1
# var: (:cont,)

julia> factor.(components(TV, ℬ, (:cont,)))
3-element Vector{Sym}:
 -(v1 - v2 - v3)/2
  (v1 - v2 + v3)/2
  (v1 + v2 - v3)/2

julia> ℬ₀ = Basis(Sym[0 1 1; 1 0 1; 1 1 1]) ;

julia> TV0 = change_tens(TV, ℬ₀)
Tens{1, 3, Sym, Vec{3, Sym}}
# data: 3-element Vec{3, Sym}:
     -v₁ + v₃
     -v₂ + v₃
 v₁ + v₂ - v₃
# basis: 3×3 Tensor{2, 3, Sym, 9}:
 0  1  1
 1  0  1
 1  1  1
# var: (:cont,)
```
"""
function change_tens(t::AbstractTens, ℬ::AbstractBasis, newvar::NTuple)
    if isequal(ℬ, get_basis(t)) && isequal(newvar, get_var(t))
        return t
    else
        return Tens(components(t, ℬ, newvar), ℬ, newvar)
    end
end

function change_tens(t::AbstractTens, newbasis::AbstractBasis)
    if isequal(newbasis, get_basis(t))
        return t
    else
        return Tens(components(t, newbasis, get_var(t)), newbasis, get_var(t))
    end
end

function change_tens(t::AbstractTens, newvar::NTuple)
    if isequal(newvar, get_var(t))
        return t
    else
        return Tens(components(t, get_basis(t), newvar), get_basis(t), newvar)
    end
end

"""
    change_tens_canon(t::AbstractTens{order,dim,T},var::NTuple{order,Symbol})

Rewrite the same tensor with components corresponding to the canonical basis

```julia
julia> ℬ = Basis(Sym[0 1 1; 1 0 1; 1 1 0]) ;

julia> TV = Tens(Tensor{1,3}(i->symbols("v\$i",real=true)), ℬ)
Tens{1, 3, Sym, Vec{3, Sym}}
# data: 3-element Vec{3, Sym}:
 v₁
 v₂
 v₃
# basis: 3×3 Tensor{2, 3, Sym, 9}:
 0  1  1
 1  0  1
 1  1  1
# var: (:cont,)

julia> TV0 = change_tens_canon(TV)
TensND.TensCanonical{1, 3, Sym, Vec{3, Sym}}
# data: 3-element Vec{3, Sym}:
      v₂ + v₃
      v₁ + v₃
 v₁ + v₂ + v₃
# basis: 3×3 TensND.LazyIdentity{3, Sym}:
 1  0  0
 0  1  0
 0  0  1
# var: (:cont,)
```
"""
change_tens_canon(t::AbstractTens) = change_tens(t, CanonicalBasis{get_dim(t), eltype(t)}())


for OP in (:(tsimplify), :(tfactor), :(tsubs), :(ttrigsimp), :(texpand_trig))
    @eval $OP(t::Tensors.AllTensors{dim, T}, args...; kwargs...) where {dim, T <: Sym} =
        Tensors.get_base(typeof(t))($OP.(Tensors.get_data(t), args...; kwargs...))
    @eval $OP(t::AbstractTens{order, dim, T}, args...; kwargs...) where {order, dim, T <: Sym} =
        Tens($OP(get_array(t), args...; kwargs...), $OP(get_basis(t), args...; kwargs...), get_var(t))
end
for OP in (:(tdiff),)
    @eval $OP(t::Tensors.AllTensors{dim, T}, args...; kwargs...) where {dim, T <: Sym} =
        Tensors.get_base(typeof(t))($OP.(Tensors.get_data(t), args...; kwargs...))
    @eval $OP(t::AbstractTens{order, dim, T}, args...; kwargs...) where {order, dim, T <: Sym} =
        Tens($OP(get_array(t), args...; kwargs...), get_basis(t), get_var(t))
end
diff_with_basis(t::AbstractTens{order, dim, T}, args...; kwargs...) where {order, dim, T <: Sym} =
    change_tens(Tens(diff(components_canon(t), args...; kwargs...)), get_basis(t), get_var(t))


for OP in (:(tsimplify), :(tsubs))
    @eval $OP(t::Tensors.AllTensors{dim, Num}, args...; kwargs...) where {dim} =
        Tensors.get_base(typeof(t))($OP.(Tensors.get_data(t), args...; kwargs...))
    @eval $OP(t::AbstractTens{order, dim, Num}, args...; kwargs...) where {order, dim} =
        Tens($OP(get_array(t), args...; kwargs...), $OP(get_basis(t), args...; kwargs...), get_var(t))
end
for OP in (:(tdiff),)
    @eval $OP(t::Tensors.AllTensors{dim, Num}, args...; kwargs...) where {dim} =
        Tensors.get_base(typeof(t))($OP.(Tensors.get_data(t), args...; kwargs...))
    @eval $OP(t::AbstractTens{order, dim, Num}, args...; kwargs...) where {order, dim} =
        Tens($OP(get_array(t), args...; kwargs...), get_basis(t), get_var(t))
end
diff_with_basis(t::AbstractTens{order, dim, Num}, args...; kwargs...) where {order, dim} =
    change_tens(Tens(diff(components_canon(t), args...; kwargs...)), get_basis(t), get_var(t))


##############
# Operations #
##############

choose_best_basis(ℬ::AbstractBasis, ::AbstractBasis) = ℬ

choose_best_basis(ℬ::OrthonormalBasis, ::AbstractBasis) = ℬ
choose_best_basis(ℬ::OrthogonalBasis, ::AbstractBasis) = ℬ
choose_best_basis(::AbstractBasis, ℬ::OrthonormalBasis) = ℬ
choose_best_basis(::AbstractBasis, ℬ::OrthogonalBasis) = ℬ

choose_best_basis(::CanonicalBasis, ℬ::OrthonormalBasis) = ℬ
choose_best_basis(::CanonicalBasis, ℬ::OrthogonalBasis) = ℬ
choose_best_basis(ℬ::OrthonormalBasis, ::CanonicalBasis) = ℬ
choose_best_basis(ℬ::OrthogonalBasis, ::CanonicalBasis) = ℬ

choose_best_basis(::OrthonormalBasis, ℬ::OrthogonalBasis) = ℬ
choose_best_basis(ℬ::OrthogonalBasis, ::OrthonormalBasis) = ℬ
choose_best_basis(ℬ::CanonicalBasis, ::CanonicalBasis) = ℬ
choose_best_basis(ℬ::OrthonormalBasis, ::OrthonormalBasis) = ℬ
choose_best_basis(ℬ::OrthogonalBasis, ::OrthogonalBasis) = ℬ

function same_basis(
        t1::AbstractTens{order1, dim},
        t2::AbstractTens{order2, dim},
    ) where {order1, order2, dim}
    ℬ = choose_best_basis(get_basis(t1), get_basis(t2))
    return change_tens(t1, ℬ), change_tens(t2, ℬ)
end

function same_basis_same_var(
        t1::AbstractTens{order1, dim},
        t2::AbstractTens{order2, dim},
    ) where {order1, order2, dim}
    ℬ = choose_best_basis(get_basis(t1), get_basis(t2))
    return change_tens(t1, ℬ, get_var(t1)), change_tens(t2, ℬ, get_var(t1))
end

# same_basis(
#     t1::AbstractTens{order1,dim},
#     t2::AbstractTens{order2,dim},
# ) where {order1,order2,dim} = t1, change_tens(t2, basis(t1))

# same_basis_same_var(
#     t1::AbstractTens{order1,dim},
#     t2::AbstractTens{order2,dim},
# ) where {order1,order2,dim} = t1, change_tens(t2, basis(t1), get_var(t1))


for OP in (:(==), :(!=), :(isequal))
    @eval @inline function Base.$OP(
            t1::AbstractTens{order, dim},
            t2::AbstractTens{order, dim},
        ) where {order, dim}
        nt1, nt2 = same_basis_same_var(t1, t2)
        return $OP(get_array(nt1), get_array(nt2))
    end
end

for OP in (:+, :-)
    @eval @inline function Base.$OP(
            t1::AbstractTens{order, dim},
            t2::AbstractTens{order, dim},
        ) where {order, dim}
        nt1, nt2 = same_basis_same_var(t1, t2)
        return Tens($OP(get_array(nt1), get_array(nt2)), get_basis(nt1), get_var(nt1))
    end
    @eval @inline function Base.$OP(
            t1::AllTensOrthogonal{order, dim, T},
            t2::UniformScaling{T},
        ) where {order, dim, T <: Sym}
        nt1 = TensOrthonormal(t1)
        return Tens($OP(get_array(nt1), t2), get_basis(nt1), get_var(nt1))
    end
    @eval @inline function Base.$OP(
            t1::UniformScaling{T},
            t2::AllTensOrthogonal{order, dim, T},
        ) where {order, dim, T <: Sym}
        nt2 = TensOrthonormal(t2)
        return Tens($OP(t1, get_array(nt2)), get_basis(nt2), get_var(nt2))
    end
end

@inline Base.:-(t::AbstractTens) = Tens(.-(get_array(t)), get_basis(t), get_var(t))
@inline Base.:*(α::Number, t::AbstractTens) = Tens(α * get_array(t), get_basis(t), get_var(t))
@inline Base.:*(t::AbstractTens, α::Number) = Tens(α * get_array(t), get_basis(t), get_var(t))
@inline Base.:/(t::AbstractTens, α::Number) = Tens(get_array(t) / α, get_basis(t), get_var(t))

@inline Base.inv(t::AbstractTens{2}) =
    Tens(inv(get_array(t)), get_basis(t), (invvar(get_var(t, 2)), invvar(get_var(t, 1))))
@inline Base.inv(t::AbstractTens{4}) = Tens(
    inv(get_array(t)),
    get_basis(t),
    (
        invvar(get_var(t, 3)),
        invvar(get_var(t, 4)),
        invvar(get_var(t, 1)),
        invvar(get_var(t, 2)),
    ),
)

"""
    KM(t::AbstractTens{order,dim}; kwargs...)
    KM(t::AbstractTens{order,dim}, var::NTuple{order,Symbol}, b::AbstractBasis{dim}; kwargs...)

Write the components of a second or fourth order tensor in Kelvin-Mandel notation

# Examples
```julia
julia> σ = Tens(SymmetricTensor{2,3}((i, j) -> symbols("σ\$i\$j", real = true))) ;

julia> KM(σ)
6-element Vector{Sym}:
         σ11
         σ22
         σ33
      √2⋅σ32
      √2⋅σ31
      √2⋅σ21

julia> C = Tens(SymmetricTensor{4,3}((i, j, k, l) -> symbols("C\$i\$j\$k\$l", real = true))) ;

julia> KM(C)
6×6 Matrix{Sym}:
         C₁₁₁₁     C₁₁₂₂     C₁₁₃₃  √2⋅C₁₁₃₂  √2⋅C₁₁₃₁  √2⋅C₁₁₂₁
         C₂₂₁₁     C₂₂₂₂     C₂₂₃₃  √2⋅C₂₂₃₂  √2⋅C₂₂₃₁  √2⋅C₂₂₂₁
         C₃₃₁₁     C₃₃₂₂     C₃₃₃₃  √2⋅C₃₃₃₂  √2⋅C₃₃₃₁  √2⋅C₃₃₂₁
      √2⋅C₃₂₁₁  √2⋅C₃₂₂₂  √2⋅C₃₂₃₃   2⋅C₃₂₃₂   2⋅C₃₂₃₁   2⋅C₃₂₂₁
      √2⋅C₃₁₁₁  √2⋅C₃₁₂₂  √2⋅C₃₁₃₃   2⋅C₃₁₃₂   2⋅C₃₁₃₁   2⋅C₃₁₂₁
      √2⋅C₂₁₁₁  √2⋅C₂₁₂₂  √2⋅C₂₁₃₃   2⋅C₂₁₃₂   2⋅C₂₁₃₁   2⋅C₂₁₂₁
```
"""
KM(t::Tensors.AllTensors; kwargs...) = tomandel(t; kwargs...)
KM(t::AbstractTens; kwargs...) = tomandel(get_array(t); kwargs...)

KM(
    t::AbstractTens{order, dim},
    b::AbstractBasis{dim},
    var::NTuple{order, Symbol},
    kwargs...,
) where {order, dim} = tomandel(tensor_or_array(components(t, b, var)); kwargs...)

KM(t::AbstractTens{order, dim}, b::AbstractBasis{dim}; kwargs...) where {order, dim} =
    tomandel(tensor_or_array(components(t, b)); kwargs...)


KM(t::AbstractArray; kwargs...) = KM(Tens(t); kwargs...)
KM(t::AbstractArray, b::AbstractBasis; kwargs...) = KM(Tens(t), b; kwargs...)

const select_type_KM = Dict(
    (6, 6) => SymmetricTensor{4, 3},
    (9, 9) => Tensor{4, 3},
    (3, 3) => SymmetricTensor{4, 2},
    (4, 4) => Tensor{4, 2},
    (6,) => SymmetricTensor{2, 3},
    (9,) => Tensor{2, 3},
    (3,) => SymmetricTensor{2, 2},
    (4,) => Tensor{2, 2},
)


"""
    inv_KM(v::AbstractVecOrMat; kwargs...)

Define a tensor from a Kelvin-Mandel vector or matrix representation
"""
inv_KM(TT::Type{<:Tensors.AllTensors}, v::AbstractVecOrMat; kwargs...) =
    Tens(frommandel(TT, v; kwargs...))
inv_KM(TT::Type{<:Tensors.AllTensors}, v::AbstractVecOrMat, b::AbstractBasis; kwargs...) =
    Tens(frommandel(TT, v; kwargs...), b)
inv_KM(v::AbstractVecOrMat; kwargs...) = inv_KM(select_type_KM[size(v)], v; kwargs...)
inv_KM(v::AbstractVecOrMat, b::AbstractBasis; kwargs...) = inv_KM(select_type_KM[size(v)], v, b; kwargs...)


"""
    otimes(t1::AbstractTens{order1,dim}, t2::AbstractTens{order2,dim})

Define a tensor product between two tensors

`(aⁱeᵢ) ⊗ (bʲeⱼ) = aⁱbʲ eᵢ⊗eⱼ`
"""
function Tensors.otimes(
        t1::AbstractTens{order1, dim},
        t2::AbstractTens{order2, dim},
    ) where {order1, order2, dim}
    if isequal(t1, t2)
        return otimes(t1)
    else
        nt1, nt2 = same_basis(t1, t2)
        data = otimes(get_array(nt1), get_array(nt2))
        var = (get_var(nt1)..., get_var(nt2)...)
        return Tens(data, get_basis(nt1), var)
    end
end

function Tensors.otimes(
        t1::TensOrthonormal{order1, dim},
        t2::TensOrthonormal{order2, dim},
    ) where {order1, order2, dim}
    if isequal(t1, t2)
        return otimes(t1)
    else
        nt1, nt2 = same_basis(t1, t2)
        data = otimes(get_array(nt1), get_array(nt2))
        return Tens(data, get_basis(nt1))
    end
end

Tensors.otimes(v::AbstractTens{1, dim}) where {dim} =
    Tens(otimes(get_array(v)), get_basis(v), (get_var(v)..., get_var(v)...))

Tensors.otimes(v::TensOrthonormal{1, dim}) where {dim} =
    Tens(otimes(get_array(v)), get_basis(v))

@inline function Tensors.otimes(S::SymmetricTensor{2, dim}) where {dim}
    return SymmetricTensor{4, dim}(
        @inline function (i, j, k, l)
            return @inbounds S[i, j] * S[k, l]
        end
    )
end

Tensors.otimes(t::AbstractTens{2, dim}) where {dim} =
    Tens(otimes(get_array(t)), get_basis(t), (get_var(t)..., get_var(t)...))

Tensors.otimes(t::TensOrthonormal{2, dim}) where {dim} =
    Tens(otimes(get_array(t)), get_basis(t))


function scontract(
        t1::AbstractArray{T1, order1},
        t2::AbstractArray{T2, order2},
    ) where {T1, T2, order1, order2}
    ec1 = ntuple(i -> i, order1)
    ec2 = ntuple(i -> order1 - 1 + i, order2)
    ec3 = (ec1[begin:(end - 1)]..., ec2[(begin + 1):end]...)
    return T1.(einsum(EinCode((ec1, ec2), ec3), (AbstractArray{T1}(t1), AbstractArray{T2}(t2))))
end

scontract(t1::AbstractArray{T1, 1}, t2::AbstractArray{T2, 1}) where {T1, T2} =
    dot(AbstractArray{T1}(t1), AbstractArray{T2}(t2))

for TT1 in (Vec, SecondOrderTensor), TT2 in (Vec, SecondOrderTensor)
    @eval scontract(S1::$TT1, S2::$TT2) = dot(S1, S2)
end

"""
    dot(t1::AbstractTens{order1,dim}, t2::AbstractTens{order2,dim})

Define a contracted product between two tensors

`a ⋅ b = aⁱbⱼ`
"""
function LinearAlgebra.dot(
        t1::AbstractTens{order1, dim},
        t2::AbstractTens{order2, dim},
    ) where {order1, order2, dim}
    nt1, nt2 = same_basis(t1, t2)
    var = (invvar(get_var(nt1)[end]), get_var(nt2)[(begin + 1):end]...)
    nt2 = change_tens(nt2, get_basis(nt2), var)
    data = scontract(get_array(nt1), get_array(nt2))
    var = (get_var(nt1)[begin:(end - 1)]..., get_var(nt2)[(begin + 1):end]...)
    return Tens(data, get_basis(nt1), var)
end

function LinearAlgebra.dot(
        t1::TensOrthonormal{order1, dim},
        t2::TensOrthonormal{order2, dim},
    ) where {order1, order2, dim}
    nt1, nt2 = same_basis(t1, t2)
    data = scontract(get_array(nt1), get_array(nt2))
    return Tens(data, get_basis(nt1))
end

function LinearAlgebra.dot(t1::AbstractTens{1, dim}, t2::AbstractTens{1, dim}) where {dim}
    nt1, nt2 = same_basis(t1, t2)
    var = (invvar(get_var(nt1)[end]), get_var(nt2)[(begin + 1):end]...)
    nt2 = change_tens(nt2, get_basis(nt2), var)
    return scontract(get_array(nt1), get_array(nt2))
end

function LinearAlgebra.dot(
        t1::TensOrthonormal{1, dim},
        t2::TensOrthonormal{1, dim},
    ) where {dim}
    nt1, nt2 = same_basis(t1, t2)
    return scontract(get_array(nt1), get_array(nt2))
end

LinearAlgebra.dot(v::TensOrthonormal{1, dim}) where {dim} = dot(get_array(v))

LinearAlgebra.dot(t::TensOrthonormal{2, dim}) where {dim} = dot(get_array(t))

LinearAlgebra.norm(u::AbstractTens{1, dim}) where {dim} = √(dot(u, u))

LinearAlgebra.norm(t::AbstractTens{2, dim}) where {dim} = √(dot(t, t))

"""
    contract(t::AbstractTens{order,dim}, i::Integer, j::Integer)

Calculate the tensor obtained after contraction with respect to the indices `i` and `j`
"""
function contract(t::AbstractTens{order, dim}, i::Integer, j::Integer) where {order, dim}
    var = ntuple(k -> k == j ? invvar(get_var(t, i)) : get_var(t, k), order)
    nt = change_tens(t, get_basis(t), var)
    data = contract(get_array(nt), i, j)
    m = min(i, j)
    M = max(i, j)
    var = (get_var(nt)[1:(m - 1)]..., get_var(nt)[(m + 1):(M - 1)]..., get_var(nt)[(M + 1):order]...)
    return Tens(data, get_basis(nt), var)
end

function contract(t::AbstractTens{2, dim}, i::Integer, j::Integer) where {dim}
    var = ntuple(k -> k == j ? invvar(get_var(t, i)) : get_var(t, k), 2)
    nt = change_tens(t, get_basis(t), var)
    return contract(get_array(nt), i, j)
end

contract(t::TensOrthonormal{order, dim}, i::Integer, j::Integer) where {order, dim} =
    Tens(contract(get_array(t), i, j), get_basis(t))

contract(t::TensOrthonormal{2, dim}, i::Integer, j::Integer) where {dim} =
    contract(get_array(t), i, j)

LinearAlgebra.tr(t::AbstractTens{2}) = contract(t, 1, 2)

"""
    dcontract(t1::AbstractTens{order1,dim}, t2::AbstractTens{order2,dim})

Define a double contracted product between two tensors

`𝛔 ⊡ 𝛆 = σⁱʲεᵢⱼ`
`𝛔 = ℂ ⊡ 𝛆`

# Examples
```julia
julia> 𝛆 = Tens(SymmetricTensor{2,3}((i, j) -> symbols("ε\$i\$j", real = true))) ;

julia> k, μ = symbols("k μ", real =true) ;

julia> ℂ = 3k * t𝕁() + 2μ * t𝕂() ;

julia> 𝛔 = ℂ ⊡ 𝛆
Tens{2, 3, Sym, Sym, SymmetricTensor{2, 3, Sym, 6}, CanonicalBasis{3, Sym}}
# data: 3×3 SymmetricTensor{2, 3, Sym, 6}:
 ε11*(k + 4*μ/3) + ε22*(k - 2*μ/3) + ε33*(k - 2*μ/3)                                              2⋅ε21⋅μ                                              2⋅ε31⋅μ
                                             2⋅ε21⋅μ  ε11*(k - 2*μ/3) + ε22*(k + 4*μ/3) + ε33*(k - 2*μ/3)                                              2⋅ε32⋅μ
                                             2⋅ε31⋅μ                                              2⋅ε32⋅μ  ε11*(k - 2*μ/3) + ε22*(k - 2*μ/3) + ε33*(k + 4*μ/3)
# var: (:cont, :cont)
# basis: 3×3 Tensor{2, 3, Sym, 9}:
 1  0  0
 0  1  0
 0  0  1
```
"""
function Tensors.dcontract(
        t1::AbstractTens{order1, dim},
        t2::AbstractTens{order2, dim},
    ) where {order1, order2, dim}
    nt1, nt2 = same_basis(t1, t2)
    var = (invvar(get_var(nt1)[end - 1]), invvar(get_var(nt1)[end]), get_var(nt2)[(begin + 2):end]...)
    nt2 = change_tens(nt2, get_basis(nt2), var)
    data = Tensors.dcontract(get_array(nt1), get_array(nt2))
    var = (get_var(nt1)[begin:(end - 2)]..., get_var(nt2)[(begin + 2):end]...)
    return Tens(data, get_basis(nt1), var)
end

function Tensors.dcontract(
        t1::TensOrthonormal{order1, dim},
        t2::TensOrthonormal{order2, dim},
    ) where {order1, order2, dim}
    nt1, nt2 = same_basis(t1, t2)
    data = Tensors.dcontract(get_array(nt1), get_array(nt2))
    return Tens(data, get_basis(nt1))
end

"""
    dotdot(v1::AbstractTens{order1,dim}, S::AbstractTens{orderS,dim}, v2::AbstractTens{order2,dim})

Define a bilinear operator `𝐯₁⋅𝕊⋅𝐯₂`

# Examples
```julia
julia> n = Tens(Sym[0, 0, 1]) ;

julia> k, μ = symbols("k μ", real =true) ;

julia> ℂ = 3k * t𝕁() + 2μ * t𝕂() ;

julia> dotdot(n,ℂ,n) # Acoustic tensor
3×3 Tens{2, 3, Sym, Sym, Tensor{2, 3, Sym, 9}, CanonicalBasis{3, Sym}}:
 μ  0          0
 0  μ          0
 0  0  k + 4*μ/3
```
"""
function Tensors.dotdot(
        v1::AbstractTens{order1, dim},
        S::AbstractTens{orderS, dim},
        v2::AbstractTens{order2, dim},
    ) where {order1, orderS, order2, dim}
    nS, nv1 = same_basis(S, v1)
    nS, nv2 = same_basis(S, v2)
    var = (invvar(get_var(nS)[begin]),)
    nv1 = change_tens(nv1, get_basis(nv1), var)
    var = (invvar(get_var(nS)[end]),)
    nv2 = change_tens(nv2, get_basis(nv2), var)
    data = dotdot(nv1.data, nS.data, nv2.data)
    var = (get_var(nS)[begin + 1], get_var(nS)[end - 1])
    return Tens(data, get_basis(nS), var)
end

function Tensors.dotdot(
        v1::TensOrthonormal{order1, dim},
        S::TensOrthonormal{orderS, dim},
        v2::TensOrthonormal{order2, dim},
    ) where {order1, orderS, order2, dim}
    nS, nv1 = same_basis(S, v1)
    nS, nv2 = same_basis(S, v2)
    data = dotdot(nv1.data, nS.data, nv2.data)
    return Tens(data, get_basis(nS))
end

"""
    qcontract(t1::AbstractTens{order1,dim}, t2::AbstractTens{order2,dim})

Define a quadruple contracted product between two tensors

`𝔸 ⊙ 𝔹 = AᵢⱼₖₗBⁱʲᵏˡ`

# Examples
```julia
julia> 𝕀 = t𝕀(Sym) ; 𝕁 = t𝕁(Sym) ; 𝕂 = t𝕂(Sym) ;

julia> 𝕀 ⊙ 𝕀
6

julia> 𝕁 ⊙ 𝕀
1

julia> 𝕂 ⊙ 𝕀
5

julia> 𝕂 ⊙ 𝕁
0
```
"""
function qcontract(
        t1::AbstractTens{order1, dim},
        t2::AbstractTens{order2, dim},
    ) where {order1, order2, dim}
    nt1, nt2 = same_basis(t1, t2)
    var = (
        invvar(get_var(nt1)[end - 3]),
        invvar(get_var(nt1)[end - 2]),
        invvar(get_var(nt1)[end - 1]),
        invvar(get_var(nt1)[end]),
        get_var(nt2)[(begin + 4):end]...,
    )
    nt2 = change_tens(nt2, get_basis(nt2), var)
    data = qcontract(get_array(nt1), get_array(nt2))
    var = (get_var(nt1)[begin:(end - 4)]..., get_var(nt2)[(begin + 4):end]...)
    return Tens(data, get_basis(nt1), var)
end

function qcontract(
        t1::TensOrthonormal{order1, dim},
        t2::TensOrthonormal{order2, dim},
    ) where {order1, order2, dim}
    nt1, nt2 = same_basis(t1, t2)
    data = qcontract(get_array(nt1), get_array(nt2))
    return Tens(data, get_basis(nt1))
end

function qcontract(t1::AbstractTens{4, dim}, t2::AbstractTens{4, dim}) where {dim}
    nt1, nt2 = same_basis(t1, t2)
    var = (
        invvar(get_var(nt1)[end - 3]),
        invvar(get_var(nt1)[end - 2]),
        invvar(get_var(nt1)[end - 1]),
        invvar(get_var(nt1)[end]),
        get_var(nt2)[(begin + 4):end]...,
    )
    nt2 = change_tens(nt2, get_basis(nt2), var)
    return qcontract(get_array(nt1), get_array(nt2))
end

function qcontract(t1::TensOrthonormal{4, dim}, t2::TensOrthonormal{4, dim}) where {dim}
    nt1, nt2 = same_basis(t1, t2)
    return qcontract(get_array(nt1), get_array(nt2))
end

"""
    otimesu(t1::AbstractTens{order1,dim}, t2::AbstractTens{order2,dim})

Define a special tensor product between two tensors of at least second order

`(𝐚 ⊠ 𝐛) ⊡ 𝐩 = 𝐚⋅𝐩⋅𝐛 = aⁱᵏbʲˡpₖₗ eᵢ⊗eⱼ`
"""
function Tensors.otimesu(
        t1::AbstractTens{order1, dim},
        t2::AbstractTens{order2, dim},
    ) where {order1, order2, dim}
    nt1, nt2 = same_basis(t1, t2)
    data = otimesu(get_array(nt1), get_array(nt2))
    var = (
        get_var(nt1)[begin:(end - 1)]...,
        get_var(nt2)[begin],
        get_var(nt1)[end],
        get_var(nt2)[(begin + 1):end]...,
    )
    return Tens(data, get_basis(nt1), var)
end

function Tensors.otimesu(
        t1::TensOrthonormal{order1, dim},
        t2::TensOrthonormal{order2, dim},
    ) where {order1, order2, dim}
    nt1, nt2 = same_basis(t1, t2)
    data = otimesu(get_array(nt1), get_array(nt2))
    return Tens(data, get_basis(nt1))
end

function Tensors.otimesl(
        t1::AbstractTens{order1, dim},
        t2::AbstractTens{order2, dim},
    ) where {order1, order2, dim}
    nt1, nt2 = same_basis(t1, t2)
    data = otimesl(get_array(nt1), get_array(nt2))
    var = (
        get_var(nt1)[begin:(end - 1)]...,
        get_var(nt2)[begin + 1],
        get_var(nt1)[end],
        get_var(nt2)[begin],
        get_var(nt2)[(begin + 2):end]...,
    )
    return Tens(data, get_basis(nt1), var)
end

function Tensors.otimesl(
        t1::TensOrthonormal{order1, dim},
        t2::TensOrthonormal{order2, dim},
    ) where {order1, order2, dim}
    nt1, nt2 = same_basis(t1, t2)
    data = otimesl(get_array(nt1), get_array(nt2))
    return Tens(data, get_basis(nt1))
end

otimesul(S1::SecondOrderTensor{dim}, S2::SecondOrderTensor{dim}) where {dim} =
    symmetric(otimesu(S1, S2))

"""
    otimesul(t1::AbstractTens{order1,dim}, t2::AbstractTens{order2,dim})

Define a special tensor product between two tensors of at least second order

`(𝐚 ⊠ˢ 𝐛) ⊡ 𝐩 = (𝐚 ⊠ 𝐛) ⊡ (𝐩 + ᵗ𝐩)/2  = 1/2(aⁱᵏbʲˡ+aⁱˡbʲᵏ) pₖₗ eᵢ⊗eⱼ`
"""
function otimesul(
        t1::AbstractTens{order1, dim},
        t2::AbstractTens{order2, dim},
    ) where {order1, order2, dim}
    nt1, nt2 = same_basis(t1, t2)
    var = (get_var(nt1)[end - 1], get_var(nt1)[end], get_var(nt2)[(begin + 2):end]...)
    nt2 = change_tens(nt2, get_basis(nt2), var)
    data = otimesul(get_array(nt1), get_array(nt2))
    var = (get_var(nt1)..., get_var(nt2)...)
    return Tens(data, get_basis(nt1), var)
end

function otimesul(
        t1::TensOrthonormal{order1, dim},
        t2::TensOrthonormal{order2, dim},
    ) where {order1, order2, dim}
    nt1, nt2 = same_basis(t1, t2)
    data = otimesul(get_array(nt1), get_array(nt2))
    return Tens(data, get_basis(nt1))
end


"""
    sotimes(t1::AbstractTens{order1,dim}, t2::AbstractTens{order2,dim})

Define a symmetric tensor product between two tensors

`(aⁱeᵢ) ⊗ˢ (bʲeⱼ) = 1/2(aⁱbʲ + aʲbⁱ) eᵢ⊗eⱼ`
"""
function sotimes(
        t1::AbstractTens{order1, dim},
        t2::AbstractTens{order2, dim},
    ) where {order1, order2, dim}
    nt1, nt2 = same_basis(t1, t2)
    var = (get_var(nt1)[end], get_var(nt2)[(begin + 1):end]...)
    nt2 = change_tens(nt2, get_basis(nt2), var)
    data = sotimes(get_array(nt1), get_array(nt2))
    var = (get_var(nt1)..., get_var(nt2)...)
    return Tens(data, get_basis(nt1), var)
end

function sotimes(
        t1::AbstractTens{1, dim},
        t2::AbstractTens{1, dim},
    ) where {dim}
    if isequal(t1, t2)
        return otimes(t1)
    else
        nt1, nt2 = same_basis(t1, t2)
        var = (get_var(nt1)[end], get_var(nt2)[(begin + 1):end]...)
        nt2 = change_tens(nt2, get_basis(nt2), var)
        data = sotimes(get_array(nt1), get_array(nt2))
        var = (get_var(nt1)..., get_var(nt2)...)
        return Tens(data, get_basis(nt1), var)
    end
end

function sotimes(
        t1::TensOrthonormal{order1, dim},
        t2::TensOrthonormal{order2, dim},
    ) where {order1, order2, dim}
    nt1, nt2 = same_basis(t1, t2)
    data = sotimes(get_array(nt1), get_array(nt2))
    return Tens(data, get_basis(nt1))
end

function sotimes(
        t1::TensOrthonormal{1, dim},
        t2::TensOrthonormal{1, dim},
    ) where {dim}
    if isequal(t1, t2)
        return otimes(t1)
    else
        nt1, nt2 = same_basis(t1, t2)
        data = sotimes(get_array(nt1), get_array(nt2))
        return Tens(data, get_basis(nt1))
    end
end

Base.transpose(t::TensArray{order, dim, T, <:SecondOrderTensor}) where {order, dim, T} =
    Tens(transpose(get_array(t)), get_basis(t), (get_var(t)[2], get_var(t)[1]))

Base.transpose(t::TensOrthonormal{order, dim, T, <:SecondOrderTensor}) where {order, dim, T} =
    Tens(transpose(get_array(t)), get_basis(t))

Base.transpose(t::TensArray{order, dim, T, <:FourthOrderTensor}) where {order, dim, T} = Tens(
    Tensors.transpose(get_array(t)),
    get_basis(t),
    (get_var(t)[2], get_var(t)[1], get_var(t)[4], get_var(t)[3]),
)

Base.transpose(t::TensOrthonormal{order, dim, T, <:FourthOrderTensor}) where {order, dim, T} =
    Tens(Tensors.transpose(get_array(t)), get_basis(t))

Tensors.majortranspose(t::TensArray{order, dim, T, <:FourthOrderTensor}) where {order, dim, T} =
    Tens(
    majortranspose(get_array(t)),
    get_basis(t),
    (get_var(t)[3], get_var(t)[4], get_var(t)[1], get_var(t)[2]),
)

Tensors.majortranspose(
    t::TensOrthonormal{order, dim, T, <:FourthOrderTensor},
) where {order, dim, T} = Tens(majortranspose(get_array(t)), get_basis(t))

Tensors.minortranspose(t::TensArray{order, dim, T, <:FourthOrderTensor}) where {order, dim, T} =
    Tens(
    minortranspose(get_array(t)),
    get_basis(t),
    (get_var(t)[2], get_var(t)[1], get_var(t)[4], get_var(t)[3]),
)

Tensors.minortranspose(
    t::TensOrthonormal{order, dim, T, <:FourthOrderTensor},
) where {order, dim, T} = Tens(minortranspose(get_array(t)), get_basis(t))

function tens_basis(ℬ::AbstractBasis{dim, T}, i::Integer, ::Val{:cov}) where {dim, T}
    t = [T(Int(j == i)) for j in 1:dim]
    return Tens(t, ℬ, (:cont,))
end
function tens_basis(ℬ::AbstractBasis{dim, T}, i::Integer, ::Val{:cont}) where {dim, T}
    t = [T(Int(j == i)) for j in 1:dim]
    return Tens(t, ℬ, (:cov,))
end
tens_basis(ℬ::AbstractBasis{dim, T}, i::Integer, var = :cov) where {dim, T} = tens_basis(ℬ, i, Val(var))
tens_basis(ℬ::AbstractBasis{dim, T}, var = :cov) where {dim, T} = ntuple(i -> tens_basis(ℬ, i, Val(var)), get_dim(ℬ))

export AbstractTens, Tens
export proj_tens, best_sym_tens
export get_order, arraytype, get_data, get_array, get_basis, get_var
export intrinsic
export components, components_canon, change_tens, change_tens_canon
export diff_with_basis
export KM, inv_KM
export get_basis, get_var
export tens_basis
