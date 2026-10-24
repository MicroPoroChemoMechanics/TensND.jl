# ── Components in another basis, with other variances ────────────────────────

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

"""
    _components_eltype(T, bases...) -> Type

Element type of the components of a tensor of element type `T` written through
the given bases: wide enough for a `Float64` tensor in a basis of duals or
symbols, which converting back to `T` cannot hold. The conversion that follows
each contraction narrows the output of the contraction to this type.

The canonical basis imposes no type: its entries are exact zeros and ones,
which every number type represents. Without that exception, the default
`CanonicalBasis()` — of element type `Sym` — would turn every `Float64`
tensor written in it into a symbolic one.
"""
_components_eltype(::Type{T}, bases::AbstractBasis...) where {T} =
    promote_type(T, map(_entry_type, bases)...)

_entry_type(::CanonicalBasis) = Union{}
_entry_type(::AbstractBasis{dim, T}) where {dim, T} = T

function components(t::Tens{order, dim, T}, var::NTuple{order, Symbol}) where {order, dim, T}
    if isequal(var, get_var(t))
        return get_array(t)
    else
        Tc = _components_eltype(T, get_basis(t))
        m = _narrow_eltype(Tc, Array(get_array(t)))
        for i in 1:order
            if !isequal(get_var(t, i), var[i])
                m = _narrow_eltype(Tc, _contract_index(m, metric(get_basis(t), var[i]), i))
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
        Tc = _components_eltype(T, get_basis(t), ℬ)
        m = _narrow_eltype(Tc, Array(get_array(t)))
        for i in 1:order
            # The matrix of the change for this index, from the variance it has
            # to the one asked: computed here rather than stored in a
            # `Dict{…, AbstractMatrix}`, whose abstract values every use read
            # at run time.
            c = transpose(vecbasis(get_basis(t), invvar(get_var(t, i)))) * vecbasis(ℬ, var[i])
            if c ≠ 1I
                m = _narrow_eltype(Tc, _contract_index(m, c, i))
            end
        end
        return m
    end
end

# A tensor whose own basis is `ℬ` is a `TensOrthogonal`, served above: here
# the basis always changes.
function components(
        t::AbstractTens{order, dim, T},
        ℬ::OrthogonalBasis{dim},
        var::NTuple{order, Symbol},
    ) where {order, dim, T}
    m0 = components(t, relevant_OrthonormalBasis(ℬ))
    m = _narrow_eltype(_components_eltype(eltype(m0), ℬ), Array(m0))
    return _scale_by_variance!(m, var, ℬ.λ, inv.(ℬ.λ))
end

components(t::AbstractTens{order, dim, T}, ℬ::AbstractBasis{dim}) where {order, dim, T} =
    components(t, ℬ, get_var(t))

components(t::AbstractTens{order, dim, T}, ℬ::OrthonormalBasis{dim}) where {order, dim, T} =
    components(t, ℬ, ntuple(_ -> :cont, Val(order)))

function components(
        t::TensOrthonormal{order, dim, T},
        ℬ::OrthonormalBasis{dim},
    ) where {order, dim, T}
    isequal(ℬ, get_basis(t)) && return get_array(t)
    return _as_array(_components_in_orthonormal(get_array(t), t, ℬ))
end

# The components of `t` in another orthonormal basis ℬ. With `R` the transpose
# of the change matrix `bb` (`bb[i, a] = eᵢ ⋅ e′ₐ`), they are `R ⋅ x ⋅ Rᵀ` at
# order 2 and `(R ⊠ R) ⊡ x ⊡ (Rᵀ ⊠ Rᵀ)` at order 4: Tensors.jl's static
# operations, about seven times faster than four contractions of an `Array`,
# and of a type known at compile time. The orders Tensors.jl has no type for
# keep the contractions.
function _components_in_orthonormal(
        x::Union{Vec{dim}, Tensor{2, dim}, SymmetricTensor{2, dim}, Tensor{4, dim}, SymmetricTensor{4, dim}},
        t::TensOrthonormal{order, dim, T},
        ℬ::OrthonormalBasis{dim},
    ) where {order, dim, T}
    Tc = _components_eltype(T, get_basis(t), ℬ)
    bb = transpose(vecbasis(get_basis(t))) * vecbasis(ℬ)
    return _rotate_components(x, Tensor{2, dim, Tc}(transpose(bb)))
end
function _components_in_orthonormal(
        x::AbstractArray, t::TensOrthonormal{order, dim, T}, ℬ::OrthonormalBasis{dim}
    ) where {order, dim, T}
    bb = transpose(vecbasis(get_basis(t))) * vecbasis(ℬ)
    Tc = _components_eltype(T, get_basis(t), ℬ)
    m = _narrow_eltype(Tc, Array(x))
    if bb ≠ 1I
        for i in 1:order
            m = _narrow_eltype(Tc, _contract_index(m, bb, i))
        end
    end
    return m
end

_rotate_components(x::Vec, R) = R ⋅ x
_rotate_components(x::Tensor{2}, R) = R ⋅ x ⋅ transpose(R)
_rotate_components(x::SymmetricTensor{2}, R) = Tensors.unsafe_symmetric(R ⋅ x ⋅ transpose(R))
_rotate_components(x::Tensor{4}, R) = otimesu(R, R) ⊡ x ⊡ otimesu(transpose(R), transpose(R))
_rotate_components(x::SymmetricTensor{4}, R) =
    Tensors.unsafe_symmetric(otimesu(R, R) ⊡ x ⊡ otimesu(transpose(R), transpose(R)))

# `components` has always returned a plain `Array` after a change of basis.
_as_array(x::Tensors.AbstractTensor) = Array(x)
_as_array(x::AbstractArray) = x

# The same change, kept in its Tensors.jl form: the type of the new tensor is
# then that of its storage, known at compile time.
function change_tens(t::TensOrthonormal{order, dim}, ℬ::OrthonormalBasis{dim}) where {order, dim}
    isequal(ℬ, get_basis(t)) && return t
    return Tens(_components_in_orthonormal(get_array(t), t, ℬ), ℬ)
end
# Covariant and contravariant components coincide in an orthonormal basis.
change_tens(t::TensOrthonormal{order, dim}, ℬ::OrthonormalBasis{dim}, ::NTuple{order, Symbol}) where {order, dim} =
    change_tens(t, ℬ)

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

Rewrite the same tensor with components corresponding to the canonical basis.
A `TensISO`, `TensTI`, `TensOrtho` or `TensCubic` has canonical components
already and is returned as it is, in its class.

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


for OP in (:(tsimplify), :(tfactor), :(tsubs), :(ttrigsimp), :(texpand_trig), :(tlimit))
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
"""
    diff_with_basis(t::AbstractTens, args...; kwargs...)

Differentiate `t` **including its basis**, then re-express the result on the
original basis and variance.

This differs from [`∂`](@ref) applied componentwise: when the basis itself
depends on the differentiation variable — as a rotated or curvilinear frame does
— the derivative of the basis vectors contributes, and dropping it is the
classic way to lose the Christoffel terms.
"""
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


# ── The common basis of two operands, and equality ───────────────────────────

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

# Two tensors in the canonical basis are in the same basis already, whatever
# the element types of their two `CanonicalBasis` objects. The generic method
# would find so at run time, through `isequal`, but return a type it cannot
# infer, and every product of a dual tensor with a plain one paid a dynamic
# dispatch for it.
same_basis(t1::TensCanonical, t2::TensCanonical) = (t1, t2)

function same_basis_same_var(
        t1::AbstractTens{order1, dim},
        t2::AbstractTens{order2, dim},
    ) where {order1, order2, dim}
    ℬ = choose_best_basis(get_basis(t1), get_basis(t2))
    return change_tens(t1, ℬ, get_var(t1)), change_tens(t2, ℬ, get_var(t1))
end


# `==` compares tensors as mathematical objects, whatever the basis and the
# variances they are written in. `isequal` is structural — same variances, same
# basis, same components — so that it agrees with `hash`, as `Dict` and `Set`
# require: the basis-independent `isequal` held two tensors equal whose hashes,
# those of their components, differed.
for OP in (:(==), :(!=))
    @eval @inline function Base.$OP(
            t1::AbstractTens{order, dim},
            t2::AbstractTens{order, dim},
        ) where {order, dim}
        nt1, nt2 = same_basis_same_var(t1, t2)
        return $OP(get_array(nt1), get_array(nt2))
    end
end

Base.isequal(t1::AbstractTens, t2::AbstractTens) =
    isequal(_var_key(t1), _var_key(t2)) && isequal(get_basis(t1), get_basis(t2)) &&
    isequal(get_array(t1), get_array(t2))

Base.hash(t::AbstractTens, h::UInt) =
    hash(get_array(t), hash(get_basis(t), hash(_var_key(t), hash(0x5d3c2b7e1a4f6089, h))))

# In an orthonormal basis covariant and contravariant components coincide: the
# variances are then no part of what the tensor is.
_var_key(t::AbstractTens) = _var_key(get_basis(t), get_var(t))
_var_key(::OrthonormalBasis, var) = nothing
_var_key(::AbstractBasis, var) = var

tens_basis(ℬ::AbstractBasis{dim, T}, i::Integer, ::Val{:cov}) where {dim, T} =
    Tens(_unit_components(Val(dim), T, i), ℬ, (:cont,))
tens_basis(ℬ::AbstractBasis{dim, T}, i::Integer, ::Val{:cont}) where {dim, T} =
    Tens(_unit_components(Val(dim), T, i), ℬ, (:cov,))

# The components of the `i`-th basis vector, in the storage `tensor_or_array`
# would choose (a `Vec` up to dimension 3) but built with their static size:
# from a `Vector`, the size is a run-time value, and so was the type of every
# basis vector, which cost a type instantiation and a dynamic dispatch per call.
_unit_components(::Val{dim}, ::Type{T}, i::Integer) where {dim, T} =
    dim ≤ 3 ? Vec{dim, T}(ntuple(j -> T(Int(j == i)), Val(dim))) : [T(Int(j == i)) for j in 1:dim]
"""
    tens_basis(ℬ::AbstractBasis, i::Integer, var = :cov) → AbstractTens{1}
    tens_basis(ℬ::AbstractBasis, var = :cov) → NTuple

The `i`-th basis vector of `ℬ` as a first-order tensor, or all of them as a
tuple. `var = :cov` gives the basis vectors `𝐞ᵢ`, `var = :cont` the dual
vectors `𝐞ⁱ`.
"""
tens_basis(ℬ::AbstractBasis, i::Integer, var::Symbol = :cov) = _by_variance(tens_basis, var, ℬ, i)
tens_basis(ℬ::AbstractBasis, var::Symbol = :cov) = _by_variance(tens_basis, var, ℬ)
tens_basis(ℬ::AbstractBasis{dim}, var::Val) where {dim} = ntuple(i -> tens_basis(ℬ, i, var), Val(dim))

export components, components_canon, change_tens, change_tens_canon
export diff_with_basis
export tens_basis
