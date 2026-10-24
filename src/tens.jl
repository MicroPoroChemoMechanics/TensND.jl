# ── Tensor types ──────────────────────────────────────────────────────────────
#
# `AbstractTens` and the general tensor types — `Tens` (any basis, with
# variances), `TensRotated`, `TensCanonical`, `TensOrthogonal` — with their
# accessors and the storage of components as Tensors.jl arrays. Displaying them
# is in `tens_display.jl`, their components in other bases in
# `tens_components.jl`, Kelvin-Mandel notation in `kelvin_mandel.jl`, and their
# algebra in `tens_algebra.jl`.

"""
    AbstractTens{order,dim,T}

Supertype of every tensor in `TensND`: a tensor of a given `order` and `dim`,
with element type `T`.

A tensor is **not** an array. It is an array of components *together with* the
basis those components refer to and the variance of each index, which is what
lets tensors expressed in different frames be combined, and what makes comparing
raw component arrays meaningless. See
[Bases and variance](@ref th-bases-variance).

The concrete subtypes fall in two families:

| Family | Types | Stores |
|:--|:--|:--|
| general | [`Tens`](@ref), `TensRotated`, `TensCanonical`, `TensOrthogonal` | a full component array |
| structured | [`TensISO`](@ref), [`TensTI`](@ref), [`TensOrtho`](@ref) | a few scalars plus an orientation |

The type follows the basis: only the non-orthonormal cases carry a variance
tuple. Accessors common to all of them are [`get_order`](@ref),
[`get_dim`](@ref), [`get_basis`](@ref), [`get_var`](@ref), [`get_data`](@ref)
and [`get_array`](@ref).
"""
abstract type AbstractTens{order, dim, T <: Number} <: AbstractArray{T, order} end

"""
    get_order(t::AbstractTens) → Int

Order (number of indices) of `t`.
"""
get_order(::AbstractTens{order, dim, T}) where {order, dim, T} = order
"""
    get_dim(t::AbstractTens) → Int
    get_dim(ℬ::AbstractBasis) → Int

Dimension of the underlying space.
"""
get_dim(::AbstractTens{order, dim, T}) where {order, dim, T} = dim
Base.eltype(::Type{AbstractTens{order, dim, T}}) where {order, dim, T} = T


"""
    Tens{order,dim,T,A<:AbstractArray,B<:Basis}

Tensor type of any order defined by
- a multidata of components (of any type inheriting from `AbstractArray`, e.g. `Tensor` or `SymmetricTensor`)
- a basis of `AbstractBasis` type, whose concrete type is the parameter `B`
- a tuple of variances (covariant `:cov` or contravariant `:cont`) of length equal to the `order` of the tensor

The variance tuple only matters when the basis is not orthonormal; on an
orthonormal basis the covariant and contravariant components coincide, which is
why `TensCanonical`, `TensRotated` and `TensOrthonormal` do not carry one. See [Bases and variance](@ref th-bases-variance).

# Examples

Storing the covariant metric of a non-orthogonal basis as a twice-covariant
tensor, and raising one index — which necessarily returns the identity, since
`gⁱᵏ gₖⱼ = δⁱⱼ`:

```julia
julia> ℬ = Basis(Sym[1 0 0; 0 1 0; 0 1 1]) ;

julia> T = Tens(metric(ℬ, :cov), ℬ, (:cov, :cov))
3×3 Tens{2, 3, Sym, SymmetricTensor{2, 3, Sym, 6}}:
 1  0  0
 0  2  1
 0  1  1

julia> components(T, (:cont, :cov))
3×3 Matrix{Sym}:
 1  0  0
 0  1  0
 0  0  1
```

See also [`components`](@ref), [`change_tens`](@ref), [`get_var`](@ref).
"""
struct Tens{order, dim, T, A <: AbstractArray, B <: Basis} <: AbstractTens{order, dim, T}
    data::A
    # A type parameter rather than an abstract field: every `get_basis`, and
    # everything dispatched on its result, would otherwise be resolved at run
    # time.
    basis::B
    var::NTuple{order, Symbol}
    function Tens(
            data::AbstractArray{T, order},
            basis::Basis{dim},
            var::NTuple{order, Symbol} = ntuple(_ -> :cont, Val(order)),
        ) where {order, dim, T}
        newdata = tensor_or_array(data)
        return new{order, dim, T, typeof(newdata), typeof(basis)}(newdata, basis, var)
    end
    Tens(data::AbstractArray, basis::RotatedBasis, args...) = TensRotated(data, basis)
    Tens(data::AbstractArray, basis::OrthogonalBasis, args...) =
        TensOrthogonal(data, basis, args...)
    Tens(
        data::AbstractArray{T},
        basis::OrthonormalBasis = CanonicalBasis{_dim_of(data), T}(),
        args...,
    ) where {T} = TensOrthonormal(data, basis)
    Tens(data::AbstractArray, var::NTuple, basis::AbstractBasis = CanonicalBasis{_dim_of(data), eltype(data)}()) = Tens(data, basis, var)
    Tens(data::AbstractArray{T, 0}, args...) where {T} = T(data[1])
    # A zero-dimensional array is a wrapped scalar and unwraps to one whatever
    # follows it. The method above says so, but it is generic in `args...` while
    # every constructor above is specific in its second argument: each pair is
    # more specific on one side and less on the other, which is an ambiguity
    # rather than a preference. Naming the crossings settles it; all of them
    # unwrap, exactly as the generic method does.
    Tens(data::AbstractArray{T, 0}, ::Basis, args...) where {T} = T(data[1])
    Tens(data::AbstractArray{T, 0}, ::RotatedBasis, args...) where {T} = T(data[1])
    Tens(data::AbstractArray{T, 0}, ::OrthogonalBasis, args...) where {T} = T(data[1])
    Tens(data::AbstractArray{T, 0}, ::OrthonormalBasis, args...) where {T} = T(data[1])
    Tens(data::AbstractArray{T, 0}, ::NTuple, args...) where {T} = T(data[1])
    # And the two three-argument forms, specific in their third argument where
    # the ones above are generic.
    Tens(data::AbstractArray{T, 0}, ::NTuple, ::AbstractBasis) where {T} = T(data[1])
    Tens(data::AbstractArray{T, 0}, ::Basis, ::NTuple{0, Symbol}) where {T} = T(data[1])
    Tens(data::T, args...) where {T} = data
end


struct TensRotated{order, dim, T, A <: AbstractArray, B <: RotatedBasis} <: AbstractTens{order, dim, T}
    data::A
    basis::B
    function TensRotated(
            data::AbstractArray{T, order},
            basis::RotatedBasis{dim},
        ) where {order, dim, T}
        newdata = tensor_or_array(data)
        return new{order, dim, T, typeof(newdata), typeof(basis)}(newdata, basis)
    end
end

# The dimension, from the type of the array when it carries one (every
# Tensors.jl array) so that the type of the tensor is inferred, from its size
# otherwise.
_dim_of(::Tensors.AbstractTensor{order, dim}) where {order, dim} = dim
_dim_of(a::AbstractArray) = size(a, 1)

struct TensCanonical{order, dim, T, A <: AbstractArray} <: AbstractTens{order, dim, T}
    data::A
    function TensCanonical(data::AbstractArray{T, order}) where {order, T}
        newdata = tensor_or_array(data)
        return new{order, _dim_of(newdata), T, typeof(newdata)}(newdata)
    end
end

struct TensOrthogonal{order, dim, T, A <: AbstractArray, B <: OrthogonalBasis} <: AbstractTens{order, dim, T}
    data::A
    basis::B
    var::NTuple{order, Symbol}
    function TensOrthogonal(
            data::AbstractArray{T, order},
            basis::OrthogonalBasis{dim},
            var::NTuple{order, Symbol} = ntuple(_ -> :cont, Val(order)),
        ) where {order, dim, T}
        newdata = tensor_or_array(data)
        return new{order, dim, T, typeof(newdata), typeof(basis)}(newdata, basis, var)
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
    ℬ = get_basis(t)
    m = _narrow_eltype(_components_eltype(T, ℬ), Array(get_array(t)))
    _scale_by_variance!(m, get_var(t), inv.(ℬ.λ), ℬ.λ)
    return Tens(m, relevant_OrthonormalBasis(ℬ))
end

# `m[i₁, …, iₙ] *= s₁[i₁] ⋯ sₙ[iₙ]`, where `sₖ` is `cov` if the k-th index is
# covariant and `cont` if it is contravariant: the scale factors between the
# components in an orthogonal basis and in its orthonormal counterpart.
function _scale_by_variance!(m::AbstractArray{<:Any, order}, var::NTuple{order, Symbol}, cov, cont) where {order}
    s = ntuple(k -> _by_variance(_pick_by_variance, var[k], cov, cont), Val(order))
    for ind in CartesianIndices(m)
        m[ind] *= prod(ntuple(k -> s[k][ind[k]], Val(order)))
    end
    return m
end
_pick_by_variance(cov, cont, ::Val{:cov}) = cov
_pick_by_variance(cov, cont, ::Val{:cont}) = cont

@inline Base.size(t::TensArray) = size(get_array(t))
@inline Base.getindex(t::TensArray, ind...) = getindex(get_array(t), ind...)

@inline Base.zero(t::AbstractTens) = Tens(zero.(get_array(t)), get_basis(t), get_var(t))

# This function aims at storing the table of components in the `Tensor` type whenever possible
# Convert a raw array into the best matching Tensors.jl type (Vec, Tensor, or
# SymmetricTensor).  Falls through to the input unchanged for types that are
# already Tensors.AllTensors, for orders not in {1,2,4}, and for dimensions
# above 3, where Tensors.jl has no type.
tensor_or_array(tab::AbstractArray{T, 1}) where {T} = size(tab, 1) ≤ 3 ? Vec{size(tab, 1)}(tab) : tab
# `dim` comes from `size(tab, 1)`, i.e. it is a RUNTIME value.  Writing
# `Tensor{order, dim}(tab)` with such a `dim` builds a non-concrete type at
# compile time, so the construction is fully dynamic: measured at 3147 ns and
# 3120 bytes for an 81-element array — versus 311 ns to produce that array in
# the first place.  Since every structured tensor reaches the generic route
# through this function (`change_tens` → `same_basis` → every binary
# operation), that cost was paid twice per `⊡` between structured operands.
#
# Funnel the runtime `dim` through `Val` so each branch constructs a concrete
# type.  Only 2 and 3 are specialized — the dimensions this package actually
# uses; anything else falls back to the previous dynamic path.
@inline function _to_tensor(tab, ::Val{dim}, ::Val{order}) where {dim, order}
    newtab = Tensor{order, dim}(tab)
    return _store_symmetric(newtab) ? Tensors.symmetric(newtab) : newtab
end

"""
    _store_symmetric(t) -> Bool

Whether a freshly built `Tensor` should be kept in its symmetric form.

`Tensors.issymmetric` compares components exactly, which makes the storage
decision hinge on the last bit. Writing a structured tensor (`TensISO`,
`TensTI`, `TensOrtho`) about a non-canonical axis leaves minor-antisymmetric
residue of order `1e-16`: a tensor whose components are of size 30 then lands
on the 81-component `Tensor` where the very same tensor written about `e₃`
lands on the 36-component `SymmetricTensor`.

The two are not interchangeable downstream. `inv` of the full fourth-order
form solves a 9×9 system whose three minor-antisymmetric directions are null,
so it throws `SingularException` — a dilute localization `inv(𝕀 + ℙ : δℂ)`
would fail or succeed depending on the last bit of an intermediate. Residue at
the scale of the tensor's own rounding error is noise, not antisymmetry, and
is treated as such.

Truly exact element types (symbolic, rational, integer) keep the exact
comparison: there is no round-off to absorb, and no tolerance to define. The
split is [`ApproxType`](@ref) — **not** `AbstractFloat`: a `ForwardDiff.Dual`
subtypes `Real` without subtyping `AbstractFloat`, yet it carries exactly the
same round-off as the float it is built on, so it belongs on the tolerant side.
Keying on `AbstractFloat` made automatic differentiation impossible through any
structured tensor written about a non-canonical axis: `mandel` of a `Dual`-valued
`TensISO` in a `RotatedBasis` returned a 9×9 matrix where the same tensor in
`Float64` returned a 6×6 one.
"""
@inline _store_symmetric(t::Tensors.AbstractTensor) = _issym(t)

@inline function _store_symmetric(
        t::Tensors.AbstractTensor{order, dim, T}
    ) where {order, dim, T <: ApproxType}
    _issym(t) && return true
    scale = maximum(abs, t)
    iszero(scale) && return true
    return maximum(abs, t - Tensor{order, dim}(Tensors.symmetric(t))) <= 16 * _approx_eps(T) * scale
end

"""
    _approx_eps(::Type{T}) -> AbstractFloat

Machine epsilon of the floating-point type underlying `T`, for the tolerance of
[`_store_symmetric`](@ref).

`eps` is already defined on a `ForwardDiff.Dual` type (it forwards to the value
type), but not on a `Complex` one — hence this thin dispatcher rather than a
bare `eps(T)`.
"""
@inline _approx_eps(::Type{T}) where {T <: Complex} = eps(real(T))
@inline _approx_eps(::Type{T}) where {T} = eps(T)

for order in (2, 4)
    @eval function tensor_or_array(tab::AbstractArray{T, $order}) where {T}
        dim = size(tab, 1)
        dim == 3 && return _to_tensor(tab, Val(3), Val($order))
        dim == 2 && return _to_tensor(tab, Val(2), Val($order))
        dim == 1 || return tab
        newtab = Tensor{$order, dim}(tab)
        return _store_symmetric(newtab) ? Tensors.symmetric(newtab) : newtab
    end
    @eval tensor_or_array(tab::Tensor{$order, dim}) where {dim} = _store_symmetric(tab) ? Tensors.symmetric(tab) : tab
end
tensor_or_array(tab::Tensors.AllTensors) = tab
tensor_or_array(tab::AbstractArray) = tab

# `tensor_or_array` for an array whose order and dimension the caller knows from
# a type: the Tensors.jl type is then built at a static dimension, with the
# same storage decision, where `tensor_or_array` reads the dimension from the
# array at run time and leaves every caller with an uninferred result.
@inline _static_tensor(a, ::Val{dim}, ::Val{order}) where {dim, order} =
    order in (2, 4) && dim in (2, 3) ? _to_tensor(a, Val(dim), Val(order)) : tensor_or_array(a)
@inline _static_tensor(t::AbstractTens{order, dim}) where {order, dim} =
    _static_tensor(get_array(t), Val(dim), Val(order))

##############################
# Utility/Accessor Functions #
##############################

get_array(t::TensArray) = t.data
"""
    get_basis(t::AbstractTens) → AbstractBasis

Basis the components of `t` refer to. Structured types report the canonical
basis, their orientation being carried separately by [`axis`](@ref) or
[`frame`](@ref).
"""
get_basis(t::TensBasis) = t.basis
get_basis(::TensCanonical{order, dim, T}) where {order, dim, T} = CanonicalBasis{dim, T}()
"""
    get_var(t::AbstractTens) → NTuple{order,Symbol}
    get_var(t::AbstractTens, i::Integer) → Symbol

Variance of each index of `t`, `:cov` (covariant) or `:cont` (contravariant).

On an orthonormal basis the two coincide, so tensors stored on one report a
constant tuple and the distinction may be ignored. See
[Bases and variance](@ref th-bases-variance).
"""
get_var(::TensOrthonormal{order}) where {order} = ntuple(_ -> :cov, Val(order))
get_var(::TensOrthonormal, i::Integer) = :cov
get_var(t::TensVar) = t.var
get_var(t::TensVar, i::Integer) = t.var[i]


export AbstractTens, Tens
export get_order, get_data, get_array, get_basis, get_var
