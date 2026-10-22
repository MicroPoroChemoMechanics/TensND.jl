# ──────────────────────────────────────────────────────────────────────────────
# Shared operations for structured tensor types (TensISO, TensTI{4}, TensTI,
# TensOrtho, TensCubic).  Each type must implement:
#   • get_data(t)        → NTuple of scalar coefficients
#   • _rebuild(t, data) → new tensor of the same kind with updated data
# ──────────────────────────────────────────────────────────────────────────────

# ── Symbolic helpers (generic T) ─────────────────────────────────────────────
# A single definition per operation, dispatched on the union of the
# structured types.

for OP in (:tsimplify, :tfactor, :tsubs, :tdiff, :ttrigsimp, :texpand_trig, :tlimit)
    @eval $OP(A::Union{TensISO, TensTI{4}, TensTI, TensOrtho, TensCubic}, args...; kwargs...) =
        _rebuild(A, $OP(get_data(A), args...; kwargs...))
end

# ── Explicit Num dispatch ────────────────────────────────────────────────────
# Resolves method ambiguity with the AbstractArray{Num} methods defined in
# array_utils.jl for Symbolics.jl support.  The element-type parameter (Num)
# sits at a different position in each type, so we enumerate them explicitly.

for OP in (:tsimplify, :tsubs, :tdiff)
    @eval $OP(A::TensISO{order, dim, Num}, args...; kwargs...) where {order, dim} =
        _rebuild(A, $OP(get_data(A), args...; kwargs...))
    @eval $OP(A::TensTI{4, Num, N}, args...; kwargs...) where {N} =
        _rebuild(A, $OP(get_data(A), args...; kwargs...))
    @eval $OP(A::TensTI{order, Num, N}, args...; kwargs...) where {order, N} =
        _rebuild(A, $OP(get_data(A), args...; kwargs...))
    @eval $OP(A::TensOrtho{Num}, args...; kwargs...) =
        _rebuild(A, $OP(get_data(A), args...; kwargs...))
    @eval $OP(A::TensCubic{Num}, args...; kwargs...) =
        _rebuild(A, $OP(get_data(A), args...; kwargs...))
end

# ── Scalar arithmetic ────────────────────────────────────────────────────────
# Negation, scalar multiplication, and scalar division are identical across all
# structured tensor types: apply the operation element-wise on get_data, then
# _rebuild.  Type promotion is handled correctly by _rebuild (which uses
# eltype(new_data)), supporting generic Number types including ForwardDiff.Dual.

for ST in (TensISO, TensTI{4}, TensTI, TensOrtho, TensCubic)
    @eval @inline Base.:-(A::$ST) = _rebuild(A, .-(get_data(A)))
    @eval @inline Base.:*(α::Number, A::$ST) = _rebuild(A, α .* get_data(A))
    @eval @inline Base.:*(A::$ST, α::Number) = _rebuild(A, get_data(A) .* α)
    @eval @inline Base.:/(A::$ST, α::Number) = _rebuild(A, get_data(A) ./ α)
end

# `Symbolics.Num <: Real`, and Symbolics divides any real array by a `Num`
# (`/(::AbstractArray{<:Real}, ::Num)`). A structured tensor with real entries
# divided by a `Num` matches that method and the one above equally well, so the
# call was refused as ambiguous; the structured division is restated on the
# intersection of the two signatures, which keeps the result structured.
for ST in (TensISO, TensTI{4}, TensTI, TensOrtho, TensCubic)
    @eval @inline Base.:/(A::$(typeintersect(ST, AbstractArray{<:Real})), α::Symbolics.Num) =
        _rebuild(A, get_data(A) ./ α)
end

# ── Reference checks for binary operations ───────────────────────────────────
# TensTI{4} and TensTI carry a symmetry axis (n); TensOrtho carries a material
# frame; TensCubic carries a cube frame, and its check is looser than equality
# for a reason given in `_same_cube_frame`.  Binary arithmetic requires the
# same reference.

# An `ArgumentError` rather than `@assert`, which may be compiled out.
@inline _check_same_reference(A::TensTI{4}, B::TensTI{4}) =
    axis(A) == axis(B) || throw(ArgumentError("TensTI{4} operation requires the same axis"))
@inline _check_same_reference(A::TensTI, B::TensTI) =
    axis(A) == axis(B) || throw(ArgumentError("TensTI operation requires the same axis"))
@inline _check_same_reference(A::TensOrtho, B::TensOrtho) =
    frame(A) == frame(B) ||
    throw(ArgumentError("TensOrtho operation requires the same material frame"))

# ── Tensors.jl form of a fourth-order structured tensor ──────────────────────
# Minor-symmetric by construction, so the `SymmetricTensor` is built from its 36
# independent components in closed form (`getindex`), instead of all 81 into an
# array whose symmetry `_store_symmetric` then checks: the same storage,
# `change_tens`, `⊡` and every path through `_static_tensor` without the detour.
@inline _static_tensor(t::Union{TensTI{4}, TensOrtho, TensCubic}) =
    SymmetricTensor{4, 3, eltype(t)}((i, j, k, l) -> @inbounds t[i, j, k, l])

# An isotropic tensor of order 2 or 4 is symmetric by construction too: built at
# the static dimension of its type, where the array route read the dimension at
# run time and left `change_tens` uninferred.
@inline _static_tensor(t::TensISO{order, dim}) where {order, dim} =
    order in (2, 4) && dim in (2, 3) ?
    SymmetricTensor{order, dim, eltype(t)}((ij...) -> @inbounds t[ij...]) :
    _static_tensor(get_array(t), Val(dim), Val(order))

# ── Real and imaginary parts ──────────────────────────────────────────────────

"""
    real(t::AbstractTens)
    imag(t::AbstractTens)

The real or imaginary part of the tensor `t`, in its class: a `TensISO`,
`TensTI`, `TensOrtho` or `TensCubic` stays one (its coefficients are linear in
the tensor, its axis or frame a real geometry), a tensor in a real basis stays in
that basis. In a complex basis the real part of the components is not the
components of the real part, so such a tensor comes back in the canonical basis.
"""
Base.real(t::AbstractTens) = _part(real, t)
Base.imag(t::AbstractTens) = _part(imag, t)
"""
    conj(t::AbstractTens)

The complex conjugate of the tensor `t`, in its class, as `real` and
`imag` are.
"""
Base.conj(t::AbstractTens) = _part(conj, t)

_part(f, t::Union{TensISO, TensOrtho, TensCubic}) = _rebuild(t, map(f, get_data(t)))
function _part(f, t::TensTI{order, T, N}) where {order, T, N}
    d = map(f, get_data(t))
    return TensTI{order, eltype(d), N}(d, map(x -> convert(eltype(d), real(x)), axis(t)))
end
_part(f, t::AbstractTens) = _part(f, t, get_basis(t))
_part(f, t::AbstractTens, ℬ::Union{CanonicalBasis, AbstractBasis{<:Any, <:Real}}) =
    Tens(f.(Array(get_array(t))), ℬ, get_var(t))
_part(f, t::AbstractTens{order, dim}, ::AbstractBasis) where {order, dim} =
    Tens(f.(Array(components_canon(t))), CanonicalBasis{dim, real(eltype(t))}())

# ── Coefficients: rebuild, norms of the coefficient basis, flat vectors ──────

"""
    rebuild(t, data)

The tensor of the class of `t` — same order, same axis or frame — whose
coefficients ([`get_data`](@ref), in the same order) are `data`, a tuple or a
vector. The coefficients are promoted to a common type; a transversely isotropic
tensor also promotes its axis, which shares their type. A `TensTI` given six or
eight coefficients becomes the six- or eight-parameter form.

```julia
julia> rebuild(TensISO{3}(3.0, 2.0), (1.0, 4.0)) == TensISO{3}(1.0, 4.0)
true
```
"""
rebuild(t::Union{TensISO, TensOrtho, TensCubic}, data) = _rebuild(t, _common_tuple(data))
function rebuild(t::TensTI, data)
    d = _common_tuple(data)
    T = promote_type(eltype(d), eltype(t))
    return _rebuild(t, map(x -> convert(T, x), d))
end

function _common_tuple(data)
    d = Tuple(data)
    T = mapreduce(typeof, promote_type, d)
    return map(x -> convert(T, x), d)
end

"""
    param_norms(t) → NTuple

The Frobenius norms of the tensors the coefficients of `t` multiply,
`‖rebuild(t, eᵢ)‖`: `(1, √5)` for the isotropic `(𝕁, 𝕂)` in 3-D, the norms of
the Walpole, orthotropic or cubic basis tensors for the other classes. These
bases are orthogonal, so `norm(t) == norm(param_norms(t) .* get_data(t))`: the
weights that make the Euclidean norm of the coefficients the norm of the tensor.
For a generic tensor, whose [`to_vector`](@ref) coordinates are orthonormal
Kelvin-Mandel components, they are all one.
"""
function param_norms(t::Union{TensISO, TensTI, TensOrtho, TensCubic})
    d = get_data(t)
    N = length(d)
    o, z = one(eltype(d)), zero(eltype(d))
    return ntuple(i -> norm(rebuild(t, ntuple(j -> j == i ? o : z, N))), N)
end
# The coordinates of `to_vector` of a generic tensor are Kelvin-Mandel components,
# orthonormal: every weight is one.
param_norms(t::AbstractTens) = ones(primal_type(real(eltype(t))), length(to_vector(t)))

"""
    to_vector(t) → Vector
    from_vector(prototype, v)

A tensor as a flat vector of real coordinates, and back: its coefficients
([`get_data`](@ref)) for a `TensISO`, `TensTI`, `TensOrtho` or `TensCubic`, the
Kelvin-Mandel components in the canonical basis otherwise (6 or 9 at order 2, 36
or 81 at order 4). `from_vector` rebuilds a tensor of the class of `prototype`,
with its axis or frame; a generic tensor comes back in the canonical basis. Both
are linear, for an ODE state, a nonlinear solver or a sequence accelerator.
"""
to_vector(t::Union{TensISO, TensTI, TensOrtho, TensCubic}) = collect(get_data(t))
to_vector(t::AbstractTens{order, dim}) where {order, dim} = vec(collect(KM(t, CanonicalBasis{dim, eltype(t)}())))

"""
    from_vector(prototype, v)

The tensor whose flat vector ([`to_vector`](@ref)) is `v`: of the class of
`prototype`, with its axis or frame, or in the canonical basis for a generic one.
"""
from_vector(prototype::Union{TensISO, TensTI, TensOrtho, TensCubic}, v) = rebuild(prototype, v)
from_vector(::AbstractTens{2}, v) = inv_KM(collect(v))
from_vector(::AbstractTens{4}, v) = (n = isqrt(length(v)); inv_KM(reshape(collect(v), n, n)))

export rebuild, param_norms, to_vector, from_vector

# A structured tensor is written with its canonical components already (its
# axis or frame is carried apart): rewriting it in the canonical basis is the
# identity, and keeps its class.
change_tens_canon(t::Union{TensISO, TensTI, TensOrtho, TensCubic}) = t

export generic_tens
