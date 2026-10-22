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
