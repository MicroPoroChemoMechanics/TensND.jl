# ──────────────────────────────────────────────────────────────────────────────
# Shared operations for structured tensor types (TensISO, TensTI{4}, TensTI,
# TensOrtho).  Each type must implement:
#   • get_data(t)        → NTuple of scalar coefficients
#   • _rebuild(t, data) → new tensor of the same kind with updated data
# ──────────────────────────────────────────────────────────────────────────────

# ── Symbolic helpers (generic T) ─────────────────────────────────────────────
# A single definition per operation, dispatched via Union, replacing the
# identical loops previously duplicated in each file.

for OP in (:tsimplify, :tfactor, :tsubs, :tdiff, :ttrigsimp, :texpand_trig)
    @eval $OP(A::Union{TensISO, TensTI{4}, TensTI, TensOrtho}, args...; kwargs...) =
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
end

# ── Scalar arithmetic ────────────────────────────────────────────────────────
# Negation, scalar multiplication, and scalar division are identical across all
# structured tensor types: apply the operation element-wise on get_data, then
# _rebuild.  Type promotion is handled correctly by _rebuild (which uses
# eltype(new_data)), supporting generic Number types including ForwardDiff.Dual.

for ST in (TensISO, TensTI{4}, TensTI, TensOrtho)
    @eval @inline Base.:-(A::$ST) = _rebuild(A, .-(get_data(A)))
    @eval @inline Base.:*(α::Number, A::$ST) = _rebuild(A, α .* get_data(A))
    @eval @inline Base.:*(A::$ST, α::Number) = _rebuild(A, get_data(A) .* α)
    @eval @inline Base.:/(A::$ST, α::Number) = _rebuild(A, get_data(A) ./ α)
end

# ── Reference checks for binary operations ───────────────────────────────────
# TensTI{4} and TensTI carry a symmetry axis (n); TensOrtho carries a material
# frame.  Binary arithmetic requires the same reference.

@inline _check_same_reference(A::TensTI{4}, B::TensTI{4}) =
    @assert axis(A) == axis(B) "TensTI{4} operation requires the same axis"
@inline _check_same_reference(A::TensTI, B::TensTI) =
    @assert axis(A) == axis(B) "TensTI operation requires the same axis"
@inline _check_same_reference(A::TensOrtho, B::TensOrtho) =
    @assert frame(A) == frame(B) "TensOrtho operation requires the same material frame"
