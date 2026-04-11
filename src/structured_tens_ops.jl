# ──────────────────────────────────────────────────────────────────────────────
# Shared operations for structured tensor types (TensISO, TensWalpole, TensTI,
# TensOrtho).  Each type must implement:
#   • getdata(t)        → NTuple of scalar coefficients
#   • _rebuild(t, data) → new tensor of the same kind with updated data
# ──────────────────────────────────────────────────────────────────────────────

# ── Symbolic helpers (generic T) ─────────────────────────────────────────────
# A single definition per operation, dispatched via Union, replacing the
# identical loops previously duplicated in each file.

for OP in (:tsimplify, :tfactor, :tsubs, :tdiff, :ttrigsimp, :texpand_trig)
    @eval $OP(A::Union{TensISO,TensWalpole,TensTI,TensOrtho}, args...; kwargs...) =
        _rebuild(A, $OP(getdata(A), args...; kwargs...))
end

# ── Explicit Num dispatch ────────────────────────────────────────────────────
# Resolves method ambiguity with the AbstractArray{Num} methods defined in
# array_utils.jl for Symbolics.jl support.  The element-type parameter (Num)
# sits at a different position in each type, so we enumerate them explicitly.

for OP in (:tsimplify, :tsubs, :tdiff)
    @eval $OP(A::TensISO{order,dim,Num}, args...; kwargs...) where {order,dim} =
        _rebuild(A, $OP(getdata(A), args...; kwargs...))
    @eval $OP(A::TensWalpole{Num,N}, args...; kwargs...) where {N} =
        _rebuild(A, $OP(getdata(A), args...; kwargs...))
    @eval $OP(A::TensTI{order,Num,N}, args...; kwargs...) where {order,N} =
        _rebuild(A, $OP(getdata(A), args...; kwargs...))
    @eval $OP(A::TensOrtho{Num}, args...; kwargs...) =
        _rebuild(A, $OP(getdata(A), args...; kwargs...))
end

# ── Scalar arithmetic ────────────────────────────────────────────────────────
# Negation, scalar multiplication, and scalar division are identical across all
# structured tensor types: apply the operation element-wise on getdata, then
# _rebuild.  Type promotion is handled correctly by _rebuild (which uses
# eltype(new_data)), supporting generic Number types including ForwardDiff.Dual.

for ST in (TensISO, TensWalpole, TensTI, TensOrtho)
    @eval @inline Base.:-(A::$ST) = _rebuild(A, .-(getdata(A)))
    @eval @inline Base.:*(α::Number, A::$ST) = _rebuild(A, α .* getdata(A))
    @eval @inline Base.:*(A::$ST, α::Number) = _rebuild(A, getdata(A) .* α)
    @eval @inline Base.:/(A::$ST, α::Number) = _rebuild(A, getdata(A) ./ α)
end

# ── Reference checks for binary operations ───────────────────────────────────
# TensWalpole and TensTI carry a symmetry axis (n); TensOrtho carries a material
# frame.  Binary arithmetic requires the same reference.

@inline _check_same_reference(A::TensWalpole, B::TensWalpole) =
    @assert getaxis(A) == getaxis(B) "TensWalpole operation requires the same axis"
@inline _check_same_reference(A::TensTI, B::TensTI) =
    @assert getaxis(A) == getaxis(B) "TensTI operation requires the same axis"
@inline _check_same_reference(A::TensOrtho, B::TensOrtho) =
    @assert getframe(A) == getframe(B) "TensOrtho operation requires the same material frame"
