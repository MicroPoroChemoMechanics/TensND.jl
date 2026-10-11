# ── Kelvin-Mandel notation ────────────────────────────────────────────────────

"""
    KM(t::AbstractTens{order,dim}; kwargs...)
    KM(t::AbstractTens{order,dim}, b::AbstractBasis{dim}; kwargs...)
    KM(t::AbstractTens{order,dim}, b::AbstractBasis{dim}, var::NTuple{order,Symbol}; kwargs...)

Write the components of a second or fourth order tensor in Kelvin-Mandel notation,
in its own basis or in the basis `b` (with the variances `var`, in the order of
[`components`](@ref)). The keyword arguments go to `Tensors.tomandel`.

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
    var::NTuple{order, Symbol};
    kwargs...,
) where {order, dim} = tomandel(tensor_or_array(components(t, b, var)); kwargs...)

KM(t::AbstractTens{order, dim}, b::AbstractBasis{dim}; kwargs...) where {order, dim} =
    tomandel(tensor_or_array(components(t, b)); kwargs...)


KM(t::AbstractArray; kwargs...) = KM(Tens(t); kwargs...)
KM(t::AbstractArray, b::AbstractBasis; kwargs...) = KM(Tens(t), b; kwargs...)


"""
    inv_KM(v::AbstractVecOrMat; kwargs...)

Define a tensor from a Kelvin-Mandel vector or matrix representation
"""
inv_KM(TT::Type{<:Tensors.AllTensors}, v::AbstractVecOrMat; kwargs...) =
    Tens(frommandel(TT, v; kwargs...))
inv_KM(TT::Type{<:Tensors.AllTensors}, v::AbstractVecOrMat, b::AbstractBasis; kwargs...) =
    Tens(frommandel(TT, v; kwargs...), b)
# The shape of `v` gives the tensor type: 6×6 → SymmetricTensor{4,3},
# 9×9 → Tensor{4,3}, 3×3 → SymmetricTensor{4,2}, 4×4 → Tensor{4,2}, and the
# vectors of the same lengths → the order-2 types. Each branch names its type
# literally at its own call site, so that each is inferred concretely and the
# return type is a union of four concrete types, which Julia splits; a type
# returned by a helper, or read from a table, would be a run-time value and
# `frommandel` would be reached by dynamic dispatch.

function inv_KM(v::AbstractMatrix; kwargs...)
    m, n = size(v)
    (m, n) == (6, 6) && return inv_KM(SymmetricTensor{4, 3}, v; kwargs...)
    (m, n) == (9, 9) && return inv_KM(Tensor{4, 3}, v; kwargs...)
    (m, n) == (3, 3) && return inv_KM(SymmetricTensor{4, 2}, v; kwargs...)
    (m, n) == (4, 4) && return inv_KM(Tensor{4, 2}, v; kwargs...)
    throw(DimensionMismatch("no Kelvin-Mandel tensor type for a $m×$n matrix"))
end

function inv_KM(v::AbstractVector; kwargs...)
    n = length(v)
    n == 6 && return inv_KM(SymmetricTensor{2, 3}, v; kwargs...)
    n == 9 && return inv_KM(Tensor{2, 3}, v; kwargs...)
    n == 3 && return inv_KM(SymmetricTensor{2, 2}, v; kwargs...)
    n == 4 && return inv_KM(Tensor{2, 2}, v; kwargs...)
    throw(DimensionMismatch("no Kelvin-Mandel tensor type for a length-$n vector"))
end

function inv_KM(v::AbstractMatrix, b::AbstractBasis; kwargs...)
    m, n = size(v)
    (m, n) == (6, 6) && return inv_KM(SymmetricTensor{4, 3}, v, b; kwargs...)
    (m, n) == (9, 9) && return inv_KM(Tensor{4, 3}, v, b; kwargs...)
    (m, n) == (3, 3) && return inv_KM(SymmetricTensor{4, 2}, v, b; kwargs...)
    (m, n) == (4, 4) && return inv_KM(Tensor{4, 2}, v, b; kwargs...)
    throw(DimensionMismatch("no Kelvin-Mandel tensor type for a $m×$n matrix"))
end

function inv_KM(v::AbstractVector, b::AbstractBasis; kwargs...)
    n = length(v)
    n == 6 && return inv_KM(SymmetricTensor{2, 3}, v, b; kwargs...)
    n == 9 && return inv_KM(Tensor{2, 3}, v, b; kwargs...)
    n == 3 && return inv_KM(SymmetricTensor{2, 2}, v, b; kwargs...)
    n == 4 && return inv_KM(Tensor{2, 2}, v, b; kwargs...)
    throw(DimensionMismatch("no Kelvin-Mandel tensor type for a length-$n vector"))
end


export KM, inv_KM
