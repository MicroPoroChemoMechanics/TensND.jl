# ── Kelvin-Mandel notation ────────────────────────────────────────────────────

"""
    mandel(t::AbstractTens{order,dim}; kwargs...)
    mandel(t::AbstractTens{order,dim}, b::AbstractBasis{dim}; kwargs...)
    mandel(t::AbstractTens{order,dim}, b::AbstractBasis{dim}, var::NTuple{order,Symbol}; kwargs...)

Write the components of a second or fourth order tensor in Kelvin-Mandel notation,
in its own basis or in the basis `b` (with the variances `var`, in the order of
[`components`](@ref)). The keyword arguments go to `Tensors.tomandel`.

# Examples
```julia
julia> σ = Tens(SymmetricTensor{2,3}((i, j) -> symbols("σ\$i\$j", real = true))) ;

julia> mandel(σ)
6-element Vector{Sym}:
         σ11
         σ22
         σ33
      √2⋅σ32
      √2⋅σ31
      √2⋅σ21

julia> C = Tens(SymmetricTensor{4,3}((i, j, k, l) -> symbols("C\$i\$j\$k\$l", real = true))) ;

julia> mandel(C)
6×6 Matrix{Sym}:
         C₁₁₁₁     C₁₁₂₂     C₁₁₃₃  √2⋅C₁₁₃₂  √2⋅C₁₁₃₁  √2⋅C₁₁₂₁
         C₂₂₁₁     C₂₂₂₂     C₂₂₃₃  √2⋅C₂₂₃₂  √2⋅C₂₂₃₁  √2⋅C₂₂₂₁
         C₃₃₁₁     C₃₃₂₂     C₃₃₃₃  √2⋅C₃₃₃₂  √2⋅C₃₃₃₁  √2⋅C₃₃₂₁
      √2⋅C₃₂₁₁  √2⋅C₃₂₂₂  √2⋅C₃₂₃₃   2⋅C₃₂₃₂   2⋅C₃₂₃₁   2⋅C₃₂₂₁
      √2⋅C₃₁₁₁  √2⋅C₃₁₂₂  √2⋅C₃₁₃₃   2⋅C₃₁₃₂   2⋅C₃₁₃₁   2⋅C₃₁₂₁
      √2⋅C₂₁₁₁  √2⋅C₂₁₂₂  √2⋅C₂₁₃₃   2⋅C₂₁₃₂   2⋅C₂₁₃₁   2⋅C₂₁₂₁
```
"""
mandel(t::Tensors.AllTensors; kwargs...) = tomandel(t; kwargs...)
mandel(t::AbstractTens; kwargs...) = tomandel(get_array(t); kwargs...)

mandel(
    t::AbstractTens{order, dim},
    b::AbstractBasis{dim},
    var::NTuple{order, Symbol};
    kwargs...,
) where {order, dim} = tomandel(_static_tensor(components(t, b, var), Val(dim), Val(order)); kwargs...)

mandel(t::AbstractTens{order, dim}, b::AbstractBasis{dim}; kwargs...) where {order, dim} =
    tomandel(_static_tensor(components(t, b), Val(dim), Val(order)); kwargs...)


mandel(t::AbstractArray; kwargs...) = mandel(Tens(t); kwargs...)
mandel(t::AbstractArray, b::AbstractBasis; kwargs...) = mandel(Tens(t), b; kwargs...)


"""
    from_mandel(v::AbstractVecOrMat[, b::AbstractBasis]; kwargs...)
    from_mandel(TT::Type{<:Tensors.AllTensors}, v[, b]; kwargs...)

The tensor whose Kelvin-Mandel vector or matrix is `v`, in the basis `b` (the
canonical one by default). The size of `v` gives the Tensors.jl type `TT` —
inferred when `v` is a static array — unless it is given.
"""
from_mandel(TT::Type{<:Tensors.AllTensors}, v::AbstractVecOrMat; kwargs...) =
    Tens(frommandel(TT, v; kwargs...))
from_mandel(TT::Type{<:Tensors.AllTensors}, v::AbstractVecOrMat, b::AbstractBasis; kwargs...) =
    Tens(frommandel(TT, v; kwargs...), b)
# The shape of `v` gives the tensor type: 6×6 → SymmetricTensor{4,3}, 9×9 →
# Tensor{4,3}, 3×3 → SymmetricTensor{4,2}, 4×4 → Tensor{4,2}, and the vectors of
# the same lengths → the order-2 types.
#
# A static array carries its size in its type, and the tensor type is chosen by
# dispatch: inferred, and the fastest route. A plain array knows its size only
# at run time, so the type is chosen by tests on it: a dispatch on
# `Val(size(v))` would be a dynamic one, measured twelve times slower. Each
# branch names its type literally, so the result is a union of four concrete
# types, which Julia splits.
from_mandel(v::StaticMatrix{6, 6}, b::AbstractBasis...; kwargs...) = from_mandel(SymmetricTensor{4, 3}, v, b...; kwargs...)
from_mandel(v::StaticMatrix{9, 9}, b::AbstractBasis...; kwargs...) = from_mandel(Tensor{4, 3}, v, b...; kwargs...)
from_mandel(v::StaticMatrix{3, 3}, b::AbstractBasis...; kwargs...) = from_mandel(SymmetricTensor{4, 2}, v, b...; kwargs...)
from_mandel(v::StaticMatrix{4, 4}, b::AbstractBasis...; kwargs...) = from_mandel(Tensor{4, 2}, v, b...; kwargs...)
from_mandel(v::StaticVector{6}, b::AbstractBasis...; kwargs...) = from_mandel(SymmetricTensor{2, 3}, v, b...; kwargs...)
from_mandel(v::StaticVector{9}, b::AbstractBasis...; kwargs...) = from_mandel(Tensor{2, 3}, v, b...; kwargs...)
from_mandel(v::StaticVector{3}, b::AbstractBasis...; kwargs...) = from_mandel(SymmetricTensor{2, 2}, v, b...; kwargs...)
from_mandel(v::StaticVector{4}, b::AbstractBasis...; kwargs...) = from_mandel(Tensor{2, 2}, v, b...; kwargs...)

function from_mandel(v::AbstractMatrix, b::AbstractBasis...; kwargs...)
    m, n = size(v)
    (m, n) == (6, 6) && return from_mandel(SymmetricTensor{4, 3}, v, b...; kwargs...)
    (m, n) == (9, 9) && return from_mandel(Tensor{4, 3}, v, b...; kwargs...)
    (m, n) == (3, 3) && return from_mandel(SymmetricTensor{4, 2}, v, b...; kwargs...)
    (m, n) == (4, 4) && return from_mandel(Tensor{4, 2}, v, b...; kwargs...)
    throw(DimensionMismatch("no Kelvin-Mandel tensor type for a $m×$n matrix"))
end

function from_mandel(v::AbstractVector, b::AbstractBasis...; kwargs...)
    n = length(v)
    n == 6 && return from_mandel(SymmetricTensor{2, 3}, v, b...; kwargs...)
    n == 9 && return from_mandel(Tensor{2, 3}, v, b...; kwargs...)
    n == 3 && return from_mandel(SymmetricTensor{2, 2}, v, b...; kwargs...)
    n == 4 && return from_mandel(Tensor{2, 2}, v, b...; kwargs...)
    throw(DimensionMismatch("no Kelvin-Mandel tensor type for a length-$n vector"))
end


export mandel, from_mandel
