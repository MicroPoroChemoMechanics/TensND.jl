# =============================================================================
#  symmetric_parts.jl — `Tensors.majorsymmetric` (order 4) and
#  `Tensors.symmetric` (order 2) on TensND tensors, the symmetry class kept
#  where the class allows it.
# =============================================================================

Tensors.majorsymmetric(t::Union{TensISO{4}, TensCubic, TensOrtho}) = t
Tensors.majorsymmetric(t::TensTI{4, T, 5}) where {T} = t
# In the Walpole basis the major transpose swaps ℓ₃ and ℓ₄ and changes the sign
# of ℓ₇ and ℓ₈, the major-antisymmetric generators of the six- and
# eight-parameter forms: the major-symmetric part is the five-parameter tensor.
function Tensors.majorsymmetric(t::TensTI{4, T}) where {T}
    ℓ₁, ℓ₂, ℓ₃, ℓ₄, ℓ₅, ℓ₆, _, _ = get_ℓ8(t)
    return TensTI{4, T, 5}((ℓ₁, ℓ₂, (ℓ₃ + ℓ₄) / 2, ℓ₅, ℓ₆), axis(t))
end

"""
    Tensors.majorsymmetric(t::AbstractTens{4}) -> AbstractTens{4}

Major-symmetric part ``\\tfrac12(t_{ijkl} + t_{klij})`` of a fourth-order
tensor, a method of the `Tensors.jl` function.

The symmetry class is kept: an isotropic, cubic or orthotropic tensor is its
own major-symmetric part, and a transversely isotropic tensor in the six- or
eight-parameter Walpole form (a product of TI tensors lands there) returns the
five-parameter form about the same axis. A general tensor stays in its basis;
in a basis that is not orthonormal, the transpose pairs components of the same
variance, so they are first written with the variance of the first pair on the
second as well.

See also `Tensors.symmetric` for a second-order tensor.
"""
function Tensors.majorsymmetric(t::AbstractTens{4})
    v = get_var(t)
    w = (v[1], v[2], v[1], v[2])
    A = components(t, w)
    return Tens((A + permutedims(A, (3, 4, 1, 2))) / 2, get_basis(t), w)
end

Tensors.symmetric(t::TensISO{2}) = t
Tensors.symmetric(t::TensTI{2, T, 2}) where {T} = t
# a·𝟏ₜ + b·n⊗n + c·w: the skew in-plane part c·w goes.
function Tensors.symmetric(t::TensTI{2, T, 3}) where {T}
    a, b, _ = get_data(t)
    return TensTI{2}(a, b, axis(t))
end

"""
    Tensors.symmetric(t::AbstractTens{2}) -> AbstractTens{2}

Symmetric part ``\\tfrac12(t_{ij} + t_{ji})`` of a second-order tensor, a
method of the `Tensors.jl` function. An isotropic or transversely isotropic
tensor keeps its class (the skew part `c·w` of a three-parameter `TensTI{2}`
goes); a general one stays in its basis, its components written with one
variance before transposing.
"""
function Tensors.symmetric(t::AbstractTens{2})
    v = get_var(t)
    w = (v[1], v[1])
    A = components(t, w)
    return Tens((A + transpose(A)) / 2, get_basis(t), w)
end
