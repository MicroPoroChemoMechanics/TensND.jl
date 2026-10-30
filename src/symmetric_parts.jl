# =============================================================================
#  symmetric_parts.jl — `Tensors.majorsymmetric` (order 4) and
#  `Tensors.symmetric` (order 2) on TensND tensors, the symmetry class kept
#  where the class allows it.
# =============================================================================

Tensors.majorsymmetric(t::Union{TensISO{4}, TensCubic, TensOrtho{4, <:Any, 9}}) = t
# The symmetric part of the normal block, the nine-constant form on the same frame.
function Tensors.majorsymmetric(t::TensOrtho{4, <:Any, 12})
    C11, C22, C33, C12, C13, C23, C44, C55, C66, C21, C31, C32 = get_data(t)
    return _rebuild(t, (C11, C22, C33, (C12 + C21) / 2, (C13 + C31) / 2, (C23 + C32) / 2, C44, C55, C66))
end
Tensors.majorsymmetric(t::TensTI{4, T, 5}) where {T} = t
# In the Walpole basis the major transpose swaps ℓ₃ and ℓ₄ and changes the sign
# of ℓ₇ and ℓ₈, the major-antisymmetric generators of the six- and
# eight-parameter forms: the major-symmetric part is the five-parameter tensor.
function Tensors.majorsymmetric(t::TensTI{4, T}) where {T}
    ℓ₁, ℓ₂, ℓ₃, ℓ₄, ℓ₅, ℓ₆, _, _ = walpole_coeffs(t, Val(8))
    return TensTI{4, T, 5}((ℓ₁, ℓ₂, (ℓ₃ + ℓ₄) / 2, ℓ₅, ℓ₆), axis(t))
end

"""
    Tensors.majorsymmetric(t::AbstractTens{4}) -> AbstractTens{4}

Major-symmetric part ``\\tfrac12(t_{ijkl} + t_{klij})`` of a fourth-order
tensor, a method of the `Tensors.jl` function.

The symmetry class is kept: an isotropic, cubic or nine-constant orthotropic
tensor is its own major-symmetric part, a twelve-constant orthotropic one (a
product of orthotropic tensors lands there) returns the nine-constant form on
the same frame, and a transversely isotropic tensor in the six- or
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

# Every structured fourth-order tensor has the minor symmetries.
Tensors.minorsymmetric(t::Union{TensISO{4}, TensTI{4}, TensOrtho{4}, TensCubic}) = t

"""
    Tensors.minorsymmetric(t::AbstractTens{4}) -> AbstractTens{4}

Minor-symmetric part ``\\tfrac14(t_{ijkl} + t_{jikl} + t_{ijlk} + t_{jilk})`` of a
fourth-order tensor, a method of the `Tensors.jl` function. A structured tensor
(isotropic, transversely isotropic, orthotropic, cubic) is its own minor-symmetric
part; a general tensor stays in its basis, stored symmetric. In a basis that is
not orthonormal, the transposed indices are first given the same variance, that
of the first index of each pair.

For instance, the half sum ``\\tfrac12\\bigl((\\underline{a}\\otimes\\underline{a})
\\stackrel{s}{\\boxtimes}\\boldsymbol{N} + \\boldsymbol{N}\\stackrel{s}{\\boxtimes}
(\\underline{a}\\otimes\\underline{a})\\bigr)``, that is
``\\underline{a}\\stackrel{s}{\\otimes}\\boldsymbol{N}\\stackrel{s}{\\otimes}\\underline{a}``
for a symmetric ``\\boldsymbol{N}``, is `minorsymmetric((a ⊗ a) ⊠ˢ N)`: one product
instead of two.

See also `Tensors.majorsymmetric`.
"""
function Tensors.minorsymmetric(t::AbstractTens{4})
    v = get_var(t)
    w = (v[1], v[1], v[3], v[3])
    return Tens(_minorsymmetric(components(t, w)), get_basis(t), w)
end
_minorsymmetric(A::Tensors.SymmetricTensor{4}) = A
_minorsymmetric(A::Tensors.Tensor{4}) = Tensors.minorsymmetric(A)
_minorsymmetric(A::AbstractArray{T, 4}) where {T} =
    ((A + permutedims(A, (2, 1, 3, 4))) + (permutedims(A, (1, 2, 4, 3)) + permutedims(A, (2, 1, 4, 3)))) / 4

Tensors.symmetric(t::TensISO{2}) = t
Tensors.symmetric(t::TensOrtho{2}) = t
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
