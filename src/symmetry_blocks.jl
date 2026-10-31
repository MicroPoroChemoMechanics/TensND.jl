# =============================================================================
#  symmetry_blocks.jl — each symmetry class as a sum of small matrix algebras.
#
#  The tensors of one symmetry class, in one frame, are closed under ⊡ and form
#  an algebra isomorphic to a direct sum of full matrix algebras: the double
#  contraction of two of them is the product of their blocks, block by block,
#  and the inverse is the inverse of each block.
#
#    class                     order 4                     order 2
#    TensISO                   (α), (β)                    (λ)
#    TensCubic                 (α), (β), (γ)               —
#    TensTI                    [ℓ₁ ℓ₃; ℓ₄ ℓ₂], (ℓ₅), (ℓ₆)   (a), (b)
#    TensOrtho                 3×3 normal block,           (a₁), (a₂), (a₃)
#                              (2C₄₄), (2C₅₅), (2C₆₆)
#
#  Whatever operates on such tensors through their products and inverses —
#  a Volterra operator whose coefficients are functions of two times, say —
#  works on the blocks, the same way for every class.
# =============================================================================


_block1(x) = SMatrix{1, 1}(x)

"""
    symmetry_blocks(A) -> Tuple of SMatrix

The blocks of `A` in the matrix algebra of its symmetry class: the small matrices whose
products and inverses, block by block, are those of the tensors,

```math
\\mathrm{blocks}(\\mathbb{A}:\\mathbb{B}) = \\mathrm{blocks}(\\mathbb{A})\\,\\mathrm{blocks}(\\mathbb{B}),
\\qquad
\\mathrm{blocks}(\\mathbb{A}^{-1}) = \\mathrm{blocks}(\\mathbb{A})^{-1}.
```

| class | order 4 | order 2 |
|:--|:--|:--|
| [`TensISO`](@ref) | ``(\\alpha)``, ``(\\beta)`` | ``(\\lambda)`` |
| [`TensCubic`](@ref) | ``(\\alpha)``, ``(\\beta)``, ``(\\gamma)`` | — |
| [`TensTI`](@ref) | ``\\begin{pmatrix}\\ell_1&\\ell_3\\\\\\ell_4&\\ell_2\\end{pmatrix}``, ``(\\ell_5)``, ``(\\ell_6)`` | ``(a)``, ``(b)`` |
| [`TensOrtho`](@ref) | the normal ``3\\times 3`` block, ``(2C_{44})``, ``(2C_{55})``, ``(2C_{66})`` | ``(a_1)``, ``(a_2)``, ``(a_3)`` |

The frame (axis, material frame) is not part of the blocks: [`from_symmetry_blocks`](@ref)
takes it from a tensor of the class.

```jldoctest
julia> A = TensTI{4}(4.0, 3.0, 1.0, 1.0, 2.0, 5.0, (0.0, 0.0, 1.0));

julia> symmetry_blocks(A)[1]
2×2 StaticArraysCore.SMatrix{2, 2, Float64, 4} with indices SOneTo(2)×SOneTo(2):
 4.0  1.0
 1.0  3.0
```
"""
symmetry_blocks(A::TensISO{2}) = (_block1(A.data[1]),)
symmetry_blocks(A::TensISO{4}) = (_block1(A.data[1]), _block1(A.data[2]))
symmetry_blocks(A::TensCubic) = map(_block1, A.data)
symmetry_blocks(A::TensTI{2, <:Any, 2}) = map(_block1, A.data)
function symmetry_blocks(A::Union{TensTI{4, <:Any, 5}, TensTI{4, <:Any, 6}})
    ℓ₁, ℓ₂, ℓ₃, ℓ₄, ℓ₅, ℓ₆ = walpole_coeffs(A)
    return (SMatrix{2, 2}(ℓ₁, ℓ₄, ℓ₃, ℓ₂), _block1(ℓ₅), _block1(ℓ₆))
end
symmetry_blocks(A::TensOrtho{2}) = map(_block1, A.data)
function symmetry_blocks(A::TensOrtho{4, <:Any, 9})
    C₁₁, C₂₂, C₃₃, C₁₂, C₁₃, C₂₃, C₄₄, C₅₅, C₆₆ = A.data
    return (SMatrix{3, 3}(C₁₁, C₁₂, C₁₃, C₁₂, C₂₂, C₂₃, C₁₃, C₂₃, C₃₃), _block1(2C₄₄), _block1(2C₅₅), _block1(2C₆₆))
end
function symmetry_blocks(A::TensOrtho{4, <:Any, 12})
    C₁₁, C₂₂, C₃₃, C₁₂, C₁₃, C₂₃, C₄₄, C₅₅, C₆₆, C₂₁, C₃₁, C₃₂ = A.data
    return (SMatrix{3, 3}(C₁₁, C₂₁, C₃₁, C₁₂, C₂₂, C₃₂, C₁₃, C₂₃, C₃₃), _block1(2C₄₄), _block1(2C₅₅), _block1(2C₆₆))
end

"""
    from_symmetry_blocks(A, blocks) -> AbstractTens

The tensor of the class and frame of `A` whose blocks are `blocks`: the inverse of
[`symmetry_blocks`](@ref). A product of two symmetric tensors of a class need not be
symmetric, so a transversely isotropic or orthotropic result takes the general shape of its
class (six Walpole coefficients, twelve orthotropic constants).
"""
from_symmetry_blocks(::TensISO{2, dim}, b) where {dim} = TensISO{dim}(b[1][1])
from_symmetry_blocks(::TensISO{4, dim}, b) where {dim} = TensISO{dim}(b[1][1], b[2][1])
from_symmetry_blocks(A::TensCubic, b) = TensCubic(b[1][1], b[2][1], b[3][1], A.frame)
from_symmetry_blocks(A::TensTI{2, <:Any, 2}, b) = TensTI{2}(b[1][1], b[2][1], A.n)
function from_symmetry_blocks(A::Union{TensTI{4, <:Any, 5}, TensTI{4, <:Any, 6}}, b)
    W = b[1]
    return TensTI{4}(W[1, 1], W[2, 2], W[1, 2], W[2, 1], b[2][1], b[3][1], A.n)
end
from_symmetry_blocks(A::TensOrtho{2}, b) =
    TensOrtho{2, _blocks_eltype(b)}((b[1][1], b[2][1], b[3][1]), A.frame)
function from_symmetry_blocks(A::TensOrtho{4}, b)
    N, s₄, s₅, s₆ = b[1], b[2][1], b[3][1], b[4][1]
    data = (N[1, 1], N[2, 2], N[3, 3], N[1, 2], N[1, 3], N[2, 3], s₄ / 2, s₅ / 2, s₆ / 2, N[2, 1], N[3, 1], N[3, 2])
    return TensOrtho{4, _blocks_eltype(b)}(data, A.frame)
end

_blocks_eltype(b) = mapreduce(eltype, promote_type, b)

export symmetry_blocks, from_symmetry_blocks

# ── The blocks of a matrix of components ───────────────────────────────────────
#
# In an orthonormal basis of the symmetric tensors adapted to the class and its
# frame (`symmetry_basis`), the Kelvin-Mandel matrix of a member of the class
# (order 4), or its matrix of components (order 2), is block diagonal: each
# block of `symmetry_blocks` sits on the diagonal once per copy listed in
# `symmetry_layout`. Read on any matrix, the first copy of each block is the
# projection a computation on matrices of components needs; the other entries
# say whether the matrix belongs to the class.

# The material-frame bases, as columns of Kelvin-Mandel vectors: the isotropic
# one (spherical, then the two normal deviators, then the three shears), and the
# transversely isotropic one about e₃ (axial, transverse spherical, then the
# two in-plane deviators, then the two out-of-plane shears).
function _iso_mandel_basis(::Type{T}) where {T}
    o, z = one(T), zero(T)
    a, b, c = inv(sqrt(3o)), inv(sqrt(2o)), inv(sqrt(6o))
    return SMatrix{6, 6, T}(
        a, a, a, z, z, z,
        b, -b, z, z, z, z,
        c, c, -2c, z, z, z,
        z, z, z, o, z, z,
        z, z, z, z, o, z,
        z, z, z, z, z, o,
    )
end
function _ti_mandel_basis(::Type{T}) where {T}
    o, z = one(T), zero(T)
    b = inv(sqrt(2o))
    return SMatrix{6, 6, T}(
        z, z, o, z, z, z,
        b, b, z, z, z, z,
        b, -b, z, z, z, z,
        z, z, z, z, z, o,
        z, z, z, o, z, z,
        z, z, z, z, o, z,
    )
end

_basis_type(A) = float(eltype(A))

"""
    symmetry_basis(A) -> SMatrix

The orthonormal basis, as columns, in which the matrix of a member of the class
and frame of `A` is block diagonal: Kelvin-Mandel vectors (``6\\times6``) at
order 4, vectors (``3\\times3``) at order 2. Each block of
[`symmetry_blocks`](@ref) sits on the diagonal of ``Q^{\\mathsf T}MQ`` at the
positions [`symmetry_layout`](@ref) lists, once per copy.

| class | columns (in the material frame) |
|:--|:--|
| ISO, CUBIC | the spherical direction, the two normal deviators, the three shears |
| TI about ``\\underline{n}`` | axial, transverse spherical, the two in-plane deviators, the two out-of-plane shears |
| ORTHO | the frame's normal and shear directions |
"""
symmetry_basis(A::TensISO{4, 3}) = _iso_mandel_basis(_basis_type(A))
symmetry_basis(A::TensCubic) = _km_congruence(_frame_matrix(A.frame)) * _iso_mandel_basis(_basis_type(A))
symmetry_basis(A::Union{TensTI{4, <:Any, 5}, TensTI{4, <:Any, 6}}) =
    _km_congruence(first(_axis_frame(A.n))) * _ti_mandel_basis(_basis_type(A))
symmetry_basis(A::TensOrtho{4}) = _km_congruence(_frame_matrix(A.frame))
symmetry_basis(A::TensISO{2, 3}) = SMatrix{3, 3, _basis_type(A)}(I)
symmetry_basis(A::TensTI{2, <:Any, 2}) = first(_axis_frame(A.n))
symmetry_basis(A::TensOrtho{2}) = _frame_matrix(A.frame)

"""
    symmetry_layout(A) -> Tuple

For each block of [`symmetry_blocks`](@ref)`(A)`, the index sets of its copies on
the diagonal of ``Q^{\\mathsf T}MQ``, ``Q`` = [`symmetry_basis`](@ref)`(A)`: the
isotropic shear modulus ``\\beta`` sits on five, the transverse shear ``\\ell_5``
on two, and so on.
"""
symmetry_layout(::TensISO{4, 3}) = ((SVector(1),), (SVector(2), SVector(3), SVector(4), SVector(5), SVector(6)))
symmetry_layout(::TensCubic) = ((SVector(1),), (SVector(2), SVector(3)), (SVector(4), SVector(5), SVector(6)))
symmetry_layout(::Union{TensTI{4, <:Any, 5}, TensTI{4, <:Any, 6}}) =
    ((SVector(1, 2),), (SVector(3), SVector(4)), (SVector(5), SVector(6)))
symmetry_layout(::TensOrtho{4}) = ((SVector(1, 2, 3),), (SVector(4),), (SVector(5),), (SVector(6),))
symmetry_layout(::TensISO{2, 3}) = ((SVector(1), SVector(2), SVector(3)),)
symmetry_layout(::TensTI{2, <:Any, 2}) = ((SVector(1), SVector(2)), (SVector(3),))
symmetry_layout(::TensOrtho{2}) = ((SVector(1),), (SVector(2),), (SVector(3),))

_matrix_size(::AbstractTens{4}) = 6
_matrix_size(::AbstractTens{2}) = 3

"""
    symmetry_blocks(A, M::AbstractMatrix) -> Tuple of SMatrix

The blocks, in the class and frame of `A`, of the matrix `M`: a ``6\\times6``
Kelvin-Mandel matrix at order 4, a ``3\\times3`` matrix of components at order 2,
in the canonical basis. Each block is read on its first copy in
``Q^{\\mathsf T}MQ``; when `M` belongs to the class,
`symmetry_matrix(A, symmetry_blocks(A, M)) ≈ M`, and the blocks are those of the
tensor `M` stands for.
"""
function symmetry_blocks(A::AbstractTens, M::AbstractMatrix)
    N = _matrix_size(A)
    size(M) == (N, N) || throw(DimensionMismatch("symmetry_blocks: expected a $(N)×$(N) matrix, got $(size(M))"))
    Q, S = symmetry_basis(A), SMatrix{N, N}(M)
    return map(copies -> (Qc = Q[:, first(copies)]; transpose(Qc) * S * Qc), symmetry_layout(A))
end

"""
    symmetry_matrix(A, blocks) -> SMatrix

The matrix (Kelvin-Mandel at order 4, of components at order 2, in the canonical
basis) of the member of the class and frame of `A` whose blocks are `blocks`:
the inverse of [`symmetry_blocks`](@ref)`(A, M)` on the class.
"""
function symmetry_matrix(A::AbstractTens, b)
    Q = symmetry_basis(A)
    # Σₖ Σ_copies Q_c bₖ Q_cᵀ, Q_c the columns of Q a copy of block k occupies.
    terms = map(symmetry_layout(A), b) do copies, bₖ
        sum(c -> (Qc = Q[:, c]; Qc * bₖ * transpose(Qc)), copies)
    end
    return reduce(+, terms)
end

export symmetry_basis, symmetry_layout, symmetry_matrix
