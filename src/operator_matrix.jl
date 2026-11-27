# =============================================================================
#  operator_matrix.jl — a tensor of order 2 or 4 as the matrix of the linear
#  map it is, and what is read on that matrix: definiteness and the departure
#  from major symmetry.
#
#  A tensor of order 4 with the minor symmetries maps symmetric tensors of
#  order 2 onto symmetric tensors of order 2 by double contraction; written on
#  the Kelvin-Mandel components (`kelvin_mandel.jl`), whose basis is
#  orthonormal, it is an n×n matrix, n = d(d+1)/2, and the algebra of the
#  tensors is that of the matrices: `⊡` is the product, the major transpose the
#  transpose, the Frobenius norm that of the tensor. A tensor of order 2 maps
#  vectors onto vectors: its d×d matrix of components.
#
#  The components are read in the canonical basis, or in an orthonormal basis
#  given (`components`), never in the basis a tensor happens to be stored in.
#  The structured classes are written from their coefficients, in closed form
#  (`_km_canonical`, `structured_tens_promotion.jl`). Dimensions 2 and 3 give
#  static matrices, the others plain ones.
# =============================================================================

# ── The matrix of a tensor ───────────────────────────────────────────────────

"""
    operator_matrix(t::AbstractTens{4}; symmetrize = true)
    operator_matrix(t::AbstractTens{4}, ℬ::OrthonormalBasis; symmetrize = true)
    operator_matrix(t::AbstractTens{2})
    operator_matrix(t::AbstractTens{2}, ℬ::OrthonormalBasis)

The matrix of the linear map the tensor `t` is, in the canonical basis or in
the orthonormal basis `ℬ`.

At order 4, ``\\mathbb{C}`` maps a symmetric tensor of order 2 onto
``\\boldsymbol{\\sigma} = \\mathbb{C} : \\boldsymbol{\\varepsilon}``. With the
Kelvin-Mandel components ``\\hat{\\varepsilon}_I = w_I\\,\\varepsilon_{ij}``
([`mandel_pairs`](@ref), [`mandel_weights`](@ref): ``w_I = 1`` on the diagonal,
``\\sqrt{2}`` off it), the map is ``\\hat{\\sigma} = M\\,\\hat{\\varepsilon}``, with
the ``n \\times n`` matrix, ``n = d(d+1)/2``,

```math
M_{IJ} = w_I\\,w_J\\,\\mathbb{C}^{s}_{ijkl},
\\qquad
\\mathbb{C}^{s}_{ijkl} = \\tfrac14\\left(\\mathbb{C}_{ijkl} + \\mathbb{C}_{jikl} + \\mathbb{C}_{ijlk} + \\mathbb{C}_{jilk}\\right),
```

``(i, j)`` and ``(k, l)`` the pairs of the slots ``I`` and ``J``. The basis of
the Kelvin-Mandel components is orthonormal, so that the matrix of
``\\mathbb{A} : \\mathbb{B}`` is the product of the matrices, that of the major
transpose the transpose, and the Frobenius norm that of the tensor.

With `symmetrize = true` the matrix is that of the minor-symmetric part
``\\mathbb{C}^{s}``, whatever the storage of `t`: it is
[`mandel`](@ref)`(t, ℬ; symmetrize = true)`, and reads a minor-symmetric tensor
exactly. With `symmetrize = false` the components must have the minor
symmetries (to the tolerance of their storage), and a tensor that has not is
refused.

At order 2 the matrix is that of the components ``t_{ij}``, the map
``\\underline{v} \\mapsto \\boldsymbol{t} \\cdot \\underline{v}`` on vectors.

In dimensions 2 and 3 the result is an `SMatrix` (`3×3` or `6×6` at order 4,
`d×d` at order 2), in other dimensions a `Matrix`. The reciprocal is
[`from_operator_matrix`](@ref).

See also [`definiteness`](@ref), [`major_asymmetry`](@ref), [`mandel`](@ref),
[`canonical_components`](@ref).
"""
operator_matrix(t::AbstractTens{4, dim}; symmetrize::Bool = true) where {dim} =
    _operator_matrix(t, CanonicalBasis{dim, eltype(t)}(), symmetrize)
operator_matrix(t::AbstractTens{4, dim}, ℬ::OrthonormalBasis{dim}; symmetrize::Bool = true) where {dim} =
    _operator_matrix(t, ℬ, symmetrize)

operator_matrix(t::AbstractTens{2, dim}) where {dim} = operator_matrix(t, CanonicalBasis{dim, eltype(t)}())
operator_matrix(t::AbstractTens{2, dim}, ℬ::OrthonormalBasis{dim}) where {dim} =
    _static_components(Val(dim), components(t, ℬ))

# Any tensor: its components in ℬ, written as the Kelvin-Mandel matrix of their
# minor-symmetric part, which reads minor-symmetric components to the last bit.
# Without symmetrization the components are checked rather than averaged.
function _operator_matrix(t::AbstractTens{4, dim}, ℬ, symmetrize::Bool) where {dim}
    c = components(t, ℬ)
    symmetrize || _minor_symmetric(c) || _throw_not_minor_symmetric()
    return _km_of_components(c, Val(dim))
end
_km_of_components(c::AbstractArray{<:Any, 4}, d::Val) = _static_km(d, _mandel_minor(c, mandel_pairs(d)))

# An isotropic tensor has the same matrix in every orthonormal basis:
# `α 𝕁 + β 𝕂 = β 𝕀 + (α − β) 𝕁`, the matrix of 𝕁 being `𝐢 𝐢ᵀ / d`, `𝐢` the
# Kelvin-Mandel vector of the identity (one on the d diagonal slots).
function _operator_matrix(t::TensISO{4, dim, T}, ::OrthonormalBasis, ::Bool) where {dim, T}
    α, β = get_data(t)
    c = (α - β) / dim
    n = mandel_length(dim)
    M = [(I == J ? β : zero(T)) + (I ≤ dim && J ≤ dim ? c : zero(T)) for I in 1:n, J in 1:n]
    return _static_km(Val(dim), M)
end

# The classes that carry a frame or an axis are minor-symmetric by construction,
# and their canonical matrix has a closed form; in another basis they take the
# general reading.
const _KMClosedForm = Union{TensOrtho{4, 2}, TensOrtho{4, 3}, TensCubic{2}, TensCubic{3}, TensTI{4}}
_operator_matrix(t::_KMClosedForm, ℬ, symmetrize::Bool) = _closed_form_km(t, ℬ, symmetrize)
_closed_form_km(t, ::CanonicalBasis, ::Bool) = _km_canonical(t)
_closed_form_km(t::AbstractTens{4, dim}, ℬ::RotatedBasis, ::Bool) where {dim} = _km_of_components(components(t, ℬ), Val(dim))

@noinline _throw_not_minor_symmetric() = throw(
    ArgumentError(
        "operator_matrix: the components do not have the minor symmetries; " *
            "`symmetrize = true` takes the matrix of their minor-symmetric part"
    )
)

# Static sizes where Tensors.jl holds the components, plain matrices beyond.
_static_km(::Val{2}, M::AbstractMatrix) = SMatrix{3, 3, eltype(M)}(M)
_static_km(::Val{3}, M::AbstractMatrix) = SMatrix{6, 6, eltype(M)}(M)
_static_km(::Val, M::AbstractMatrix) = Matrix(M)

# Entry by entry: a `SymmetricTensor` stores only its independent components,
# which a linear index would read.
_static_components(::Val{2}, c::AbstractMatrix) = SMatrix{2, 2, eltype(c)}(c[i, j] for i in 1:2, j in 1:2)
_static_components(::Val{3}, c::AbstractMatrix) = SMatrix{3, 3, eltype(c)}(c[i, j] for i in 1:3, j in 1:3)
_static_components(::Val, c::AbstractMatrix) = Matrix(c)

"""
    from_operator_matrix(M::AbstractMatrix, Val(4)[, ℬ::OrthonormalBasis])
    from_operator_matrix(M::AbstractMatrix, Val(2)[, ℬ::OrthonormalBasis])

The tensor of order 4 or 2 whose matrix ([`operator_matrix`](@ref)) is `M`, in
the canonical basis or in the orthonormal basis `ℬ`: the reciprocal of
[`operator_matrix`](@ref) on the tensors with the minor symmetries.

At order 4, `M` is the ``n \\times n`` Kelvin-Mandel matrix of dimension ``d``,
``n = d(d+1)/2`` ([`mandel_dimension`](@ref)), and the tensor is
[`from_mandel`](@ref)`(M)`; a size that is that of no dimension is refused,
where `from_mandel` would read it as the form without symmetry. At order 2,
`M` is the ``d \\times d`` matrix of the components.
"""
from_operator_matrix(M::AbstractMatrix, ::Val{4}) = from_mandel(_checked_km(M))
from_operator_matrix(M::AbstractMatrix, ::Val{4}, ℬ::OrthonormalBasis) = from_mandel(_checked_km(M), ℬ)
from_operator_matrix(M::AbstractMatrix, ::Val{2}) = Tens(_checked_square(M))
from_operator_matrix(M::AbstractMatrix, ::Val{2}, ℬ::OrthonormalBasis) = Tens(_checked_square(M), ℬ)

_checked_square(M::AbstractMatrix) = (LinearAlgebra.checksquare(M); M)
_checked_km(M::AbstractMatrix) = (mandel_dimension(LinearAlgebra.checksquare(M)); M)

# ── What the matrix tells ────────────────────────────────────────────────────

"""
    definiteness(t::AbstractTens; rtol = 1.0e-12) → Symbol

The sign of the quadratic form of the tensor `t` of order 2 or 4,
``\\boldsymbol{\\varepsilon} : \\mathbb{C} : \\boldsymbol{\\varepsilon}`` on the
symmetric tensors of order 2 (``\\underline{v} \\cdot \\boldsymbol{t} \\cdot
\\underline{v}`` on the vectors at order 2). With the matrix ``M`` of `t`
([`operator_matrix`](@ref)), the form is that of its symmetric part,

```math
\\boldsymbol{\\varepsilon} : \\mathbb{C} : \\boldsymbol{\\varepsilon}
= \\hat{\\varepsilon} \\cdot M\\,\\hat{\\varepsilon}
= \\hat{\\varepsilon} \\cdot M^{s}\\,\\hat{\\varepsilon},
\\qquad M^{s} = \\tfrac12\\left(M + M^{\\mathsf{T}}\\right),
```

the matrix of the major-symmetric part of ``\\mathbb{C}`` at order 4, and its
sign is that of the eigenvalues ``\\lambda_a`` of ``M^{s}``. An eigenvalue with
``|\\lambda_a| \\le`` `rtol` ``\\max_b |\\lambda_b|`` is taken as zero, and the
result is

| result | eigenvalues |
|:--|:--|
| `:positive` | all ``> 0`` |
| `:negative` | all ``< 0`` |
| `:positive_semidefinite` | all ``\\ge 0``, one at least ``= 0`` |
| `:negative_semidefinite` | all ``\\le 0``, one at least ``= 0`` |
| `:indefinite` | some ``> 0``, some ``< 0`` |
| `:null` | all ``= 0`` (``M^{s} = 0``) |

The eigenvalues are those of the values: a dual number is read through
[`primal`](@ref). Comparing eigenvalues needs numbers: a symbolic tensor is
refused.

See also [`major_asymmetry`](@ref).
"""
function definiteness(t::Union{AbstractTens{2}, AbstractTens{4}}; rtol::Real = 1.0e-12)
    M = operator_matrix(t)
    return _definiteness(Val(is_hard_numeric(eltype(M))), M, rtol)
end

# A sign read on eigenvalues: numbers and duals only.
_definiteness(::Val{false}, M, rtol) = throw(
    ArgumentError(
        "definiteness compares eigenvalues, which symbolic entries of type $(eltype(M)) " *
            "do not allow; substitute numbers first"
    )
)
function _definiteness(::Val{true}, M, rtol)
    P = float.(primal.(M))
    λ = eigvals(Symmetric((P + transpose(P)) / 2))
    tol = rtol * maximum(abs, λ)
    positive, negative = any(>(tol), λ), any(<(-tol), λ)
    zero_found = any(x -> abs(x) ≤ tol, λ)
    positive && negative && return :indefinite
    positive && return zero_found ? :positive_semidefinite : :positive
    negative && return zero_found ? :negative_semidefinite : :negative
    return :null
end

"""
    major_asymmetry(t::AbstractTens{4}) → Real

The relative distance of the minor-symmetric part of `t` to its major-symmetric
part, read on its matrix ``M`` ([`operator_matrix`](@ref)):

```math
\\frac{\\lVert M - M^{\\mathsf{T}} \\rVert}{\\lVert M + M^{\\mathsf{T}} \\rVert},
```

the Frobenius norm of the skew part over that of the symmetric part: zero for a
major-symmetric tensor, one when the symmetric and skew parts have the same
norm. It is zero for the null tensor. Apart from that case it involves no
comparison, so that dual numbers carry its derivative.

See also [`definiteness`](@ref).
"""
function major_asymmetry(t::AbstractTens{4})
    M = operator_matrix(t)
    Mᵀ = transpose(M)
    s = norm(M + Mᵀ)
    return iszero(s) ? zero(s) : norm(M - Mᵀ) / s
end

export operator_matrix, from_operator_matrix, definiteness, major_asymmetry
