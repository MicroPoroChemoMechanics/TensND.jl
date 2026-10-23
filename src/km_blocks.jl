# ── Kelvin-Mandel block algebra about a normal ───────────────────────────────
#
# In an orthonormal frame (ℓ, m, n) whose third axis is a normal n, the
# Kelvin-Mandel order of TensND and Tensors.jl, (11, 22, 33, √2·23, √2·13, √2·12),
# reads (ℓℓ, mm, nn, √2·mn, √2·ℓn, √2·ℓm): the in-plane symmetric tensors are the
# slots KM_IP = (1, 2, 6), the out-of-plane ones the slots KM_OP = (3, 4, 5), with
# their √2 weights already carried by the Mandel basis. The two sets interleave,
# which costs nothing: every map below indexes by slot.
#
# The out-of-plane slots (3, 4, 5) are the index pairs (33, 23, 13), so they
# correspond to the acoustic indices in the REVERSED order (3, 2, 1):
# KM(ℂ)[OP_a, OP_b] = s_a s_b K[π(a), π(b)], with π = (3, 2, 1) and s = (1, √2, √2).
# The reversal is invisible on an isotropic tensor (K is then diagonal); only a
# triclinic one tests it.
#
# Everything is written on static matrices with no pivoting and no SVD, so that it
# holds for duals and for symbols: a 3×3 inverse by cofactors (`_inv3`), the
# 6×6 inverse by a Schur complement on the partition.

"""
    KM_IP
    KM_OP

The Kelvin-Mandel slots of the in-plane symmetric tensors, `(1, 2, 6)` for
`(ℓ⊗ℓ, m⊗m, √2 ℓ⊗ˢm)`, and of the out-of-plane ones, `(3, 4, 5)` for
`(n⊗n, √2 m⊗ˢn, √2 ℓ⊗ˢn)`, in an orthonormal frame `(ℓ, m, n)` whose third axis
is the normal `n`. Across the interfaces of a laminate normal to `n`, the
in-plane strain and the out-of-plane stress are continuous.
"""
const KM_IP = SVector(1, 2, 6)
"The Kelvin-Mandel slots of the out-of-plane symmetric tensors; see [`KM_IP`](@ref)."
const KM_OP = SVector(3, 4, 5)

@inline _op_slot(i::Int) = i == 3 ? 1 : (i == 4 ? 2 : (i == 5 ? 3 : 0))
@inline _ip_slot(i::Int) = i == 1 ? 1 : (i == 2 ? 2 : (i == 6 ? 3 : 0))

"""
    km_ip_block(M) → SMatrix{3,3}
    km_op_block(M) → SMatrix{3,3}

The in-plane or out-of-plane block of a 6×6 Kelvin-Mandel matrix, in slot order
([`KM_IP`](@ref), [`KM_OP`](@ref)), with its Mandel weights: no rescaling and no
permutation ([`acoustic_tensor`](@ref) reads the physical acoustic tensor).
"""
@inline km_ip_block(M::AbstractMatrix) = SMatrix{3, 3}(M[KM_IP, KM_IP])
"The out-of-plane block of a 6×6 Kelvin-Mandel matrix; see [`km_ip_block`](@ref)."
@inline km_op_block(M::AbstractMatrix) = SMatrix{3, 3}(M[KM_OP, KM_OP])

"""
    km_ip_embed(B) → SMatrix{6,6}
    km_op_embed(B) → SMatrix{6,6}

A slot-ordered 3×3 block scattered into a 6×6 Kelvin-Mandel matrix, zero
everywhere else: the right inverse of [`km_ip_block`](@ref) or
[`km_op_block`](@ref).
"""
@inline km_ip_embed(B::AbstractMatrix) = _km_embed(_ip_slot, B)
"An out-of-plane block scattered into a 6×6 Kelvin-Mandel matrix; see [`km_ip_embed`](@ref)."
@inline km_op_embed(B::AbstractMatrix) = _km_embed(_op_slot, B)

@inline function _km_embed(slot::F, B::AbstractMatrix{T}) where {F, T}
    return SMatrix{6, 6, T}(
        ntuple(Val(36)) do k
            a, b = slot(mod1(k, 6)), slot((k - 1) ÷ 6 + 1)
            (a == 0 || b == 0) ? zero(T) : @inbounds(B[a, b])
        end
    )
end

"""
    km_from_blocks(X_II, X_IO, X_OI, X_OO) → SMatrix{6,6}

The 6×6 Kelvin-Mandel matrix of its four slot-ordered 3×3 blocks on the
in-plane / out-of-plane partition ([`KM_IP`](@ref), [`KM_OP`](@ref)).
"""
@inline function km_from_blocks(
        X11::AbstractMatrix, X12::AbstractMatrix, X21::AbstractMatrix, X22::AbstractMatrix
    )
    T = promote_type(eltype(X11), eltype(X12), eltype(X21), eltype(X22))
    return SMatrix{6, 6, T}(
        ntuple(Val(36)) do k
            i, j = mod1(k, 6), (k - 1) ÷ 6 + 1
            ai, aj, bi, bj = _ip_slot(i), _ip_slot(j), _op_slot(i), _op_slot(j)
            @inbounds if ai != 0 && aj != 0
                X11[ai, aj]
            elseif ai != 0
                X12[ai, bj]
            elseif aj != 0
                X21[bi, aj]
            else
                X22[bi, bj]
            end
        end
    )
end

# The inverse of a 3×3 matrix by its cofactors and one determinant: no pivot, and
# on symbols an expression of the size of the cofactors (the generic inverse of
# StaticArrays swells symbolic results by orders of magnitude).
@inline function _inv3(K::AbstractMatrix{T}) where {T}
    return @inbounds begin
        a11, a12, a13 = K[1, 1], K[1, 2], K[1, 3]
        a21, a22, a23 = K[2, 1], K[2, 2], K[2, 3]
        a31, a32, a33 = K[3, 1], K[3, 2], K[3, 3]
        c11, c12, c13 = a22 * a33 - a23 * a32, a13 * a32 - a12 * a33, a12 * a23 - a13 * a22
        c21, c22, c23 = a23 * a31 - a21 * a33, a11 * a33 - a13 * a31, a13 * a21 - a11 * a23
        c31, c32, c33 = a21 * a32 - a22 * a31, a12 * a31 - a11 * a32, a11 * a22 - a12 * a21
        d = inv(a11 * c11 + a12 * c21 + a13 * c31)
        SMatrix{3, 3}(c11 * d, c21 * d, c31 * d, c12 * d, c22 * d, c32 * d, c13 * d, c23 * d, c33 * d)
    end
end

"""
    km_inv_schur(M) → SMatrix{6,6}

The inverse of a 6×6 Kelvin-Mandel matrix by block elimination on the in-plane /
out-of-plane partition: two 3×3 inverses by cofactors and a Schur complement,
with no pivoting and no factorization, so that it holds for duals and for
symbols, where an LU factorization pivots on values. Needs the out-of-plane block
and the in-plane Schur complement invertible, as they are for any stiffness whose
acoustic tensor about `n` is positive definite.
"""
@inline function km_inv_schur(M::AbstractMatrix)
    A, D = km_ip_block(M), km_op_block(M)
    B = SMatrix{3, 3}(M[KM_IP, KM_OP])
    C = SMatrix{3, 3}(M[KM_OP, KM_IP])
    Di = _inv3(D)
    BDi, DiC = B * Di, Di * C
    S = _inv3(A - BDi * C)
    return km_from_blocks(S, -S * BDi, -DiC * S, Di + DiC * S * BDi)
end

"""
    plane_pinv(M) → SMatrix{6,6}
    plane_pinv2(K) → SMatrix{3,3}

The Moore-Penrose pseudo-inverse of the out-of-plane part of a 6×6
Kelvin-Mandel matrix (`plane_pinv`): the inverse of its [`km_op_block`](@ref),
embedded back. For the Kelvin-Mandel matrix of a stiffness ℂ it is the Hill
tensor of a flat inclusion of normal `n`, `n ⊗ˢ (n⋅ℂ⋅n)⁻¹ ⊗ˢ n`; for a matrix
supported on the out-of-plane subspace it is its pseudo-inverse. `plane_pinv2` is
the order-2 counterpart for a 3×3 matrix written in the frame `(ℓ, m, n)`:
`(1/K₃₃) n⊗n`. No SVD, so duals and symbols go through.
"""
@inline plane_pinv(M::AbstractMatrix) = km_op_embed(_inv3(km_op_block(M)))
"The order-2 out-of-plane pseudo-inverse `(1/K₃₃) n⊗n`; see [`plane_pinv`](@ref)."
@inline function plane_pinv2(K3::AbstractMatrix{T}) where {T}
    z = zero(T)
    return SMatrix{3, 3, T}(z, z, z, z, z, z, z, z, inv(@inbounds K3[3, 3]))
end

"""
    acoustic_tensor(M::AbstractMatrix) → SMatrix{3,3}
    acoustic_tensor(ℂ::AbstractTens{4,3}, n) → Tens{2,3}

The acoustic tensor `K = n⋅ℂ⋅n`, `Kᵢⱼ = nₖ ℂₖᵢₗⱼ nₗ`: read off the out-of-plane
block of a 6×6 Kelvin-Mandel matrix written in a frame whose third axis is `n`
(its Mandel weights and the (3, 2, 1) order of its slots undone), or computed
from a tensor and a direction.
"""
@inline function acoustic_tensor(M::AbstractMatrix{T}) where {T}
    B = km_op_block(M)
    r2 = sqrt(2 * one(T))
    K11, K12, K13 = B[3, 3] / 2, B[3, 2] / 2, B[3, 1] / r2
    K22, K23, K33 = B[2, 2] / 2, B[2, 1] / r2, B[1, 1]
    return SMatrix{3, 3, T}(K11, K12, K13, K12, K22, K23, K13, K23, K33)
end
acoustic_tensor(ℂ::AbstractTens{4, 3}, n) = dotdot(Tens(Vec{3}(_extract_vec(n))), ℂ, Tens(Vec{3}(_extract_vec(n))))

"""
    compliance_op_block(X) → SMatrix{3,3}

The slot-ordered out-of-plane block ([`km_op_block`](@ref)) of the fourth-order
compliance `X ⊗ˢ n ⊗ˢ n` — the strain `(X⋅(σ⋅n)) ⊗ˢ n` — from the symmetric 3×3
matrix `X` written in the frame `(ℓ, m, n)`: the map
[`acoustic_tensor`](@ref) inverts on the compliance side,
`Ŝₐᵦ = X_{π(a)π(b)} / (sₐ sᵦ)`, `π = (3, 2, 1)`, `s = (1, √2, √2)`. For a spring
interface of compliance `X`, `km_op_embed(compliance_op_block(X))` is its
Kelvin-Mandel matrix.
"""
@inline function compliance_op_block(X::AbstractMatrix{T}) where {T}
    r2 = sqrt(2 * one(T))
    S11, S12, S13 = X[3, 3], X[3, 2] / r2, X[3, 1] / r2
    S22, S23, S33 = X[2, 2] / 2, X[2, 1] / 2, X[1, 1] / 2
    return SMatrix{3, 3, T}(S11, S12, S13, S12, S22, S23, S13, S23, S33)
end

export KM_IP, KM_OP, km_ip_block, km_op_block, km_ip_embed, km_op_embed, km_from_blocks
export km_inv_schur, plane_pinv, plane_pinv2, acoustic_tensor, compliance_op_block
