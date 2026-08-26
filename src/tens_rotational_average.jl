# =============================================================================
#  tens_rotational_average.jl — EXACT rotation-group averages of 2nd/4th-order
#  tensors with minor symmetry only (no major symmetry assumed).
#
#  Two operations, mirroring a runtime tensor symmetrization:
#
#  * `isotropify(t)`            — exact average over SO(3):
#        4th order → TensISO{4} (the isotropic subspace of minor-symmetric
#        tensors is 2-dimensional {𝕁, 𝕂}, even without major symmetry)
#        2nd order → TensISO{2}
#
#  * `transverse_isotropify(t, n)` — exact average over rotations about `n`:
#        4th order → TensTI{4,T,8} (full 8-dim commutant of the SO(2) action
#        on the Kelvin-Mandel space: Walpole W₁..W₆ + antisymmetric couplings
#        W₇ (m=1) and W₈ (m=2))
#        2nd order → TensTI{2,T,3} (a·nT + b·nₙ + c·w, w the in-plane
#        rotation generator w·p = n×p)
#
#  These are the operators to use on CONCENTRATION / CONTRIBUTION tensors —
#  such tensors generally lack major symmetry, and the average preserves that
#  (ℓ₃ ≠ ℓ₄, ℓ₇, ℓ₈). For parameter extraction or reporting, use the best-fit
#  projections instead (`best_fit_ti`, `best_fit_iso`, `best_fit_ortho`).
#
#  The two are NOT the same operation and must never be conflated: an average
#  is exact and lossless on the invariant subspace, a projection is a
#  least-squares fit that discards whatever does not fit.
#
#  Kelvin-Mandel convention (Tensors.jl): index order (11, 22, 33, 23, 13, 12)
#  with weight √2 on the shear slots. The closed-form azimuthal average about
#  e₃ is the orthogonal projection onto the commutant algebra:
#      m=0 block {33, (11+22)/√2}      → kept in full (ℓ₁, ℓ₂, ℓ₃, ℓ₄)
#      m=1 doublet {23, 13}            → aI + bG → (ℓ₆, ℓ₇)
#      m=2 doublet {(11−22)/√2, 12}    → aI + bG → (ℓ₅, ℓ₈)
# =============================================================================

# ── Kelvin-Mandel helpers (minor-symmetric, Tensors.jl index convention) ─────

const _MANDEL_IDX = ((1, 1), (2, 2), (3, 3), (2, 3), (1, 3), (1, 2))

"""
    mandel66_minor(arr::AbstractArray{T,4}) → Matrix{T} (6×6)

Kelvin-Mandel 6×6 matrix of a 4th-order array, minor-symmetrizing on read
(exactly, so tensors that are minor-symmetric up to round-off are cleanly
projected).  Index convention (11, 22, 33, 23, 13, 12), weights √2.

Unlike [`KM`](@ref) this makes no assumption of major symmetry and never goes
through a `Tensors.SymmetricTensor`, so it accepts a concentration tensor.
"""
function mandel66_minor(arr::AbstractArray{T, 4}) where {T}
    sq2 = sqrt(T(2))
    Λ(I) = I ≤ 3 ? one(T) : sq2
    M = Matrix{T}(undef, 6, 6)
    @inbounds for I in 1:6, J in 1:6
        (i, j) = _MANDEL_IDX[I]
        (k, l) = _MANDEL_IDX[J]
        v = (arr[i, j, k, l] + arr[j, i, k, l] + arr[i, j, l, k] + arr[j, i, l, k]) / 4
        M[I, J] = Λ(I) * Λ(J) * v
    end
    return M
end

"""
    array_from_mandel66(M::AbstractMatrix{T}) → Array{T,4}

Inverse of [`mandel66_minor`](@ref): rebuild the (exactly minor-symmetric)
3×3×3×3 array from a 6×6 Kelvin-Mandel matrix.
"""
function array_from_mandel66(M::AbstractMatrix{T}) where {T}
    sq2 = sqrt(T(2))
    Λ(I) = I ≤ 3 ? one(T) : sq2
    arr = Array{T, 4}(undef, 3, 3, 3, 3)
    @inbounds for I in 1:6, J in 1:6
        (i, j) = _MANDEL_IDX[I]
        (k, l) = _MANDEL_IDX[J]
        v = M[I, J] / (Λ(I) * Λ(J))
        arr[i, j, k, l] = v
        arr[j, i, k, l] = v
        arr[i, j, l, k] = v
        arr[j, i, l, k] = v
    end
    return arr
end

"""
    _axis_frame(n) → (R, n̂)

Return an orthonormal, right-handed frame `(e₁′, e₂′, n̂)` as the columns of a
3×3 matrix, with `n̂` the normalized axis.  The in-plane vectors are chosen
deterministically; the azimuthal average is independent of that choice.

For a comparable element type the reference vector is the canonical axis least
aligned with `n`, which is the best-conditioned choice. A symbolic axis answers
no comparison (`is_hard_numeric` is false), so the reference is then picked
structurally, by the first canonical axis whose cross product with `n` is not
*identically* zero — no `argmin`, no tolerance.
"""
function _axis_frame(n)
    T = promote_type(typeof(n[1]), typeof(n[2]), typeof(n[3]))
    nn = sqrt(n[1]^2 + n[2]^2 + n[3]^2)
    n̂ = (T(n[1] / nn), T(n[2] / nn), T(n[3] / nn))
    h = _inplane_reference(n̂, T)
    # u = h × n̂ (in-plane), v = n̂ × u  →  (u, v, n̂) right-handed
    u1 = h[2] * n̂[3] - h[3] * n̂[2]
    u2 = h[3] * n̂[1] - h[1] * n̂[3]
    u3 = h[1] * n̂[2] - h[2] * n̂[1]
    un = sqrt(u1^2 + u2^2 + u3^2)
    u1, u2, u3 = u1 / un, u2 / un, u3 / un
    v1 = n̂[2] * u3 - n̂[3] * u2
    v2 = n̂[3] * u1 - n̂[1] * u3
    v3 = n̂[1] * u2 - n̂[2] * u1
    R = Matrix{T}(undef, 3, 3)
    R[1, 1], R[2, 1], R[3, 1] = u1, u2, u3
    R[1, 2], R[2, 2], R[3, 2] = v1, v2, v3
    R[1, 3], R[2, 3], R[3, 3] = n̂[1], n̂[2], n̂[3]
    return R, n̂
end

function _inplane_reference(n̂::NTuple{3, T}, ::Type{T}) where {T}
    if is_hard_numeric(T)
        a1, a2, a3 = abs(n̂[1]), abs(n̂[2]), abs(n̂[3])
        return a1 ≤ a2 ? (a1 ≤ a3 ? (one(T), zero(T), zero(T)) : (zero(T), zero(T), one(T))) :
            (a2 ≤ a3 ? (zero(T), one(T), zero(T)) : (zero(T), zero(T), one(T)))
    end
    # Symbolic: `h × n̂ = 0` is decided structurally (both other components
    # identically zero), never by a tolerance. Falls back to e₃, which is the
    # right answer whenever `n` is an in-plane symbolic direction.
    iszero(n̂[2]) && iszero(n̂[3]) && return (zero(T), one(T), zero(T))
    return (one(T), zero(T), zero(T))
end

# Components of a 4th-order array in the frame whose COLUMNS are given by R:
# arr′[a,b,c,d] = R[i,a] R[j,b] R[k,c] R[l,d] arr[i,j,k,l], via four successive
# single-index contractions (cost 4·3⁵ instead of 3⁸).
function _rotate4(arr::AbstractArray{TA, 4}, R::AbstractMatrix{TR}) where {TA, TR}
    T = promote_type(TA, TR)
    t1 = zeros(T, 3, 3, 3, 3)
    @inbounds for a in 1:3, j in 1:3, k in 1:3, l in 1:3
        s = zero(T)
        for i in 1:3
            s += R[i, a] * arr[i, j, k, l]
        end
        t1[a, j, k, l] = s
    end
    t2 = zeros(T, 3, 3, 3, 3)
    @inbounds for a in 1:3, b in 1:3, k in 1:3, l in 1:3
        s = zero(T)
        for j in 1:3
            s += R[j, b] * t1[a, j, k, l]
        end
        t2[a, b, k, l] = s
    end
    t1 .= 0
    @inbounds for a in 1:3, b in 1:3, c in 1:3, l in 1:3
        s = zero(T)
        for k in 1:3
            s += R[k, c] * t2[a, b, k, l]
        end
        t1[a, b, c, l] = s
    end
    t2 .= 0
    @inbounds for a in 1:3, b in 1:3, c in 1:3, d in 1:3
        s = zero(T)
        for l in 1:3
            s += R[l, d] * t1[a, b, c, l]
        end
        t2[a, b, c, d] = s
    end
    return t2
end

# ── Closed-form azimuthal average about e₃ (Kelvin-Mandel, Tensors order) ────

"""
    ti8_params_from_KM(M::AbstractMatrix) → NTuple{8}

Coefficients `(ℓ₁, …, ℓ₈)` of the exact azimuthal average about `e₃` of a
minor-symmetric tensor given by its 6×6 Kelvin-Mandel matrix (Tensors.jl index
order 11, 22, 33, 23, 13, 12). Orthogonal projection onto the commutant
algebra; TensND Walpole convention (ℓ₃ ↔ the `C₃₃₁₁` side, ℓ₄ ↔ the `C₁₁₃₃`
side).

Two distinct uses, and the difference matters:

- as an **average**, on a matrix that is not axially invariant — this is what
  [`transverse_isotropify`](@ref) calls;
- as an exact **read-off**, on a matrix already known to be axially invariant
  about `e₃`, where the projection is the identity and the eight coefficients
  are simply recovered. This is the non-major-symmetric counterpart of
  [`ti_params_from_KM`](@ref), which projects onto the 5-coefficient
  major-symmetric span and would silently discard `ℓ₃ ≠ ℓ₄`, `ℓ₇` and `ℓ₈` —
  precisely the content a concentration tensor carries.

Reciprocal: [`KM_from_ti8_params`](@ref).
"""
function ti8_params_from_KM(M::AbstractMatrix{T}) where {T}
    sq2 = sqrt(T(2))
    ℓ₁ = M[3, 3]
    ℓ₂ = (M[1, 1] + M[2, 2] + M[1, 2] + M[2, 1]) / 2
    ℓ₃ = (M[3, 1] + M[3, 2]) / sq2
    ℓ₄ = (M[1, 3] + M[2, 3]) / sq2
    ℓ₅ = (M[1, 1] + M[2, 2] - M[1, 2] - M[2, 1]) / 4 + M[6, 6] / 2
    ℓ₆ = (M[4, 4] + M[5, 5]) / 2
    ℓ₇ = (M[5, 4] - M[4, 5]) / 2
    ℓ₈ = (M[6, 1] - M[6, 2] - M[1, 6] + M[2, 6]) / (2 * sq2)
    return (ℓ₁, ℓ₂, ℓ₃, ℓ₄, ℓ₅, ℓ₆, ℓ₇, ℓ₈)
end

"""
    KM_from_ti8_params(p::NTuple{8,T}) → Matrix{T} (6×6)

Kelvin-Mandel matrix (Tensors order, axis `e₃`) of the axially-invariant tensor
with coefficients `(ℓ₁,…,ℓ₈)`. Inverse of [`ti8_params_from_KM`](@ref) on the
commutant subspace, and the 8-coefficient counterpart of
[`KM_from_ti_params`](@ref).
"""
function KM_from_ti8_params(p::NTuple{8, T}) where {T}
    ℓ₁, ℓ₂, ℓ₃, ℓ₄, ℓ₅, ℓ₆, ℓ₇, ℓ₈ = p
    sq2 = sqrt(T(2))
    M = zeros(T, 6, 6)
    M[1, 1] = M[2, 2] = (ℓ₂ + ℓ₅) / 2
    M[1, 2] = M[2, 1] = (ℓ₂ - ℓ₅) / 2
    M[3, 3] = ℓ₁
    M[3, 1] = M[3, 2] = ℓ₃ / sq2
    M[1, 3] = M[2, 3] = ℓ₄ / sq2
    M[4, 4] = M[5, 5] = ℓ₆
    M[4, 5] = -ℓ₇
    M[5, 4] = ℓ₇
    M[6, 6] = ℓ₅
    M[6, 1] = ℓ₈ / sq2
    M[6, 2] = -ℓ₈ / sq2
    M[1, 6] = -ℓ₈ / sq2
    M[2, 6] = ℓ₈ / sq2
    return M
end

# ── Public API — 4th order ───────────────────────────────────────────────────

"""
    isotropify(t::AbstractTens{4,3}) → TensISO{4}

Exact average of `t` over SO(3): `α = T_iijj/3`, `β = (T_ijij − α)/5` → `α𝕁 + β𝕂`.
Valid for minor-symmetric tensors with or without major symmetry (the isotropic
subspace is `{𝕁, 𝕂}` in both cases), which is what distinguishes it from the
`AbstractArray` method above — the latter is a Frobenius projection and the two
coincide exactly on minor-symmetric input.
"""
function isotropify(t::AbstractTens{4, 3})
    arr = get_array(t)
    α = sum(arr[i, i, j, j] for i in 1:3, j in 1:3) / 3
    full_trace = sum(arr[i, j, i, j] for i in 1:3, j in 1:3)
    β = (full_trace - α) / 5
    return TensISO{3}(α, β)
end

# Closed form on the Walpole coefficients, avoiding `get_array`.
#
# With `nₙ = n⊗n`, `nT = 𝟏 − nₙ` and the Walpole basis (W₁ = nₙ⊗nₙ,
# W₂ = (nT⊗nT)/2, W₃ = (nₙ⊗nT)/√2, W₄ = (nT⊗nₙ)/√2, W₅ = nT⊠ˢnT − (nT⊗nT)/2,
# W₆ = nT⊠ˢnₙ + nₙ⊠ˢnT), the two invariants read
#
#     T_iijj = ℓ₁ + 2ℓ₂ + √2 (ℓ₃ + ℓ₄)
#     T_ijij = ℓ₁ +  ℓ₂ + 2ℓ₅ + 2ℓ₆
#
# Checked on the two tensors whose isotropic average is known by inspection:
# 𝕀 = W₁+W₂+W₅+W₆ gives α = β = 1, and 𝕁 = (W₁ + 2W₂ + √2W₃ + √2W₄)/3 gives
# α = 1, β = 0.
#
# `ℓ₇` and `ℓ₈` (the antisymmetric azimuthal couplings of the `N = 8` form)
# contribute to neither invariant — W₇ and W₈ annihilate `𝟏` and are purely
# off-diagonal in Kelvin-Mandel, hence traceless — so `get_ℓ`, which drops
# them, is the right accessor here.
#
# Motivation: `get_array` rebuilds all 81 components from the Walpole basis,
# which on a symbolic element type produces expressions an order of magnitude
# larger than the answer. This path also speeds up a numeric polycrystal
# self-consistent scheme, where the average runs once per phase per iteration.
function isotropify(t::TensTI{4})
    ℓ₁, ℓ₂, ℓ₃, ℓ₄, ℓ₅, ℓ₆ = get_ℓ(t)
    r2 = sqrt(2 * one(eltype(t)))
    α = (ℓ₁ + 2ℓ₂ + r2 * (ℓ₃ + ℓ₄)) / 3
    β = (ℓ₁ + ℓ₂ + 2ℓ₅ + 2ℓ₆ - α) / 5
    return TensISO{3}(α, β)
end

# An isotropic tensor is its own SO(3) average, at either order.
isotropify(t::TensISO{4, 3}) = t
isotropify(t::TensISO{2, 3}) = t

# Order 2, same idea: `t = a nT + b nₙ (+ c w)` has trace `2a + b`, the
# antisymmetric in-plane part `w` being traceless. No `get_array`.
function isotropify(t::TensTI{2})
    d = get_data(t)
    return TensISO{3}((2 * d[1] + d[2]) / 3)
end

"""
    isotropify(t::AbstractTens{2,3}) → TensISO{2}

Exact SO(3) average of a 2nd-order tensor: `(tr t / 3) 𝟏`.
"""
function isotropify(t::AbstractTens{2, 3})
    arr = get_array(t)
    λ = (arr[1, 1] + arr[2, 2] + arr[3, 3]) / 3
    return TensISO{3}(λ)
end

"""
    transverse_isotropify(t::AbstractTens{4,3}, n) → TensTI{4,T,8}

Exact average of `t` over all rotations about the axis `n`
(`(1/2π)∫ R_φ ⋆ t dφ`), for a minor-symmetric `t` with or without major
symmetry.  The result lives in the full 8-dimensional axially-invariant space:
the non-major-symmetric components (ℓ₃ ≠ ℓ₄) and the antisymmetric azimuthal
couplings (ℓ₇, ℓ₈) are preserved — they are dropped by naive symmetric TI
projections but present in e.g. averaged strain-concentration tensors.
"""
function transverse_isotropify(t::AbstractTens{4, 3}, n)
    R, n̂ = _axis_frame(n)
    # Structured fast paths (already axially invariant about n̂)
    if t isa TensISO{4, 3}
        return _lift_walpole_N8(fromISO(t, n̂))
    elseif t isa TensTI{4} && axis(t) == n̂
        return _lift_walpole_N8(t)
    end
    arr = get_array(t)
    arr_ez = _rotate4(arr, R)
    p = ti8_params_from_KM(mandel66_minor(arr_ez))
    return TensTI{4}(p..., n̂)
end

"""
    transverse_isotropify(t::AbstractTens{2,3}, n) → TensTI{2,T,3}

Exact azimuthal average of a 2nd-order tensor about `n`:
`a·nT + b·nₙ + c·w` with `b = n̂ᵀ t n̂`, `a = (tr t − b)/2` and `c = (w : t)/2`
(`w` the in-plane rotation generator `w·p = n̂ × p`).  The antisymmetric
in-plane part `c` is preserved (a symmetric TI parametrization would silently
drop it).
"""
function transverse_isotropify(t::AbstractTens{2, 3}, n)
    arr = get_array(t)
    T0 = eltype(arr)
    nn = sqrt(n[1]^2 + n[2]^2 + n[3]^2)
    T = promote_type(T0, typeof(nn))
    n̂ = (T(n[1] / nn), T(n[2] / nn), T(n[3] / nn))
    b = zero(T)
    @inbounds for i in 1:3, j in 1:3
        b += n̂[i] * arr[i, j] * n̂[j]
    end
    a = ((arr[1, 1] + arr[2, 2] + arr[3, 3]) - b) / 2
    # c = (1/2) Σ w[i,j] t[i,j],  w[i,j] = ε[i,k,j] n̂[k]
    c = (
        n̂[1] * (arr[3, 2] - arr[2, 3]) +
            n̂[2] * (arr[1, 3] - arr[3, 1]) +
            n̂[3] * (arr[2, 1] - arr[1, 2])
    ) / 2
    return TensTI{2}(a, b, c, n̂)
end

# ── Kelvin-Mandel block helpers ─────────────────────────────────────────────
#
# For callers working directly on 6×6 blocks rather than on tensors — the
# ageing-viscoelastic Volterra matrices of MeanFieldHomogenization.jl are the
# motivating case.

"""
    ti_average_mandel66(M::AbstractMatrix, n) → Matrix (6×6)

Exact azimuthal average about `n` of the minor-symmetric tensor whose 6×6
Kelvin-Mandel matrix (Tensors.jl order) is `M`, returned as a 6×6 matrix.
"""
function ti_average_mandel66(M::AbstractMatrix, n)
    R, n̂ = _axis_frame(n)
    arr_ez = _rotate4(array_from_mandel66(M), R)
    p = ti8_params_from_KM(mandel66_minor(arr_ez))
    M_ez = KM_from_ti8_params(p)
    # rotate back to the canonical frame: columns of R are (e₁′, e₂′, n̂), so
    # the inverse basis change uses Rᵀ.
    Rt = permutedims(R)
    return mandel66_minor(_rotate4(array_from_mandel66(M_ez), Rt))
end

"""
    iso_average_mandel66(M::AbstractMatrix{T}) → (α, β)

Exact SO(3) average of the minor-symmetric tensor whose Kelvin-Mandel matrix is
`M`: `α = (1/3)Σ_{i,j≤3} M[i,j]`, `β = (tr M − α)/5`.  Returns the `(α, β)`
coefficients of `α𝕁 + β𝕂`.
"""
function iso_average_mandel66(M::AbstractMatrix{T}) where {T}
    α = (
        M[1, 1] + M[1, 2] + M[1, 3] + M[2, 1] + M[2, 2] +
            M[2, 3] + M[3, 1] + M[3, 2] + M[3, 3]
    ) / 3
    β = (M[1, 1] + M[2, 2] + M[3, 3] + M[4, 4] + M[5, 5] + M[6, 6] - α) / 5
    return (α, β)
end

# ── Exports ──────────────────────────────────────────────────────────────────

export transverse_isotropify
export mandel66_minor, array_from_mandel66
export ti8_params_from_KM, KM_from_ti8_params
export ti_average_mandel66, iso_average_mandel66
