##############################################################################
# Tensor projection onto symmetry subspaces (TI, ORTHO)                     #
#                                                                            #
# Extends proj_tens (ISO already in tens_isotropic.jl) to:                   #
#   - TI  (transverse isotropy)  — fixed axis or optimized                   #
#   - ORTHO (orthotropy)         — fixed frame or optimized                  #
# for both 4th-order and 2nd-order tensors.                                  #
#                                                                            #
# Inspired by ECHOES C++ (tensor_approx.h, tensor_ti.h, tensor_ortho.h).    #
##############################################################################

using StaticArrays

# ── Kelvin-Mandel index couples ──────────────────────────────────────────────
# Kelvin-Mandel (KM) ordering maps symmetric 3×3 tensor index pairs (i,j)
# to a single index k ∈ 1:6.  The convention is:
#   k=1↔(1,1), k=2↔(2,2), k=3↔(3,3), k=4↔(2,3), k=5↔(1,3), k=6↔(1,2)
# Off-diagonal components are scaled by √2 so that the Frobenius norm is
# preserved (‖A‖² = Aₖₘ Aₖₘ).  This differs from Voigt notation which
# uses engineering shear (factor 2 instead of √2).
const _KM_COUPLES = ((1, 1), (2, 2), (3, 3), (2, 3), (1, 3), (1, 2))

# ── ForwardDiff-compatible rotation helpers ──────────────────────────────────

"""
    _rot3_raw(θ, ϕ, ψ) → SMatrix{3,3}

Rotation matrix (RotZYZ convention: R = Rz(ϕ)·Ry(θ)·Rz(ψ)) built from
explicit cos/sin, compatible with ForwardDiff Dual numbers.

The third column of the matrix is the direction `(sinθ cosϕ, sinθ sinϕ, cosθ)`.

# Examples
```julia
julia> R = _rot3_raw(0.3, 0.5, 0.0)
3×3 SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
 ...
```
"""
function _rot3_raw(θ, ϕ, ψ)
    cθ, sθ = cos(θ), sin(θ)
    cϕ, sϕ = cos(ϕ), sin(ϕ)
    cψ, sψ = cos(ψ), sin(ψ)
    T = promote_type(typeof(cθ), typeof(cϕ), typeof(cψ))
    return SMatrix{3, 3, T}(
        cθ * cψ * cϕ - sψ * sϕ, cθ * cψ * sϕ + cϕ * sψ, -cψ * sθ,
        -cθ * cϕ * sψ - cψ * sϕ, -cθ * sψ * sϕ + cψ * cϕ, sθ * sψ,
        cϕ * sθ, sθ * sϕ, cθ
    )
end

"""
    _KM_rotation(θ, ϕ, ψ) → SMatrix{6,6}

6×6 Kelvin-Mandel (Bond) rotation matrix from Euler angles (θ, ϕ, ψ).
Compatible with ForwardDiff.

Transforms a 6×6 KM matrix as: `C' = Q' · C · Q` where `Q = _KM_rotation(...)`.

Uses the standard Bond matrix construction from `R = _rot3_raw(θ, ϕ, ψ)`.

# Examples
```julia
julia> Q = _KM_rotation(0.3, 0.5, 0.0);

julia> size(Q)
(6, 6)
```
"""
function _KM_rotation(θ, ϕ, ψ)
    R = _rot3_raw(θ, ϕ, ψ)
    T = eltype(R)
    sq2 = sqrt(T(2))
    Q = MMatrix{6, 6, T}(undef)

    # (I,J) both in 1:3 — diagonal block: Q[I,J] = R[I,J]²
    for I in 1:3, J in 1:3
        Q[I, J] = R[I, J]^2
    end

    # I ∈ 1:3, J ∈ 4:6 — normal-shear cross
    for J in 4:6
        k, l = _KM_COUPLES[J]
        for I in 1:3
            Q[I, J] = sq2 * R[I, k] * R[I, l]
        end
    end

    # I ∈ 4:6, J ∈ 1:3 — shear-normal cross
    for I in 4:6
        i, j = _KM_COUPLES[I]
        for J in 1:3
            Q[I, J] = sq2 * R[i, J] * R[j, J]
        end
    end

    # (I,J) both in 4:6 — shear-shear block
    for I in 4:6
        i, j = _KM_COUPLES[I]
        for J in 4:6
            k, l = _KM_COUPLES[J]
            Q[I, J] = R[i, k] * R[j, l] + R[i, l] * R[j, k]
        end
    end

    return SMatrix{6, 6, T}(Q)
end

"""
    _KM_of_array(A::AbstractArray{T,4}) → Matrix{T}

Compute the 6×6 Kelvin-Mandel matrix of a 3×3×3×3 array.
"""
function _KM_of_array(A::AbstractArray{T, 4}) where {T}
    return tomandel(tensor_or_array(A))
end

"""
    _KM_of_array(A::AbstractArray{T,2}) → Matrix{T}

Compute the 6-vector (or just return the 3×3 matrix) Kelvin-Mandel for a 2nd-order tensor.
"""
function _KM_of_array(A::AbstractArray{T, 2}) where {T}
    return A
end

# ── TI projection helpers ────────────────────────────────────────────────────

"""
    _project_TI_KM(C::AbstractMatrix) → (ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆)

Project a 6×6 KM matrix (assumed to be in a frame where `e₃` is the TI symmetry
axis) onto the major-symmetric TI subspace. Returns the 5 Walpole coefficients.

Formulas (from ECHOES `tensor_ti.h`):
- `c = (C[1,1] + C[2,2]) / 2`
- `d = (C[1,2] + C[2,1]) / 2`
- `ℓ₁ = C[3,3]`
- `ℓ₂ = c + d`
- `ℓ₃ = (C[1,3] + C[2,3] + C[3,1] + C[3,2]) / (2√2)`
- `ℓ₅ = (c − d + C[6,6]) / 2`
- `ℓ₆ = (C[4,4] + C[5,5]) / 2`

# Examples
```julia
julia> C = Float64[10 3 2.5 0 0 0; 3 10 2.5 0 0 0; 2.5 2.5 12 0 0 0;
                   0 0 0 4 0 0; 0 0 0 0 4 0; 0 0 0 0 0 7];

julia> _project_TI_KM(C)
(12.0, 13.0, 3.5355339059327378, 7.0, 4.0)
```
"""
function _project_TI_KM(C)
    T = eltype(C)
    sq2 = sqrt(T(2))
    c = (C[1, 1] + C[2, 2]) / 2
    d = (C[1, 2] + C[2, 1]) / 2
    ℓ₁ = C[3, 3]
    ℓ₂ = c + d
    ℓ₃ = (C[1, 3] + C[2, 3] + C[3, 1] + C[3, 2]) / (2 * sq2)
    ℓ₅ = (c - d + C[6, 6]) / 2
    ℓ₆ = (C[4, 4] + C[5, 5]) / 2
    return (ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆)
end

"""
    _build_TI_KM(ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆) → SMatrix{6,6}

Build a 6×6 KM matrix (in the frame where `e₃` is the TI axis) from
5 major-symmetric Walpole coefficients.

# Examples
```julia
julia> B = _build_TI_KM(12.0, 13.0, sqrt(2)*2.5, 7.0, 4.0);

julia> B[1,1]  # = (ℓ₂ + ℓ₅)/2 = (13+7)/2 = 10
10.0
```
"""
function _build_TI_KM(ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆)
    T = promote_type(typeof(ℓ₁), typeof(ℓ₂), typeof(ℓ₃), typeof(ℓ₅), typeof(ℓ₆))
    sq2 = sqrt(T(2))
    C1111 = (ℓ₂ + ℓ₅) / 2
    C1122 = (ℓ₂ - ℓ₅) / 2
    C1133 = ℓ₃ / sq2
    C3333 = ℓ₁
    C2323 = ℓ₆ / 2
    C1212 = ℓ₅ / 2
    z = zero(T)
    return SMatrix{6, 6, T}(
        C1111, C1122, C1133, z, z, z,
        C1122, C1111, C1133, z, z, z,
        C1133, C1133, C3333, z, z, z,
        z, z, z, 2C2323, z, z,
        z, z, z, z, 2C2323, z,
        z, z, z, z, z, 2C1212
    )
end

# ── ORTHO projection helpers ─────────────────────────────────────────────────

"""
    _project_ORTHO_KM(C::AbstractMatrix) → NTuple{9}

Extract 9 orthotropic parameters from a 6×6 KM matrix in the material frame.
Returns `(C₁₁, C₂₂, C₃₃, C₁₂, C₁₃, C₂₃, C₄₄, C₅₅, C₆₆)`.

For non-symmetric KM matrices (non-major-symmetric tensors), the off-diagonal
entries are averaged: `Cᵢⱼ = (C[i,j] + C[j,i]) / 2`.

# Examples
```julia
julia> C = Float64[10 3 2.5 0 0 0; 3 8 1.5 0 0 0; 2.5 1.5 12 0 0 0;
                   0 0 0 4 0 0; 0 0 0 0 6 0; 0 0 0 0 0 7];

julia> _project_ORTHO_KM(C)
(10.0, 8.0, 12.0, 3.0, 2.5, 1.5, 2.0, 3.0, 3.5)
```
"""
function _project_ORTHO_KM(C)
    C₁₁ = C[1, 1]
    C₂₂ = C[2, 2]
    C₃₃ = C[3, 3]
    C₁₂ = (C[1, 2] + C[2, 1]) / 2
    C₁₃ = (C[1, 3] + C[3, 1]) / 2
    C₂₃ = (C[2, 3] + C[3, 2]) / 2
    C₄₄ = C[4, 4] / 2
    C₅₅ = C[5, 5] / 2
    C₆₆ = C[6, 6] / 2
    return (C₁₁, C₂₂, C₃₃, C₁₂, C₁₃, C₂₃, C₄₄, C₅₅, C₆₆)
end

"""
    _build_ORTHO_KM(C₁₁, C₂₂, C₃₃, C₁₂, C₁₃, C₂₃, C₄₄, C₅₅, C₆₆) → SMatrix{6,6}

Build a 6×6 KM matrix in the material frame from 9 orthotropic parameters.

# Examples
```julia
julia> B = _build_ORTHO_KM(10.0, 8.0, 12.0, 3.0, 2.5, 1.5, 2.0, 3.0, 3.5);

julia> B[4,4]  # = 2*C₄₄ = 4
4.0
```
"""
function _build_ORTHO_KM(C₁₁, C₂₂, C₃₃, C₁₂, C₁₃, C₂₃, C₄₄, C₅₅, C₆₆)
    T = promote_type(
        typeof(C₁₁), typeof(C₂₂), typeof(C₃₃), typeof(C₁₂),
        typeof(C₁₃), typeof(C₂₃), typeof(C₄₄), typeof(C₅₅), typeof(C₆₆)
    )
    z = zero(T)
    return SMatrix{6, 6, T}(
        C₁₁, C₁₂, C₁₃, z, z, z,
        C₁₂, C₂₂, C₂₃, z, z, z,
        C₁₃, C₂₃, C₃₃, z, z, z,
        z, z, z, 2 * C₄₄, z, z,
        z, z, z, z, 2 * C₅₅, z,
        z, z, z, z, z, 2 * C₆₆
    )
end

# ── Norm helpers ─────────────────────────────────────────────────────────────

_frobenius(A::AbstractArray) = sqrt(sum(x -> x^2, A))

# ── Angles from a unit vector ────────────────────────────────────────────────

"""
    _angles_from_n(n) → (θ, ϕ)

Extract spherical angles (θ, ϕ) from a unit vector `n`:
- `θ = atan(√(n₁² + n₂²), n₃)` (polar angle from z-axis)
- `ϕ = atan(n₂, n₁)` (azimuthal angle)

# Examples
```julia
julia> _angles_from_n((0.0, 0.0, 1.0))
(0.0, 0.0)

julia> _angles_from_n((1.0, 0.0, 0.0))
(1.5707963267948966, 0.0)
```
"""
function _angles_from_n(n)
    nv = _extract_vec(n)
    θ = atan(sqrt(nv[1]^2 + nv[2]^2), nv[3])
    ϕ = atan(nv[2], nv[1])
    return (θ, ϕ)
end

"""
    _n_from_angles(θ, ϕ) → NTuple{3}

Build a unit vector from spherical angles:
`n = (sinθ·cosϕ, sinθ·sinϕ, cosθ)`.

# Examples
```julia
julia> _n_from_angles(0.0, 0.0)
(0.0, 0.0, 1.0)
```
"""
function _n_from_angles(θ, ϕ)
    sθ = sin(θ)
    return (sθ * cos(ϕ), sθ * sin(ϕ), cos(θ))
end

# ── Candidate axis / frame from eigenstructure (cheap, O(1)) ─────────────────
#
# For a 4th-order elasticity tensor `C`, the 2nd-order trace
#   dᵢⱼ = Cᵢₖⱼₖ = ∑ₖ C[i,k,j,k]
# captures the principal-direction structure (Cowin-Mehrabadi 1990).  Its
# eigenvectors are candidate principal axes for TI / ORTHO detection without
# running any optimisation.

"""
    _trace_tensor(C::AbstractArray{T,4}) → SMatrix{3,3,T,9}

Return the 2nd-order "trace" tensor `dᵢⱼ = ∑ₖ C[i,k,j,k]` — its eigenstructure
exposes the principal directions of `C` cheaply (no optimisation required).

Used to propose candidate TI axes and candidate orthotropic frames when the
user asks for symmetry detection without angle optimisation.
"""
function _trace_tensor(C::AbstractArray{T, 4}) where {T}
    d = MMatrix{3, 3, T}(undef)
    for i in 1:3, j in 1:3
        s = zero(T)
        for k in 1:3
            s += C[i, k, j, k]
        end
        d[i, j] = s
    end
    return SMatrix{3, 3, T}((d + d') / 2)  # symmetrise to guard numerical noise
end

"""
    _candidate_TI_axis(C::AbstractArray{T,4}; rtol=1e-8) → NTuple{3,T}

Propose a candidate TI symmetry axis for a 4th-order tensor `C` by looking
at the eigenvalues of its trace tensor `dᵢⱼ = Cᵢₖⱼₖ`.  If two eigenvalues
coincide (relative tolerance `rtol`), the third one's eigenvector is the
candidate axis — this is exact for a genuinely TI tensor.  Otherwise falls
back to `e₃ = (0, 0, 1)`.

Only meaningful for numeric element types: for symbolic `T`, the eigendecomposition is skipped and `e₃` is returned.
"""
function _candidate_TI_axis(C::AbstractArray{T, 4}; rtol = 1.0e-8) where {T <: AbstractFloat}
    d = _trace_tensor(C)
    eig = eigen(Symmetric(Matrix(d)))
    λ = eig.values
    V = eig.vectors
    span = max(abs(λ[3] - λ[1]), one(T))
    close12 = abs(λ[1] - λ[2]) < rtol * span
    close23 = abs(λ[2] - λ[3]) < rtol * span
    # Identify which eigenvalue is isolated (the TI axis direction):
    if close12 && !close23
        v = V[:, 3]        # λ₁ ≈ λ₂ → axis along v₃
    elseif !close12 && close23
        v = V[:, 1]        # λ₂ ≈ λ₃ → axis along v₁
    else
        return (zero(T), zero(T), one(T))   # ambiguous — fallback to e₃
    end
    return (T(v[1]), T(v[2]), T(v[3]))
end

# Non-floating-point fallback (symbolic, Int, etc.): return canonical e₃.
_candidate_TI_axis(::AbstractArray{T, 4}; kwargs...) where {T} = (zero(T), zero(T), one(T))

"""
    _candidate_TI_axis(M::AbstractMatrix{T}; rtol=1e-8) → NTuple{3,T}

For a 2nd-order tensor `M`, propose its unique/isolated eigenvector as the
candidate TI axis.  Falls back to `e₃` when the three eigenvalues are
indistinguishable.
"""
function _candidate_TI_axis(M::AbstractArray{T, 2}; rtol = 1.0e-8) where {T <: AbstractFloat}
    eig = eigen(Symmetric((M + M') / 2))
    λ = eig.values
    V = eig.vectors
    span = max(abs(λ[3] - λ[1]), one(T))
    close12 = abs(λ[1] - λ[2]) < rtol * span
    close23 = abs(λ[2] - λ[3]) < rtol * span
    if close12 && !close23
        v = V[:, 3]
    elseif !close12 && close23
        v = V[:, 1]
    else
        return (zero(T), zero(T), one(T))
    end
    return (T(v[1]), T(v[2]), T(v[3]))
end

_candidate_TI_axis(::AbstractArray{T, 2}; kwargs...) where {T} = (zero(T), zero(T), one(T))

"""
    _candidate_ORTHO_frame(C::AbstractArray{T,4}) → OrthonormalBasis{3,T}

Propose a candidate orthotropic material frame for `C` by taking the three
eigenvectors of its trace tensor `dᵢⱼ = Cᵢₖⱼₖ`.  Returns a
`RotatedBasis{3,T}` (or `CanonicalBasis{3,T}` if the eigenvectors coincide
with the canonical frame).

Only meaningful for numeric element types: for symbolic `T`, the canonical
frame is returned.
"""
function _candidate_ORTHO_frame(C::AbstractArray{T, 4}) where {T <: AbstractFloat}
    d = _trace_tensor(C)
    eig = eigen(Symmetric(Matrix(d)))
    R = Matrix{T}(eig.vectors)
    # Ensure a right-handed basis (det = +1, not −1)
    if det(R) < 0
        R[:, 1] = -R[:, 1]
    end
    return RotatedBasis(R)
end

_candidate_ORTHO_frame(::AbstractArray{T, 4}) where {T} = CanonicalBasis{3, T}()

"""
    _candidate_ORTHO_frame(M::AbstractMatrix{T}) → OrthonormalBasis{3,T}

For a 2nd-order tensor, the ORTHO frame candidate is directly built from its
eigenvectors.
"""
function _candidate_ORTHO_frame(M::AbstractArray{T, 2}) where {T <: AbstractFloat}
    eig = eigen(Symmetric((M + M') / 2))
    R = Matrix{T}(eig.vectors)
    if det(R) < 0
        R[:, 1] = -R[:, 1]
    end
    return RotatedBasis(R)
end

_candidate_ORTHO_frame(::AbstractArray{T, 2}) where {T} = CanonicalBasis{3, T}()

# ── proj_tens : TI, order 4, fixed axis ──────────────────────────────────────

"""
    proj_tens(::Val{:TI}, A::AbstractArray{T,4}, n) → (TensTI{4, T, 5}, d, drel)

Project a 4th-order tensor `A` (3×3×3×3) onto the transversely isotropic subspace
with fixed symmetry axis `n`. Returns a major-symmetric `TensTI{4, T, 5}`.

The projection minimises the Frobenius distance `‖B − A‖` over all TI tensors `B`
with axis `n`.

Returns a 3-tuple `(B, d, drel)`:
- `B`: the projected `TensTI{4, T, 5}`
- `d`: absolute Frobenius distance `‖B − A‖`
- `drel`: relative distance `d / ‖A‖`

# Examples
```julia
julia> n = [0., 0., 1.];

julia> C = tens_TI(10., 3., 2.5, 12., 2., n);

julia> B, d, drel = proj_tens(:TI, get_array(C), n);

julia> d < 1e-12
true

julia> arg_TI(B) == arg_TI(C)
true
```
"""
function proj_tens(::Val{:TI}, A::AbstractArray{T, 4}, n) where {T}
    nA = _frobenius(A)
    if nA ≈ zero(T)
        z = zero(T)
        return TensTI{4}(z, z, z, z, z, n), z, z
    end

    # Compute KM of A and rotate to frame where e₃ = n
    C_KM = _KM_of_array(A)
    θ, ϕ = _angles_from_n(n)
    P₆ = _KM_rotation(θ, ϕ, zero(T))
    C_rot = P₆' * C_KM * P₆

    # Project in the rotated frame
    ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆ = _project_TI_KM(C_rot)

    # Build projected tensor
    B = TensTI{4}(ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆, n)

    # Compute distances
    d = _frobenius(get_array(B) - A)
    return B, d, d / nA
end

# ── proj_tens : TI, order 2, fixed axis ──────────────────────────────────────

"""
    proj_tens(::Val{:TI}, A::AbstractArray{T,2}, n) → (TensTI{2,T,2}, d, drel)

Project a 2nd-order tensor `A` (3×3) onto the transversely isotropic subspace
with fixed symmetry axis `n`. Returns a `TensTI{2,T,2}`.

In the rotated frame where `e₃ = n`:
- `a = (M[1,1] + M[2,2]) / 2` (transverse)
- `b = M[3,3]` (axial)

# Examples
```julia
julia> n = [0., 0., 1.];

julia> A = [5. 0 0; 0 5 0; 0 0 8];

julia> B, d, drel = proj_tens(:TI, A, n);

julia> B.data
(5.0, 8.0)
```
"""
function proj_tens(::Val{:TI}, A::AbstractArray{T, 2}, n) where {T}
    nA = _frobenius(A)
    if nA ≈ zero(T)
        z = zero(T)
        return TensTI{2}(z, z, n), z, z
    end

    # Rotate to frame where e₃ = n
    θ, ϕ = _angles_from_n(n)
    R = _rot3_raw(θ, ϕ, zero(T))
    M_rot = R' * A * R

    # TI projection: average transverse, keep axial
    a = (M_rot[1, 1] + M_rot[2, 2]) / 2
    b = M_rot[3, 3]

    B = TensTI{2}(a, b, n)

    d = _frobenius(get_array(B) - A)
    return B, d, d / nA
end

# ── proj_tens : ORTHO, order 4, fixed frame ──────────────────────────────────

"""
    proj_tens(::Val{:ORTHO}, A::AbstractArray{T,4}, frame::OrthonormalBasis{3}) → (TensOrtho{T}, d, drel)

Project a 4th-order tensor `A` (3×3×3×3) onto the orthotropic subspace with
fixed material frame `frame`. Returns a `TensOrtho{T}`.

# Examples
```julia
julia> frame = CanonicalBasis{3,Float64}();

julia> t = TensOrtho(10., 8., 12., 3., 2.5, 1.5, 2., 3., 3.5, frame);

julia> B, d, drel = proj_tens(:ORTHO, get_array(t), frame);

julia> d < 1e-12
true
```
"""
function proj_tens(::Val{:ORTHO}, A::AbstractArray{T, 4}, frame::OrthonormalBasis{3}) where {T}
    nA = _frobenius(A)
    if nA ≈ zero(T)
        z = zero(T)
        B = TensOrtho(z, z, z, z, z, z, z, z, z, frame)
        return B, z, z
    end

    # Compute KM of A and rotate to material frame
    C_KM = _KM_of_array(A)
    angs = angles(Matrix{T}(vecbasis(frame, :cov)), Val(3))
    P₆ = _KM_rotation(angs.θ, angs.ϕ, angs.ψ)
    C_rot = P₆' * C_KM * P₆

    # Project
    C₁₁, C₂₂, C₃₃, C₁₂, C₁₃, C₂₃, C₄₄, C₅₅, C₆₆ = _project_ORTHO_KM(C_rot)
    B = TensOrtho(C₁₁, C₂₂, C₃₃, C₁₂, C₁₃, C₂₃, C₄₄, C₅₅, C₆₆, frame)

    d = _frobenius(get_array(B) - A)
    return B, d, d / nA
end

# ── proj_tens : ORTHO, order 2, fixed frame ──────────────────────────────────

"""
    proj_tens(::Val{:ORTHO}, A::AbstractArray{T,2}, frame::OrthonormalBasis{3}) → (Array{T,2}, d, drel)

Project a 2nd-order tensor `A` (3×3) onto the orthotropic subspace with
fixed material frame `frame`. The projection is `diag(M₁₁, M₂₂, M₃₃)` in the
material frame.

# Examples
```julia
julia> frame = CanonicalBasis{3,Float64}();

julia> A = [5. 1 2; 1 8 3; 2 3 12];

julia> B, d, drel = proj_tens(:ORTHO, A, frame);

julia> B ≈ diagm([5., 8., 12.])
true
```
"""
function proj_tens(::Val{:ORTHO}, A::AbstractArray{T, 2}, frame::OrthonormalBasis{3}) where {T}
    nA = _frobenius(A)
    if nA ≈ zero(T)
        z = zero(T)
        B = zeros(T, 3, 3)
        return B, z, z
    end

    # Rotate to material frame
    angs = angles(Matrix{T}(vecbasis(frame, :cov)), Val(3))
    R = _rot3_raw(angs.θ, angs.ϕ, angs.ψ)
    M_rot = R' * A * R

    # ORTHO projection: keep diagonal only
    B_rot = zeros(T, 3, 3)
    for i in 1:3
        B_rot[i, i] = M_rot[i, i]
    end

    # Rotate back to canonical frame
    B = R * B_rot * R'

    d = _frobenius(B - A)
    return B, d, d / nA
end

# ── Symbol dispatch with extra argument ──────────────────────────────────────

"""
    proj_tens(sym::Symbol, A::AbstractArray, arg) → (projected, d, drel)

Convenience dispatch: `proj_tens(:TI, A, n)` calls `proj_tens(Val(:TI), A, n)`.
"""
proj_tens(sym::Symbol, A::AbstractArray, arg) = proj_tens(Val(sym), A, arg)

# ── Fallbacks for optimized versions (require NLopt) ─────────────────────────
#
# `TensNDNLoptExt` supplies the real optimizers.  It must NOT define methods
# with the same signatures as the fallbacks below: overwriting a method is
# forbidden during module precompilation, so the extension would fail to
# precompile and be rebuilt (with a warning) on every session.
#
# Hence the indirection through `_proj_tens_opt`: the fallback here is the
# catch-all `Val{S}`/`AbstractArray` method, and the extension registers
# strictly *more specific* methods, which is plain dispatch — no overwrite.

"""
    _proj_tens_opt(::Val{S}, A) -> (B, d, drel)

Internal hook for the rotation-optimized projections.  The fallback throws;
`TensNDNLoptExt` adds the concrete methods when NLopt is loaded.
"""
function _proj_tens_opt(::Val{S}, A::AbstractArray) where {S}
    fixed = S === :ORTHO ? "proj_tens(:ORTHO, A, frame)" : "proj_tens(:TI, A, n)"
    return error(
        "NLopt.jl is required for rotation-optimized $S projection. " *
            "Run `using NLopt` or add NLopt to your project. " *
            "For fixed-axis/frame projection, use $fixed."
    )
end

"""
    proj_tens(::Val{:TI}, A::AbstractArray{T,4}) where {T<:AbstractFloat}

Find the best TI approximation of `A` by optimizing over all possible
symmetry axes. Requires the NLopt package: `using NLopt`.

See also [`proj_tens(::Val{:TI}, A, n)`](@ref) for fixed-axis projection.
"""
proj_tens(::Val{:TI}, A::AbstractArray{T, 4}) where {T <: AbstractFloat} =
    _proj_tens_opt(Val(:TI), A)

"""
    proj_tens(::Val{:TI}, A::AbstractArray{T,2}) where {T<:AbstractFloat}

Find the best TI approximation of a 2nd-order tensor `A` by optimizing
the symmetry axis. Requires the NLopt package: `using NLopt`.
"""
proj_tens(::Val{:TI}, A::AbstractArray{T, 2}) where {T <: AbstractFloat} =
    _proj_tens_opt(Val(:TI), A)

"""
    proj_tens(::Val{:ORTHO}, A::AbstractArray{T,4}) where {T<:AbstractFloat}

Find the best orthotropic approximation of `A` by optimizing over all
possible material frames. Requires the NLopt package: `using NLopt`.
"""
proj_tens(::Val{:ORTHO}, A::AbstractArray{T, 4}) where {T <: AbstractFloat} =
    _proj_tens_opt(Val(:ORTHO), A)

"""
    proj_tens(::Val{:ORTHO}, A::AbstractArray{T,2}) where {T<:AbstractFloat}

Find the best orthotropic approximation of a 2nd-order tensor `A` by
optimizing the material frame. Requires the NLopt package: `using NLopt`.
"""
proj_tens(::Val{:ORTHO}, A::AbstractArray{T, 2}) where {T <: AbstractFloat} =
    _proj_tens_opt(Val(:ORTHO), A)

# ──────────────────────────────────────────────────────────────────────────────
# Value-level symmetry predicates on raw arrays
# ──────────────────────────────────────────────────────────────────────────────
#
# These predicates complement the type-level predicates defined in
# `tens_isotropic.jl` and `tens_walpole.jl`: they answer "do the numerical
# components satisfy this symmetry to tolerance ε?" rather than "does this
# container impose this symmetry?".  The default cheap path uses the
# Kelvin-Mandel eigenstructure to propose a candidate axis/frame without any
# optimisation (O(1)); the expensive path with `optimize_angles=true` routes
# to the NLopt-backed `proj_tens(Val(:TI|:ORTHO), A)` fallbacks.

"""
    is_ISO(A::AbstractArray; ε=1e-6) → Bool

Return `true` when the components of `A` satisfy material isotropy up to
relative tolerance `ε` (closed-form projection: no optimisation).

The value-level predicate complements the type-level `is_ISO(::TensISO) = true`
defined in `tens_isotropic.jl`: the type version asks whether the
container *imposes* isotropy, this version asks whether the numerical
components *satisfy* it.

# Examples
```julia
julia> C = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0];

julia> is_ISO(C)
true

julia> is_ISO([1.0 2.0 0.0; 2.0 3.0 0.0; 0.0 0.0 1.0])
false
```
"""
function is_ISO(A::AbstractArray; ε = 1.0e-6)
    _, _, drel = proj_tens(Val(:ISO), A)
    return drel < ε
end

"""
    is_TI(A::AbstractArray, n; ε=1e-6) → Bool

Return `true` when `A` is transversely isotropic about the given axis `n`
up to relative tolerance `ε` (closed-form projection on the fixed axis).

# Examples
```julia
julia> n = [0.0, 0.0, 1.0];

julia> C = tens_TI(10., 3., 2.5, 12., 2., n);

julia> is_TI(get_array(C), n)
true

julia> is_TI(get_array(C), [1.0, 0.0, 0.0])
false
```
"""
function is_TI(A::AbstractArray, n; ε = 1.0e-6)
    _, _, drel = proj_tens(Val(:TI), A, n)
    return drel < ε
end

"""
    is_TI(A::AbstractArray; ε=1e-6, optimize_angles=false) → Bool

Cheap default: propose a candidate TI axis from the Kelvin-Mandel
eigendecomposition of `A` and check the residual.  With
`optimize_angles=true`, runs the NLopt-backed axis search — requires
`using NLopt`.

See also [`is_TI(A, n)`](@ref), [`_candidate_TI_axis`](@ref).
"""
function is_TI(A::AbstractArray; ε = 1.0e-6, optimize_angles::Bool = false)
    if optimize_angles
        _, _, drel = proj_tens(Val(:TI), A)
    else
        n = _candidate_TI_axis(A)
        _, _, drel = proj_tens(Val(:TI), A, n)
    end
    return drel < ε
end

"""
    is_ORTHO(A::AbstractArray, frame::OrthonormalBasis{3}; ε=1e-6) → Bool

Return `true` when `A` is orthotropic in the given material frame up to
relative tolerance `ε` (closed-form projection on the fixed frame).
"""
function is_ORTHO(A::AbstractArray, frame::OrthonormalBasis{3}; ε = 1.0e-6)
    _, _, drel = proj_tens(Val(:ORTHO), A, frame)
    return drel < ε
end

"""
    is_ORTHO(A::AbstractArray; ε=1e-6, optimize_angles=false) → Bool

Cheap default: propose a candidate orthotropic frame from the Kelvin-Mandel
eigendecomposition and check the residual.  With `optimize_angles=true`,
runs the NLopt-backed frame search — requires `using NLopt`.
"""
function is_ORTHO(A::AbstractArray; ε = 1.0e-6, optimize_angles::Bool = false)
    if optimize_angles
        _, _, drel = proj_tens(Val(:ORTHO), A)
    else
        frame = _candidate_ORTHO_frame(A)
        _, _, drel = proj_tens(Val(:ORTHO), A, frame)
    end
    return drel < ε
end

# ──────────────────────────────────────────────────────────────────────────────
# AbstractTens-dispatching convenience wrappers
# ──────────────────────────────────────────────────────────────────────────────
#
# Every `proj_tens`/`is_ISO`/`is_TI`/`is_ORTHO` method above operates on a raw
# `AbstractArray` (the caller must `get_array(t)` first). These thin wrappers
# let them be called directly on any `TensND.AbstractTens` (order 2 or 4) —
# the analogue of echoes' `.paramsym(sym)` attached to the tensor itself.
# They dispatch below the type-specific `is_ISO(::TensISO)`, `is_ISO(::TensTI)`,
# `is_ISO(::TensOrtho)`, … predicates already defined in `tens_isotropic.jl` /
# `tens_walpole.jl` (those remain the more specific, and cheaper, match for
# structured types) — in practice these wrappers fire for the fully generic
# `Tens{4,3}` / `Tens{2,3}` (no assumed symmetry), which is exactly the case
# that needs a best-fit projection rather than a type-level answer.

proj_tens(v::Val{:ISO}, t::AbstractTens) = proj_tens(v, get_array(t))
proj_tens(v::Val{:TI}, t::AbstractTens{4, dim, T}, n) where {dim, T} = proj_tens(v, get_array(t), n)
proj_tens(v::Val{:TI}, t::AbstractTens{2, dim, T}, n) where {dim, T} = proj_tens(v, get_array(t), n)
proj_tens(v::Val{:ORTHO}, t::AbstractTens{4, dim, T}, frame::OrthonormalBasis{3}) where {dim, T} =
    proj_tens(v, get_array(t), frame)
proj_tens(v::Val{:ORTHO}, t::AbstractTens{2, dim, T}, frame::OrthonormalBasis{3}) where {dim, T} =
    proj_tens(v, get_array(t), frame)
proj_tens(sym::Symbol, t::AbstractTens, args...) = proj_tens(Val(sym), t, args...)

is_ISO(t::AbstractTens; kwargs...) = is_ISO(get_array(t); kwargs...)
is_TI(t::AbstractTens, n; kwargs...) = is_TI(get_array(t), n; kwargs...)
is_TI(t::AbstractTens; kwargs...) = is_TI(get_array(t); kwargs...)
is_ORTHO(t::AbstractTens, frame; kwargs...) = is_ORTHO(get_array(t), frame; kwargs...)
is_ORTHO(t::AbstractTens; kwargs...) = is_ORTHO(get_array(t); kwargs...)

# ──────────────────────────────────────────────────────────────────────────────
# Public aliases for the per-block Kelvin-Mandel ↔ symmetry-parameter
# conversions (previously private `_project_TI_KM`/`_build_TI_KM`/
# `_project_ORTHO_KM`/`_build_ORTHO_KM`, used internally by `proj_tens`).
# Exposed publicly so that callers working directly on 6×6 Kelvin-Mandel
# blocks (e.g. MeanFieldHom's per-timestep ALV Volterra matrices) can reuse
# the same closed-form projection without going through a full `AbstractTens`.
# ──────────────────────────────────────────────────────────────────────────────

"""
    ti_params_from_KM(C::AbstractMatrix) → (ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆)

Project a 6×6 Kelvin-Mandel matrix (in a frame where `e₃` is the TI symmetry
axis) onto the major-symmetric TI subspace, returning the 5 Walpole
coefficients. See [`proj_tens`](@ref) for the tensor-level equivalent.
"""
const ti_params_from_KM = _project_TI_KM

"""
    KM_from_ti_params(ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆) → SMatrix{6,6}

Build a 6×6 Kelvin-Mandel matrix (in the frame where `e₃` is the TI axis)
from 5 major-symmetric Walpole coefficients — the reciprocal of
[`ti_params_from_KM`](@ref).
"""
const KM_from_ti_params = _build_TI_KM

"""
    ortho_params_from_KM(C::AbstractMatrix) → NTuple{9}

Extract the 9 orthotropic parameters `(C₁₁,C₂₂,C₃₃,C₁₂,C₁₃,C₂₃,C₄₄,C₅₅,C₆₆)`
from a 6×6 Kelvin-Mandel matrix in the material frame. See
[`proj_tens`](@ref) for the tensor-level equivalent.
"""
const ortho_params_from_KM = _project_ORTHO_KM

"""
    KM_from_ortho_params(C₁₁, C₂₂, C₃₃, C₁₂, C₁₃, C₂₃, C₄₄, C₅₅, C₆₆) → SMatrix{6,6}

Build a 6×6 Kelvin-Mandel matrix in the material frame from 9 orthotropic
parameters — the reciprocal of [`ortho_params_from_KM`](@ref).
"""
const KM_from_ortho_params = _build_ORTHO_KM

# ── Exports ──────────────────────────────────────────────────────────────────

export proj_tens
export ti_params_from_KM, KM_from_ti_params, ortho_params_from_KM, KM_from_ortho_params
