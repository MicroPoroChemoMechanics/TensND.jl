##############################################################################
# Tensor projection onto symmetry subspaces (TI, ORTHO)                     #
#                                                                            #
# Extends project (ISO already in tens_isotropic.jl) to:                   #
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
    mandel_rotation(θ, ϕ, ψ) → SMatrix{6,6}
    mandel_rotation(R::AbstractMatrix) → SMatrix{6,6}

6×6 Kelvin-Mandel (Bond) rotation matrix, from ZYZ Euler angles `(θ, ϕ, ψ)`
or directly from a 3×3 rotation matrix `R` whose **columns** are the new
basis vectors.

Transforms a 6×6 KM matrix as `C' = Q' · C · Q` with `Q = mandel_rotation(…)`;
`Q` is orthogonal.

Built entirely with `SMatrix` / `ntuple`, so it is compatible with
`ForwardDiff.Dual` **and** with non-`isbits` element types such as
`SymPyPythonCall.Sym` or `Symbolics.Num`, for which an `MMatrix{6,6,T}(undef)` cannot
even be constructed.

# Examples
```julia
julia> Q = mandel_rotation(0.3, 0.5, 0.0);

julia> size(Q)
(6, 6)
```
"""
mandel_rotation(θ, ϕ, ψ) = mandel_rotation(_rot3_raw(θ, ϕ, ψ))

export mandel_rotation

function mandel_rotation(R::AbstractMatrix{T}) where {T}
    sq2 = sqrt(2 * one(T))
    return SMatrix{6, 6, T}(
        ntuple(Val(36)) do n
            # column-major: element n is (row I, column J)
            I = mod1(n, 6)
            J = (n - 1) ÷ 6 + 1
            i, j = _KM_COUPLES[I]
            k, l = _KM_COUPLES[J]
            @inbounds if I ≤ 3 && J ≤ 3
                R[i, k]^2                                   # normal-normal
            elseif I ≤ 3
                sq2 * R[i, k] * R[i, l]                     # normal-shear
            elseif J ≤ 3
                sq2 * R[i, k] * R[j, k]                     # shear-normal
            else
                R[i, k] * R[j, l] + R[i, l] * R[j, k]       # shear-shear
            end
        end
    )
end

"""
    _KM_of_array(A::AbstractArray{T,4}) → Matrix{T}

Compute the 6×6 Kelvin-Mandel matrix of a 3×3×3×3 array.
"""
function _KM_of_array(A::AbstractArray{T, 4}) where {T}
    return _km66(tensor_or_array(A), A)
end

# The 6×6 matrix whether or not `tensor_or_array` recognized the minor
# symmetries: past its tolerance it returns a full `Tensor`, whose `tomandel` is
# 9×9 and broke every 6×6 product of the projections. The minor-symmetric part
# is what those projections act on.
_km66(t::Tensors.SymmetricTensor{4}, A) = tomandel(t)
_km66(t, A) = _mandel_minor(A)

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

# ── CUBIC projection helpers ─────────────────────────────────────────────────

"""
    _project_CUBIC_KM(C::AbstractMatrix) → NTuple{3}

Extract the three cubic coefficients `(α, β, γ)` on `(𝕁, 𝔼, 𝕋)` from a 6×6
Kelvin-Mandel matrix expressed **in the cube frame**.

The three projectors being mutually orthogonal and idempotent with traces
`1, 2, 3`, the orthogonal projection is read off by three traces:

```math
\\alpha = \\operatorname{tr}(\\mathbb C\\,\\mathbb J), \\qquad
\\beta  = \\tfrac12\\operatorname{tr}(\\mathbb C\\,\\mathbb E), \\qquad
\\gamma = \\tfrac13\\operatorname{tr}(\\mathbb C\\,\\mathbb T),
```

which is what makes this an averaging of the entries rather than a fit: `α` is
the mean of the whole upper 3×3 block, `β` the mean of its diagonal minus its
off-diagonal, `γ` the mean of the shear diagonal. A non-major-symmetric input
is handled by the same traces, which symmetrize it implicitly.

# Examples
```julia
julia> C = Float64[10 4 4 0 0 0; 4 10 4 0 0 0; 4 4 10 0 0 0;
                   0 0 0 4 0 0; 0 0 0 0 4 0; 0 0 0 0 0 4];

julia> _project_CUBIC_KM(C)             # (C₁₁+2C₁₂, C₁₁-C₁₂, 2C₄₄)
(18.0, 6.0, 4.0)
```
"""
function _project_CUBIC_KM(C)
    # tr(C 𝕁) with 𝕁 = v vᵀ/3: the sum of the whole upper block over 3.
    α = (
        C[1, 1] + C[1, 2] + C[1, 3] + C[2, 1] + C[2, 2] +
            C[2, 3] + C[3, 1] + C[3, 2] + C[3, 3]
    ) / 3
    # tr(C 𝔼)/2 with 𝔼 = upper block minus its own 𝕁 part.
    β = (C[1, 1] + C[2, 2] + C[3, 3] - α) / 2
    γ = (C[4, 4] + C[5, 5] + C[6, 6]) / 3
    return (α, β, γ)
end

"""
    _build_CUBIC_KM(α, β, γ) → SMatrix{6,6}

Build the 6×6 Kelvin-Mandel matrix in the cube frame from the three projector
coefficients — the reciprocal of [`_project_CUBIC_KM`](@ref).
"""
function _build_CUBIC_KM(α, β, γ)
    T = promote_type(typeof(α), typeof(β), typeof(γ))
    C₁₁ = α / 3 + 2β / 3
    C₁₂ = α / 3 - β / 3
    z = zero(T)
    return SMatrix{6, 6, T}(
        C₁₁, C₁₂, C₁₂, z, z, z,
        C₁₂, C₁₁, C₁₂, z, z, z,
        C₁₂, C₁₂, C₁₁, z, z, z,
        z, z, z, γ, z, z,
        z, z, z, z, γ, z,
        z, z, z, z, z, γ
    )
end

# ── Why cubic has no eigenstructure candidate ────────────────────────────────
#
#  `_candidate_TI_axis` and `_candidate_ORTHO_frame` both work by diagonalizing
#  a second-order contraction of the tensor.  For a cubic tensor that cannot
#  work, and the reason is not a numerical one: **every second-order
#  contraction of a cubic tensor is isotropic**.  Measured on
#  `cubic_from_constants(10, 4, 2, frame)`, in any frame,
#
#      C_iikl = 18 δ_kl ,      C_ikil = 14 δ_kl ,
#
#  so no eigenframe carries any information about where the cube points; the
#  cube axes are a genuinely *fourth-order* feature.  Handing
#  `_candidate_ORTHO_frame` a cubic tensor returns an arbitrary frame — measured
#  at a relative residual of 0.099 where the true frame gives 3.5e-16.
#
#  So there is no cheap exact candidate to start from, and the free-orientation
#  cubic projection optimizes from a fixed angular grid alone.  It is the one
#  place where this class is *less* well served than its neighbors, and it is a
#  property of the symmetry rather than of the implementation.

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
# running any optimization.

"""
    _trace_tensor(C::AbstractArray{T,4}) → SMatrix{3,3,T,9}

Return the 2nd-order "trace" tensor `dᵢⱼ = ∑ₖ C[i,k,j,k]` — its eigenstructure
exposes the principal directions of `C` cheaply (no optimization required).

Used to propose candidate TI axes and candidate orthotropic frames when the
user asks for symmetry detection without angle optimization.
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
    return SMatrix{3, 3, T}((d + d') / 2)  # symmetrize to guard numerical noise
end

"""
    _sym_eigen(d::AbstractMatrix; rtol = 1e-8) -> (λ, V)

Eigenvalues, ascending, and orthonormal eigenvectors of the real symmetric
matrix `d`, as `eigen(Symmetric(d))` returns them.

`eigen` goes through LAPACK and stops at `AbstractFloat`. For
`ForwardDiff.Dual` entries the values are those of the value matrix, and each
derivative comes from first-order perturbation theory, one nesting level at a
time:

    λ̇ᵢ = vᵢ·ḋ·vᵢ,    v̇ᵢ = Σⱼ≠ᵢ (vⱼ·ḋ·vᵢ) / (λᵢ − λⱼ) vⱼ.

Two eigenvalues closer than `rtol · max(λₙ − λ₁, 1)` — the closeness by which
`_candidate_TI_axis` recognizes a double eigenvalue — are treated as one:
within their plane no direction is preferred, so neither contributes to the
derivative of the other's vector.
"""
_sym_eigen(d::AbstractMatrix{<:AbstractFloat}; rtol = 1.0e-8) =
    (e = eigen(Symmetric(Matrix(d))); (e.values, e.vectors))

function _sym_eigen(d::AbstractMatrix{D}; rtol = 1.0e-8) where {D <: ForwardDiff.Dual}
    λ₀, V₀ = _sym_eigen(ForwardDiff.value.(d); rtol)
    n = length(λ₀)
    span = max(λ₀[n] - λ₀[1], one(eltype(λ₀)))
    # Each derivative of `d` in the eigenbasis: entry (j, i) is vⱼ·ḋ·vᵢ.
    ḋ = ntuple(k -> V₀' * ForwardDiff.partials.(d, k) * V₀, Val(ForwardDiff.npartials(D)))
    function v̇(m, a, i)
        s = zero(eltype(m))
        for j in 1:n
            (j == i || abs(λ₀[i] - λ₀[j]) < rtol * span) && continue
            s += m[j, i] / (λ₀[i] - λ₀[j]) * V₀[a, j]
        end
        return s
    end
    λ = [D(λ₀[i], ForwardDiff.Partials(map(m -> m[i, i], ḋ))) for i in 1:n]
    V = [D(V₀[a, i], ForwardDiff.Partials(map(m -> v̇(m, a, i), ḋ))) for a in 1:n, i in 1:n]
    return λ, V
end

# The element types whose candidates come from an eigendecomposition.
const _EigenType = Union{AbstractFloat, ForwardDiff.Dual}

"""
    _candidate_TI_axis(C::AbstractArray{T,4}; rtol=1e-8) → NTuple{3,T}

Propose a candidate TI symmetry axis for a 4th-order tensor `C` by looking
at the eigenvalues of its trace tensor `dᵢⱼ = Cᵢₖⱼₖ`.  If two eigenvalues
coincide (relative tolerance `rtol`), the third one's eigenvector is the
candidate axis — this is exact for a genuinely TI tensor.  Otherwise falls
back to `e₃ = (0, 0, 1)`.

Only meaningful for numeric element types: for symbolic `T`, the eigendecomposition is skipped and `e₃` is returned.
"""
function _candidate_TI_axis(C::AbstractArray{T, 4}; rtol = 1.0e-8) where {T <: _EigenType}
    λ, V = _sym_eigen(_trace_tensor(C); rtol)
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

# Fallback for the other element types (symbolic, integer…): the canonical e₃.
_candidate_TI_axis(::AbstractArray{T, 4}; kwargs...) where {T} = (zero(T), zero(T), one(T))

"""
    _candidate_TI_axis(M::AbstractMatrix{T}; rtol=1e-8) → NTuple{3,T}

For a 2nd-order tensor `M`, propose its unique/isolated eigenvector as the
candidate TI axis.  Falls back to `e₃` when the three eigenvalues are
indistinguishable.
"""
function _candidate_TI_axis(M::AbstractArray{T, 2}; rtol = 1.0e-8) where {T <: _EigenType}
    λ, V = _sym_eigen((M + M') / 2; rtol)
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
function _candidate_ORTHO_frame(C::AbstractArray{T, 4}) where {T <: _EigenType}
    R = Matrix{T}(last(_sym_eigen(_trace_tensor(C))))
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
function _candidate_ORTHO_frame(M::AbstractArray{T, 2}) where {T <: _EigenType}
    R = Matrix{T}(last(_sym_eigen((M + M') / 2)))
    if det(R) < 0
        R[:, 1] = -R[:, 1]
    end
    return RotatedBasis(R)
end

_candidate_ORTHO_frame(::AbstractArray{T, 2}) where {T} = CanonicalBasis{3, T}()

# ── project : TI, order 4, fixed axis ──────────────────────────────────────

"""
    project(A::AbstractArray{T,4}, ::Val{:TI}, n) → (TensTI{4, T, 5}, d, drel)

Project a 4th-order tensor `A` (3×3×3×3) onto the transversely isotropic subspace
with fixed symmetry axis `n`. Returns a major-symmetric `TensTI{4, T, 5}`.

The projection minimizes the Frobenius distance `‖B − A‖` over all TI tensors `B`
with axis `n`.

Returns a 3-tuple `(B, d, drel)`:
- `B`: the projected `TensTI{4, T, 5}`
- `d`: absolute Frobenius distance `‖B − A‖`
- `drel`: relative distance `d / ‖A‖`

# Examples
```julia
julia> n = [0., 0., 1.];

julia> C = ti_from_constants(10., 3., 2.5, 12., 2., n);

julia> B, d, drel = project(get_array(C), :TI, n);

julia> d < 1e-12
true

julia> ti_constants(B) == ti_constants(C)
true
```
"""
function project(A::AbstractArray{T, 4}, ::Val{:TI}, n) where {T}
    nA = _frobenius(A)
    if nA ≈ zero(T)
        z = zero(T)
        return TensTI{4}(z, z, z, z, z, n), z, z
    end

    # The Kelvin-Mandel matrix of A in a frame whose third axis is n, by the 6×6
    # congruence of that frame: no Euler angles, whose inverse trigonometry has an
    # infinite derivative at n = e₃ (a NaN through ForwardDiff), and no `'`, which
    # conjugates a symbol not declared real. The projection does not depend on
    # the in-plane axes the frame picks.
    C_rot = _km_in_frame_of(_KM_of_array(A), first(_axis_frame(n)))

    # Project in the rotated frame
    ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆ = _project_TI_KM(C_rot)

    # Build projected tensor
    B = TensTI{4}(ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆, n)

    # Compute distances
    d = _frobenius(get_array(B) - A)
    return B, d, d / nA
end

# ── project : TI, order 2, fixed axis ──────────────────────────────────────

"""
    project(A::AbstractArray{T,2}, ::Val{:TI}, n) → (TensTI{2,T,2}, d, drel)

Project a 2nd-order tensor `A` (3×3) onto the transversely isotropic subspace
with fixed symmetry axis `n`. Returns a `TensTI{2,T,2}`.

In the rotated frame where `e₃ = n`:
- `a = (M[1,1] + M[2,2]) / 2` (transverse)
- `b = M[3,3]` (axial)

# Examples
```julia
julia> n = [0., 0., 1.];

julia> A = [5. 0 0; 0 5 0; 0 0 8];

julia> B, d, drel = project(A, :TI, n);

julia> B.data
(5.0, 8.0)
```
"""
function project(A::AbstractArray{T, 2}, ::Val{:TI}, n) where {T}
    nA = _frobenius(A)
    if nA ≈ zero(T)
        z = zero(T)
        return TensTI{2}(z, z, n), z, z
    end

    # In a frame whose third axis is n (no Euler angles: their derivative is
    # infinite at n = e₃; `transpose`: `'` conjugates a symbol), the transverse
    # average and the axial component.
    R = first(_axis_frame(n))
    a, b, _ = ti2_params_from_matrix(transpose(R) * A * R)

    B = TensTI{2}(a, b, n)

    d = _frobenius(get_array(B) - A)
    return B, d, d / nA
end

# ── project : ORTHO, order 4, fixed frame ──────────────────────────────────

"""
    project(A::AbstractArray{T,4}, ::Val{:ORTHO}, frame::OrthonormalBasis{3}) → (TensOrtho{4,T,9}, d, drel)

Project a 4th-order tensor `A` (3×3×3×3) onto the orthotropic subspace with
fixed material frame `frame`. Returns a `TensOrtho{4,T,9}`.

# Examples
```julia
julia> frame = CanonicalBasis{3,Float64}();

julia> t = TensOrtho(10., 8., 12., 3., 2.5, 1.5, 2., 3., 3.5, frame);

julia> B, d, drel = project(get_array(t), :ORTHO, frame);

julia> d < 1e-12
true
```
"""
function project(A::AbstractArray{T, 4}, ::Val{:ORTHO}, frame::OrthonormalBasis{3}) where {T}
    nA = _frobenius(A)
    if nA ≈ zero(T)
        z = zero(T)
        B = TensOrtho(z, z, z, z, z, z, z, z, z, frame)
        return B, z, z
    end

    # The Kelvin-Mandel matrix of A in the material frame, by the 6×6 congruence
    # of the frame matrix itself. Recovering Euler angles from the frame went
    # through inverse trigonometry with an infinite derivative at the canonical
    # frame (gimbal lock), a NaN under ForwardDiff for a frame that depends on a
    # parameter, and `'` conjugated a symbol not declared real.
    C_rot = _km_in_frame_of(_KM_of_array(A), Matrix(vecbasis(frame, :cov)))

    # Project
    C₁₁, C₂₂, C₃₃, C₁₂, C₁₃, C₂₃, C₄₄, C₅₅, C₆₆ = _project_ORTHO_KM(C_rot)
    B = TensOrtho(C₁₁, C₂₂, C₃₃, C₁₂, C₁₃, C₂₃, C₄₄, C₅₅, C₆₆, frame)

    d = _frobenius(get_array(B) - A)
    return B, d, d / nA
end

# ── project : CUBIC, order 4, fixed cube frame ─────────────────────────────

"""
    project(A::AbstractArray{T,4}, ::Val{:CUBIC}, frame::OrthonormalBasis{3}) → (TensCubic{T}, d, drel)

Project a 4th-order tensor onto the cubic subspace with cube axes `frame`.

**Read `drel` rather than trusting the projection.** Where the morphology and
the medium leave the octahedral group invariant, the answer belongs to the
class by group theory, so its distance to it is discretization error and
nothing else — an error estimate that costs nothing and involves no reference
solution. Where they do not, `drel` says how much was discarded, and the
projection is a fit like any other.

# Examples
```julia
julia> frame = CanonicalBasis{3,Float64}();

julia> t = cubic_from_constants(10.0, 4.0, 2.0, frame);

julia> B, d, drel = project(get_array(t), :CUBIC, frame);

julia> d < 1e-12
true
```

See also [`best_fit_cubic`](@ref), [`cubic_anisotropy`](@ref).
"""
function project(A::AbstractArray{T, 4}, ::Val{:CUBIC}, frame::OrthonormalBasis{3}) where {T}
    nA = _frobenius(A)
    if nA ≈ zero(T)
        z = zero(T)
        return TensCubic{T}((z, z, z), frame), z, z
    end
    # The material frame by its 6×6 congruence, as for the ORTHO projection.
    C_rot = _km_in_frame_of(_KM_of_array(A), Matrix(vecbasis(frame, :cov)))
    α, β, γ = _project_CUBIC_KM(C_rot)
    B = TensCubic(α, β, γ, frame)
    d = _frobenius(get_array(B) - A)
    return B, d, d / nA
end

"""
    project(A::AbstractArray{T,2}, ::Val{:CUBIC}, frame) → (TensISO{2,3,T}, d, drel)

**At order two the cubic class is the isotropic class**, the octahedral group
leaving no second-order tensor invariant but a multiple of the identity. The
frame is accepted and ignored, and the isotropic projection is returned.

That is not a shortcut: it is why a cube-symmetric pore has a single scalar
resistivity contribution while its compliance contribution needs three
constants, and why a conduction computation on such a shape carries no
anisotropy signal whatever the shape does in elasticity.
"""
project(A::AbstractArray{T, 2}, ::Val{:CUBIC}, ::OrthonormalBasis{3}) where {T} =
    project(A, Val(:ISO))

# ── project : ORTHO, order 2, fixed frame ──────────────────────────────────

"""
    project(A::AbstractArray{T,2}, ::Val{:ORTHO}, frame::OrthonormalBasis{3}) → (TensOrtho{2,T,3}, d, drel)

Project a 2nd-order tensor `A` (3×3) onto the orthotropic subspace with
fixed material frame `frame`. The projection is `diag(M₁₁, M₂₂, M₃₃)` in the
material frame, returned as a `TensOrtho{2}` on that frame (as the TI projection
returns a `TensTI{2}`).

# Examples
```julia
julia> frame = CanonicalBasis{3,Float64}();

julia> A = [5. 1 2; 1 8 3; 2 3 12];

julia> B, d, drel = project(A, :ORTHO, frame);

julia> get_array(B) ≈ diagm([5., 8., 12.])
true
```
"""
function project(A::AbstractArray{T, 2}, ::Val{:ORTHO}, frame::OrthonormalBasis{3}) where {T}
    nA = _frobenius(A)
    if nA ≈ zero(T)
        z = zero(T)
        return TensOrtho{2}(z, z, z, frame), z, z
    end

    # The frame's own matrix (columns: its vectors in canonical components):
    # exact, in the frame's element type, and with no round trip through Euler
    # angles, whose inverse trigonometry has an infinite derivative at the
    # gimbal-lock point (the canonical frame itself).
    R = Matrix(vecbasis(frame, :cov))
    M_rot = transpose(R) * A * R
    # ORTHO projection: the diagonal in the material frame.
    B = TensOrtho{2}(M_rot[1, 1], M_rot[2, 2], M_rot[3, 3], frame)

    d = _frobenius(get_array(B) - A)
    return B, d, d / nA
end

# ── Symbol dispatch with extra argument ──────────────────────────────────────

"""
    project(A::AbstractArray, sym::Symbol, arg) → (projected, d, drel)

Convenience dispatch: `project(A, :TI, n)` calls `project(A, Val(:TI), n)`.
"""
project(A::AbstractArray, sym::Symbol, arg) = project(A, Val(sym), arg)

# ── Projections with a free orientation ──────────────────────────────────────
#
# The orientation is searched by `_proj_tens_opt` (`orientation_search.jl`),
# for floating-point arrays.

# The orientation found by the optimizer carries no derivative: it would take
# the implicit function theorem on its optimality condition. Refused, rather
# than a projection whose derivative ignores how the orientation moves.
for (S, N) in ((:TI, 4), (:TI, 2), (:ORTHO, 4), (:ORTHO, 2), (:CUBIC, 4))
    msg = "the rotation-optimized :$S projection does not differentiate the " *
        "orientation it finds; project on a given axis or frame, or call " *
        "best_sym_tens with optimize_angles = false, whose candidates carry derivatives."
    @eval project(::AbstractArray{<:ForwardDiff.Dual, $N}, ::Val{$(QuoteNode(S))}) =
        throw(ArgumentError($msg))
end

"""
    project(A::AbstractArray{T,4}, ::Val{:TI}) where {T<:AbstractFloat}

Find the best TI approximation of `A` over all symmetry axes. The axis is
searched — the objective screened on 200 axes, Newton's method on the angles
from the best ones, from the eigenstructure candidate and from `e₃` — with no
random step, so a call is repeatable to the bit.

See also [`project(A, ::Val{:TI}, n)`](@ref) for fixed-axis projection.
"""
project(A::AbstractArray{T, 4}, ::Val{:TI}) where {T <: AbstractFloat} =
    _proj_tens_opt(Val(:TI), A)

"""
    project(A::AbstractArray{T,2}, ::Val{:TI}) where {T<:AbstractFloat}

Find the best TI approximation of a 2nd-order tensor `A` over all symmetry
axes. No search is needed: the axis is an eigenvector of the smallest or the
largest eigenvalue of the symmetric part of `A`, whichever gives the closer
tensor.
"""
project(A::AbstractArray{T, 2}, ::Val{:TI}) where {T <: AbstractFloat} =
    _proj_tens_opt(Val(:TI), A)

"""
    project(A::AbstractArray{T,4}, ::Val{:ORTHO}) where {T<:AbstractFloat}

Find the best orthotropic approximation of `A` over all material frames. The
frame is searched as the TI axis is (400 frames screened, Newton's method on
the best ones, on the eigenstructure candidate and on the canonical frame).
"""
project(A::AbstractArray{T, 4}, ::Val{:ORTHO}) where {T <: AbstractFloat} =
    _proj_tens_opt(Val(:ORTHO), A)

"""
    project(A::AbstractArray{T,4}, ::Val{:CUBIC}) where {T<:AbstractFloat}

Find the best cubic approximation of `A` over all **cube orientations**, searched
as an orthotropic frame is, but with no eigenstructure candidate, there being
none for a cubic tensor.

The octahedral group is discrete, but the *orientation of the cube* is not: it
is an ordinary rotation, so the objective is smooth in the Euler angles exactly
as it is for a TI axis or an orthotropic frame. The only trace the group leaves
is that the optimum is 24-fold degenerate — the optimizer returns one of the 24
equivalent frames, and they all describe the same tensor with the same
coefficients (see `TensCubic`).

See also [`project(A, ::Val{:CUBIC}, frame)`](@ref) for fixed-frame
projection.
"""
project(A::AbstractArray{T, 4}, ::Val{:CUBIC}) where {T <: AbstractFloat} =
    _proj_tens_opt(Val(:CUBIC), A)

"""
    project(A::AbstractArray{T,2}, ::Val{:CUBIC})

At order two the cubic class is the isotropic class, which has no orientation
to optimize, so this is the isotropic projection, for any element type the
isotropic one takes.
"""
project(A::AbstractArray{T, 2}, ::Val{:CUBIC}) where {T} = project(A, Val(:ISO))

"""
    project(A::AbstractArray{T,2}, ::Val{:ORTHO}) where {T<:AbstractFloat}

Find the best orthotropic approximation of a 2nd-order tensor `A` over all
material frames. No search is needed: it is the eigenframe of the symmetric
part of `A`, in which the projection keeps all of that part.
"""
project(A::AbstractArray{T, 2}, ::Val{:ORTHO}) where {T <: AbstractFloat} =
    _proj_tens_opt(Val(:ORTHO), A)

# ──────────────────────────────────────────────────────────────────────────────
# Value-level symmetry predicates on raw arrays
# ──────────────────────────────────────────────────────────────────────────────
#
# These predicates complement the type-level predicates defined in
# `tens_isotropic.jl` and `tens_anisotropic.jl`: they answer "do the numerical
# components satisfy this symmetry to tolerance ε?" rather than "does this
# container impose this symmetry?".  The default cheap path uses the
# Kelvin-Mandel eigenstructure to propose a candidate axis/frame without any
# optimization (O(1)); the path with `optimize_angles=true` searches the
# orientation with `project(A, Val(:TI|:ORTHO))`.

"""
    isisotropic(A::AbstractArray; ε=1e-6) → Bool

Return `true` when the components of `A` satisfy material isotropy up to
relative tolerance `ε` (closed-form projection: no optimization).

The value-level predicate complements the type-level `isisotropic(::TensISO) = true`
defined in `tens_isotropic.jl`: the type version asks whether the
container *imposes* isotropy, this version asks whether the numerical
components *satisfy* it.

# Examples
```julia
julia> C = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0];

julia> isisotropic(C)
true

julia> isisotropic([1.0 2.0 0.0; 2.0 3.0 0.0; 0.0 0.0 1.0])
false
```
"""
function isisotropic(A::AbstractArray; ε = 1.0e-6)
    _, _, drel = project(A, Val(:ISO))
    return drel < ε
end

"""
    istransverselyisotropic(A::AbstractArray, n; ε=1e-6) → Bool

Return `true` when `A` is transversely isotropic about the given axis `n`
up to relative tolerance `ε` (closed-form projection on the fixed axis).

# Examples
```julia
julia> n = [0.0, 0.0, 1.0];

julia> C = ti_from_constants(10., 3., 2.5, 12., 2., n);

julia> istransverselyisotropic(get_array(C), n)
true

julia> istransverselyisotropic(get_array(C), [1.0, 0.0, 0.0])
false
```
"""
function istransverselyisotropic(A::AbstractArray, n; ε = 1.0e-6)
    _, _, drel = project(A, Val(:TI), n)
    return drel < ε
end

"""
    istransverselyisotropic(A::AbstractArray; ε=1e-6, optimize_angles=false) → Bool

Cheap default: propose a candidate TI axis from the Kelvin-Mandel
eigendecomposition of `A` and check the residual.  With
`optimize_angles=true`, searches the axis (`project(A, :TI)`).

See also [`istransverselyisotropic(A, n)`](@ref), `_candidate_TI_axis`.
"""
function istransverselyisotropic(A::AbstractArray; ε = 1.0e-6, optimize_angles::Bool = false)
    if optimize_angles
        _, _, drel = project(A, Val(:TI))
    else
        n = _candidate_TI_axis(A)
        _, _, drel = project(A, Val(:TI), n)
    end
    return drel < ε
end

"""
    isorthotropic(A::AbstractArray, frame::OrthonormalBasis{3}; ε=1e-6) → Bool

Return `true` when `A` is orthotropic in the given material frame up to
relative tolerance `ε` (closed-form projection on the fixed frame).
"""
function isorthotropic(A::AbstractArray, frame::OrthonormalBasis{3}; ε = 1.0e-6)
    _, _, drel = project(A, Val(:ORTHO), frame)
    return drel < ε
end

"""
    isorthotropic(A::AbstractArray; ε=1e-6, optimize_angles=false) → Bool

Cheap default: propose a candidate orthotropic frame from the Kelvin-Mandel
eigendecomposition and check the residual.  With `optimize_angles=true`,
searches the frame (`project(A, :ORTHO)`).
"""
function isorthotropic(A::AbstractArray; ε = 1.0e-6, optimize_angles::Bool = false)
    if optimize_angles
        _, _, drel = project(A, Val(:ORTHO))
    else
        frame = _candidate_ORTHO_frame(A)
        _, _, drel = project(A, Val(:ORTHO), frame)
    end
    return drel < ε
end

# ──────────────────────────────────────────────────────────────────────────────
# AbstractTens-dispatching convenience wrappers
# ──────────────────────────────────────────────────────────────────────────────
#
# The kernels above work on an array and read an axis or a frame in the
# coordinates of that array. A tensor is projected in the orthonormal frame its
# components are written in: an axis or a frame given in canonical coordinates
# is brought into that frame (a 3-vector or a 3×3 product), never the 81
# components out of it, and the axis or frame of the result is taken back to
# canonical coordinates. A tensor written in the canonical basis pays nothing
# for this, and the structured types keep their cheaper type-level predicates
# (`isisotropic(::TensISO)`, …), which are more specific than these wrappers.
#
# `get_array(t)` alone, read as canonical, was wrong for every tensor written in
# a rotated basis: the projection about a canonical axis was computed about that
# axis expressed in the wrong frame.

"""
    _orthonormal_components(t) -> (A, Q)

Components `A` of `t` in an orthonormal frame, and the rotation `Q` from that
frame to the canonical one (its columns are the frame's vectors in canonical
coordinates), or `nothing` when the frame is canonical.

A tensor written in an orthonormal basis is read as it is: covariant and
contravariant components coincide there, so its stored array is the answer.
Only an orthogonal or a general basis changes the components, through the
metric and the variances of each index (`components`), to its orthonormal
parent or to the canonical basis.
"""
_orthonormal_components(t::AbstractTens) = _orthonormal_components(t, get_basis(t))
_orthonormal_components(t::AbstractTens, ::CanonicalBasis) = (get_array(t), nothing)
_orthonormal_components(t::AbstractTens, ℬ::OrthonormalBasis) = (get_array(t), _rotation_to_canon(ℬ))
function _orthonormal_components(t::AbstractTens, ℬ::AbstractBasis)
    ℬo = relevant_OrthonormalBasis(ℬ)
    return components(t, ℬo), _rotation_to_canon(ℬo)
end

_rotation_to_canon(::CanonicalBasis) = nothing
_rotation_to_canon(ℬ::OrthonormalBasis) = Matrix(vecbasis(ℬ))

# The frame vectors as the columns of a static matrix: read in hot loops
# (`symmetry_basis`), so never a heap copy.
_frame_matrix(::CanonicalBasis{3, T}) where {T} = SMatrix{3, 3, T}(I)
_frame_matrix(ℬ::OrthonormalBasis{3}) = SMatrix{3, 3}(vecbasis(ℬ, Val(:cov)))

# A canonical axis or frame, in the frame of `Q`.
_to_frame(::Nothing, x) = x
_to_frame(Q::AbstractMatrix, n) = transpose(Q) * collect(_extract_vec(n))
_to_frame(Q::AbstractMatrix, f::OrthonormalBasis{3}) =
    RotatedBasis(transpose(Q) * _frame_matrix(f))

# A projection computed in the frame of `Q`, with its axis or frame in canonical
# coordinates. The coefficients do not depend on the frame.
_from_frame(::Nothing, P) = P
_from_frame(::AbstractMatrix, P::TensISO) = P
function _from_frame(Q::AbstractMatrix, P::TensTI{order, T, N}) where {order, T, N}
    T2 = promote_type(T, eltype(Q))
    n = Q * collect(axis(P))
    return TensTI{order, T2, N}(map(T2, get_data(P)), (T2(n[1]), T2(n[2]), T2(n[3])))
end
_from_frame(Q::AbstractMatrix, P::TensOrtho) =
    _rebuild_frame(P, RotatedBasis(Q * _frame_matrix(frame(P))))
_from_frame(Q::AbstractMatrix, P::TensCubic) =
    TensCubic{eltype(P)}(get_data(P), RotatedBasis(Q * _frame_matrix(frame(P))))

_from_frame(Q::AbstractMatrix, (P, d, drel)::Tuple) = (_from_frame(Q, P), d, drel)

function _proj_in_frame(v::Val, t::AbstractTens, args...)
    A, Q = _orthonormal_components(t)
    return _from_frame(Q, project(A, v, map(x -> _to_frame(Q, x), args)...))
end

project(t::AbstractTens, v::Val{:ISO}) = project(first(_orthonormal_components(t)), v)
project(t::AbstractTens{4, dim, T}, v::Val{:TI}, n) where {dim, T} = _proj_in_frame(v, t, n)
project(t::AbstractTens{2, dim, T}, v::Val{:TI}, n) where {dim, T} = _proj_in_frame(v, t, n)
project(t::AbstractTens{4, dim, T}, v::Val{:ORTHO}, frame::OrthonormalBasis{3}) where {dim, T} =
    _proj_in_frame(v, t, frame)
project(t::AbstractTens{2, dim, T}, v::Val{:ORTHO}, frame::OrthonormalBasis{3}) where {dim, T} =
    _proj_in_frame(v, t, frame)
project(t::AbstractTens{4, dim, T}, v::Val{:CUBIC}, frame::OrthonormalBasis{3}) where {dim, T} =
    _proj_in_frame(v, t, frame)
project(t::AbstractTens{2, dim, T}, v::Val{:CUBIC}, frame::OrthonormalBasis{3}) where {dim, T} =
    _proj_in_frame(v, t, frame)
# The searches, axis or frame free. Without these, `project(t, :TI)` reached
# the array kernel with `t` itself, read in its own basis as if canonical. The
# element-type bound matches the kernels', so each wrapper is strictly more
# specific than the kernel it unwraps for.
project(t::AbstractTens{4, dim, T}, v::Val{:TI}) where {dim, T <: AbstractFloat} = _proj_in_frame(v, t)
project(t::AbstractTens{2, dim, T}, v::Val{:TI}) where {dim, T <: AbstractFloat} = _proj_in_frame(v, t)
project(t::AbstractTens{4, dim, T}, v::Val{:ORTHO}) where {dim, T <: AbstractFloat} = _proj_in_frame(v, t)
project(t::AbstractTens{2, dim, T}, v::Val{:ORTHO}) where {dim, T <: AbstractFloat} = _proj_in_frame(v, t)
project(t::AbstractTens{4, dim, T}, v::Val{:CUBIC}) where {dim, T <: AbstractFloat} = _proj_in_frame(v, t)
project(t::AbstractTens{2}, ::Val{:CUBIC}) = project(t, Val(:ISO))
project(t::AbstractTens, sym::Symbol, args...) = project(t, Val(sym), args...)

"""
    iscubic(A::AbstractArray, frame::OrthonormalBasis{3}; ε = 1.0e-6) -> Bool

Whether `A` is cubic about `frame` to within `ε`, measured as the relative
Frobenius distance to the class.

The **value-level** companion of the type-level `iscubic(::TensCubic)`, exactly
as for [`isisotropic`](@ref), [`istransverselyisotropic`](@ref) and [`isorthotropic`](@ref): the first asks
what the components say, the second what the container guarantees.
"""
function iscubic(A::AbstractArray, frame::OrthonormalBasis{3}; ε = 1.0e-6)
    _, _, drel = project(A, Val(:CUBIC), frame)
    return drel ≤ ε
end

"""
    iscubic(A::AbstractArray; ε = 1.0e-6) -> Bool

Whether `A` is cubic about **some** cube, the cube orientation being searched
(`project(A, :CUBIC)`). Unlike [`istransverselyisotropic`](@ref) and [`isorthotropic`](@ref) it has
no cheap alternative: every second-order contraction of a cubic tensor is
isotropic, so no eigenframe reveals where the cube points and there is nothing
to start from but a search.

Give the frame — `iscubic(A, frame)` — whenever it is known; that path is exact
and costs one projection.
"""
function iscubic(A::AbstractArray; ε = 1.0e-6)
    _, _, drel = project(A, Val(:CUBIC))
    return drel ≤ ε
end

# The predicates: a distance, which no rotation changes, so only the axis or the
# frame is brought into the tensor's frame.
function _in_frame(t::AbstractTens, f::F) where {F}
    A, Q = _orthonormal_components(t)
    return f(A, Q)
end
isisotropic(t::AbstractTens; kwargs...) = isisotropic(first(_orthonormal_components(t)); kwargs...)
istransverselyisotropic(t::AbstractTens, n; kwargs...) = _in_frame(t, (A, Q) -> istransverselyisotropic(A, _to_frame(Q, n); kwargs...))
istransverselyisotropic(t::AbstractTens; kwargs...) = istransverselyisotropic(first(_orthonormal_components(t)); kwargs...)
isorthotropic(t::AbstractTens, frame; kwargs...) =
    _in_frame(t, (A, Q) -> isorthotropic(A, _to_frame(Q, frame); kwargs...))
# `AbstractTens <: AbstractArray`, so this wrapper and the closed-form worker
# taking an `OrthonormalBasis{3}` are each more specific on one argument and
# less on the other. Unwrapping first is what the wrapper is for; saying it for
# this combination resolves the pair.
isorthotropic(t::AbstractTens, frame::OrthonormalBasis{3}; kwargs...) =
    _in_frame(t, (A, Q) -> isorthotropic(A, _to_frame(Q, frame); kwargs...))
isorthotropic(t::AbstractTens; kwargs...) = isorthotropic(first(_orthonormal_components(t)); kwargs...)
iscubic(t::AbstractTens, frame; kwargs...) =
    _in_frame(t, (A, Q) -> iscubic(A, _to_frame(Q, frame); kwargs...))
# Same disambiguation as `isorthotropic` just above: `AbstractTens <: AbstractArray`,
# so the wrapper and the worker are each more specific on one argument.
iscubic(t::AbstractTens, frame::OrthonormalBasis{3}; kwargs...) =
    _in_frame(t, (A, Q) -> iscubic(A, _to_frame(Q, frame); kwargs...))
iscubic(t::AbstractTens; kwargs...) = iscubic(first(_orthonormal_components(t)); kwargs...)

# ──────────────────────────────────────────────────────────────────────────────
# Public aliases for the per-block Kelvin-Mandel ↔ symmetry-parameter
# conversions `project` uses internally (`_project_TI_KM`, `_build_TI_KM`,
# `_project_ORTHO_KM`, `_build_ORTHO_KM`), so that callers working directly on
# 6×6 Kelvin-Mandel blocks (e.g. MeanFieldHomogenization's per-timestep ALV
# Volterra matrices) can reuse the same closed-form projection without going
# through a full `AbstractTens`.
# ──────────────────────────────────────────────────────────────────────────────

"""
    ti_params_from_mandel(C::AbstractMatrix) → (ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆)

Project a 6×6 Kelvin-Mandel matrix (in a frame where `e₃` is the TI symmetry
axis) onto the major-symmetric TI subspace, returning the 5 Walpole
coefficients. See [`project`](@ref) for the tensor-level equivalent.
"""
const ti_params_from_mandel = _project_TI_KM

"""
    ti2_params_from_matrix(M::AbstractMatrix) → (a, b, c)

The transversely isotropic part of a 3×3 matrix written in a frame whose third
axis is the symmetry axis: the transverse coefficient `a = (M₁₁ + M₂₂)/2`, the
axial one `b = M₃₃` and the in-plane antisymmetric one `c = (M₂₁ − M₁₂)/2`, those
of `a (𝟏 − n⊗n) + b n⊗n + c [n]ₓ` (see `TensTI{2}`). For a symmetric matrix
`c = 0`. The order-2 counterpart of [`ti_params_from_mandel`](@ref).
"""
ti2_params_from_matrix(M::AbstractMatrix) = ((M[1, 1] + M[2, 2]) / 2, M[3, 3], (M[2, 1] - M[1, 2]) / 2)

"""
    mandel_from_ti_params(ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆) → SMatrix{6,6}

Build a 6×6 Kelvin-Mandel matrix (in the frame where `e₃` is the TI axis)
from 5 major-symmetric Walpole coefficients — the reciprocal of
[`ti_params_from_mandel`](@ref).
"""
const mandel_from_ti_params = _build_TI_KM

"""
    ortho_params_from_mandel(C::AbstractMatrix) → NTuple{9}

Extract the 9 orthotropic parameters `(C₁₁,C₂₂,C₃₃,C₁₂,C₁₃,C₂₃,C₄₄,C₅₅,C₆₆)`
from a 6×6 Kelvin-Mandel matrix in the material frame. See
[`project`](@ref) for the tensor-level equivalent.
"""
const ortho_params_from_mandel = _project_ORTHO_KM

"""
    cubic_params_from_mandel(C::AbstractMatrix) → (α, β, γ)

Extract the three cubic coefficients on `(𝕁, 𝔼, 𝕋)` from a 6×6 Kelvin-Mandel
matrix in the cube frame. See [`project`](@ref) for the tensor-level
equivalent, and [`cubic_constants`](@ref) for the engineering constants.
"""
const cubic_params_from_mandel = _project_CUBIC_KM

"""
    mandel_from_cubic_params(α, β, γ) → SMatrix{6,6}

Build a 6×6 Kelvin-Mandel matrix in the cube frame from the three cubic
coefficients — the reciprocal of [`cubic_params_from_mandel`](@ref).
"""
const mandel_from_cubic_params = _build_CUBIC_KM

"""
    mandel_from_ortho_params(C₁₁, C₂₂, C₃₃, C₁₂, C₁₃, C₂₃, C₄₄, C₅₅, C₆₆) → SMatrix{6,6}

Build a 6×6 Kelvin-Mandel matrix in the material frame from 9 orthotropic
parameters — the reciprocal of [`ortho_params_from_mandel`](@ref).
"""
const mandel_from_ortho_params = _build_ORTHO_KM

# ──────────────────────────────────────────────────────────────────────────────
# Best-fit projections — the "give me the closest tensor of that class" façade
# ──────────────────────────────────────────────────────────────────────────────
#
# `project` returns the triple `(projection, distance, relative distance)`,
# which is what a diagnostic wants. These three names return the projection
# alone, which is what a caller extracting parameters wants, and read as the
# echoes `.paramsym(sym = …)` they mirror.
#
# They are PROJECTIONS, not averages. Whenever the question is "what does this
# tensor become once its orientation is averaged out", the answer is
# `isotropify` / `transverse_isotropify`, which are exact and lossless on the
# invariant subspace. See `tens_rotational_average.jl`.

"""
    best_fit_iso(t::AbstractTens) -> TensISO

Orthogonal (Frobenius) projection of `t` onto `span(𝕁, 𝕂)` at order 4, onto
`span(𝟏)` at order 2 — the projection component of
[`project`](@ref)`(Val(:ISO), t)`.

For a minor-symmetric tensor this coincides with the exact SO(3) average
[`isotropify`](@ref); for a general one it does not, and the average is the
meaningful object.
"""
best_fit_iso(t::AbstractTens) = project(t, Val(:ISO))[1]

"""
    best_fit_ti(t::AbstractTens{4,3}, axis) -> TensTI{4,T,5}
    best_fit_ti(t::AbstractTens{2,3}, axis) -> TensTI{2,T,2}

Orthogonal (Frobenius) projection of `t` onto the **major-symmetric** Walpole
span about `axis` — the projection component of [`project`](@ref)`(Val(:TI),
t, axis)`.

!!! warning "A projection, not an orientation average"
    It forces major symmetry, `(ℓ₃+ℓ₄)/2`, and drops the antisymmetric
    azimuthal couplings `ℓ₇`, `ℓ₈`. On a concentration tensor that content is
    physical and its loss is silent. Use [`transverse_isotropify`](@ref) to
    average, and this only to report parameters.
"""
best_fit_ti(t::AbstractTens{4, 3}, axis) = project(t, Val(:TI), axis)[1]
best_fit_ti(t::AbstractTens{2, 3}, axis) = project(t, Val(:TI), axis)[1]

"""
    best_fit_ortho(t::AbstractTens{4,3}, frame) -> TensOrtho
    best_fit_ortho(t::AbstractTens{2,3}, frame) -> TensOrtho{2}

Orthogonal (Frobenius) projection of `t` onto the orthotropic span in the given
material `frame` — the projection component of [`project`](@ref)`(
Val(:ORTHO), t, frame)`. Same caveat as [`best_fit_ti`](@ref): a fit, not an
average.
"""
best_fit_ortho(t::AbstractTens{4, 3}, frame) = project(t, Val(:ORTHO), frame)[1]
best_fit_ortho(t::AbstractTens{2, 3}, frame) = project(t, Val(:ORTHO), frame)[1]

"""
    best_fit_cubic(t::AbstractTens{4,3}, frame) -> TensCubic
    best_fit_cubic(t::AbstractTens{2,3}, frame) -> TensISO{2,3}

Orthogonal (Frobenius) projection of `t` onto the cubic span with cube axes
`frame` — the projection component of [`project`](@ref)`(Val(:CUBIC), t,
frame)`.

Unlike [`best_fit_ti`](@ref), this loses nothing that a *major-symmetric*
result would have kept: the cubic class is automatically major-symmetric, so
there is no antisymmetric content to drop silently. What it does discard is
everything outside the class, and `project` reports how much.

At order two it returns an isotropic tensor, the two classes coinciding there.
"""
best_fit_cubic(t::AbstractTens{4, 3}, frame) = project(t, Val(:CUBIC), frame)[1]
best_fit_cubic(t::AbstractTens{2, 3}, frame) = project(t, Val(:CUBIC), frame)[1]

# ── Exports ──────────────────────────────────────────────────────────────────

export project
export ti_params_from_mandel, mandel_from_ti_params, ortho_params_from_mandel, mandel_from_ortho_params, ti2_params_from_matrix
export cubic_params_from_mandel, mandel_from_cubic_params
export best_fit_iso, best_fit_ti, best_fit_ortho, best_fit_cubic
export iscubic
