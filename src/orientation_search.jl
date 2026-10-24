# ── Orientation search for the symmetry projections ──────────────────────────
#
# The projections of `tens_projection.jl` onto the transversely isotropic,
# orthotropic and cubic classes take the axis or the frame as given. Finding it
# is minimizing, over two or three Euler angles, the relative distance of the
# tensor to the class turned by those angles,
#
#     j(angles) = 1 − ‖P(angles) C‖² / ‖C‖²,
#
# P being the orthogonal projection onto the turned class, so that
# j = (‖C − P C‖ / ‖C‖)². j is smooth and periodic in the angles.
#
# Strategy, for a fourth-order tensor: the objective is screened on a fixed,
# quasi-uniform set of orientations — cheap, an evaluation costs about a third
# of a microsecond — and Newton's method on the angles refines the few best
# points, the eigenstructure candidate (`_candidate_TI_axis` /
# `_candidate_ORTHO_frame`, exact whenever the tensor has the symmetry) and the
# canonical orientation. The Hessian is exact (ForwardDiff, a 2×2 or 3×3
# matrix), so the refinement converges quadratically near a minimum, and a
# backtracking line search keeps it descending away from one. The best
# objective over every refined point is kept. No random search is involved:
# the same tensor gives the same orientation on every call. On 150 random
# anisotropic stiffnesses the result matched, to rounding, the minimum of a
# dense search (3,000 axes for TI, 7,200 frames for the orthotropic and cubic
# classes), which a multi-start from a dozen points missed for one tensor in
# twenty on the TI axis.
#
# For a second-order tensor no search is needed: the best orthotropic frame is
# the eigenframe of the symmetric part, and the best TI axis an eigenvector of
# one of its two extreme eigenvalues (`_best_TI2_axis`).

# ── Objectives ───────────────────────────────────────────────────────────────

"""
    _obj_TI4(x, C_KM, sqnorm_C) → j

Objective of the TI projection of a fourth-order tensor:
`j(θ,ϕ) = 1 − ‖B‖²_KM / ‖C‖²_KM`, with `x = (θ, ϕ)` and `C_KM` the 6×6
Kelvin-Mandel matrix of the tensor.
"""
function _obj_TI4(x, C_KM, sqnorm_C)
    θ, ϕ = x[1], x[2]
    P₆ = mandel_rotation(θ, ϕ, zero(eltype(x)))
    C_rot = P₆' * C_KM * P₆
    ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆ = _project_TI_KM(C_rot)
    B_KM = _build_TI_KM(ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆)
    return one(eltype(x)) - sum(x -> x^2, B_KM) / sqnorm_C
end

"""
    _obj_ORTHO4(x, C_KM, sqnorm_C) → j

Objective of the orthotropic projection of a fourth-order tensor:
`j(θ,ϕ,ψ) = 1 − ‖B‖²_KM / ‖C‖²_KM`.
"""
function _obj_ORTHO4(x, C_KM, sqnorm_C)
    θ, ϕ, ψ = x[1], x[2], x[3]
    P₆ = mandel_rotation(θ, ϕ, ψ)
    C_rot = P₆' * C_KM * P₆
    B_KM = _build_ORTHO_KM(_project_ORTHO_KM(C_rot)...)
    return one(eltype(x)) - sum(x -> x^2, B_KM) / sqnorm_C
end

"""
    _obj_CUBIC4(x, C_KM, sqnorm_C) → j

Objective of the cubic projection with a free cube orientation:
`j(θ,ϕ,ψ) = 1 − ‖B‖²_KM / ‖C‖²_KM`, the same form as the TI and orthotropic
objectives, since all three are orthogonal projections onto a turned subspace.

Smooth in the angles: the octahedral group is discrete, but the orientation of
the cube is an ordinary rotation. The group shows up only as a 24-fold
degeneracy of the minimum, which costs the search nothing — every one of the
24 frames describes the same tensor with the same coefficients.
"""
function _obj_CUBIC4(x, C_KM, sqnorm_C)
    θ, ϕ, ψ = x[1], x[2], x[3]
    P₆ = mandel_rotation(θ, ϕ, ψ)
    C_rot = P₆' * C_KM * P₆
    B_KM = _build_CUBIC_KM(_project_CUBIC_KM(C_rot)...)
    return one(eltype(x)) - sum(x -> x^2, B_KM) / sqnorm_C
end

# ── Local refinement: Newton's method on the angles ──────────────────────────

const _NEWTON_MAXITER = 50

"""
    _refine_angles(obj, x0::SVector) → SVector

Minimize `obj` from `x0` by Newton's method with the exact Hessian
(`ForwardDiff`) and a backtracking (Armijo) line search.

Where the Hessian is not positive definite — near a saddle or a maximum of
`obj` — the step divides each eigen-component of the gradient by the absolute
value of its curvature, floored at a small fraction of the largest one: a
descent direction in every case, and the Newton step wherever the Hessian is
positive definite. The iteration stops when the gradient vanishes to rounding,
when a step no longer moves the angles or lowers the objective, or after
`_NEWTON_MAXITER` steps. The caller compares the objectives of every start and
refined point, so a refinement that stops early can only miss an improvement.
"""
function _refine_angles(obj, x0::SVector{N, Float64}) where {N}
    x = x0
    f = obj(x)
    for _ in 1:_NEWTON_MAXITER
        g = ForwardDiff.gradient(obj, x)
        maximum(abs, g) ≤ 1.0e-14 && break
        λ, V = eigen(Symmetric(ForwardDiff.hessian(obj, x)))
        floor = 1.0e-8 * max(maximum(abs, λ), 1.0)
        d = -(V * ((V' * g) ./ max.(abs.(λ), floor)))
        slope = dot(g, d)
        t = 1.0
        xt = x + d
        ft = obj(xt)
        while ft > f + 1.0e-4 * t * slope
            t /= 2
            t < 1.0e-12 && return x
            xt = x + t * d
            ft = obj(xt)
        end
        moved = maximum(abs, xt - x)
        lowered = f - ft
        x, f = xt, ft
        (moved ≤ 1.0e-13 || lowered ≤ 1.0e-17) && break
    end
    return x
end

"""
    _optimize_angles(obj, starts, grid, k) → SVector

Minimize `obj` over the angles: refine every point of `starts`, then the `k`
points of `grid` where `obj` is lowest, and keep the best result.

`obj` is evaluated at each point *before* refining it, so the returned point is
never worse than any of them — in particular never worse than the
eigenstructure candidate, which is exact for a tensor that genuinely has the
symmetry. The points are scanned in order and a challenger must beat the
incumbent by `_IMPROVE_TOL`, so the earliest one wins whenever the difference
is rounding noise: callers put their most trustworthy candidate first.
"""
function _optimize_angles(obj, starts, grid, k)
    values = map(obj, grid)
    screened = grid[partialsortperm(values, 1:k)]
    best_x = first(starts)
    best_j = obj(best_x)
    for x0 in Iterators.flatten((starts, screened))
        j0 = obj(x0)
        if j0 < best_j - _IMPROVE_TOL
            best_j, best_x = j0, x0
        end
        x = _refine_angles(obj, x0)
        j = obj(x)
        if j < best_j - _IMPROVE_TOL
            best_j, best_x = j, x
        end
    end
    return best_x
end

const _IMPROVE_TOL = 1.0e-14

# ── Orientations screened ────────────────────────────────────────────────────

# `n` axes spread over the upper hemisphere along a golden-angle spiral, as
# polar angles `(θ, ϕ)`: an axis and its opposite give the same class, so the
# hemisphere covers every axis once.
_hemisphere_spiral(n) = [
    SVector(acos((k + 0.5) / n), mod(k * π * (3 - sqrt(5)), 2π)) for k in 0:(n - 1)
]

# 200 axes for the TI class, and for the frames of the orthotropic and cubic
# classes 100 axes for the third vector times 4 turns ψ about it: a quarter turn
# about any axis of the frame maps it to an equivalent one, so ψ ∈ [0, π/2)
# suffices. The numbers are those for which the search matched the dense
# reference on every test tensor, with a margin (half as many TI axes missed
# once).
const _TI_GRID = _hemisphere_spiral(200)
const _FRAME_GRID = [SVector(a[1], a[2], ψ) for a in _hemisphere_spiral(100) for ψ in (0.0, π / 8, π / 4, 3π / 8)]
const _TI_SCREENED = 3
const _FRAME_SCREENED = 4

"""
    _ti_starts(A) → Vector{SVector{2,Float64}}

Points refined whatever the screening gives, the TI axis: the eigenstructure
candidate (`_candidate_TI_axis`, exact for a genuinely TI tensor), then `e₃`.
"""
function _ti_starts(A)
    θ, ϕ = _angles_from_n(_candidate_TI_axis(A))
    return [SVector(Float64(θ), Float64(ϕ)), SVector(0.0, 0.0)]
end

"""
    _ortho_starts(A) → Vector{SVector{3,Float64}}

Points refined whatever the screening gives, the orthotropic frame: the
eigenstructure candidate (`_candidate_ORTHO_frame`, exact for a genuinely
orthotropic tensor), then the canonical frame.
"""
function _ortho_starts(A)
    a = angles(Matrix(vecbasis(_candidate_ORTHO_frame(A), :cov)))
    return [SVector(Float64(a.θ), Float64(a.ϕ), Float64(a.ψ)), SVector(0.0, 0.0, 0.0)]
end

"""
    _cubic_starts(A) → Vector{SVector{3,Float64}}

Points refined whatever the screening gives, the cube orientation: the
canonical frame, and **no** eigenstructure candidate. Every second-order
contraction of a cubic tensor is isotropic, so no eigenframe carries any
information about where the cube points — handing `_candidate_ORTHO_frame` a
cubic tensor returns an arbitrary frame, measured at a relative residual of
0.099 against 3.5e-16 for the true one.
"""
_cubic_starts(_A) = [SVector(0.0, 0.0, 0.0)]

"""
    _best_TI2_axis(A) → NTuple{3}

Axis of the best TI approximation of the second-order tensor `A`, in closed
form. About the axis `n`, the projection is `a (𝟏 − n⊗n) + b n⊗n` with
`b = q = n·A·n` and `a = (tr A − q)/2`, so its squared norm
`(tr A − q)²/2 + q²` is a convex function of `q`, and `q` ranges over the
eigenvalues of the symmetric part of `A`: the maximum is at the smallest or at
the largest one, and the axis is the corresponding eigenvector.
"""
function _best_TI2_axis(A)
    λ, V = eigen(Symmetric((A + A') / 2))
    t = λ[1] + λ[2] + λ[3]
    f(q) = (t - q)^2 / 2 + q^2
    k = f(λ[1]) > f(λ[3]) ? 1 : 3
    return (V[1, k], V[2, k], V[3, k])
end

# ── The searches ─────────────────────────────────────────────────────────────

function _proj_tens_opt(::Val{:TI}, A::AbstractArray{T, 4}) where {T <: AbstractFloat}
    C_KM = _KM_of_array(A)
    sqnorm_C = sum(x -> x^2, C_KM)
    if sqnorm_C ≈ zero(T)
        z = zero(T)
        return TensTI{4}(z, z, z, z, z, (z, z, one(T))), z, z
    end
    x = _optimize_angles(x -> _obj_TI4(x, C_KM, sqnorm_C), _ti_starts(A), _TI_GRID, _TI_SCREENED)
    return proj_tens(Val(:TI), A, _n_from_angles(x[1], x[2]))
end

function _proj_tens_opt(::Val{:TI}, A::AbstractArray{T, 2}) where {T <: AbstractFloat}
    sqnorm_A = sum(x -> x^2, A)
    if sqnorm_A ≈ zero(T)
        z = zero(T)
        return TensTI{2}(z, z, (z, z, one(T))), z, z
    end
    return proj_tens(Val(:TI), A, _best_TI2_axis(A))
end

function _proj_tens_opt(::Val{:ORTHO}, A::AbstractArray{T, 4}) where {T <: AbstractFloat}
    C_KM = _KM_of_array(A)
    sqnorm_C = sum(x -> x^2, C_KM)
    if sqnorm_C ≈ zero(T)
        z = zero(T)
        return TensOrtho(z, z, z, z, z, z, z, z, z, CanonicalBasis{3, T}()), z, z
    end
    x = _optimize_angles(x -> _obj_ORTHO4(x, C_KM, sqnorm_C), _ortho_starts(A), _FRAME_GRID, _FRAME_SCREENED)
    return proj_tens(Val(:ORTHO), A, RotatedBasis(x[1], x[2], x[3]))
end

function _proj_tens_opt(::Val{:ORTHO}, A::AbstractArray{T, 2}) where {T <: AbstractFloat}
    sqnorm_A = sum(x -> x^2, A)
    if sqnorm_A ≈ zero(T)
        z = zero(T)
        return TensOrtho{2}(z, z, z, CanonicalBasis{3, T}()), z, z
    end
    # The eigenframe of the symmetric part: in it the projection keeps all of
    # that part, the most any frame can.
    return proj_tens(Val(:ORTHO), A, _candidate_ORTHO_frame(A))
end

function _proj_tens_opt(::Val{:CUBIC}, A::AbstractArray{T, 4}) where {T <: AbstractFloat}
    C_KM = _KM_of_array(A)
    sqnorm_C = sum(x -> x^2, C_KM)
    if sqnorm_C ≈ zero(T)
        z = zero(T)
        return TensCubic{T}((z, z, z), CanonicalBasis{3, T}()), z, z
    end
    x = _optimize_angles(x -> _obj_CUBIC4(x, C_KM, sqnorm_C), _cubic_starts(A), _FRAME_GRID, _FRAME_SCREENED)
    return proj_tens(Val(:CUBIC), A, RotatedBasis(x[1], x[2], x[3]))
end
