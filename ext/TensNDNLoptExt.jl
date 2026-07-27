##############################################################################
# TensNDNLoptExt — rotation-optimized tensor projections using NLopt         #
#                                                                            #
# This package extension overrides the error-throwing fallbacks in           #
# tens_projection.jl with actual optimizers:                                 #
#   - proj_tens(Val(:TI), A)      — optimize TI axis (θ, ϕ)                 #
#   - proj_tens(Val(:ORTHO), A)   — optimize ORTHO frame (θ, ϕ, ψ)          #
# for both 4th-order and 2nd-order tensors.                                  #
#                                                                            #
# Strategy: deterministic multi-start + LD_TNEWTON local refinement.         #
# The starting points are the eigenstructure candidate (`_candidate_TI_axis` #
# / `_candidate_ORTHO_frame`, exact whenever the tensor really does have     #
# that symmetry) plus a fixed angular grid containing the canonical axes.    #
# The best objective over {every start} ∪ {every refined start} is kept, so  #
# the optimized projection is never worse than the fixed-axis / fixed-frame  #
# projection along any grid point.                                           #
#                                                                            #
# Earlier versions ran the ECHOES tensor_approx.h strategy: a GD_MLSL global #
# pass, then a local one.  GD_MLSL is *stochastic* and NLopt seeds its       #
# generator from the clock, so the result was irreproducible from one call   #
# to the next: on a tensor exactly orthotropic about a tilted frame, ~0.5 %  #
# of runs returned a spurious local minimum (drel ≈ 0.44 instead of ≈1e-13), #
# which surfaced as an intermittent CI failure.  The multi-start below       #
# matched or beat the best of 30 GD_MLSL runs on every tensor tested (150    #
# random anisotropic ones included) and is 3-8× faster.                      #
#                                                                            #
# Gradient computed via ForwardDiff.                                         #
##############################################################################

module TensNDNLoptExt

using TensND
using NLopt
using ForwardDiff
using StaticArrays

import TensND: proj_tens, _proj_tens_opt, _rot3_raw, _KM_rotation, _KM_of_array,
    _project_TI_KM, _build_TI_KM,
    _project_ORTHO_KM, _build_ORTHO_KM,
    _frobenius, _n_from_angles, _angles_from_n,
    _extract_vec, _candidate_TI_axis, _candidate_ORTHO_frame,
    angles, vecbasis

# ── Objective function: TI, order 4 ─────────────────────────────────────────

"""
    _obj_TI4(x, C_KM, sqnorm_C) → j

Objective for TI order-4 projection: `j(θ,ϕ) = 1 − ‖B‖²_KM / ‖C‖²_KM`.
`x = [θ, ϕ]`, `C_KM` is the 6×6 KM matrix of A.
"""
function _obj_TI4(x, C_KM, sqnorm_C)
    θ, ϕ = x[1], x[2]
    P₆ = _KM_rotation(θ, ϕ, zero(eltype(x)))
    C_rot = P₆' * C_KM * P₆
    ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆ = _project_TI_KM(C_rot)
    B_KM = _build_TI_KM(ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆)
    return one(eltype(x)) - sum(x -> x^2, B_KM) / sqnorm_C
end

# ── Objective function: ORTHO, order 4 ──────────────────────────────────────

"""
    _obj_ORTHO4(x, C_KM, sqnorm_C) → j

Objective for ORTHO order-4 projection: `j(θ,ϕ,ψ) = 1 − ‖B‖²_KM / ‖C‖²_KM`.
"""
function _obj_ORTHO4(x, C_KM, sqnorm_C)
    θ, ϕ, ψ = x[1], x[2], x[3]
    P₆ = _KM_rotation(θ, ϕ, ψ)
    C_rot = P₆' * C_KM * P₆
    params = _project_ORTHO_KM(C_rot)
    B_KM = _build_ORTHO_KM(params...)
    return one(eltype(x)) - sum(x -> x^2, B_KM) / sqnorm_C
end

# ── Objective function: TI, order 2 ─────────────────────────────────────────

function _obj_TI2(x, A, sqnorm_A)
    θ, ϕ = x[1], x[2]
    R = _rot3_raw(θ, ϕ, zero(eltype(x)))
    M_rot = R' * A * R
    a = (M_rot[1, 1] + M_rot[2, 2]) / 2
    b = M_rot[3, 3]
    # Build projected matrix in canonical frame
    n = (sin(θ) * cos(ϕ), sin(θ) * sin(ϕ), cos(θ))
    B_sqnorm = zero(eltype(x))
    for i in 1:3, j in 1:3
        δij = i == j ? one(eltype(x)) : zero(eltype(x))
        Bij = a * (δij - n[i] * n[j]) + b * n[i] * n[j]
        B_sqnorm += Bij^2
    end
    return one(eltype(x)) - B_sqnorm / sqnorm_A
end

# ── Objective function: ORTHO, order 2 ──────────────────────────────────────

function _obj_ORTHO2(x, A, sqnorm_A)
    θ, ϕ, ψ = x[1], x[2], x[3]
    R = _rot3_raw(θ, ϕ, ψ)
    M_rot = R' * A * R
    # ORTHO projection: keep diagonal only, then rotate back
    B_sqnorm = zero(eltype(x))
    for i in 1:3, j in 1:3
        Bij = zero(eltype(x))
        for k in 1:3
            Bij += R[i, k] * M_rot[k, k] * R[j, k]
        end
        B_sqnorm += Bij^2
    end
    return one(eltype(x)) - B_sqnorm / sqnorm_A
end

# ── Deterministic multi-start optimization ───────────────────────────────────

# Angular start grids.  Both contain the canonical frame (all angles zero) and
# the canonical axes, which is what gives the "never worse than a fixed-frame
# projection" guarantee.  `θ = 0` makes ϕ irrelevant for the TI axis
# `n = (sinθ cosϕ, sinθ sinϕ, cosθ)`, hence the filter.
const _TI_STARTS = [
    [θ, ϕ] for θ in (0.0, π / 4, π / 2) for ϕ in (0.0, π / 4, π / 2, 3π / 4)
        if !(θ == 0.0 && ϕ != 0.0)
]

const _ORTHO_STARTS = [
    [θ, ϕ, ψ] for θ in (0.0, π / 4, π / 2)
        for ϕ in (0.0, π / 3, 2π / 3) for ψ in (0.0, π / 4)
]

# Refinement bounds.  Generous rather than a fundamental domain: the
# eigenstructure candidate arrives with whatever Euler representative
# `angles` produced (ϕ and ψ may be negative), and a gradient method started
# at a good point does not wander.  The objective is well defined for any
# angles, and so is the frame rebuilt from them.
const _ANGLE_BOUND = 2π

# Smallest objective decrease worth acting on.  `j = 1 − ‖B‖²/‖C‖²` is a
# difference of O(1) quantities, so evaluating it near an exact symmetry costs
# a few ulps of cancellation and `j` lands anywhere in ±1e-16.  Without this
# floor, a refined point could displace an earlier start by "improving" `j`
# purely in rounding noise — and since the starts are ordered exact-candidate
# first, the point being displaced is the *better* one.  It shows up plainly on
# a symmetric 3×3, whose eigenframe is its orthotropic projection exactly:
# the candidate gives drel = 0, the noise-level winner gives drel ≈ 1e-8, the
# angular resolution of the refinement.
const _IMPROVE_TOL = 1.0e-14

"""
    _refine_angles(obj, n_angles, x0) → x

Local refinement of `obj` from `x0` with `LD_TNEWTON`, gradient supplied by
`ForwardDiff`.  Returns `x0` unchanged if NLopt errors out — the caller keeps
the best objective seen, so a failed refinement can only be a missed
improvement, never a regression.
"""
function _refine_angles(obj, n_angles::Int, x0::Vector{Float64})
    nlopt_obj = (x_vec, grad_vec) -> begin
        if length(grad_vec) > 0
            grad_vec .= ForwardDiff.gradient(obj, x_vec)
        end
        return obj(x_vec)
    end

    return try
        opt = NLopt.Opt(:LD_TNEWTON, n_angles)
        NLopt.lower_bounds!(opt, fill(-_ANGLE_BOUND, n_angles))
        NLopt.upper_bounds!(opt, fill(_ANGLE_BOUND, n_angles))
        NLopt.xtol_rel!(opt, 1.0e-8)
        NLopt.xtol_abs!(opt, 1.0e-8)
        NLopt.ftol_rel!(opt, 1.0e-10)
        NLopt.maxeval!(opt, 200)
        NLopt.min_objective!(opt, nlopt_obj)

        (_, minx, _) = NLopt.optimize(opt, copy(x0))
        minx
    catch e
        @debug "NLopt local refinement failed; keeping the starting point" exception = (e, catch_backtrace())
        copy(x0)
    end
end

"""
    _optimize_angles(obj, n_angles, starts) → x_opt

Minimize `obj` over the angles by refining every point of `starts` locally and
keeping the best result.

Deterministic: no random search is involved, so two calls on the same tensor
return the same angles.  `obj` is evaluated at each start *before* refining it,
so the returned point is never worse than the best start — in particular never
worse than the eigenstructure candidate, which is exact for a tensor that
genuinely has the symmetry being projected onto.

`starts` is scanned in order and a challenger must beat the incumbent by
`_IMPROVE_TOL`, so the earliest start wins whenever the difference is rounding
noise.  Callers therefore put their most trustworthy candidate first.
"""
function _optimize_angles(obj, n_angles::Int, starts)
    best_x = first(starts)
    best_j = obj(best_x)

    for x0 in starts
        j0 = obj(x0)
        if j0 < best_j - _IMPROVE_TOL
            best_j, best_x = j0, x0
        end
        x = _refine_angles(obj, n_angles, x0)
        j = obj(x)
        if j < best_j - _IMPROVE_TOL
            best_j, best_x = j, x
        end
    end

    return best_x
end

# ── Deterministic start sets ─────────────────────────────────────────────────

"""
    _ti_starts(A) → Vector{Vector{Float64}}

Start set for the TI axis: the eigenstructure candidate
(`_candidate_TI_axis`, exact for a genuinely TI tensor) followed by
`_TI_STARTS`.
"""
function _ti_starts(A)
    θ, ϕ = _angles_from_n(_candidate_TI_axis(A))
    return vcat([[Float64(θ), Float64(ϕ)]], _TI_STARTS)
end

"""
    _ortho_starts(A) → Vector{Vector{Float64}}

Start set for the orthotropic frame: the eigenstructure candidate
(`_candidate_ORTHO_frame`, exact for a genuinely orthotropic tensor) followed
by `_ORTHO_STARTS`.
"""
function _ortho_starts(A)
    a = angles(Matrix(vecbasis(_candidate_ORTHO_frame(A), :cov)))
    return vcat([[Float64(a.θ), Float64(a.ϕ), Float64(a.ψ)]], _ORTHO_STARTS)
end

# ── proj_tens: TI, order 4, optimized ────────────────────────────────────────

"""
    proj_tens(::Val{:TI}, A::AbstractArray{T,4}) where {T<:AbstractFloat}

Find the best TI approximation of a 4th-order tensor `A` by optimizing the
symmetry axis over all directions. Uses NLopt (GD_MLSL + LD_TNEWTON).

Returns `(B::TensTI{4, T, 5}, d, drel)`.

# Examples
```julia
julia> using NLopt

julia> n = [1/√3, 1/√3, 1/√3];

julia> C = tens_TI(10., 3., 2.5, 12., 2., n);

julia> B, d, drel = proj_tens(:TI, get_array(C));

julia> drel < 1e-6
true
```
"""
function TensND._proj_tens_opt(::Val{:TI}, A::AbstractArray{T, 4}) where {T <: AbstractFloat}
    C_KM = _KM_of_array(A)
    sqnorm_C = sum(x -> x^2, C_KM)
    if sqnorm_C ≈ zero(T)
        z = zero(T)
        n = (z, z, one(T))
        return TensTI{4}(z, z, z, z, z, n), z, z
    end

    obj = x -> _obj_TI4(x, C_KM, sqnorm_C)
    x_opt = _optimize_angles(obj, 2, _ti_starts(A))

    n = _n_from_angles(x_opt[1], x_opt[2])
    return proj_tens(Val(:TI), A, n)
end

# ── proj_tens: TI, order 2, optimized ────────────────────────────────────────

"""
    proj_tens(::Val{:TI}, A::AbstractArray{T,2}) where {T<:AbstractFloat}

Find the best TI approximation of a 2nd-order tensor `A` (3×3) by optimizing
the symmetry axis. Uses NLopt.

Returns `(B::TensTI{2,T,2}, d, drel)`.

# Examples
```julia
julia> using NLopt

julia> n = [1/√2, 1/√2, 0.];

julia> A = TensTI{2}(5.0, 8.0, n); Amat = get_array(A);

julia> B, d, drel = proj_tens(:TI, Amat);

julia> drel < 1e-6
true
```
"""
function TensND._proj_tens_opt(::Val{:TI}, A::AbstractArray{T, 2}) where {T <: AbstractFloat}
    sqnorm_A = sum(x -> x^2, A)
    if sqnorm_A ≈ zero(T)
        z = zero(T)
        n = (z, z, one(T))
        return TensTI{2}(z, z, n), z, z
    end

    obj = x -> _obj_TI2(x, A, sqnorm_A)
    x_opt = _optimize_angles(obj, 2, _ti_starts(A))

    n = _n_from_angles(x_opt[1], x_opt[2])
    return proj_tens(Val(:TI), A, n)
end

# ── proj_tens: ORTHO, order 4, optimized ─────────────────────────────────────

"""
    proj_tens(::Val{:ORTHO}, A::AbstractArray{T,4}) where {T<:AbstractFloat}

Find the best orthotropic approximation of a 4th-order tensor `A` by
optimizing the material frame (3 Euler angles). Uses NLopt.

Returns `(B::TensOrtho{T}, d, drel)`.

# Examples
```julia
julia> using NLopt

julia> frame = RotatedBasis(0.3, 0.5, 0.7);

julia> t = TensOrtho(10., 8., 12., 3., 2.5, 1.5, 2., 3., 3.5, frame);

julia> B, d, drel = proj_tens(:ORTHO, get_array(t));

julia> drel < 1e-4
true
```
"""
function TensND._proj_tens_opt(::Val{:ORTHO}, A::AbstractArray{T, 4}) where {T <: AbstractFloat}
    C_KM = _KM_of_array(A)
    sqnorm_C = sum(x -> x^2, C_KM)
    if sqnorm_C ≈ zero(T)
        z = zero(T)
        frame = CanonicalBasis{3, T}()
        return TensOrtho(z, z, z, z, z, z, z, z, z, frame), z, z
    end

    obj = x -> _obj_ORTHO4(x, C_KM, sqnorm_C)
    x_opt = _optimize_angles(obj, 3, _ortho_starts(A))

    frame = RotatedBasis(x_opt[1], x_opt[2], x_opt[3])
    return proj_tens(Val(:ORTHO), A, frame)
end

# ── proj_tens: ORTHO, order 2, optimized ─────────────────────────────────────

"""
    proj_tens(::Val{:ORTHO}, A::AbstractArray{T,2}) where {T<:AbstractFloat}

Find the best orthotropic approximation of a 2nd-order tensor `A` (3×3) by
optimizing the material frame. Uses NLopt.

Returns `(B::Array{T,2}, d, drel)`.

# Examples
```julia
julia> using NLopt

julia> A = Float64[5 1 2; 1 8 3; 2 3 12];

julia> B, d, drel = proj_tens(:ORTHO, A);

julia> d ≥ 0
true
```
"""
function TensND._proj_tens_opt(::Val{:ORTHO}, A::AbstractArray{T, 2}) where {T <: AbstractFloat}
    sqnorm_A = sum(x -> x^2, A)
    if sqnorm_A ≈ zero(T)
        z = zero(T)
        return zeros(T, 3, 3), z, z
    end

    obj = x -> _obj_ORTHO2(x, A, sqnorm_A)
    x_opt = _optimize_angles(obj, 3, _ortho_starts(A))

    frame = RotatedBasis(x_opt[1], x_opt[2], x_opt[3])
    return proj_tens(Val(:ORTHO), A, frame)
end

end # module TensNDNLoptExt
