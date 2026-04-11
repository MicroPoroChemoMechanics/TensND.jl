##############################################################################
# TensNDNLoptExt — rotation-optimized tensor projections using NLopt         #
#                                                                            #
# This package extension overrides the error-throwing fallbacks in           #
# tens_projection.jl with actual optimizers:                                 #
#   - proj_tens(Val(:TI), A)      — optimize TI axis (θ, ϕ)                 #
#   - proj_tens(Val(:ORTHO), A)   — optimize ORTHO frame (θ, ϕ, ψ)          #
# for both 4th-order and 2nd-order tensors.                                  #
#                                                                            #
# Strategy (from ECHOES tensor_approx.h):                                    #
#   Pass 1: GD_MLSL (global) + LD_TNEWTON (local), coarse tolerances        #
#   Pass 2: LD_TNEWTON (local), fine tolerances                              #
#                                                                            #
# Gradient computed via ForwardDiff.                                         #
##############################################################################

module TensNDNLoptExt

using TensND
using NLopt
using ForwardDiff
using StaticArrays

import TensND: proj_tens, _rot3_raw, _KM_rotation, _KM_of_array,
    _project_TI_KM, _build_TI_KM,
    _project_ORTHO_KM, _build_ORTHO_KM,
    _frobenius, _n_from_angles, _angles_from_n,
    _extract_vec

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

# ── Two-pass NLopt optimization ──────────────────────────────────────────────

"""
    _optimize_angles(obj, n_angles, x0; kwargs...) → (x_opt, j_opt)

Two-pass NLopt optimization matching ECHOES tensor_approx.h strategy:
- Pass 1: GD_MLSL (global) + LD_TNEWTON (local), coarse tolerances
- Pass 2: LD_TNEWTON (local), fine tolerances

Bounds: θ ∈ [0, π/2], ϕ ∈ [0, 2π] (TI) or [0, π] (ORTHO), ψ ∈ [0, π/2].
"""
function _optimize_angles(obj, n_angles::Int, x0::Vector{Float64})
    lb = zeros(n_angles)
    ub = fill(π / 2, n_angles)
    if n_angles >= 2
        ub[2] = n_angles > 2 ? π : 2π   # ϕ bound: 2π for TI (2 angles), π for ORTHO (3 angles)
    end

    x = copy(x0)

    # NLopt objective with ForwardDiff gradient (shared by both passes)
    nlopt_obj = (x_vec, grad_vec) -> begin
        if length(grad_vec) > 0
            grad_vec .= ForwardDiff.gradient(obj, x_vec)
        end
        return obj(x_vec)
    end

    # ── Pass 1: Global search ──
    try
        opt1 = NLopt.Opt(:GD_MLSL, n_angles)
        NLopt.lower_bounds!(opt1, lb)
        NLopt.upper_bounds!(opt1, ub)
        NLopt.xtol_rel!(opt1, 1.0e-2)
        NLopt.xtol_abs!(opt1, 1.0e-2)
        NLopt.ftol_rel!(opt1, 1.0e-3)
        NLopt.maxeval!(opt1, 1000)

        local_opt = NLopt.Opt(:LD_TNEWTON, n_angles)
        NLopt.lower_bounds!(local_opt, lb)
        NLopt.upper_bounds!(local_opt, ub)
        NLopt.xtol_rel!(local_opt, 1.0e-3)
        NLopt.xtol_abs!(local_opt, 1.0e-3)
        NLopt.ftol_rel!(local_opt, 1.0e-3)
        NLopt.maxeval!(local_opt, 1000)
        NLopt.local_optimizer!(opt1, local_opt)

        NLopt.min_objective!(opt1, nlopt_obj)

        (minf, minx, ret) = NLopt.optimize(opt1, x)
        x = minx
    catch e
        @debug "NLopt global optimizer failed; proceeding to local refinement" exception = (e, catch_backtrace())
    end

    # ── Pass 2: Local refinement ──
    try
        opt2 = NLopt.Opt(:LD_TNEWTON, n_angles)
        NLopt.lower_bounds!(opt2, lb)
        NLopt.upper_bounds!(opt2, ub)
        NLopt.xtol_rel!(opt2, 1.0e-6)
        NLopt.xtol_abs!(opt2, 1.0e-6)
        NLopt.ftol_rel!(opt2, 1.0e-6)
        NLopt.maxeval!(opt2, 100)

        NLopt.min_objective!(opt2, nlopt_obj)

        (minf, minx, ret) = NLopt.optimize(opt2, x)
        x = minx
    catch e
        @debug "NLopt local optimizer failed; returning best x found so far" exception = (e, catch_backtrace())
    end

    return x
end

# ── proj_tens: TI, order 4, optimized ────────────────────────────────────────

"""
    proj_tens(::Val{:TI}, A::AbstractArray{T,4}) where {T<:AbstractFloat}

Find the best TI approximation of a 4th-order tensor `A` by optimizing the
symmetry axis over all directions. Uses NLopt (GD_MLSL + LD_TNEWTON).

Returns `(B::TensWalpole{T,5}, d, drel)`.

# Examples
```julia
julia> using NLopt

julia> n = [1/√3, 1/√3, 1/√3];

julia> C = tensTI(10., 3., 2.5, 12., 2., n);

julia> B, d, drel = proj_tens(:TI, getarray(C));

julia> drel < 1e-6
true
```
"""
function TensND.proj_tens(::Val{:TI}, A::AbstractArray{T, 4}) where {T <: AbstractFloat}
    C_KM = _KM_of_array(A)
    sqnorm_C = sum(x -> x^2, C_KM)
    if sqnorm_C ≈ zero(T)
        z = zero(T)
        n = (z, z, one(T))
        return TensWalpole(z, z, z, z, z, n), z, z
    end

    obj = x -> _obj_TI4(x, C_KM, sqnorm_C)
    x0 = [T(π / 4), T(π / 4)]
    x_opt = _optimize_angles(obj, 2, Float64.(x0))

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

julia> A = TensTI{2}(5.0, 8.0, n); Amat = getarray(A);

julia> B, d, drel = proj_tens(:TI, Amat);

julia> drel < 1e-6
true
```
"""
function TensND.proj_tens(::Val{:TI}, A::AbstractArray{T, 2}) where {T <: AbstractFloat}
    sqnorm_A = sum(x -> x^2, A)
    if sqnorm_A ≈ zero(T)
        z = zero(T)
        n = (z, z, one(T))
        return TensTI{2}(z, z, n), z, z
    end

    obj = x -> _obj_TI2(x, A, sqnorm_A)
    x0 = [T(π / 4), T(π / 4)]
    x_opt = _optimize_angles(obj, 2, Float64.(x0))

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

julia> B, d, drel = proj_tens(:ORTHO, getarray(t));

julia> drel < 1e-4
true
```
"""
function TensND.proj_tens(::Val{:ORTHO}, A::AbstractArray{T, 4}) where {T <: AbstractFloat}
    C_KM = _KM_of_array(A)
    sqnorm_C = sum(x -> x^2, C_KM)
    if sqnorm_C ≈ zero(T)
        z = zero(T)
        frame = CanonicalBasis{3, T}()
        return TensOrtho(z, z, z, z, z, z, z, z, z, frame), z, z
    end

    obj = x -> _obj_ORTHO4(x, C_KM, sqnorm_C)
    x0 = [T(π / 4), T(π / 4), T(π / 4)]
    x_opt = _optimize_angles(obj, 3, Float64.(x0))

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
function TensND.proj_tens(::Val{:ORTHO}, A::AbstractArray{T, 2}) where {T <: AbstractFloat}
    sqnorm_A = sum(x -> x^2, A)
    if sqnorm_A ≈ zero(T)
        z = zero(T)
        return zeros(T, 3, 3), z, z
    end

    obj = x -> _obj_ORTHO2(x, A, sqnorm_A)
    x0 = [T(π / 4), T(π / 4), T(π / 4)]
    x_opt = _optimize_angles(obj, 3, Float64.(x0))

    frame = RotatedBasis(x_opt[1], x_opt[2], x_opt[3])
    return proj_tens(Val(:ORTHO), A, frame)
end

end # module TensNDNLoptExt
