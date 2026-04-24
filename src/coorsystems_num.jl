using ForwardDiff

"""
    CoorSystemNum{dim,T<:Real} <: AbstractCoorSystem{dim,T}

Numerical coordinate system for automatic-differentiation-based differential operators.

Stores three point-wise functions:
- `χ_func(x)` : Lamé coefficients at point `x`
- `R_func(x)`  : rotation matrix (columns = unit vectors of the normalized basis) at `x`
- `Γ_func(x)`  : Christoffel symbols at `x` → `Array{T,3}` with `Γ[i,j,k] = Γᵏᵢⱼ`

All three functions accept an `AbstractVector` of coordinates and are called at each
evaluation point. Predefined constructors are available:

    coorsys_cartesian_num(T=Float64)
    coorsys_polar_num(T=Float64)
    coorsys_cylindrical_num(T=Float64)
    coorsys_spherical_num(T=Float64)

A generic constructor from an `OM::Function` (position vector) is also provided:

    CoorSystemNum(OM::Function, dim::Integer, T::Type=Float64)

Unlike `CoorSystemSym`, the normalized basis, natural basis vectors, Lamé coefficients,
and Christoffel symbols all depend on the evaluation point `x₀` and are accessed via
point-wise accessors `normalized_basis(CS, x₀)`, `natural_basis(CS, x₀)`, etc.
"""
struct CoorSystemNum{dim, T <: Real} <: AbstractCoorSystem{dim, T}
    χ_func::Function   # x -> Vector of Lamé coefficients
    R_func::Function   # x -> Matrix  (columns = unit vectors eᵢ in Cartesian frame)
    Γ_func::Function   # x -> Array{T,3}  Γ[i,j,k] = Γᵏᵢⱼ
end

@pure get_dim(::CoorSystemNum{dim}) where {dim} = dim

# ─────────────────────────────────────────────────────────────────
# Christoffel symbol helpers
# ─────────────────────────────────────────────────────────────────

"""
Compute Christoffel symbols for an **orthogonal** system from the Jacobian of
the Lamé-coefficient function.

Convention: `Γ[i,j,k] = Γᵏᵢⱼ`.
"""
function _christoffel_orthogonal(χ::AbstractVector, ∂χ::AbstractMatrix)
    # ∂χ[k,i] = ∂χₖ/∂xⁱ  (Jacobian: row=Lamé index, col=coordinate index)
    dim = length(χ)
    T = eltype(∂χ)
    Γ = zeros(T, dim, dim, dim)
    for i in 1:dim, j in 1:dim, k in 1:dim
        if i == j == k
            Γ[i, j, k] = ∂χ[k, i] / χ[k]
        elseif i == j && i != k
            Γ[i, j, k] = -(χ[i] / χ[k]^2) * ∂χ[i, k]
        elseif i != j && k == j
            Γ[i, j, k] = ∂χ[j, i] / χ[j]
        elseif i != j && k == i
            Γ[i, j, k] = ∂χ[i, j] / χ[i]
            # else 0
        end
    end
    return Γ
end

"""
Compute Christoffel symbols from position vector function `OM` via nested AD.
`OM_func(x)` returns a vector of Cartesian coordinates.
"""
function _christoffel_from_OM(OM_func::Function, x::AbstractVector{T}) where {T}
    dim = length(x)
    # Jacobian J[α,i] = ∂OM[α]/∂x[i]
    J = ForwardDiff.jacobian(OM_func, x)
    # Metric function g(y) = J(y)'J(y)
    g_func(y) = let Jy = ForwardDiff.jacobian(OM_func, y)
        Jy' * Jy
    end
    # Reshape flat Jacobian of metric into 3-array:
    #   dg[a,b,c] = ∂g[a,b]/∂x[c]
    # (column-major: vec(g)[(b-1)*dim+a] = g[a,b], so dg_flat[(b-1)*dim+a, c] = ∂g[a,b]/∂x[c])
    dg = reshape(ForwardDiff.jacobian(y -> vec(g_func(y)), x), dim, dim, dim)
    # Metric and its inverse at x
    g = J' * J
    ginv = inv(g)
    Γ = zeros(T, dim, dim, dim)
    for i in 1:dim, j in 1:dim, k in 1:dim
        # Γᵏᵢⱼ = ½ gᵏˡ (∂ᵢgⱼₗ + ∂ⱼgᵢₗ - ∂ₗgᵢⱼ)
        #       = ½ gᵏˡ (dg[l,j,i] + dg[l,i,j] - dg[j,i,l])
        val = zero(T)
        for l in 1:dim
            val += ginv[k, l] * (dg[l, j, i] + dg[l, i, j] - dg[j, i, l])
        end
        Γ[i, j, k] = val / 2
    end
    return Γ
end

# ─────────────────────────────────────────────────────────────────
# Covariant derivative correction: ∇ᵢtʲ… += Γᵢₖʲ tᵏ…
# ─────────────────────────────────────────────────────────────────

"""
Apply the Christoffel correction to a rank-`order` tensor in natural contravariant components.
`t_arr` : array of natural contravariant components (size dim^order)
`Γᵢ`    : Γ[i,:,:] slice (dim×dim matrix, entry [k,j] = Γᵏᵢₖ ... = Γᵏᵢⱼ)
"""
function _apply_christoffel(t_arr::AbstractArray, Γᵢ::AbstractMatrix)
    dim = size(Γᵢ, 1)
    order = ndims(t_arr)
    result = zero(t_arr)
    for o in 1:order
        for idx in CartesianIndices(t_arr)
            j = idx[o]
            for k in 1:dim
                new_idx = CartesianIndex(ntuple(s -> s == o ? k : idx[s], order))
                result[idx] += Γᵢ[k, j] * t_arr[new_idx]
            end
        end
    end
    return result
end

_apply_christoffel(t_arr::Number, ::AbstractMatrix) = zero(t_arr)

# ─────────────────────────────────────────────────────────────────
# χ outer product helper: builds χ^⊗order for natural-frame conversions
# ─────────────────────────────────────────────────────────────────

"""
Return an array of shape (dim, dim, ...) (`order` times) where entry [i₁,…,iₙ] = χᵢ₁ ⋯ χᵢₙ.
"""
function _χ_outer(χ::AbstractVector, order::Integer)
    dim = length(χ)
    result = ones(eltype(χ), ntuple(_ -> dim, order)...)
    for o in 1:order
        sz = ntuple(k -> k == o ? dim : 1, order)
        result = result .* reshape(χ, sz)
    end
    return result
end

# ─────────────────────────────────────────────────────────────────
# Point-wise covariant derivative ∂
# Returns a plain scalar or Array (ForwardDiff-compatible)
# ─────────────────────────────────────────────────────────────────

"""
    ∂(f::Function, i::Integer, CS::CoorSystemNum, x₀)

Covariant derivative of `f` with respect to the `i`-th coordinate, evaluated at `x₀`.

`f` must accept an `AbstractVector` of length `dim` and return either a scalar or a
plain `Array` of **physical (normalized-frame) components** of the tensor field.
Returns a plain scalar or plain `Array` (suitable for further `ForwardDiff` differentiation).

See also `GRAD`, `DIV`, `LAPLACE`, `HESS`, `SYMGRAD` for high-level operators that wrap
results in `AbstractTens`.
"""
function ∂(f::Function, i::Integer, CS::CoorSystemNum{dim, T}, x₀::AbstractVector) where {dim, T}
    # Ensure f_arr returns a plain array (handle AbstractTens-returning functions)
    f_arr = x -> begin
        v = f(x); v isa AbstractTens ? Array(v) : v
    end
    val = f_arr(x₀)

    if val isa Number
        # Scalar: plain partial derivative (no Christoffel correction needed)
        Jac = ForwardDiff.jacobian(x -> [f_arr(x)], collect(x₀))
        return Jac[1, i]
    else
        χ₀ = CS.χ_func(x₀)
        order = ndims(val)
        # Convert to natural contravariant components: t_nat = t_phys / (χᵢ₁ ⋯ χᵢₙ)
        f_nat = x -> f_arr(x) ./ _χ_outer(CS.χ_func(x), order)
        Jac_nat = ForwardDiff.jacobian(x -> vec(f_nat(x)), collect(x₀))
        ∂_nat = reshape(Jac_nat[:, i], size(val))
        # Christoffel correction in natural contravariant frame
        val_nat = val ./ _χ_outer(χ₀, order)
        Γᵢ = CS.Γ_func(x₀)[i, :, :]
        correction_nat = _apply_christoffel(Array(val_nat), Γᵢ)
        # Convert back to physical (normalized) components
        return (∂_nat .+ correction_nat) .* _χ_outer(χ₀, order)
    end
end

# ─────────────────────────────────────────────────────────────────
# Point-wise basis accessors (analogues of CoorSystemSym fields)
# ─────────────────────────────────────────────────────────────────

"""
    normalized_basis(CS::CoorSystemNum, x₀) → AbstractBasis

Return the normalized (orthonormal) basis at point `x₀` as a `RotatedBasis` or
`CanonicalBasis` (when the system is Cartesian at that point).
"""
function normalized_basis(CS::CoorSystemNum{dim, T}, x₀::AbstractVector) where {dim, T}
    R = CS.R_func(x₀)
    ET = eltype(R)
    return RotatedBasis(Matrix{ET}(R))
end

"""
    natural_basis(CS::CoorSystemNum, x₀) → AbstractBasis

Return the natural (non-normalized) basis at point `x₀` as an `OrthogonalBasis`
(i.e. a `RotatedBasis` scaled by the Lamé coefficients).
"""
function natural_basis(CS::CoorSystemNum{dim, T}, x₀::AbstractVector) where {dim, T}
    ℬnorm = normalized_basis(CS, x₀)
    χ = CS.χ_func(x₀)
    ET = eltype(χ)
    return Basis(ℬnorm, Vector{ET}(χ))
end

"""
    natvec(CS::CoorSystemNum, x₀, i, var=:cov) → AbstractTens

Return the `i`-th natural basis vector at point `x₀`, either covariant (`var=:cov`,
proportional to `χᵢ`) or contravariant (`var=:cont`, proportional to `1/χᵢ`).
"""
function natvec(CS::CoorSystemNum{dim, T}, x₀::AbstractVector, i::Integer, var::Symbol = :cov) where {dim, T}
    ℬnorm = normalized_basis(CS, x₀)
    χ = CS.χ_func(x₀)
    ET = eltype(χ)
    if var == :cov
        data = Vec{dim}(j -> j == i ? χ[i] : zero(ET))
        return Tens(data, ℬnorm, (:cov,))
    else  # :cont
        data = Vec{dim}(j -> j == i ? one(ET) / χ[i] : zero(ET))
        return Tens(data, ℬnorm, (:cont,))
    end
end

natvec(CS::CoorSystemNum{dim, T}, x₀::AbstractVector, i::Integer, var::Val{:cov}) where {dim, T} =
    natvec(CS, x₀, i, :cov)
natvec(CS::CoorSystemNum{dim, T}, x₀::AbstractVector, i::Integer, var::Val{:cont}) where {dim, T} =
    natvec(CS, x₀, i, :cont)

"""
    unitvec(CS::CoorSystemNum, x₀, i) → AbstractTens

Return the `i`-th unit vector of the normalized basis at point `x₀`.
"""
function unitvec(CS::CoorSystemNum{dim, T}, x₀::AbstractVector, i::Integer) where {dim, T}
    ℬnorm = normalized_basis(CS, x₀)
    ET = eltype(CS.χ_func(x₀))
    data = Vec{dim}(j -> j == i ? one(ET) : zero(ET))
    return Tens(data, ℬnorm, (:cov,))
end

# ─────────────────────────────────────────────────────────────────
# Point-wise accessors (override AbstractCoorSystem field access)
# ─────────────────────────────────────────────────────────────────

Lame(CS::CoorSystemNum, x::AbstractVector) = CS.χ_func(x)
Christoffel(CS::CoorSystemNum, x::AbstractVector) = CS.Γ_func(x)

# ─────────────────────────────────────────────────────────────────
# Internal raw operators — return plain Arrays (ForwardDiff-safe)
# Used internally by LAPLACE and HESS to avoid creating bases with Dual numbers
# ─────────────────────────────────────────────────────────────────

function _grad_raw(f::Function, CS::CoorSystemNum{dim}) where {dim}
    return function (x₀::AbstractVector)
        f_arr = x -> begin
            v = f(x); v isa AbstractTens ? Array(v) : v
        end
        val = f_arr(x₀)
        χ = CS.χ_func(x₀)
        ET = eltype(χ)
        if val isa Number
            Jac = ForwardDiff.jacobian(x -> [f_arr(x)], collect(x₀))
            return [Jac[1, i] / χ[i] for i in 1:dim]
        else
            order = ndims(val)
            result = zeros(ET, size(val)..., dim)
            for i in 1:dim
                result[fill(:, order)..., i] = ∂(f_arr, i, CS, x₀) ./ χ[i]
            end
            return result
        end
    end
end

function _div_raw(f::Function, CS::CoorSystemNum{dim}) where {dim}
    return function (x₀::AbstractVector)
        f_arr = x -> begin
            v = f(x); v isa AbstractTens ? Array(v) : v
        end
        val = f_arr(x₀)
        χ = CS.χ_func(x₀)
        order = ndims(val)
        result_size = size(val)[1:(end - 1)]
        result = zeros(eltype(χ), result_size...)
        for i in 1:dim
            cov_∂ = ∂(f_arr, i, CS, x₀)
            sl = selectdim(cov_∂, order, i)
            result .+= sl ./ χ[i]
        end
        return ndims(result) == 0 ? result[] : result
    end
end

# ─────────────────────────────────────────────────────────────────
# Public differential operators — return AbstractTens or scalar
# ─────────────────────────────────────────────────────────────────

"""
    GRAD(f::Function, CS::CoorSystemNum) -> Function

Return a function `x₀ -> gradient` as an `AbstractTens` in the normalized basis at `x₀`.
- If `f` returns a scalar, the gradient is a rank-1 `AbstractTens`.
- If `f` returns a rank-`n` `AbstractTens` or array, the gradient is a rank-`n+1` `AbstractTens`.
"""
function GRAD(f::Function, CS::CoorSystemNum{dim}) where {dim}
    return function (x₀::AbstractVector)
        f_arr = x -> begin
            v = f(x); v isa AbstractTens ? Array(v) : v
        end
        val = f_arr(x₀)
        G = _grad_raw(f_arr, CS)(x₀)
        ℬnorm = normalized_basis(CS, x₀)
        if val isa Number
            return Tens(Vec{dim}(i -> G[i]), ℬnorm, (:cont,))
        else
            order = ndims(val)
            return Tens(G, ℬnorm, ntuple(_ -> :cont, order + 1))
        end
    end
end

"""
    SYMGRAD(f::Function, CS::CoorSystemNum) -> Function

Return a function `x₀ -> symmetric gradient` as an `AbstractTens{2}` in the normalized
basis at `x₀`. Applies to vector-valued functions `f`.
"""
function SYMGRAD(f::Function, CS::CoorSystemNum{dim}) where {dim}
    return function (x₀::AbstractVector)
        f_arr = x -> begin
            v = f(x); v isa AbstractTens ? Array(v) : v
        end
        G = _grad_raw(f_arr, CS)(x₀)   # dim×dim plain matrix
        S = (G + permutedims(G, (2, 1))) / 2
        ℬnorm = normalized_basis(CS, x₀)
        return Tens(S, ℬnorm, (:cont, :cont))
    end
end

"""
    DIV(f::Function, CS::CoorSystemNum) -> Function

Return a function `x₀ -> divergence` as a scalar (if `f` is a vector field) or
`AbstractTens` (if `f` is a tensor field of order ≥ 2).

`f` may return a plain `Array` or an `AbstractTens` (components are extracted via `Array`).
"""
function DIV(f::Function, CS::CoorSystemNum{dim}) where {dim}
    return function (x₀::AbstractVector)
        d = _div_raw(f, CS)(x₀)
        if d isa Number
            return d
        else
            ℬnorm = normalized_basis(CS, x₀)
            return Tens(d, ℬnorm, ntuple(_ -> :cont, ndims(d)))
        end
    end
end

"""
    LAPLACE(f::Function, CS::CoorSystemNum) -> Function

Return a function `x₀ -> Laplacian` (scalar). Uses the raw (plain-array) gradient and
divergence internally, so it is fully compatible with `ForwardDiff` differentiation.
"""
function LAPLACE(f::Function, CS::CoorSystemNum{dim}) where {dim}
    f_arr = x -> begin
        v = f(x); v isa AbstractTens ? Array(v) : v
    end
    return x₀ -> _div_raw(_grad_raw(f_arr, CS), CS)(x₀)
end

"""
    HESS(f::Function, CS::CoorSystemNum) -> Function

Return a function `x₀ -> Hessian` as an `AbstractTens{2}` in the normalized basis.
Uses the raw (plain-array) double gradient internally.
"""
function HESS(f::Function, CS::CoorSystemNum{dim}) where {dim}
    f_arr = x -> begin
        v = f(x); v isa AbstractTens ? Array(v) : v
    end
    return function (x₀::AbstractVector)
        H = _grad_raw(_grad_raw(f_arr, CS), CS)(x₀)
        ℬnorm = normalized_basis(CS, x₀)
        return Tens(H, ℬnorm, (:cont, :cont))
    end
end

# ─────────────────────────────────────────────────────────────────
# Predefined factory functions
# ─────────────────────────────────────────────────────────────────

"""
    coorsys_cartesian_num(T=Float64)

Numerical Cartesian coordinate system in 3D. Coordinates: (x, y, z).
"""
function coorsys_cartesian_num(::Type{T} = Float64) where {T <: Real}
    dim = 3
    χ_func = x -> fill(one(x[1]), dim)
    R_func = x -> [i == j ? one(x[1]) : zero(x[1]) for i in 1:dim, j in 1:dim]
    Γ_func = x -> zeros(typeof(x[1]), dim, dim, dim)
    return CoorSystemNum{dim, T}(χ_func, R_func, Γ_func)
end

"""
    coorsys_polar_num(T=Float64)

Numerical polar coordinate system. Coordinates: (r, θ). Lamé coefficients: (1, r).
"""
function coorsys_polar_num(::Type{T} = Float64) where {T <: Real}
    dim = 2
    χ_func = (x::AbstractVector{<:Number}) -> [one(x[1]), x[1]]          # (1, r)
    R_func = (x::AbstractVector{<:Number}) -> begin
        θ = x[2]
        [cos(θ) -sin(θ); sin(θ) cos(θ)]
    end
    Γ_func = (x::AbstractVector{<:Number}) -> begin
        χ = χ_func(x)
        ∂χ = ForwardDiff.jacobian(y -> χ_func(y), collect(x))
        _christoffel_orthogonal(χ, ∂χ)
    end
    return CoorSystemNum{dim, T}(χ_func, R_func, Γ_func)
end

"""
    coorsys_cylindrical_num(T=Float64)

Numerical cylindrical coordinate system. Coordinates: (r, θ, z). Lamé coefficients: (1, r, 1).
"""
function coorsys_cylindrical_num(::Type{T} = Float64) where {T <: Real}
    dim = 3
    χ_func = (x::AbstractVector{<:Number}) -> [one(x[1]), x[1], one(x[1])]  # (1, r, 1)
    R_func = (x::AbstractVector{<:Number}) -> begin
        θ = x[2]
        z = zero(θ); o = one(θ)
        [cos(θ) -sin(θ) z; sin(θ) cos(θ) z; z z o]
    end
    Γ_func = (x::AbstractVector{<:Number}) -> begin
        χ = χ_func(x)
        ∂χ = ForwardDiff.jacobian(y -> χ_func(y), collect(x))
        _christoffel_orthogonal(χ, ∂χ)
    end
    return CoorSystemNum{dim, T}(χ_func, R_func, Γ_func)
end

"""
    coorsys_spherical_num(T=Float64)

Numerical spherical coordinate system. Coordinates: (θ, ϕ, r) — same convention as
`coorsys_spherical()`. Lamé coefficients: (r, r·sin(θ), 1).
"""
function coorsys_spherical_num(::Type{T} = Float64) where {T <: Real}
    dim = 3
    # coords = (θ, ϕ, r)
    χ_func = (x::AbstractVector{<:Number}) -> [x[3], x[3] * sin(x[1]), one(x[1])]
    R_func = (x::AbstractVector{<:Number}) -> begin
        θ, ϕ = x[1], x[2]
        z = θ * 0
        # Columns = 𝐞ᶿ, 𝐞ᵠ, 𝐞ʳ expressed in Cartesian
        [
            cos(θ) * cos(ϕ)  -sin(ϕ)  sin(θ) * cos(ϕ);
            sin(θ) * cos(ϕ)   cos(ϕ)  sin(θ) * sin(ϕ);
            -sin(θ)          z        cos(θ)
        ]
    end
    Γ_func = (x::AbstractVector{<:Number}) -> begin
        χ = χ_func(x)
        ∂χ = ForwardDiff.jacobian(y -> χ_func(y), collect(x))
        _christoffel_orthogonal(χ, ∂χ)
    end
    return CoorSystemNum{dim, T}(χ_func, R_func, Γ_func)
end

"""
    CoorSystemNum(OM_func::Function, dim::Integer, T::Type=Float64)

Generic constructor: build a `CoorSystemNum` from a position-vector function.
`OM_func(x)` must return a `Vector` of Cartesian components given a coordinate vector `x`
of length `dim`. Christoffel symbols are computed by nested automatic differentiation
of the metric tensor `g = J'J` where `J = ∂OM/∂x`.
"""
function CoorSystemNum(OM_func::Function, dim::Integer, ::Type{T} = Float64) where {T <: Real}
    χ_func = x -> begin
        J = ForwardDiff.jacobian(OM_func, collect(x))
        g = J' * J
        [sqrt(g[i, i]) for i in 1:dim]
    end
    R_func = x -> begin
        J = ForwardDiff.jacobian(OM_func, collect(x))
        g = J' * J
        # Normalized columns: eᵢ = aᵢ/χᵢ  where aᵢ = J[:,i]
        hcat([J[:, i] / sqrt(g[i, i]) for i in 1:dim]...)
    end
    Γ_func = x -> _christoffel_from_OM(OM_func, collect(x))
    return CoorSystemNum{dim, T}(χ_func, R_func, Γ_func)
end

export CoorSystemNum
export coorsys_cartesian_num, coorsys_polar_num, coorsys_cylindrical_num, coorsys_spherical_num
