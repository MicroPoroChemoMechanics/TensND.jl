using ForwardDiff

"""
    PointwiseGeometry{X,R,G} <: AbstractGeometry

The geometry of a chart as three functions of the point `x` (an `AbstractVector`
of coordinates), for differential operators by automatic differentiation:

- `χ_func(x)` : the Lamé coefficients at `x`;
- `R_func(x)` : the rotation matrix whose columns are the unit vectors of the
  normalized basis at `x`, in Cartesian components;
- `Γ_func(x)` : the Christoffel symbols at `x`, `Γ[i,j,k] = Γᵏᵢⱼ`.

A [`CoordinateSystem`](@ref) of this geometry takes fields as functions of the
point and returns functions of the point: `grad(f, CS)` is `x₀ -> ∇f(x₀)`. The
bases and the Lamé coefficients depend on the point too, and are read with it:
`normalized_basis(CS, x₀)`, `lame_coefficients(CS, x₀)`, ….

Predefined: `coorsys_cartesian(PointwiseGeometry)`, `coorsys_polar(PointwiseGeometry)`,
`coorsys_cylindrical(PointwiseGeometry)`, `coorsys_spherical(PointwiseGeometry)`,
each with an optional element type (`Float64` by default); any chart from its
position vector with `CoordinateSystem(OM_func, dim, T)`.
"""
struct PointwiseGeometry{X, R, G} <: AbstractGeometry
    # The three functions are type parameters rather than `::Function` fields,
    # so that what they return is inferred.
    χ_func::X   # x -> Vector of Lamé coefficients
    R_func::R   # x -> Matrix  (columns = unit vectors eᵢ in Cartesian frame)
    Γ_func::G   # x -> Array{T,3}  Γ[i,j,k] = Γᵏᵢⱼ
end

# The pointwise systems, and how a dimension and an element type make one.
const _PointwiseSystem{dim, T} = CoordinateSystem{dim, T, <:PointwiseGeometry}
CoordinateSystem{dim, T}(g::PointwiseGeometry) where {dim, T <: Real} =
    CoordinateSystem{dim, T, typeof(g)}(g)
_pointwise_system(::Val{dim}, ::Type{T}, χ_func, R_func, Γ_func) where {dim, T} =
    CoordinateSystem{dim, T}(PointwiseGeometry(χ_func, R_func, Γ_func))

# The three functions of a pointwise system.
_χ_func(CS::_PointwiseSystem) = getfield(CS, :geometry).χ_func
_R_func(CS::_PointwiseSystem) = getfield(CS, :geometry).R_func
_Γ_func(CS::_PointwiseSystem) = getfield(CS, :geometry).Γ_func

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

# A field as a function of plain values: the components of a tensor, a number
# or an array as it is.
_plain_value(v::AbstractTens) = Array(v)
_plain_value(v) = v
_plain_field(f) = x -> _plain_value(f(x))

"""
    ∂(f::Function, i::Integer, CS::_PointwiseSystem, x₀)

Covariant derivative of `f` with respect to the `i`-th coordinate, evaluated at `x₀`.

`f` must accept an `AbstractVector` of length `dim` and return either a scalar or a
plain `Array` of **physical (normalized-frame) components** of the tensor field.
Returns a plain scalar or plain `Array` (suitable for further `ForwardDiff` differentiation).

See also `grad`, `divergence`, `laplacian`, `hess`, `symgrad` for high-level operators that wrap
results in `AbstractTens`.
"""
function ∂(f::Function, i::Integer, CS::_PointwiseSystem{dim, T}, x₀::AbstractVector) where {dim, T}
    f_arr = _plain_field(f)
    return _∂(f_arr(x₀), f_arr, i, CS, x₀)
end

# A scalar field: the plain partial derivative, no Christoffel correction.
function _∂(::Number, f_arr, i, CS, x₀)
    Jac = ForwardDiff.jacobian(x -> [f_arr(x)], collect(x₀))
    return Jac[1, i]
end

function _∂(val::AbstractArray, f_arr, i, CS, x₀)
    χ₀ = _χ_func(CS)(x₀)
    order = ndims(val)
    # Convert to natural contravariant components: t_nat = t_phys / (χᵢ₁ ⋯ χᵢₙ)
    f_nat = x -> f_arr(x) ./ _χ_outer(_χ_func(CS)(x), order)
    Jac_nat = ForwardDiff.jacobian(x -> vec(f_nat(x)), collect(x₀))
    ∂_nat = reshape(Jac_nat[:, i], size(val))
    # Christoffel correction in natural contravariant frame
    val_nat = val ./ _χ_outer(χ₀, order)
    Γᵢ = _Γ_func(CS)(x₀)[i, :, :]
    correction_nat = _apply_christoffel(Array(val_nat), Γᵢ)
    # Convert back to physical (normalized) components
    return (∂_nat .+ correction_nat) .* _χ_outer(χ₀, order)
end

# ─────────────────────────────────────────────────────────────────
# Point-wise basis accessors (the exact geometry stores them as fields)
# ─────────────────────────────────────────────────────────────────

"""
    normalized_basis(CS::_PointwiseSystem, x₀) → AbstractBasis

Return the normalized (orthonormal) basis at point `x₀` as a `RotatedBasis` or
`CanonicalBasis` (when the system is Cartesian at that point).
"""
function normalized_basis(CS::_PointwiseSystem{dim, T}, x₀::AbstractVector) where {dim, T}
    R = _R_func(CS)(x₀)
    ET = eltype(R)
    return RotatedBasis(Matrix{ET}(R))
end

"""
    natural_basis(CS::_PointwiseSystem, x₀) → AbstractBasis

Return the natural (non-normalized) basis at point `x₀` as an `OrthogonalBasis`
(i.e. a `RotatedBasis` scaled by the Lamé coefficients).
"""
function natural_basis(CS::_PointwiseSystem{dim, T}, x₀::AbstractVector) where {dim, T}
    ℬnorm = normalized_basis(CS, x₀)
    χ = _χ_func(CS)(x₀)
    ET = eltype(χ)
    return Basis(ℬnorm, Vector{ET}(χ))
end

"""
    natvec(CS::_PointwiseSystem, x₀, i, var=:cov) → AbstractTens

Return the `i`-th natural basis vector at point `x₀`, either covariant (`var=:cov`,
proportional to `χᵢ`) or contravariant (`var=:cont`, proportional to `1/χᵢ`).
"""
natvec(CS::_PointwiseSystem, x₀::AbstractVector, i::Integer, var::Symbol = :cov) =
    _by_variance(natvec, var, CS, x₀, i)

function natvec(CS::_PointwiseSystem{dim}, x₀::AbstractVector, i::Integer, ::Val{:cov}) where {dim}
    χ = _χ_func(CS)(x₀)
    ET = eltype(χ)
    data = Vec{dim}(j -> j == i ? χ[i] : zero(ET))
    return Tens(data, normalized_basis(CS, x₀), (:cov,))
end

function natvec(CS::_PointwiseSystem{dim}, x₀::AbstractVector, i::Integer, ::Val{:cont}) where {dim}
    χ = _χ_func(CS)(x₀)
    ET = eltype(χ)
    data = Vec{dim}(j -> j == i ? one(ET) / χ[i] : zero(ET))
    return Tens(data, normalized_basis(CS, x₀), (:cont,))
end

"""
    unitvec(CS::_PointwiseSystem, x₀, i) → AbstractTens

Return the `i`-th unit vector of the normalized basis at point `x₀`.
"""
function unitvec(CS::_PointwiseSystem{dim, T}, x₀::AbstractVector, i::Integer) where {dim, T}
    ℬnorm = normalized_basis(CS, x₀)
    ET = eltype(_χ_func(CS)(x₀))
    data = Vec{dim}(j -> j == i ? one(ET) : zero(ET))
    return Tens(data, ℬnorm, (:cov,))
end

# ─────────────────────────────────────────────────────────────────
# Point-wise accessors (override AbstractCoorSystem field access)
# ─────────────────────────────────────────────────────────────────

lame_coefficients(CS::_PointwiseSystem, x::AbstractVector) = _χ_func(CS)(x)
christoffel(CS::_PointwiseSystem, x::AbstractVector) = _Γ_func(CS)(x)

# ─────────────────────────────────────────────────────────────────
# Internal raw operators — return plain Arrays (ForwardDiff-safe)
# Used internally by laplacian and hess to avoid creating bases with Dual numbers
# ─────────────────────────────────────────────────────────────────

function _grad_raw(f::Function, CS::_PointwiseSystem{dim}) where {dim}
    return function (x₀::AbstractVector)
        f_arr = _plain_field(f)
        val = f_arr(x₀)
        χ = _χ_func(CS)(x₀)
        # The buffer must hold both the geometry (χ, from the coordinate system)
        # and the field value. Taking `eltype(χ)` alone makes the operators
        # undifferentiable with respect to a parameter carried by the *field*:
        # the system is built at a `Float64` point, so a `Dual`-valued field had
        # nowhere to go and the call failed converting a `Dual` to a `Float64`.
        # Promoting costs nothing when both are `Float64`.
        return _grad_at(val, f_arr, CS, χ, promote_type(eltype(χ), eltype(val)), x₀)
    end
end

function _grad_at(::Number, f_arr, CS::_PointwiseSystem{dim}, χ, ET, x₀) where {dim}
    Jac = ForwardDiff.jacobian(x -> [f_arr(x)], collect(x₀))
    return [Jac[1, i] / χ[i] for i in 1:dim]
end

function _grad_at(val::AbstractArray, f_arr, CS::_PointwiseSystem{dim}, χ, ET, x₀) where {dim}
    order = ndims(val)
    result = zeros(ET, size(val)..., dim)
    for i in 1:dim
        result[fill(:, order)..., i] = ∂(f_arr, i, CS, x₀) ./ χ[i]
    end
    return result
end

function _div_raw(f::Function, CS::_PointwiseSystem{dim}) where {dim}
    return function (x₀::AbstractVector)
        f_arr = _plain_field(f)
        val = f_arr(x₀)
        χ = _χ_func(CS)(x₀)
        order = ndims(val)
        result_size = size(val)[1:(end - 1)]
        # Same promotion as in `_grad_raw`: the buffer holds geometry *and*
        # field, so a field carrying `Dual` numbers stays differentiable.
        result = zeros(promote_type(eltype(χ), eltype(val)), result_size...)
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
    grad(f::Function, CS::_PointwiseSystem) -> Function

Return a function `x₀ -> gradient` as an `AbstractTens` in the normalized basis at `x₀`.
- If `f` returns a scalar, the gradient is a rank-1 `AbstractTens`.
- If `f` returns a rank-`n` `AbstractTens` or array, the gradient is a rank-`n+1` `AbstractTens`.
"""
function grad(f::Function, CS::_PointwiseSystem{dim}) where {dim}
    return function (x₀::AbstractVector)
        f_arr = _plain_field(f)
        val = f_arr(x₀)
        G = _grad_raw(f_arr, CS)(x₀)
        return _grad_tens(val, G, normalized_basis(CS, x₀))
    end
end

# The gradient of a scalar field is a vector; of an order-n field, order n + 1.
_grad_tens(::Number, G, ℬ::AbstractBasis{dim}) where {dim} = Tens(Vec{dim}(i -> G[i]), ℬ, (:cont,))
_grad_tens(val::AbstractArray, G, ℬ) = Tens(G, ℬ, ntuple(_ -> :cont, ndims(val) + 1))

"""
    symgrad(f::Function, CS::_PointwiseSystem) -> Function

Return a function `x₀ -> symmetric gradient` as an `AbstractTens{2}` in the normalized
basis at `x₀`. Applies to vector-valued functions `f`.
"""
function symgrad(f::Function, CS::_PointwiseSystem{dim}) where {dim}
    return function (x₀::AbstractVector)
        f_arr = _plain_field(f)
        G = _grad_raw(f_arr, CS)(x₀)   # dim×dim plain matrix
        S = (G + permutedims(G, (2, 1))) / 2
        ℬnorm = normalized_basis(CS, x₀)
        return Tens(S, ℬnorm, (:cont, :cont))
    end
end

"""
    divergence(f::Function, CS::_PointwiseSystem) -> Function

Return a function `x₀ -> divergence` as a scalar (if `f` is a vector field) or
`AbstractTens` (if `f` is a tensor field of order ≥ 2).

`f` may return a plain `Array` or an `AbstractTens` (components are extracted via `Array`).
"""
function divergence(f::Function, CS::_PointwiseSystem{dim}) where {dim}
    return function (x₀::AbstractVector)
        return _div_tens(_div_raw(f, CS)(x₀), CS, x₀)
    end
end

# The divergence of a vector field is a number; of an order-n field, order n − 1.
_div_tens(d::Number, CS, x₀) = d
_div_tens(d::AbstractArray, CS, x₀) = Tens(d, normalized_basis(CS, x₀), ntuple(_ -> :cont, ndims(d)))

"""
    laplacian(f::Function, CS::_PointwiseSystem) -> Function

Return a function `x₀ -> Laplacian` (scalar). Uses the raw (plain-array) gradient and
divergence internally, so it is fully compatible with `ForwardDiff` differentiation.
"""
function laplacian(f::Function, CS::_PointwiseSystem{dim}) where {dim}
    f_arr = _plain_field(f)
    return x₀ -> _div_raw(_grad_raw(f_arr, CS), CS)(x₀)
end

"""
    hess(f::Function, CS::_PointwiseSystem) -> Function

Return a function `x₀ -> Hessian` as an `AbstractTens{2}` in the normalized basis.
Uses the raw (plain-array) double gradient internally.
"""
function hess(f::Function, CS::_PointwiseSystem{dim}) where {dim}
    f_arr = _plain_field(f)
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
    coorsys_cartesian(PointwiseGeometry, T = Float64)

Numerical Cartesian coordinate system in 3D. Coordinates: (x, y, z).
"""
function coorsys_cartesian(::Type{PointwiseGeometry}, ::Type{T} = Float64) where {T <: Real}
    dim = 3
    χ_func = x -> fill(one(x[1]), dim)
    R_func = x -> [i == j ? one(x[1]) : zero(x[1]) for i in 1:dim, j in 1:dim]
    Γ_func = x -> zeros(typeof(x[1]), dim, dim, dim)
    return _PointwiseSystem{dim, T}(χ_func, R_func, Γ_func)
end

"""
    coorsys_polar(PointwiseGeometry, T = Float64)

Numerical polar coordinate system. Coordinates: (r, θ). Lamé coefficients: (1, r).
"""
function coorsys_polar(::Type{PointwiseGeometry}, ::Type{T} = Float64) where {T <: Real}
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
    return _PointwiseSystem{dim, T}(χ_func, R_func, Γ_func)
end

"""
    coorsys_cylindrical(PointwiseGeometry, T = Float64)

Numerical cylindrical coordinate system. Coordinates: (r, θ, z). Lamé coefficients: (1, r, 1).
"""
function coorsys_cylindrical(::Type{PointwiseGeometry}, ::Type{T} = Float64) where {T <: Real}
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
    return _PointwiseSystem{dim, T}(χ_func, R_func, Γ_func)
end

"""
    coorsys_spherical(PointwiseGeometry, T = Float64)

Numerical spherical coordinate system. Coordinates: (θ, ϕ, r) — same convention as
`coorsys_spherical()`. Lamé coefficients: (r, r·sin(θ), 1).
"""
function coorsys_spherical(::Type{PointwiseGeometry}, ::Type{T} = Float64) where {T <: Real}
    dim = 3
    # coords = (θ, ϕ, r)
    χ_func = (x::AbstractVector{<:Number}) -> [x[3], x[3] * sin(x[1]), one(x[1])]
    R_func = (x::AbstractVector{<:Number}) -> begin
        θ, ϕ = x[1], x[2]
        z = θ * 0
        # Columns = 𝐞ᶿ, 𝐞ᵠ, 𝐞ʳ expressed in Cartesian. 𝐞ᶿ = (cos θ cos ϕ,
        # cos θ sin ϕ, −sin θ): its second entry read sin θ cos ϕ, so the
        # "normalized" basis was not even normed (|𝐞ᶿ|² = cos²ϕ + sin²θ).
        [
            cos(θ) * cos(ϕ)  -sin(ϕ)  sin(θ) * cos(ϕ);
            cos(θ) * sin(ϕ)   cos(ϕ)  sin(θ) * sin(ϕ);
            -sin(θ)          z        cos(θ)
        ]
    end
    Γ_func = (x::AbstractVector{<:Number}) -> begin
        χ = χ_func(x)
        ∂χ = ForwardDiff.jacobian(y -> χ_func(y), collect(x))
        _christoffel_orthogonal(χ, ∂χ)
    end
    return _PointwiseSystem{dim, T}(χ_func, R_func, Γ_func)
end

"""
    CoordinateSystem(OM_func::Function, dim::Integer, T::Type = Float64)

A pointwise system from its position vector: `OM_func(x)` returns the Cartesian
components of the point of coordinates `x` (a vector of length `dim`). The Lamé
coefficients and the normalized basis come from the Jacobian of `OM_func`, the
Christoffel symbols from the derivatives of the metric `g = JᵀJ`, all by
automatic differentiation.
"""
function CoordinateSystem(OM_func::Function, dim::Integer, ::Type{T} = Float64) where {T <: Real}
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
    return CoordinateSystem{dim, T}(PointwiseGeometry(χ_func, R_func, Γ_func))
end
