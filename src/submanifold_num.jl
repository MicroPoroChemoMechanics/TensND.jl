# =============================================================================
#  submanifold_num.jl — a hypersurface whose geometry is computed at a point.
#
#  The pointwise counterpart of `submanifold.jl`: the surface is given by a
#  function `OM_func(u)` of its `dim-1` coordinates, and its tangent frame, unit
#  normal, fundamental forms and Gauss–Weingarten coefficients are evaluated at
#  a point `u` by automatic differentiation of that function. It is a
#  `Submanifold` of `PointwiseGeometry`, so that the operators of the pointwise
#  coordinate systems (`coorsystems_num.jl`) apply to it unchanged, along its
#  `dim-1` surface coordinates. The parametrization must be orthogonal
#  (`∂_α OM ⋅ ∂_β OM = 0` for `α ≠ β`), as for `CoordinateSystem(OM_func, dim)`.
# =============================================================================

# The tangent vectors `𝐚_α = ∂_α OM` (columns) and the second derivatives
# `∂_α ∂_β OM` (`H[:, α, β]`) at `u`.
_tangents(OM_func, u) = ForwardDiff.jacobian(OM_func, u)
function _second_derivatives(OM_func, u)
    n = length(u)
    dJ = ForwardDiff.jacobian(v -> vec(ForwardDiff.jacobian(OM_func, v)), u)
    return reshape(dJ, :, n, n)                 # H[i, α, β] = ∂_β ∂_α OMᵢ
end

# The unit normal: the generalized cross product of the normalized tangent
# vectors, `nᵢ = det(𝐞₁, …, 𝐞_{dim-1}, 𝐞ᵢ)`, normalized, as in the exact case.
function _unit_normal(E::AbstractMatrix)
    dim = size(E, 1)
    n = [det(hcat(E, [j == i ? one(eltype(E)) : zero(eltype(E)) for j in 1:dim])) for i in 1:dim]
    return n / norm(n)
end

# Lamé coefficients of the surface, completed by 1 for the normal.
function _surface_lame(OM_func, u)
    J = _tangents(OM_func, u)
    return [[norm(view(J, :, α)) for α in axes(J, 2)]; one(eltype(J))]
end

# The orthonormal frame `(𝐞_α, 𝐧)` as columns, in Cartesian components.
function _surface_frame(OM_func, u)
    J = _tangents(OM_func, u)
    E = reduce(hcat, [view(J, :, α) / norm(view(J, :, α)) for α in axes(J, 2)])
    return hcat(E, _unit_normal(E))
end

# The second fundamental form in natural covariant components,
# `b_αβ = 𝐧 ⋅ ∂_α ∂_β OM`.
function _second_form(OM_func, u)
    n = _surface_frame(OM_func, u)[:, end]
    H = _second_derivatives(OM_func, u)
    return [transpose(n) * view(H, :, α, β) for α in axes(H, 2), β in axes(H, 3)]
end

# The Gauss–Weingarten coefficients `Γ[α, j, k]`, `α` a surface direction, `j`, `k`
# ambient indices, laid out as in the exact case: the intrinsic Christoffel
# symbols `Γᵞ_αβ` of the induced metric, `Γ[α, β, dim] = b_αβ` (Gauss) and
# `Γ[α, dim, γ] = -b_α^γ = -b_αγ / χ_γ²` (Weingarten, the metric being
# diagonal).
function _gauss_weingarten(OM_func, u)
    m = length(u)
    dim = m + 1
    χ = _surface_lame(OM_func, u)
    ∂χ = ForwardDiff.jacobian(v -> _surface_lame(OM_func, v)[1:m], u)
    b = _second_form(OM_func, u)
    T = promote_type(eltype(∂χ), eltype(b))
    Γ = zeros(T, m, dim, dim)
    Γ[:, 1:m, 1:m] .= _christoffel_orthogonal(χ[1:m], ∂χ)
    for α in 1:m, β in 1:m
        Γ[α, β, dim] = b[α, β]
        Γ[α, dim, β] = -b[α, β] / χ[β]^2
    end
    return Γ
end

"""
    Submanifold(OM_func::Function, dim::Integer, T::Type = Float64)

A pointwise hypersurface of the `dim`-dimensional space: `OM_func(u)` returns
the Cartesian components of the point of surface coordinates `u`, a vector of
length `dim-1`. The tangent frame, the unit normal, the fundamental forms and
the Gauss–Weingarten coefficients are computed at a point by automatic
differentiation of `OM_func`; the parametrization must be orthogonal.

The operators of the pointwise coordinate systems apply to it, along the
surface coordinates: [`grad`](@ref) of a function of `u` is the surface
gradient, [`laplacian`](@ref) the Laplace–Beltrami operator. The quantities
of the exact [`Submanifold`](@ref) are read at a point:
[`normal`](@ref)`(SM, u)`, [`submetric`](@ref)`(SM, u)`,
[`curvature`](@ref)`(SM, u)`, [`connection`](@ref)`(SM, u)`.

# Examples
```jldoctest
julia> R = 2.0 ;

julia> Sphere = Submanifold(u -> R * [sin(u[1])cos(u[2]), sin(u[1])sin(u[2]), cos(u[1])], 3) ;

julia> u = [0.7, 0.3] ;

julia> normal(Sphere, u) ≈ Tens([sin(0.7)cos(0.3), sin(0.7)sin(0.3), cos(0.7)])
true

julia> laplacian(v -> cos(v[1]), Sphere)(u) ≈ -2cos(0.7) / R^2    # Laplace–Beltrami
true
```
"""
function Submanifold(OM_func::F, dim::Integer, ::Type{T} = Float64) where {F <: Function, T <: Real}
    dim ≥ 2 || throw(ArgumentError("Submanifold: a hypersurface lives in a space of dimension 2 at least, got $dim"))
    χ_func = u -> _surface_lame(OM_func, collect(u))
    R_func = u -> _surface_frame(OM_func, collect(u))
    Γ_func = u -> _gauss_weingarten(OM_func, collect(u))
    return Submanifold{dim, T}(PointwiseGeometry(χ_func, R_func, Γ_func))
end

# ── The quantities of the surface at a point ─────────────────────────────────

const _PointwiseSubmanifold{dim, T} = Submanifold{dim, T, <:PointwiseGeometry}

normal(SM::_PointwiseSubmanifold{dim}, u::AbstractVector) where {dim} = unitvec(SM, u, dim)

# A form of the tangent block, stored in `dim` dimensions with a vanishing last
# row and column, in the natural basis at `u`, as in the exact case.
function _tangent_form(SM::_PointwiseSubmanifold{dim}, u, f) where {dim}
    S = SymmetricTensor{2, dim}((i, j) -> i < dim && j < dim ? f(i, j) : zero(f(1, 1)))
    return Tens(S, natural_basis(SM, u), (:cov, :cov))
end

"""
    submetric(SM::Submanifold, u) → AbstractTens{2,dim}

The first fundamental form of a pointwise submanifold at the point `u`; see the
exact method.
"""
function submetric(SM::_PointwiseSubmanifold, u::AbstractVector)
    χ = lame_coefficients(SM, u)
    return _tangent_form(SM, u, (α, β) -> α == β ? χ[α]^2 : zero(χ[α]))
end

"""
    curvature(SM::Submanifold, u) → AbstractTens{2,dim}

The second fundamental form of a pointwise submanifold at the point `u`; see
the exact method.
"""
function curvature(SM::_PointwiseSubmanifold, u::AbstractVector)
    b = christoffel(SM, u)
    return _tangent_form(SM, u, (α, β) -> b[α, β, end])
end

"""
    connection(SM::Submanifold, u) → Array{T,3}

The intrinsic connection coefficients `Γᵞ_αβ` of a pointwise submanifold at the
point `u`; see the exact method.
"""
connection(SM::_PointwiseSubmanifold{dim}, u::AbstractVector) where {dim} =
    christoffel(SM, u)[:, 1:(dim - 1), 1:(dim - 1)]
