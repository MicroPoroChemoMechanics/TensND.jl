# ── Charts: a curvilinear system as a map between coordinates and points ─────
#
# A chart holds the parameters of a coordinate system (a semi-focal distance,
# none for the spherical one) and answers, at a point, with numbers: the
# Cartesian components of the point of given coordinates, the coordinates of a
# given point, the Lamé coefficients and the unit vectors. The same closed
# forms, written once here, make the exact systems (`coorsys_spheroidal`,
# `coorsys_ellipsoidal`, …, called with symbols) and the pointwise ones.

"""
    AbstractChart{dim}

A curvilinear coordinate system of dimension `dim` as a map: the Cartesian
components of the point of coordinates `ξ` are `coords_to_cartesian(chart, ξ)`,
and the coordinates of the point `x` are `cartesian_to_coords(chart, x)`.

Predefined: [`PolarChart`](@ref), [`CylindricalChart`](@ref),
[`SphericalChart`](@ref), [`ProlateSpheroidalChart`](@ref),
[`OblateSpheroidalChart`](@ref), [`EllipsoidalChart`](@ref).
"""
abstract type AbstractChart{dim} end

"""
    PolarChart()

Polar coordinates `(r, θ)`: `x₁ = r cos θ`, `x₂ = r sin θ`.
"""
struct PolarChart <: AbstractChart{2} end

"""
    CylindricalChart()

Cylindrical coordinates `(r, θ, z)`: `x₁ = r cos θ`, `x₂ = r sin θ`, `x₃ = z`.
"""
struct CylindricalChart <: AbstractChart{3} end

"""
    SphericalChart()

Spherical coordinates in the order `(θ, ϕ, r)` of [`coorsys_spherical`](@ref):
`x₁ = r sin θ cos ϕ`, `x₂ = r sin θ sin ϕ`, `x₃ = r cos θ`.
"""
struct SphericalChart <: AbstractChart{3} end

# A parameter is checked where a comparison means something: a symbol is the
# caller's responsibility, as for the bases.
_check_positive(x, name) = _check_positive(x, name, Val(is_hard_numeric(typeof(x))))
_check_positive(x, name, ::Val{false}) = nothing
function _check_positive(x, name, ::Val{true})
    x > zero(x) || throw(DomainError(x, "the $name must be positive"))
    return nothing
end

"""
    ProlateSpheroidalChart(c)

Prolate spheroidal coordinates `(ϕ, p, q)` of axis `𝐞₃` and semi-focal distance
`c > 0`, with `-1 ≤ p ≤ 1` and `q ≥ 1`:

    x₁ = c √(1-p²) √(q²-1) cos ϕ,   x₂ = c √(1-p²) √(q²-1) sin ϕ,   x₃ = c p q.

The surfaces `q = const` are the confocal prolate spheroids of semi-axes `c q`
along `𝐞₃` and `c √(q²-1)` across; the foci are the points `± c 𝐞₃`.

See also [`OblateSpheroidalChart`](@ref), [`coorsys_spheroidal`](@ref).
"""
struct ProlateSpheroidalChart{T <: Number} <: AbstractChart{3}
    c::T
    function ProlateSpheroidalChart(c::T) where {T <: Number}
        _check_positive(c, "semi-focal distance")
        return new{T}(c)
    end
end

"""
    OblateSpheroidalChart(c)

Oblate spheroidal coordinates `(ϕ, p, τ)` of axis `𝐞₃` and semi-focal distance
`c > 0`, with `-1 ≤ p ≤ 1` and `τ ≥ 0`:

    x₁ = c √(1-p²) √(τ²+1) cos ϕ,   x₂ = c √(1-p²) √(τ²+1) sin ϕ,   x₃ = c p τ.

The surfaces `τ = const` are the confocal oblate spheroids of semi-axes `c τ`
along `𝐞₃` and `c √(τ²+1)` across; the focal circle has radius `c` in the plane
`x₃ = 0`. The chart is the prolate one continued to `c → -i c`, `q → i τ`.

See also [`ProlateSpheroidalChart`](@ref), [`coorsys_spheroidal_oblate`](@ref).
"""
struct OblateSpheroidalChart{T <: Number} <: AbstractChart{3}
    c::T
    function OblateSpheroidalChart(c::T) where {T <: Number}
        _check_positive(c, "semi-focal distance")
        return new{T}(c)
    end
end

const _SpheroidalChart{T} = Union{ProlateSpheroidalChart{T}, OblateSpheroidalChart{T}}

# `q² - s` and `q² - s p²` with `s = 1` (prolate) or `s = -1` (oblate): the two
# charts share every formula up to this sign.
_spheroidal_sign(::ProlateSpheroidalChart) = 1
_spheroidal_sign(::OblateSpheroidalChart) = -1

"""
    EllipsoidalChart(h₂, h₃)

Ellipsoidal coordinates `(ρ, μ, ν)` of semi-focal distances `h₂ > h₃ > 0`, with
`h₁² = h₂² - h₃²` and `ν² ≤ h₃² ≤ μ² ≤ h₂² ≤ ρ²`:

    x₁ = ρ μ ν / (h₂ h₃),
    x₂ = √(ρ²-h₃²) √(μ²-h₃²) √(h₃²-ν²) / (h₁ h₃),
    x₃ = √(ρ²-h₂²) √(h₂²-μ²) √(h₂²-ν²) / (h₁ h₂).

The surfaces `ρ = const` are the confocal ellipsoids of semi-axes `ρ`,
`√(ρ²-h₃²)`, `√(ρ²-h₂²)` along `𝐞₁`, `𝐞₂`, `𝐞₃`; `μ = const` are hyperboloids of
one sheet and `ν = const` hyperboloids of two sheets. For an ellipsoid of
semi-axes `a₁ > a₂ > a₃`, `h₂² = a₁² - a₃²` and `h₃² = a₁² - a₂²`, and its
surface is `ρ = a₁`.

With the positive square roots, the chart covers the quarter space `x₂ ≥ 0`,
`x₃ ≥ 0`, the sign of `ν` being that of `x₁`. The unit vectors
`(𝐞ᵨ, 𝐞ᵤ, 𝐞ᵥ)`, taken in the order of the coordinates, form a left-handed
triad; `(𝐞ᵨ, 𝐞ᵥ, 𝐞ᵤ)` is right-handed.

See also [`coorsys_ellipsoidal`](@ref).
"""
struct EllipsoidalChart{T <: Number} <: AbstractChart{3}
    h₂::T
    h₃::T
    function EllipsoidalChart(h₂::T, h₃::T) where {T <: Number}
        _check_focal_order(h₂, h₃, Val(is_hard_numeric(T)))
        return new{T}(h₂, h₃)
    end
end
EllipsoidalChart(h₂::Number, h₃::Number) = EllipsoidalChart(promote(h₂, h₃)...)

_check_focal_order(h₂, h₃, ::Val{false}) = nothing
function _check_focal_order(h₂, h₃, ::Val{true})
    h₂ > h₃ > zero(h₃) || throw(
        DomainError(
            (h₂, h₃),
            "the semi-focal distances must satisfy h₂ > h₃ > 0; for a spheroid " *
                "(h₃ = 0 or h₃ = h₂) use a spheroidal chart",
        )
    )
    return nothing
end

"""
    EllipsoidalChart(a₁, a₂, a₃)

The ellipsoidal chart whose ellipsoid `ρ = a₁` has the semi-axes
`a₁ > a₂ > a₃ > 0` along `𝐞₁`, `𝐞₂`, `𝐞₃`.
"""
EllipsoidalChart(a₁::Number, a₂::Number, a₃::Number) =
    EllipsoidalChart(sqrt(a₁^2 - a₃^2), sqrt(a₁^2 - a₂^2))

_three(ξ) = (ξ[1], ξ[2], ξ[3])

# ── Coordinates → point ──────────────────────────────────────────────────────

"""
    coords_to_cartesian(chart::AbstractChart, ξ) -> SVector

Cartesian components of the point whose coordinates in `chart` are `ξ` (a
tuple or a vector). The result carries the element type of the coordinates and
of the parameters of the chart together: duals and symbols pass through.

# Examples
```jldoctest
julia> coords_to_cartesian(ProlateSpheroidalChart(1.0), (0.0, 0.0, 2.0)) ≈ [√3, 0, 0]
true
```

See also [`cartesian_to_coords`](@ref), [`change_coords`](@ref).
"""
function coords_to_cartesian(::PolarChart, ξ)
    r, θ = ξ[1], ξ[2]
    return SVector(r * cos(θ), r * sin(θ))
end
function coords_to_cartesian(::CylindricalChart, ξ)
    r, θ, z = _three(ξ)
    return SVector(promote(r * cos(θ), r * sin(θ), z)...)
end
function coords_to_cartesian(::SphericalChart, ξ)
    θ, ϕ, r = _three(ξ)
    return SVector(r * sin(θ) * cos(ϕ), r * sin(θ) * sin(ϕ), r * cos(θ))
end
function coords_to_cartesian(chart::_SpheroidalChart, ξ)
    ϕ, p, q = _three(ξ)
    ϱ = chart.c * sqrt(1 - p^2) * sqrt(q^2 - _spheroidal_sign(chart))
    return SVector(promote(ϱ * cos(ϕ), ϱ * sin(ϕ), chart.c * p * q)...)
end
function coords_to_cartesian(chart::EllipsoidalChart, ξ)
    ρ, μ, ν = _three(ξ)
    h₂, h₃ = chart.h₂, chart.h₃
    h₁ = sqrt(h₂^2 - h₃^2)
    return SVector(
        promote(
            ρ * μ * ν / (h₂ * h₃),
            sqrt(ρ^2 - h₃^2) * sqrt(μ^2 - h₃^2) * sqrt(h₃^2 - ν^2) / (h₁ * h₃),
            sqrt(ρ^2 - h₂^2) * sqrt(h₂^2 - μ^2) * sqrt(h₂^2 - ν^2) / (h₁ * h₂),
        )...
    )
end

# ── Lamé coefficients and unit vectors ───────────────────────────────────────

# `χᵢ = ‖∂𝐱/∂ξᵢ‖`, in the order of the coordinates.
function _lame(chart::_SpheroidalChart, ξ)
    _, p, q = _three(ξ)
    c, s = chart.c, _spheroidal_sign(chart)
    return (
        c * sqrt(1 - p^2) * sqrt(q^2 - s),
        c * sqrt(q^2 - s * p^2) / sqrt(1 - p^2),
        c * sqrt(q^2 - s * p^2) / sqrt(q^2 - s),
    )
end
function _lame(chart::EllipsoidalChart, ξ)
    ρ, μ, ν = _three(ξ)
    h₂, h₃ = chart.h₂, chart.h₃
    return (
        sqrt(ρ^2 - μ^2) * sqrt(ρ^2 - ν^2) / (sqrt(ρ^2 - h₃^2) * sqrt(ρ^2 - h₂^2)),
        sqrt(ρ^2 - μ^2) * sqrt(μ^2 - ν^2) / (sqrt(μ^2 - h₃^2) * sqrt(h₂^2 - μ^2)),
        sqrt(ρ^2 - ν^2) * sqrt(μ^2 - ν^2) / (sqrt(h₃^2 - ν^2) * sqrt(h₂^2 - ν^2)),
    )
end

# The unit vectors `𝐞ᵢ = (∂𝐱/∂ξᵢ)/χᵢ` as the columns of a matrix, in Cartesian
# components.
function _frame(chart::_SpheroidalChart, ξ)
    ϕ, p, q = _three(ξ)
    s = _spheroidal_sign(chart)
    p̄, q̄, d = sqrt(1 - p^2), sqrt(q^2 - s), sqrt(q^2 - s * p^2)
    T = typeof(q̄ / d)
    return SMatrix{3, 3, T}(
        -sin(ϕ), cos(ϕ), zero(T),
        -p * q̄ * cos(ϕ) / d, -p * q̄ * sin(ϕ) / d, q * p̄ / d,
        q * p̄ * cos(ϕ) / d, q * p̄ * sin(ϕ) / d, p * q̄ / d,
    )
end
function _frame(chart::EllipsoidalChart, ξ)
    ρ, μ, ν = _three(ξ)
    h₂, h₃ = chart.h₂, chart.h₃
    x₁, x₂, x₃ = coords_to_cartesian(chart, ξ)
    χᵨ, χᵤ, χᵥ = _lame(chart, ξ)
    # `∂xᵢ/∂ξ = xᵢ ∂ln xᵢ/∂ξ`, each coordinate entering one factor of each `xᵢ`.
    # Entry by entry: Symbolics takes an `hcat` of its own numbers as a
    # symbolic array.
    T = typeof(x₁ / χᵨ)
    return SMatrix{3, 3, T}(
        x₁ / (ρ * χᵨ), x₂ * ρ / ((ρ^2 - h₃^2) * χᵨ), x₃ * ρ / ((ρ^2 - h₂^2) * χᵨ),
        x₁ / (μ * χᵤ), x₂ * μ / ((μ^2 - h₃^2) * χᵤ), -x₃ * μ / ((h₂^2 - μ^2) * χᵤ),
        x₁ / (ν * χᵥ), -x₂ * ν / ((h₃^2 - ν^2) * χᵥ), -x₃ * ν / ((h₂^2 - ν^2) * χᵥ),
    )
end

# ── Point → coordinates ──────────────────────────────────────────────────────
#
# Every inversion below branches on values (a sign, an ordering of roots): it
# takes numbers, duals included, not symbols.

"""
    cartesian_to_coords(chart::AbstractChart, x) -> SVector

Coordinates in `chart` of the point of Cartesian components `x`, the inverse of
[`coords_to_cartesian`](@ref). The inversion is in closed form and
differentiable (ForwardDiff duals in `x` and in the parameters of the chart),
away from the singular set of the chart (axis, foci, focal curves).

For an [`EllipsoidalChart`](@ref), the coordinates returned are those of the
mirror image of `x` in the quarter space `x₂ ≥ 0`, `x₃ ≥ 0` that the chart
covers: `ν` has the sign of `x₁`, and
`coords_to_cartesian(chart, cartesian_to_coords(chart, x)) ≈ (x₁, |x₂|, |x₃|)`.

# Examples
```jldoctest
julia> chart = OblateSpheroidalChart(2.0); ξ = (0.3, -0.4, 1.5);

julia> cartesian_to_coords(chart, coords_to_cartesian(chart, ξ)) ≈ [ξ...]
true
```
"""
function cartesian_to_coords(::PolarChart, x)
    return SVector(hypot(x[1], x[2]), atan(x[2], x[1]))
end
function cartesian_to_coords(::CylindricalChart, x)
    return SVector(promote(hypot(x[1], x[2]), atan(x[2], x[1]), x[3])...)
end
function cartesian_to_coords(::SphericalChart, x)
    ϱ = hypot(x[1], x[2])
    return SVector(atan(ϱ, x[3]), atan(x[2], x[1]), hypot(ϱ, x[3]))
end
function cartesian_to_coords(chart::ProlateSpheroidalChart, x)
    c = chart.c
    ϱ = hypot(x[1], x[2])
    # `2 c q` is the sum of the distances to the two foci; `p` follows from
    # `x₃ = c p q`, which avoids the difference of the two distances.
    q = (hypot(ϱ, x[3] - c) + hypot(ϱ, x[3] + c)) / (2c)
    # On the axis `|x₃| = c q` up to rounding: `p` is kept in its interval.
    p = clamp(x[3] / (c * q), -one(q), one(q))
    return SVector(promote(atan(x[2], x[1]), p, q)...)
end
function cartesian_to_coords(chart::OblateSpheroidalChart, x)
    c = chart.c
    ϱ² = x[1]^2 + x[2]^2
    z = x[3]
    # `τ²` and `-p²` are the two roots of `c² u² - d u - z² = 0`. The root that
    # adds two terms of the same sign is formed first, the other follows from
    # the product of the roots: no difference of close numbers.
    d = ϱ² + z^2 - c^2
    Δ = sqrt(d^2 + 4 * c^2 * z^2)
    large = (abs(d) + Δ) / (2 * c^2)
    # On the focal circle both roots vanish.
    small = iszero(large) ? large : z^2 / (c^2 * large)
    τ², p² = d ≥ zero(d) ? (large, small) : (small, large)
    p = z < zero(z) ? -sqrt(p²) : sqrt(p²)
    p = clamp(p, -one(p), one(p))            # `|p| = 1` on the axis, up to rounding
    return SVector(promote(atan(x[2], x[1]), p, sqrt(τ²))...)
end
function cartesian_to_coords(chart::EllipsoidalChart, x)
    h₂², h₃² = chart.h₂^2, chart.h₃^2
    x₁², x₂², x₃² = x[1]^2, x[2]^2, x[3]^2
    # `ρ²`, `μ²`, `ν²` are the three real roots of the cubic
    # `κ³ + a κ² + b κ + c = 0` that the three confocal quadrics through `x`
    # share, one root in each of `[0, h₃²]`, `[h₃², h₂²]`, `[h₂², ∞)`.
    a = -(x₁² + x₂² + x₃² + h₃² + h₂²)
    b = (h₃² + h₂²) * x₁² + h₂² * x₂² + h₃² * x₃² + h₃² * h₂²
    # Trigonometric form of the roots of the depressed cubic `t³ + P t + Q = 0`,
    # `κ = t - a/3`; three real roots make `P ≤ 0` and the argument of `acos`
    # lie in `[-1, 1]`, up to rounding on the focal curves, hence the clamp.
    P = b - a^2 / 3
    Q = 2 * a^3 / 27 - a * b / 3 - h₃² * h₂² * x₁²
    m = 2 * sqrt(-P / 3)
    θ = acos(clamp(3 * Q / (P * m), -one(Q), one(Q))) / 3
    # Each root is kept in its interval, which rounding may leave by a few
    # units of the last place on a coordinate plane.
    ρ² = max(m * cos(θ) - a / 3, h₂²)
    μ² = clamp(m * cos(θ - 2 * oftype(θ, π) / 3) - a / 3, h₃², h₂²)
    # The smallest root from the product of the three, `ρ² μ² ν² = h₂² h₃² x₁²`,
    # rather than as a difference: it keeps its relative accuracy near `x₁ = 0`.
    ν = clamp(chart.h₂ * chart.h₃ * x[1] / sqrt(ρ² * μ²), -chart.h₃, chart.h₃)
    return SVector(promote(sqrt(ρ²), sqrt(μ²), ν)...)
end

"""
    change_coords(to::AbstractChart, from::AbstractChart, ξ) -> SVector

Coordinates in the chart `to` of the point whose coordinates in the chart
`from` are `ξ`. Between two spheroidal charts of the same axis and different
semi-focal distances, this is the change of coordinates that accompanies the
re-expansion of spheroidal harmonics.

# Examples
```jldoctest
julia> ξ₂ = change_coords(ProlateSpheroidalChart(2.0), ProlateSpheroidalChart(1.0), (0.3, 0.5, 3.0));

julia> coords_to_cartesian(ProlateSpheroidalChart(2.0), ξ₂) ≈ coords_to_cartesian(ProlateSpheroidalChart(1.0), (0.3, 0.5, 3.0))
true
```
"""
change_coords(to::AbstractChart{dim}, from::AbstractChart{dim}, ξ) where {dim} =
    cartesian_to_coords(to, coords_to_cartesian(from, ξ))

export AbstractChart, PolarChart, CylindricalChart, SphericalChart
export ProlateSpheroidalChart, OblateSpheroidalChart, EllipsoidalChart
export coords_to_cartesian, cartesian_to_coords, change_coords
