# Charts: a curvilinear system as a map between coordinates and points.
# References are independent of the code under test: the focal distances of the
# prolate system, the equation of the confocal quadrics, the Jacobian of the
# map by ForwardDiff, the implicit-function theorem for the parameters, and
# extended precision.

import Symbolics
using StaticArrays: SVector

# Points of the three charts, fixed once: every octant (and several azimuths)
# is represented.
const _ch_prolate_ξ = [
    (0.3, 0.4, 1.7), (-2.1, -0.6, 2.5), (1.2, 0.9, 1.05), (-0.7, -0.2, 4.0),
    (2.8, 0.55, 1.3), (-1.9, 0.1, 3.3), (0.9, -0.95, 1.8), (-3.0, -0.35, 1.2),
]
const _ch_oblate_ξ = [
    (0.3, 0.4, 1.7), (-2.1, -0.6, 0.5), (1.2, 0.9, 0.05), (-0.7, -0.2, 4.0),
    (2.8, 0.55, 1.3), (-1.9, 0.1, 0.3), (0.9, -0.95, 0.8), (-3.0, -0.35, 2.2),
]
# (ρ, μ, ν) with h₂ = 2, h₃ = 1: ρ > 2, 1 < μ < 2, |ν| < 1.
const _ch_ellipsoidal_ξ = [
    (2.5, 1.3, 0.4), (3.1, 1.8, -0.7), (2.05, 1.05, 0.95), (4.0, 1.5, -0.1),
    (2.2, 1.95, 0.6), (5.0, 1.1, -0.9), (2.8, 1.6, 0.25), (3.7, 1.4, -0.5),
]
# Cartesian points covering the eight octants.
const _ch_x = [
    SVector(s₁ * 0.7, s₂ * 1.9, s₃ * 1.3) for s₁ in (-1, 1) for s₂ in (-1, 1) for s₃ in (-1, 1)
]

_ch_cases = (
    (ProlateSpheroidalChart(1.5), _ch_prolate_ξ),
    (OblateSpheroidalChart(1.5), _ch_oblate_ξ),
    (EllipsoidalChart(2.0, 1.0), _ch_ellipsoidal_ξ),
)

# The orderings of the ellipsoidal coordinates.
_ch_ordered(h₂, h₃, ξ) = ξ[3]^2 <= h₃^2 <= ξ[2]^2 <= h₂^2 <= ξ[1]^2

# Allocation of a call, measured inside a function after warm-up.
function _ch_alloc(chart, ξ, x)
    coords_to_cartesian(chart, ξ)
    cartesian_to_coords(chart, x)
    a = @allocated coords_to_cartesian(chart, ξ)
    b = @allocated cartesian_to_coords(chart, x)
    return a, b
end

@testsection "Charts" begin

    @testsection "Round trips" begin
        for (chart, ξs) in _ch_cases, ξ in ξs
            x = coords_to_cartesian(chart, ξ)
            @test cartesian_to_coords(chart, x) ≈ SVector(ξ) rtol = 1.0e-12
        end
        for chart in (
                    PolarChart(), CylindricalChart(), SphericalChart(),
                    ProlateSpheroidalChart(1.5), OblateSpheroidalChart(1.5),
                ), x in _ch_x
            n = chart isa PolarChart ? 2 : 3
            xs = SVector{n}(x[1:n])
            @test coords_to_cartesian(chart, cartesian_to_coords(chart, xs)) ≈ xs rtol = 1.0e-12
        end
        @test cartesian_to_coords(PolarChart(), SVector(0.0, 2.0)) ≈ [2.0, π / 2]
        @test cartesian_to_coords(SphericalChart(), SVector(0.0, 0.0, -3.0)) ≈ [π, 0.0, 3.0]
        @test coords_to_cartesian(SphericalChart(), (π / 2, π / 2, 2.0)) ≈ [0.0, 2.0, 0.0] atol = 1.0e-15
        @test coords_to_cartesian(CylindricalChart(), (2.0, π, 5.0)) ≈ [-2.0, 0.0, 5.0] atol = 1.0e-15
        @test coords_to_cartesian(PolarChart(), (2.0, π / 2)) ≈ [0.0, 2.0] atol = 1.0e-15
        # the docstring examples
        @test coords_to_cartesian(ProlateSpheroidalChart(1.0), (0.0, 0.0, 2.0)) ≈ [√3, 0, 0]
    end

    @testsection "Ellipsoidal contract: mirror image in the quarter space" begin
        chart = EllipsoidalChart(2.0, 1.0)
        for x in _ch_x
            ξ = cartesian_to_coords(chart, x)
            @test coords_to_cartesian(chart, ξ) ≈ SVector(x[1], abs(x[2]), abs(x[3])) rtol = 1.0e-12
            @test sign(ξ[3]) == sign(x[1])
            @test _ch_ordered(2.0, 1.0, ξ)
        end
        for ξ in _ch_ellipsoidal_ξ
            @test _ch_ordered(2.0, 1.0, ξ)
        end
        # semi-axes constructor: h₂² = a₁² - a₃², h₃² = a₁² - a₂²
        @test EllipsoidalChart(3.0, 2.0, 1.0) == EllipsoidalChart(√8, √5)
        # EllipsoidalChart(h₂, h₃) promotes mixed arguments
        @test EllipsoidalChart(2, 1.0) isa EllipsoidalChart{Float64}
    end

    @testsection "Independent references" begin
        # prolate: 2cq = sum of the distances to the foci ±c 𝐞₃, 2cp = difference
        c = 1.5
        chart = ProlateSpheroidalChart(c)
        for ξ in _ch_prolate_ξ
            x = coords_to_cartesian(chart, ξ)
            d₊ = norm(x - SVector(0, 0, c))
            d₋ = norm(x + SVector(0, 0, c))
            @test d₊ + d₋ ≈ 2c * ξ[3] rtol = 1.0e-13
            @test d₋ - d₊ ≈ 2c * ξ[2] rtol = 1.0e-13
            # on the spheroid of semi-axes c q (axis) and c √(q²-1)
            @test (x[1]^2 + x[2]^2) / (c^2 * (ξ[3]^2 - 1)) + x[3]^2 / (c * ξ[3])^2 ≈ 1 rtol = 1.0e-13
            # and cartesian_to_coords finds these distances
            ξ′ = cartesian_to_coords(chart, x)
            @test d₊ + d₋ ≈ 2c * ξ′[3] rtol = 1.0e-13
        end
        # oblate: on the spheroid ϱ²/(c²(τ²+1)) + z²/(c²τ²) = 1
        chart = OblateSpheroidalChart(c)
        for ξ in _ch_oblate_ξ
            x = coords_to_cartesian(chart, ξ)
            ϱ² = x[1]^2 + x[2]^2
            @test ϱ² / (c^2 * (ξ[3]^2 + 1)) + x[3]^2 / (c^2 * ξ[3]^2) ≈ 1 rtol = 1.0e-13
            ξ′ = cartesian_to_coords(chart, x)
            @test ϱ² / (c^2 * (ξ′[3]^2 + 1)) + x[3]^2 / (c^2 * ξ′[3]^2) ≈ 1 rtol = 1.0e-13
            # the hyperboloid p = const: ϱ²/(c²(1-p²)) - z²/(c²p²) = 1
            @test ϱ² / (c^2 * (1 - ξ′[2]^2)) - x[3]^2 / (c^2 * ξ′[2]^2) ≈ 1 rtol = 1.0e-12
        end
        # ellipsoidal: the three quadrics of parameters ρ², μ², ν²
        h₂, h₃ = 2.0, 1.0
        chart = EllipsoidalChart(h₂, h₃)
        quadric(x, κ) = x[1]^2 / κ + x[2]^2 / (κ - h₃^2) + x[3]^2 / (κ - h₂^2)
        for x in _ch_x
            ξ = cartesian_to_coords(chart, x)
            for κ in ξ .^ 2
                @test quadric(x, κ) ≈ 1 rtol = 1.0e-9
            end
        end
        for ξ in _ch_ellipsoidal_ξ
            x = coords_to_cartesian(chart, ξ)
            for κ in ξ .^ 2
                @test quadric(x, κ) ≈ 1 rtol = 1.0e-12
            end
        end
        # a point of the ellipsoid of semi-axes a has ρ = a₁
        a₁, a₂, a₃ = 3.0, 2.0, 1.0
        chart = EllipsoidalChart(a₁, a₂, a₃)
        for (u, v) in ((0.3, 0.8), (1.1, 0.2), (2.0, 0.9), (0.6, 0.5))
            x = SVector(a₁ * sin(u) * cos(v), a₂ * sin(u) * sin(v), a₃ * cos(u))
            @test cartesian_to_coords(chart, x)[1] ≈ a₁ rtol = 1.0e-12
        end
    end

    @testsection "Degenerate points" begin
        # prolate: on the axis beyond the foci, on the focal segment, at a focus
        c = 1.5
        chart = ProlateSpheroidalChart(c)
        for z in (-4.0, -2.0, 2.0, 5.0)
            ξ = cartesian_to_coords(chart, SVector(0.0, 0.0, z))
            @test all(isfinite, ξ)
            @test ξ[3] ≈ abs(z) / c && ξ[2] ≈ sign(z)
            @test coords_to_cartesian(chart, ξ) ≈ [0, 0, z] atol = 1.0e-14
        end
        for z in (-1.4, -0.5, 0.0, 0.9, 1.4)
            ξ = cartesian_to_coords(chart, SVector(0.0, 0.0, z))
            @test ξ[3] ≈ 1 atol = 1.0e-15
            @test ξ[2] ≈ z / c
            @test coords_to_cartesian(chart, ξ) ≈ [0, 0, z] atol = 1.0e-14
        end
        # a point of the focal segment off the axis is not on it: the equator point
        ξ = cartesian_to_coords(chart, SVector(0.0, 0.0, c))
        @test all(isfinite, ξ) && ξ[3] ≈ 1 && ξ[2] ≈ 1

        # oblate: focal disk (z = 0, ϱ < c), focal circle, axis, equatorial plane beyond
        c = 1.5
        chart = OblateSpheroidalChart(c)
        for ϱ in (0.0, 0.4, 1.0, 1.49)
            x = SVector(ϱ * cos(0.6), ϱ * sin(0.6), 0.0)
            ξ = cartesian_to_coords(chart, x)
            @test all(isfinite, ξ)
            @test ξ[3] == 0
            @test ξ[2] ≈ sqrt(1 - ϱ^2 / c^2)
            @test coords_to_cartesian(chart, ξ) ≈ x atol = 1.0e-14
        end
        x = SVector(c * cos(0.6), c * sin(0.6), 0.0)           # focal circle
        ξ = cartesian_to_coords(chart, x)
        # τ = √(rounding of order 1e-16): absolute accuracy of order 1e-8 on the circle
        @test ξ[2] == 0 && abs(ξ[3]) < 1.0e-7 && ξ[1] ≈ 0.6
        @test coords_to_cartesian(chart, ξ) ≈ x atol = 1.0e-14
        for z in (-3.0, -0.2, 0.2, 2.0)                        # axis
            ξ = cartesian_to_coords(chart, SVector(0.0, 0.0, z))
            @test all(isfinite, ξ)
            @test ξ[2] ≈ sign(z) && ξ[3] ≈ abs(z) / c
            @test coords_to_cartesian(chart, ξ) ≈ [0, 0, z] atol = 1.0e-14
        end
        for ϱ in (1.6, 3.0)                                    # equatorial plane beyond the circle
            x = SVector(ϱ, 0.0, 0.0)
            ξ = cartesian_to_coords(chart, x)
            @test ξ[2] == 0 && ξ[3] ≈ sqrt(ϱ^2 / c^2 - 1)
            @test coords_to_cartesian(chart, ξ) ≈ x atol = 1.0e-14
        end

        # ellipsoidal: the three coordinate planes, and the axes
        chart = EllipsoidalChart(2.0, 1.0)
        planes = (
            SVector(0.0, 1.3, 0.8), SVector(0.0, 0.0, 0.8), SVector(0.0, 1.3, 0.0),     # x₁ = 0
            SVector(0.9, 0.0, 0.7), SVector(2.6, 0.0, 1.1), SVector(0.4, 0.0, 1.9),     # x₂ = 0
            SVector(0.9, 1.4, 0.0), SVector(0.2, 2.3, 0.0), SVector(3.1, 0.6, 0.0),     # x₃ = 0
            SVector(1.2, 0.0, 0.0), SVector(0.0, 2.0, 0.0), SVector(0.0, 0.0, 1.5),     # axes
        )
        # On these planes the roots of the cubic meet the bounds `h₃²`, `h₂²`:
        # the coordinates stay in their intervals and map back to the point.
        in_range(ξ) = ξ[1] ≥ 2.0 && 1.0 ≤ ξ[2] ≤ 2.0 && -1.0 ≤ ξ[3] ≤ 1.0
        for x in planes
            ξ = cartesian_to_coords(chart, x)
            @test all(isfinite, ξ)
            @test in_range(ξ)
            @test coords_to_cartesian(chart, ξ) ≈ x rtol = 1.0e-6 atol = 1.0e-6
        end
        # A point at a distance `ε` of a plane is recovered to the accuracy of
        # the square root of a rounded difference, about `1e-8` in absolute terms.
        for ε in (1.0e-3, 1.0e-6, 1.0e-9)
            for x in (SVector(ε, 1.3, 0.8), SVector(0.9, ε, 0.7), SVector(0.9, 1.4, ε))
                ξ = cartesian_to_coords(chart, x)
                @test in_range(ξ)
                @test coords_to_cartesian(chart, ξ) ≈ x rtol = 1.0e-6 atol = 1.0e-7
            end
        end
        # On the axis of a spheroidal chart, beyond the foci, `|p| = 1` up to
        # rounding and never more, so that it stays an argument of `legendre_P`.
        for chart in (ProlateSpheroidalChart(1.3), OblateSpheroidalChart(1.3)), z in (-7.3, -1.9, 2.2, 1.0e6 / 3, 3.0e-7 + 1.3)
            p = cartesian_to_coords(chart, SVector(0.0, 0.0, z))[2]
            @test abs(p) ≤ 1.0
            @test abs(p) ≈ 1.0
            @test all(isfinite, legendre_P(4, 2, p))
        end
    end

    @testsection "Lamé coefficients and unit vectors" begin
        for (chart, ξs) in _ch_cases, ξ in ξs
            ξv = SVector(ξ)
            J = ForwardDiff.jacobian(u -> coords_to_cartesian(chart, u), ξv)
            χ = TensND._lame(chart, ξ)
            R = TensND._frame(chart, ξ)
            @test collect(χ) ≈ [norm(J[:, i]) for i in 1:3] rtol = 1.0e-12
            @test Matrix(R) ≈ J ./ SVector(χ...)' rtol = 1.0e-12
            @test R' * R ≈ I atol = 1.0e-13
            @test det(R) ≈ (chart isa EllipsoidalChart ? -1 : 1) atol = 1.0e-13
        end
    end

    @testsection "ForwardDiff: Jacobian of the inverse" begin
        # J = R diag(χ) ⇒ J⁻¹ = diag(1/χ) Rᵀ
        for (chart, ξs) in _ch_cases, ξ in ξs
            x = coords_to_cartesian(chart, ξ)
            Jinv = ForwardDiff.jacobian(y -> cartesian_to_coords(chart, y), x)
            ξ′ = cartesian_to_coords(chart, x)
            χ = TensND._lame(chart, ξ′)
            R = TensND._frame(chart, ξ′)
            @test Jinv ≈ (R ./ SVector(χ...)')' rtol = 1.0e-8
        end
        # a Jacobian at points of every sign for the spheroidal charts
        for chart in (ProlateSpheroidalChart(1.5), OblateSpheroidalChart(1.5)), x in _ch_x
            Jinv = ForwardDiff.jacobian(y -> cartesian_to_coords(chart, y), x)
            ξ′ = cartesian_to_coords(chart, x)
            χ = TensND._lame(chart, ξ′)
            R = TensND._frame(chart, ξ′)
            @test Jinv ≈ (R ./ SVector(χ...)')' rtol = 1.0e-8
        end
    end

    @testsection "ForwardDiff: parameters (implicit-function theorem)" begin
        # x = X(ξ, c) with ξ = ξ(x, c): ∂ξ/∂c = -J⁻¹ ∂X/∂c
        x = SVector(0.7, 1.9, 1.3)
        for (mk, c₀) in (
                (c -> ProlateSpheroidalChart(c), 1.5),
                (c -> OblateSpheroidalChart(c), 1.5),
            )
            ξ = cartesian_to_coords(mk(c₀), x)
            dξ = ForwardDiff.derivative(c -> cartesian_to_coords(mk(c), x), c₀)
            dX = ForwardDiff.derivative(c -> coords_to_cartesian(mk(c), ξ), c₀)
            J = ForwardDiff.jacobian(u -> coords_to_cartesian(mk(c₀), u), ξ)
            @test dξ ≈ -(J \ dX) rtol = 1.0e-8
        end
        # ellipsoidal: h₂ and h₃ separately
        h₂, h₃ = 2.0, 1.0
        ξ = cartesian_to_coords(EllipsoidalChart(h₂, h₃), x)
        J = ForwardDiff.jacobian(u -> coords_to_cartesian(EllipsoidalChart(h₂, h₃), u), ξ)
        dξ₂ = ForwardDiff.derivative(h -> cartesian_to_coords(EllipsoidalChart(h, h₃), x), h₂)
        dX₂ = ForwardDiff.derivative(h -> coords_to_cartesian(EllipsoidalChart(h, h₃), ξ), h₂)
        @test dξ₂ ≈ -(J \ dX₂) rtol = 1.0e-8
        dξ₃ = ForwardDiff.derivative(h -> cartesian_to_coords(EllipsoidalChart(h₂, h), x), h₃)
        dX₃ = ForwardDiff.derivative(h -> coords_to_cartesian(EllipsoidalChart(h₂, h), ξ), h₃)
        @test dξ₃ ≈ -(J \ dX₃) rtol = 1.0e-8
        # the three semi-axes constructor is differentiable too
        dξ = ForwardDiff.derivative(a -> cartesian_to_coords(EllipsoidalChart(3.0, 2.0, a), x), 1.0)
        @test all(isfinite, dξ)
        # the coordinates and the parameter dual together
        c₀ = 1.5
        f(c, y) = cartesian_to_coords(ProlateSpheroidalChart(c), y)
        g = ForwardDiff.gradient(v -> f(v[1], SVector(v[2], v[3], v[4]))[3], [c₀, x...])
        Jinv = ForwardDiff.jacobian(y -> f(c₀, y), x)
        @test g[2:4] ≈ Jinv[3, :] rtol = 1.0e-10
        @test g[1] ≈ ForwardDiff.derivative(c -> f(c, x)[3], c₀) rtol = 1.0e-10
    end

    @testsection "ForwardDiff: nested derivatives" begin
        # X(ξ(x, c), c) = x for every c: all its derivatives vanish, the second
        # one through the nested duals of cartesian_to_coords.
        x = SVector(0.7, 1.9, 1.3)
        for mk in (
                c -> ProlateSpheroidalChart(c), c -> OblateSpheroidalChart(c),
                c -> EllipsoidalChart(2.0 + c, 1.0),
            )
            id(c) = coords_to_cartesian(mk(c), cartesian_to_coords(mk(c), x))
            @test ForwardDiff.derivative(id, 1.0) ≈ zeros(3) atol = 1.0e-10
            @test ForwardDiff.derivative(c -> ForwardDiff.derivative(id, c), 1.0) ≈ zeros(3) atol = 1.0e-9
            # mixed second derivative with respect to a coordinate and the parameter
            idx(c, y) = coords_to_cartesian(mk(c), cartesian_to_coords(mk(c), y))
            H = ForwardDiff.derivative(c -> ForwardDiff.jacobian(y -> idx(c, y), x), 1.0)
            @test H ≈ zeros(3, 3) atol = 1.0e-9
        end
        # the Hessian of the coordinate q of the prolate chart against the
        # closed form of its Laplacian-free part: for q = (d₊ + d₋)/(2c),
        # ∇q·∇q = (q² - 1)/(c²(q² - p²)) (= 1/χ_q²).
        c = 1.5
        chart = ProlateSpheroidalChart(c)
        q(y) = cartesian_to_coords(chart, y)[3]
        ξ = cartesian_to_coords(chart, x)
        gq = ForwardDiff.gradient(q, x)
        @test dot(gq, gq) ≈ (ξ[3]^2 - 1) / (c^2 * (ξ[3]^2 - ξ[2]^2)) rtol = 1.0e-10
        # the second derivative of q along 𝐞₃ against the second derivative of
        # the distances d± to the foci: ∂²d/∂z² = (ϱ²)/d³
        ϱ² = x[1]^2 + x[2]^2
        d₊, d₋ = norm(x - SVector(0, 0, c)), norm(x + SVector(0, 0, c))
        h₃₃ = ForwardDiff.hessian(q, x)[3, 3]
        @test h₃₃ ≈ (ϱ² / d₊^3 + ϱ² / d₋^3) / (2c) rtol = 1.0e-9
    end

    @testsection "BigFloat round trip" begin
        setprecision(BigFloat, 256) do
            for (chart, ξs) in (
                        (ProlateSpheroidalChart(big"1.5"), _ch_prolate_ξ),
                        (OblateSpheroidalChart(big"1.5"), _ch_oblate_ξ),
                        (EllipsoidalChart(big"2.0", big"1.0"), _ch_ellipsoidal_ξ),
                    ), ξ in ξs
                ξb = SVector(BigFloat.(ξ)...)
                x = coords_to_cartesian(chart, ξb)
                @test eltype(x) == BigFloat
                @test maximum(abs.(cartesian_to_coords(chart, x) - ξb)) < big"1.0e-60"
            end
        end
    end

    @testsection "Type stability and allocations" begin
        Dual = ForwardDiff.Dual{Nothing, Float64, 3}
        for (chart, ξs) in _ch_cases
            ξ = SVector(ξs[1])
            x = coords_to_cartesian(chart, ξ)
            @inferred coords_to_cartesian(chart, ξ)
            @inferred coords_to_cartesian(chart, ξs[1])
            @inferred cartesian_to_coords(chart, x)
            @inferred change_coords(chart, chart, ξ)
            ξd = SVector(ForwardDiff.Dual(ξ[1], 1.0, 0.0, 0.0), ForwardDiff.Dual(ξ[2], 0.0, 1.0, 0.0), ForwardDiff.Dual(ξ[3], 0.0, 0.0, 1.0))
            xd = coords_to_cartesian(chart, ξd)
            @test eltype(@inferred coords_to_cartesian(chart, ξd)) <: ForwardDiff.Dual
            @test eltype(@inferred cartesian_to_coords(chart, xd)) <: ForwardDiff.Dual
            @test eltype(@inferred change_coords(chart, chart, ξd)) <: ForwardDiff.Dual
            # a dual parameter
            cd = chart isa EllipsoidalChart ? EllipsoidalChart(ForwardDiff.Dual(2.0, 1.0), ForwardDiff.Dual(1.0, 0.0)) :
                chart isa ProlateSpheroidalChart ? ProlateSpheroidalChart(ForwardDiff.Dual(1.5, 1.0)) :
                OblateSpheroidalChart(ForwardDiff.Dual(1.5, 1.0))
            @inferred coords_to_cartesian(cd, ξ)
            @inferred cartesian_to_coords(cd, x)
            a, b = _ch_alloc(chart, ξ, x)
            @test a == 0
            @test b == 0
        end
        for chart in (PolarChart(), CylindricalChart(), SphericalChart())
            n = chart isa PolarChart ? 2 : 3
            ξ = SVector{n}(ntuple(i -> 0.3i, n))
            @inferred coords_to_cartesian(chart, ξ)
            @inferred cartesian_to_coords(chart, SVector{n}(ntuple(i -> 0.4i, n)))
        end
    end

    @testsection "change_coords" begin
        ξ = (0.4, 0.3, 2.5)
        P₁, P₂ = ProlateSpheroidalChart(1.0), ProlateSpheroidalChart(2.5)
        x = coords_to_cartesian(P₁, ξ)
        ξ₂ = change_coords(P₂, P₁, ξ)
        @test coords_to_cartesian(P₂, ξ₂) ≈ x rtol = 1.0e-12
        @test change_coords(P₁, P₂, ξ₂) ≈ [ξ...] rtol = 1.0e-12
        @test change_coords(P₁, P₁, ξ) ≈ [ξ...] rtol = 1.0e-12
        O₁, O₂ = OblateSpheroidalChart(0.8), OblateSpheroidalChart(1.7)
        ξo = (-1.3, -0.5, 0.9)
        @test coords_to_cartesian(O₂, change_coords(O₂, O₁, ξo)) ≈ coords_to_cartesian(O₁, ξo) rtol = 1.0e-12
        # prolate ↔ oblate
        @test coords_to_cartesian(O₁, change_coords(O₁, P₁, ξ)) ≈ x rtol = 1.0e-12
        @test coords_to_cartesian(P₂, change_coords(P₂, O₁, ξo)) ≈ coords_to_cartesian(O₁, ξo) rtol = 1.0e-12
        # spheroidal ↔ spherical: r = |x|, cos θ = x₃/r, ϕ unchanged
        θ, ϕ, r = change_coords(SphericalChart(), P₁, ξ)
        @test r ≈ norm(x) && cos(θ) ≈ x[3] / norm(x) && ϕ ≈ ξ[1]
        @test coords_to_cartesian(P₂, change_coords(P₂, SphericalChart(), (0.7, -1.1, 3.0))) ≈
            coords_to_cartesian(SphericalChart(), (0.7, -1.1, 3.0)) rtol = 1.0e-12
        # cylindrical ↔ spheroidal
        @test coords_to_cartesian(CylindricalChart(), change_coords(CylindricalChart(), O₁, ξo)) ≈
            coords_to_cartesian(O₁, ξo) rtol = 1.0e-12
        # ellipsoidal ↔ prolate (a point of the quarter space)
        E = EllipsoidalChart(2.0, 1.0)
        ξe = (2.5, 1.3, 0.4)
        @test coords_to_cartesian(P₁, change_coords(P₁, E, ξe)) ≈ coords_to_cartesian(E, ξe) rtol = 1.0e-12
        # dimensions must agree
        @test_throws MethodError change_coords(PolarChart(), P₁, ξ)
    end

    @testsection "Limits" begin
        # a large q: the prolate chart approaches the spherical one
        c = 1.5
        chart = ProlateSpheroidalChart(c)
        for p in (-0.8, 0.1, 0.9)
            q = 1.0e6
            x = coords_to_cartesian(chart, (0.7, p, q))
            θ, ϕ, r = cartesian_to_coords(SphericalChart(), x)
            @test r ≈ c * q rtol = 1.0e-9
            @test cos(θ) ≈ p rtol = 1.0e-9
            @test ϕ ≈ 0.7
        end
        # the oblate chart too, with p = cos θ
        chart = OblateSpheroidalChart(c)
        for p in (-0.8, 0.1, 0.9)
            τ = 1.0e6
            θ, ϕ, r = cartesian_to_coords(SphericalChart(), coords_to_cartesian(chart, (0.7, p, τ)))
            @test r ≈ c * τ rtol = 1.0e-9
            @test cos(θ) ≈ p rtol = 1.0e-9
        end
        # q → 1: the segment between the foci; p = x₃/c
        x = coords_to_cartesian(ProlateSpheroidalChart(c), (0.2, 0.6, 1.0))
        @test x ≈ [0, 0, 0.6c] atol = 1.0e-15
        # τ → 0: the focal disk, ϱ = c √(1 - p²)
        x = coords_to_cartesian(OblateSpheroidalChart(c), (0.0, 0.6, 0.0))
        @test x ≈ [0.8c, 0, 0] atol = 1.0e-15
        # ellipsoidal, ρ large: x ≈ ρ times a unit vector (confocal ellipsoids become spheres)
        chart = EllipsoidalChart(2.0, 1.0)
        ρ = 1.0e6
        x = coords_to_cartesian(chart, (ρ, 1.3, 0.4))
        @test norm(x) ≈ ρ rtol = 1.0e-6
    end

    @testsection "Parameters: errors and symbols" begin
        @test_throws DomainError ProlateSpheroidalChart(-1.0)
        @test_throws DomainError ProlateSpheroidalChart(0.0)
        @test_throws DomainError OblateSpheroidalChart(0.0)
        @test_throws DomainError OblateSpheroidalChart(-2)
        @test_throws DomainError EllipsoidalChart(0.5, 0.7)
        @test_throws DomainError EllipsoidalChart(1.0, 0.0)
        @test_throws DomainError EllipsoidalChart(1.0, 1.0)
        @test_throws DomainError EllipsoidalChart(1.0, 2.0, 3.0)
        @test ProlateSpheroidalChart(2) isa ProlateSpheroidalChart{Int}
        @test ProlateSpheroidalChart(1.5f0) isa ProlateSpheroidalChart{Float32}
        @test eltype(coords_to_cartesian(ProlateSpheroidalChart(1.5f0), (0.3f0, 0.2f0, 2.0f0))) == Float32

        # SymPy
        cs, h₂s, h₃s = symbols("c h2 h3", positive = true)
        ϕs, ps, qs = symbols("ϕ p q", real = true)
        for (chart, ξ, ξ₀) in (
                (ProlateSpheroidalChart(cs), (ϕs, ps, qs), (0.3, 0.4, 1.7)),
                (OblateSpheroidalChart(cs), (ϕs, ps, qs), (0.3, 0.4, 1.7)),
            )
            @test chart.c == cs
            X = coords_to_cartesian(chart, ξ)
            Xnum = coords_to_cartesian(chart isa ProlateSpheroidalChart ? ProlateSpheroidalChart(1.5) : OblateSpheroidalChart(1.5), ξ₀)
            at = Dict(cs => 1.5, ϕs => ξ₀[1], ps => ξ₀[2], qs => ξ₀[3])
            @test [Float64(X[i].subs(at)) for i in 1:3] ≈ Xnum rtol = 1.0e-12
        end
        E = EllipsoidalChart(h₂s, h₃s)
        ρs, μs, νs = symbols("ρ μ", positive = true)..., symbols("ν", real = true)
        X = coords_to_cartesian(E, (ρs, μs, νs))
        at = Dict(h₂s => 2.0, h₃s => 1.0, ρs => 2.5, μs => 1.3, νs => 0.4)
        @test [Float64(X[i].subs(at)) for i in 1:3] ≈ coords_to_cartesian(EllipsoidalChart(2.0, 1.0), (2.5, 1.3, 0.4)) rtol = 1.0e-12

        # Symbolics
        Symbolics.@variables c ϕ p q ρ μ ν h₂ h₃
        val(e, d) = Float64(Symbolics.value(Symbolics.substitute(e, d; fold = Val(true))))
        P = ProlateSpheroidalChart(c)
        X = coords_to_cartesian(P, (ϕ, p, q))
        d = Dict(c => 1.5, ϕ => 0.3, p => 0.4, q => 1.7)
        @test [val(X[i], d) for i in 1:3] ≈ coords_to_cartesian(ProlateSpheroidalChart(1.5), (0.3, 0.4, 1.7)) rtol = 1.0e-12
        O = OblateSpheroidalChart(c)
        X = coords_to_cartesian(O, (ϕ, p, q))
        @test [val(X[i], d) for i in 1:3] ≈ coords_to_cartesian(OblateSpheroidalChart(1.5), (0.3, 0.4, 1.7)) rtol = 1.0e-12
        E = EllipsoidalChart(h₂, h₃)
        X = coords_to_cartesian(E, (ρ, μ, ν))
        d = Dict(h₂ => 2.0, h₃ => 1.0, ρ => 2.5, μ => 1.3, ν => 0.4)
        @test [val(X[i], d) for i in 1:3] ≈ coords_to_cartesian(EllipsoidalChart(2.0, 1.0), (2.5, 1.3, 0.4)) rtol = 1.0e-12
        # the Lamé coefficients and the frame accept symbols too
        χ = TensND._lame(E, (ρ, μ, ν))
        @test [val(χ[i], d) for i in 1:3] ≈ collect(TensND._lame(EllipsoidalChart(2.0, 1.0), (2.5, 1.3, 0.4))) rtol = 1.0e-12
    end
end
