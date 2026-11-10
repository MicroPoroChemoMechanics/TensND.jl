# The confocal systems built from a chart: prolate and oblate spheroidal,
# ellipsoidal. Exact (SymPy, Symbolics) and pointwise versions, against closed
# forms, against each other, and against the exact Cartesian operators.

import Symbolics
using StaticArrays: SVector

# Numeric value of a SymPy expression at a point.
_cf_val(e, d) = Float64(e.subs(d))

# The Cartesian components of the unit vectors are those of the chart.
_cf_local(chart, ξ, V) = TensND._frame(chart, ξ)' * V(coords_to_cartesian(chart, ξ))

@testsection "Confocal coordinate systems" begin

    c = symbols("c", positive = true)
    τ = symbols("τ", positive = true)
    (ϕ, p, q) = get_coords(coorsys_spheroidal())
    s∂² = tsimplify ∘ laplacian

    @testsection "Prolate spheroidal (unchanged behavior)" begin
        S = coorsys_spheroidal()
        @test get_coords(S) == (ϕ, p, q)
        OM = position_vector(S)
        χ = lame_coefficients(S)
        expected = (c^2 * (1 - p^2) * (q^2 - 1), c^2 * (q^2 - p^2) / (1 - p^2), c^2 * (q^2 - p^2) / (q^2 - 1))
        for i in 1:3
            @test iszero(tsimplify(χ[i]^2 - expected[i]))
        end
        @test tsimplify(laplacian(OM[1]^2, S)) == 2
        @test iszero(tsimplify(laplacian(OM[1]^2 - OM[3]^2, S)))
        @test iszero(tsimplify(laplacian(OM[1] * OM[3], S)))
        @test iszero(tsimplify(laplacian(OM[2] * OM[3], S)))
        # the Laplacian of a function of p alone
        f = SymFunction("f", real = true)(p)
        @test iszero(
            tsimplify(laplacian(f, S) - ∂((1 - p^2) * ∂(f, p), p) / (c^2 * (q^2 - p^2)))
        )
        # the position vector has the components of the chart
        ξ₀ = (0.3, 0.4, 1.7)
        p̄, q̄ = symbols("p̄ q̄", positive = true)
        at = Dict(c => 1.5, ϕ => ξ₀[1], p => ξ₀[2], q => ξ₀[3], p̄ => sqrt(1 - ξ₀[2]^2), q̄ => sqrt(ξ₀[3]^2 - 1))
        @test [_cf_val(OM[i], at) for i in 1:3] ≈ coords_to_cartesian(ProlateSpheroidalChart(1.5), ξ₀) rtol = 1.0e-12
    end

    @testsection "Oblate spheroidal" begin
        O = coorsys_spheroidal_oblate()
        @test get_coords(O) == (ϕ, p, τ)
        OM = position_vector(O)
        χ = lame_coefficients(O)
        closed = (
            c * sqrt(1 - p^2) * sqrt(τ^2 + 1), c * sqrt(τ^2 + p^2) / sqrt(1 - p^2),
            c * sqrt(τ^2 + p^2) / sqrt(τ^2 + 1),
        )
        for i in 1:3
            @test iszero(tsimplify(χ[i]^2 - closed[i]^2))
        end
        @test tsimplify(laplacian(OM[1]^2, O)) == 2
        @test iszero(tsimplify(laplacian(OM[1]^2 - OM[3]^2, O)))
        @test iszero(tsimplify(laplacian(OM[1] * OM[3], O)))
        f = SymFunction("f", real = true)(p)
        @test iszero(
            tsimplify(laplacian(f, O) - ∂((1 - p^2) * ∂(f, p), p) / (c^2 * (τ^2 + p^2)))
        )
        p̄, τ̄ = symbols("p̄ τ̄", positive = true)
        ξ₀ = (0.3, 0.4, 0.7)
        at = Dict(c => 1.5, ϕ => ξ₀[1], p => ξ₀[2], τ => ξ₀[3], p̄ => sqrt(1 - ξ₀[2]^2), τ̄ => sqrt(ξ₀[3]^2 + 1))
        @test [_cf_val(OM[i], at) for i in 1:3] ≈ coords_to_cartesian(OblateSpheroidalChart(1.5), ξ₀) rtol = 1.0e-12

        # the oblate system is the prolate one continued to c → -i c, q → i τ
        S = coorsys_spheroidal()
        IM = Sym(im)
        continued(e) = tsubs(e, c => -IM * c, q => IM * τ)
        χˢ = lame_coefficients(S)
        # dq = i dτ: the third Lamé coefficient changes sign under the continuation
        for (i, σ) in zip(1:3, (1, 1, -1))
            @test iszero(tsimplify(σ * continued(χˢ[i]^2) - χ[i]^2))
        end
        @test iszero(tsimplify(continued(tsimplify(laplacian(f, S))) - tsimplify(laplacian(f, O))))
    end

    @testsection "Ellipsoidal" begin
        E = coorsys_ellipsoidal()
        ρ, μ = symbols("ρ μ", positive = true)
        ν = symbols("ν", real = true)
        h₂, h₃ = symbols("h₂ h₃", positive = true)
        ρ̄₂, ρ̄₃, μ̄₂, μ̄₃, ν̄₂, ν̄₃ = symbols("ρ̄₂ ρ̄₃ μ̄₂ μ̄₃ ν̄₂ ν̄₃", positive = true)
        @test get_coords(E) == (ρ, μ, ν)
        OM = position_vector(E)
        # position vector: components of the chart, at a point
        ξ₀ = (2.5, 1.3, -0.4)
        d = Dict(h₂ => 2.0, h₃ => 1.0, ρ => ξ₀[1], μ => ξ₀[2], ν => ξ₀[3])
        full = merge(
            d, Dict(
                ρ̄₂ => sqrt(ξ₀[1]^2 - 4), ρ̄₃ => sqrt(ξ₀[1]^2 - 1), μ̄₂ => sqrt(4 - ξ₀[2]^2),
                μ̄₃ => sqrt(ξ₀[2]^2 - 1), ν̄₂ => sqrt(4 - ξ₀[3]^2), ν̄₃ => sqrt(1 - ξ₀[3]^2),
            )
        )
        chart = EllipsoidalChart(2.0, 1.0)
        @test [_cf_val(OM[i], full) for i in 1:3] ≈ coords_to_cartesian(chart, ξ₀) rtol = 1.0e-12
        # the three confocal quadrics through the point
        x = [_cf_val(OM[i], full) for i in 1:3]
        for κ in ξ₀ .^ 2
            @test x[1]^2 / κ + x[2]^2 / (κ - 1) + x[3]^2 / (κ - 4) ≈ 1 rtol = 1.0e-12
        end

        # Lamé coefficients
        χ = lame_coefficients(E)
        closed = (
            (ρ^2 - μ^2) * (ρ^2 - ν^2) / ((ρ^2 - h₃^2) * (ρ^2 - h₂^2)),
            (ρ^2 - μ^2) * (μ^2 - ν^2) / ((μ^2 - h₃^2) * (h₂^2 - μ^2)),
            (ρ^2 - ν^2) * (μ^2 - ν^2) / ((h₃^2 - ν^2) * (h₂^2 - ν^2)),
        )
        for i in 1:3
            @test iszero(tsimplify(χ[i]^2 - closed[i]))
        end

        # Laplacians of polynomials of the Cartesian coordinates
        @test tsimplify(laplacian(OM[1]^2, E)) == 2
        @test tsimplify(laplacian(OM[2]^2, E)) == 2
        @test iszero(tsimplify(laplacian(OM[1]^2 - OM[3]^2, E)))
        # `tsimplify` of the next two takes about 20 s each: checked at points
        # (the auxiliary symbols replaced by their values) instead
        for ξ in ((2.5, 1.3, -0.4), (3.3, 1.7, 0.8))
            pt = Dict(
                h₂ => 2.0, h₃ => 1.0, ρ => ξ[1], μ => ξ[2], ν => ξ[3],
                ρ̄₂ => sqrt(ξ[1]^2 - 4), ρ̄₃ => sqrt(ξ[1]^2 - 1), μ̄₂ => sqrt(4 - ξ[2]^2),
                μ̄₃ => sqrt(ξ[2]^2 - 1), ν̄₂ => sqrt(4 - ξ[3]^2), ν̄₃ => sqrt(1 - ξ[3]^2),
            )
            @test abs(_cf_val(laplacian(OM[2] * OM[3], E), pt)) < 1.0e-10
            @test abs(_cf_val(laplacian(OM[1] * OM[2] * OM[3], E), pt)) < 1.0e-10
            @test _cf_val(laplacian(OM[1]^2 + 3 * OM[3]^2, E), pt) ≈ 8 rtol = 1.0e-10
        end

        # a function of ρ alone: the Laplacian in closed form
        f = SymFunction("f", real = true)(ρ)
        r₂, r₃ = sqrt(ρ^2 - h₂^2), sqrt(ρ^2 - h₃^2)
        Δρ = r₃ * r₂ * ∂(r₃ * r₂ * ∂(f, ρ), ρ) / ((ρ^2 - μ^2) * (ρ^2 - ν^2))
        @test iszero(tsimplify(laplacian(f, E) - Δρ))
    end

    @testsection "Christoffel symbols: exact against pointwise" begin
        p̄, q̄ = symbols("p̄ q̄", positive = true)
        p̄o, τ̄ = symbols("p̄ τ̄", positive = true)
        ρ, μ = symbols("ρ μ", positive = true)
        ν = symbols("ν", real = true)
        h₂, h₃ = symbols("h₂ h₃", positive = true)
        ρ̄₂, ρ̄₃, μ̄₂, μ̄₃, ν̄₂, ν̄₃ = symbols("ρ̄₂ ρ̄₃ μ̄₂ μ̄₃ ν̄₂ ν̄₃", positive = true)

        ξ₀ = [0.3, 0.4, 1.7]
        Γsym = map(γ -> _cf_val(γ, Dict(c => 1.5, ϕ => ξ₀[1], p => ξ₀[2], q => ξ₀[3], p̄ => sqrt(1 - ξ₀[2]^2), q̄ => sqrt(ξ₀[3]^2 - 1))), christoffel(coorsys_spheroidal()))
        @test christoffel(coorsys_spheroidal(PointwiseGeometry, 1.5), ξ₀) ≈ Γsym rtol = 1.0e-10 atol = 1.0e-12

        ξ₀ = [0.3, 0.4, 0.7]
        at = Dict(c => 1.5, ϕ => ξ₀[1], p => ξ₀[2], τ => ξ₀[3], p̄o => sqrt(1 - ξ₀[2]^2), τ̄ => sqrt(ξ₀[3]^2 + 1))
        Γsym = map(γ -> _cf_val(γ, at), christoffel(coorsys_spheroidal_oblate()))
        @test christoffel(coorsys_spheroidal_oblate(PointwiseGeometry, 1.5), ξ₀) ≈ Γsym rtol = 1.0e-10 atol = 1.0e-12

        ξ₀ = [2.5, 1.3, -0.4]
        at = Dict(
            h₂ => 2.0, h₃ => 1.0, ρ => ξ₀[1], μ => ξ₀[2], ν => ξ₀[3],
            ρ̄₂ => sqrt(ξ₀[1]^2 - 4), ρ̄₃ => sqrt(ξ₀[1]^2 - 1), μ̄₂ => sqrt(4 - ξ₀[2]^2),
            μ̄₃ => sqrt(ξ₀[2]^2 - 1), ν̄₂ => sqrt(4 - ξ₀[3]^2), ν̄₃ => sqrt(1 - ξ₀[3]^2),
        )
        Γsym = map(γ -> _cf_val(γ, at), christoffel(coorsys_ellipsoidal()))
        @test christoffel(coorsys_ellipsoidal(PointwiseGeometry, 2.0, 1.0), ξ₀) ≈ Γsym rtol = 1.0e-10 atol = 1.0e-12
        # the Lamé coefficients too
        χsym = [_cf_val(χ, at) for χ in lame_coefficients(coorsys_ellipsoidal())]
        @test collect(lame_coefficients(coorsys_ellipsoidal(PointwiseGeometry, 2.0, 1.0), ξ₀)) ≈ χsym rtol = 1.0e-12
    end

    @testsection "Pointwise systems: operators against the Cartesian ones" begin
        # a scalar field and a vector field, polynomials of the Cartesian components
        f(x) = x[1]^3 + x[1] * x[2] * x[3] + 2 * x[3]^2 + x[2]^2
        ∇f(x) = [3x[1]^2 + x[2] * x[3], x[1] * x[3] + 2x[2], x[1] * x[2] + 4x[3]]
        Δf(x) = 6x[1] + 2 + 4
        V(x) = [x[1]^2, x[2] * x[3], x[1] * x[3]^2]
        divV(x) = 2x[1] + x[3] + 2x[1] * x[3]
        # harmonic polynomial
        h(x) = x[1]^2 - x[3]^2 + x[2] * x[3]

        cases = (
            (ProlateSpheroidalChart(1.5), coorsys_spheroidal(PointwiseGeometry, 1.5), [0.3, 0.4, 1.7]),
            (ProlateSpheroidalChart(1.5), coorsys_spheroidal(PointwiseGeometry, 1.5), [-2.1, -0.6, 2.5]),
            (OblateSpheroidalChart(1.5), coorsys_spheroidal_oblate(PointwiseGeometry, 1.5), [0.3, 0.4, 0.7]),
            (OblateSpheroidalChart(1.5), coorsys_spheroidal_oblate(PointwiseGeometry, 1.5), [-1.9, -0.3, 1.2]),
            (EllipsoidalChart(2.0, 1.0), coorsys_ellipsoidal(PointwiseGeometry, 2.0, 1.0), [2.5, 1.3, -0.4]),
            (EllipsoidalChart(2.0, 1.0), coorsys_ellipsoidal(PointwiseGeometry, 2.0, 1.0), [3.1, 1.8, 0.7]),
        )
        for (chart, CS, ξ) in cases
            x = coords_to_cartesian(chart, ξ)
            F = ζ -> f(coords_to_cartesian(chart, ζ))
            @test laplacian(F, CS)(ξ) ≈ Δf(x) rtol = 1.0e-9
            @test laplacian(ζ -> h(coords_to_cartesian(chart, ζ)), CS)(ξ) ≈ 0 atol = 1.0e-9
            g = grad(F, CS)(ξ)
            @test g isa AbstractTens
            @test collect(canonical_components(g)) ≈ ∇f(x) rtol = 1.0e-9
            # a vector field is given by its components on the unit vectors
            W = ζ -> _cf_local(chart, ζ, V)
            @test divergence(W, CS)(ξ) ≈ divV(x) rtol = 1.0e-9
            # the generic system built from the position vector agrees
            CSg = CoordinateSystem(ζ -> Vector(coords_to_cartesian(chart, ζ)), 3)
            @test laplacian(F, CSg)(ξ) ≈ laplacian(F, CS)(ξ) rtol = 1.0e-9
            @test collect(canonical_components(grad(F, CSg)(ξ))) ≈ collect(canonical_components(g)) rtol = 1.0e-9
            @test divergence(W, CSg)(ξ) ≈ divergence(W, CS)(ξ) rtol = 1.0e-9
            @test collect(lame_coefficients(CSg, ξ)) ≈ collect(lame_coefficients(CS, ξ)) rtol = 1.0e-10
            @test christoffel(CSg, ξ) ≈ christoffel(CS, ξ) rtol = 1.0e-7 atol = 1.0e-9
            @test Array(normalized_basis(CS, ξ)) ≈ Matrix(TensND._frame(chart, ξ)) rtol = 1.0e-12
        end

        # derivative of a pointwise Laplacian with respect to the semi-focal distance:
        # at fixed coordinates the point is x(ξ, c) = c x(ξ, 1), so d/dc of
        # Δf = 6x₁ + 6 is 6x₁/c for the cubic below.
        g₃(x) = x[1]^3 + x[1] * x[3]^2 - 3 * x[1] * x[2]^2    # Δ = 6x₁ + 2x₁ - 6x₁ = 2x₁
        for (mkchart, mkCS, ξ) in (
                (c -> ProlateSpheroidalChart(c), c -> coorsys_spheroidal(PointwiseGeometry, c), [0.3, 0.4, 1.7]),
                (c -> OblateSpheroidalChart(c), c -> coorsys_spheroidal_oblate(PointwiseGeometry, c), [0.3, 0.4, 0.7]),
            )
            lap(c) = laplacian(ζ -> g₃(coords_to_cartesian(mkchart(c), ζ)), mkCS(c))(ξ)
            x₁(c) = coords_to_cartesian(mkchart(c), ξ)[1]
            c₀ = 1.5
            @test lap(c₀) ≈ 2 * x₁(c₀) rtol = 1.0e-9
            @test ForwardDiff.derivative(lap, c₀) ≈ 2 * x₁(c₀) / c₀ rtol = 1.0e-8
        end
        # ellipsoidal: x₁ = ρμν/(h₂h₃) so ∂x₁/∂h₂ = -x₁/h₂ and ∂x₁/∂h₃ = -x₁/h₃
        ξ = [2.5, 1.3, -0.4]
        lapE(h₂, h₃) = laplacian(
            ζ -> g₃(coords_to_cartesian(EllipsoidalChart(h₂, h₃), ζ)), coorsys_ellipsoidal(PointwiseGeometry, h₂, h₃)
        )(ξ)
        x₁E = ξ[1] * ξ[2] * ξ[3] / 2.0
        @test lapE(2.0, 1.0) ≈ 2 * x₁E rtol = 1.0e-9
        @test ForwardDiff.derivative(h -> lapE(h, 1.0), 2.0) ≈ -2 * x₁E / 2.0 rtol = 1.0e-8
        @test ForwardDiff.derivative(h -> lapE(2.0, h), 1.0) ≈ -2 * x₁E / 1.0 rtol = 1.0e-8

        # element types and errors
        @test coorsys_spheroidal(PointwiseGeometry, 2) isa CoordinateSystem{3, Float64, <:PointwiseGeometry}
        @test coorsys_spheroidal_oblate(PointwiseGeometry, 1.5f0) isa CoordinateSystem{3, Float32, <:PointwiseGeometry}
        @test_throws DomainError coorsys_spheroidal(PointwiseGeometry, -1.0)
        @test_throws DomainError coorsys_spheroidal_oblate(PointwiseGeometry, 0.0)
        @test_throws DomainError coorsys_ellipsoidal(PointwiseGeometry, 1.0, 2.0)
    end

    @testsection "Symbolics" begin
        for (mk, cs, ξ₀) in (
                (coorsys_spheroidal, :prolate, (0.3, 0.4, 1.7)),
                (coorsys_spheroidal_oblate, :oblate, (0.3, 0.4, 0.7)),
            )
            CS = mk(Symbolics.Num)
            @test CS isa CoordinateSystem{3, Symbolics.Num}
            ϕn, pn, qn = get_coords(CS)
            OM = position_vector(CS)
            𝐟 = OM[1]^2 - OM[3]^2 + OM[2] * OM[3]
            Δ = divergence(grad(𝐟, CS), CS)
            cn = Symbolics.@variables(c)[1]
            d = Dict(cn => 1.5, ϕn => ξ₀[1], pn => ξ₀[2], qn => ξ₀[3])
            val = Symbolics.value(Symbolics.substitute(Δ, d; fold = Val(true)))
            @test abs(Float64(val)) < 1.0e-10
            # ∇²(x₁²) = 2
            val = Symbolics.value(Symbolics.substitute(divergence(grad(OM[1]^2, CS), CS), d; fold = Val(true)))
            @test Float64(val) ≈ 2 rtol = 1.0e-10
        end
        CS = coorsys_ellipsoidal(Symbolics.Num)
        @test CS isa CoordinateSystem{3, Symbolics.Num}
        ρn, μn, νn = get_coords(CS)
        OM = position_vector(CS)
        Symbolics.@variables h₂ h₃
        d = Dict(h₂ => 2.0, h₃ => 1.0, ρn => 2.5, μn => 1.3, νn => -0.4)
        Δ = divergence(grad(OM[1]^2 - OM[3]^2, CS), CS)
        @test abs(Float64(Symbolics.value(Symbolics.substitute(Δ, d; fold = Val(true))))) < 1.0e-9
        Δ = divergence(grad(OM[2] * OM[3], CS), CS)
        @test abs(Float64(Symbolics.value(Symbolics.substitute(Δ, d; fold = Val(true))))) < 1.0e-9
    end

    @testsection "tsubs on a RotatedBasis" begin
        # SymPy
        S = coorsys_spheroidal()
        ℬ = normalized_basis(S)
        @test ℬ isa RotatedBasis
        ℬ₀ = tsubs(ℬ, ϕ => 0.3, p => 0.4, q => 1.7)
        @test ℬ₀ isa RotatedBasis
        R = TensND._frame(ProlateSpheroidalChart(1.5), (0.3, 0.4, 1.7))
        @test [Float64(ℬ₀.eᵢ[i, j]) for i in 1:3, j in 1:3] ≈ Matrix(R) rtol = 1.0e-12
        Rn = [Float64(ℬ₀.eᵢ[i, j]) for i in 1:3, j in 1:3]
        @test Rn' * Rn ≈ I(3)
        # the identity gives back the canonical basis
        Sp = coorsys_spherical()
        θ, φ, r = get_coords(Sp)
        ℬˢ = tsubs(normalized_basis(Sp), θ => 0, φ => 0)
        @test ℬˢ isa CanonicalBasis
        @test tsubs(normalized_basis(Sp), θ => 0.7, φ => 0.2) isa RotatedBasis

        # Symbolics
        Sn = coorsys_spheroidal(Symbolics.Num)
        ℬn = normalized_basis(Sn)
        ϕn, pn, qn = get_coords(Sn)
        ℬn₀ = tsubs(ℬn, ϕn => 0.3, pn => 0.4, qn => 1.7)
        @test ℬn₀ isa RotatedBasis
        fold(x) = Float64(Symbolics.value(Symbolics.substitute(x, Dict(); fold = Val(true))))
        @test [fold(ℬn₀.eᵢ[i, j]) for i in 1:3, j in 1:3] ≈ Matrix(R) rtol = 1.0e-12
        # No simplification of `cos(0)` on a `Num`: the identity is not recognized
        # as such, but it is the identity once evaluated.
        Spn = coorsys_spherical(Symbolics.Num)
        θn, φn, rn = get_coords(Spn)
        ℬˢn = tsubs(normalized_basis(Spn), θn => 0, φn => 0)
        @test ℬˢn isa RotatedBasis
        @test [fold(ℬˢn.eᵢ[i, j]) for i in 1:3, j in 1:3] ≈ I(3) atol = 1.0e-14
        @test tsubs(normalized_basis(Spn), θn => 0.7, φn => 0.2) isa RotatedBasis
    end
end
