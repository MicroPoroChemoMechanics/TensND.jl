# The toroidal system built from its chart. Exact (SymPy, Symbolics) and
# pointwise versions, against the closed form of the Laplacian in toroidal
# coordinates, against each other, and against the exact Cartesian operators.

import Symbolics
using StaticArrays: SVector

# Numeric value of a SymPy expression at a point.
_tor_val(e, d) = Float64(e.subs(d))

# Components of a Cartesian vector field on the unit vectors of the chart.
_tor_local(chart, ξ, V) = TensND._frame(chart, ξ)' * V(coords_to_cartesian(chart, ξ))

@testsection "Toroidal coordinate system" begin

    a = symbols("a", positive = true)
    d̄ = symbols("d̄", positive = true)

    @testsection "Exact system (SymPy)" begin
        T = coorsys_toroidal()
        ϕ, α, β = get_coords(T)
        OM = position_vector(T)
        χ = lame_coefficients(T)
        d = cosh(α) - cos(β)
        closed = (a * sinh(α) / d, a / d, a / d)
        for i in 1:3
            @test iszero(tsimplify(χ[i] - closed[i]))
        end
        # harmonic polynomials of the Cartesian components
        @test tsimplify(laplacian(OM[1]^2, T)) == 2
        @test iszero(tsimplify(laplacian(OM[1]^2 - OM[3]^2, T)))
        @test iszero(tsimplify(laplacian(OM[1] * OM[3], T)))
        @test iszero(tsimplify(laplacian(OM[2] * OM[3], T)))
        # the Laplacian of a function of α alone, and of β alone:
        # Δf = d³/(a² sinh α) [∂α(sinh α/d ∂αf) + ∂β(sinh α/d ∂βf)]
        f = SymFunction("f", real = true)(α)
        @test iszero(
            tsimplify(laplacian(f, T) - d^3 / (a^2 * sinh(α)) * ∂(sinh(α) / d * ∂(f, α), α))
        )
        g = SymFunction("g", real = true)(β)
        @test iszero(
            tsimplify(laplacian(g, T) - d^3 / (a^2 * sinh(α)) * ∂(sinh(α) / d * ∂(g, β), β))
        )
        # 1/‖x‖ is harmonic, and ‖x‖² = a² (cosh α + cos β)/(cosh α - cos β)
        r² = a^2 * (cosh(α) + cos(β)) / d
        @test iszero(tsimplify((OM[1]^2 + OM[2]^2 + OM[3]^2 - r²).subs(d̄, d)))
        @test iszero(tsimplify(laplacian(1 / sqrt(r²), T)))
        # the position vector has the components of the chart
        ξ₀ = (0.3, 0.4, 1.7)
        at = Dict(a => 1.5, ϕ => ξ₀[1], α => ξ₀[2], β => ξ₀[3], d̄ => cosh(ξ₀[2]) - cos(ξ₀[3]))
        @test [_tor_val(OM[i], at) for i in 1:3] ≈ coords_to_cartesian(ToroidalChart(1.5), ξ₀) rtol = 1.0e-12
        # the unit vectors are those of the chart
        ℬ = normalized_basis(T)
        @test [_tor_val(ℬ[i, j], at) for i in 1:3, j in 1:3] ≈ Matrix(TensND._frame(ToroidalChart(1.5), ξ₀)) rtol = 1.0e-12
        # coordinates given by the caller
        u, v, w = symbols("u", real = true), symbols("v", positive = true), symbols("w", real = true)
        b = symbols("b", positive = true)
        T′ = coorsys_toroidal((u, v, w), b)
        @test get_coords(T′) == (u, v, w)
        @test iszero(tsimplify(lame_coefficients(T′)[2] - b / (cosh(v) - cos(w))))
        @test coorsys_toroidal(Sym) isa CoordinateSystem{3, <:Sym}
    end

    @testsection "Christoffel symbols and Lamé coefficients: exact against pointwise" begin
        T = coorsys_toroidal()
        ϕ, α, β = get_coords(T)
        for ξ₀ in ([0.3, 0.4, 1.7], [-2.1, 1.6, -2.5])
            at = Dict(a => 1.5, ϕ => ξ₀[1], α => ξ₀[2], β => ξ₀[3], d̄ => cosh(ξ₀[2]) - cos(ξ₀[3]))
            Γsym = map(γ -> _tor_val(γ, at), christoffel(T))
            CS = coorsys_toroidal(PointwiseGeometry, 1.5)
            @test christoffel(CS, ξ₀) ≈ Γsym rtol = 1.0e-10 atol = 1.0e-12
            χsym = [_tor_val(χ, at) for χ in lame_coefficients(T)]
            @test collect(lame_coefficients(CS, ξ₀)) ≈ χsym rtol = 1.0e-12
        end
    end

    @testsection "Pointwise system: operators against the Cartesian ones" begin
        f(x) = x[1]^3 + x[1] * x[2] * x[3] + 2 * x[3]^2 + x[2]^2
        ∇f(x) = [3x[1]^2 + x[2] * x[3], x[1] * x[3] + 2x[2], x[1] * x[2] + 4x[3]]
        Δf(x) = 6x[1] + 2 + 4
        V(x) = [x[1]^2, x[2] * x[3], x[1] * x[3]^2]
        divV(x) = 2x[1] + x[3] + 2x[1] * x[3]
        h(x) = x[1]^2 - x[3]^2 + x[2] * x[3]

        chart = ToroidalChart(1.5)
        CS = coorsys_toroidal(PointwiseGeometry, 1.5)
        for ξ in ([0.3, 0.4, 1.7], [-2.1, 1.6, -2.5], [1.2, 0.9, 0.6], [2.8, 2.2, -1.1])
            x = coords_to_cartesian(chart, ξ)
            F = ζ -> f(coords_to_cartesian(chart, ζ))
            @test laplacian(F, CS)(ξ) ≈ Δf(x) rtol = 1.0e-8
            @test laplacian(ζ -> h(coords_to_cartesian(chart, ζ)), CS)(ξ) ≈ 0 atol = 1.0e-8
            g = grad(F, CS)(ξ)
            @test g isa AbstractTens
            @test collect(canonical_components(g)) ≈ ∇f(x) rtol = 1.0e-9
            W = ζ -> _tor_local(chart, ζ, V)
            @test divergence(W, CS)(ξ) ≈ divV(x) rtol = 1.0e-8
            # the generic system built from the position vector agrees
            CSg = CoordinateSystem(ζ -> Vector(coords_to_cartesian(chart, ζ)), 3)
            @test laplacian(F, CSg)(ξ) ≈ laplacian(F, CS)(ξ) rtol = 1.0e-8
            @test collect(canonical_components(grad(F, CSg)(ξ))) ≈ collect(canonical_components(g)) rtol = 1.0e-9
            @test divergence(W, CSg)(ξ) ≈ divergence(W, CS)(ξ) rtol = 1.0e-8
            @test collect(lame_coefficients(CSg, ξ)) ≈ collect(lame_coefficients(CS, ξ)) rtol = 1.0e-10
            @test christoffel(CSg, ξ) ≈ christoffel(CS, ξ) rtol = 1.0e-7 atol = 1.0e-9
            @test Array(normalized_basis(CS, ξ)) ≈ Matrix(TensND._frame(chart, ξ)) rtol = 1.0e-12
        end

        # derivative of a pointwise Laplacian with respect to the focal radius:
        # at fixed coordinates the point is x(ξ, a) = a x(ξ, 1), so d/da of
        # Δg₃ = 2x₁ is 2x₁/a.
        g₃(x) = x[1]^3 + x[1] * x[3]^2 - 3 * x[1] * x[2]^2
        ξ = [0.3, 0.4, 1.7]
        lap(a) = laplacian(ζ -> g₃(coords_to_cartesian(ToroidalChart(a), ζ)), coorsys_toroidal(PointwiseGeometry, a))(ξ)
        x₁(a) = coords_to_cartesian(ToroidalChart(a), ξ)[1]
        @test lap(1.5) ≈ 2 * x₁(1.5) rtol = 1.0e-8
        @test ForwardDiff.derivative(lap, 1.5) ≈ 2 * x₁(1.5) / 1.5 rtol = 1.0e-7

        # element types and errors
        @test coorsys_toroidal(PointwiseGeometry, 2) isa CoordinateSystem{3, Float64, <:PointwiseGeometry}
        @test coorsys_toroidal(PointwiseGeometry, 1.5f0) isa CoordinateSystem{3, Float32, <:PointwiseGeometry}
        @test_throws DomainError coorsys_toroidal(PointwiseGeometry, -1.0)
        @test_throws DomainError coorsys_toroidal(PointwiseGeometry, 0.0)
    end

    @testsection "Symbolics" begin
        CS = coorsys_toroidal(Symbolics.Num)
        @test CS isa CoordinateSystem{3, Symbolics.Num}
        ϕn, αn, βn = get_coords(CS)
        OM = position_vector(CS)
        an = Symbolics.@variables(a)[1]
        d = Dict(an => 1.5, ϕn => 0.3, αn => 0.4, βn => 1.7)
        val(e) = Float64(Symbolics.value(Symbolics.substitute(e, d; fold = Val(true))))
        𝐟 = OM[1]^2 - OM[3]^2 + OM[2] * OM[3]
        @test abs(val(divergence(grad(𝐟, CS), CS))) < 1.0e-9
        @test val(divergence(grad(OM[1]^2, CS), CS)) ≈ 2 rtol = 1.0e-9
        χ = lame_coefficients(CS)
        @test [val(χ[i]) for i in 1:3] ≈ collect(TensND._lame(ToroidalChart(1.5), (0.3, 0.4, 1.7))) rtol = 1.0e-12
    end
end
