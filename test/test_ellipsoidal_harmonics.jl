# Lamé functions and ellipsoidal harmonics. References: the closed forms of the
# low degrees (Dassios 2012, Eqs. (3.37), (3.38), (3.95), (3.96)), harmonicity
# checked by automatic differentiation, the degree 1 in closed form, and the
# normalization constants of Kushch (2020), Table 1, p. 8 (ellipsoid of
# semi-axes (1.5, 1, 0.5)), computed here by a quadrature of the surface
# harmonics independent of the code under test.

# The chart of the ellipsoid of semi-axes a₁ > a₂ > a₃.
_eh_chart(a₁, a₂, a₃) = EllipsoidalChart(sqrt(a₁^2 - a₃^2), sqrt(a₁^2 - a₂^2))

# Kushch (2020), Table 1: (n, γ) for the ellipsoid (1.5, 1, 0.5).
const _EH_KUSHCH_2020 = ((1, 10.471975), (1, 6.2831853), (5, 0.021145271), (5, 0.050203007), (10, 4.1852909e-5), (10, 2.1196723e-4))

@testsection "Ellipsoidal harmonics" begin
    a = (1.3, 0.9, 0.5)
    chart = _eh_chart(a...)
    h₂, h₃ = chart.h₂, chart.h₃
    h₁ = sqrt(h₂^2 - h₃^2)

    @testsection "Lamé functions: classes and closed forms" begin
        for n in 0:6
            E = lame_functions(chart, n)
            @test length(E) == 2n + 1
            @test all(e -> e.coefficients[1] == 1, E)
        end
        # degree 2, class K: x² + Λ - a₁², Λ the roots of Σ 1/(Λ - aᵢ²) = 0
        a² = a .^ 2
        s₁, s₂ = sum(a²), a²[1] * a²[2] + a²[2] * a²[3] + a²[3] * a²[1]
        Λ = sort([(s₁ + sqrt(s₁^2 - 3s₂)) / 3, (s₁ - sqrt(s₁^2 - 3s₂)) / 3])
        @test sort([E.coefficients[2] for E in lame_functions(chart, 2)[1:2]]) ≈ Λ .- a²[1] rtol = 1.0e-13
        # class N: p = 1 at degree 2 and 4 at degree 3 (Dassios, § 3.4)
        @test only(E.p for E in lame_functions(chart, 2) if E.a == E.b == 1) ≈ 1 rtol = 1.0e-13
        @test only(E.p for E in lame_functions(chart, 3) if E.a == E.b == 1) ≈ 4 rtol = 1.0e-13
        @test_throws DomainError lame_functions(chart, -1)
        @test_throws DomainError lame_value(first(lame_functions(chart, 1)), ellipsoidal_signed_coordinates(chart, [0.1, 0.2, 0.3]), 4)
    end

    @testsection "Harmonic inside and outside" begin
        nodes = TensND._gauss_legendre_01(30)
        interior(E, x) = (c = ellipsoidal_signed_coordinates(chart, x); prod(first(lame_value(E, c, i)) for i in 1:3))
        exterior(E, x) = (c = ellipsoidal_signed_coordinates(chart, x); first(lame_second_kind(chart, E, c.ρ, nodes)) * first(lame_value(E, c, 2)) * first(lame_value(E, c, 3)))
        # points in several octants
        for x in ([0.7, -0.4, 0.3], [-0.5, 0.6, -0.2]), n in 0:3, E in lame_functions(chart, n)
            @test abs(tr(ForwardDiff.hessian(z -> interior(E, z), x))) ≤ 1.0e-11
        end
        for y in ([1.9, 0.8, -1.2], [-1.4, -1.1, 0.9]), n in 0:3, E in lame_functions(chart, n)
            @test abs(tr(ForwardDiff.hessian(z -> exterior(E, z), y))) ≤ 1.0e-12
        end
        # degree 1: the three coordinates, x₁ = ρμν/(h₂h₃), and their decaying counterparts
        x = [0.7, -0.4, 0.3]
        K1, L1, M1 = lame_functions(chart, 1)
        @test interior(K1, x) ≈ h₂ * h₃ * x[1] rtol = 1.0e-12
        @test interior(L1, x) ≈ h₁ * h₃ * x[2] rtol = 1.0e-12
        @test interior(M1, x) ≈ h₁ * h₂ * x[3] rtol = 1.0e-12
        far = 300 .* [0.3, 0.5, -0.8]
        @test exterior(K1, far) ≈ h₂ * h₃ * far[1] / norm(far)^3 rtol = 1.0e-4
        # the derivative of the second kind against the derivative of its integral
        E = lame_functions(chart, 2)[1]
        @test ForwardDiff.derivative(ρ -> first(lame_second_kind(chart, E, ρ)), 1.7) ≈ last(lame_second_kind(chart, E, 1.7)) rtol = 1.0e-10
        # the default nodes, and the rule itself
        xs, ws = TensND._gauss_legendre_01(12)
        @test sum(ws) ≈ 1 && sum(ws .* xs .^ 11) ≈ 1 / 12
    end

    @testsection "Normalization: Kushch (2020), Table 1" begin
        # γ = ∫ (E(μ)E(ν))² dΩ on the surface, dΩ = (μ²-ν²) dμ dν / [√(μ²-h₃²)√(h₂²-μ²)√(h₃²-ν²)√(h₂²-ν²)],
        # in μ² = h₃² + (h₂²-h₃²) sin²α, α ∈ [0, π], and ν = h₃ sin β, β ∈ [0, 2π),
        # where dΩ = (μ²-ν²)/(μ √(h₂²-ν²)) dα dβ and the roots take the signs of
        # cos α and cos β
        c₂ = _eh_chart(1.5, 1.0, 0.5)
        H₂, H₃ = c₂.h₂, c₂.h₃
        H₁ = sqrt(H₂^2 - H₃^2)
        α, wα = TensND._gauss_legendre_01(80)
        nβ = 160
        γ(E) = sum(
            begin
                A, B = π * α[i], 2π * (j - 1) / nβ
                μ, ν = sqrt(H₃^2 + H₁^2 * sin(A)^2), H₃ * sin(B)
                c = (; μ, ν, μ3 = H₁ * sin(A), μ2 = H₁ * cos(A), ν3 = H₃ * cos(B), ν2 = sqrt(H₂^2 - ν^2))
                S = first(lame_value(E, c, 2)) * first(lame_value(E, c, 3))
                π * wα[i] * (2π / nβ) * S^2 * (μ^2 - ν^2) / (μ * sqrt(H₂^2 - ν^2))
            end for i in eachindex(α), j in 1:nβ
        )
        g₁ = γ.(lame_functions(c₂, 1))
        @test g₁ ≈ (4π / 3) .* [(H₂ * H₃)^2, (H₁ * H₃)^2, (H₁ * H₂)^2] rtol = 1.0e-12
        for (n, ref) in _EH_KUSHCH_2020
            @test minimum(abs.(γ.(lame_functions(c₂, n)) .- ref)) ≤ 1.0e-6 * ref
        end
    end
end
