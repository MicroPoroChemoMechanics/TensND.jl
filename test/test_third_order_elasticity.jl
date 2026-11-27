# =============================================================================
#  test_third_order_elasticity.jl — the solid isotropic at orders two and
#  three (`ThirdOrderIsotropic`, `strain_energy`, `cauchy_stress`,
#  `acoustoelastic_tensor`, `prestressed_waves`, `third_order_constants`).
#
#  References, independent of the closed forms under test:
#  - the tangent A = ∂²W/∂F∂F and the first Piola-Kirchhoff stress ∂W/∂F by
#    ForwardDiff on the energy written as a function of F;
#  - the third derivatives of the energy with respect to the strain, by nested
#    ForwardDiff, for the constants C₁₁₁, C₁₁₂, C₁₂₃;
#  - the formulas of Hughes & Kelly (1953), Phys. Rev. 92(5), for the true
#    velocities under a hydrostatic pressure p, with K = λ + 2μ/3:
#        ρ₀ v_P² = λ + 2μ − p (6l + 4m + 7λ + 10μ)/(3K) ,
#        ρ₀ v_S² = μ − p (3m − n/2 + 3λ + 6μ)/(3K) ,
#    the true velocity being (1 + e) times the natural one for F = (1 + e)𝟏,
#    and p = −3Ke at the first order;
#  - the identity g·N = V between natural group and phase velocities.
# =============================================================================

using StaticArrays: SMatrix

const TOE_I3 = Matrix(1.0I, 3, 3)
# Lamé's constants and Murnaghan's (GPa), of the order of those of a sandstone.
const TOE_CONSTANTS = (8.3, 11.9, -3.4e3, -6.8e3, -6.7e3)
const TOE_F = TOE_I3 + 1.0e-3 * [0.3 -0.2 0.5; 0.1 0.4 -0.3; -0.6 0.2 0.1]

toe_energy(c, F) = strain_energy(ThirdOrderIsotropic(c...), (transpose(F) * F - I) / 2)

@testsection "third-order elasticity" begin
    s = ThirdOrderIsotropic(TOE_CONSTANTS...)
    λ, μ, l, m, n = TOE_CONSTANTS

    @testsection "constructor and notations" begin
        @test ThirdOrderIsotropic(1, 2.0, 3, 4, 5) isa ThirdOrderIsotropic{Float64}
        c = third_order_constants(s)
        W(e) = strain_energy(s, [e[1] e[6] e[5]; e[6] e[2] e[4]; e[5] e[4] e[3]])
        d3(i, j, k) = ForwardDiff.derivative(
            a -> ForwardDiff.derivative(
                b -> ForwardDiff.derivative(t -> W([(q == i) * t + (q == j) * b + (q == k) * a for q in 1:6]), 0.0), 0.0
            ), 0.0
        )
        @test c.C111 ≈ d3(1, 1, 1)
        @test c.C112 ≈ d3(1, 1, 2)
        @test c.C123 ≈ d3(1, 2, 3)
        # Landau's energy μ tr E² + ½λ (tr E)² + ⅓A tr E³ + B tr E tr E² + ⅓C (tr E)³
        E = (transpose(TOE_F) * TOE_F - I) / 2
        @test strain_energy(s, E) ≈ μ * tr(E^2) + λ * tr(E)^2 / 2 + c.A * tr(E^3) / 3 + c.B * tr(E) * tr(E^2) + c.C * tr(E)^3 / 3
    end

    @testsection "stress and tangent against the derivatives of the energy" begin
        A = acoustoelastic_tensor(s, TOE_F)
        H = ForwardDiff.hessian(f -> toe_energy(TOE_CONSTANTS, reshape(f, 3, 3)), vec(TOE_F))
        @test reshape(A, 9, 9) ≈ H rtol = 1.0e-12
        @test A ≈ permutedims(A, (3, 4, 1, 2))
        @test !(A ≈ permutedims(A, (2, 1, 3, 4)))
        P = reshape(ForwardDiff.gradient(f -> toe_energy(TOE_CONSTANTS, reshape(f, 3, 3)), vec(TOE_F)), 3, 3)
        σ = cauchy_stress(s, TOE_F)
        @test σ ≈ P * transpose(TOE_F) / det(TOE_F) rtol = 1.0e-12
        @test σ ≈ transpose(σ)
        @test @inferred(cauchy_stress(s, SMatrix{3, 3}(TOE_F))) ≈ σ
        # the tangent of the second Piola-Kirchhoff stress: Hessian of the energy in
        # the six strains (engineering shears), and its part linear in the strain
        E = (transpose(TOE_F) * TOE_F - I) / 2
        strain(e) = [e[1] e[6] / 2 e[5] / 2; e[6] / 2 e[2] e[4] / 2; e[5] / 2 e[4] / 2 e[3]]
        H₆ = ForwardDiff.hessian(e -> strain_energy(s, strain(e)), [E[1, 1], E[2, 2], E[3, 3], 2E[2, 3], 2E[1, 3], 2E[1, 2]])
        @test voigt_stiffness(@inferred(tangent_stiffness(s, E))) ≈ H₆ rtol = 1.0e-12
        Δ(X) = canonical_components(tangent_stiffness(s, X)) - canonical_components(tangent_stiffness(s, zero(X)))
        @test Δ(2.5 * E) ≈ 2.5 * Δ(E)
        @test canonical_components(tangent_stiffness(s, zero(E))) ≈ get_array(TensISO{3}(3λ + 2μ, 2μ))
        ℂ = TensISO{3}(3λ + 2μ, 2μ)
        @test acoustoelastic_tensor(s, TOE_I3) ≈ get_array(ℂ)
        @test iszero(cauchy_stress(s, TOE_I3))
    end

    @testsection "waves under a hydrostatic pressure (Hughes & Kelly)" begin
        K = λ + 2μ / 3
        ρ₀ = 2.3
        N = [0.3, -0.5, 0.7]
        true_moduli(e) = (w = prestressed_waves(s, ρ₀, (1 + e) * TOE_I3, N); ρ₀ * (1 + e)^2 * w.velocity .^ 2)
        M₀ = true_moduli(0.0)
        @test M₀ ≈ [λ + 2μ, μ, μ]
        dM = ForwardDiff.derivative(true_moduli, 0.0)          # dp = −3K de
        @test dM[1] ≈ 6l + 4m + 7λ + 10μ
        @test dM[2] ≈ 3m - n / 2 + 3λ + 6μ
        @test dM[3] ≈ dM[2]
        # a finite pressure: the first-order formula holds to the second order
        e = -1.0e-5
        @test true_moduli(e)[1] ≈ λ + 2μ + e * (6l + 4m + 7λ + 10μ) rtol = 1.0e-5
    end

    @testsection "any deformation: eigenpairs, group velocity, derivatives" begin
        ρ₀ = 2.3
        N = [0.2, -0.5, 0.8]
        w = prestressed_waves(s, ρ₀, TOE_F, N)
        A = acoustoelastic_tensor(s, TOE_F)
        ν = N / norm(N)
        Γ = [sum(ν[j] * A[i, j, k, q] * ν[q] for j in 1:3, q in 1:3) for i in 1:3, k in 1:3]
        @test issorted(w.velocity; rev = true)
        for k in 1:3
            @test Γ * w.polarization[:, k] ≈ ρ₀ * w.velocity[k]^2 * w.polarization[:, k]
            @test dot(w.group[:, k], w.direction) ≈ w.velocity[k]
        end
        # the group velocity is the gradient of the frequency ω(k) = |k| V(k/|k|)
        F₁ = TOE_I3 + 2.0e-3 * [0.3 -0.2 0.5; 0.1 0.4 -0.3; -0.6 0.2 0.1]     # the stress breaks the minor symmetries
        w₁ = prestressed_waves(s, ρ₀, F₁, N)
        for k in 1:3
            ω(κ) = norm(κ) * prestressed_waves(s, ρ₀, F₁, κ).velocity[k]
            @test ForwardDiff.gradient(ω, ν) ≈ w₁.group[:, k] rtol = 1.0e-9
        end
        # derivative of the fastest velocity with respect to l, by the Rayleigh quotient
        dV = ForwardDiff.derivative(x -> prestressed_waves(ThirdOrderIsotropic(λ, μ, x, m, n), ρ₀, TOE_F, N).velocity[1], l)
        dA = (acoustoelastic_tensor(ThirdOrderIsotropic(0.0, 0.0, 1.0, 0.0, 0.0), TOE_F))
        U = w.polarization[:, 1]
        @test dV ≈ sum(U[i] * ν[j] * dA[i, j, k, q] * U[k] * ν[q] for i in 1:3, j in 1:3, k in 1:3, q in 1:3) / (2ρ₀ * w.velocity[1])
        # nested: second derivative with respect to the stretch of a hydrostatic state
        V(e) = prestressed_waves(s, ρ₀, (1 + e) * TOE_I3, [1.0, 0.0, 0.0]).velocity[1]
        A1111(e) = acoustoelastic_tensor(s, (1 + e) * TOE_I3)[1, 1, 1, 1]
        d2(f) = ForwardDiff.derivative(a -> ForwardDiff.derivative(f, a), 0.0)
        @test d2(e -> ρ₀ * V(e)^2) ≈ d2(A1111)
        # a deformation that makes the tangent lose its definiteness
        @test_throws DomainError prestressed_waves(s, ρ₀, 1.05 * TOE_I3, N)
    end

    @testsection "uniaxial stress: closed-form coefficients" begin
        # (K, μ, l, m, n) in GPa and the predicted (b_S, a_P) in GPa⁻¹ of Table 1 of
        # Johnson & Rasolofosaon (1996), J. Geophys. Res. 101(B2), 3113-3124, printed
        # with two or three digits.
        table = (
            (164.0, 82.0, -348.0, -1030.0, 1100.0, 0.027, -0.032),      # Armco iron
            (31.8, 27.5, 14.0, 92.0, 420.0, 0.088, 0.092),              # Pyrex
            (13.8, 18.2, -3371.0, -6742.0, -6600.0, -2.47, -9.67),      # Barre granite
        )
        e₁, e₃ = [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]
        for (K, G, tl, tm, tn, b_S, a_P) in table
            L = K - 2G / 3
            t = ThirdOrderIsotropic(L, G, tl, tm, tn)
            c = acoustoelastic_coefficients(t)
            @test c.b_S ≈ b_S atol = 6.0e-4 + 5.0e-3 * abs(b_S)
            @test c.a_P ≈ a_P atol = 6.0e-4 + 5.0e-3 * abs(a_P)
            # a compression σ along e₃ at the first order: F = 𝟏 + 𝕊 : (−σ e₃ ⊗ e₃)
            E, ν = 9K * G / (3K + G), (3K - 2G) / (2 * (3K + G))
            F(σ) = [1 + ν * σ / E 0 0; 0 1 + ν * σ / E 0; 0 0 1 - σ / E]
            # ρ₀v² of the wave of direction N polarized along U at σ = 0 (true velocity)
            function slope(N, U)
                ForwardDiff.derivative(0.0) do σ
                    A = acoustoelastic_tensor(t, F(σ))
                    sum(U[i] * N[j] * A[i, j, k, q] * U[k] * N[q] for i in 1:3, j in 1:3, k in 1:3, q in 1:3) * norm(F(σ) * N)^2
                end
            end
            @test slope(e₃, e₃) ≈ c.P_parallel
            @test slope(e₃, e₁) ≈ c.S_parallel
            @test slope(e₁, e₁) ≈ c.P_perpendicular
            @test slope(e₁, e₃) ≈ c.S_perpendicular_parallel
            @test slope(e₁, [0.0, 1.0, 0.0]) ≈ c.S_perpendicular_perpendicular
            @test c.a_P ≈ (c.P_perpendicular - c.P_parallel) / (2 * (L + 2G))
            @test c.b_S ≈ (c.S_perpendicular_perpendicular - c.S_perpendicular_parallel) / (2G)
            @test c.P_hydrostatic ≈ -(6tl + 4tm + 7L + 10G) / (3K)
            @test c.S_hydrostatic ≈ -(3tm - tn / 2 + 3L + 6G) / (3K)
            # the anisotropy of P at an angle θ from the stress goes as sin²θ
            θ = 0.6
            N = [sin(θ), 0.0, cos(θ)]
            vP(σ, N) = norm(F(σ) * N) * prestressed_waves(t, 1.0, F(σ), N).velocity[1]
            A_P = ForwardDiff.derivative(σ -> (vP(σ, e₃) - vP(σ, N)) / vP(σ, e₃), 0.0)
            @test A_P ≈ -c.a_P * sin(θ)^2
        end
    end

    @testsection "symbols" begin
        @syms λs μs ls ms ns es
        ss = ThirdOrderIsotropic(λs, μs, ls, ms, ns)
        Fs = [i == j ? 1 + es : zero(es) for i in 1:3, j in 1:3]
        A = acoustoelastic_tensor(ss, Fs)
        # natural moduli A₁₁₁₁ and A₂₁₂₁: true ones after the factor (1 + e)²
        dP = diff((1 + es)^2 * A[1, 1, 1, 1], es).subs(es, 0)
        dS = diff((1 + es)^2 * A[2, 1, 2, 1], es).subs(es, 0)
        @test iszero(simplify(dP - (6ls + 4ms + 7λs + 10μs)))
        @test iszero(simplify(dS - (3ms - ns / 2 + 3λs + 6μs)))
        σ = cauchy_stress(ss, Fs)
        @test iszero(simplify(diff(σ[1, 1], es).subs(es, 0) - (3λs + 2μs)))
        k = acoustoelastic_coefficients(ss)
        @test iszero(simplify(k.P_hydrostatic + (6ls + 4ms + 7λs + 10μs) / (3λs + 2μs)))
        @test iszero(simplify(k.S_hydrostatic + (3ms - ns / 2 + 3λs + 6μs) / (3λs + 2μs)))
        @test iszero(simplify(k.b_S - (k.S_perpendicular_perpendicular - k.S_perpendicular_parallel) / (2μs)))
        c = third_order_constants(ss)
        @test iszero(simplify(c.C111 - c.C112 - 4ms))
    end
end
