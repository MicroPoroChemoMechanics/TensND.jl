# Toroidal functions: Legendre functions of half-integer degree and argument
# cosh α. References independent of the code under test: Laplace's integral
# for the first kind, Heine's for the second, both by the trapezoidal rule on a
# smooth periodic integrand; mpmath; the Casoratian; derivatives by ForwardDiff;
# the Laplace equation in the toroidal chart.

import Symbolics
using SymPyPythonCall: PythonCall

# The trapezoidal rule on [0, π] for an even 2π-periodic integrand converges
# faster than any power of the step.
_th_trapz(f, n = 4000) = (π / n) * (sum(f(π * k / n) for k in 1:(n - 1)) + (f(0.0) + f(π)) / 2)

# Laplace: P_ν(cosh α) = (1/π) ∫₀^π (cosh α + sinh α cos t)^ν dt.
_th_P(n, α) = _th_trapz(t -> (cosh(α) + sinh(α) * cos(t))^(n - 1 / 2)) / π
# Heine: Q_{n-1/2}(cosh α) = ∫₀^π cos(n t) / √(2 cosh α - 2 cos t) dt.
_th_Q(n, α) = _th_trapz(t -> cos(n * t) / sqrt(2 * cosh(α) - 2 * cos(t)), 20000)

@testsection "Toroidal functions" begin

    αs = (0.05, 0.3, 1.0, 2.5, 5.0)

    @testsection "Order 0 against the integrals of Laplace and Heine" begin
        N = 12
        for α in αs
            P = toroidal_P(N, 0, α)
            Q = toroidal_Q(N, 0, α)
            @test length(P) == N + 1
            @test length(Q) == N + 1
            for n in 0:N
                @test P[n + 1] ≈ _th_P(n, α) rtol = 1.0e-11
                # the integrand of Heine's integral is sharply peaked at small α
                α ≥ 0.3 && @test Q[n + 1] ≈ _th_Q(n, α) rtol = 1.0e-10 atol = 1.0e-13
            end
        end
        # a single degree
        @test toroidal_P(0, 0, 1.0) ≈ [_th_P(0, 1.0)] rtol = 1.0e-11
        @test toroidal_Q(0, 0, 1.0) ≈ [_th_Q(0, 1.0)] rtol = 1.0e-10
    end

    @testsection "Order 0 against mpmath" begin
        mp = PythonCall.pyimport("mpmath")
        mp.mp.dps = 30
        ref(f, n, α) = PythonCall.pyconvert(Float64, mp.re(f(mp.mpf(n) - mp.mpf(1) / 2, 0, mp.cosh(mp.mpf(α)), type = 3)))
        for α in (0.02, 0.05, 0.3, 1.0, 2.5, 5.0, 9.0), n in (0, 1, 2, 5, 11)
            @test toroidal_P(n, 0, α)[n + 1] ≈ ref(mp.legenp, n, α) rtol = 1.0e-12
            @test toroidal_Q(n, 0, α)[n + 1] ≈ ref(mp.legenq, n, α) rtol = 1.0e-12
        end
    end

    @testsection "Casoratian and limits" begin
        N = 30
        for α in αs
            P = toroidal_P(N, 0, α)
            Q = toroidal_Q(N, 0, α)
            # P_{ν+1} Q_ν - P_ν Q_{ν+1} = 1/(ν+1), ν = n - 1/2
            for n in 0:(N - 1)
                @test P[n + 2] * Q[n + 1] - P[n + 1] * Q[n + 2] ≈ 1 / (n + 1 / 2) rtol = 1.0e-9
            end
        end
        # on the axis every function of the first kind equals 1
        @test toroidal_P(5, 0, 0.0) ≈ ones(6)
        # far along α: Q_{-1/2} → π e^{-α/2}, P_{-1/2} → (2/π)(α + 2 ln 2) e^{-α/2}
        α = 30.0
        @test toroidal_Q(0, 0, α)[1] ≈ π * exp(-α / 2) rtol = 1.0e-12
        @test toroidal_P(0, 0, α)[1] ≈ 2 / π * (α + 2 * log(2)) * exp(-α / 2) rtol = 1.0e-10
        # the ratio of consecutive degrees tends to e^{-α}
        Q = toroidal_Q(80, 0, 2.0)
        @test Q[81] / Q[80] ≈ exp(-2.0) rtol = 2.0e-2
    end

    @testsection "Rounding: upward and downward runs against extended precision" begin
        setprecision(BigFloat, 256) do
            for α in (0.01, 0.02, 0.05, 0.5, 3.0, 8.0), m in (0, 1, 2)
                N = 60
                Qb = toroidal_Q(N, m, big(α))
                Pb = toroidal_P(N, m, big(α))
                @test eltype(Qb) == BigFloat
                @test maximum(abs.(toroidal_Q(N, m, α) ./ Float64.(Qb) .- 1)) < 1.0e-11
                # the order is raised by a quotient by sinh α: its rounding grows as 1/α²
                @test maximum(abs.(toroidal_P(N, m, α) ./ Float64.(Pb) .- 1)) < 1.0e-13 / min(α, 1)^(2 * min(m, 1) + 1)
            end
        end
        # both runs give the same values where the switch occurs
        α = 2 * log(2) / 20
        @test toroidal_Q(20, 0, α) ≈ toroidal_Q(21, 0, α)[1:21] rtol = 1.0e-12
    end

    @testsection "Orders: derivatives of the order 0 by ForwardDiff" begin
        # P_ν^m = wᵐ dᵐP_ν/dxᵐ with x = cosh α and w = sinh α
        N = 6
        for α in (0.3, 1.0, 2.5)
            x, w = cosh(α), sinh(α)
            for (f, name) in ((toroidal_P, :P), (toroidal_Q, :Q))
                g(y) = f(N, 0, acosh(y))
                d1 = ForwardDiff.derivative(g, x)
                d2 = ForwardDiff.derivative(y -> ForwardDiff.derivative(g, y), x)
                d3 = ForwardDiff.derivative(z -> ForwardDiff.derivative(y -> ForwardDiff.derivative(g, y), z), x)
                @test f(N, 1, α) ≈ w .* d1 rtol = 1.0e-9
                @test f(N, 2, α) ≈ w^2 .* d2 rtol = 1.0e-8
                @test f(N, 3, α) ≈ w^3 .* d3 rtol = 1.0e-7
            end
        end
        # the recurrence in the degree at the order m: (ν-m+1) F_{ν+1} = (2ν+1) x F_ν - (ν+m) F_{ν-1}
        for α in (0.3, 2.5), m in (1, 2, 3), f in (toroidal_P, toroidal_Q)
            F = f(10, m, α)
            for n in 1:9
                ν = n - 1 / 2
                @test (ν - m + 1) * F[n + 2] ≈ 2n * cosh(α) * F[n + 1] - (ν + m) * F[n] rtol = 1.0e-9
            end
        end
    end

    @testsection "Derivatives with respect to α" begin
        for α in (0.3, 1.0, 2.5), m in (0, 1, 2)
            @test toroidal_dP(8, m, α) ≈ ForwardDiff.derivative(a -> toroidal_P(8, m, a), α) rtol = 1.0e-10
            @test toroidal_dQ(8, m, α) ≈ ForwardDiff.derivative(a -> toroidal_Q(8, m, a), α) rtol = 1.0e-10
            # second derivatives, through nested duals
            @test ForwardDiff.derivative(a -> toroidal_dP(4, m, a), α) ≈
                ForwardDiff.derivative(b -> ForwardDiff.derivative(a -> toroidal_P(4, m, a), b), α) rtol = 1.0e-9
        end
        @test length(toroidal_dP(0, 0, 1.0)) == 1
        @test length(toroidal_dQ(0, 1, 1.0)) == 1
        # Legendre's equation in α: F″ + coth α F′ - [ν(ν+1) + m²/sinh²α] F = 0
        for f in (toroidal_P, toroidal_Q), m in (0, 2), α in (0.4, 1.7)
            F = f(5, m, α)
            dF = ForwardDiff.derivative(a -> f(5, m, a), α)
            d²F = ForwardDiff.derivative(b -> ForwardDiff.derivative(a -> f(5, m, a), b), α)
            for n in 0:5
                ν = n - 1 / 2
                @test d²F[n + 1] + coth(α) * dF[n + 1] ≈ (ν * (ν + 1) + m^2 / sinh(α)^2) * F[n + 1] rtol = 1.0e-8
            end
        end
    end

    @testsection "Harmonic functions of the toroidal chart" begin
        chart = ToroidalChart(1.3)
        CS = coorsys_toroidal(PointwiseGeometry, 1.3)
        for (kind, n, m) in ((toroidal_P, 0, 0), (toroidal_P, 2, 1), (toroidal_Q, 1, 0), (toroidal_Q, 3, 2))
            h = ζ -> sqrt(cosh(ζ[2]) - cos(ζ[3])) * kind(n, m, ζ[2])[n + 1] * cos(n * ζ[3] + 0.4) * cos(m * ζ[1] - 0.2)
            for ξ in ([0.3, 0.6, 1.1], [-1.2, 1.4, -2.0])
                scale = abs(h(ξ)) / lame_coefficients(chart, ξ)[2]^2
                @test abs(laplacian(h, CS)(ξ)) < 1.0e-7 * scale
            end
        end
    end

    @testsection "Types and errors" begin
        @inferred toroidal_P(4, 1, 0.7)
        @inferred toroidal_Q(4, 1, 0.7)
        @inferred toroidal_dP(4, 1, 0.7)
        @inferred toroidal_dQ(4, 1, 0.7)
        @inferred toroidal_Q(40, 2, 3.0)
        @test eltype(toroidal_P(3, 0, 0.7f0)) == Float32
        @test eltype(toroidal_Q(3, 1, 0.7f0)) == Float32
        @test eltype(toroidal_Q(3, 0, ForwardDiff.Dual(0.7, 1.0))) <: ForwardDiff.Dual
        @test_throws DomainError toroidal_P(3, 0, -0.1)
        @test_throws DomainError toroidal_P(3, 1, 0.0)
        @test_throws DomainError toroidal_dP(3, 0, 0.0)
        @test_throws DomainError toroidal_Q(3, 0, 0.0)
        @test_throws DomainError toroidal_dQ(3, 0, 0.0)
        @test_throws DomainError toroidal_Q(3, 0, -1.0)
        @test_throws DomainError toroidal_P(-1, 0, 1.0)
        @test_throws DomainError toroidal_Q(2, -1, 1.0)
        # an iteration that stops on a value: symbols are refused
        Symbolics.@variables αs
        @test_throws ArgumentError toroidal_P(2, 0, αs)
        @test_throws ArgumentError toroidal_Q(2, 0, αs)
        @test_throws ArgumentError toroidal_dP(2, 0, symbols("α", positive = true))
        @test_throws ArgumentError toroidal_dQ(2, 0, αs)
        # the complete elliptic integrals behind the seeds, at a known value
        K, KE = TensND._complete_elliptic(sqrt(0.5), sqrt(0.5))
        @test K ≈ 1.8540746773013719 rtol = 1.0e-14
        @test K - KE ≈ 1.3506438810476755 rtol = 1.0e-14
        @test_throws DomainError TensND._complete_elliptic(1.0, 0.0)
    end
end
