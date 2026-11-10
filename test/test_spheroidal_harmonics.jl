# Spheroidal harmonics and their re-expansion. References are independent of
# the code under test: closed forms of the Legendre functions, Gauss-Legendre
# quadrature built here from the Golub-Welsch eigenproblem, Neumann's integral
# for Q_n, the Casoratian of the two radial solutions, the terminating
# hypergeometric sum, the identities between the two sides of every
# re-expansion, and extended precision.

import Symbolics
using StaticArrays: SVector

# Gauss-Legendre rule of `n` nodes: eigenvalues and first eigenvector components
# of the Jacobi matrix.
function _sh_gauss(n::Integer)
    β = [k / sqrt(4k^2 - 1) for k in 1:(n - 1)]
    E = eigen(SymTridiagonal(zeros(n), β))
    return E.values, 2 .* E.vectors[1, :] .^ 2
end

# P_n(t) by Bonnet's recurrence, for real or complex t.
function _sh_legendre(n::Integer, t)
    a, b = one(t), t
    n == 0 && return a
    for k in 1:(n - 1)
        a, b = b, ((2k + 1) * t * b - k * a) / (k + 1)
    end
    return b
end

# Meridian-plane point of a chart, from the angular and the radial coordinates.
_sh_point(::SphericalChart, p, ζ) = SVector(acos(p), 0.6, ζ)
_sh_point(::Union{ProlateSpheroidalChart, OblateSpheroidalChart}, p, ζ) = SVector(0.6, p, ζ)

_sh_charts = (
    ProlateSpheroidalChart(1.0), ProlateSpheroidalChart(1.7),
    OblateSpheroidalChart(0.8), OblateSpheroidalChart(2.0), SphericalChart(),
)
_sh_pairs = [
    (to, from) for to in _sh_charts for from in _sh_charts
        if !(to isa SphericalChart && from isa SphericalChart)
]

# Largest |lhs - rhs| / (|lhs| + Σ|terms|) of the interior re-expansion over all
# degrees n ≤ N and orders m ≤ n at the Cartesian point `x`.
function _sh_interior_error(to, from, x, N)
    ξt, ξf = cartesian_to_coords(to, x), cartesian_to_coords(from, x)
    err = zero(eltype(x))
    for n in 0:N, m in 0:n
        L = interior_reexpansion(to, from, n, m)
        Ht = interior_harmonics(to, n, m, ξt)
        lhs = interior_harmonics(from, n, m, ξf)[n + 1]
        terms = [L[k + 1] * Ht[n - 2k + 1] for k in 0:((n - m) ÷ 2)]
        den = abs(lhs) + sum(abs, terms)
        err = max(err, iszero(den) ? zero(err) : abs(lhs - sum(terms)) / den)
    end
    return err
end

@testsection "Spheroidal harmonics" begin

    @testsection "legendre_P" begin
        for x in (-0.8, -0.1, 0.35, 0.9)
            s = sqrt(1 - x^2)
            @test legendre_P(3, 1, x) ≈ [0, -s, -3x * s, -1.5 * (5x^2 - 1) * s] rtol = 1.0e-13
            @test legendre_P(3, 2, x) ≈ [0, 0, 3 * (1 - x^2), 15x * (1 - x^2)] rtol = 1.0e-13
            @test legendre_P(3, 3, x) ≈ [0, 0, 0, -15 * s^3] rtol = 1.0e-13
            @test legendre_P(3, 0, x) ≈ [1, x, (3x^2 - 1) / 2, (5x^3 - 3x) / 2] rtol = 1.0e-13
            @test legendre_P(1, 2, x) == [0, 0]            # n < m
        end
        @test legendre_P(2, 1, 0.5) ≈ [0, -√0.75, -1.5 * √0.75]
        # orthogonality ∫P_n^m P_r^m = 2/(2n+1) (n+m)!/(n-m)! δ_nr
        t, w = _sh_gauss(24)
        for m in 0:6
            vals = [legendre_P(6, m, tj) for tj in t]
            nrm(n) = 2 / (2n + 1) * Float64(factorial(big(n + m)) / factorial(big(n - m)))
            for n in m:6, r in m:6
                I = sum(w[j] * vals[j][n + 1] * vals[j][r + 1] for j in eachindex(t))
                @test I ≈ (n == r ? nrm(n) : 0.0) atol = 1.0e-13 * sqrt(nrm(n) * nrm(r))
            end
        end
        # BigFloat agrees with Float64
        @test Float64.(legendre_P(6, 2, big"0.37")) ≈ legendre_P(6, 2, 0.37) rtol = 1.0e-13
    end

    @testsection "radial_P and radial_Q: closed forms" begin
        for q in (1.3, 2.5, 7.0)
            P = radial_P(ProlateSpheroidalChart(1.0), 2, 0, q)
            @test P ≈ [1, q, (3q^2 - 1) / 2] rtol = 1.0e-13
            @test radial_P(ProlateSpheroidalChart(1.0), 1, 1, q) ≈ [0, sqrt(q^2 - 1)] rtol = 1.0e-13
            @test radial_P(ProlateSpheroidalChart(1.0), 2, 2, q)[3] ≈ 3 * (q^2 - 1) rtol = 1.0e-13
            Q₀ = atanh(1 / q)                              # arccoth q
            Q = radial_Q(ProlateSpheroidalChart(1.0), 1, 0, q)
            @test Q ≈ [Q₀, q * Q₀ - 1] rtol = 1.0e-12
            Q₁¹ = radial_Q(ProlateSpheroidalChart(1.0), 1, 1, q)
            @test Q₁¹[1] == 0
            @test Q₁¹[2] ≈ sqrt(q^2 - 1) * (Q₀ - q / (q^2 - 1)) rtol = 1.0e-12
        end
        for τ in (0.4, 1.5, 6.0)
            O = OblateSpheroidalChart(1.0)
            @test radial_P(O, 2, 0, τ) ≈ [1, τ, (3τ^2 + 1) / 2] rtol = 1.0e-13
            Q₀ = atan(1 / τ)                               # arccot τ
            @test radial_Q(O, 1, 0, τ) ≈ [Q₀, 1 - τ * Q₀] rtol = 1.0e-12
        end
        for r in (0.4, 2.0)
            @test radial_P(SphericalChart(), 4, 1, r) ≈ [0, r, r^2, r^3, r^4]
            @test radial_Q(SphericalChart(), 3, 0, r) ≈ [1 / r, 1 / r^2, 1 / r^3, 1 / r^4]
            @test radial_Q(SphericalChart(), 3, 2, r)[1:2] == [0, 0]
        end
        # the semi-focal distance does not enter
        @test radial_P(ProlateSpheroidalChart(3.0), 3, 1, 2.0) == radial_P(ProlateSpheroidalChart(1.0), 3, 1, 2.0)
    end

    @testsection "Casoratian" begin
        # 𝖯_{n+1}𝖰_n - s 𝖯_n𝖰_{n+1} = (-1)^m (n+m)!/(n-m+1)!
        for (chart, s, ζs) in (
                    (ProlateSpheroidalChart(1.0), 1, (1.3, 2.5, 7.0)),
                    (OblateSpheroidalChart(1.0), -1, (0.4, 1.5, 4.0)),
                ), ζ in ζs
            N = 8
            for m in 0:5
                P, Q = radial_P(chart, N, m, ζ), radial_Q(chart, N, m, ζ)
                for n in m:(N - 1)
                    W = P[n + 2] * Q[n + 1] - s * P[n + 1] * Q[n + 2]
                    ref = (-1)^m * Float64(factorial(big(n + m)) / factorial(big(n - m + 1)))
                    @test W ≈ ref rtol = 1.0e-10
                end
            end
        end
    end

    @testsection "radial_Q: Neumann's integral" begin
        t, w = _sh_gauss(60)
        for q in (1.3, 2.5), n in 0:5
            Qref = sum(w[j] * _sh_legendre(n, t[j]) / (q - t[j]) for j in eachindex(t)) / 2
            @test radial_Q(ProlateSpheroidalChart(1.0), n, 0, q)[n + 1] ≈ Qref rtol = 1.0e-10
        end
        # 𝒬_n(τ) = iⁿ⁺¹ Q_n(iτ)
        for τ in (0.4, 1.5), n in 0:5
            Qiτ = sum(w[j] * _sh_legendre(n, t[j]) / (im * τ - t[j]) for j in eachindex(t)) / 2
            ref = im^(n + 1) * Qiτ
            @test abs(imag(ref)) < 1.0e-12
            @test radial_Q(OblateSpheroidalChart(1.0), n, 0, τ)[n + 1] ≈ real(ref) rtol = 1.0e-11
        end
    end

    @testsection "High degree and extended precision" begin
        for (chart, ζ, s) in (
                (ProlateSpheroidalChart(1.0), 1.5, 1), (OblateSpheroidalChart(1.0), 0.8, -1),
                (ProlateSpheroidalChart(1.0), 1.02, 1), (OblateSpheroidalChart(1.0), 0.05, -1),
            )
            Q = radial_Q(chart, 200, 0, ζ)
            @test all(isfinite, Q) && all(>(0), Q)
            @test all(diff(Q) .< 0)
            # the ratio of consecutive degrees tends to e^{-ξ}
            ξ = s == 1 ? acosh(ζ) : asinh(ζ)
            @test Q[201] / Q[200] ≈ exp(-ξ) rtol = 2.0e-2
            # the Casoratian at degree 199 (m = 0): (n+m)!/(n-m+1)! = 1/200
            P = radial_P(chart, 200, 0, ζ)
            @test P[201] * Q[200] - s * P[200] * Q[201] ≈ 1 / 200 rtol = 1.0e-8
        end
        for (c64, cbig, ζ) in (
                    (ProlateSpheroidalChart(1.0), ProlateSpheroidalChart(big"1.0"), 1.3),
                    (OblateSpheroidalChart(1.0), OblateSpheroidalChart(big"1.0"), 0.7),
                ), m in (0, 2)
            Q64 = radial_Q(c64, 12, m, ζ)
            Qbig = radial_Q(cbig, 12, m, BigFloat(ζ))
            @test eltype(Qbig) == BigFloat
            @test Float64.(Qbig) ≈ Q64 rtol = 1.0e-13
            @test Float64.(radial_P(cbig, 12, m, BigFloat(ζ))) ≈ radial_P(c64, 12, m, ζ) rtol = 1.0e-13
        end
    end

    @testsection "Each harmonic is harmonic" begin
        # the Laplacian of the harmonic times cos mϕ, by the pointwise system
        for (chart, CS, ξ) in (
                    (ProlateSpheroidalChart(1.2), coorsys_spheroidal(PointwiseGeometry, 1.2), [0.7, 0.3, 1.8]),
                    (OblateSpheroidalChart(1.2), coorsys_spheroidal_oblate(PointwiseGeometry, 1.2), [0.7, 0.3, 0.9]),
                ), (n, m) in ((0, 0), (1, 0), (2, 1), (3, 2), (4, 0), (5, 3))
            Fi = ζ -> interior_harmonics(chart, n, m, ζ)[n + 1] * cos(m * ζ[1])
            Fe = ζ -> exterior_harmonics(chart, n, m, ζ)[n + 1] * cos(m * ζ[1])
            @test laplacian(Fi, CS)(ξ) ≈ 0 atol = 1.0e-9 * (1 + abs(Fi(ξ)))
            @test laplacian(Fe, CS)(ξ) ≈ 0 atol = 1.0e-9 * (1 + abs(Fe(ξ)))
        end
        # and the harmonics are not trivially zero
        @test abs(interior_harmonics(ProlateSpheroidalChart(1.2), 3, 2, [0.7, 0.3, 1.8])[4]) > 1.0e-3
        # spherical solid harmonics
        CS = coorsys_spherical(PointwiseGeometry)
        ξ = [0.8, 0.7, 1.4]
        for (n, m) in ((2, 1), (4, 0))
            Fi = ζ -> interior_harmonics(SphericalChart(), n, m, ζ)[n + 1] * cos(m * ζ[2])
            Fe = ζ -> exterior_harmonics(SphericalChart(), n, m, ζ)[n + 1] * cos(m * ζ[2])
            @test laplacian(Fi, CS)(ξ) ≈ 0 atol = 1.0e-9
            @test laplacian(Fe, CS)(ξ) ≈ 0 atol = 1.0e-9
        end
    end

    @testsection "reexpansion_polynomials" begin
        # Φᵢ⁽ⁿ⁾(x) = (2n+4i+1)/(2n+1) (n+½)ᵢ/i! ₂F₁(-i, n+i+½; n+3/2; x), terminating
        function Φref(n, i, x::T) where {T}
            a, c = T(2n + 1) / 2, T(2n + 3) / 2
            S, term = zero(T), one(T)
            for k in 0:i
                S += term * x^k
                term *= (-i + k) * (n + i + T(1) / 2 + k) / ((c + k) * (k + 1))
            end
            pre = T(2n + 4i + 1) / (2n + 1) * prod((a + j) / (j + 1) for j in 0:(i - 1); init = one(T))
            return pre * S
        end
        for n in 0:4
            for x in (0.2, 0.5, 0.8, -0.7, 1.9)
                Φ = reexpansion_polynomials(n, 8, x)
                @test Φ[1] == 1
                @test Φ[2] ≈ (2n + 5) * (1 - x) / 2 rtol = 1.0e-13
                @test Φ ≈ [Float64(Φref(n, i, big(x))) for i in 0:8] rtol = 1.0e-10
            end
            @test all(≈(0; atol = 1.0e-12), reexpansion_polynomials(n, 8, 1.0)[2:end])
            Φb = reexpansion_polynomials(n, 8, big"0.3")
            @test maximum(abs.(Φb - [Φref(n, i, big"0.3") for i in 0:8]) ./ abs.(Φb)) < big"1.0e-60"
        end
        @test reexpansion_polynomials(2, 0, 0.5) == [1.0]
    end

    @testsection "interior_reexpansion: exact identity" begin
        pts = (
            SVector(0.7, -0.4, 0.3), SVector(0.3, 0.1, 1.5), SVector(2.5, -1.0, 1.9),
            SVector(0.5, 0.2, 0.0), SVector(0.2, 0.1, 0.0), SVector(-0.9, 0.6, -0.4),
        )
        # On the plane z = 0 the angle θ of a spherical chart is π/2 up to
        # rounding, and the harmonics that vanish there are not zero to the
        # precision of the others: those points are not used with it.
        for (to, from) in _sh_pairs, x in pts
            (iszero(x[3]) && (to isa SphericalChart || from isa SphericalChart)) && continue
            @test _sh_interior_error(to, from, x, 6) < 1.0e-11
        end
        # the same in extended precision
        bigcharts = (
            ProlateSpheroidalChart(big"1.0"), ProlateSpheroidalChart(big"1.7"),
            OblateSpheroidalChart(big"0.8"), OblateSpheroidalChart(big"2.0"), SphericalChart(),
        )
        for to in bigcharts, from in bigcharts
            (to isa SphericalChart && from isa SphericalChart) && continue
            for x in (SVector(big"0.7", big"-0.4", big"0.3"), SVector(big"0.5", big"0.2", big"0.01"), SVector(big"0.5", big"0.2", big"0.0"))
                (iszero(x[3]) && (to isa SphericalChart || from isa SphericalChart)) && continue
                @test _sh_interior_error(to, from, x, 6) < big"1.0e-60"
            end
        end
        # the same chart: the identity
        for chart in _sh_charts[1:4], n in 0:6, m in 0:n
            L = interior_reexpansion(chart, chart, n, m)
            @test length(L) == (n - m) ÷ 2 + 1
            @test L ≈ [1.0; zeros(length(L) - 1)] atol = 1.0e-13
        end
    end

    @testsection "exterior_reexpansion" begin
        imax = 80
        bounds = Dict()
        for (to, from) in _sh_pairs
            b = to isa SphericalChart ? reexpansion_bound(to, from) :
                from isa SphericalChart ? reexpansion_bound(to, from) : reexpansion_bound(to, from)
            ζs = (b + 0.5, b + 2.0)
            worst = 0.0
            for ζ in ζs, p in (0.35, -0.7), (n, m) in ((0, 0), (1, 1), (2, 0), (3, 2), (5, 1))
                ξt = _sh_point(to, p, ζ)
                x = coords_to_cartesian(to, ξt)
                ξf = cartesian_to_coords(from, x)
                K = exterior_reexpansion(to, from, n, m, imax)
                Ht = exterior_harmonics(to, n + 2imax, m, ξt)
                lhs = exterior_harmonics(from, n, m, ξf)[n + 1]
                terms = [K[i + 1] * Ht[n + 2i + 1] for i in 0:imax]
                worst = max(worst, abs(lhs - sum(terms)) / (abs(lhs) + sum(abs, terms)))
            end
            @test worst < 1.0e-10
        end
        # K₀ = (c₁/c₂)ⁿ⁺¹ between spheroidal charts
        for (to, from) in _sh_pairs
            (to isa SphericalChart || from isa SphericalChart) && continue
            for n in 0:4, m in 0:n
                @test exterior_reexpansion(to, from, n, m, 3)[1] ≈ (from.c / to.c)^(n + 1) rtol = 1.0e-13
            end
        end
        # the terms decay as e^{-2(ξ-ξ*)}: radial terms Kᵢ Q_{n+2i}(q₂), prolate → prolate
        from, to = ProlateSpheroidalChart(2.0), ProlateSpheroidalChart(1.0)
        q₂ = 3.0
        b = reexpansion_bound(to, from)
        @test b == 2.0
        K = exterior_reexpansion(to, from, 0, 0, 90)
        Q = radial_Q(to, 180, 0, q₂)
        t = [abs(K[i + 1] * Q[2i + 1]) for i in 0:90]
        @test all(isfinite, t)
        pred = exp(-2 * (acosh(q₂) - acosh(b)))
        for i in (50, 70, 85)
            @test t[i + 2] / t[i + 1] ≈ pred rtol = 5.0e-2
        end
        # table of reexpansion_bound
        P(c) = ProlateSpheroidalChart(c)
        O(c) = OblateSpheroidalChart(c)
        @test reexpansion_bound(P(1.0), P(2.0)) == 2.0
        @test reexpansion_bound(P(2.0), P(1.0)) == 1.0
        @test reexpansion_bound(O(1.0), O(2.0)) ≈ √3
        @test reexpansion_bound(O(2.0), O(1.0)) == 0.0
        @test reexpansion_bound(O(1.0), P(1.5)) == 1.5
        @test reexpansion_bound(P(1.0), O(1.5)) ≈ sqrt(1 + 1.5^2)
        @test reexpansion_bound(SphericalChart(), P(1.3)) == 1.3
        @test reexpansion_bound(SphericalChart(), O(0.7)) == 0.7
        @test reexpansion_bound(P(2.0), SphericalChart()) == 1.0
        @test reexpansion_bound(O(2.0), SphericalChart()) == 0.0
    end

    @testsection "ForwardDiff" begin
        # d/dq of Q₁ and Q₀
        for q in (1.3, 2.5)
            Q₀ = atanh(1 / q)
            @test ForwardDiff.derivative(q -> radial_Q(ProlateSpheroidalChart(1.0), 1, 0, q)[2], q) ≈
                Q₀ - q / (q^2 - 1) rtol = 1.0e-10
            @test ForwardDiff.derivative(q -> radial_Q(ProlateSpheroidalChart(1.0), 1, 0, q)[1], q) ≈
                -1 / (q^2 - 1) rtol = 1.0e-10
            # P: derivative of (3q²-1)/2 is 3q
            @test ForwardDiff.derivative(q -> radial_P(ProlateSpheroidalChart(1.0), 2, 0, q)[3], q) ≈ 3q
        end
        # oblate: d/dτ 𝒬₁ = -acot τ + τ/(τ²+1)
        for τ in (0.4, 1.5)
            @test ForwardDiff.derivative(τ -> radial_Q(OblateSpheroidalChart(1.0), 1, 0, τ)[2], τ) ≈
                -atan(1 / τ) + τ / (τ^2 + 1) rtol = 1.0e-10
        end
        # d/dc₁ of a coefficient of interior_reexpansion, (n, m) = (2, 0):
        #   L₀ = (c₂/c₁)²,  L₁ = -s₁/2 + s₂ c₂²/(2c₁²)
        for (mkfrom, s₁) in ((c -> ProlateSpheroidalChart(c), 1), (c -> OblateSpheroidalChart(c), -1)),
                (to, s₂, c₂) in (
                    (ProlateSpheroidalChart(1.7), 1, 1.7), (OblateSpheroidalChart(0.8), -1, 0.8),
                )
            c₁ = 1.3
            L(c) = interior_reexpansion(to, mkfrom(c), 2, 0)
            @test L(c₁) ≈ [(c₂ / c₁)^2, -s₁ / 2 + s₂ * c₂^2 / (2c₁^2)] rtol = 1.0e-12
            dL = ForwardDiff.derivative(L, c₁)
            @test dL ≈ [-2c₂^2 / c₁^3, -s₂ * c₂^2 / c₁^3] rtol = 1.0e-10
            # nested: second derivatives 6c₂²/c₁⁴ and 3 s₂ c₂²/c₁⁴
            d2L = ForwardDiff.derivative(c -> ForwardDiff.derivative(L, c), c₁)
            @test d2L ≈ [6c₂^2 / c₁^4, 3 * s₂ * c₂^2 / c₁^4] rtol = 1.0e-10
        end
        # a dual in the harmonic: gradient of a harmonic at a point, against the
        # derivative of the product rule of the closed forms (n = m = 0, 1)
        ξ = [0.7, 0.3, 1.8]
        g = ForwardDiff.gradient(ζ -> interior_harmonics(ProlateSpheroidalChart(1.2), 1, 0, ζ)[2], ξ)
        @test g ≈ [0, 1.8, 0.3] rtol = 1.0e-12          # P₁(p)P₁(q) = p q
    end

    @testsection "Type stability" begin
        P, O, S = ProlateSpheroidalChart(1.2), OblateSpheroidalChart(0.8), SphericalChart()
        @inferred legendre_P(4, 1, 0.3)
        for chart in (P, O, S)
            @inferred radial_P(chart, 4, 1, 1.5)
            @inferred radial_Q(chart, 4, 1, 1.5)
            @inferred interior_harmonics(chart, 4, 1, SVector(0.5, 0.3, 1.5))
            @inferred exterior_harmonics(chart, 4, 1, SVector(0.5, 0.3, 1.5))
        end
        @inferred reexpansion_polynomials(2, 5, 0.4)
        for (to, from) in ((P, O), (O, P), (S, P), (O, S))
            @inferred interior_reexpansion(to, from, 4, 1)
            @inferred exterior_reexpansion(to, from, 4, 1, 10)
            @inferred reexpansion_bound(to, from)
        end
        dq = ForwardDiff.Dual(1.5, 1.0)
        @test eltype(@inferred radial_Q(P, 4, 1, dq)) <: ForwardDiff.Dual
        @test eltype(@inferred radial_P(P, 4, 1, dq)) <: ForwardDiff.Dual
        @test eltype(@inferred legendre_P(4, 1, ForwardDiff.Dual(0.3, 1.0))) <: ForwardDiff.Dual
    end

    @testsection "Symbols" begin
        # SymPy
        xs, qs, τs = symbols("x q τ", real = true)
        c₁s, c₂s = symbols("c1 c2", positive = true)
        v = legendre_P(4, 1, xs)
        @test [Float64(v[i].subs(Dict(xs => 0.3))) for i in 1:5] ≈ legendre_P(4, 1, 0.3) rtol = 1.0e-12
        v = radial_P(ProlateSpheroidalChart(1.0), 4, 1, qs)
        @test [Float64(v[i].subs(Dict(qs => 1.7))) for i in 1:5] ≈ radial_P(ProlateSpheroidalChart(1.0), 4, 1, 1.7) rtol = 1.0e-12
        v = radial_P(OblateSpheroidalChart(1.0), 4, 1, τs)
        @test [Float64(v[i].subs(Dict(τs => 0.9))) for i in 1:5] ≈ radial_P(OblateSpheroidalChart(1.0), 4, 1, 0.9) rtol = 1.0e-12
        for (mk, s₁, s₂, to₀, from₀) in (
                (ProlateSpheroidalChart, 1, 1, ProlateSpheroidalChart(1.7), ProlateSpheroidalChart(1.1)),
                (OblateSpheroidalChart, -1, -1, OblateSpheroidalChart(1.7), OblateSpheroidalChart(1.1)),
            )
            L = interior_reexpansion(mk(c₂s), mk(c₁s), 4, 1)
            at = Dict(c₁s => 1.1, c₂s => 1.7)
            @test [Float64(l.subs(at)) for l in L] ≈ interior_reexpansion(to₀, from₀, 4, 1) rtol = 1.0e-12
        end

        # Symbolics
        Symbolics.@variables x q τ c₁ c₂
        val(e, d) = Float64(Symbolics.value(Symbolics.substitute(e, d; fold = Val(true))))
        v = legendre_P(4, 1, x)
        @test [val(v[i], Dict(x => 0.3)) for i in 1:5] ≈ legendre_P(4, 1, 0.3) rtol = 1.0e-12
        v = radial_P(ProlateSpheroidalChart(1.0), 4, 1, q)
        @test [val(v[i], Dict(q => 1.7)) for i in 1:5] ≈ radial_P(ProlateSpheroidalChart(1.0), 4, 1, 1.7) rtol = 1.0e-12
        v = radial_P(OblateSpheroidalChart(1.0), 4, 1, τ)
        @test [val(v[i], Dict(τ => 0.9)) for i in 1:5] ≈ radial_P(OblateSpheroidalChart(1.0), 4, 1, 0.9) rtol = 1.0e-12
        L = interior_reexpansion(ProlateSpheroidalChart(c₂), OblateSpheroidalChart(c₁), 4, 1)
        @test [val(l, Dict(c₁ => 1.1, c₂ => 1.7)) for l in L] ≈
            interior_reexpansion(ProlateSpheroidalChart(1.7), OblateSpheroidalChart(1.1), 4, 1) rtol = 1.0e-12
    end

    @testsection "Errors" begin
        P, O = ProlateSpheroidalChart(1.0), OblateSpheroidalChart(1.0)
        @test_throws DomainError legendre_P(3, -1, 0.3)
        @test_throws DomainError legendre_P(-1, 0, 0.3)
        @test_throws DomainError radial_P(P, 3, -1, 1.5)
        @test_throws DomainError radial_Q(O, 3, -1, 1.5)
        @test_throws DomainError radial_P(SphericalChart(), 3, -1, 1.5)
        @test_throws DomainError radial_Q(SphericalChart(), -2, 0, 1.5)
        @test_throws DomainError interior_reexpansion(P, O, 1, 2)
        @test_throws DomainError exterior_reexpansion(P, O, 1, 2, 5)
        @test_throws DomainError interior_reexpansion(SphericalChart(), O, 1, 2)
        @test_throws DomainError exterior_reexpansion(O, SphericalChart(), 1, 2, 3)
        @test_throws DomainError interior_reexpansion(P, O, 2, -1)
        # the focal set
        @test_throws DomainError radial_Q(P, 3, 0, 1.0)
        @test_throws DomainError radial_Q(O, 3, 0, 0.0)
    end
end
