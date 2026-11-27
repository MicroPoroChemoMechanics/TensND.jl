# =============================================================================
#  test_waves_thomsen.jl — Thomsen's parameters of a transversely isotropic
#  medium (`thomsen_parameters`, the form `:thomsen` of `ti_constants` and
#  `ti_from_constants`, `anellipticity`) and the velocities of weak
#  anisotropy (`weak_wave_modes`, `squared_moveout_velocities`).
#
#  References, independent of the code under test:
#  - the definitions of Thomsen (1986), Geophysics 51(10), written from the
#    five components, and his identity δ* = (1 − V_S0²/V_P0²)(2δ − ε) between
#    δ and δ* = [2(C₁₃ + C₄₄)² − (C₃₃ − C₄₄)(C₁₁ + C₃₃ − 2C₄₄)]/(2C₃₃²);
#  - an elliptical medium, (C₁₃ + C₄₄)² = (C₁₁ − C₄₄)(C₃₃ − C₄₄), that is
#    δ = ε: η = 0 and ρV_P² = C₁₁ sin²θ + C₃₃ cos²θ exactly;
#  - the first-order Taylor polynomial, in the strength of the anisotropy, of
#    the exact velocities of `wave_modes`, taken symbolically (SymPy);
#  - the moveout velocity as V(0)² [1 + V″(0)/V(0)], the second derivative of
#    the exact velocity with respect to the angle by nested ForwardDiff;
#  - the analytic derivatives of the definitions with respect to C₁₃ and C₁₁.
#
#  The media are the Mesaverde sandstone and clayshale (stiffness in GPa,
#  density in g/cm³, velocities in km/s).
# =============================================================================

using StaticArrays: SVector

const WTH_SEED = 20261128
const WTH_E3 = [0.0, 0.0, 1.0]
# (C₁₁₁₁, C₁₁₂₂, C₁₁₃₃, C₃₃₃₃, C₂₃₂₃), i.e. (C₁₁, C₁₂, C₁₃, C₃₃, C₄₄) in two-index notation.
const WTH_SANDSTONE = (52.9, 15.7, 21.8, 47.8, 16.4)
const WTH_RHO_SANDSTONE = 2.45
const WTH_CLAYSHALE = (66.7, 32.2, 39.4, 40.0, 10.9)
const WTH_RHO_CLAYSHALE = 2.59

"The TI stiffness of the five components `c`, of axis `a`."
wth_ti(c, a = WTH_E3) = ti_from_constants(c..., a; form = :components, role = :stiffness)

"Thomsen's parameters (V_P0, V_S0, ε, δ, γ) from the components, by their definitions."
function wth_thomsen(c, ρ)
    C11, C12, C13, C33, C44 = c
    C66 = (C11 - C12) / 2
    ε = (C11 - C33) / (2 * C33)
    γ = (C66 - C44) / (2 * C44)
    δ = ((C13 + C44)^2 - (C33 - C44)^2) / (2 * C33 * (C33 - C44))
    return (sqrt(C33 / ρ), sqrt(C44 / ρ), ε, δ, γ)
end

"Thomsen's δ* (1986) from the components."
function wth_delta_star(c)
    C11, _, C13, C33, C44 = c
    return (2 * (C13 + C44)^2 - (C33 - C44) * (C11 + C33 - 2 * C44)) / (2 * C33^2)
end

"The exact velocity of the mode `k` along (sin θ, 0, cos θ)."
wth_velocity(t, ρ, k, θ) = wave_modes(t, ρ, [sin(θ), 0, cos(θ)])[k][1]

"Evaluate a SymPy expression at the substitutions `s` (a vector of pairs)."
wth_num(x::Sym, s) = Float64(N(tsubs(x, s...)))
wth_num(x::Real, s) = Float64(x)

@testsection "Thomsen's parameters" begin
    @testsection "definitions and published media" begin
        for (c, ρ) in ((WTH_SANDSTONE, WTH_RHO_SANDSTONE), (WTH_CLAYSHALE, WTH_RHO_CLAYSHALE))
            t = wth_ti(c)
            p = thomsen_parameters(t, ρ)
            @test keys(p) == (:V_P0, :V_S0, :ε, :δ, :γ)
            @test all(isapprox.(values(p), wth_thomsen(c, ρ); rtol = 1.0e-13))
            # Thomsen's identity between δ and δ*.
            @test isapprox(wth_delta_star(c), (1 - p.V_S0^2 / p.V_P0^2) * (2 * p.δ - p.ε); rtol = 1.0e-13)
            # The velocities along the axis are those of the Christoffel problem.
            @test isapprox(phase_velocities(t, ρ, WTH_E3), [p.V_P0, p.V_S0, p.V_S0]; rtol = 1.0e-14)
            @test ti_constants(t; form = :thomsen, density = ρ) == Tuple(p)
        end
        p = thomsen_parameters(wth_ti(WTH_SANDSTONE), WTH_RHO_SANDSTONE)
        @test all(
            isapprox.(
                values(p), (4.417035666782945, 2.58725289661069, 0.053347280334728055, 0.15766329984276323, 0.0670731707317074);
                rtol = 1.0e-12,
            )
        )
        p = thomsen_parameters(wth_ti(WTH_CLAYSHALE), WTH_RHO_CLAYSHALE)
        @test all(
            isapprox.(
                values(p), (3.9298874594592967, 2.0514614811139418, 0.33375, 0.7230584192439862, 0.29128440366972475);
                rtol = 1.0e-12,
            )
        )
    end

    @testsection "anellipticity" begin
        _, _, ε, δ, _ = wth_thomsen(WTH_SANDSTONE, 1.0)
        @test isapprox(anellipticity(wth_ti(WTH_SANDSTONE)), (ε - δ) / (1 + 2 * δ); rtol = 1.0e-13)
        @test isapprox(anellipticity(wth_ti(WTH_SANDSTONE)), -0.07930807415662047; rtol = 1.0e-12)
        # δ = ε: an elliptical medium, ρV_P² = C₁₁ sin²θ + C₃₃ cos²θ in every direction.
        ρ = 2.0
        t = ti_from_constants(3.0, 1.5, 0.2, 0.2, 0.1, WTH_E3; form = :thomsen, density = ρ)
        @test abs(anellipticity(t)) ≤ 1.0e-14
        C11, _, _, C33, _ = ti_constants(t; form = :components)
        for θ in range(0, π / 2; length = 7)
            @test isapprox(ρ * wth_velocity(t, ρ, 1, θ)^2, C11 * sin(θ)^2 + C33 * cos(θ)^2; rtol = 1.0e-13)
        end
    end

    @testsection "the form :thomsen of ti_constants and ti_from_constants" begin
        rng = Xoshiro(WTH_SEED)
        a = normalize(randn(rng, 3))
        for (c, ρ) in ((WTH_SANDSTONE, WTH_RHO_SANDSTONE), (WTH_CLAYSHALE, WTH_RHO_CLAYSHALE)), axis in (WTH_E3, a)
            t = wth_ti(c, axis)
            p = thomsen_parameters(t, ρ)
            t′ = ti_from_constants(values(p)..., axis; form = :thomsen, density = ρ)
            A, A′ = Array(canonical_components(t)), Array(canonical_components(t′))
            @test norm(A′ - A) ≤ 1.0e-13 * norm(A)
        end
        t = wth_ti(WTH_SANDSTONE)
        p = thomsen_parameters(t, WTH_RHO_SANDSTONE)
        # Two velocities among the parameters: no density, no answer.
        @test_throws ArgumentError ti_constants(t; form = :thomsen)
        @test_throws ArgumentError ti_from_constants(values(p)..., WTH_E3; form = :thomsen)
        # An unknown form still lists the forms, this one included.
        @test_throws ArgumentError ti_constants(t; form = :thompson)
        @test_throws r":thomsen" ti_constants(t; form = :thompson)
        @test_throws r":thomsen" ti_from_constants(values(p)..., WTH_E3; form = :thompson)
    end

    @testsection "derivatives" begin
        C11, C12, C13, C33, C44 = WTH_SANDSTONE
        ρ = WTH_RHO_SANDSTONE
        dp = ForwardDiff.derivative(x -> collect(values(thomsen_parameters(wth_ti((C11, C12, x, C33, C44)), ρ))), C13)
        # ∂δ/∂C₁₃ = (C₁₃ + C₄₄)/(C₃₃(C₃₃ − C₄₄)); C₁₃ enters nothing else.
        @test isapprox(dp[4], (C13 + C44) / (C33 * (C33 - C44)); rtol = 1.0e-12)
        @test all(iszero, (dp[1], dp[2], dp[3], dp[5]))
        dp = ForwardDiff.derivative(x -> collect(values(thomsen_parameters(wth_ti((x, C12, C13, C33, C44)), ρ))), C11)
        # ∂ε/∂C₁₁ = 1/(2C₃₃), ∂γ/∂C₁₁ = 1/(4C₄₄).
        @test isapprox(dp[3], 1 / (2 * C33); rtol = 1.0e-12)
        @test isapprox(dp[5], 1 / (4 * C44); rtol = 1.0e-12)
        # V_P0 ∝ ρ^(−1/2).
        p = thomsen_parameters(wth_ti(WTH_SANDSTONE), ρ)
        @test isapprox(ForwardDiff.derivative(ϱ -> thomsen_parameters(wth_ti(WTH_SANDSTONE), ϱ).V_P0, ρ), -p.V_P0 / (2 * ρ); rtol = 1.0e-12)
        @test (@inferred thomsen_parameters(wth_ti(WTH_SANDSTONE), ρ)) isa NamedTuple
    end
end

@testsection "Weak anisotropy of a TI medium" begin
    @testsection "along the axis and under rotation" begin
        t = wth_ti(WTH_SANDSTONE)
        ρ = WTH_RHO_SANDSTONE
        p = thomsen_parameters(t, ρ)
        # Along the axis the first-order velocities are exact.
        w = weak_wave_modes(t, ρ, WTH_E3)
        @test isapprox(SVector(w.qP, w.qSV, w.SH), SVector(p.V_P0, p.V_S0, p.V_S0); rtol = 1.0e-14)
        # They depend on the angle to the axis only.
        R = Matrix(RotXYZ(0.4, -0.9, 1.3))
        n = normalize([0.3, -0.5, 0.8])
        wR = weak_wave_modes(wth_ti(WTH_SANDSTONE, R * WTH_E3), ρ, R * n)
        @test all(isapprox.(values(wR), values(weak_wave_modes(t, ρ, n)); rtol = 1.0e-13))
        @test (@inferred weak_wave_modes(t, ρ, n)) isa NamedTuple
    end

    @testsection "an error of second order" begin
        # The sandstone's parameters scaled by τ: the error to the exact velocities
        # goes as τ², a factor 100 between τ = 10⁻² and 10⁻³ (10 at first order).
        ρ = WTH_RHO_SANDSTONE
        p = thomsen_parameters(wth_ti(WTH_SANDSTONE), ρ)
        θ = 0.6
        n = [sin(θ), 0.0, cos(θ)]
        err(τ) = let t = ti_from_constants(p.V_P0, p.V_S0, τ * p.ε, τ * p.δ, τ * p.γ, WTH_E3; form = :thomsen, density = ρ)
            abs.(SVector(values(weak_wave_modes(t, ρ, n))...) .- SVector(map(first, values(wave_modes(t, ρ, n)))...))
        end
        ratio = err(1.0e-2) ./ err(1.0e-3)
        @test all(r -> 50 < r < 200, ratio)
    end

    @testsection "the first-order series of the exact velocities, with symbols" begin
        @syms vs::positive dv::positive e::real d::real g::real r::positive τ::real θ::real
        vp = vs + dv                                       # vp > vs, for SymPy to simplify
        Cs(x) = ti_from_constants(vp, vs, x * e, x * d, x * g, [0, 0, 1]; form = :thomsen, density = r)
        n = [sin(θ), 0, cos(θ)]
        exact = wave_modes(Cs(τ), r, n)
        weak = weak_wave_modes(Cs(τ), r, n)
        rng = Xoshiro(WTH_SEED + 1)
        pts = [
            [
                θ => 0.1 + 1.3 * rand(rng), vs => 1 + rand(rng), dv => 0.5 + 1.5 * rand(rng),
                e => 0.6 * rand(rng) - 0.3, d => 0.6 * rand(rng) - 0.3, g => 0.6 * rand(rng) - 0.3,
                r => 1 + 2 * rand(rng), τ => 2 * rand(rng) - 1,
            ] for _ in 1:5
        ]
        for k in (:qP, :qSV, :SH)
            V = exact[k][1]
            taylor = tsubs(V, τ => 0) + τ * tsubs(tdiff(V, τ), τ => 0)
            Δ = taylor - weak[k]
            @test all(s -> abs(wth_num(Δ, s)) ≤ 1.0e-12, pts)
            # The weak velocities are linear in the parameters.
            @test tsimplify(tdiff(weak[k], τ, τ)) == 0
        end
        V = exact[:SH][1]
        @test tsimplify(tsubs(V, τ => 0) + τ * tsubs(tdiff(V, τ), τ => 0) - weak[:SH]) == 0
    end

    @testsection "moveout velocities" begin
        t = wth_ti(WTH_SANDSTONE)
        ρ = WTH_RHO_SANDSTONE
        mo = squared_moveout_velocities(t, ρ)
        @test keys(mo) == (:qP, :qSV, :SH)
        for k in 1:3
            V(θ) = wth_velocity(t, ρ, k, θ)
            V″ = ForwardDiff.derivative(θ -> ForwardDiff.derivative(V, θ), 0.0)
            V₀ = V(0.0)
            @test isapprox(V₀^2 * (1 + V″ / V₀), mo[k]; rtol = 1.0e-12)
        end
        V_P0, V_S0, ε, δ, γ = wth_thomsen(WTH_SANDSTONE, ρ)
        @test isapprox(mo.qP, V_P0^2 * (1 + 2 * δ); rtol = 1.0e-13)
        @test isapprox(mo.qSV, V_S0^2 * (1 + 2 * (V_P0 / V_S0)^2 * (ε - δ)); rtol = 1.0e-13)
        @test isapprox(mo.SH, V_S0^2 * (1 + 2 * γ); rtol = 1.0e-13)
        @test all(isapprox.(sqrt.(values(mo)), (5.065796126362031, 1.6196987040340953, 2.7553287888551297); rtol = 1.0e-12))
        # the clayshale: the qSV square is negative, and still the curvature of the exact velocity
        t₂, ρ₂ = wth_ti(WTH_CLAYSHALE), WTH_RHO_CLAYSHALE
        m₂ = squared_moveout_velocities(t₂, ρ₂)
        V₂(θ) = wth_velocity(t₂, ρ₂, 2, θ)
        @test m₂.qSV < 0
        @test isapprox(V₂(0.0)^2 + V₂(0.0) * ForwardDiff.derivative(θ -> ForwardDiff.derivative(V₂, θ), 0.0), m₂.qSV; rtol = 1.0e-10)
    end
end
