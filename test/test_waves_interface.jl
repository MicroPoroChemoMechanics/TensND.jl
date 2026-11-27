# =============================================================================
#  test_waves_interface.jl — plane waves at a plane interface
#  (`interface_waves`, `reflection_transmission`, `InterfaceWave`).
#
#  References, independent of the code under test:
#  - the Christoffel problem (Γ(s) − ρ𝟏)·U = 0 with Γ(s)_ik = s_j C_jikl s_l
#    and the traction (ν·ℂ·s)·U, both as explicit index sums on the canonical
#    components, and the flux Re(τ·Ū);
#  - at normal incidence on isotropic media, the impedance formulas
#    R = (Z₂ − Z₁)/(Z₁ + Z₂) and T = 2Z₁/(Z₁ + Z₂), Z = ρV_P;
#  - at oblique incidence on isotropic media, the law of Snell-Descartes
#    sin θ/V constant, and the Zoeppritz problem (continuity of u_x, u_z, σ_zz
#    and σ_xz for P-SV) written here with displacement amplitudes and the
#    isotropic traction λ(s·d)ν + μ[s(d·ν) + (s·ν)d], its energy fractions
#    from the flux ρV²η of each wave (η its vertical slowness);
#  - beyond the critical angle, the evanescent normal slowness
#    i√(p² − 1/V_P2²);
#  - for an SH wave on isotropic media, R = (μ₁η₁ − μ₂η₂)/(μ₁η₁ + μ₂η₂) and
#    the transmitted energy 1 − R²;
#  - for a linear slip κ/(ρV_S) between identical media at normal incidence,
#    |R|² = (κ/2)²/(1 + (κ/2)²) and |T|² = 1/(1 + (κ/2)²);
#  - reciprocity: the P→P energy transmission is the same from either side;
#  - with dual numbers, the derivatives of 2Z₁/(Z₁ + Z₂) with respect to a
#    bulk modulus and of the normal slowness √(1/V_P² − s²) with respect to
#    the tangential one.
#
#  Moduli in GPa, densities in g/cm³, slownesses in s/km. Phases of the
#  polarizations are arbitrary: only moduli, displacement vectors and energy
#  fractions are compared.
# =============================================================================

using StaticArrays: SVector

const WIF_NU = [0.0, 0.0, 1.0]
const WIF_M1 = (k = 20.0, μ = 12.0, ρ = 2.5)
const WIF_M2 = (k = 45.0, μ = 25.0, ρ = 2.8)

wif_stiffness(m) = iso_stiffness(m.k, m.μ)
wif_vp(m) = sqrt((m.k + 4 * m.μ / 3) / m.ρ)
wif_vs(m) = sqrt(m.μ / m.ρ)
"Lamé's λ, μ and the density of an isotropic medium."
wif_lame(m) = (m.k - 2 * m.μ / 3, m.μ, m.ρ)

"The forward waves of the isotropic medium `m`."
wif_forward(m, ν, sₜ) = first(interface_waves(wif_stiffness(m), m.ρ, ν, sₜ))

"Reflection and transmission of `incident` from the isotropic medium `m₁` into `m₂`."
wif_rt(m₁, m₂, ν, incident; kw...) = reflection_transmission(wif_stiffness(m₁), m₁.ρ, wif_stiffness(m₂), m₂.ρ, ν, incident; kw...)

"The displacement Σ a_m U_m of the waves `idx` of a group, with amplitudes `a`."
wif_field(a, waves, idx) = sum(a[m] * waves[m].polarization for m in idx)

"The isotropic traction (ν·ℂ·s)·d = λ(s·d)ν + μ[s(d·ν) + (s·ν)d], without conjugation."
wif_iso_traction(λ, μ, ν, s, d) = λ * sum(s .* d) * ν + μ * (s * sum(d .* ν) + sum(s .* ν) * d)

"Γ(s)_ik = s_j C_jikl s_l by an explicit index sum (s may be complex)."
wif_gamma(A, s) = [sum(s[j] * A[j, i, k, l] * s[l] for j in 1:3, l in 1:3) for i in 1:3, k in 1:3]

"(ν·ℂ·s)_ik = ν_j C_jikl s_l by an explicit index sum."
wif_nus(A, ν, s) = [sum(ν[j] * A[j, i, k, l] * s[l] for j in 1:3, l in 1:3) for i in 1:3, k in 1:3]

"""
The Zoeppritz problem of a P wave of horizontal slowness `p` (along e₁) coming
from the isotropic medium `m₁` (z < 0) onto `m₂`, interface z = 0: amplitudes
(R_P, R_S, T_P, T_S) of unit displacement polarizations and their energy
fractions. Pre-critical incidence.
"""
function wif_zoeppritz(m₁, m₂, p)
    (λ₁, μ₁, ρ₁), (λ₂, μ₂, ρ₂) = wif_lame(m₁), wif_lame(m₂)
    α₁, β₁, α₂, β₂ = wif_vp(m₁), wif_vs(m₁), wif_vp(m₂), wif_vs(m₂)
    η(V) = sqrt(1 / V^2 - p^2)
    ηP₁, ηS₁, ηP₂, ηS₂ = η(α₁), η(β₁), η(α₂), η(β₂)
    ν = WIF_NU
    # (slowness, unit polarization) of each wave; SV polarizations normal to s in the xz-plane.
    inc = ([p, 0, ηP₁], α₁ * [p, 0, ηP₁])
    rP = ([p, 0, -ηP₁], α₁ * [p, 0, -ηP₁])
    rS = ([p, 0, -ηS₁], β₁ * [ηS₁, 0, p])
    tP = ([p, 0, ηP₂], α₂ * [p, 0, ηP₂])
    tS = ([p, 0, ηS₂], β₂ * [ηS₂, 0, -p])
    # rows: u_x, u_z, τ_x, τ_z
    row(w, λ, μ) = (t = wif_iso_traction(λ, μ, ν, w[1], w[2]); [w[2][1], w[2][3], t[1], t[3]])
    A = hcat(row(rP, λ₁, μ₁), row(rS, λ₁, μ₁), -row(tP, λ₂, μ₂), -row(tS, λ₂, μ₂))
    x = A \ (-row(inc, λ₁, μ₁))
    # energy flux through the interface of a unit-amplitude wave: ρV²η
    E₀ = ρ₁ * α₁^2 * ηP₁
    energy = (
        abs2(x[1]) * ρ₁ * α₁^2 * ηP₁, abs2(x[2]) * ρ₁ * β₁^2 * ηS₁,
        abs2(x[3]) * ρ₂ * α₂^2 * ηP₂, abs2(x[4]) * ρ₂ * β₂^2 * ηS₂,
    ) ./ E₀
    return (amplitude = x, energy = energy)
end

"Sine of the angle between a real slowness vector and the normal."
wif_sin(s, ν) = norm(cross(s, ν)) / norm(s)

@testsection "Waves at an interface" begin
    @testsection "the six waves solve the Christoffel problem" begin
        G₁ = from_mandel(
            let M = randn(Xoshiro(1), 6, 6)
                M * transpose(M) + 6I
            end
        )
        ν = normalize([0.2, -0.3, 0.9])
        sₜ = [0.05, 0.08, 0.0]
        sₚ = sₜ - dot(sₜ, ν) * ν                                    # the tangential part kept
        ρ = 2.0
        A = Array(canonical_components(G₁))
        forward, backward = interface_waves(G₁, ρ, ν, sₜ)
        @test forward isa NTuple{3, InterfaceWave{Float64}}
        @test backward isa NTuple{3, InterfaceWave{Float64}}
        for w in (forward..., backward...)
            s, U = w.slowness, w.polarization
            @test w.homogeneous
            @test norm((wif_gamma(A, s) - ρ * I) * U) ≤ 1.0e-13 * ρ
            @test norm(cross(ν, s - sₚ)) ≤ 1.0e-14
            @test s ≈ sₚ + w.normal_slowness * ν
            τ = wif_nus(A, ν, s) * U
            @test norm(w.traction - τ) ≤ 1.0e-13 * norm(τ)
            @test w.flux ≈ real(sum(τ .* conj.(U)))
            @test sum(abs2, U) ≈ 1
        end
        @test all(w -> w.flux > 0, forward)
        @test all(w -> w.flux < 0, backward)
        # No mirror symmetry about the interface: the two families differ.
        @test maximum(abs(f.normal_slowness + b.normal_slowness) for (f, b) in zip(forward, backward)) > 1.0e-3
        # With it (isotropic), ξ and −ξ.
        f, b = interface_waves(wif_stiffness(WIF_M1), WIF_M1.ρ, WIF_NU, [0.1, 0.05, 0.0])
        for m in 1:3
            @test isapprox(f[m].normal_slowness, -b[m].normal_slowness; atol = 1.0e-14)
        end
    end

    @testsection "normal incidence of a P wave" begin
        incident = wif_forward(WIF_M1, WIF_NU, zeros(3))[1]
        r = wif_rt(WIF_M1, WIF_M2, WIF_NU, incident)
        Z₁, Z₂ = WIF_M1.ρ * wif_vp(WIF_M1), WIF_M2.ρ * wif_vp(WIF_M2)
        @test isapprox(abs(r.R[1]), abs(Z₂ - Z₁) / (Z₁ + Z₂); rtol = 1.0e-13)
        @test isapprox(abs(r.T[1]), 2Z₁ / (Z₁ + Z₂); rtol = 1.0e-13)
        @test isapprox(abs(r.R[1]), 0.219085922925755; rtol = 1.0e-13)
        @test isapprox(abs(r.T[1]), 0.7809140770742451; rtol = 1.0e-13)
        for a in (r.R[2], r.R[3], r.T[2], r.T[3])
            @test abs(a) ≤ 1.0e-13
        end
        @test isapprox(r.energy[1][1], ((Z₂ - Z₁) / (Z₁ + Z₂))^2; rtol = 1.0e-12)
        @test isapprox(r.energy[2][1], 4Z₁ * Z₂ / (Z₁ + Z₂)^2; rtol = 1.0e-12)
        @test isapprox(sum(r.energy[1]) + sum(r.energy[2]), 1; atol = 1.0e-14)
    end

    @testsection "P wave at 25°: Snell-Descartes and Zoeppritz" begin
        θ₁ = deg2rad(25)
        p = sin(θ₁) / wif_vp(WIF_M1)
        incident = wif_forward(WIF_M1, WIF_NU, [p, 0.0, 0.0])[1]
        @test isapprox(incident.slowness, [p, 0, sqrt(1 / wif_vp(WIF_M1)^2 - p^2)]; rtol = 1.0e-13)
        r = wif_rt(WIF_M1, WIF_M2, WIF_NU, incident)
        # Each wave on its slowness sphere, with sin θ / V = p.
        for (w, V) in zip((r.reflected..., r.transmitted...), (wif_vp(WIF_M1), wif_vs(WIF_M1), wif_vs(WIF_M1), wif_vp(WIF_M2), wif_vs(WIF_M2), wif_vs(WIF_M2)))
            s = real.(w.slowness)
            @test w.homogeneous
            @test isapprox(norm(s), 1 / V; rtol = 1.0e-13)
            @test isapprox(wif_sin(s, WIF_NU) / V, sin(θ₁) / wif_vp(WIF_M1); rtol = 1.0e-13)
        end
        z = wif_zoeppritz(WIF_M1, WIF_M2, p)
        @test isapprox(sum(z.energy), 1; atol = 1.0e-14)              # the reference itself
        @test isapprox(abs(r.R[1]), abs(z.amplitude[1]); rtol = 1.0e-12)
        @test isapprox(norm(wif_field(r.R, r.reflected, 2:3)), abs(z.amplitude[2]); rtol = 1.0e-12)
        @test isapprox(abs(r.T[1]), abs(z.amplitude[3]); rtol = 1.0e-12)
        @test isapprox(norm(wif_field(r.T, r.transmitted, 2:3)), abs(z.amplitude[4]); rtol = 1.0e-12)
        @test isapprox(r.energy[1][1], z.energy[1]; rtol = 1.0e-12)
        @test isapprox(r.energy[1][2] + r.energy[1][3], z.energy[2]; rtol = 1.0e-12)
        @test isapprox(r.energy[2][1], z.energy[3]; rtol = 1.0e-12)
        @test isapprox(r.energy[2][2] + r.energy[2][3], z.energy[4]; rtol = 1.0e-12)
        @test isapprox(sum(r.energy[1]) + sum(r.energy[2]), 1; atol = 1.0e-14)
        # P-SV only: no displacement along e₂.
        @test abs(wif_field(r.R, r.reflected, 1:3)[2]) ≤ 1.0e-13
        @test abs(wif_field(r.T, r.transmitted, 1:3)[2]) ≤ 1.0e-13
    end

    @testsection "beyond the critical angle of P" begin
        θ₁ = asin(wif_vp(WIF_M1) / wif_vp(WIF_M2)) + 0.1
        p = sin(θ₁) / wif_vp(WIF_M1)
        incident = wif_forward(WIF_M1, WIF_NU, [p, 0.0, 0.0])[1]
        r = wif_rt(WIF_M1, WIF_M2, WIF_NU, incident)
        @test all(w -> w.homogeneous, r.reflected)
        evanescent = findall(w -> !w.homogeneous, r.transmitted)
        @test length(evanescent) == 1
        w = r.transmitted[only(evanescent)]
        # The transmitted P wave, decaying towards +ν.
        @test isapprox(w.normal_slowness, im * sqrt(p^2 - 1 / wif_vp(WIF_M2)^2); rtol = 1.0e-12)
        @test abs(r.energy[2][only(evanescent)]) ≤ 1.0e-14
        @test isapprox(sum(r.energy[1]) + sum(r.energy[2]), 1; atol = 1.0e-13)
    end

    @testsection "SH wave on isotropic media" begin
        θ₁ = deg2rad(30)
        p = sin(θ₁) / wif_vs(WIF_M1)
        λ₁, μ₁, _ = wif_lame(WIF_M1)
        # The S waves of an isotropic medium share one normal slowness, so the
        # SH wave is built here: polarization e₂, traction from the isotropic
        # formula.
        S = wif_forward(WIF_M1, WIF_NU, [p, 0.0, 0.0])[2]
        ξ = S.normal_slowness
        s = SVector{3, ComplexF64}(p, 0, ξ)
        U = SVector{3, ComplexF64}(0, 1, 0)
        τ = SVector{3, ComplexF64}(wif_iso_traction(λ₁, μ₁, WIF_NU, s, U))
        incident = InterfaceWave{Float64}(ξ, s, U, τ, real(sum(τ .* conj.(U))), true)
        η₁, η₂ = sqrt(1 / wif_vs(WIF_M1)^2 - p^2), sqrt(1 / wif_vs(WIF_M2)^2 - p^2)
        @test isapprox(real(ξ), η₁; rtol = 1.0e-13)
        # Z cos θ = ρV cos θ = μη
        a₁, a₂ = WIF_M1.μ * η₁, WIF_M2.μ * η₂
        R = (a₁ - a₂) / (a₁ + a₂)
        r = wif_rt(WIF_M1, WIF_M2, WIF_NU, incident)
        uR, uT = wif_field(r.R, r.reflected, 1:3), wif_field(r.T, r.transmitted, 1:3)
        @test isapprox(norm(uR), abs(R); rtol = 1.0e-12)
        @test isapprox(norm(uT), abs(2a₁ / (a₁ + a₂)); rtol = 1.0e-12)
        @test isapprox(abs(uR[2]), norm(uR); rtol = 1.0e-12)          # SH stays SH
        @test isapprox(sum(r.energy[1]), R^2; rtol = 1.0e-12)
        @test isapprox(sum(r.energy[2]), 1 - R^2; rtol = 1.0e-12)
    end

    @testsection "linear slip" begin
        incident = wif_forward(WIF_M1, WIF_NU, zeros(3))[2]          # an S wave
        welded = wif_rt(WIF_M1, WIF_M1, WIF_NU, incident)
        @test norm(welded.R) ≤ 1.0e-14
        @test isapprox(norm(wif_field(welded.T, welded.transmitted, 1:3)), 1; rtol = 1.0e-13)
        Z = WIF_M1.ρ * wif_vs(WIF_M1)
        for κ in (0.0, 0.5, 5.0)
            r = wif_rt(WIF_M1, WIF_M1, WIF_NU, incident; slip = κ / Z * I(3))
            h = (κ / 2)^2
            @test isapprox(norm(wif_field(r.R, r.reflected, 1:3))^2, h / (1 + h); atol = 1.0e-13)
            @test isapprox(norm(wif_field(r.T, r.transmitted, 1:3))^2, 1 / (1 + h); rtol = 1.0e-13)
            @test isapprox(sum(r.energy[1]), h / (1 + h); atol = 1.0e-13)
            @test isapprox(sum(r.energy[1]) + sum(r.energy[2]), 1; atol = 1.0e-14)
            if iszero(κ)
                @test isapprox(r.R, welded.R; atol = 1.0e-14)
                @test isapprox(r.T, welded.T; atol = 1.0e-14)
            end
        end
    end

    @testsection "energy balance between two anisotropic media" begin
        G₁ = from_mandel(
            let M = randn(Xoshiro(1), 6, 6)
                M * transpose(M) + 6I
            end
        )
        G₂ = from_mandel(
            let M = randn(Xoshiro(2), 6, 6)
                M * transpose(M) + 9I
            end
        )
        ν = normalize([0.2, -0.3, 0.9])
        sₜ = [0.05, 0.08, 0.0]
        for incident in first(interface_waves(G₁, 2.0, ν, sₜ))
            r = reflection_transmission(G₁, 2.0, G₂, 3.0, ν, incident)
            @test all(w -> w.homogeneous, (r.reflected..., r.transmitted...))
            @test all(≥(-1.0e-14), (r.energy[1]..., r.energy[2]...))
            @test isapprox(sum(r.energy[1]) + sum(r.energy[2]), 1; atol = 1.0e-12)
        end
    end

    @testsection "coinciding shear waves in any frame" begin
        # The two shear waves of an isotropic medium share their slowness: their
        # polarizations must be chosen so that the fluxes add, whatever the frame.
        rng = Xoshiro(20261134)
        for trial in 1:4
            Q = Matrix(rand(rng, RotMatrix{3}))
            ν, t = Q * [0.0, 0.0, 1.0], Q * [cos(trial), sin(trial), 0.0]
            forward = wif_forward(WIF_M1, ν, 0.12 * t)
            @test abs(real(sum(forward[2].polarization .* forward[3].polarization))) ≤ 1.0e-12
            @test forward[2].normal_slowness ≈ forward[3].normal_slowness
            for incident in forward
                r = wif_rt(WIF_M1, WIF_M2, ν, incident)
                @test sum(r.energy[1]) + sum(r.energy[2]) ≈ 1 rtol = 1.0e-12
            end
            # derivative of the balance with respect to a modulus of the second medium
            balance(μ) = (r = wif_rt(WIF_M1, (k = WIF_M2.k, μ = μ, ρ = WIF_M2.ρ), ν, forward[2]); sum(r.energy[1]) + sum(r.energy[2]))
            @test abs(ForwardDiff.derivative(balance, WIF_M2.μ)) ≤ 1.0e-12
        end
    end

    @testsection "reciprocity of the P-P transmission" begin
        p = sin(deg2rad(25)) / wif_vp(WIF_M1)
        sₜ = [p, 0.0, 0.0]
        E₁₂ = wif_rt(WIF_M1, WIF_M2, WIF_NU, wif_forward(WIF_M1, WIF_NU, sₜ)[1]).energy[2][1]
        E₂₁ = wif_rt(WIF_M2, WIF_M1, -WIF_NU, wif_forward(WIF_M2, -WIF_NU, sₜ)[1]).energy[2][1]
        @test isapprox(E₁₂, E₂₁; rtol = 1.0e-12)
        # The same from the Zoeppritz problem, each way.
        @test isapprox(E₁₂, wif_zoeppritz(WIF_M1, WIF_M2, p).energy[3]; rtol = 1.0e-12)
        @test isapprox(E₂₁, wif_zoeppritz(WIF_M2, WIF_M1, p).energy[3]; rtol = 1.0e-12)
    end

    @testsection "refusals" begin
        @syms k::positive
        @test_throws ArgumentError interface_waves(iso_stiffness(k, 12), 2.5, WIF_NU, [0.1, 0, 0])
        @test_throws ArgumentError interface_waves(wif_stiffness(WIF_M1), WIF_M1.ρ, zeros(3), [0.1, 0, 0])
    end

    @testsection "derivatives" begin
        # |T_P| = 2Z₁/(Z₁ + Z₂) at normal incidence, Z₂ = √(ρ₂(k₂ + 4μ₂/3)):
        # ∂Z₂/∂k₂ = ρ₂/(2Z₂).
        incident = wif_forward(WIF_M1, WIF_NU, zeros(3))[1]
        transmitted_p(k) = abs(wif_rt(WIF_M1, (k = k, μ = WIF_M2.μ, ρ = WIF_M2.ρ), WIF_NU, incident).T[1])
        Z₁, Z₂ = WIF_M1.ρ * wif_vp(WIF_M1), WIF_M2.ρ * wif_vp(WIF_M2)
        dT = -2Z₁ / (Z₁ + Z₂)^2 * WIF_M2.ρ / (2Z₂)
        @test isapprox(ForwardDiff.derivative(transmitted_p, WIF_M2.k), dT; rtol = 1.0e-10)
        # The normal slowness ξ = √(1/V_P² − s²) of the transmitted P wave: ∂ξ/∂s = −s/ξ.
        ξ(x) = real(wif_forward(WIF_M2, WIF_NU, [x, 0.0, 0.0])[1].normal_slowness)
        s = 0.1
        @test isapprox(ξ(s), sqrt(1 / wif_vp(WIF_M2)^2 - s^2); rtol = 1.0e-13)
        @test isapprox(ForwardDiff.derivative(ξ, s), -s / ξ(s); rtol = 1.0e-10)
        # The two shear waves of an isotropic medium share the normal slowness
        # ξ_S = √(ρ/μ − s²): ∂ξ_S/∂μ = −ρ/(2μ²ξ_S), at normal incidence too.
        for x in (0.0, 0.1), m in 2:3
            ξS(μ) = real(wif_forward((k = WIF_M2.k, μ = μ, ρ = WIF_M2.ρ), WIF_NU, [x, 0.0, 0.0])[m].normal_slowness)
            ξ₀ = sqrt(WIF_M2.ρ / WIF_M2.μ - x^2)
            @test isapprox(ForwardDiff.derivative(ξS, WIF_M2.μ), -WIF_M2.ρ / (2 * WIF_M2.μ^2 * ξ₀); rtol = 1.0e-10)
        end
        # Second derivative of |T_P| at normal incidence, nested duals.
        d²T = 4Z₁ / (Z₁ + Z₂)^3 * (WIF_M2.ρ / (2Z₂))^2 + 2Z₁ / (Z₁ + Z₂)^2 * WIF_M2.ρ^2 / (4Z₂^3)
        @test isapprox(ForwardDiff.derivative(k -> ForwardDiff.derivative(transmitted_p, k), WIF_M2.k), d²T; rtol = 1.0e-8)
        # Oblique incidence: the energy fractions against those of the Zoeppritz
        # system, both differentiated with respect to the shear modulus of the
        # second medium. The shear waves coincide, so their energies are summed.
        p = 0.05
        inc = wif_forward(WIF_M1, WIF_NU, [p, 0.0, 0.0])[1]
        m₂(μ) = (k = WIF_M2.k, μ = μ, ρ = WIF_M2.ρ)
        energies(μ) = (r = wif_rt(WIF_M1, m₂(μ), WIF_NU, inc); [r.energy[1][1], r.energy[1][2] + r.energy[1][3], r.energy[2][1], r.energy[2][2] + r.energy[2][3]])
        reference(μ) = collect(wif_zoeppritz(WIF_M1, m₂(μ), p).energy)
        @test isapprox(energies(WIF_M2.μ), reference(WIF_M2.μ); rtol = 1.0e-11)
        @test isapprox(ForwardDiff.derivative(energies, WIF_M2.μ), ForwardDiff.derivative(reference, WIF_M2.μ); rtol = 1.0e-8)
        @test abs(sum(ForwardDiff.derivative(energies, WIF_M2.μ))) ≤ 1.0e-12
    end
end

# The stiffnesses, the normal and the tangential slowness written in bases
# drawn at random: the waves and the coefficients do not change.
@testsection "Interface waves: bases drawn at random, inference" begin
    rng = Xoshiro(20261142)
    ℂ₁ = ti_from_constants(52.9, 15.7, 21.8, 47.8, 16.4, normalize([1.0, 2.0, 2.0]); form = :components, role = :stiffness)
    ℂ₂ = wif_stiffness(WIF_M2)
    ν, sₜ = WIF_NU, [0.12, 0.05, 0.0]
    f₀, b₀ = interface_waves(ℂ₁, 2.45, ν, sₜ)
    r₀ = reflection_transmission(ℂ₁, 2.45, ℂ₂, WIF_M2.ρ, ν, f₀[1])
    for _ in 1:3
        ℂ₁ᵦ, ℂ₂ᵦ = express(ℂ₁, random_basis(rng)), express(ℂ₂, random_basis(rng))
        νᵦ, sᵦ = express(Tens(ν), random_basis(rng)), express(Tens(sₜ), random_basis(rng))
        f, b = interface_waves(ℂ₁ᵦ, 2.45, νᵦ, sᵦ)
        @test [w.normal_slowness for w in f] ≈ [w.normal_slowness for w in f₀] rtol = 1.0e-12
        @test [w.normal_slowness for w in b] ≈ [w.normal_slowness for w in b₀] rtol = 1.0e-12
        r = reflection_transmission(ℂ₁ᵦ, 2.45, ℂ₂ᵦ, WIF_M2.ρ, νᵦ, f[1])
        # the reflected waves are distinct in the TI medium; the two transmitted
        # shear waves of the isotropic one are defined together only
        @test abs.(r.R) ≈ abs.(r₀.R) rtol = 1.0e-10
        @test collect(r.energy[1]) ≈ collect(r₀.energy[1]) atol = 1.0e-12
        @test abs(r.T[1]) ≈ abs(r₀.T[1]) rtol = 1.0e-10
        @test r.energy[2][1] ≈ r₀.energy[2][1] rtol = 1.0e-10
        @test r.energy[2][2] + r.energy[2][3] ≈ r₀.energy[2][2] + r₀.energy[2][3] rtol = 1.0e-10
        field(q) = sum(q.T[m] * q.transmitted[m].polarization for m in 2:3)
        phase(q) = q.T[1] * q.transmitted[1].polarization          # the P wave fixes the common phase
        @test field(r) / (phase(r) ⋅ conj.(phase(r₀))) ≈ field(r₀) / (phase(r₀) ⋅ conj.(phase(r₀))) rtol = 1.0e-9
    end
    @inferred interface_waves(ℂ₁, 2.45, ν, sₜ)
    @inferred reflection_transmission(ℂ₁, 2.45, ℂ₂, WIF_M2.ρ, ν, f₀[1])
    @inferred reflection_transmission(ℂ₁, 2.45, ℂ₂, WIF_M2.ρ, ν, f₀[1]; slip = 0.1 * I(3))
    D = ForwardDiff.Dual{:wif}(2.45, 1.0)
    @inferred interface_waves(ℂ₁, D, ν, sₜ)
end
