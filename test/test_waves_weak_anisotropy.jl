# =============================================================================
#  test_waves_weak_anisotropy.jl — plane waves to first order in the deviation
#  from isotropy (`weak_plane_waves`) and shear-wave splitting
#  (`shear_splitting`).
#
#  References, independent of the code under test:
#  - the exact waves of `plane_waves` (themselves checked against the
#    Christoffel equation in `test_waves_christoffel.jl`): the error of the
#    first-order waves is of second order in the deviation, a factor 100
#    between two deviations in the ratio 10;
#  - the derivative of the exact P velocity with respect to the deviation, at
#    the isotropic tensor, equal to that of the first-order one (ForwardDiff);
#  - an isotropic tensor, for which the first-order waves are exact;
#  - with symbols, ρV_P² = (n⊗n):ℂ:(n⊗n) written out for a TI medium of axis
#    e₃, C₁₁ sin⁴θ + C₃₃ cos⁴θ + (2C₁₃ + 4C₄₄) sin²θ cos²θ;
#  - the S velocities normal to the axis of a TI medium, √(C₆₆/ρ) (SH, along
#    e₂) and √(C₄₄/ρ), and their equality along the axis and in an isotropic
#    medium; the delay ∝ ρ^(1/2) and the splitting independent of ρ.
#
#  The transversely isotropic medium is the Mesaverde sandstone (stiffness in
#  GPa, density in g/cm³, velocities in km/s).
# =============================================================================

using StaticArrays: SVector

const WWA_SEED = 20261129
const WWA_N = normalize([0.3, -0.5, 0.8])
const WWA_E3 = [0.0, 0.0, 1.0]
const WWA_RHO = 2.5
# An isotropic stiffness and a symmetric deviation of any symmetry, as canonical arrays.
const WWA_C0 = Array(canonical_components(iso_stiffness(20.0, 12.0)))
const WWA_G = Array(
    canonical_components(
        from_mandel(
            let M = randn(Xoshiro(3), 6, 6)
                (M + transpose(M)) / 2
            end
        )
    )
)
const WWA_SANDSTONE = (52.9, 15.7, 21.8, 47.8, 16.4)        # (C₁₁, C₁₂, C₁₃, C₃₃, C₄₄), GPa
const WWA_RHO_SANDSTONE = 2.45

"The isotropic stiffness deviated by `t` times the deviation."
wwa_stiffness(t) = Tens(WWA_C0 .+ t .* WWA_G)

"The three velocities of the first-order waves, as a vector."
wwa_velocities(w) = SVector(w.qP[1], w.qS1[1], w.qS2[1])

"Errors of the first-order velocities to the exact ones."
function wwa_errors(t, n = WWA_N)
    C = wwa_stiffness(t)
    return abs.(wwa_velocities(weak_plane_waves(C, WWA_RHO, n)) .- plane_waves(C, WWA_RHO, n).velocity)
end

wwa_ti(c, a = WWA_E3) = ti_from_constants(c..., a; form = :components, role = :stiffness)

"Evaluate a SymPy expression at the substitutions `s` (a vector of pairs)."
wwa_num(x::Sym, s) = Float64(N(tsubs(x, s...)))
wwa_num(x::Real, s) = Float64(x)

@testsection "Weak anisotropy of any symmetry" begin
    @testsection "an error of second order" begin
        e₂, e₃ = wwa_errors(1.0e-2), wwa_errors(1.0e-3)
        @test all(<(1.0e-7), e₂)
        @test all(r -> isapprox(r, 100; rtol = 0.02), e₂ ./ e₃)
    end

    @testsection "polarizations" begin
        C = wwa_stiffness(1.0e-2)
        w = weak_plane_waves(C, WWA_RHO, WWA_N)
        exact = plane_waves(C, WWA_RHO, WWA_N)
        @test w.qP[2] ≈ WWA_N
        u₁, u₂ = normalize(w.qS1[2]), normalize(w.qS2[2])
        @test abs(dot(u₁, WWA_N)) ≤ 1.0e-14
        @test abs(dot(u₂, WWA_N)) ≤ 1.0e-14
        @test abs(dot(u₁, u₂)) ≤ 1.0e-11
        @test 1 - abs(dot(u₁, exact.polarization[:, 2])) ≤ 5.0e-8
    end

    @testsection "exact for an isotropic tensor" begin
        k, μ, ρ = 20.0, 12.0, WWA_RHO
        w = weak_plane_waves(iso_stiffness(k, μ), ρ, WWA_N)
        VP, VS = sqrt((k + 4 * μ / 3) / ρ), sqrt(μ / ρ)
        @test isapprox(wwa_velocities(w), [VP, VS, VS]; rtol = 1.0e-14)
    end

    @testsection "the same first derivative as the exact P wave" begin
        # At t = 0 the S waves are degenerate and the square root of the 2×2
        # problem is not differentiable: the P wave only.
        dweak = ForwardDiff.derivative(t -> weak_plane_waves(wwa_stiffness(t), WWA_RHO, WWA_N).qP[1], 0.0)
        dexact = ForwardDiff.derivative(t -> plane_waves(wwa_stiffness(t), WWA_RHO, WWA_N).velocity[1], 0.0)
        @test isapprox(dweak, dexact; rtol = 1.0e-10)
    end

    @testsection "a TI medium of moderate anisotropy" begin
        C = wwa_ti(WWA_SANDSTONE)
        ρ = WWA_RHO_SANDSTONE
        maxerr = maximum(range(0, π / 2; length = 91)) do θ
            n = [sin(θ), 0.0, cos(θ)]
            V = plane_waves(C, ρ, n).velocity
            maximum(abs.(wwa_velocities(weak_plane_waves(C, ρ, n)) .- V) ./ V)
        end
        @test 1.0e-3 < maxerr < 5.0e-3
    end

    @testsection "symbols: the P wave of a TI medium" begin
        @syms c11::positive c12::positive c13::positive c33::positive c44::positive r::positive θ::real
        C = ti_from_constants(c11, c12, c13, c33, c44, [0, 0, 1]; form = :components)
        w = weak_plane_waves(C, r, [sin(θ), 0, cos(θ)])
        ref = c11 * sin(θ)^4 + c33 * cos(θ)^4 + (2 * c13 + 4 * c44) * sin(θ)^2 * cos(θ)^2
        Δ = r * w.qP[1]^2 - ref
        rng = Xoshiro(WWA_SEED)
        for _ in 1:5
            s = [
                c11 => 50 + 5 * rand(rng), c12 => 15 + 3 * rand(rng), c13 => 20 + 3 * rand(rng),
                c33 => 45 + 5 * rand(rng), c44 => 15 + 3 * rand(rng), r => 2 + rand(rng),
                θ => 0.1 + 1.3 * rand(rng),
            ]
            @test abs(wwa_num(Δ, s)) ≤ 1.0e-12 * wwa_num(c11, s)
        end
    end
end

@testsection "Shear-wave splitting" begin
    @testsection "read on the exact waves" begin
        C = wwa_stiffness(0.3)
        w = plane_waves(C, WWA_RHO, WWA_N)
        s = shear_splitting(C, WWA_RHO, WWA_N)
        V₁, V₂ = w.velocity[2], w.velocity[3]
        @test s.splitting ≈ (V₁ - V₂) / V₁
        @test s.delay ≈ 1 / V₂ - 1 / V₁
        @test s.fast == w.polarization[:, 2]
        @test 0 < s.splitting < 1
    end

    @testsection "no splitting in an isotropic medium and along a TI axis" begin
        s = shear_splitting(iso_stiffness(5.0, 2.0), 2.0, WWA_N)
        @test abs(s.splitting) ≤ 1.0e-14
        @test abs(s.delay) ≤ 1.0e-14
        s = shear_splitting(wwa_ti(WWA_SANDSTONE), WWA_RHO_SANDSTONE, WWA_E3)
        @test abs(s.splitting) ≤ 1.0e-14
        @test abs(s.delay) ≤ 1.0e-14
    end

    @testsection "normal to the axis of a TI medium" begin
        C11, C12, _, _, C44 = WWA_SANDSTONE
        C66 = (C11 - C12) / 2                                     # 18.6 > C₄₄: SH is the fast wave
        ρ = WWA_RHO_SANDSTONE
        s = shear_splitting(wwa_ti(WWA_SANDSTONE), ρ, (1, 0, 0))
        @test isapprox(s.splitting, 1 - sqrt(C44 / C66); rtol = 1.0e-13)
        @test isapprox(s.delay, sqrt(ρ / C44) - sqrt(ρ / C66); rtol = 1.0e-13)
        @test isapprox(abs(s.fast[2]), 1; atol = 1.0e-14)
        # The splitting does not depend on the density, the delay goes as √ρ.
        S(ϱ) = shear_splitting(wwa_ti(WWA_SANDSTONE), ϱ, (1, 0, 0))
        @test abs(ForwardDiff.derivative(ϱ -> S(ϱ).splitting, ρ)) ≤ 1.0e-14
        @test isapprox(ForwardDiff.derivative(ϱ -> S(ϱ).delay, ρ), s.delay / (2 * ρ); rtol = 1.0e-12)
    end
end

# The stiffness and the direction written in two bases drawn at random, each
# of its own.
@testsection "Weak anisotropy and splitting: bases drawn at random, inference" begin
    rng = Xoshiro(WWA_SEED + 11)
    ℂ = wwa_stiffness(0.05)
    v₀ = wwa_velocities(weak_plane_waves(ℂ, WWA_RHO, WWA_N))
    s₀ = shear_splitting(ℂ, WWA_RHO, WWA_N)
    for _ in 1:3
        ℂᵦ = express(ℂ, random_basis(rng))
        nᵦ = express(Tens(WWA_N), random_basis(rng))
        @test wwa_velocities(weak_plane_waves(ℂᵦ, WWA_RHO, nᵦ)) ≈ v₀ rtol = 1.0e-13
        s = shear_splitting(ℂᵦ, WWA_RHO, nᵦ)
        @test s.splitting ≈ s₀.splitting rtol = 1.0e-10
        @test s.delay ≈ s₀.delay rtol = 1.0e-10
        @test min(norm(s.fast - s₀.fast), norm(s.fast + s₀.fast)) ≤ 1.0e-10
    end
    @inferred weak_plane_waves(ℂ, WWA_RHO, WWA_N)
    @inferred shear_splitting(ℂ, WWA_RHO, WWA_N)
end
