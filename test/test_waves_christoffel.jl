# =============================================================================
#  test_waves_christoffel.jl — the plane waves of a homogeneous elastic medium:
#  `plane_waves`, `phase_velocities`, `slowness`, `christoffel_invariants` and
#  the closed forms of `wave_modes`.
#
#  References, independent of the code under test:
#  - the Christoffel equation (Γ/ρ − V² 𝟏)·U = 0, the acoustic tensor
#    Γ = n·ℂ·n written as an explicit index sum (OMEinsum);
#  - the group velocity as the gradient, with respect to the wave vector k, of
#    the frequency k ↦ |k| V(k/|k|), homogeneous of degree one (ForwardDiff);
#  - the closed forms along and normal to the axis of a transversely isotropic
#    medium (√(C₃₃/ρ), √(C₄₄/ρ), √(C₁₁/ρ), √(C₆₆/ρ)) and those of an isotropic
#    one (√((k + 4μ/3)/ρ), √(μ/ρ), group velocity V n);
#  - the invariants of Γ/ρ as the symmetric functions of the squared
#    velocities, and, with symbols, as those of the acoustic tensor of a TI
#    medium of axis e₃ written out;
#  - the closed forms of `wave_modes` checked as eigenpairs of Γ/ρ;
#  - invariance under a change of the basis the tensor is written in, and
#    covariance under a rotation of the medium and of the direction, the
#    components rotated by explicit index sums;
#  - the derivative of a velocity with respect to a stiffness component by the
#    Hellmann–Feynman formula U·(∂Γ/∂c)·U/(2ρV), and with respect to the
#    density from V ∝ ρ^(−1/2).
#
#  The transversely isotropic media are the Mesaverde sandstone and clayshale
#  (stiffness in GPa, density in g/cm³, velocities in km/s).
# =============================================================================

using StaticArrays: SVector

const WCH_SEED = 20261127
const WCH_N = normalize([0.3, -0.5, 0.8])
const WCH_E3 = [0.0, 0.0, 1.0]

"Canonical components of a tensor, as a plain array."
wch_arr(t) = Array(canonical_components(t))

"Γ(n) = n·ℂ·n of a canonical component array, by an explicit index sum."
wch_gamma(A, n) = ein"k,kijl,l->ij"(collect(n), A, collect(n))

"The smaller of |u − v| and |u + v|: the distance between two unit vectors up to their sign."
wch_signless(u, v) = min(norm(u - v), norm(u + v))

"Mesaverde sandstone, (C₁₁₁₁, C₁₁₂₂, C₁₁₃₃, C₃₃₃₃, C₂₃₂₃) in GPa, TI of axis `a`."
wch_sandstone(a = WCH_E3) = ti_from_constants(52.9, 15.7, 21.8, 47.8, 16.4, a; form = :components, role = :stiffness)
const WCH_RHO_SANDSTONE = 2.45
"Mesaverde clayshale, same layout."
wch_clayshale(a = WCH_E3) = ti_from_constants(66.7, 32.2, 39.4, 40.0, 10.9, a; form = :components, role = :stiffness)
const WCH_RHO_CLAYSHALE = 2.59

"A positive definite stiffness without any symmetry class."
wch_general() = from_mandel(
    let M = randn(Xoshiro(1), 6, 6)
        M * transpose(M) + 3I
    end
)

"The three media of the generic checks, with their densities."
wch_media() = (
    ("isotropic", iso_stiffness(5.0, 2.0), 2.0),
    ("sandstone", wch_sandstone(), WCH_RHO_SANDSTONE),
    ("general", wch_general(), 2.0),
)

"Evaluate a SymPy expression at the substitutions `s` (a vector of pairs)."
wch_num(x::Sym, s) = Float64(N(tsubs(x, s...)))
wch_num(x::Real, s) = Float64(x)

@testsection "Plane waves: the Christoffel problem" begin
    @testsection "eigenpairs of the acoustic tensor" begin
        for (name, C, ρ) in wch_media()
            w = plane_waves(C, ρ, WCH_N)
            @test w isa PlaneWaves{Float64}
            Γ = wch_gamma(wch_arr(C), WCH_N) / ρ
            scale = w.velocity[1]^2
            for m in 1:3
                U = w.polarization[:, m]
                @test norm(Γ * U - w.velocity[m]^2 * U) ≤ 1.0e-13 * scale
                # The projection of the group velocity on the direction is the phase velocity.
                @test abs(dot(w.group[:, m], WCH_N) - w.velocity[m]) ≤ 1.0e-14 * w.velocity[1]
            end
            @test isapprox(transpose(w.polarization) * w.polarization, Matrix(1.0I, 3, 3); atol = 1.0e-14)
            @test issorted(w.velocity; rev = true)
            @test w.direction ≈ WCH_N
            @test phase_velocities(C, ρ, WCH_N) == w.velocity
            # Any nonzero vector is a direction: only its line counts.
            @test plane_waves(C, ρ, 3.7 * WCH_N).velocity ≈ w.velocity
        end
    end

    @testsection "group velocity as the gradient of the frequency" begin
        n = collect(WCH_N)
        for (name, C, ρ) in wch_media()
            w = plane_waves(C, ρ, WCH_N)
            ω(k, m) = norm(k) * plane_waves(C, ρ, k).velocity[m]
            tol = 1.0e-13 * w.velocity[1]
            if name == "isotropic"
                # The two S waves are degenerate: only their sum is differentiable.
                @test norm(ForwardDiff.gradient(k -> ω(k, 1), n) - w.group[:, 1]) ≤ tol
                ∇S = ForwardDiff.gradient(k -> ω(k, 2) + ω(k, 3), n)
                @test norm(∇S - (w.group[:, 2] + w.group[:, 3])) ≤ tol
            else
                for m in 1:3
                    @test norm(ForwardDiff.gradient(k -> ω(k, m), n) - w.group[:, m]) ≤ tol
                end
            end
        end
    end

    @testsection "closed forms: isotropic medium" begin
        ρ = 2.0
        C = iso_stiffness(5.0, 2.0)
        VP, VS = sqrt((5.0 + 4 * 2.0 / 3) / ρ), sqrt(2.0 / ρ)
        w = plane_waves(C, ρ, WCH_N)
        @test isapprox(w.velocity, [VP, VS, VS]; rtol = 1.0e-14)
        for m in 1:3
            # Isotropic: the energy goes along the direction, at the phase velocity.
            @test norm(w.group[:, m] - w.velocity[m] * WCH_N) ≤ 1.0e-14 * VP
        end
        @test wch_signless(w.polarization[:, 1], WCH_N) ≤ 1.0e-14
        modes = wave_modes(C, ρ, WCH_N)
        @test modes.P[1] ≈ VP
        @test modes.S[1] ≈ VS
        @test modes.P[2] == Tuple(WCH_N)
        @test modes.S[2] === nothing
    end

    @testsection "closed forms: TI medium along and normal to its axis" begin
        for (C, c, ρ) in (
                (wch_sandstone(), (52.9, 15.7, 21.8, 47.8, 16.4), WCH_RHO_SANDSTONE),
                (wch_clayshale(), (66.7, 32.2, 39.4, 40.0, 10.9), WCH_RHO_CLAYSHALE),
            )
            C11, C12, C13, C33, C44 = c
            C66 = (C11 - C12) / 2
            # Along the axis: P at √(C₃₃/ρ), the two S at √(C₄₄/ρ), the energy along the axis.
            w = plane_waves(C, ρ, (0, 0, 1))
            @test isapprox(w.velocity, [sqrt(C33 / ρ), sqrt(C44 / ρ), sqrt(C44 / ρ)]; rtol = 1.0e-14)
            for m in 1:3
                @test norm(w.group[:, m] - w.velocity[m] * WCH_E3) ≤ 1.0e-14 * w.velocity[1]
            end
            # Normal to it, along e₁: P at √(C₁₁/ρ), SH (along e₂) at √(C₆₆/ρ), SV (along e₃) at √(C₄₄/ρ).
            @test C66 > C44           # both media: SH faster than SV normal to the axis
            w = plane_waves(C, ρ, (1, 0, 0))
            @test isapprox(w.velocity, [sqrt(C11 / ρ), sqrt(C66 / ρ), sqrt(C44 / ρ)]; rtol = 1.0e-14)
            @test wch_signless(w.polarization[:, 1], [1.0, 0.0, 0.0]) ≤ 1.0e-14
            @test wch_signless(w.polarization[:, 2], [0.0, 1.0, 0.0]) ≤ 1.0e-14
            @test wch_signless(w.polarization[:, 3], [0.0, 0.0, 1.0]) ≤ 1.0e-14
        end
        # The values measured once along the axis of the sandstone.
        w = plane_waves(wch_sandstone(), WCH_RHO_SANDSTONE, WCH_E3)
        @test isapprox(w.velocity, [4.417035666782945, 2.58725289661069, 2.58725289661069]; rtol = 1.0e-14)
    end

    @testsection "slowness" begin
        w = plane_waves(wch_general(), 2.0, WCH_N)
        S = slowness(w)
        for m in 1:3
            @test S[:, m] ≈ WCH_N / w.velocity[m]
            @test dot(S[:, m], w.group[:, m]) ≈ 1     # s·g = 1: the group velocity projects on n at V
        end
    end

    @testsection "Christoffel invariants" begin
        for (name, C, ρ) in wch_media()
            v = plane_waves(C, ρ, WCH_N).velocity .^ 2
            I₁, I₂, I₃ = christoffel_invariants(C, ρ, WCH_N)
            @test isapprox(I₁, sum(v); rtol = 1.0e-11)
            @test isapprox(I₂, v[1] * v[2] + v[1] * v[3] + v[2] * v[3]; rtol = 1.0e-11)
            @test isapprox(I₃, prod(v); rtol = 1.0e-11)
        end
    end

    @testsection "wave_modes of a TI medium" begin
        ρ = WCH_RHO_SANDSTONE
        rng = Xoshiro(WCH_SEED)
        # Axis e₃ with n in (e₁, e₃), then a random axis and a random direction.
        a₀ = normalize(randn(rng, 3))
        n₀ = normalize(randn(rng, 3))
        n₀ = dot(n₀, a₀) ≥ 0 ? n₀ : -n₀       # the coupling term of `wave_modes` assumes cos θ ≥ 0
        θ = 0.7
        for (a, n) in ((WCH_E3, [sin(θ), 0.0, cos(θ)]), (a₀, n₀))
            C = wch_sandstone(a)
            modes = wave_modes(C, ρ, n)
            w = plane_waves(C, ρ, n)
            @test isapprox(sort([modes.qP[1], modes.qSV[1], modes.SH[1]]; rev = true), w.velocity; rtol = 1.0e-14)
            @test modes.qP[1] ≈ w.velocity[1]
            Γ = wch_gamma(wch_arr(C), n) / ρ
            tol = 1.0e-13 * w.velocity[1]^2
            # qP and qSV in the plane of the axis and the direction, components (transverse, axial).
            eₜ = normalize(n - dot(n, a) * a)
            for (V, p) in (modes.qP, modes.qSV)
                u = p[1] * eₜ + p[2] * a
                @test norm(Γ * u - V^2 * u) ≤ tol * norm(u)
            end
            # SH along a × n.
            V, u = modes.SH
            @test u ≈ cross(a, n)
            @test norm(Γ * u - V^2 * u) ≤ tol * norm(u)
            # The direction needs no normalization.
            @test wave_modes(C, ρ, 2.5 * n).SH[1] ≈ V
            @test wave_modes(C, ρ, 2.5 * n).qSV[1] ≈ modes.qSV[1]
        end
    end

    @testsection "symbols: closed forms and invariants" begin
        @syms k::positive μ::positive r::positive
        m = wave_modes(iso_stiffness(k, μ), r, (0, 0, 1))
        @test tsimplify(m.P[1]^2 - (k + 4 * μ / 3) / r) == 0
        @test tsimplify(m.S[1]^2 - μ / r) == 0

        @syms c11::positive c12::positive c13::positive c33::positive c44::positive θ::real
        C = ti_from_constants(c11, c12, c13, c33, c44, [0, 0, 1]; form = :components)
        n = [sin(θ), 0, cos(θ)]
        modes = wave_modes(C, r, n)
        I₁, I₂, I₃ = christoffel_invariants(C, r, n)
        # The acoustic tensor of a TI medium of axis e₃ along (sin θ, 0, cos θ), written out.
        s, c = sin(θ), cos(θ)
        Γ = [
            c11 * s^2 + c44 * c^2   0                                   (c13 + c44) * s * c
            0                       (c11 - c12) / 2 * s^2 + c44 * c^2   0
            (c13 + c44) * s * c     0                                   c44 * s^2 + c33 * c^2
        ]
        rng = Xoshiro(WCH_SEED + 1)
        for _ in 1:3
            # A stable medium: c₁₁ > |c₁₂|, (c₁₁ + c₁₂) c₃₃ > 2c₁₃², c₄₄ > 0.
            pts = [
                c11 => 50 + 5 * rand(rng), c12 => 15 + 3 * rand(rng), c13 => 20 + 3 * rand(rng),
                c33 => 45 + 5 * rand(rng), c44 => 15 + 3 * rand(rng), r => 2 + rand(rng),
                θ => 0.1 + 1.3 * rand(rng),
            ]
            G = [wch_num(x, pts) for x in Γ] / wch_num(r, pts)
            λ = sort(eigvals(Symmetric(G)); rev = true)
            v = sort([wch_num(modes.qP[1]^2, pts), wch_num(modes.qSV[1]^2, pts), wch_num(modes.SH[1]^2, pts)]; rev = true)
            @test isapprox(v, λ; rtol = 1.0e-12)
            @test isapprox(wch_num(I₁, pts), tr(G); rtol = 1.0e-12)
            @test isapprox(wch_num(I₂, pts), (tr(G)^2 - tr(G * G)) / 2; rtol = 1.0e-12)
            @test isapprox(wch_num(I₃, pts), det(G); rtol = 1.0e-12)
        end
    end

    @testsection "independent of the basis, covariant under rotations" begin
        ρ = 2.0
        C = wch_general()
        w = plane_waves(C, ρ, WCH_N)
        # The same tensor written in another basis.
        wℬ = plane_waves(express(C, RotatedBasis(0.3, 0.7, 1.1)), ρ, WCH_N)
        @test isapprox(wℬ.velocity, w.velocity; rtol = 1.0e-12)
        @test isapprox(wℬ.group, w.group; rtol = 1.0e-12)
        @test all(m -> wch_signless(wℬ.polarization[:, m], w.polarization[:, m]) ≤ 1.0e-12, 1:3)

        # The medium and the direction turned by R: velocities kept, vectors turned.
        R = Matrix(RotXYZ(0.4, -0.9, 1.3))
        CR = Tens(ein"ia,jb,kc,ld,abcd->ijkl"(R, R, R, R, wch_arr(C)))
        for (C₀, CR₀, ρ₀) in ((C, CR, ρ), (wch_sandstone(), wch_sandstone(R * WCH_E3), WCH_RHO_SANDSTONE))
            w₀ = plane_waves(C₀, ρ₀, WCH_N)
            wR = plane_waves(CR₀, ρ₀, R * WCH_N)
            @test isapprox(wR.velocity, w₀.velocity; rtol = 1.0e-12)
            @test isapprox(wR.group, R * w₀.group; atol = 1.0e-12 * w₀.velocity[1])
            @test all(m -> wch_signless(wR.polarization[:, m], R * w₀.polarization[:, m]) ≤ 1.0e-10, 1:3)
        end
    end

    @testsection "derivatives" begin
        ρ = 2.0
        A = wch_arr(wch_general())
        G = wch_arr(
            from_mandel(
                let M = randn(Xoshiro(3), 6, 6)
                    (M + transpose(M)) / 2
                end
            )
        )
        w = plane_waves(Tens(A), ρ, WCH_N)
        Γ′ = wch_gamma(G, WCH_N)
        for m in 1:3
            # Hellmann–Feynman: ∂(ρV²)/∂c = U·(∂Γ/∂c)·U for a simple eigenvalue.
            dV = ForwardDiff.derivative(c -> phase_velocities(Tens(A .+ c .* G), ρ, WCH_N)[m], 0.0)
            U = w.polarization[:, m]
            @test isapprox(dV, dot(U, Γ′ * U) / (2 * ρ * w.velocity[m]); rtol = 1.0e-10)
            # V ∝ ρ^(−1/2): V′ = −V/(2ρ), V″ = 3V/(4ρ²).
            V(ϱ) = phase_velocities(Tens(A), ϱ, WCH_N)[m]
            @test isapprox(ForwardDiff.derivative(V, ρ), -w.velocity[m] / (2 * ρ); rtol = 1.0e-12)
            V″ = ForwardDiff.derivative(ϱ -> ForwardDiff.derivative(V, ϱ), ρ)
            @test isapprox(V″, 3 * w.velocity[m] / (4 * ρ^2); rtol = 1.0e-12)
        end
        @test plane_waves(Tens(A), ForwardDiff.Dual(ρ, 1.0), WCH_N) isa PlaneWaves{<:ForwardDiff.Dual}
    end

    @testsection "types and refusals" begin
        C = wch_sandstone()
        @test (@inferred plane_waves(C, WCH_RHO_SANDSTONE, SVector(0.3, 0.1, 0.9))) isa PlaneWaves{Float64}
        @test (@inferred plane_waves(wch_general(), 2.0, SVector(0.3, 0.1, 0.9))) isa PlaneWaves{Float64}
        @test (@inferred phase_velocities(C, WCH_RHO_SANDSTONE, SVector(0.3, 0.1, 0.9))) isa SVector{3, Float64}
        @test occursin("PlaneWaves along", sprint(show, plane_waves(C, WCH_RHO_SANDSTONE, WCH_E3)))

        @test_throws ArgumentError plane_waves(C, WCH_RHO_SANDSTONE, (0, 0, 0))
        @test_throws DimensionMismatch plane_waves(C, WCH_RHO_SANDSTONE, [1.0, 0.0])
        @test_throws DimensionMismatch plane_waves(C, WCH_RHO_SANDSTONE, [1.0, 0.0, 0.0, 0.0])

        # C₁₂₁₂ = −1/2: along e₁ the transverse wave polarized along e₂ has ρV² = −1/2;
        # along e₃ the acoustic tensor is diag(1/2, 1/2, 1), positive.
        Cbad = from_mandel(Matrix(Diagonal([1.0, 1, 1, 1, 1, -1])))
        @test_throws DomainError plane_waves(Cbad, 1.0, (1, 0, 0))
        @test isapprox(plane_waves(Cbad, 1.0, (0, 0, 1)).velocity, [1, sqrt(0.5), sqrt(0.5)]; rtol = 1.0e-14)
        @test_throws DomainError plane_waves(iso_stiffness(5.0, 2.0), -1.0, (0, 0, 1))

        # Sorting the waves compares velocities, which symbols do not take.
        @syms k::positive μ::positive
        @test_throws ArgumentError plane_waves(iso_stiffness(k, μ), 1, (0, 0, 1))
    end
end
