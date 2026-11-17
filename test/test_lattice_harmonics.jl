# Solid harmonics of a periodic array of sources (Ewald sums). References: the
# constant of the periodic potential of the simple cubic lattice,
# -2.837297479480620/L; the lattice sums of Kushch (2013), Table A.1, p. 465, for
# the unit cube; direct sums where they converge absolutely; the independence
# of the splitting parameter; extended precision.

# Kushch (2013), Table A.1: (t, s, Y_ts(0)/N_ts) of the simple cubic lattice of
# unit period, N_ts = √((t-s)!(t+s)!), for the degrees above 2.
const _LH_KUSHCH_A1 = (
    (4, 0, 3.10822668269944), (4, 4, 1.85752072772946), (6, 0, 0.573329289434499),
    (6, 4, -1.07260088543324), (8, 0, 3.25929309334983), (8, 4, 1.22565950318498),
    (8, 8, 1.86744362272512),
)

@testsection "Lattice solid harmonics" begin
    x = [0.2, -0.1, 0.35]
    L = (1.0, 1.3, 0.8)

    @testsection "Simple cubic values" begin
        h = lattice_solid_harmonics(zeros(3), (1.0, 1.0, 1.0), 8)
        @test real(h[0, 0]) ≈ -2.83729747948062 rtol = 1.0e-13
        # cubic symmetry leaves no harmonic of degree 2, and none of odd degree
        @test abs(h[2, 0]) < 1.0e-13 && abs(h[2, 2]) < 1.0e-13 && abs(h[3, 2]) < 1.0e-12
        for (t, s, ref) in _LH_KUSHCH_A1
            @test real(h[t, s]) / solid_harmonic_norm(t, s) ≈ ref rtol = 1.0e-12
            @test abs(imag(h[t, s])) < 1.0e-9
        end
        # the sums scale as L^-(t+1)
        h2 = lattice_solid_harmonics(zeros(3), (2.0, 2.0, 2.0), 4)
        @test h2[4, 0] ≈ h[4, 0] / 2^5 rtol = 1.0e-12
        @test h2[0, 0] ≈ h[0, 0] / 2 rtol = 1.0e-12
        # in extended precision, with a cutoff that follows it
        Lb = (big"1.0", big"1.0", big"1.0")
        hb = lattice_solid_harmonics(zeros(3), Lb, 4)
        @test eltype(hb.values) == Complex{BigFloat}
        @test Float64(real(hb[0, 0])) ≈ real(h[0, 0]) rtol = 1.0e-14
        hb₂ = lattice_solid_harmonics(zeros(3), Lb, 4; α = big"2.5")
        @test maximum(abs, hb.values .- hb₂.values) < big"1.0e-25"
    end

    @testsection "Direct sums and periodicity" begin
        h = lattice_solid_harmonics(x, L, 8)
        direct = sum(irregular_solid_harmonics(x .+ [i * L[1], j * L[2], k * L[3]], 8)[8, 3] for i in -9:9, j in -9:9, k in -9:9)
        @test h[8, 3] ≈ direct rtol = 1.0e-6
        hp = lattice_solid_harmonics(x .+ [L[1], -2 * L[2], L[3]], L, 4)
        for t in 0:4, s in 0:t
            @test hp[t, s] ≈ h[t, s] rtol = 1.0e-10 atol = 1.0e-11
        end
    end

    @testsection "A lattice point and its neighborhood" begin
        # near a lattice point: the harmonic of that point, the regular field of
        # the others by the translation formula, and on the degree 1 the uniform
        # field -(4π/3V) x of the background
        L₂ = (1.0, 1.2, 0.9)
        V = prod(L₂)
        p = [0.04, -0.03, 0.05]
        kmax = 10
        on = lattice_solid_harmonics(zeros(3), L₂, kmax + 4)
        off = lattice_solid_harmonics(p, L₂, 4)
        Y = irregular_solid_harmonics(p, 4)
        y = regular_solid_harmonics(p, kmax)
        regular(t, s) = sum((-1)^(k + l) * on[t + k, s - l] * y[k, l] for k in 0:kmax for l in (-k):k)
        for (t, s) in ((2, 0), (2, 1), (3, 2), (4, 0))
            @test off[t, s] ≈ Y[t, s] + regular(t, s) rtol = 1.0e-9
        end
        @test off[1, 0] ≈ Y[1, 0] + regular(1, 0) - 4π / (3V) * y[1, 0] rtol = 1.0e-9
        @test off[1, 1] ≈ Y[1, 1] + regular(1, 1) - 4π / (3V) * 2 * y[1, 1] rtol = 1.0e-9
    end

    @testsection "Derivatives and refusals" begin
        h = lattice_solid_harmonics(x, L, 4)
        d = ForwardDiff.derivative(u -> real.(lattice_solid_harmonics(x .+ u .* [0, 0, 1], L, 3).values), 0.0)
        k = 0
        for t in 0:3, s in 0:t
            k += 1
            # ∂₃ Y*_t^s = -Y*_{t+1}^s, but for the dipole along 𝐞₃, a second
            # derivative of S whose trace 4π/V the degree 2 does not hold
            trace = (t, s) == (1, 0) ? 4π / (3 * prod(L)) : 0.0
            @test d[k] ≈ -real(h[t + 1, s]) - trace rtol = 1.0e-9 atol = 1.0e-10
        end
        c = real(lattice_solid_harmonics(zeros(3), (1.0, 1.0, 1.0), 4)[4, 0])
        @test ForwardDiff.derivative(a -> real(lattice_solid_harmonics(zeros(3), (a, a, a), 4)[4, 0]), 1.5) ≈ -5 * c / 1.5^6 rtol = 1.0e-10
        @test @inferred(lattice_solid_harmonics(x, L, 2)) isa SolidHarmonicTable{ComplexF64}
        @test_throws DimensionMismatch lattice_solid_harmonics([0.1, 0.2], L, 2)
        @test_throws DimensionMismatch lattice_solid_harmonics(x, (1.0, 1.0), 2)
        @test_throws DomainError lattice_solid_biharmonics(x, L, -1)
    end

    @testsection "Biharmonics: splitting parameter, direct sums, derivatives" begin
        for p in (x, zeros(3))
            S₁, B₁ = lattice_solid_biharmonics(p, L, 6; α = 2.0)
            S₂, B₂ = lattice_solid_biharmonics(p, L, 6; α = 3.5)
            @test maximum(abs, B₁.values .- B₂.values) < 1.0e-10 * maximum(abs, B₁.values)
            @test maximum(abs, S₁.values .- S₂.values) < 1.0e-10 * maximum(abs, S₁.values)
        end
        # B*_n^m = Σ -r² Y_n^m/(2n-1) over the images, absolutely convergent at n = 7
        S, B = lattice_solid_biharmonics(x, L, 8)
        @test S.values == lattice_solid_harmonics(x, L, 8).values
        direct = zero(ComplexF64)
        for a in -12:12, b in -12:12, c in -12:12
            xn = x .+ (a * L[1], b * L[2], c * L[3])
            direct -= sum(abs2, xn) * irregular_solid_harmonics(xn, 7)[7, 3] / 13
        end
        @test B[7, 3] ≈ direct rtol = 1.0e-4
        # ∂₃ B*_n^m = -B*_{n+1}^m - 2(n² - m²) Y*_{n-1}^m / [(2n+1)(2n-1)] (Δρ = 2S)
        d = ForwardDiff.derivative(u -> lattice_solid_biharmonics(x .+ u .* [0, 0, 1], L, 3)[2].values, 0.0)
        S, B = lattice_solid_biharmonics(x, L, 4)
        k = 0
        for n in 0:3, m in 0:n
            k += 1
            n == 0 && continue
            @test d[k] ≈ -B[n + 1, m] - 2 * (n^2 - m^2) * S[n - 1, m] / ((2n + 1) * (2n - 1)) rtol = 1.0e-9
        end
        @test @inferred(lattice_solid_biharmonics(x, L, 2)) isa Tuple{SolidHarmonicTable{ComplexF64}, SolidHarmonicTable{ComplexF64}}
    end

    @testsection "Excluded images" begin
        out = [(0, 0, 1), (-1, 0, 0), (1, 1, -1)]
        S, B = lattice_solid_biharmonics(x, L, 6; α = 1.2, exclude = out)
        Sf, Bf = lattice_solid_biharmonics(x, L, 6)
        Sd, Bd = zero(S.values), zero(B.values)
        for n in out
            xn = x .+ n .* L
            Y = irregular_solid_harmonics(xn, 6)
            r² = sum(abs2, xn)
            k = 0
            for N in 0:6, M in 0:N
                k += 1
                Sd[k] += Y[N, M]
                Bd[k] += N == 0 ? sqrt(r²) : -r² * Y[N, M] / (2N - 1)
            end
        end
        @test S.values ≈ Sf.values .- Sd rtol = 1.0e-10
        @test B.values ≈ Bf.values .- Bd rtol = 1.0e-10
        @test lattice_solid_harmonics(x, L, 6; α = 1.2, exclude = out).values == S.values
    end
end
