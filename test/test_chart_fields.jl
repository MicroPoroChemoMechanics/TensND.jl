# Tensor fields of constant components in the local frame of a chart.
# References independent of the code under test: the moduli of a transversely
# isotropic tensor read by contraction with the radial and the tangential
# directions built from the point; the covariance of the field under a rotation
# of the point (any rotation for the spherical chart, about 𝐞₃ for the
# cylindrical one); its invariance along a ray; the symbolic field contracted
# with the unit vectors of the exact system.

using StaticArrays: SVector

@testsection "ChartTens" begin
    hill = (n = 12.0, l = 5.0, k = 19.0, m = 11.0, p = 4.0)
    x = [0.3, -0.5, 0.8]
    contraction(C, a, b, c, d) = sum(C[i, j, k, l] * a[i] * b[j] * c[k] * d[l] for i in 1:3, j in 1:3, k in 1:3, l in 1:3)
    # (radial direction, two directions across it) at the point x.
    frames = (
        (SphericalChart(), x / norm(x), normalize(cross(x, [0.0, 0.0, 1.0])), nothing),
        (CylindricalChart(), [x[1], x[2], 0.0] / hypot(x[1], x[2]), [0.0, 0.0, 1.0], nothing),
    )
    for (chart, r, t, _) in frames
        F = radial_ti(chart; hill...)
        T = tensor_at(F, x)
        @test T isa TensTI{4, Float64, 5}
        C = canonical_components(T)
        t′ = cross(r, t)
        @test contraction(C, r, r, r, r) ≈ hill.n
        @test contraction(C, r, r, t, t) ≈ hill.l
        @test (contraction(C, t, t, t, t) + contraction(C, t, t, t′, t′)) / 2 ≈ hill.k
        @test (contraction(C, t, t, t, t) - contraction(C, t, t, t′, t′)) / 2 ≈ hill.m
        @test contraction(C, r, t, r, t) ≈ hill.p
        @test contraction(C, r, t′, r, t′) ≈ hill.p
        # The same field from its Walpole coefficients, and from a pattern of no class.
        @test radial_ti(chart, hill.n, 2 * hill.k, √2 * hill.l, 2 * hill.m, 2 * hill.p) == F
        @test canonical_components(tensor_at(ChartTens(Tens(canonical_components(F.pattern)), chart), x)) ≈ C
        # Constant along a ray from the origin (spherical) or from the axis (cylindrical).
        ray(s) = chart isa SphericalChart ? (1 + s) .* x : [(1 + s) * x[1], (1 + s) * x[2], x[3] + s]
        @test norm(ForwardDiff.derivative(s -> vec(canonical_components(tensor_at(F, ray(s)))), 0.0)) < 1.0e-12
        # Order 2: aᵣ along the radius, aₜ across it.
        K = canonical_components(tensor_at(radial_ti(chart, 2.0, 5.0), x))
        @test r' * K * r ≈ 5.0 && t' * K * t ≈ 2.0 && abs(r' * K * t) < 1.0e-14
        # An isotropic tensor is the same field in every frame.
        @test tensor_at(radial_ti(chart, TensISO{3}(30.0, 12.0)), x) ≈ TensISO{3}(30.0, 12.0)
        @test canonical_components(tensor_at(radial_ti(chart, TensISO{3}(7.0)), x)) ≈ 7I
        @test tensor_at(ChartTens(TensISO{3}(30.0, 12.0), chart), x) === TensISO{3}(30.0, 12.0)
    end

    @testset "covariance under a rotation of the point" begin
        F = radial_ti(SphericalChart(); hill...)
        Q = collect(rand(RotMatrix{3}))
        rotated(C) = [sum(Q[i, a] * Q[j, b] * Q[k, c] * Q[l, d] * C[a, b, c, d] for a in 1:3, b in 1:3, c in 1:3, d in 1:3) for i in 1:3, j in 1:3, k in 1:3, l in 1:3]
        @test canonical_components(tensor_at(F, Q * x)) ≈ rotated(canonical_components(tensor_at(F, x)))
        G = radial_ti(CylindricalChart(); hill...)
        Q = [cos(0.7) -sin(0.7) 0; sin(0.7) cos(0.7) 0; 0 0 1]
        @test canonical_components(tensor_at(G, Q * x)) ≈ rotated(canonical_components(tensor_at(G, x)))
    end

    @testset "algebra" begin
        F = radial_ti(SphericalChart(); hill...)
        G = radial_ti(SphericalChart(), TensISO{3}(30.0, 12.0))
        @test canonical_components(tensor_at(2F + G / 2 - F, x)) ≈ canonical_components(tensor_at(F, x)) + canonical_components(TensISO{3}(30.0, 12.0)) / 2
        @test (F * 2).pattern == (2F).pattern && (-F).pattern == (-1 * F).pattern
        @test canonical_components(tensor_at(inv(F), x) ⊡ tensor_at(F, x)) ≈ canonical_components(TensISO{3}(1.0, 1.0))
        @test eltype(F) == Float64 && eltype(typeof(F)) == Float64
        @test F == radial_ti(SphericalChart(); hill...) && F != G
        @test_throws ArgumentError F + radial_ti(CylindricalChart(); hill...)
        @test occursin("SphericalChart", sprint(show, F))
        @test radial_axis(SphericalChart()) == (0, 0, 1) && radial_axis(CylindricalChart()) == (1, 0, 0)
    end

    @testset "dual and symbolic coefficients" begin
        d = ForwardDiff.derivative(n -> canonical_components(tensor_at(radial_ti(SphericalChart(); hill..., n = n), x))[3, 3, 3, 3], 12.0)
        @test d ≈ (x[3] / norm(x))^4
        a, b = symbols("a b", positive = true)
        for (chart, CS, i) in ((SphericalChart(), coorsys_spherical(), 3), (CylindricalChart(), coorsys_cylindrical(), 1))
            K = tensor_at(radial_ti(chart, a, b), CS)
            e = unitvec(CS)
            @test iszero(simplify((e[i] ⋅ K) ⋅ e[i] - b))
            @test all(iszero(simplify((e[j] ⋅ K) ⋅ e[j] - a)) for j in 1:3 if j != i)
        end
        # Symbolic moduli stay exact: √2 is the symbolic one.
        n, l, k, m, p = symbols("n l k m p", positive = true)
        @test get_data(radial_ti(SphericalChart(); n, l, k, m, p).pattern) == (n, 2k, sqrt(Sym(2)) * l, 2m, 2p)
        # Symbolic coefficients at a numeric point: the frame is made of numbers.
        K = canonical_components(tensor_at(radial_ti(SphericalChart(), a, b), (0.0, 2.0, 0.0)))
        @test Float64.(N.(subs.(K, a => 1, b => 5))) ≈ [1 0 0; 0 5 0; 0 0 1]
    end
end
