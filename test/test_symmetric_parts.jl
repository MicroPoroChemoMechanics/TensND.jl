# Major-symmetric and symmetric parts of TensND tensors: the class is kept where
# it can be, and the components in the canonical basis are those of the
# symmetrized canonical array, whatever basis a tensor is written in.

@testsection "Symmetric parts" begin
    rng = Xoshiro(20261023)
    R = Matrix(qr(randn(rng, 3, 3)).Q)
    det(R) < 0 && (R[:, 1] .*= -1)
    ℬr = RotatedBasis(R)
    ℬg = Basis([1.0 0.3 0.1; 0.0 1.2 0.2; 0.1 0.0 0.9])     # neither orthogonal nor normed
    msym(A) = (A + permutedims(A, (3, 4, 1, 2))) / 2
    n = (0.36, 0.48, 0.8)

    # Structured classes: the class is kept.
    for t in (
            TensISO{3}(3.0, 2.0), TensCubic(9.0, 4.0, 3.0, ℬr),
            TensOrtho(9.0, 8.0, 7.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, ℬr),
            TensTI{4, Float64, 5}((1.0, 2.0, 0.3, 4.0, 5.0), n),
        )
        @test majorsymmetric(t) === t
    end
    for N in (6, 8)
        data = N == 6 ? (1.0, 2.0, 0.3, 0.7, 4.0, 5.0) : (1.0, 2.0, 0.3, 0.7, 4.0, 5.0, 0.2, 0.1)
        t = TensTI{4, Float64, N}(data, n)
        s = majorsymmetric(t)
        @test s isa TensTI{4, Float64, 5}
        @test canonical_components(s) ≈ msym(canonical_components(t))
    end
    # A general tensor: in a rotated basis and in a general one, with variances.
    A = randn(rng, 3, 3, 3, 3)
    A = (A + permutedims(A, (2, 1, 3, 4)) + permutedims(A, (1, 2, 4, 3)) + permutedims(A, (2, 1, 4, 3))) / 4
    for t in (
            Tens(A), express(Tens(A), ℬr), express(Tens(A), ℬg),
            express(Tens(A), ℬg, (:cov, :cont, :cont, :cov)),
        )
        s = majorsymmetric(t)
        @test canonical_components(s) ≈ msym(A)
        @test get_basis(s) == get_basis(t)
    end

    # Order 2.
    @test symmetric(TensISO{3}(2.0)) === TensISO{3}(2.0)
    @test symmetric(TensTI{2}(1.0, 2.0, n)) === TensTI{2}(1.0, 2.0, n)
    P = TensTI{2}(2.0, 5.0, 0.7, n)
    @test symmetric(P) isa TensTI{2, Float64, 2}
    @test canonical_components(symmetric(P)) ≈ (canonical_components(P) + transpose(canonical_components(P))) / 2
    B = randn(rng, 3, 3)
    for t in (Tens(B), express(Tens(B), ℬr), express(Tens(B), ℬg, (:cov, :cont)))
        @test canonical_components(symmetric(t)) ≈ (B + transpose(B)) / 2
    end

    # Derivatives and symbols go through.
    f(x) = canonical_components(majorsymmetric(TensTI{4, typeof(x), 8}((x, 2.0, 0.3, 0.7, 4.0, 5.0, x^2, 0.1), n)))[1, 1, 1, 1]
    @test ForwardDiff.derivative(f, 0.4) ≈ ForwardDiff.derivative(x -> canonical_components(TensTI{4, typeof(x), 8}((x, 2.0, 0.3, 0.7, 4.0, 5.0, x^2, 0.1), n))[1, 1, 1, 1], 0.4)
    TensND.Symbolics.@variables a b c
    Ps = symmetric(TensTI{2}(a, b, c, (0, 0, 1)))
    @test isequal(get_data(Ps), (a, b))
end
