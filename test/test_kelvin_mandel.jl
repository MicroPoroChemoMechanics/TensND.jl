# The Kelvin-Mandel notation in any dimension.
# References independent of the code under test:
#  - Tensors.jl's `tomandel` in dimensions 2 and 3, on the same components
#    given as plain arrays to the kernels written for dimension 4 and above;
#  - in dimensions 4 and 5, the two properties that define the notation, with
#    the contractions written as explicit sums over the indices: the vector has
#    the norm of the tensor, and the matrix of a tensor of order 4 is that of
#    the linear map it is on symmetric tensors of order 2.

using LinearAlgebra: Symmetric

km_contract(C, ε) = [sum(C[i, j, k, l] * ε[k, l] for k in axes(C, 3), l in axes(C, 4)) for i in axes(C, 1), j in axes(C, 2)]
km_symmetric(rng, d) = (a = randn(rng, d, d); (a + transpose(a)) / 2)
km_minor_symmetric_part(A) = [((A[i, j, k, l] + A[j, i, k, l]) + (A[i, j, l, k] + A[j, i, l, k])) / 4 for i in axes(A, 1), j in axes(A, 2), k in axes(A, 3), l in axes(A, 4)]
km_minor_symmetric(rng, d) = km_minor_symmetric_part(randn(rng, d, d, d, d))
km_minor_form(A) = mandel(A; symmetrize = true)        # the keyword written at the call, so that it is propagated

@testsection "Kelvin-Mandel: the index pairs of a dimension" begin
    @test TensND._mandel_pairs(Val(2)) === ((1, 1), (2, 2), (1, 2))
    @test TensND._mandel_pairs(Val(3)) === ((1, 1), (2, 2), (3, 3), (2, 3), (1, 3), (1, 2))
    @test TensND._mandel_pairs(Val(4)) === ((1, 1), (2, 2), (3, 3), (4, 4), (3, 4), (2, 4), (1, 4), (2, 3), (1, 3), (1, 2))
    for d in 1:6
        P = @inferred TensND._mandel_pairs(Val(d))
        @test length(P) == TensND._mandel_length(d) == d * (d + 1) ÷ 2
        @test allunique(P) && all(p -> 1 <= p[1] <= p[2] <= d, P)
        @test TensND._mandel_pairs(d) == collect(P)                # the run-time form, as a vector
        @test [TensND._mandel_weight(Float64, p) for p in P] == [p[1] == p[2] ? 1.0 : sqrt(2.0) for p in P]
    end
    # The tables of the 3-D kernels are these, read once.
    @test TensND._KM_COUPLES === TensND._mandel_pairs(Val(3))
    @test TensND._KM_PAIR === ((1, 1), (2, 2), (3, 3), (3, 2), (3, 1), (2, 1))
end

@testsection "Kelvin-Mandel: the kernels of plain arrays against Tensors.jl, dimensions 2 and 3" begin
    rng = Xoshiro(20261127)
    for d in (2, 3)
        ε, C = km_symmetric(rng, d), km_minor_symmetric(rng, d)
        @test TensND._mandel_array(ε) ≈ tomandel(SymmetricTensor{2, d}(ε))
        @test TensND._mandel_array(C) ≈ tomandel(SymmetricTensor{4, d}(C))
        a, A = randn(rng, d, d), randn(rng, d, d, d, d)                      # no symmetry
        @test TensND._mandel_array(a) == tomandel(Tensor{2, d}(a))
        @test TensND._mandel_array(A) == tomandel(Tensor{4, d}(A))
        # … and back.
        @test TensND._from_mandel_array(tomandel(SymmetricTensor{2, d}(ε))) ≈ ε
        @test TensND._from_mandel_array(tomandel(SymmetricTensor{4, d}(C))) ≈ C
        @test TensND._from_mandel_array(tomandel(Tensor{2, d}(a))) == a
        @test TensND._from_mandel_array(tomandel(Tensor{4, d}(A))) == A
        # The public functions keep the Tensors.jl route there.
        @test mandel(Tens(C)) == tomandel(SymmetricTensor{4, d}(C))
        @test get_array(from_mandel(mandel(Tens(C)))) isa SymmetricTensor{4, d}
    end
end

@testsection "Kelvin-Mandel: dimensions 4 and 5" begin
    rng = Xoshiro(20261128)
    for d in (4, 5)
        n = d * (d + 1) ÷ 2
        ε, C = km_symmetric(rng, d), km_minor_symmetric(rng, d)
        tε, tC = Tens(ε), Tens(C)
        @test get_dim(tC) == d && get_array(tC) isa Array{Float64, 4}
        v, M = mandel(tε), mandel(tC)
        @test size(v) == (n,) && size(M) == (n, n)
        @test norm(v) ≈ norm(ε)                                    # an isometry
        @test norm(M) ≈ norm(C)
        @test M * v ≈ mandel(Tens(km_contract(C, ε)))              # the matrix of the map
        @test mandel(C) == M && mandel(ε) == v                     # arrays as well
        @test mandel(tC, CanonicalBasis{d, Float64}()) ≈ M
        @test mandel(C; symmetrize = true) ≈ M
        @test (@inferred km_minor_form(C)) isa Matrix{Float64}
        @test (@inferred TensND._mandel_array(C)) isa Matrix{Float64}
        @test (@inferred TensND._from_mandel_array(M)) isa Array{Float64, 4}
        # Back, by the size alone.
        @test get_array(from_mandel(M)) ≈ C && get_array(from_mandel(v)) ≈ ε
        @test get_dim(from_mandel(M)) == d
        # The inverse of the tensor has the inverse matrix, and is its inverse
        # on the symmetric tensors.
        S = km_minor_symmetric(rng, d)
        P = [S[i, j, k, l] + S[k, l, i, j] for i in 1:d, j in 1:d, k in 1:d, l in 1:d]     # major-symmetric too
        K = Tens(from_mandel(mandel(Tens(P)) * mandel(Tens(P)) + I))                        # positive definite
        @test mandel(inv(K)) ≈ inv(mandel(K))
        @test km_contract(get_array(inv(K)), km_contract(get_array(K), ε)) ≈ ε
        # Without the minor symmetries: the d²×d² form of the bare components.
        a, A = randn(rng, d, d), randn(rng, d, d, d, d)
        @test size(mandel(Tens(a))) == (d^2,) && size(mandel(Tens(A))) == (d^2, d^2)
        @test mandel(Tens(A)) * mandel(Tens(a)) ≈ mandel(Tens(km_contract(A, a)))
        @test get_array(from_mandel(mandel(Tens(A)); dim = d)) == A
        @test get_array(from_mandel(mandel(Tens(a)); dim = d)) == a
        @test mandel(Tens(A); symmetrize = true) ≈ mandel(Tens(km_minor_symmetric_part(A)))
        # A component that breaks the symmetry by more than rounding changes the form.
        C′ = copy(C); C′[1, 2, 1, 1] += 1.0e-6
        @test size(mandel(Tens(C′))) == (d^2, d^2)
        C″ = copy(C); C″[1, 2, 1, 1] += 4 * eps() * maximum(abs, C)
        @test size(mandel(Tens(C″))) == (n, n)
        # The display of a tensor of order 4 goes through its matrix.
        @test occursin("KM", sprint(show, MIME"text/plain"(), tC))
        # Distances and functions of a tensor are those of its matrix.
        @test to_vector(tC) ≈ vec(M)
    end
end

@testsection "Kelvin-Mandel: sizes that two dimensions share, and sizes of none" begin
    # 36 = 8·9/2 = 6²: the symmetric form of dimension 8 unless the dimension is given.
    M = randn(Xoshiro(1), 36, 36)
    @test get_dim(from_mandel(M)) == 8
    @test get_dim(from_mandel(M; dim = 6)) == 6
    @test get_dim(from_mandel(randn(Xoshiro(2), 36))) == 8
    @test get_dim(from_mandel(randn(Xoshiro(2), 36); dim = 6)) == 6
    @test_throws DimensionMismatch from_mandel(randn(7, 7))
    @test_throws "dimension 4" from_mandel(randn(15, 15); dim = 4)
    @test_throws DimensionMismatch from_mandel(randn(10, 15))
    @test TensND._mandel_dimension(10, nothing) == (4, true)
    @test TensND._mandel_dimension(16, nothing) == (4, false)
    @test TensND._mandel_dimension(16, 4) == (4, false)
end

@testsection "Kelvin-Mandel: dual numbers and symbols in dimension 4" begin
    rng = Xoshiro(20261129)
    d = 4
    C, ε = km_minor_symmetric(rng, d), km_symmetric(rng, d)
    # d/dx of ‖mandel(x C)‖² is 2x ‖C‖², first then second order.
    sq(x) = sum(abs2, mandel(Tens(x .* C)))
    @test ForwardDiff.derivative(sq, 1.5) ≈ 2 * 1.5 * norm(C)^2
    @test ForwardDiff.derivative(x -> ForwardDiff.derivative(sq, x), 1.5) ≈ 2 * norm(C)^2
    # The inverse of x K is K⁻¹ / x.
    K = get_array(from_mandel(mandel(Tens(C)) * transpose(mandel(Tens(C))) + I))
    inv11(x) = get_array(inv(Tens(x .* K)))[1, 1, 1, 1]
    @test ForwardDiff.derivative(inv11, 2.0) ≈ -get_array(inv(Tens(K)))[1, 1, 1, 1] / 4
    @test size(mandel(Tens(ForwardDiff.Dual.(C, 1.0)))) == (10, 10)          # the form is that of the values
    # Exactly, on symbols: the √2 is the symbolic one, and the round trip is the identity.
    s = [i <= j ? symbols("s$i$j", real = true) : symbols("s$j$i", real = true) for i in 1:d, j in 1:d]
    v = mandel(Tens(s))
    @test length(v) == 10 && v[5] == sqrt(Sym(2)) * s[3, 4]
    @test all(iszero, simplify.(get_array(from_mandel(v)) .- s))
    @test iszero(simplify(sum(x -> x^2, v) - sum(x -> x^2, s)))
    a = [symbols("a$i$j", real = true) for i in 1:d, j in 1:d]                # no symmetry
    @test length(mandel(Tens(a))) == 16 && get_array(from_mandel(mandel(Tens(a)); dim = 4)) == a
end
