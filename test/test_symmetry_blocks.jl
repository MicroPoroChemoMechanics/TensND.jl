# The blocks of each symmetry class (`symmetry_blocks`, `from_symmetry_blocks`): products
# and inverses, block by block, against the double contraction and the inverse of the same
# tensors taken through their generic form, which knows nothing of the blocks.

import Symbolics

@testsection "Symmetry blocks: products and inverses block by block" begin
    rng = Xoshiro(20261031)
    frame = random_basis(rng)
    axis = Tuple(normalize(randn(rng, 3)))
    # Positive definite members of each class, two per class, in one frame.
    spd2() = (s = rand(rng, 2, 2) .+ 0.1; s * transpose(s) + I)
    function ortho9()
        N = rand(rng, 3, 3); N = N * transpose(N) + 2I
        return (N[1, 1], N[2, 2], N[3, 3], N[1, 2], N[1, 3], N[2, 3], 1 + rand(rng), 1 + rand(rng), 1 + rand(rng))
    end
    function ti6()
        W = spd2()
        return TensTI{4}(W[1, 1], W[2, 2], W[1, 2], W[2, 1], 1 + rand(rng), 1 + rand(rng), axis)
    end
    pairs = [
        "ISO order 4" => () -> TensISO{3}(2 + rand(rng), 1 + rand(rng)),
        "ISO order 2" => () -> TensISO{3}(1 + rand(rng)),
        "cubic" => () -> TensCubic(2 + rand(rng), 1 + rand(rng), 1 + rand(rng), frame),
        "TI order 4" => ti6,
        "TI order 4, symmetric" => () -> TensTI{4}(3 + rand(rng), 2 + rand(rng), rand(rng), 1 + rand(rng), 1 + rand(rng), axis),
        "TI order 2" => () -> TensTI{2}(1 + rand(rng), 2 + rand(rng), axis),
        "ortho order 4" => () -> TensOrtho{4, Float64}(ortho9(), frame),
        "ortho order 2" => () -> TensOrtho{2, Float64}((1 + rand(rng), 2 + rand(rng), 3 + rand(rng)), frame),
    ]
    dense(X) = canonical_components(generic_tens(X))
    matrix(X) = get_order(X) == 4 ? mandel(generic_tens(X)) : Matrix(dense(X))
    # The product of the algebra: ⊡ at order 4, ⋅ at order 2.
    prod(X, Y) = get_order(X) == 4 ? X ⊡ Y : X ⋅ Y
    for (name, draw) in pairs
        A, B = draw(), draw()
        bA, bB = symmetry_blocks(A), symmetry_blocks(B)
        # The round trip, and the product and the inverse block by block.
        @test dense(from_symmetry_blocks(A, bA)) ≈ dense(A) rtol = 1.0e-13
        @test dense(from_symmetry_blocks(A, map(*, bA, bB))) ≈ dense(prod(generic_tens(A), generic_tens(B))) rtol = 1.0e-12
        @test dense(from_symmetry_blocks(A, map(inv, bA))) ≈ dense(inv(generic_tens(A))) rtol = 1.0e-12
        # The structured product and inverse of the class have the same blocks.
        @test all(map((x, y) -> isapprox(x, y; rtol = 1.0e-12), symmetry_blocks(prod(A, B)), map(*, bA, bB)))
        @test all(map((x, y) -> isapprox(x, y; rtol = 1.0e-12), symmetry_blocks(inv(A)), map(inv, bA)))
        @inferred symmetry_blocks(A)
        # The matrix of components (Kelvin-Mandel at order 4): its blocks in the adapted
        # basis are those of the tensor, and the matrix is rebuilt from them.
        M = matrix(A)
        Q = symmetry_basis(A)
        @test transpose(Q) * Q ≈ I
        @test all(map((x, y) -> isapprox(x, y; rtol = 1.0e-12), symmetry_blocks(A, M), bA))
        @test symmetry_matrix(A, bA) ≈ M rtol = 1.0e-12
        @inferred symmetry_blocks(A, M)
        @inferred symmetry_matrix(A, bA)
    end
    # A matrix outside the class has a residual off the blocks.
    ti = TensTI{4}(4.0, 3.0, 1.0, 2.0, 5.0, axis)
    M = matrix(TensOrtho{4, Float64}(ortho9(), frame))
    @test !isapprox(symmetry_matrix(ti, symmetry_blocks(ti, M)), M; rtol = 1.0e-3)
    @test_throws DimensionMismatch symmetry_blocks(ti, zeros(3, 3))
end

@testsection "Symmetry blocks: duals and symbols" begin
    # The derivative of a block with respect to a modulus, against its closed form.
    d = ForwardDiff.derivative(k -> symmetry_blocks(TensISO{3}(3k, 4.0))[1][1], 2.0)
    @test d == 3.0
    W(x) = TensTI{4}(x, 2x, x^2, 1.0, 3.0, 4.0, (0.0, 0.0, 1.0))
    dW = ForwardDiff.derivative(x -> symmetry_blocks(W(x))[1], 1.5)
    # ℓ₁ = x, ℓ₂ = 2x, ℓ₃ = x², ℓ₄ constant: the block [ℓ₁ ℓ₃; ℓ₄ ℓ₂] has derivative [1 2x; 0 2].
    @test dW ≈ [1.0 3.0; 0.0 2.0]
    # Symbols: the blocks of the product are the products of the blocks, exactly.
    Symbolics.@variables α β γ δ
    A, B = TensISO{3}(α, β), TensISO{3}(γ, δ)
    bp = symmetry_blocks(A ⊡ B)
    @test isequal(Symbolics.simplify(bp[1][1] - α * γ), 0)
    @test isequal(Symbolics.simplify(bp[2][1] - β * δ), 0)
    @test isequal(from_symmetry_blocks(A, symmetry_blocks(A)).data, A.data)
    # Symbols through the matrix of components: the blocks of the Kelvin-Mandel
    # matrix of an isotropic tensor are its moduli.
    bm = symmetry_blocks(A, symmetry_matrix(A, symmetry_blocks(A)))
    @test isequal(Symbolics.simplify(bm[1][1] - α; expand = true), 0)
    @test isequal(Symbolics.simplify(bm[2][1] - β; expand = true), 0)
    # Duals through the matrix of components, in a rotated frame.
    f = TensOrtho{2, Float64}((1.0, 2.0, 3.0), random_basis(Xoshiro(7)))
    fx(x) = symmetry_blocks(TensOrtho{2, typeof(x)}((x, x^2, 2.0), f.frame))
    dm = ForwardDiff.derivative(x -> symmetry_blocks(f, symmetry_matrix(f, fx(x)))[2][1], 3.0)
    @test dm ≈ 6.0
end
