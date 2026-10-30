# A tensor stored with its symmetries keeps them through a change of basis, and a
# fourth-order tensor whose components have the minor symmetries is inverted on
# the symmetric tensors whatever its storage. From a non-orthonormal basis the
# round-off of the new components exceeded the tolerance of the storage
# decision: a minor-symmetric stiffness came back in full storage, and the
# inverse of 𝕀 + ℙ : δℂ, a 9×9 matrix with three null directions, failed in
# 2 of 200 random bases (seen in MeanFieldHomogenization's dilute schemes).

using TensND, Test, LinearAlgebra, Random, ForwardDiff
const _SK = TensND.Tensors

@testset "symmetries kept through a change of basis" begin
    rng = Xoshiro(20261030)
    I4 = tens_Id4(Val(3), Val(Float64))
    for _ in 1:200
        # Two stiffnesses in a random, ill-conditioned basis, brought back to the
        # canonical one by the product: every intermediate keeps the symmetric storage.
        B = Basis(Matrix(random_rotation(rng)) + randn(rng, 3, 3), Val(:cov))
        P = rand(rng, Tens{4, 3}; definite = true)
        δC = change_tens(rand(rng, Tens{4, 3}; definite = true), B)
        @test get_array(change_tens(δC, CanonicalBasis{3, Float64}())) isa _SK.SymmetricTensor{4, 3}
        M = I4 + P ⊡ δC
        @test get_array(M) isa _SK.SymmetricTensor{4, 3}
        # The inverse is the inverse: M⁻¹ : M is the symmetric identity.
        @test components_canon(inv(M) ⊡ M) ≈ components_canon(I4) atol = 1.0e-10
    end
    # At order 2, a symmetric tensor stays symmetric, components and storage.
    for _ in 1:20
        B = Basis(Matrix(random_rotation(rng)) + randn(rng, 3, 3), Val(:cov))
        K = rand(rng, Tens{2, 3}; definite = true)
        KB = change_tens(K, B)
        @test get_array(KB) isa _SK.SymmetricTensor{2, 3}
        @test components_canon(KB) ≈ components_canon(K)
    end
    # A pair of different variances is not symmetrized.
    T2 = change_tens(rand(rng, Tens{2, 3}; definite = true), Basis(Matrix(random_rotation(rng)) + 0.3 * randn(rng, 3, 3), Val(:cov)), (:cont, :cov))
    @test get_var(T2) == (:cont, :cov)
end

@testset "inverse of a minor-symmetric tensor in full storage" begin
    rng = Xoshiro(3)
    S = get_array(rand(rng, Tens{4, 3}; definite = true))
    full = _SK.Tensor{4, 3}(S)                        # the same components, 81 of them
    @test_throws SingularException inv(full)          # as a 9×9 matrix: singular
    @test TensND._inv_array(full) ≈ _SK.Tensor{4, 3}(inv(S))
    # Derivatives through it: d(A⁻¹)/dx = -A⁻¹ (dA/dx) A⁻¹, here A = x S.
    d = ForwardDiff.derivative(x -> TensND._inv_array(_SK.Tensor{4, 3}(x * S))[1, 1, 1, 1], 1.0)
    @test d ≈ -inv(S)[1, 1, 1, 1]
    # A tensor without the minor symmetries keeps the 9×9 inverse.
    G = _SK.Tensor{4, 3}(randn(rng, 3, 3, 3, 3))
    @test TensND._inv_array(G) == inv(G)
    dG = ForwardDiff.derivative(x -> TensND._inv_array(x * G)[1, 2, 3, 1], 1.0)
    @test dG ≈ -inv(G)[1, 2, 3, 1]
end
