# Tensors.jl bridge and Voigt forms: components always canonical, whatever the
# basis a tensor is written in.

@testsection "Tensors.jl bridge and Voigt forms" begin
    rng = Xoshiro(20261022)
    R = Matrix(qr(randn(rng, 3, 3)).Q)
    det(R) < 0 && (R[:, 1] .*= -1)
    ℬ = RotatedBasis(R)

    # A tensor written in a rotated basis: its own components differ from the
    # canonical ones, and the bridge must give the latter.
    S = randn(rng, 3, 3); S = S + transpose(S)
    σ = change_tens(Tens(S), ℬ)
    @test !(get_array(σ) ≈ S)
    @test Array(to_tensors(σ)) ≈ S
    @test to_tensors(σ) isa SymmetricTensor{2, 3}
    C = change_tens(Tens(get_array(TensCubic(3.0, 1.0, 0.5, CanonicalBasis{3, Float64}()))), ℬ)
    @test Array(to_tensors(C)) ≈ components_canon(C)
    @test to_tensors(TensISO{3}(3.0, 2.0)) ≈ SymmetricTensor{4, 3}((i, j, k, l) -> get_array(TensISO{3}(3.0, 2.0))[i, j, k, l])
    @test to_tensors(Tens([1.0 2.0; 2.0 4.0])) isa SymmetricTensor{2, 2}
    # A non-symmetric second-order tensor is symmetrized.
    @test Array(to_tensors(Tens([0.0 2.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]))) ≈ [0 1 0; 1 0 0; 0 0 0]

    T2 = SymmetricTensor{2, 3}((i, j) -> 1.0 * i + 2.0 * j)
    @test to_tensors(from_tensors(T2)) ≈ T2
    @test components_canon(from_tensors(T2)) ≈ Array(T2)

    # Voigt: canonical components, engineering shears only on the strain.
    v = voigt_strain(σ)
    w = voigt_stress(σ)
    @test collect(w) ≈ [S[1, 1], S[2, 2], S[3, 3], S[1, 2], S[1, 3], S[2, 3]]
    @test collect(v) ≈ [S[1, 1], S[2, 2], S[3, 3], 2S[1, 2], 2S[1, 3], 2S[2, 3]]
    @test components_canon(strain_from_voigt(v)) ≈ S
    @test components_canon(stress_from_voigt(w)) ≈ S
    @test components_canon(stress_from_voigt((1, 2, 3, 4, 5, 6))) == [1 4 5; 4 2 6; 5 6 3]
    @test_throws DimensionMismatch stress_from_voigt(ones(5))

    # Derivatives go through, under nested duals too.
    f(x) = voigt_strain(change_tens(Tens(x * S), ℬ))[4]
    @test ForwardDiff.derivative(f, 0.7) ≈ 2S[1, 2]
    @test ForwardDiff.derivative(y -> ForwardDiff.derivative(x -> x * f(x), y), 0.7) ≈ 4S[1, 2]
    @test ForwardDiff.derivative(x -> components_canon(strain_from_voigt(ntuple(i -> x * i, 6)))[1, 2], 1.0) ≈ 2

    # Symbolic components stay symbolic.
    TensND.Symbolics.@variables a b
    sym = stress_from_voigt((a, b, a, b, a, b))
    @test isequal(TensND.Symbolics.simplify(components_canon(sym)[1, 2] - b), 0)
    @test isequal(TensND.Symbolics.simplify(voigt_strain(sym)[6] - 2b), 0)
end
