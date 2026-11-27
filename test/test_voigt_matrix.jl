# =============================================================================
#  test_voigt_matrix.jl — the 6×6 matrices of a stiffness and of a compliance
#  in the two-index (Voigt) notation (`voigt_stiffness`, `voigt_compliance`,
#  `stiffness_from_voigt`, `compliance_from_voigt`).
#
#  References, independent of the code under test:
#  - the published two-index matrix of the Mesaverde sandstone, a TI medium of
#    axis e₃, with C₆₆ = (C₁₁ − C₁₂)/2;
#  - the canonical components C_ijkl of a general tensor, read entry by entry
#    (C_IJ = C_ijkl, S_IJ = S_ijkl times 1, 2 or 4);
#  - the meaning of the matrices: σ = ℂ:ε computed by an explicit index sum on
#    a random symmetric strain, against the matrix product with engineering
#    shears (2ε₂₃, 2ε₁₃, 2ε₁₂), and the same for a compliance;
#  - the inverse of a stiffness, whose compliance matrix is the inverse matrix;
#  - with symbols, the entries C₄₄ and C₆₆ = (C₁₁ − C₁₂)/2 of a TI tensor; with
#    dual numbers, the derivative 1 of an entry of the round trip.
# =============================================================================

using StaticArrays: SMatrix

const VGT_SEED = 20261130
# Two-index order of the pairs: 11, 22, 33, 23, 13, 12.
const VGT_PAIRS = ((1, 1), (2, 2), (3, 3), (2, 3), (1, 3), (1, 2))
const VGT_SANDSTONE = (52.9, 15.7, 21.8, 47.8, 16.4)        # (C₁₁, C₁₂, C₁₃, C₃₃, C₄₄), GPa

vgt_ti(c) = ti_from_constants(c..., [0.0, 0.0, 1.0]; form = :components, role = :stiffness)

"A positive definite tensor without any symmetry class, written in a rotated basis."
vgt_general() = express(
    from_mandel(
        let M = randn(Xoshiro(1), 6, 6)
            M * transpose(M) + 3I
        end
    ), RotatedBasis(0.3, 0.7, 1.1)
)

"σ_ij = A_ijkl ε_kl by an explicit index sum."
vgt_contract(A, ε) = [sum(A[i, j, k, l] * ε[k, l] for k in 1:3, l in 1:3) for i in 1:3, j in 1:3]

"The six components of a symmetric 3×3 matrix in the two-index order, shears weighted by `w`."
vgt_vector(a, w) = [a[i, j] * (i == j ? 1 : w) for (i, j) in VGT_PAIRS]

@testsection "Two-index (Voigt) matrices" begin
    @testsection "a published TI matrix" begin
        C11, C12, C13, C33, C44 = VGT_SANDSTONE
        C66 = (C11 - C12) / 2
        ref = [
            C11 C12 C13 0 0 0
            C12 C11 C13 0 0 0
            C13 C13 C33 0 0 0
            0 0 0 C44 0 0
            0 0 0 0 C44 0
            0 0 0 0 0 C66
        ]
        C = vgt_ti(VGT_SANDSTONE)
        @test isapprox(voigt_stiffness(C), ref; atol = 1.0e-13)
        @test (@inferred voigt_stiffness(C)) isa SMatrix{6, 6, Float64}
        @test (@inferred voigt_compliance(inv(C))) isa SMatrix{6, 6, Float64}
        # The compliance matrix of the inverse tensor is the inverse matrix.
        @test isapprox(voigt_compliance(inv(C)) * voigt_stiffness(C), I(6); atol = 1.0e-15)
    end

    @testsection "entries of a general tensor in a rotated basis" begin
        G = vgt_general()
        A = Array(canonical_components(G))
        V, S = voigt_stiffness(G), voigt_compliance(G)
        for (p, (i, j)) in enumerate(VGT_PAIRS), (q, (k, l)) in enumerate(VGT_PAIRS)
            nshear = (p > 3) + (q > 3)
            @test isapprox(V[p, q], A[i, j, k, l]; atol = 1.0e-13)
            @test isapprox(S[p, q], 2^nshear * A[i, j, k, l]; atol = 1.0e-13)
        end
        @test V[4, 4] ≈ A[2, 3, 2, 3]
        @test V[1, 6] ≈ A[1, 1, 1, 2]
        @test V[5, 6] ≈ A[1, 3, 1, 2]
        @test S[4, 4] ≈ 4 * A[2, 3, 2, 3]
        @test S[1, 5] ≈ 2 * A[1, 1, 1, 3]
        @test S[1, 2] ≈ A[1, 1, 2, 2]
    end

    @testsection "stress and engineering strain" begin
        G = vgt_general()
        A = Array(canonical_components(G))
        rng = Xoshiro(VGT_SEED)
        for _ in 1:3
            ε = let e = randn(rng, 3, 3)
                (e + transpose(e)) / 2
            end
            σ = vgt_contract(A, ε)
            # σ_I = C_IJ ε_J with ε_J the engineering strains.
            @test isapprox(vgt_vector(σ, 1), voigt_stiffness(G) * vgt_vector(ε, 2); rtol = 1.0e-13)
            # ε_I = S_IJ σ_J: G read as a compliance maps a stress σ to the
            # strain G:σ, whose engineering components come out.
            τ = ε                                                 # a stress, this time
            @test isapprox(vgt_vector(vgt_contract(A, τ), 2), voigt_compliance(G) * vgt_vector(τ, 1); rtol = 1.0e-13)
        end
    end

    @testsection "round trips" begin
        G = vgt_general()
        A = Array(canonical_components(G))
        V, S = voigt_stiffness(G), voigt_compliance(G)
        @test isapprox(Array(canonical_components(stiffness_from_voigt(V))), A; atol = 1.0e-14 * maximum(abs, A))
        @test isapprox(Array(canonical_components(compliance_from_voigt(S))), A; atol = 1.0e-14 * maximum(abs, A))
        M = let m = randn(Xoshiro(VGT_SEED + 1), 6, 6)
            m + transpose(m)
        end
        @test isapprox(voigt_stiffness(stiffness_from_voigt(M)), M; atol = 1.0e-14 * maximum(abs, M))
        @test isapprox(voigt_compliance(compliance_from_voigt(M)), M; atol = 1.0e-14 * maximum(abs, M))
    end

    @testsection "a matrix that is not 6×6" begin
        @test_throws DimensionMismatch stiffness_from_voigt(zeros(5, 5))
        @test_throws DimensionMismatch compliance_from_voigt(zeros(6, 5))
    end

    @testsection "symbols" begin
        @syms c11::positive c12::positive c13::positive c33::positive c44::positive
        C = ti_from_constants(c11, c12, c13, c33, c44, [0, 0, 1]; form = :components)
        V = voigt_stiffness(C)
        @test tsimplify(V[4, 4] - c44) == 0
        @test tsimplify(V[6, 6] - (c11 - c12) / 2) == 0
        @test tsimplify(V[1, 3] - c13) == 0
        @test tsimplify(V[1, 4]) == 0
    end

    @testsection "dual numbers" begin
        M₀ = let m = randn(Xoshiro(VGT_SEED + 2), 6, 6)
            m * transpose(m) + 6I
        end
        # A symmetric perturbation of the entry (p, q): its own derivative is 1,
        # the others 0.
        for (p, q) in ((1, 1), (4, 6), (2, 5))
            E = zeros(6, 6)
            E[p, q] = E[q, p] = 1
            d = ForwardDiff.derivative(x -> voigt_stiffness(stiffness_from_voigt(M₀ + x * E)), 0.0)
            @test isapprox(d, E; atol = 1.0e-14)
            d = ForwardDiff.derivative(x -> voigt_compliance(compliance_from_voigt(M₀ + x * E)), 0.0)
            @test isapprox(d, E; atol = 1.0e-14)
        end
    end
end
