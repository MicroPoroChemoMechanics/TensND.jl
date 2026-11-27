# =============================================================================
#  test_ortho_constants.jl — the orthotropic tensor of the engineering
#  constants (`ortho_from_constants`, `ortho_constants`), in space and in the
#  plane.
#
#  References, independent of the code under test:
#  - the definitions read on a test: under a unit tension along the material
#    axis e_i, ε_ii = 1/E_i and ε_jj = −ν_ij/E_i; under a unit shear stress
#    σ_ij = σ_ji = 1, the engineering shear 2ε_ij = 1/G_ij. The strains are
#    computed from canonical components, the material frame being rotated;
#  - the transversely isotropic and isotropic tensors of the same constants
#    (`ti_from_constants`, `iso_compliance`);
#  - analytic derivatives of the compliance, and C₂₃₂₃ = G₂₃ for the stiffness.
# =============================================================================

const OCS_CONSTANTS = (12.0, 7.0, 4.5, 0.31, 0.22, 0.27, 1.9, 2.3, 3.1)   # E₁ E₂ E₃ ν₁₂ ν₁₃ ν₂₃ G₂₃ G₁₃ G₁₂
const OCS_PLANE = (12.0, 5.0, 0.07, 0.7)                                  # E₁ E₂ ν₁₂ G₁₂

"The strain of the unit stress e_i ⊗ e_j (symmetrized) in the material frame `e`, as material components."
function ocs_strain(S, e, i, j)
    A = canonical_components(S)
    d = size(e, 1)
    σ = (e[:, i] * e[:, j]' + e[:, j] * e[:, i]') / (i == j ? 2 : 1)
    ε = [sum(A[a, b, c, f] * σ[c, f] for c in 1:d, f in 1:d) for a in 1:d, b in 1:d]
    return e' * ε * e
end

@testsection "orthotropic engineering constants" begin
    @testsection "space: the constants read on tests" begin
        E₁, E₂, E₃, ν₁₂, ν₁₃, ν₂₃, G₂₃, G₁₃, G₁₂ = OCS_CONSTANTS
        for frame in (CanonicalBasis{3, Float64}(), RotatedBasis(0.3, 0.7, 1.1))
            e = Matrix(vecbasis(frame))
            S = @inferred ortho_from_constants(OCS_CONSTANTS..., frame)
            @test S isa TensOrtho{4, 3, Float64, 9}
            ε = ocs_strain(S, e, 1, 1)
            @test [ε[1, 1], ε[2, 2], ε[3, 3]] ≈ [1 / E₁, -ν₁₂ / E₁, -ν₁₃ / E₁]
            ε = ocs_strain(S, e, 2, 2)
            @test [ε[1, 1], ε[2, 2], ε[3, 3]] ≈ [-ν₁₂ / E₁, 1 / E₂, -ν₂₃ / E₂]
            ε = ocs_strain(S, e, 3, 3)
            @test ε[3, 3] ≈ 1 / E₃
            @test 2 * ocs_strain(S, e, 2, 3)[2, 3] ≈ 1 / G₂₃
            @test 2 * ocs_strain(S, e, 1, 3)[1, 3] ≈ 1 / G₁₃
            @test 2 * ocs_strain(S, e, 1, 2)[1, 2] ≈ 1 / G₁₂
            @test collect(@inferred ortho_constants(S)) ≈ collect(OCS_CONSTANTS)
            C = @inferred ortho_from_constants(OCS_CONSTANTS..., frame; role = :stiffness)
            @test mandel(C ⊡ S, CanonicalBasis{3, Float64}()) ≈ I
            @test collect(ortho_constants(C; role = :stiffness)) ≈ collect(OCS_CONSTANTS)
        end
    end

    @testsection "space: transversely isotropic and isotropic limits" begin
        E₁, E₃, ν₁₂, ν₃₁, G₃₁ = 9.0, 14.0, 0.28, 0.33, 2.7
        G₁₂ = E₁ / (2 * (1 + ν₁₂))
        frame = CanonicalBasis{3, Float64}()
        for role in (:compliance, :stiffness)
            O = ortho_from_constants(E₁, E₁, E₃, ν₁₂, ν₃₁ * E₁ / E₃, ν₃₁ * E₁ / E₃, G₃₁, G₃₁, G₁₂, frame; role)
            Tᵢ = ti_from_constants(E₁, E₃, ν₁₂, ν₃₁, G₃₁, [0.0, 0.0, 1.0]; form = :engineering, role)
            @test canonical_components(O) ≈ canonical_components(Tᵢ)
        end
        E, ν = 10.0, 0.25
        G = E / (2 * (1 + ν))
        O = ortho_from_constants(E, E, E, ν, ν, ν, G, G, G, frame)
        @test canonical_components(O) ≈ canonical_components(inv(iso_stiffness(E / (3 * (1 - 2ν)), G)))
    end

    @testsection "plane" begin
        E₁, E₂, ν₁₂, G₁₂ = OCS_PLANE
        for frame in (CanonicalBasis{2, Float64}(), Basis(0.37))
            e = Matrix(vecbasis(frame))
            S = @inferred ortho_from_constants(OCS_PLANE..., frame)
            @test S isa TensOrtho{4, 2, Float64, 4}
            ε = ocs_strain(S, e, 1, 1)
            @test [ε[1, 1], ε[2, 2]] ≈ [1 / E₁, -ν₁₂ / E₁]
            @test ocs_strain(S, e, 2, 2)[2, 2] ≈ 1 / E₂
            @test 2 * ocs_strain(S, e, 1, 2)[1, 2] ≈ 1 / G₁₂
            @test collect(@inferred ortho_constants(S)) ≈ collect(OCS_PLANE)
            C = ortho_from_constants(OCS_PLANE..., frame; role = :stiffness)
            @test collect(ortho_constants(C; role = :stiffness)) ≈ collect(OCS_PLANE)
            # Lekhnitskii's parameters of the sheet: p = −√(E₁/E₂) in its axes
            frame isa CanonicalBasis && @test characteristic_pair(C)[2] ≈ -sqrt(E₁ / E₂)
        end
    end

    @testsection "derivatives and symbols" begin
        frame = CanonicalBasis{3, Float64}()
        c = OCS_CONSTANTS
        S1111(E) = canonical_components(ortho_from_constants(E, c[2:end]..., frame))[1, 1, 1, 1]
        @test ForwardDiff.derivative(S1111, c[1]) ≈ -1 / c[1]^2
        @test ForwardDiff.derivative(E -> ForwardDiff.derivative(S1111, E), c[1]) ≈ 2 / c[1]^3
        C2323(G) = canonical_components(ortho_from_constants(c[1:6]..., G, c[8:9]..., frame; role = :stiffness))[2, 3, 2, 3]
        @test ForwardDiff.derivative(C2323, c[7]) ≈ 1
        E1_back(E) = ortho_constants(ortho_from_constants(E, c[2:end]..., frame; role = :stiffness); role = :stiffness)[1]
        @test ForwardDiff.derivative(E1_back, c[1]) ≈ 1
        @syms E₁::positive E₂::positive E₃::positive ν₁₂::real ν₁₃::real ν₂₃::real G₂₃::positive G₁₃::positive G₁₂::positive
        sym = (E₁, E₂, E₃, ν₁₂, ν₁₃, ν₂₃, G₂₃, G₁₃, G₁₂)
        S = ortho_from_constants(sym..., CanonicalBasis{3, Sym}())
        @test all(iszero ∘ simplify, collect(ortho_constants(S)) .- collect(sym))
        C = ortho_from_constants(sym..., CanonicalBasis{3, Sym}(); role = :stiffness)
        @test all(iszero ∘ simplify, collect(ortho_constants(C; role = :stiffness)) .- collect(sym))
        @test iszero(simplify(canonical_components(C)[1, 2, 1, 2] - G₁₂))
    end

    @testsection "refusals" begin
        frame = CanonicalBasis{3, Float64}()
        @test_throws ArgumentError ortho_from_constants(OCS_CONSTANTS..., frame; role = :modulus)
        general = TensOrtho(10.0, 9.0, 8.0, 3.0, 2.0, 1.0, 4.0, 5.0, 6.0, 3.5, 2.5, 1.5, frame)
        @test_throws ArgumentError ortho_constants(general)
    end
end
