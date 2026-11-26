# =============================================================================
#  test_iso_moduli_dimensions.jl — the isotropic moduli in any dimension, and
#  the two conventions of the plane.
#
#  The reference is the definition of each modulus, read on the full arrays of
#  the stiffness λ 𝟏⊗𝟏 + 2μ 𝕀 and of its inverse (computed here, from the
#  Kelvin-Mandel matrix): E = 1/S₁₁₁₁ and ν = −S₁₁₂₂/S₁₁₁₁ on a uniaxial
#  stress, k from tr σ = d k tr ε, λ = C₁₁₂₂, μ = C₁₂₁₂, M = C₁₁₁₁.
# =============================================================================

const IMD_NAMES = (:bulk, :young, :lame, :shear, :poisson, :pwave)
const IMD_READERS = (bulk_modulus, young_modulus, lame_lambda, shear_modulus, poisson_ratio, pwave_modulus)

"The six moduli of the medium of dimension `d` with Lamé coefficients `(λ, μ)`, from their definitions."
function imd_moduli(d, λ, μ)
    δ(i, j) = i == j ? 1 : 0
    C = [λ * δ(i, j) * δ(k, l) + μ * (δ(i, k) * δ(j, l) + δ(i, l) * δ(j, k)) for i in 1:d, j in 1:d, k in 1:d, l in 1:d]
    S = get_array(from_mandel(inv(mandel(Tens(C)))))
    return (
        bulk = sum(C[i, i, j, j] for i in 1:d, j in 1:d) / d^2, young = 1 / S[1, 1, 1, 1], lame = C[1, 1, 2, 2],
        shear = C[1, 2, 1, 2], poisson = -S[1, 1, 2, 2] / S[1, 1, 1, 1], pwave = C[1, 1, 1, 1],
    )
end
imd_stiffness(d, λ, μ) = TensISO{d}(d * λ + 2μ, 2μ)
imd_builder(kind, p, q) = getfield(TensND, (p, q) == (:bulk, :shear) ? Symbol(:iso_, kind) : Symbol(:iso_, kind, :_, p, :_, q))
# (k, μ) is the same in both conventions: its builder takes no keyword.
imd_build(kind, p, q, v, a, b, convention) =
    (p, q) == (:bulk, :shear) ? imd_builder(kind, p, q)(v, a, b) : imd_builder(kind, p, q)(v, a, b; convention)

imd_inferred_build(v, E, ν) = @inferred iso_stiffness_young_poisson(v, E, ν; convention = :intrinsic)
imd_inferred_default(v, E, ν) = @inferred iso_stiffness_bulk_young(v, E, ν)
imd_inferred_read(C) = @inferred young_poisson(C; convention = :intrinsic)
imd_inferred_compliance(S) = @inferred young_modulus(S; role = :compliance, convention = :intrinsic)

@testsection "Isotropic moduli in any dimension" begin

    @testsection "the medium of dimension d: every pair, against the definitions" begin
        for d in 2:5, ν in (-0.8, -0.3, 0.0, 0.1, 0.2)
            v, μ = Val(d), 7.0
            λ = 2μ * ν / (1 - (d - 1) * ν)
            ref = imd_moduli(d, λ, μ)
            @test ref.poisson ≈ ν atol = 1.0e-12
            C = imd_stiffness(d, λ, μ)
            S = inv(C)
            for (n, r) in zip(IMD_NAMES, IMD_READERS)
                @test r(C; convention = :intrinsic) ≈ ref[n] atol = 1.0e-11 * μ
                @test r(S; role = :compliance, convention = :intrinsic) ≈ ref[n] atol = 1.0e-11 * μ
            end
            for i in 1:6, j in (i + 1):6
                p, q = IMD_NAMES[i], IMD_NAMES[j]
                pair = getfield(TensND, Symbol(p, :_, q))
                @test collect(pair(C; convention = :intrinsic)) ≈ [ref[p], ref[q]] atol = 1.0e-11 * μ
                @test collect(pair(S; role = :compliance, convention = :intrinsic)) ≈ [ref[p], ref[q]] atol = 1.0e-11 * μ
                # (λ, ν) does not give μ at ν = 0; (E, M) keeps the root ν ≥ 0.
                (p, q) == (:lame, :poisson) && ν == 0 && continue
                (p, q) == (:young, :pwave) && ν < 0 && continue
                Cb = imd_build(:stiffness, p, q, v, ref[p], ref[q], :intrinsic)
                @test Cb isa TensISO{4, d, Float64}
                @test collect(get_data(Cb)) ≈ collect(get_data(C)) rtol = 1.0e-9
                Sb = imd_build(:compliance, p, q, v, ref[p], ref[q], :intrinsic)
                @test collect(get_data(Sb)) ≈ collect(get_data(S)) rtol = 1.0e-9
            end
        end
    end

    @testsection "the two conventions" begin
        # Space: one convention, and the keyword changes nothing, to the bit.
        E, ν = 70.0, 0.33
        @test iso_stiffness_young_poisson(Val(3), E, ν; convention = :intrinsic) === iso_stiffness_young_poisson(E, ν)
        @test iso_compliance_young_poisson(Val(3), E, ν; convention = :intrinsic) === iso_compliance_young_poisson(E, ν)
        C₃ = iso_stiffness_young_poisson(E, ν)
        for i in 1:6, j in (i + 1):6
            p, q = IMD_NAMES[i], IMD_NAMES[j]
            pair = getfield(TensND, Symbol(p, :_, q))
            @test pair(C₃; convention = :intrinsic) === pair(C₃)
            a, b = pair(C₃)
            @test imd_build(:stiffness, p, q, Val(3), a, b, :intrinsic) === imd_builder(:stiffness, p, q)(Val(3), a, b)
            @test imd_build(:compliance, p, q, Val(3), a, b, :intrinsic) === imd_builder(:compliance, p, q)(Val(3), a, b)
        end
        # The plane: plane strain by default, to the bit; the plane medium on demand.
        C₂ = iso_stiffness_young_poisson(Val(2), E, ν)
        @test iso_stiffness_young_poisson(Val(2), E, ν; convention = :plane_strain) === C₂
        @test young_poisson(C₂; convention = :plane_strain) === young_poisson(C₂)
        @test collect(young_poisson(C₂)) ≈ [E, ν]
        # E₂ and ν₂ of the plane medium, from those of the 3-D material in plane strain.
        E₂, ν₂ = young_poisson(C₂; convention = :intrinsic)
        @test E₂ ≈ E / (1 - ν^2) && ν₂ ≈ ν / (1 - ν)
        @test collect(plane_strain_to_plane_stress(E, ν)) ≈ [E₂, ν₂]
        k, μ = bulk_shear(C₂)
        @test E₂ ≈ 4k * μ / (k + μ) && ν₂ ≈ (k - μ) / (k + μ)
        @test iso_stiffness_young_poisson(Val(2), E₂, ν₂; convention = :intrinsic) ≈ C₂
        @test !(iso_stiffness_young_poisson(Val(2), E, ν; convention = :intrinsic) ≈ C₂)
        # k, λ, μ and M do not see the convention.
        for r in (bulk_modulus, lame_lambda, shear_modulus, pwave_modulus)
            @test r(C₂; convention = :intrinsic) === r(C₂)
        end
        @test iso_stiffness_lame_shear(Val(2), 3.0, 2.0; convention = :intrinsic) === iso_stiffness_lame_shear(Val(2), 3.0, 2.0)
        # Beyond dimension 3 the medium is the default.
        C₄ = iso_stiffness_young_poisson(Val(4), E, 0.2)
        @test C₄ === iso_stiffness_young_poisson(Val(4), E, 0.2; convention = :intrinsic)
        @test collect(young_poisson(C₄)) ≈ [E, 0.2]
    end

    @testsection "refusals" begin
        for d in (3, 4)
            C = TensISO{d}(3.0, 2.0)
            @test_throws ArgumentError young_modulus(C; convention = :plane_strain)
            @test_throws ArgumentError bulk_modulus(C; convention = :plane_strain)
            @test_throws ArgumentError young_poisson(C; convention = :plane_strain)
            @test_throws ArgumentError iso_stiffness_young_poisson(Val(d), 1.0, 0.2; convention = :plane_strain)
            @test_throws ArgumentError iso_compliance_young_poisson(Val(d), 1.0, 0.2; convention = :plane_strain)
            @test_throws ArgumentError iso_stiffness_bulk_young(Val(d), 1.0, 0.2; convention = :plane_strain)
            @test_throws ArgumentError iso_stiffness_lame_shear(Val(d), 1.0, 0.2; convention = :plane_strain)
        end
        @test_throws ArgumentError young_modulus(TensISO{2}(3.0, 2.0); convention = :plane_stress)
        @test_throws ArgumentError iso_stiffness_young_shear(Val(2), 1.0, 0.4; convention = :plane)
        @test occursin(
            "intrinsic", sprint(
                showerror, try
                    young_modulus(TensISO{4}(3.0, 2.0); convention = :plane_strain)
                catch e
                    e
                end
            )
        )
    end

    @testsection "inference, derivatives, symbols" begin
        for v in (Val(2), Val(4))
            @test imd_inferred_build(v, 70.0, 0.2) isa TensISO{4}
            @test imd_inferred_default(v, 5.0, 3.0) isa TensISO{4}
            @test imd_inferred_read(imd_inferred_build(v, 70.0, 0.2)) isa NTuple{2, Float64}
            @test imd_inferred_compliance(inv(imd_inferred_build(v, 70.0, 0.2))) ≈ 70.0
        end
        # First and second derivatives of E_d(k, μ) = 2μ d²k / (2μ + d(d−1)k).
        d = 4
        E_closed(x) = 2x[2] * d^2 * x[1] / (2x[2] + d * (d - 1) * x[1])
        E_read(x) = young_modulus(iso_stiffness(Val(d), x[1], x[2]))
        x₀ = [5.0, 3.0]
        @test ForwardDiff.gradient(E_read, x₀) ≈ ForwardDiff.gradient(E_closed, x₀)
        @test ForwardDiff.hessian(E_read, x₀) ≈ ForwardDiff.hessian(E_closed, x₀)
        # A stiffness is linear in E at fixed ν.
        C = iso_stiffness_young_poisson(Val(d), 70.0, 0.2)
        dC = ForwardDiff.derivative(E -> collect(get_data(iso_stiffness_young_poisson(Val(d), E, 0.2))), 70.0)
        @test dC ≈ collect(get_data(C)) ./ 70.0
        # Through the two quadratic conversions.
        E₀, λ₀, M₀ = young_modulus(C), lame_lambda(C), pwave_modulus(C)
        @test ForwardDiff.derivative(t -> shear_modulus(iso_stiffness_young_lame(Val(d), t * E₀, t * λ₀)), 1.0) ≈ shear_modulus(C)
        @test ForwardDiff.derivative(t -> shear_modulus(iso_stiffness_young_pwave(Val(d), t * E₀, t * M₀)), 1.0) ≈ shear_modulus(C)

        # With symbols: the compliance built from (E, ν) has S₁₁₁₁ = 1/E and S₁₁₂₂ = −ν/E.
        E, ν = symbols("E ν", positive = true)
        for (v, convention) in ((Val(2), :intrinsic), (Val(4), :intrinsic), (Val(5), :intrinsic), (Val(3), :intrinsic))
            S = get_array(iso_compliance_young_poisson(v, E, ν; convention))
            @test simplify(S[1, 1, 1, 1] - 1 / E) == 0
            @test simplify(S[1, 1, 2, 2] + ν / E) == 0
            C = iso_stiffness_young_poisson(v, E, ν; convention)
            @test all(iszero, simplify.(collect(young_poisson(C; convention)) .- [E, ν]))
            @test all(iszero, simplify.(collect(young_poisson(inv(C); role = :compliance, convention)) .- [E, ν]))
        end
        # In plane strain it is S₁₁₁₁ = (1 − ν²)/E: the strain e₃₃ is blocked.
        S₂ = get_array(iso_compliance_young_poisson(Val(2), E, ν))
        @test simplify(S₂[1, 1, 1, 1] - (1 - ν^2) / E) == 0
        k, μ = symbols("k μ", positive = true)
        C₂ = iso_stiffness(Val(2), k, μ)
        @test simplify(young_modulus(C₂; convention = :intrinsic) - 4k * μ / (k + μ)) == 0
        @test simplify(poisson_ratio(C₂; convention = :intrinsic) - (k - μ) / (k + μ)) == 0
        C₄ = iso_stiffness(Val(4), k, μ)
        @test simplify(poisson_ratio(C₄) - (4k - 2μ) / (12k + 2μ)) == 0
        @test simplify(lame_lambda(C₄) - (k - μ / 2)) == 0
        # Symbolics.Num goes through the rational conversions.
        Symbolics.@variables En νn
        @test iso_stiffness_young_poisson(Val(4), En, νn) isa TensISO{4, 4, Symbolics.Num}
        @test young_poisson(iso_stiffness_young_poisson(Val(2), En, νn; convention = :intrinsic); convention = :intrinsic) isa NTuple{2, Symbolics.Num}
    end
end
