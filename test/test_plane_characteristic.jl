# =============================================================================
#  test_plane_characteristic.jl — the characteristic roots of a stiffness of
#  the plane (`engineering_compliances`, `characteristic_pair`,
#  `characteristic_roots`).
#
#  References, independent of the code under test:
#  - the equilibrium of u = U f(x₁ + μ x₂), written with the stiffness:
#    det[Q + μ(R + Rᵀ) + μ²T] = 0, Q = e₁·ℂ·e₁, R = e₁·ℂ·e₂, T = e₂·ℂ·e₂,
#    by sums over the indices (the code works on the compliance);
#  - the compliances read as strains under the three unit stresses;
#  - the closed form in the axes of an orthotropic sheet of engineering
#    constants (E₁, E₂, G, ν₁): s = i n, p = −k with k = √(E₁/E₂),
#    n = √(2k + E₁/G − 2ν₁) (Lekhnitskii 1981, Eq. (31.1));
#  - the closed forms of the orthotropic and square classes against the
#    numeric route on their general copies, values and derivatives.
# =============================================================================

const PCH_SY = TensND.Symbolics
const PCH_SEED = 20261127

pch_general() = from_mandel([30.0 4.0 3.0; 4.0 9.0 -2.0; 3.0 -2.0 7.0])
pch_ortho(β = 0.37, C₁₁ = 40.0) = TensOrtho(C₁₁, 6.0, 2.5, 3.0, Basis(β))
pch_square(β = 1.1) = cubic_from_constants(20.0, 12.0, 9.0, Basis(β))
pch_iso() = iso_stiffness(Val(2), 5.0, 2.0)
# The same tensor as a general one: the numeric route.
pch_copy(C) = Tens(Array(canonical_components(C)))

"det[Q + μ(R + Rᵀ) + μ²T] from the canonical components of the stiffness."
function pch_equilibrium(C, μ)
    A = canonical_components(C)
    d = (1, μ)
    M = [sum(d[j] * A[i, j, k, l] * d[l] for j in 1:2, l in 1:2) for i in 1:2, k in 1:2]
    return M[1, 1] * M[2, 2] - M[1, 2] * M[2, 1]
end

@testsection "characteristic roots of a plane stiffness" begin
    rng = Random.Xoshiro(PCH_SEED)

    @testsection "engineering compliances" begin
        C = pch_general()
        a₁₁, a₂₂, a₁₂, a₆₆, a₁₆, a₂₆ = engineering_compliances(C)
        strain(S) = canonical_components(inv(C) ⊡ Tens(S))
        ε = strain([1.0 0.0; 0.0 0.0])
        @test [ε[1, 1], ε[2, 2], 2ε[1, 2]] ≈ [a₁₁, a₁₂, a₁₆] rtol = 1.0e-12
        ε = strain([0.0 0.0; 0.0 1.0])
        @test [ε[1, 1], ε[2, 2], 2ε[1, 2]] ≈ [a₁₂, a₂₂, a₂₆] rtol = 1.0e-12
        ε = strain([0.0 1.0; 1.0 0.0])
        @test [ε[1, 1], ε[2, 2], 2ε[1, 2]] ≈ [a₁₆, a₂₆, a₆₆] rtol = 1.0e-12
    end

    @testsection "the roots solve the equilibrium" begin
        for C in (pch_general(), pch_ortho(), pch_square(), pch_iso(), pch_ortho(2π * rand(rng)), pch_square(2π * rand(rng)))
            s, p = characteristic_pair(C)
            μ₁, μ₂ = characteristic_roots(C)
            @test μ₁ + μ₂ ≈ s && μ₁ * μ₂ ≈ p
            @test imag(μ₁) > 0 && imag(μ₂) > 0
            scale = abs(pch_equilibrium(C, 1.0im))
            # the quartic is proportional to (μ² − sμ + p)(μ² − s̄μ + p̄)
            z = randn(rng, ComplexF64, 4)
            q(μ) = (μ^2 - s * μ + p) * (μ^2 - conj(s) * μ + conj(p))
            @test [pch_equilibrium(C, μ) / pch_equilibrium(C, z[1]) for μ in z] ≈ [q(μ) / q(z[1]) for μ in z] rtol = 1.0e-10
            C isa TensISO || @test abs(pch_equilibrium(C, μ₁)) + abs(pch_equilibrium(C, μ₂)) ≤ 1.0e-10 * scale
        end
    end

    @testsection "closed forms" begin
        @test characteristic_pair(pch_iso()) === (complex(0.0, 2.0), complex(-1.0))
        @test all(≈(1.0im), characteristic_roots(pch_iso()))
        for β in (0.0, 2π * rand(rng), 2π * rand(rng))
            for C in (pch_ortho(β), pch_square(β))
                s, p = characteristic_pair(C)
                s₀, p₀ = characteristic_pair(pch_copy(C))
                @test s ≈ s₀ rtol = 1.0e-12
                @test p ≈ p₀ rtol = 1.0e-12
            end
        end
        # an orthotropic sheet in its axes, from its engineering constants
        E₁, E₂, G, ν₁ = 12.0, 5.0, 0.7, 0.07
        S = [1 / E₁ -ν₁ / E₁ 0; -ν₁ / E₁ 1 / E₂ 0; 0 0 1 / (2G)]
        k, n = sqrt(E₁ / E₂), sqrt(2 * sqrt(E₁ / E₂) + E₁ / G - 2ν₁)
        for C in (inv(from_mandel(S)), TensOrtho(inv(S)[1, 1], inv(S)[2, 2], inv(S)[1, 2], G, Basis(0.0)))
            s, p = characteristic_pair(C)
            @test s ≈ im * n rtol = 1.0e-12
            @test p ≈ -k rtol = 1.0e-12
        end
    end

    @testsection "derivatives" begin
        # numeric route against the closed form, with respect to a modulus and to the angle
        for x₀ in (40.0, 55.0)
            closed(x) = collect(reim(characteristic_pair(pch_ortho(0.37, x))[1]))
            numeric(x) = collect(reim(characteristic_pair(pch_copy(pch_ortho(0.37, x)))[1]))
            @test ForwardDiff.derivative(numeric, x₀) ≈ ForwardDiff.derivative(closed, x₀) rtol = 1.0e-9
            d²(f) = ForwardDiff.derivative(x -> ForwardDiff.derivative(f, x), x₀)
            @test d²(numeric) ≈ d²(closed) rtol = 1.0e-7
        end
        # in its axes p = −√(a₂₂/a₁₁): ∂p/∂a₂₂ = p/(2a₂₂)
        p_of(a₂₂) = real(characteristic_pair(inv(from_mandel([0.1 -0.01 0; -0.01 a₂₂ 0; 0 0 0.6])))[2])
        @test ForwardDiff.derivative(p_of, 0.3) ≈ -sqrt(0.3 / 0.1) / (2 * 0.3) rtol = 1.0e-10
        # through isotropy, where the two roots merge: the pair stays smooth
        iso_pair(x) = collect(reim(characteristic_pair(pch_copy(TensOrtho(x, 9.0, 5.0, 2.0, Basis(0.0))))[2]))
        ortho_pair(x) = collect(reim(characteristic_pair(TensOrtho(x, 9.0, 5.0, 2.0, Basis(0.0)))[2]))
        @test ForwardDiff.derivative(iso_pair, 9.0) ≈ ForwardDiff.derivative(ortho_pair, 9.0) rtol = 1.0e-9
    end

    @testsection "symbols" begin
        PCH_SY.@variables c₁₁ c₂₂ c₁₂ c₆₆ k μ
        at = Dict(c₁₁ => 40.0, c₂₂ => 6.0, c₁₂ => 2.5, c₆₆ => 3.0)
        number(x) = Float64(PCH_SY.value(PCH_SY.substitute(x, at; fold = Val(true))))
        number(z::Complex) = complex(number(real(z)), number(imag(z)))
        for (frame, numeric_frame) in ((CanonicalBasis{2, PCH_SY.Num}(), CanonicalBasis{2, Float64}()), (Basis(0.37), Basis(0.37)))
            s, p = characteristic_pair(TensOrtho(c₁₁, c₂₂, c₁₂, c₆₆, frame))
            s₀, p₀ = characteristic_pair(pch_copy(TensOrtho(40.0, 6.0, 2.5, 3.0, numeric_frame)))
            @test number(s) ≈ s₀ rtol = 1.0e-12
            @test number(p) ≈ p₀ rtol = 1.0e-12
        end
        s, p = characteristic_pair(iso_stiffness(Val(2), k, μ))
        @test isequal(s, complex(PCH_SY.Num(0), PCH_SY.Num(2))) && isequal(p, complex(PCH_SY.Num(-1)))
        general = Tens(canonical_components(TensOrtho(c₁₁, c₂₂, c₁₂, c₆₆, CanonicalBasis{2, PCH_SY.Num}())))
        @test_throws ArgumentError characteristic_pair(general)
    end

    @testsection "inference" begin
        @inferred characteristic_pair(pch_general())
        @inferred characteristic_pair(pch_ortho())
        @inferred characteristic_pair(pch_iso())
        @inferred characteristic_roots(pch_general())
        @inferred engineering_compliances(pch_general())
    end

    @testsection "refusals" begin
        A = Array(canonical_components(pch_general()))
        A[1, 1, 2, 2] += 1.5
        @test_throws ArgumentError characteristic_pair(Tens(A))
        # μ⁴ + 2μ² − 1 has two real roots
        @test_throws DomainError characteristic_pair(from_mandel([1.0 0 0; 0 -1.0 0; 0 0 1.0]))
    end
end
