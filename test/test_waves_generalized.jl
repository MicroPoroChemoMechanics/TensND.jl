# =============================================================================
#  test_waves_generalized.jl — the generalized Thomsen parameters of a
#  stiffness of any symmetry (`generalized_thomsen_parameters`) and the unit
#  direction of two angles (`direction_from_angles`).
#
#  References, independent of the code under test:
#  - the Lavoux limestone of Mensch & Rasolofosaon (1997, Geophys. J. Int.
#    128, Appendix B), a triclinic stiffness in two-index notation (GPa), and
#    its seventeen parameters recomputed apart, to four digits;
#  - the two-index entries C_IJ read on the canonical components C_ijkl, and
#    the definitions of the parameters written on them;
#  - the meaning of the parameters: the P modulus n·Γ(n)·n = C_ijkl n_i n_j n_k n_l,
#    an explicit index sum, written as the exact polynomial
#        C_ijkl n_i n_j n_k n_l / C₃₃ = 1 + 2ε_y n₁⁴ + 2ε_x n₂⁴
#            + 2δ_z n₁²n₂² + 2δ_y n₁²n₃² + 2δ_x n₂²n₃²
#            + 4(ε₁₅ n₁³n₃ + ε₁₆ n₁³n₂ + ε₂₄ n₂³n₃ + ε₂₆ n₂³n₁ + ε₃₄ n₃³n₂ + ε₃₅ n₃³n₁)
#            + 4(χ_x n₁²n₂n₃ + χ_y n₁n₂²n₃ + χ_z n₁n₂n₃²)
#    for any unit n (derived below from (n₁² + n₂² + n₃²)² = 1), against the
#    P wave of `weak_plane_waves`, ρV_P² = n·Γ(n)·n;
#  - for a TI medium of axis e₃, Thomsen's (1986) ε, γ and δ written from the
#    components, and the exact identity δ = δ_x + δ_x²/(2(1 − C₄₄/C₃₃))
#    (Sayers) between his δ and δ_x;
#  - for `direction_from_angles`, the spherical coordinates written out, the
#    converse `angles`, and the derivatives of sin and cos.
#
#  The polynomial. With m = (n₁², n₂², n₃², 2n₂n₃, 2n₁n₃, 2n₁n₂) the
#  engineering form of n⊗n, C_ijkl n_i n_j n_k n_l = mᵀ C m, that is
#      C₁₁n₁⁴ + C₂₂n₂⁴ + C₃₃n₃⁴ + 2(C₁₂ + 2C₆₆)n₁²n₂² + 2(C₁₃ + 2C₅₅)n₁²n₃²
#      + 2(C₂₃ + 2C₄₄)n₂²n₃² + 4C₁₅n₁³n₃ + 4C₁₆n₁³n₂ + 4C₂₄n₂³n₃ + 4C₂₆n₂³n₁
#      + 4C₃₄n₃³n₂ + 4C₃₅n₃³n₁ + 4(C₁₄ + 2C₅₆)n₁²n₂n₃ + 4(C₂₅ + 2C₄₆)n₁n₂²n₃
#      + 4(C₃₆ + 2C₄₅)n₁n₂n₃² ;
#  subtracting C₃₃ (n₁⁴ + n₂⁴ + n₃⁴ + 2n₁²n₂² + 2n₁²n₃² + 2n₂²n₃²) = C₃₃ and
#  dividing by C₃₃ gives the polynomial of the parameters above.
# =============================================================================

using StaticArrays: SVector

const WGN_SEED = 20261131
const WGN_PAIRS = ((1, 1), (2, 2), (3, 3), (2, 3), (1, 3), (1, 2))
# Lavoux limestone, Mensch & Rasolofosaon (1997), Appendix B, GPa.
const WGN_LAVOUX = [
    24.1 11.7 13.5 -0.5 0.6 -0.5
    11.7 23.5 12.6 -0.2 0.0 1.1
    13.5 12.6 26.5 0.1 -0.4 -0.2
    -0.5 -0.2 0.1 5.9 0.0 0.2
    0.6 0.0 -0.4 0.0 5.2 -0.3
    -0.5 1.1 -0.2 0.2 -0.3 5.5
]
# Its parameters, recomputed apart and rounded to four digits.
const WGN_LAVOUX_PARAMETERS = (
    εx = -0.0566, εy = -0.0453, δx = -0.0792, δy = -0.0981, δz = -0.1434,
    χx = -0.0415, χy = 0.0151, χz = -0.0075,
    ε15 = 0.0226, ε16 = -0.0189, ε24 = -0.0075, ε26 = 0.0415, ε34 = 0.0038, ε35 = -0.0151,
    γx = 0.0673, γy = -0.0339, γz = -0.0593,
)
const WGN_RHO = 2.19                 # Lavoux limestone, g/cm³
const WGN_SANDSTONE = (52.9, 15.7, 21.8, 47.8, 16.4)        # (C₁₁, C₁₂, C₁₃, C₃₃, C₄₄), GPa

wgn_ti(c) = ti_from_constants(c..., [0.0, 0.0, 1.0]; form = :components, role = :stiffness)

"The two-index entries C_IJ = C_ijkl, read on the canonical components."
wgn_entries(t) = (A = Array(canonical_components(t)); [A[i, j, k, l] for (i, j) in WGN_PAIRS, (k, l) in WGN_PAIRS])

"The parameters from their definitions on the two-index entries `C`."
function wgn_definitions(C)
    M = C[3, 3]
    return (
        εx = (C[2, 2] - M) / (2M), εy = (C[1, 1] - M) / (2M),
        δx = (C[2, 3] + 2C[4, 4] - M) / M, δy = (C[1, 3] + 2C[5, 5] - M) / M, δz = (C[1, 2] + 2C[6, 6] - M) / M,
        χx = (C[1, 4] + 2C[5, 6]) / M, χy = (C[2, 5] + 2C[4, 6]) / M, χz = (C[3, 6] + 2C[4, 5]) / M,
        ε15 = C[1, 5] / M, ε16 = C[1, 6] / M, ε24 = C[2, 4] / M, ε26 = C[2, 6] / M, ε34 = C[3, 4] / M, ε35 = C[3, 5] / M,
        γx = (C[4, 4] - C[5, 5]) / (2C[5, 5]), γy = (C[6, 6] - C[4, 4]) / (2C[4, 4]), γz = (C[5, 5] - C[4, 4]) / (2C[4, 4]),
    )
end

"C_ijkl n_i n_j n_k n_l by an explicit index sum."
wgn_quartic(A, n) = sum(A[i, j, k, l] * n[i] * n[j] * n[k] * n[l] for i in 1:3, j in 1:3, k in 1:3, l in 1:3)

"The same quartic form written with the two-index entries (header)."
function wgn_voigt_quartic(C, n)
    n₁, n₂, n₃ = n
    return C[1, 1] * n₁^4 + C[2, 2] * n₂^4 + C[3, 3] * n₃^4 +
        2 * (C[1, 2] + 2C[6, 6]) * n₁^2 * n₂^2 + 2 * (C[1, 3] + 2C[5, 5]) * n₁^2 * n₃^2 +
        2 * (C[2, 3] + 2C[4, 4]) * n₂^2 * n₃^2 +
        4 * (C[1, 5] * n₁^3 * n₃ + C[1, 6] * n₁^3 * n₂ + C[2, 4] * n₂^3 * n₃ + C[2, 6] * n₂^3 * n₁ + C[3, 4] * n₃^3 * n₂ + C[3, 5] * n₃^3 * n₁) +
        4 * ((C[1, 4] + 2C[5, 6]) * n₁^2 * n₂ * n₃ + (C[2, 5] + 2C[4, 6]) * n₁ * n₂^2 * n₃ + (C[3, 6] + 2C[4, 5]) * n₁ * n₂ * n₃^2)
end

"The quartic form divided by C₃₃, written with the parameters `p` (header)."
function wgn_parameter_quartic(p, n)
    n₁, n₂, n₃ = n
    return 1 + 2 * p.εy * n₁^4 + 2 * p.εx * n₂^4 +
        2 * p.δz * n₁^2 * n₂^2 + 2 * p.δy * n₁^2 * n₃^2 + 2 * p.δx * n₂^2 * n₃^2 +
        4 * (p.ε15 * n₁^3 * n₃ + p.ε16 * n₁^3 * n₂ + p.ε24 * n₂^3 * n₃ + p.ε26 * n₂^3 * n₁ + p.ε34 * n₃^3 * n₂ + p.ε35 * n₃^3 * n₁) +
        4 * (p.χx * n₁^2 * n₂ * n₃ + p.χy * n₁ * n₂^2 * n₃ + p.χz * n₁ * n₂ * n₃^2)
end

"A unit direction of colatitude θ and azimuth λ, written out."
wgn_direction(θ, λ) = [cos(λ) * sin(θ), sin(λ) * sin(θ), cos(θ)]

"Random unit directions."
wgn_directions(rng, m) = [wgn_direction(π * rand(rng), 2π * rand(rng)) for _ in 1:m]

"Checks the P-wave polynomial of the parameters on the stiffness `C`, at random directions."
function wgn_check_polynomial(C, rng)
    A = Array(canonical_components(C))
    V = wgn_entries(C)
    p = generalized_thomsen_parameters(C)
    for n in wgn_directions(rng, 8)
        q = wgn_quartic(A, n)
        # The two-index polynomial is the index sum (the derivation of the header).
        @test isapprox(wgn_voigt_quartic(V, n), q; rtol = 1.0e-13)
        # The parameters give the same, and so does the P wave of weak anisotropy.
        @test isapprox(V[3, 3] * wgn_parameter_quartic(p, n), q; rtol = 1.0e-13)
        @test isapprox(WGN_RHO * weak_plane_waves(C, WGN_RHO, n).qP[1]^2, q; rtol = 1.0e-13)
    end
    return nothing
end

@testsection "Generalized Thomsen parameters" begin
    @testsection "the Lavoux limestone" begin
        C = stiffness_from_voigt(WGN_LAVOUX)
        p = @inferred generalized_thomsen_parameters(C)
        @test keys(p) == keys(WGN_LAVOUX_PARAMETERS)
        for k in keys(p)
            @test isapprox(p[k], WGN_LAVOUX_PARAMETERS[k]; atol = 5.0e-5)
        end
        # From the definitions, on the published matrix and on the entries read
        # on the canonical components.
        for ref in (wgn_definitions(WGN_LAVOUX), wgn_definitions(wgn_entries(C)))
            for k in keys(p)
                @test isapprox(p[k], ref[k]; atol = 1.0e-14)
            end
        end
    end

    @testsection "the P wave of a triclinic medium" begin
        wgn_check_polynomial(stiffness_from_voigt(WGN_LAVOUX), Xoshiro(WGN_SEED))
    end

    @testsection "an orthotropic medium in its axes" begin
        C = TensOrtho(30.0, 26.0, 22.0, 9.0, 8.0, 7.5, 6.0, 7.0, 8.0, CanonicalBasis{3, Float64}())
        p = generalized_thomsen_parameters(C)
        # No odd term: the P modulus is even in each component of n.
        for k in (:χx, :χy, :χz, :ε15, :ε16, :ε24, :ε26, :ε34, :ε35)
            @test abs(p[k]) ≤ 1.0e-15
        end
        wgn_check_polynomial(C, Xoshiro(WGN_SEED + 1))
    end

    @testsection "a TI medium: Thomsen's parameters" begin
        C11, C12, C13, C33, C44 = WGN_SANDSTONE
        C66 = (C11 - C12) / 2
        ε = (C11 - C33) / (2C33)
        γ = (C66 - C44) / (2C44)
        δ = ((C13 + C44)^2 - (C33 - C44)^2) / (2C33 * (C33 - C44))
        p = generalized_thomsen_parameters(wgn_ti(WGN_SANDSTONE))
        @test p.εx ≈ ε
        @test p.εy ≈ ε
        @test p.γy ≈ γ
        @test abs(p.γx) ≤ 1.0e-15
        @test abs(p.γz) ≤ 1.0e-15
        @test p.δx ≈ p.δy
        @test p.δx ≈ (C13 + 2C44 - C33) / C33
        @test isapprox(p.δx, 0.14225941422594138; rtol = 1.0e-14)
        # Thomsen's δ = δ_x to first order, by an exact identity.
        @test isapprox(δ, p.δx + p.δx^2 / (2 * (1 - C44 / C33)); rtol = 1.0e-14)
        for k in (:χx, :χy, :χz, :ε15, :ε16, :ε24, :ε26, :ε34, :ε35)
            @test abs(p[k]) ≤ 1.0e-15
        end
    end

    @testsection "symbols" begin
        @syms c11::positive c12::positive c13::positive c33::positive c44::positive
        C = ti_from_constants(c11, c12, c13, c33, c44, [0, 0, 1]; form = :components)
        p = generalized_thomsen_parameters(C)
        @test tsimplify(p.εy - (c11 - c33) / (2 * c33)) == 0
        @test tsimplify(p.εx - p.εy) == 0
        @test tsimplify(p.δx - (c13 + 2 * c44 - c33) / c33) == 0
        @test tsimplify(p.γy - ((c11 - c12) / 2 - c44) / (2 * c44)) == 0
        @test tsimplify(p.χx) == 0
    end

    @testsection "dual numbers" begin
        E(p, q) = (e = zeros(6, 6); e[p, q] = e[q, p] = 1; e)
        params(x, e) = generalized_thomsen_parameters(stiffness_from_voigt(WGN_LAVOUX + x * e))
        M, C11, C23, C44 = WGN_LAVOUX[3, 3], WGN_LAVOUX[1, 1], WGN_LAVOUX[2, 3], WGN_LAVOUX[4, 4]
        # εy = (C₁₁ − M)/(2M): ∂/∂C₁₁ = 1/(2M), ∂/∂M = −C₁₁/(2M²); δx: ∂/∂M = −(C₂₃ + 2C₄₄)/M².
        @test isapprox(ForwardDiff.derivative(x -> params(x, E(1, 1)).εy, 0.0), 1 / (2M); rtol = 1.0e-13)
        @test isapprox(ForwardDiff.derivative(x -> params(x, E(3, 3)).εy, 0.0), -C11 / (2M^2); rtol = 1.0e-13)
        @test isapprox(ForwardDiff.derivative(x -> params(x, E(3, 3)).δx, 0.0), -(C23 + 2C44) / M^2; rtol = 1.0e-13)
    end
end

@testsection "Direction from angles" begin
    @testsection "numbers and dual numbers" begin
        rng = Xoshiro(WGN_SEED + 2)
        for _ in 1:5
            θ, ϕ = π * (0.05 + 0.9 * rand(rng)), π * (1.9 * rand(rng) - 0.95)
            n = @inferred direction_from_angles(θ, ϕ)
            @test n isa SVector{3, Float64}
            @test n ≈ [sin(θ) * cos(ϕ), sin(θ) * sin(ϕ), cos(θ)]
            @test norm(n) ≈ 1
            a = angles(Vector(n))
            @test a.θ ≈ θ
            @test a.ϕ ≈ ϕ
            m = @inferred direction_from_angles(ϕ)
            @test m isa SVector{2, Float64}
            @test m ≈ [cos(ϕ), sin(ϕ)]
            @test angles(Vector(m)).θ ≈ ϕ
            # ∂n/∂θ is the unit vector of the meridian.
            @test ForwardDiff.derivative(t -> direction_from_angles(t, ϕ), θ) ≈ [cos(θ) * cos(ϕ), cos(θ) * sin(ϕ), -sin(θ)]
        end
    end
    @testsection "symbols" begin
        @syms θ::real ϕ::real
        n = direction_from_angles(θ, ϕ)
        @test tsimplify(n[1] - sin(θ) * cos(ϕ)) == 0
        @test tsimplify(n[3] - cos(θ)) == 0
        @test tsimplify(sum(x -> x^2, n) - 1) == 0
        @test tsimplify(sum(x -> x^2, direction_from_angles(θ)) - 1) == 0
    end
end
