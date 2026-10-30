# The isotropic elastic moduli: the stiffness built from each of the fifteen
# pairs of (k, E, λ, μ, ν, M) and every modulus read back, in 3-D and in plane
# strain, against the closed forms written from (λ, μ) alone; the physical
# roots of the two quadratic conversions over Poisson's ratios down to −1;
# derivatives; symbolic moduli; the aliases.

using TensND, Test, LinearAlgebra, ForwardDiff
import Symbolics
using SymPyPythonCall: Sym, symbols, simplify

# The six moduli of the 3-D material (λ, μ), with `k` the bulk modulus of the
# dimension (in 2-D, the in-plane λ + μ), and its stiffness written directly.
_im_moduli(::Val{3}, λ, μ) = (bulk = λ + 2μ / 3, young = μ * (3λ + 2μ) / (λ + μ), lame = λ, shear = μ, poisson = λ / (2(λ + μ)), pwave = λ + 2μ)
_im_moduli(::Val{2}, λ, μ) = merge(_im_moduli(Val(3), λ, μ), (bulk = λ + μ,))
_im_stiffness(::Val{dim}, λ, μ) where {dim} = TensISO{dim}(dim * λ + 2μ, 2μ)

const _IM_NAMES = (:bulk, :young, :lame, :shear, :poisson, :pwave)
const _IM_READERS = (bulk_modulus, young_modulus, lame_lambda, shear_modulus, poisson_ratio, pwave_modulus)
_im_builder(p, q) = getfield(TensND, (p, q) == (:bulk, :shear) ? :iso_stiffness : Symbol(:iso_stiffness_, p, :_, q))
_im_pair_reader(p, q) = getfield(TensND, Symbol(p, :_, q))

@testset "isotropic moduli: every pair, against the closed forms" begin
    for v in (Val(3), Val(2)), ν in (-0.9, -0.5, -0.2, 0.0, 0.1, 0.25, 0.3, 0.45, 0.49)
        μ = 7.0
        λ = 2μ * ν / (1 - 2ν)
        ref = _im_moduli(v, λ, μ)
        C = _im_stiffness(v, λ, μ)
        for (n, r) in zip(_IM_NAMES, _IM_READERS)
            @test r(C) ≈ ref[n] atol = 1.0e-12 * μ
        end
        for i in 1:6, j in (i + 1):6
            p, q = _IM_NAMES[i], _IM_NAMES[j]
            @test all(_im_pair_reader(p, q)(C) .≈ (ref[p], ref[q]))
            # The 2-D (k, E) conversion takes the root ν ≥ −1/4; the other root
            # is the material of Poisson's ratio −1/2 − ν, with the same k and E.
            (v === Val(2) && (p, q) == (:bulk, :young) && ν < -0.25) && continue
            # (E, M) is shared by a material of ν ≥ 0 and one of ν ≤ 0: the first is kept.
            ((p, q) == (:young, :pwave) && ν < 0) && continue
            # λ = 0 exactly when ν = 0: the pair (λ, ν) then says nothing of μ.
            ((p, q) == (:lame, :poisson) && iszero(ν)) && continue
            B = _im_builder(p, q)(v, ref[p], ref[q])
            @test B isa TensISO{4, v === Val(3) ? 3 : 2}
            @test all(isapprox.(get_data(B), get_data(C); atol = 1.0e-10 * μ))
            v === Val(3) && @test get_data(_im_builder(p, q)(ref[p], ref[q])) == get_data(B)
        end
    end
    # The other root of the 2-D (k, E) conversion.
    μ, ν = 7.0, -0.4
    λ = 2μ * ν / (1 - 2ν)
    ref = _im_moduli(Val(2), λ, μ)
    @test poisson_ratio(iso_stiffness_bulk_young(Val(2), ref.bulk, ref.young)) ≈ -0.5 - ν
    # The other material of the same (E, M), of ν ≤ 0.
    μ, ν = 7.0, -0.3
    λ = 2μ * ν / (1 - 2ν)
    ref = _im_moduli(Val(3), λ, μ)
    B = iso_stiffness_young_pwave(ref.young, ref.pwave)
    @test poisson_ratio(B) > 0 && young_modulus(B) ≈ ref.young && pwave_modulus(B) ≈ ref.pwave
    # The builders of 0.10 keep their results.
    @test get_data(iso_stiffness(5.0, 2.0)) == (15.0, 4.0)
    @test get_data(iso_stiffness_young_poisson(200.0, 0.3)) == get_data(iso_stiffness(200.0 / (3 * 0.4), 200.0 / 2.6))
end

@testset "isotropic moduli: derivatives and inference" begin
    λ, μ = 3.0, 2.0
    g = ForwardDiff.gradient(x -> young_modulus(iso_stiffness_lame_shear(x[1], x[2])), [λ, μ])
    @test g ≈ [μ^2, 3λ^2 + 4λ * μ + 2μ^2] ./ (λ + μ)^2
    # Through a quadratic conversion: E read back from (E, λ), λ unchanged.
    E = young_modulus(iso_stiffness_lame_shear(λ, μ))
    @test ForwardDiff.derivative(e -> young_modulus(iso_stiffness_young_lame(e, λ)), E) ≈ 1
    @test abs(ForwardDiff.derivative(e -> lame_lambda(iso_stiffness_young_lame(e, λ)), E)) < 1.0e-12
    M = λ + 2μ
    @test ForwardDiff.derivative(m -> pwave_modulus(iso_stiffness_young_pwave(E, m)), M) ≈ 1
    @test ForwardDiff.derivative(k -> bulk_modulus(iso_stiffness_bulk_young(Val(2), k, E)), λ + μ) ≈ 1
    # Nested, under two tags.
    d2 = ForwardDiff.derivative(a -> ForwardDiff.derivative(b -> shear_modulus(iso_stiffness_young_poisson(a, b)), 0.3), 200.0)
    @test d2 ≈ -1 / (2 * 1.3^2)
    @test (@inferred iso_stiffness_young_lame(E, λ)) isa TensISO{4, 3, Float64}
    @test (@inferred iso_stiffness_bulk_poisson(Val(2), 5.0, 0.3)) isa TensISO{4, 2, Float64}
    @test (@inferred young_shear(iso_stiffness(5.0, 2.0))) isa NTuple{2, Float64}
end

@testset "isotropic moduli: symbolic" begin
    Symbolics.@variables λ μ
    vals = Dict(λ => 3.0, μ => 2.0)
    num(x) = Symbolics.value(Symbolics.substitute(x, vals; fold = Val(true)))
    ref = _im_moduli(Val(3), 3.0, 2.0)
    C = iso_stiffness_lame_shear(λ, μ)
    for (n, r) in zip(_IM_NAMES, _IM_READERS)
        @test num(r(C)) ≈ ref[n]
    end
    # A stiffness from symbolic (E, ν), read back.
    Symbolics.@variables E ν
    B = iso_stiffness_young_poisson(Val(2), E, ν)
    @test num(Symbolics.substitute(lame_lambda(B), Dict(E => ref.young, ν => ref.poisson))) ≈ 3.0
    # SymPy: the closed forms, exactly.
    l, m = symbols("λ μ", positive = true)
    S = iso_stiffness_lame_shear(l, m)
    @test simplify(young_modulus(S) - m * (3l + 2m) / (l + m)) == 0
    @test simplify(poisson_ratio(S) - l / (2(l + m))) == 0
    @test simplify(lame_lambda(iso_stiffness_shear_poisson(m, poisson_ratio(S))) - l) == 0
end

@testset "isotropic moduli: aliases" begin
    @test k_mu === k_μ === bulk_shear
    @test E_nu === E_ν === young_poisson
    @test mu_nu === μ_ν === shear_poisson
    @test lambda_mu === λ_μ === lame_shear
    @test E_lambda === E_λ === young_lame
    @test iso_stiffness_E_nu === iso_stiffness_E_ν === iso_stiffness_young_poisson
    @test iso_stiffness_k_mu === iso_stiffness_k_μ === iso_stiffness_bulk_shear === iso_stiffness
    @test iso_stiffness_λ_ν === iso_stiffness_lambda_nu === iso_stiffness_lame_poisson
    # No short alias for the P-wave modulus, nor mixed ASCII and Greek ones.
    @test !any(n -> isdefined(TensND, n), (:k_M, :E_M, :mu_M, :μ_M, :nu_M, :lambda_μ, :λ_mu))
end

@testset "plane strain and plane stress" begin
    E, ν = 200.0, 0.3
    C = iso_stiffness_young_poisson(E, ν)
    λ, μ = lame_shear(C)
    # Plane strain: the in-plane block of the 3-D stiffness.
    P = plane_strain_stiffness(C)
    @test get_array(P) ≈ get_array(C)[1:2, 1:2, 1:2, 1:2]
    @test all(young_poisson(P) .≈ (E, ν))
    # Plane stress: ε₃₃ condensed out of the 3-D Hooke law, σ₃₃ = 0.
    S = plane_stress_stiffness(C)
    A = get_array(C)
    cond = [A[i, j, k, l] - A[i, j, 3, 3] * A[3, 3, k, l] / A[3, 3, 3, 3] for i in 1:2, j in 1:2, k in 1:2, l in 1:2]
    @test get_array(S) ≈ cond
    @test lame_lambda(S) ≈ E * ν / (1 - ν^2)
    # The equivalent moduli, both ways.
    @test all(young_poisson(S) .≈ plane_stress_to_plane_strain(E, ν))
    @test all(get_data(iso_stiffness_young_poisson(Val(2), plane_stress_to_plane_strain(E, ν)...)) .≈ get_data(S))
    E′, ν′ = plane_strain_to_plane_stress(E, ν)
    @test all(get_data(plane_stress_stiffness(iso_stiffness_young_poisson(E′, ν′))) .≈ get_data(P))
    @test all(plane_stress_to_plane_strain(plane_strain_to_plane_stress(E, ν)...) .≈ (E, ν))
    # Symbolic, exactly.
    e, n = symbols("E ν", positive = true)
    Ss = plane_stress_stiffness(iso_stiffness_young_poisson(e, n))
    @test simplify(lame_lambda(Ss) - e * n / (1 - n^2)) == 0
    @test ForwardDiff.derivative(x -> first(plane_strain_to_plane_stress(x, ν)), E) ≈ 1 / (1 - ν^2)
end
