# The names of TensND 0.7, renamed in 0.8 to Julia's convention: each former
# name still works, warns naming its successor, and returns what the successor
# returns; the moduli keep their symbols as exact aliases.

@testsection "Former names: deprecated, same results" begin
    C = TensISO{3}(3 * 20.0, 2 * 12.0)
    ℬ = RotatedBasis(0.3, 0.5, 0.2)
    O = TensOrtho(20.0, 25.0, 30.0, -4.0, 1.0, 3.0, 6.0, 7.0, 8.0, ℬ)
    n = (0.0, 0.0, 1.0)
    T = tens_ti(10.0, 3.0, 2.5, 12.0, 2.0, n)
    M = mandel(C)
    CS = coorsys_spherical()

    # (former call, call in 0.8), one per renamed function
    pairs = [
        (() -> KM(O), () -> mandel(O)),
        (() -> inv_KM(M), () -> from_mandel(M)),
        (() -> KM_material(O), () -> mandel_material(O)),
        (() -> KM_rotation(0.3, 0.5, 0.2), () -> mandel_rotation(0.3, 0.5, 0.2)),
        (() -> ti_params_from_KM(M), () -> ti_params_from_mandel(M)),
        (() -> KM_from_ti_params(ti_params_from_mandel(M)...), () -> mandel_from_ti_params(ti_params_from_mandel(M)...)),
        (() -> ti8_params_from_KM(M), () -> ti8_params_from_mandel(M)),
        (() -> KM_from_ti8_params(ti8_params_from_mandel(M)), () -> mandel_from_ti8_params(ti8_params_from_mandel(M))),
        (() -> ortho_params_from_KM(M), () -> ortho_params_from_mandel(M)),
        (() -> KM_from_ortho_params(ortho_params_from_mandel(M)...), () -> mandel_from_ortho_params(ortho_params_from_mandel(M)...)),
        (() -> cubic_params_from_KM(M), () -> cubic_params_from_mandel(M)),
        (() -> KM_from_cubic_params(cubic_params_from_mandel(M)...), () -> mandel_from_cubic_params(cubic_params_from_mandel(M)...)),
        (() -> km_ip_block(M), () -> mandel_ip_block(M)),
        (() -> km_op_block(M), () -> mandel_op_block(M)),
        (() -> km_ip_embed(mandel_ip_block(M)), () -> mandel_ip_embed(mandel_ip_block(M))),
        (() -> km_op_embed(mandel_op_block(M)), () -> mandel_op_embed(mandel_op_block(M))),
        (() -> km_inv_schur(M), () -> mandel_inv_schur(M)),
        (() -> is_ISO(C), () -> isisotropic(C)),
        (() -> is_TI(T), () -> istransverselyisotropic(T)),
        (() -> is_ORTHO(O), () -> isorthotropic(O)),
        (() -> is_CUBIC(C), () -> iscubic(C)),
        (() -> tens_TI(10.0, 3.0, 2.5, 12.0, 2.0, n), () -> tens_ti(10.0, 3.0, 2.5, 12.0, 2.0, n)),
        (() -> arg_TI(T), () -> ti_args(T)),
        (() -> arg_TI_eng(inv(T)), () -> ti_args_eng(inv(T))),
        (() -> arg_TI_Hoenig(inv(T)), () -> ti_args_hoenig(inv(T))),
        (() -> fromISO(C, n), () -> iso_to_ti(C, n)),
        (() -> ISO(Val(3), Val(Float64)), () -> iso_projectors(Val(3), Val(Float64))),
        (() -> Walpole(n), () -> walpole_basis(n)),
        (() -> Walpole(n; sym = true), () -> walpole_basis_sym(n)),
        (() -> LeviCivita(Float64), () -> levi_civita(Float64)),
        (() -> Christoffel(CS), () -> christoffel(CS)),
        (() -> Lame(CS), () -> lame_coefficients(CS)),
        (() -> getOM(CS), () -> position_vector(CS)),
    ]
    for (old, new) in pairs
        @test isequal(@test_deprecated(old()), new())
    end
    let E₁ = 10.0, E₃ = 3.0
        S = tens_ti_eng(E₁, E₃, 0.2, 0.3, 2.0, n)
        @test isequal(@test_deprecated(tens_TI_eng(E₁, E₃, 0.2, 0.3, 2.0, n)), S)
        @test isequal(@test_deprecated(tens_TI_Hoenig(E₁, 0.2, 0.3, 1.0, 2.0, n)), tens_ti_hoenig(E₁, 0.2, 0.3, 1.0, 2.0, n))
    end

    # The constants, and the differential operators on a scalar field.
    # A deprecated binding warns where the name is resolved, not when it runs.
    @test Base.isdeprecated(TensND, :KM_IP) && TensND.KM_IP === MANDEL_IP
    @test Base.isdeprecated(TensND, :KM_OP) && TensND.KM_OP === MANDEL_OP
    r, θ, ϕ = getcoords(CS)
    f = r^2 * cos(θ)
    for (old, new) in ((GRAD, grad), (LAPLACE, laplacian), (HESS, hess))
        @test isequal(@test_deprecated(old(f, CS)), new(f, CS))
    end
    v = grad(f, CS)
    @test isequal(@test_deprecated(SYMGRAD(v, CS)), symgrad(v, CS))
    @test isequal(@test_deprecated(DIV(v, CS)), divergence(v, CS))
    # `divergence` is Tensors.jl's function, extended: one name for both.
    @test divergence === TensND.Tensors.divergence
end

@testsection "Moduli: the symbols are aliases of the words" begin
    @test E_nu === young_poisson
    @test k_mu === bulk_shear
    @test mu_nu === shear_poisson
    @test iso_stiffness_E_nu === iso_stiffness_young_poisson
    C = iso_stiffness_young_poisson(200.0, 0.3)
    E, ν = young_poisson(C)
    @test E ≈ 200.0 && ν ≈ 0.3
    @test all(bulk_shear(C) .≈ (200.0 / (3 * (1 - 0.6)), 200.0 / 2.6))
end
