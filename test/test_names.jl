# TensND 0.11 names each operation after what it does: the names of 0.10 it
# replaced, and those 0.8 had deprecated, are gone. The moduli keep their
# symbols as exact aliases.

@testsection "One name per operation: the former names are gone" begin
    former = (
        # renamed in 0.11
        :change_tens, :change_tens_canon, :components_canon, :mandel66, :mandel66_minor, :mandel33_minor,
        :array_from_mandel66, :array_from_mandel33, :ti_args, :ti_args_eng,
        :ti_args_hoenig, :tens_ti, :tens_ti_eng, :tens_ti_hoenig, :tens_cubic, :arg_cubic,
        :get_ℓ, :get_ℓ8, :getcoords, :rot2, :rot3, :rot6, :proj_tens,
        :ti_average_mandel66, :iso_average_mandel66, :ti_average_mat33,
        # deprecated in 0.8
        :KM, :inv_KM, :KM_material, :KM_rotation, :KM_IP, :KM_OP, :ti_params_from_KM,
        :is_ISO, :is_TI, :is_ORTHO, :is_CUBIC, :tens_TI, :arg_TI, :fromISO, :ISO,
        :Walpole, :LeviCivita, :Christoffel, :Lame, :getOM, :GRAD, :SYMGRAD, :DIV,
        :LAPLACE, :HESS, :CoorSystemSym, :CoorSystemNum, :SubManifoldSym,
        :coorsys_spherical_num, :intrinsic, :print_tensor, :Riemann,
    )
    @test all(n -> !isdefined(TensND, n), former)
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

@testsection "Coordinate systems: one type, two geometries" begin
    S = coorsys_spherical()
    P = coorsys_spherical(PointwiseGeometry)
    @test S isa CoordinateSystem{3, <:Any, <:ExactGeometry}
    @test P isa CoordinateSystem{3, Float64, <:PointwiseGeometry}
    @test coorsys_polar(PointwiseGeometry, Float32) isa CoordinateSystem{2, Float32, <:PointwiseGeometry}
end
