@testsection "Tensor projections (TI, ORTHO)" begin

    atol_num = 1e-10

    # ── helpers ──────────────────────────────────────────────────────────────────
    n_e3 = [0., 0., 1.]
    n_e1 = [1., 0., 0.]
    frame_canon = CanonicalBasis{3,Float64}()

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensTI{2} — construction & getarray" begin
        A = TensTI{2}(5.0, 8.0, n_e3)
        @test A isa TensTI{2,Float64,2}
        @test A.data == (5.0, 8.0)
        @test A.n == (0.0, 0.0, 1.0)

        M = getarray(A)
        @test size(M) == (3, 3)
        @test M[1,1] ≈ 5.0  atol=atol_num
        @test M[2,2] ≈ 5.0  atol=atol_num
        @test M[3,3] ≈ 8.0  atol=atol_num
        @test M[1,2] ≈ 0.0  atol=atol_num
        @test M[1,3] ≈ 0.0  atol=atol_num
    end

    @testsection "TensTI{2} — traits" begin
        A = TensTI{2}(5.0, 8.0, n_e3)
        @test tr(A) ≈ 18.0  atol=atol_num
        @test issymmetric(A)
        @test isTI(A)
        @test !isISO(A)
        @test !isOrtho(A)

        B = TensTI{2}(5.0, 5.0, n_e3)
        @test isISO(B)
        @test tr(B) ≈ 15.0  atol=atol_num
    end

    @testsection "TensTI{2} — arithmetic" begin
        A = TensTI{2}(5.0, 8.0, n_e3)
        B = TensTI{2}(3.0, 2.0, n_e3)

        C = A + B
        @test C.data == (8.0, 10.0)

        D = A - B
        @test D.data == (2.0, 6.0)

        E = 2.0 * A
        @test E.data == (10.0, 16.0)

        F = A / 2.0
        @test F.data == (2.5, 4.0)

        G = -A
        @test G.data == (-5.0, -8.0)

        H = inv(A)
        @test H.data[1] ≈ 0.2    atol=atol_num
        @test H.data[2] ≈ 0.125  atol=atol_num
    end

    @testsection "TensTI{2} — rotated axis" begin
        # TI tensor with axis along e₁
        A = TensTI{2}(5.0, 8.0, n_e1)
        M = getarray(A)
        @test M[1,1] ≈ 8.0  atol=atol_num   # axial
        @test M[2,2] ≈ 5.0  atol=atol_num   # transverse
        @test M[3,3] ≈ 5.0  atol=atol_num   # transverse

        # General axis
        n45 = [1/√2, 0., 1/√2]
        A45 = TensTI{2}(3.0, 7.0, n45)
        M45 = getarray(A45)
        # diagonal entries:
        # M[1,1] = 3*(1-0.5) + 7*0.5 = 1.5 + 3.5 = 5.0
        @test M45[1,1] ≈ 5.0  atol=atol_num
        @test M45[2,2] ≈ 3.0  atol=atol_num  # n₂=0 → pure transverse
        @test M45[3,3] ≈ 5.0  atol=atol_num
        @test tr(A45) ≈ 13.0  atol=atol_num  # 2*3 + 7
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "proj_tens :TI order 4 — round-trip (n=e₃)" begin
        C = tensTI(10., 3., 2.5, 12., 2., n_e3)
        A = getarray(C)
        B, d, drel = proj_tens(:TI, A, n_e3)
        @test d < atol_num
        @test drel < atol_num
        @test B isa TensWalpole
        # Walpole coefficients should match
        @test collect(argTI(B)) ≈ collect(argTI(C))  atol=atol_num
    end

    @testsection "proj_tens :TI order 4 — isotropic tensor" begin
        𝕀, 𝕁, 𝕂 = ISO(Val(3), Val(Float64))
        k, μ = 10., 5.
        C_iso = 3k * 𝕁 + 2μ * 𝕂
        A = getarray(C_iso)
        # Any axis should give distance ≈ 0 for an isotropic tensor
        B, d, drel = proj_tens(:TI, A, n_e3)
        @test drel < 1e-10
    end

    @testsection "proj_tens :TI order 4 — rotated axis" begin
        n45 = [1/√2, 0., 1/√2]
        C = tensTI(10., 3., 2.5, 12., 2., n45)
        A = getarray(C)
        B, d, drel = proj_tens(:TI, A, n45)
        @test drel < 1e-8
    end

    @testsection "proj_tens :TI order 4 — non-trivial projection" begin
        # Random anisotropic tensor → project to TI → should have TI structure
        C = tensTI(10., 3., 2.5, 12., 2., n_e3)
        A_TI = getarray(C)
        # Add a perturbation that breaks TI symmetry
        A_pert = copy(A_TI)
        A_pert[1,1,1,1] += 1.0
        B, d, drel = proj_tens(:TI, A_pert, n_e3)
        @test d > 0  # non-zero distance (perturbation was applied)
        @test drel < 0.1  # but not too far from TI
        # Re-projecting should give distance 0
        B2, d2, _ = proj_tens(:TI, getarray(B), n_e3)
        @test d2 < atol_num
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "proj_tens :TI order 2 — round-trip (n=e₃)" begin
        A_TI = TensTI{2}(5.0, 8.0, n_e3)
        M = getarray(A_TI)
        B, d, drel = proj_tens(:TI, M, n_e3)
        @test d < atol_num
        @test B isa TensTI{2}
        @test B.data[1] ≈ 5.0  atol=atol_num
        @test B.data[2] ≈ 8.0  atol=atol_num
    end

    @testsection "proj_tens :TI order 2 — rotated axis" begin
        n45 = [1/√2, 0., 1/√2]
        A_TI = TensTI{2}(5.0, 8.0, n45)
        M = getarray(A_TI)
        B, d, drel = proj_tens(:TI, M, n45)
        @test drel < 1e-8
        @test B.data[1] ≈ 5.0  atol=1e-8
        @test B.data[2] ≈ 8.0  atol=1e-8
    end

    @testsection "proj_tens :TI order 2 — isotropic → TI" begin
        M_iso = 7.0 * I(3) |> Matrix{Float64}
        B, d, drel = proj_tens(:TI, M_iso, n_e3)
        @test drel < atol_num
        @test isISO(B)
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "proj_tens :ORTHO order 4 — round-trip" begin
        t = TensOrtho(10., 8., 12., 3., 2.5, 1.5, 2., 3., 3.5, frame_canon)
        A = getarray(t)
        B, d, drel = proj_tens(:ORTHO, A, frame_canon)
        @test drel < atol_num
        @test B isa TensOrtho
        @test collect(getdata(B)) ≈ collect(getdata(t))  atol=atol_num
    end

    @testsection "proj_tens :ORTHO order 4 — TI tensor → ORTHO distance ≈ 0" begin
        C = tensTI(10., 3., 2.5, 12., 2., n_e3)
        A = getarray(C)
        B, d, drel = proj_tens(:ORTHO, A, frame_canon)
        @test drel < 1e-8
    end

    @testsection "proj_tens :ORTHO order 4 — rotated frame" begin
        frame_rot = RotatedBasis(0.3, 0.5, 0.7)
        t = TensOrtho(10., 8., 12., 3., 2.5, 1.5, 2., 3., 3.5, frame_rot)
        A = getarray(t)
        B, d, drel = proj_tens(:ORTHO, A, frame_rot)
        @test drel < 1e-8
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "proj_tens :ORTHO order 2 — round-trip" begin
        M = diagm([5., 8., 12.])
        B, d, drel = proj_tens(:ORTHO, M, frame_canon)
        @test drel < atol_num
        @test B ≈ M  atol=atol_num
    end

    @testsection "proj_tens :ORTHO order 2 — non-diagonal" begin
        M = Float64[5 1 2; 1 8 3; 2 3 12]
        B, d, drel = proj_tens(:ORTHO, M, frame_canon)
        @test B ≈ diagm([5., 8., 12.])  atol=atol_num
        @test d > 0  # off-diagonal terms removed
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "best_sym_tens — fixed basis" begin
        # TI tensor → should detect :TI
        # (skip :ISO in proj due to pre-existing simplify ambiguity in tens_isotropic.jl)
        C_ti = tensTI(10., 3., 2.5, 12., 2., n_e3)
        _, _, _, sym_ti = best_sym_tens(C_ti, n_e3; proj = (:TI, :ORTHO))
        @test sym_ti == :TI

        # ORTHO tensor → should detect :ORTHO
        t_ortho = TensOrtho(10., 8., 12., 3., 2.5, 1.5, 2., 3., 3.5, frame_canon)
        _, _, _, sym_ortho = best_sym_tens(t_ortho, frame_canon; proj = (:TI, :ORTHO))
        @test sym_ortho == :ORTHO
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "Projection helpers — _rot3_raw" begin
        using TensND: _rot3_raw
        R0 = _rot3_raw(0., 0., 0.)
        @test R0 ≈ I(3)  atol=atol_num

        # Third column should be (sinθ cosϕ, sinθ sinϕ, cosθ)
        θ, ϕ = 0.5, 0.8
        R = _rot3_raw(θ, ϕ, 0.)
        # Third column = (sinθ·cosϕ, sinθ·sinϕ, cosθ)
        @test R[1,3] ≈ sin(θ)*cos(ϕ)  atol=atol_num
        @test R[2,3] ≈ sin(θ)*sin(ϕ)  atol=atol_num
        @test R[3,3] ≈ cos(θ)  atol=atol_num
        # R should be orthogonal
        @test R' * R ≈ I(3)  atol=atol_num
    end

    @testsection "Projection helpers — _KM_rotation" begin
        using TensND: _KM_rotation
        Q0 = _KM_rotation(0., 0., 0.)
        @test Q0 ≈ I(6)  atol=atol_num

        # Q should be orthogonal for any angles
        Q = _KM_rotation(0.3, 0.5, 0.7)
        @test Q' * Q ≈ I(6)  atol=atol_num
        @test Q * Q' ≈ I(6)  atol=atol_num
    end

    @testsection "Projection helpers — _project_TI_KM round-trip" begin
        using TensND: _project_TI_KM, _build_TI_KM
        # Build a TI matrix, project it, rebuild → should match
        ℓ = (12., 13., sqrt(2)*2.5, 7., 4.)
        B = _build_TI_KM(ℓ...)
        ℓ2 = _project_TI_KM(B)
        @test collect(ℓ2) ≈ collect(ℓ)  atol=atol_num
        B2 = _build_TI_KM(ℓ2...)
        @test B2 ≈ B  atol=atol_num
    end

    @testsection "Projection helpers — _project_ORTHO_KM round-trip" begin
        using TensND: _project_ORTHO_KM, _build_ORTHO_KM
        params = (10., 8., 12., 3., 2.5, 1.5, 2., 3., 3.5)
        B = _build_ORTHO_KM(params...)
        params2 = _project_ORTHO_KM(B)
        @test collect(params2) ≈ collect(params)  atol=atol_num
    end

end
