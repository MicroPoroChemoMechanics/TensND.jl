@testsection "Tensor projections (TI, ORTHO)" begin

    atol_num = 1.0e-10

    # ── helpers ──────────────────────────────────────────────────────────────────
    n_e3 = [0.0, 0.0, 1.0]
    n_e1 = [1.0, 0.0, 0.0]
    frame_canon = CanonicalBasis{3, Float64}()

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensTI{2} — construction & get_array" begin
        A = TensTI{2}(5.0, 8.0, n_e3)
        @test A isa TensTI{2, Float64, 2}
        @test A.data == (5.0, 8.0)
        @test A.n == (0.0, 0.0, 1.0)

        M = get_array(A)
        @test size(M) == (3, 3)
        @test M[1, 1] ≈ 5.0  atol = atol_num
        @test M[2, 2] ≈ 5.0  atol = atol_num
        @test M[3, 3] ≈ 8.0  atol = atol_num
        @test M[1, 2] ≈ 0.0  atol = atol_num
        @test M[1, 3] ≈ 0.0  atol = atol_num
    end

    @testsection "TensTI{2} — traits" begin
        A = TensTI{2}(5.0, 8.0, n_e3)
        @test tr(A) ≈ 18.0  atol = atol_num
        @test issymmetric(A)
        @test is_TI(A)
        @test !is_ISO(A)
        @test !is_ORTHO(A)

        B = TensTI{2}(5.0, 5.0, n_e3)
        @test is_ISO(B)
        @test tr(B) ≈ 15.0  atol = atol_num
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
        @test H.data[1] ≈ 0.2    atol = atol_num
        @test H.data[2] ≈ 0.125  atol = atol_num
    end

    @testsection "TensTI{2} — rotated axis" begin
        # TI tensor with axis along e₁
        A = TensTI{2}(5.0, 8.0, n_e1)
        M = get_array(A)
        @test M[1, 1] ≈ 8.0  atol = atol_num   # axial
        @test M[2, 2] ≈ 5.0  atol = atol_num   # transverse
        @test M[3, 3] ≈ 5.0  atol = atol_num   # transverse

        # General axis
        n45 = [1 / √2, 0.0, 1 / √2]
        A45 = TensTI{2}(3.0, 7.0, n45)
        M45 = get_array(A45)
        # diagonal entries:
        # M[1,1] = 3*(1-0.5) + 7*0.5 = 1.5 + 3.5 = 5.0
        @test M45[1, 1] ≈ 5.0  atol = atol_num
        @test M45[2, 2] ≈ 3.0  atol = atol_num  # n₂=0 → pure transverse
        @test M45[3, 3] ≈ 5.0  atol = atol_num
        @test tr(A45) ≈ 13.0  atol = atol_num  # 2*3 + 7
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "proj_tens :TI order 4 — round-trip (n=e₃)" begin
        C = tens_TI(10.0, 3.0, 2.5, 12.0, 2.0, n_e3)
        A = get_array(C)
        B, d, drel = proj_tens(:TI, A, n_e3)
        @test d < atol_num
        @test drel < atol_num
        @test B isa TensTI{4}
        # Walpole coefficients should match
        @test collect(arg_TI(B)) ≈ collect(arg_TI(C))  atol = atol_num
    end

    @testsection "proj_tens :TI order 4 — isotropic tensor" begin
        𝕀, 𝕁, 𝕂 = ISO(Val(3), Val(Float64))
        k, μ = 10.0, 5.0
        C_iso = 3k * 𝕁 + 2μ * 𝕂
        A = get_array(C_iso)
        # Any axis should give distance ≈ 0 for an isotropic tensor
        B, d, drel = proj_tens(:TI, A, n_e3)
        @test drel < 1.0e-10
    end

    @testsection "proj_tens :TI order 4 — rotated axis" begin
        n45 = [1 / √2, 0.0, 1 / √2]
        C = tens_TI(10.0, 3.0, 2.5, 12.0, 2.0, n45)
        A = get_array(C)
        B, d, drel = proj_tens(:TI, A, n45)
        @test drel < 1.0e-8
    end

    @testsection "proj_tens :TI order 4 — non-trivial projection" begin
        # Random anisotropic tensor → project to TI → should have TI structure
        C = tens_TI(10.0, 3.0, 2.5, 12.0, 2.0, n_e3)
        A_TI = get_array(C)
        # Add a perturbation that breaks TI symmetry
        A_pert = copy(A_TI)
        A_pert[1, 1, 1, 1] += 1.0
        B, d, drel = proj_tens(:TI, A_pert, n_e3)
        @test d > 0  # non-zero distance (perturbation was applied)
        @test drel < 0.1  # but not too far from TI
        # Re-projecting should give distance 0
        B2, d2, _ = proj_tens(:TI, get_array(B), n_e3)
        @test d2 < atol_num
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "proj_tens :TI order 2 — round-trip (n=e₃)" begin
        A_TI = TensTI{2}(5.0, 8.0, n_e3)
        M = get_array(A_TI)
        B, d, drel = proj_tens(:TI, M, n_e3)
        @test d < atol_num
        @test B isa TensTI{2}
        @test B.data[1] ≈ 5.0  atol = atol_num
        @test B.data[2] ≈ 8.0  atol = atol_num
    end

    @testsection "proj_tens :TI order 2 — rotated axis" begin
        n45 = [1 / √2, 0.0, 1 / √2]
        A_TI = TensTI{2}(5.0, 8.0, n45)
        M = get_array(A_TI)
        B, d, drel = proj_tens(:TI, M, n45)
        @test drel < 1.0e-8
        @test B.data[1] ≈ 5.0  atol = 1.0e-8
        @test B.data[2] ≈ 8.0  atol = 1.0e-8
    end

    @testsection "proj_tens :TI order 2 — isotropic → TI" begin
        M_iso = 7.0 * I(3) |> Matrix{Float64}
        B, d, drel = proj_tens(:TI, M_iso, n_e3)
        @test drel < atol_num
        @test is_ISO(B)
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "proj_tens :ORTHO order 4 — round-trip" begin
        t = TensOrtho(10.0, 8.0, 12.0, 3.0, 2.5, 1.5, 2.0, 3.0, 3.5, frame_canon)
        A = get_array(t)
        B, d, drel = proj_tens(:ORTHO, A, frame_canon)
        @test drel < atol_num
        @test B isa TensOrtho
        @test collect(get_data(B)) ≈ collect(get_data(t))  atol = atol_num
    end

    @testsection "proj_tens :ORTHO order 4 — TI tensor → ORTHO distance ≈ 0" begin
        C = tens_TI(10.0, 3.0, 2.5, 12.0, 2.0, n_e3)
        A = get_array(C)
        B, d, drel = proj_tens(:ORTHO, A, frame_canon)
        @test drel < 1.0e-8
    end

    @testsection "proj_tens :ORTHO order 4 — rotated frame" begin
        frame_rot = RotatedBasis(0.3, 0.5, 0.7)
        t = TensOrtho(10.0, 8.0, 12.0, 3.0, 2.5, 1.5, 2.0, 3.0, 3.5, frame_rot)
        A = get_array(t)
        B, d, drel = proj_tens(:ORTHO, A, frame_rot)
        @test drel < 1.0e-8
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "proj_tens :ORTHO order 2 — round-trip" begin
        M = diagm([5.0, 8.0, 12.0])
        B, d, drel = proj_tens(:ORTHO, M, frame_canon)
        @test drel < atol_num
        @test B ≈ M  atol = atol_num
    end

    @testsection "proj_tens :ORTHO order 2 — non-diagonal" begin
        M = Float64[5 1 2; 1 8 3; 2 3 12]
        B, d, drel = proj_tens(:ORTHO, M, frame_canon)
        @test B ≈ diagm([5.0, 8.0, 12.0])  atol = atol_num
        @test d > 0  # off-diagonal terms removed
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "best_sym_tens — fixed basis" begin
        # TI tensor → should detect :TI
        # (skip :ISO in proj due to pre-existing simplify ambiguity in tens_isotropic.jl)
        C_ti = tens_TI(10.0, 3.0, 2.5, 12.0, 2.0, n_e3)
        _, _, _, sym_ti = best_sym_tens(C_ti, n_e3; proj = (:TI, :ORTHO))
        @test sym_ti == :TI

        # ORTHO tensor → should detect :ORTHO
        t_ortho = TensOrtho(10.0, 8.0, 12.0, 3.0, 2.5, 1.5, 2.0, 3.0, 3.5, frame_canon)
        _, _, _, sym_ortho = best_sym_tens(t_ortho, frame_canon; proj = (:TI, :ORTHO))
        @test sym_ortho == :ORTHO
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "Projection helpers — _rot3_raw" begin
        using TensND: _rot3_raw
        R0 = _rot3_raw(0.0, 0.0, 0.0)
        @test R0 ≈ I(3)  atol = atol_num

        # Third column should be (sinθ cosϕ, sinθ sinϕ, cosθ)
        θ, ϕ = 0.5, 0.8
        R = _rot3_raw(θ, ϕ, 0.0)
        # Third column = (sinθ·cosϕ, sinθ·sinϕ, cosθ)
        @test R[1, 3] ≈ sin(θ) * cos(ϕ)  atol = atol_num
        @test R[2, 3] ≈ sin(θ) * sin(ϕ)  atol = atol_num
        @test R[3, 3] ≈ cos(θ)  atol = atol_num
        # R should be orthogonal
        @test R' * R ≈ I(3)  atol = atol_num
    end

    @testsection "Projection helpers — _KM_rotation" begin
        using TensND: _KM_rotation
        Q0 = _KM_rotation(0.0, 0.0, 0.0)
        @test Q0 ≈ I(6)  atol = atol_num

        # Q should be orthogonal for any angles
        Q = _KM_rotation(0.3, 0.5, 0.7)
        @test Q' * Q ≈ I(6)  atol = atol_num
        @test Q * Q' ≈ I(6)  atol = atol_num
    end

    @testsection "Projection helpers — _project_TI_KM round-trip" begin
        using TensND: _project_TI_KM, _build_TI_KM
        # Build a TI matrix, project it, rebuild → should match
        ℓ = (12.0, 13.0, sqrt(2) * 2.5, 7.0, 4.0)
        B = _build_TI_KM(ℓ...)
        ℓ2 = _project_TI_KM(B)
        @test collect(ℓ2) ≈ collect(ℓ)  atol = atol_num
        B2 = _build_TI_KM(ℓ2...)
        @test B2 ≈ B  atol = atol_num
    end

    @testsection "Projection helpers — _project_ORTHO_KM round-trip" begin
        using TensND: _project_ORTHO_KM, _build_ORTHO_KM
        params = (10.0, 8.0, 12.0, 3.0, 2.5, 1.5, 2.0, 3.0, 3.5)
        B = _build_ORTHO_KM(params...)
        params2 = _project_ORTHO_KM(B)
        @test collect(params2) ≈ collect(params)  atol = atol_num
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "Value-level predicates — is_ISO / is_TI / is_ORTHO" begin
        # ── is_ISO on arrays
        I4 = TensISO{3}(2.0, 3.0)
        @test is_ISO(get_array(I4))
        @test !is_ISO(Float64[10 3 2.5; 3 8 1.5; 2.5 1.5 12])

        # ── is_TI (fixed axis)
        C = tens_TI(10.0, 3.0, 2.5, 12.0, 2.0, n_e3)
        @test is_TI(get_array(C), n_e3)
        @test !is_TI(get_array(C), n_e1)

        # ── is_ORTHO (fixed frame)
        O = TensOrtho(10.0, 8.0, 12.0, 3.0, 2.5, 1.5, 2.0, 3.0, 3.5, frame_canon)
        @test is_ORTHO(get_array(O), frame_canon)
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "best_sym_tens — cheap path (no NLopt)" begin
        # ── ISO detection is always cheap, always succeeds
        Iso = TensISO{3}(2.0, 3.0)
        _, _, _, sym = best_sym_tens(Iso)
        @test sym === :ISO

        # ── TI axis = e₃: candidate axis matches, cheap path succeeds
        C = tens_TI(10.0, 3.0, 2.5, 12.0, 2.0, n_e3)
        _, _, drel, sym = best_sym_tens(C)
        @test sym === :TI
        @test drel < 1.0e-10

        # ── TI axis = (1,1,1)/√3: KM eigendecomposition recovers the axis
        n_tilt = [1.0, 1.0, 1.0] ./ sqrt(3.0)
        Ct = tens_TI(10.0, 3.0, 2.5, 12.0, 2.0, n_tilt)
        _, _, drel_t, sym_t = best_sym_tens(Ct)
        @test sym_t === :TI
        @test drel_t < 1.0e-8

        # ── ORTHO in a non-canonical frame: KM eigendecomposition recovers frame
        ϕ = π / 6
        R = [cos(ϕ) -sin(ϕ) 0.0; sin(ϕ) cos(ϕ) 0.0; 0.0 0.0 1.0]
        rotframe = RotatedBasis(R)
        Orot = TensOrtho(10.0, 8.0, 12.0, 3.0, 2.5, 1.5, 2.0, 3.0, 3.5, rotframe)
        _, _, drel_o, sym_o = best_sym_tens(Orot)
        @test sym_o ∈ (:TI, :ORTHO)
        @test drel_o < 1.0e-8

        # ── Truly anisotropic: no erroneous match
        aniso_KM = Float64[
            10 3 2 1 0 0
            3 8 1 0 1 0
            2 1 12 0 0 1
            1 0 0 4 0 0
            0 1 0 0 5 0
            0 0 1 0 0 6
        ]
        aniso_arr = frommandel(SymmetricTensor{4, 3}, aniso_KM)
        _, _, _, sym_a = best_sym_tens(Tens(Array(aniso_arr)))
        @test sym_a === :ANISO
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "best_sym_tens — fixed-axis/frame path" begin
        # Given axis, TI projection must be exact for a TI tensor
        C = tens_TI(10.0, 3.0, 2.5, 12.0, 2.0, n_e3)
        _, _, drel, sym = best_sym_tens(C, n_e3; proj = (:ISO, :TI))
        @test sym === :TI
        @test drel < 1.0e-12

        # Given frame, ORTHO round-trip on an Ortho tensor
        O = TensOrtho(10.0, 8.0, 12.0, 3.0, 2.5, 1.5, 2.0, 3.0, 3.5, frame_canon)
        _, _, drel_o, sym_o = best_sym_tens(O, frame_canon; proj = (:ISO, :ORTHO))
        @test sym_o === :ORTHO
        @test drel_o < 1.0e-12
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "best_sym_tens — no NLopt needed by default" begin
        # Regression: the default no-argument form must not throw when NLopt
        # is absent from the session.  (NLopt is not loaded in the base test env.)
        aniso_KM = Float64[
            10 3 2 1 0 0
            3 8 1 0 1 0
            2 1 12 0 0 1
            1 0 0 4 0 0
            0 1 0 0 5 0
            0 0 1 0 0 6
        ]
        aniso = Tens(Array(frommandel(SymmetricTensor{4, 3}, aniso_KM)))
        # Must succeed — returns :ANISO without raising.
        res = best_sym_tens(aniso)
        @test res isa Tuple
        @test res[end] === :ANISO
    end

end
