@testsection "Walpole & Ortho tensors" begin

    # ─── helpers shared across sub-sections ───────────────────────────────────
    n3 = 𝐞(Val(3), Val(3), Val(Float64))   # e₃ as Float64 Vec
    n3s = 𝐞(Val(3), Val(3), Val(Sym))       # e₃ as Sym Vec
    atol_num = 1.0e-12

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensTI{4} — construction & traits" begin
        W1, W2, W3, W4, W5, W6 = Walpole(n3)
        @test W1 isa TensTI{4, Float64, 6}
        @test W2 isa TensTI{4, Float64, 6}
        @test size(W1) == (3, 3, 3, 3)
        @test get_basis(W1) isa CanonicalBasis{3, Float64}
        @test get_var(W1) == (:cont, :cont, :cont, :cont)

        # N=5 (symmetric) basis
        W1s, W2s, W3s, W4s, W5s = Walpole(n3; sym = true)
        @test W1s isa TensTI{4, Float64, 5}
        @test size(W1s) == (3, 3, 3, 3)

        # individual constructors
        @test tens_W1(n3) isa TensTI{4, Float64, 6}
        @test tens_W3(n3) isa TensTI{4, Float64, 6}

        # axis accessor
        @test axis(W1) == (0.0, 0.0, 1.0)
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensTI{4} — get_array & KM (n=e₃, Float64)" begin
        W1, W2, W3, W4, W5, W6 = Walpole(n3)
        sq2 = sqrt(2.0)

        # W₁ = nₙ⊗nₙ  →  only entry [3,3,3,3] = 1
        A1 = get_array(W1)
        @test A1[3, 3, 3, 3] ≈ 1.0   atol = atol_num
        @test all(
            abs.(A1[i, j, k, l]) < atol_num
                for i in 1:3, j in 1:3, k in 1:3, l in 1:3
                if !(i == 3 && j == 3 && k == 3 && l == 3)
        )

        # W₂ = (nT⊗nT)/2 with n=e₃ → nT = diag(1,1,0)
        # so W₂[i,j,k,l] = δᵢⱼ(1-δᵢ₃)δₖₗ(1-δₖ₃)/2
        A2 = get_array(W2)
        @test A2[1, 1, 1, 1] ≈ 0.5   atol = atol_num
        @test A2[1, 1, 2, 2] ≈ 0.5   atol = atol_num
        @test A2[2, 2, 2, 2] ≈ 0.5   atol = atol_num
        @test A2[3, 3, 3, 3] ≈ 0.0   atol = atol_num

        # W₃[3,3,1,1] = 1/√2
        A3 = get_array(W3)
        @test A3[3, 3, 1, 1] ≈ 1.0 / sq2   atol = atol_num
        @test A3[3, 3, 2, 2] ≈ 1.0 / sq2   atol = atol_num
        @test A3[3, 3, 3, 3] ≈ 0.0       atol = atol_num

        # W₄[1,1,3,3] = 1/√2
        A4 = get_array(W4)
        @test A4[1, 1, 3, 3] ≈ 1.0 / sq2   atol = atol_num

        # W₅: shear in transverse plane (e₁,e₂)
        # W₅[1,2,1,2] = W₅[2,1,1,2] = W₅[1,2,2,1] = W₅[2,1,2,1] = 1/2
        A5 = get_array(W5)
        @test A5[1, 2, 1, 2] ≈ 0.5   atol = atol_num
        @test A5[1, 1, 1, 1] ≈ 0.5   atol = atol_num
        @test A5[3, 3, 3, 3] ≈ 0.0   atol = atol_num

        # W₆: shear between transverse and axial
        A6 = get_array(W6)
        @test A6[1, 3, 1, 3] ≈ 0.5   atol = atol_num
        @test A6[2, 3, 2, 3] ≈ 0.5   atol = atol_num
        @test A6[1, 1, 1, 1] ≈ 0.0   atol = atol_num
        @test A6[3, 3, 3, 3] ≈ 0.0   atol = atol_num

        # KM structure for n=e₃: blocks should separate
        L = TensTI{4}(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, n3)
        Km = KM(L)
        @test size(Km) == (6, 6)
        # Off-diagonal shear coupling should be zero for axis n=e₃
        @test abs(Km[1, 4]) < atol_num
        @test abs(Km[1, 5]) < atol_num
        @test abs(Km[4, 6]) < atol_num
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensTI{4} — Walpole product rule (Float64)" begin
        W1, W2, W3, W4, W5, W6 = Walpole(n3)

        # Idempotents: W₁⊡W₁=W₁, W₂⊡W₂=W₂, W₅⊡W₅=W₅, W₆⊡W₆=W₆
        @test opequal(get_array(W1 ⊡ W1), get_array(W1))
        @test opequal(get_array(W2 ⊡ W2), get_array(W2))
        @test opequal(get_array(W5 ⊡ W5), get_array(W5))
        @test opequal(get_array(W6 ⊡ W6), get_array(W6))

        # Cross products
        @test opequal(get_array(W3 ⊡ W4), get_array(W1))
        @test opequal(get_array(W4 ⊡ W3), get_array(W2))

        # Zero cross products between incompatible blocks
        zero4 = zeros(3, 3, 3, 3)
        @test opequal(get_array(W1 ⊡ W2), zero4)
        @test opequal(get_array(W1 ⊡ W5), zero4)
        @test opequal(get_array(W5 ⊡ W6), zero4)
        @test opequal(get_array(W6 ⊡ W5), zero4)

        # General product: Walpole vs direct array contraction
        L = TensTI{4}(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, n3)
        M = TensTI{4}(0.5, 1.5, 2.0, 0.3, 0.8, 1.2, n3)
        LM_walpole = get_array(L ⊡ M)
        LM_direct = Tensor{4, 3}(get_array(L)) ⊡ Tensor{4, 3}(get_array(M))
        @test opequal(LM_walpole, Array(LM_direct))
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensTI{4} — inverse (Float64)" begin
        𝕀 = tens_Id4(Val(3), Val(Float64))
        n = [1.0 / sqrt(3.0), 1.0 / sqrt(3.0), 1.0 / sqrt(3.0)]

        # N=5 (symmetric)
        L5 = TensTI{4}(2.0, 3.0, 1.0, 4.0, 5.0, n)
        Li5 = inv(L5)
        @test Li5 isa TensTI{4, Float64, 5}
        LLi5 = get_array(L5 ⊡ Li5)
        𝕀arr = get_array(𝕀)
        @test opequal(LLi5, 𝕀arr)

        # N=6 (general)
        L6 = TensTI{4}(2.0, 3.0, 1.0, 0.5, 4.0, 5.0, n)
        Li6 = inv(L6)
        @test Li6 isa TensTI{4, Float64, 6}
        LLi6 = get_array(L6 ⊡ Li6)
        @test opequal(LLi6, 𝕀arr)
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensTI{4} — fromISO" begin
        𝕀, 𝕁, 𝕂 = ISO(Val(3), Val(Float64))
        α, β = 3.0, 2.0
        ℂiso = α * 𝕁 + β * 𝕂

        # For any axis, fromISO should give the same array as the isotropic tensor
        for n in (
                [0.0, 0.0, 1.0], [1.0, 0.0, 0.0],
                [1.0 / sqrt(3), 1.0 / sqrt(3), 1.0 / sqrt(3)],
            )
            Wiso = fromISO(ℂiso, n)
            @test Wiso isa TensTI{4, Float64, 5}
            @test is_TI(Wiso)
            @test opequal(get_array(Wiso), get_array(ℂiso))
        end

        # Symbolic — use Sym ISO tensors to avoid Float64/Sym residuals
        αs, βs = symbols("α β", real = true)
        𝕁s, 𝕂s = tens_J4(Val(3), Val(Sym)), tens_K4(Val(3), Val(Sym))
        ℂisos = αs * 𝕁s + βs * 𝕂s
        Wisos = fromISO(ℂisos, n3s)
        @test Wisos isa TensTI{4, <:Any, 5}
        for i in 1:3, j in 1:3, k in 1:3, l in 1:3
            @test simplify(Wisos[i, j, k, l] - ℂisos[i, j, k, l]) == 0
        end
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensTI{4} — is_ISO / is_TI / is_ORTHO" begin
        W = tens_W1(n3)
        @test !is_ISO(W)
        @test  is_TI(W)
        @test !is_ORTHO(W)
        𝕀, 𝕁, 𝕂 = ISO(Val(3), Val(Float64))
        @test !is_TI(𝕁)
        @test !is_ORTHO(𝕁)
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensOrtho — is_ISO / is_TI / is_ORTHO" begin
        frame3 = CanonicalBasis{3, Float64}()
        t = TensOrtho(10.0, 8.0, 9.0, 3.0, 2.0, 4.0, 2.5, 3.0, 1.5, frame3)
        @test !is_ISO(t)
        @test !is_TI(t)
        @test  is_ORTHO(t)
        # universal fallback
        @test !is_ORTHO(42)
        @test !is_ORTHO("string")
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensTI{4} — show" begin
        L5 = TensTI{4}(1.0, 2.0, 0.5, 3.0, 4.0, n3)   # N=5
        L6 = TensTI{4}(1.0, 2.0, 0.5, 0.3, 3.0, 4.0, n3)   # N=6
        buf5 = IOBuffer()
        show(buf5, L5)
        s5 = String(take!(buf5))
        @test contains(s5, "W") && contains(s5, "axis")

        buf6 = IOBuffer()
        show(buf6, L6)
        s6 = String(take!(buf6))
        @test contains(s6, "W") && contains(s6, "axis")
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensOrtho — show" begin
        frame3 = CanonicalBasis{3, Float64}()
        t = TensOrtho(10.0, 8.0, 9.0, 3.0, 2.0, 4.0, 2.5, 3.0, 1.5, frame3)
        buf = IOBuffer()
        show(buf, t)
        s = String(take!(buf))
        @test contains(s, "P₁⊗P₁") && contains(s, "frame")
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensTI{4} — tsimplify (symbolic, _rebuild)" begin
        ℓ₁, ℓ₂, ℓ₃ = symbols("ℓ₁ ℓ₂ ℓ₃", real = true)
        L = TensTI{4}(ℓ₁, ℓ₂, ℓ₃, ℓ₁ + ℓ₂, ℓ₂ + ℓ₃, n3s)   # N=6
        Ls = tsimplify(L)
        @test Ls isa TensTI{4}   # _rebuild preserves type
        @test get_ℓ(Ls)[1] == ℓ₁   # simplification is a no-op here
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensOrtho — tsimplify (_rebuild)" begin
        # TensOrtho frames are always numeric; test _rebuild with Float64 data.
        # tsimplify on a non-symbolic NTuple is a no-op, but _rebuild must still
        # reconstruct the TensOrtho, verifying the dispatch path.
        frame3 = CanonicalBasis{3, Float64}()
        t = TensOrtho(10.0, 8.0, 9.0, 3.0, 2.0, 4.0, 2.5, 3.0, 1.5, frame3)
        ts = tsimplify(t)
        @test ts isa TensOrtho   # _rebuild preserves type
        @test get_data(ts)[1] ≈ 10.0
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensTI{4} — arithmetic" begin
        W1, W2, W3, W4, W5, W6 = Walpole(n3)
        L = TensTI{4}(1.0, 2.0, 0.5, 0.5, 3.0, 4.0, n3)   # N=6
        M = TensTI{4}(0.5, 1.0, 0.25, 0.25, 1.5, 2.0, n3)

        @test opequal(get_array(L + M), get_array(L) .+ get_array(M))
        @test opequal(get_array(L - M), get_array(L) .- get_array(M))
        @test opequal(get_array(2.0 * L), 2.0 .* get_array(L))
        @test opequal(get_array(-L), .-get_array(L))

        # Symmetric N=5
        Ls = TensTI{4}(1.0, 2.0, 0.5, 3.0, 4.0, n3)
        Ms = TensTI{4}(0.5, 1.0, 0.25, 1.5, 2.0, n3)
        @test (Ls + Ms) isa TensTI{4, Float64, 5}
        @test opequal(get_array(Ls + Ms), get_array(Ls) .+ get_array(Ms))
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensTI{4} — symbolic inverse" begin
        ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆ = symbols("ℓ₁ ℓ₂ ℓ₃ ℓ₅ ℓ₆", real = true)
        L = TensTI{4}(ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆, n3s)
        Li = inv(L)
        @test Li isa TensTI{4, <:Any, 5}
        # L⊡inv(L) should be identity
        𝕀 = tens_Id4(Val(3), Val(Sym))
        prod = L ⊡ Li
        for i in 1:3, j in 1:3, k in 1:3, l in 1:3
            @test simplify(prod[i, j, k, l] - 𝕀[i, j, k, l]) == 0
        end
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensOrtho — construction & traits" begin
        frame3 = CanonicalBasis{3, Float64}()
        C11, C22, C33 = 10.0, 8.0, 12.0
        C12, C13, C23 = 3.0, 4.0, 2.5
        C44, C55, C66 = 2.0, 3.0, 1.5
        t = TensOrtho(C11, C22, C33, C12, C13, C23, C44, C55, C66, frame3)
        @test t isa TensOrtho{Float64}
        @test size(t) == (3, 3, 3, 3)
        @test get_basis(t) isa CanonicalBasis{3, Float64}
        @test get_var(t) == (:cont, :cont, :cont, :cont)
        @test frame(t) === frame3
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensOrtho — KM_material (canonical frame)" begin
        frame3 = CanonicalBasis{3, Float64}()
        C11, C22, C33 = 10.0, 8.0, 12.0
        C12, C13, C23 = 3.0, 4.0, 2.5
        C44, C55, C66 = 2.0, 3.0, 1.5
        t = TensOrtho(C11, C22, C33, C12, C13, C23, C44, C55, C66, frame3)
        Km = KM_material(t)
        @test size(Km) == (6, 6)
        # Diagonal blocks
        @test Km[1, 1] ≈ C11  atol = atol_num
        @test Km[2, 2] ≈ C22  atol = atol_num
        @test Km[3, 3] ≈ C33  atol = atol_num
        @test Km[4, 4] ≈ 2 * C44 atol = atol_num
        @test Km[5, 5] ≈ 2 * C55 atol = atol_num
        @test Km[6, 6] ≈ 2 * C66 atol = atol_num
        # Off-diagonal within normal block
        @test Km[1, 2] ≈ C12  atol = atol_num
        @test Km[1, 3] ≈ C13  atol = atol_num
        @test Km[2, 3] ≈ C23  atol = atol_num
        # Zeros between normal and shear blocks
        @test Km[1, 4] ≈ 0.0  atol = atol_num
        @test Km[2, 5] ≈ 0.0  atol = atol_num
        @test Km[3, 6] ≈ 0.0  atol = atol_num
        @test Km[4, 5] ≈ 0.0  atol = atol_num
        @test Km[4, 6] ≈ 0.0  atol = atol_num
        @test Km[5, 6] ≈ 0.0  atol = atol_num
        # Symmetry
        @test Km ≈ Km'  atol = atol_num
        # KM in canonical frame should match
        Km2 = KM(t)
        @test Km2 ≈ Km  atol = atol_num
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensOrtho — from KM matrix" begin
        frame3 = CanonicalBasis{3, Float64}()
        C11, C22, C33 = 10.0, 8.0, 12.0
        C12, C13, C23 = 3.0, 4.0, 2.5
        C44, C55, C66 = 2.0, 3.0, 1.5
        t1 = TensOrtho(C11, C22, C33, C12, C13, C23, C44, C55, C66, frame3)
        Km = KM_material(t1)
        t2 = TensOrtho(Km, frame3)
        @test opequal(get_array(t1), get_array(t2))
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensOrtho — inverse (canonical frame)" begin
        frame3 = CanonicalBasis{3, Float64}()
        C11, C22, C33 = 10.0, 8.0, 12.0
        C12, C13, C23 = 3.0, 4.0, 2.5
        C44, C55, C66 = 2.0, 3.0, 1.5
        t = TensOrtho(C11, C22, C33, C12, C13, C23, C44, C55, C66, frame3)
        ti = inv(t)
        @test ti isa TensOrtho{Float64}
        𝕀 = tens_Id4(Val(3), Val(Float64))
        A = Tensor{4, 3}(get_array(t))
        Ai = Tensor{4, 3}(get_array(ti))
        prod = Array(A ⊡ Ai)
        𝕀arr = get_array(𝕀)
        @test opequal(prod, 𝕀arr)
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensOrtho — rotated frame" begin
        # Rotate material frame: e₁→e₂, e₂→e₃, e₃→e₁ (cyclic permutation)
        R = Float64[0 0 1; 1 0 0; 0 1 0]
        frame_rot = RotatedBasis(R)
        C11, C22, C33 = 10.0, 10.0, 12.0
        C12, C13, C23 = 3.0, 3.0, 3.0
        C44, C55, C66 = 2.0, 2.0, 2.0
        t = TensOrtho(C11, C22, C33, C12, C13, C23, C44, C55, C66, frame_rot)
        # KM_material should still have block-diagonal structure
        Km = KM_material(t)
        @test Km[1, 4] ≈ 0.0  atol = atol_num
        @test Km[4, 5] ≈ 0.0  atol = atol_num
        @test Km[5, 6] ≈ 0.0  atol = atol_num
        # change_tens to canonical basis should give a valid Tens
        tc = change_tens(t, CanonicalBasis{3, Float64}())
        @test tc isa AbstractTens
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensOrtho — TI consistency" begin
        # A TI tensor (C11=C22, C13=C23, C44=C55, C66=(C11-C12)/2)
        # should give same array as corresponding TensTI{4} (with n=e₃)
        frame3 = CanonicalBasis{3, Float64}()
        C11 = 10.0; C33 = 12.0; C12 = 3.0; C13 = 2.5; C44 = 2.0
        C22 = C11; C23 = C13; C55 = C44; C66 = (C11 - C12) / 2
        to = TensOrtho(C11, C22, C33, C12, C13, C23, C44, C55, C66, frame3)

        # Build equivalent TensTI{4} (n=e₃) via the engineering constants.
        # For a TI material with n=e₃:
        #   C₁₁=C₂₂=(ℓ₂+ℓ₅)/2, C₃₃=ℓ₁, C₁₂=(ℓ₂-ℓ₅)/2, C₁₃=C₂₃=ℓ₃/√2, C₄₄=ℓ₆/2, C₆₆=ℓ₅/2
        # Inverting: ℓ₂=C₁₁+C₁₂, ℓ₅=C₁₁−C₁₂=2C₆₆, ℓ₃=C₁₃√2, ℓ₆=2C₄₄
        sq2 = sqrt(2.0)
        ℓ₁ = C33
        ℓ₅ = C66 * 2       # = C11 - C12
        ℓ₂ = C11 + C12     # NOT (C11+C12)/2
        ℓ₃ = C13 * sq2
        ℓ₆ = C44 * 2
        tw = TensTI{4}(ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆, n3)
        @test opequal(get_array(to), get_array(tw))
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensTI{4} — dcontract with TensISO" begin
        𝕀, 𝕁, 𝕂 = ISO(Val(3), Val(Float64))
        α, β = 3.0, 2.0
        ℂiso = α * 𝕁 + β * 𝕂
        L = TensTI{4}(2.0, 3.0, 1.0, 4.0, 5.0, n3)   # N=5

        # L⊡ℂiso via Walpole product rule (convert ISO first)
        res_w = L ⊡ ℂiso
        @test res_w isa TensTI{4}
        # Compare with direct array contraction
        res_direct = Tensor{4, 3}(get_array(L)) ⊡ Tensor{4, 3}(get_array(ℂiso))
        @test opequal(get_array(res_w), Array(res_direct))
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "tens_TI & arg_TI — numeric round-trip" begin
        C₁₁₁₁, C₁₁₂₂, C₁₁₃₃, C₃₃₃₃, C₂₃₂₃ = 10.0, 3.0, 2.5, 12.0, 2.0
        ℂ = tens_TI(C₁₁₁₁, C₁₁₂₂, C₁₁₃₃, C₃₃₃₃, C₂₃₂₃, n3)
        @test ℂ isa TensTI{4, Float64, 5}

        c₁₁₁₁, c₁₁₂₂, c₁₁₃₃, c₃₃₃₃, c₂₃₂₃ = arg_TI(ℂ)
        @test c₁₁₁₁ ≈ C₁₁₁₁  atol = atol_num
        @test c₁₁₂₂ ≈ C₁₁₂₂  atol = atol_num
        @test c₁₁₃₃ ≈ C₁₁₃₃  atol = atol_num
        @test c₃₃₃₃ ≈ C₃₃₃₃  atol = atol_num
        @test c₂₃₂₃ ≈ C₂₃₂₃  atol = atol_num
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "tens_TI — consistency with TensOrtho" begin
        frame3 = CanonicalBasis{3, Float64}()
        C₁₁₁₁ = 10.0; C₃₃₃₃ = 12.0; C₁₁₂₂ = 3.0; C₁₁₃₃ = 2.5; C₂₃₂₃ = 2.0
        C₆₆ = (C₁₁₁₁ - C₁₁₂₂) / 2
        to = TensOrtho(
            C₁₁₁₁, C₁₁₁₁, C₃₃₃₃, C₁₁₂₂, C₁₁₃₃, C₁₁₃₃,
            C₂₃₂₃, C₂₃₂₃, C₆₆, frame3
        )
        tw = tens_TI(C₁₁₁₁, C₁₁₂₂, C₁₁₃₃, C₃₃₃₃, C₂₃₂₃, n3)
        @test opequal(get_array(to), get_array(tw))
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "tens_TI_eng & arg_TI_eng — numeric round-trip" begin
        E₁, E₃, ν₁₂, ν₃₁, G₃₁ = 100.0, 200.0, 0.25, 0.15, 40.0
        𝕊 = tens_TI_eng(E₁, E₃, ν₁₂, ν₃₁, G₃₁, n3)
        @test 𝕊 isa TensTI{4, Float64, 5}

        e₁, e₃, n₁₂, n₃₁, g₃₁ = arg_TI_eng(𝕊)
        @test e₁ ≈ E₁   atol = atol_num
        @test e₃ ≈ E₃   atol = atol_num
        @test n₁₂ ≈ ν₁₂  atol = atol_num
        @test n₃₁ ≈ ν₃₁  atol = atol_num
        @test g₃₁ ≈ G₃₁  atol = atol_num
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "tens_TI_eng — compliance components check" begin
        E₁, E₃, ν₁₂, ν₃₁, G₃₁ = 100.0, 200.0, 0.25, 0.15, 40.0
        𝕊 = tens_TI_eng(E₁, E₃, ν₁₂, ν₃₁, G₃₁, n3)

        S₁₁₁₁, S₁₁₂₂, S₁₁₃₃, S₃₃₃₃, S₂₃₂₃ = arg_TI(𝕊)
        @test S₁₁₁₁ ≈ 1 / E₁          atol = atol_num
        @test S₃₃₃₃ ≈ 1 / E₃          atol = atol_num
        @test S₁₁₂₂ ≈ -ν₁₂ / E₁      atol = atol_num
        @test S₁₁₃₃ ≈ -ν₃₁ / E₃      atol = atol_num
        @test S₂₃₂₃ ≈ 1 / (4 * G₃₁)    atol = atol_num

        # inv(𝕊) ⊡ 𝕊 = 𝕀
        ℂ = inv(𝕊)
        𝕀 = tens_Id4(Val(3), Val(Float64))
        prod = ℂ ⊡ 𝕊
        @test opequal(get_array(prod), get_array(𝕀))
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "tens_TI & arg_TI — symbolic round-trip" begin
        C₁, C₂, C₃, C₄, C₅ = symbols("C₁₁₁₁ C₁₁₂₂ C₁₁₃₃ C₃₃₃₃ C₂₃₂₃", real = true)
        ℂ = tens_TI(C₁, C₂, C₃, C₄, C₅, n3s)
        @test ℂ isa TensTI{4, <:Any, 5}
        c₁, c₂, c₃, c₄, c₅ = arg_TI(ℂ)
        @test simplify(c₁ - C₁) == 0
        @test simplify(c₂ - C₂) == 0
        @test simplify(c₃ - C₃) == 0
        @test simplify(c₄ - C₄) == 0
        @test simplify(c₅ - C₅) == 0
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "tens_TI_Hoenig & arg_TI_Hoenig — numeric round-trip" begin
        E, ν₁, ν₂, H, Γ = 100.0, 0.25, 0.15, 2.0, 3.0
        𝕊 = tens_TI_Hoenig(E, ν₁, ν₂, H, Γ, n3)
        @test 𝕊 isa TensTI{4, Float64, 5}

        e, n1, n2, h, g = arg_TI_Hoenig(𝕊)
        @test e ≈ E   atol = atol_num
        @test n1 ≈ ν₁  atol = atol_num
        @test n2 ≈ ν₂  atol = atol_num
        @test h ≈ H   atol = atol_num
        @test g ≈ Γ   atol = atol_num
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "tens_TI_Hoenig — consistency with tens_TI_eng" begin
        E, ν₁, ν₂, H, Γ = 100.0, 0.25, 0.15, 2.0, 3.0
        𝕊h = tens_TI_Hoenig(E, ν₁, ν₂, H, Γ, n3)

        # Convert Hoenig → standard engineering constants
        E₁ = E
        E₃ = E * H
        ν₁₂ = ν₁
        ν₃₁ = H * ν₂
        G₃₁ = E * Γ / (2 * (1 + ν₁))
        𝕊e = tens_TI_eng(E₁, E₃, ν₁₂, ν₃₁, G₃₁, n3)

        @test opequal(get_array(𝕊h), get_array(𝕊e))
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "tens_TI_Hoenig — isotropic limit" begin
        # Setting ν₁=ν₂=ν, H=Γ=1 should give the isotropic compliance tensor
        E, ν = 100.0, 0.3
        𝕊 = tens_TI_Hoenig(E, ν, ν, 1.0, 1.0, n3)
        k = E / (3 * (1 - 2ν))
        μ = E / (2 * (1 + ν))
        𝕁, 𝕂 = tens_J4(Val(3), Val(Float64)), tens_K4(Val(3), Val(Float64))
        𝕊iso = (1 / (3k)) * 𝕁 + (1 / (2μ)) * 𝕂
        @test opequal(get_array(𝕊), get_array(𝕊iso))
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "otimes — TensTI{2} self-product" begin
        n = [0.0, 0.0, 1.0]
        A = TensTI{2}(5.0, 8.0, n)
        R = otimes(A)
        @test R isa TensTI{4, Float64, 5}
        sq2 = sqrt(2.0)
        @test get_data(R) == (64.0, 50.0, sq2 * 40.0, 0.0, 0.0)
        @test opequal(get_array(R), otimes(get_array(A), get_array(A)))
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "otimes — TensTI{2} × TensTI{2} (same axis)" begin
        n = [0.0, 0.0, 1.0]
        A = TensTI{2}(5.0, 8.0, n)
        B = TensTI{2}(3.0, 2.0, n)
        R = A ⊗ B
        @test R isa TensTI{4, Float64, 6}
        sq2 = sqrt(2.0)
        # ℓ₁=b₁b₂, ℓ₂=2a₁a₂, ℓ₃=√2·b₁a₂, ℓ₄=√2·a₁b₂
        @test R.data[1] ≈ 16.0     atol = atol_num   # 8*2
        @test R.data[2] ≈ 30.0     atol = atol_num   # 2*5*3
        @test R.data[3] ≈ sq2 * 24.0 atol = atol_num   # √2*8*3
        @test R.data[4] ≈ sq2 * 10.0 atol = atol_num   # √2*5*2
        @test R.data[5] ≈ 0.0      atol = atol_num
        @test R.data[6] ≈ 0.0      atol = atol_num
        @test opequal(get_array(R), otimes(get_array(A), get_array(B)))
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "otimes — TensISO{2,3} × TensTI{2} and reverse" begin
        n = [0.0, 0.0, 1.0]
        I2 = tens_Id2(Val(3), Val(Float64))   # TensISO{2,3}(1.0)
        B = TensTI{2}(5.0, 8.0, n)
        R1 = I2 ⊗ B
        @test R1 isa TensTI{4, Float64, 6}
        @test opequal(get_array(R1), otimes(get_array(I2), get_array(B)))
        R2 = B ⊗ I2
        @test R2 isa TensTI{4, Float64, 6}
        @test opequal(get_array(R2), otimes(get_array(B), get_array(I2)))
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "otimes — TensTI{2} isotropic limit" begin
        n = [0.0, 0.0, 1.0]
        λ = 3.0
        A = TensTI{2}(λ, λ, n)        # isotropic: a == b
        Rw = otimes(A)                  # → TensTI{4, T, 5}
        Riso = otimes(TensISO{3}(λ))   # → TensISO{4,3}
        @test opequal(get_array(Rw), get_array(Riso))
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "otimes — TensTI{2} non-canonical axis" begin
        n = [1.0, 1.0, 1.0] / sqrt(3.0)
        A = TensTI{2}(4.0, 7.0, n)
        B = TensTI{2}(2.0, 5.0, n)
        R = A ⊗ B
        @test R isa TensTI{4, Float64, 6}
        @test opequal(get_array(R), otimes(get_array(A), get_array(B)))
        Rs = otimes(A)
        @test Rs isa TensTI{4, Float64, 5}
        @test opequal(get_array(Rs), otimes(get_array(A), get_array(A)))
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "otimes — TensTI{2} different axes (fallback)" begin
        n1 = [0.0, 0.0, 1.0]
        n2 = [1.0, 0.0, 0.0]
        A = TensTI{2}(5.0, 8.0, n1)
        B = TensTI{2}(3.0, 2.0, n2)
        R = A ⊗ B
        # Different axes → generic fallback → Tens (not TensTI{4})
        @test !(R isa TensTI{4})
        @test opequal(get_array(R), otimes(get_array(A), get_array(B)))
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "Cross-type dispatch — promotion output types" begin
        frame = CanonicalBasis{3, Float64}()
        n = (0.0, 0.0, 1.0)

        I4 = TensISO{3}(2.0, 3.0)   # α=2, β=3
        W = tens_TI(10.0, 3.0, 2.5, 12.0, 2.0, n)                 # TensTI{4, T, 5}
        O = TensOrtho(10.0, 8.0, 12.0, 3.0, 2.5, 1.5, 2.0, 3.0, 3.5, frame)

        # ── + / − between ISO and Ortho ────────────────────────────────────
        @test I4 + O isa TensOrtho{Float64}
        @test O + I4 isa TensOrtho{Float64}
        @test I4 - O isa TensOrtho{Float64}
        @test O - I4 isa TensOrtho{Float64}

        # Numerical agreement with generic array arithmetic
        @test opequal(get_array(I4 + O), get_array(I4) + get_array(O))
        @test opequal(get_array(O - I4), get_array(O) - get_array(I4))

        # ── + / − between Walpole(N=5) and Ortho (aligned axis) ─────────────
        @test W + O isa TensOrtho{Float64}
        @test O - W isa TensOrtho{Float64}
        @test opequal(get_array(W + O), get_array(W) + get_array(O))
        @test opequal(get_array(O - W), get_array(O) - get_array(W))

        # Reference mismatch → assertion
        n_mis = (1.0, 1.0, 1.0) ./ sqrt(3.0)
        Wmis = tens_TI(10.0, 3.0, 2.5, 12.0, 2.0, n_mis)
        @test_throws AssertionError Wmis + O

        # ── TensTI{4} N=5 + N=6 (same axis) → N=6 ─────────────────────────
        W5 = tens_TI(10.0, 3.0, 2.5, 12.0, 2.0, n)          # N=5
        W6 = TensTI{4}(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, n)  # N=6
        @test (W5 + W6) isa TensTI{4, Float64, 6}
        @test (W6 - W5) isa TensTI{4, Float64, 6}
        @test opequal(get_array(W5 + W6), get_array(W5) + get_array(W6))

        # ── TensISO{2,3} ± TensTI{2} ────────────────────────────────────────
        I2 = TensISO{3}(5.0)
        T2 = TensTI{2}(3.0, 7.0, n)
        @test (I2 + T2) isa TensTI{2, Float64, 2}
        @test (T2 - I2) isa TensTI{2, Float64, 2}
        @test opequal(get_array(I2 + T2), get_array(I2) + get_array(T2))
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "Cross-type dispatch — cross-order double contraction" begin
        n = (0.0, 0.0, 1.0)
        I4 = TensISO{3}(2.0, 3.0)
        W = tens_TI(10.0, 3.0, 2.5, 12.0, 2.0, n)
        T2 = TensTI{2}(4.0, 7.0, n)

        # ── TensISO{4} ⊡ TensTI{2} / reverse ─────────────────────────────────
        R1 = I4 ⊡ T2
        @test R1 isa TensTI{2, Float64, 2}
        @test opequal(get_array(R1), dcontract(get_array(I4), get_array(T2)))

        R2 = T2 ⊡ I4
        @test R2 isa TensTI{2, Float64, 2}
        @test opequal(get_array(R2), dcontract(get_array(T2), get_array(I4)))

        # ── TensTI{4} ⊡ TensTI{2} (same axis) ──────────────────────────────
        R3 = W ⊡ T2
        @test R3 isa TensTI{2, Float64, 2}
        @test opequal(get_array(R3), dcontract(get_array(W), get_array(T2)))

        R4 = T2 ⊡ W
        @test R4 isa TensTI{2, Float64, 2}
        @test opequal(get_array(R4), dcontract(get_array(T2), get_array(W)))

        # Axis mismatch → generic fallback (no assertion), values exact
        n2 = (1.0, 0.0, 0.0)
        Wm = tens_TI(10.0, 3.0, 2.5, 12.0, 2.0, n2)
        T2m = TensTI{2}(4.0, 7.0, n)
        Rm = Wm ⊡ T2m
        @test !(Rm isa TensTI)
        @test opequal(Array(get_array(Rm)), dcontract(get_array(Wm), get_array(T2m)))
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "Cross-type dispatch — single contraction (·)" begin
        n = (0.0, 0.0, 1.0)
        A = TensTI{2}(4.0, 7.0, n)
        B = TensTI{2}(2.0, 5.0, n)

        # TensTI{2} · TensTI{2} (same axis)
        R = A ⋅ B
        @test R isa TensTI{2, Float64, 2}
        @test opequal(get_array(R), get_array(A) * get_array(B))

        # TensTI{2} · TensISO{2,3}
        I2 = TensISO{3}(3.0)
        @test (A ⋅ I2) isa TensTI{2, Float64, 2}
        @test (I2 ⋅ A) isa TensTI{2, Float64, 2}
        @test opequal(get_array(A ⋅ I2), get_array(A) * get_array(I2))
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "Cross-type dispatch — inverse consistency" begin
        n = (0.0, 0.0, 1.0)
        frame = CanonicalBasis{3, Float64}()
        Id4 = tens_Id4(Val(3), Val(Float64))

        # TensTI{4}
        W = tens_TI(100.0, 30.0, 25.0, 120.0, 20.0, n)
        @test opequal(get_array(inv(W) ⊡ W), get_array(Id4))

        # TensOrtho
        O = TensOrtho(10.0, 8.0, 12.0, 3.0, 2.5, 1.5, 2.0, 3.0, 3.5, frame)
        @test opequal(get_array(inv(O)) ⊡ get_array(O), get_array(Id4))

        # TensISO{4}
        I4 = TensISO{3}(2.0, 3.0)
        @test opequal(get_array(inv(I4) ⊡ I4), get_array(Id4))
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "Promotion helpers — iso_to_ortho / walpole_to_ortho" begin
        frame = CanonicalBasis{3, Float64}()

        # ── iso_to_ortho: ISO tensor rebuilt as TensOrtho yields identical array
        I4 = TensISO{3}(2.5, 4.0)
        O_from_iso = iso_to_ortho(I4, frame)
        @test O_from_iso isa TensOrtho{Float64}
        @test opequal(get_array(I4), get_array(O_from_iso))

        # ── walpole_to_ortho (aligned axis)
        n = (0.0, 0.0, 1.0)
        W = tens_TI(10.0, 3.0, 2.5, 12.0, 2.0, n)
        k = TensND._axis_on_frame_index(n, frame)
        @test k == 3
        O_from_w = walpole_to_ortho(W, frame, k)
        @test O_from_w isa TensOrtho{Float64}
        @test opequal(get_array(W), get_array(O_from_w))

        # Axis along e₁
        n1 = (1.0, 0.0, 0.0)
        W1 = tens_TI(10.0, 3.0, 2.5, 12.0, 2.0, n1)
        O_from_w1 = walpole_to_ortho(W1, frame, 1)
        @test opequal(get_array(W1), get_array(O_from_w1))
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensTI{4,T,8} — full axially-invariant algebra" begin
        ez = (0.0, 0.0, 1.0)
        nax = (0.6, 0.0, 0.8)   # unit, non-canonical axis

        # Mandel patterns of the antisymmetric generators (axis = e₃)
        M7 = KM(tens_W7(ez))
        M8 = KM(tens_W8(ez))
        s2 = 1 / sqrt(2.0)
        @test M7[4, 5] ≈ -1.0 atol = atol_num
        @test M7[5, 4] ≈ 1.0 atol = atol_num
        @test sum(abs, M7) ≈ 2.0 atol = atol_num
        @test M8[6, 1] ≈ s2 atol = atol_num
        @test M8[6, 2] ≈ -s2 atol = atol_num
        @test M8[1, 6] ≈ -s2 atol = atol_num
        @test M8[2, 6] ≈ s2 atol = atol_num
        @test sum(abs, M8) ≈ 4s2 atol = atol_num

        # constructor / accessors / lifts
        A8 = TensTI{4}(1.0, 2.0, 0.3, -0.4, 1.5, 2.5, 0.7, -0.2, nax)
        @test A8 isa TensTI{4, Float64, 8}
        @test get_ℓ8(A8) == (1.0, 2.0, 0.3, -0.4, 1.5, 2.5, 0.7, -0.2)
        A5 = TensTI{4}(1.0, 2.0, 0.5, 3.0, 4.0, nax)
        @test get_ℓ8(A5) == (1.0, 2.0, 0.5, 0.5, 3.0, 4.0, 0.0, 0.0)
        @test TensND._lift_walpole_N8(A5) isa TensTI{4, Float64, 8}
        @test !issymmetric(A8)
        @test !Tensors.ismajorsymmetric(A8)
        @test issymmetric(TensTI{4}(1.0, 2.0, 0.5, 0.5, 3.0, 4.0, 0.0, 0.0, nax))

        # algebraic closure: product rule vs array dcontract, inverse
        B8 = TensTI{4}(0.8, -1.2, 0.1, 0.9, 2.0, 0.4, -0.6, 1.1, nax)
        C8 = A8 ⊡ B8
        @test C8 isa TensTI{4, Float64, 8}
        Ca = Tensor{4, 3}(get_array(A8)) ⊡ Tensor{4, 3}(get_array(B8))
        @test maximum(abs, get_array(C8) .- Array(Ca)) < 1.0e-12
        Id = A8 ⊡ inv(A8)
        @test maximum(abs, KM(Id) .- I(6)) < 1.0e-12

        # minor symmetry of the array is exact (6×6 Mandel route preserved)
        @test Tensors.issymmetric(Tensor{4, 3}(get_array(A8)))

        # mixed-N ± lift
        S = A5 + A8
        @test S isa TensTI{4, Float64, 8}
        @test maximum(abs, get_array(S) .- (get_array(A5) .+ get_array(A8))) < 1.0e-12
        ISO4 = TensISO{3}(2.0, 3.0)
        SI = ISO4 + A8
        @test SI isa TensTI{4, Float64, 8}
        @test maximum(abs, get_array(SI) .- (get_array(ISO4) .+ get_array(A8))) < 1.0e-12

        # mixed-axis fallback (no assertion, generic result, values exact)
        B8b = TensTI{4}(0.8, -1.2, 0.1, 0.9, 2.0, 0.4, -0.6, 1.1, ez)
        S2 = A8 + B8b
        @test !(S2 isa TensTI)
        @test maximum(abs, Array(get_array(S2)) .- (get_array(A8) .+ get_array(B8b))) < 1.0e-12
        P2 = A8 ⊡ B8b
        Pref = Tensor{4, 3}(get_array(A8)) ⊡ Tensor{4, 3}(get_array(B8b))
        @test maximum(abs, Array(get_array(P2)) .- Array(Pref)) < 1.0e-12

        # W₇, W₈ annihilate symmetric 2nd-order tensors
        tsym = TensTI{2}(1.3, -0.4, nax)
        @test maximum(abs, get_array(tens_W7(nax) ⊡ tsym)) < 1.0e-13
        @test maximum(abs, get_array(tens_W8(nax) ⊡ tsym)) < 1.0e-13

        # ForwardDiff genericity through the closed algebra
        FD = TensND.ForwardDiff
        f = x -> get_ℓ8(inv(TensTI{4}(x, 2.0, 0.3, -0.4, 1.5, 2.5, 0.7, -0.2, ez)))[1]
        g = FD.derivative(f, 1.0)
        h = 1.0e-6
        @test g ≈ (f(1.0 + h) - f(1.0 - h)) / 2h atol = 1.0e-6
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensTI{2,T,3} — axially-invariant 2nd order with rotation part" begin
        ez = (0.0, 0.0, 1.0)
        nax = (0.6, 0.0, 0.8)

        t3 = TensTI{2}(2.0, 5.0, 0.7, nax)
        @test t3 isa TensTI{2, Float64, 3}
        @test !issymmetric(t3)
        @test tr(t3) ≈ 9.0

        # w convention: w·p = n×p → arr[1,2] = −c, arr[2,1] = +c for n = e₃
        arr = get_array(TensTI{2}(0.0, 0.0, 1.0, ez))
        @test arr[1, 2] ≈ -1.0
        @test arr[2, 1] ≈ 1.0

        # dot = complex product in the plane ⊕ scalar on the axis
        t3b = TensTI{2}(1.5, -1.0, -0.3, nax)
        d = t3 ⋅ t3b
        @test d isa TensTI{2, Float64, 3}
        @test maximum(abs, get_array(d) .- (get_array(t3) * get_array(t3b))) < 1.0e-12

        # inverse
        @test maximum(abs, (get_array(t3) * get_array(inv(t3))) .- I(3)) < 1.0e-12

        # lifts and ISO interplay
        t2 = TensTI{2}(1.0, 2.0, nax)
        s = t2 + t3
        @test s isa TensTI{2, Float64, 3}
        @test maximum(abs, get_array(s) .- (get_array(t2) .+ get_array(t3))) < 1.0e-12
        iso2 = TensISO{3}(4.0)
        s2 = iso2 + t3
        @test s2 isa TensTI{2, Float64, 3}
        @test maximum(abs, get_array(s2) .- (get_array(iso2) .+ get_array(t3))) < 1.0e-12

        # 4th ⊡ 2nd drops the antisymmetric part
        A5 = TensTI{4}(1.0, 2.0, 0.5, 3.0, 4.0, nax)
        r = A5 ⊡ t3
        rsym = A5 ⊡ TensTI{2}(2.0, 5.0, nax)
        @test maximum(abs, get_array(r) .- get_array(rsym)) < 1.0e-13

        # mixed-axis dot falls back generically
        t3c = TensTI{2}(1.5, -1.0, -0.3, ez)
        dm = t3 ⋅ t3c
        @test !(dm isa TensTI)
        @test maximum(abs, Array(get_array(dm)) .- (get_array(t3) * get_array(t3c))) < 1.0e-12
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensOrtho — ForwardDiff compatibility (frame eltype decoupled)" begin
        FD = TensND.ForwardDiff
        frame = CanonicalBasis{3, Float64}()
        # build TensOrtho with a Dual elastic constant against a Float64 frame
        O = TensOrtho(
            FD.Dual(250.0e3, 1.0), 250.0e3, 250.0e3, 100.0e3, 100.0e3, 100.0e3,
            80.0e3, 80.0e3, 80.0e3, frame
        )
        @test eltype(O) <: FD.Dual
        @test frame isa TensND.OrthonormalBasis{3, Float64}   # frame stays Float64
        # ∂C₁₁₁₁/∂C₁₁ = 1
        g = FD.derivative(
            c11 -> get_array(
                TensOrtho(c11, 250.0e3, 250.0e3, 100.0e3, 100.0e3, 100.0e3, 80.0e3, 80.0e3, 80.0e3, frame)
            )[1, 1, 1, 1],
            250.0e3
        )
        @test g ≈ 1.0
        # inverse is also differentiable (∂ over the KM block inverse)
        gi = FD.derivative(
            c11 -> get_array(
                inv(TensOrtho(c11, 250.0e3, 250.0e3, 100.0e3, 100.0e3, 100.0e3, 80.0e3, 80.0e3, 80.0e3, frame))
            )[1, 1, 1, 1],
            250.0e3
        )
        h = 1.0
        cc(x) = get_array(
            inv(TensOrtho(x, 250.0e3, 250.0e3, 100.0e3, 100.0e3, 100.0e3, 80.0e3, 80.0e3, 80.0e3, frame))
        )[1, 1, 1, 1]
        @test gi ≈ (cc(250.0e3 + h) - cc(250.0e3 - h)) / 2h atol = 1.0e-9
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensOrtho — concrete frame field, closed-form inv (v0.2.4 regression)" begin
        frame = CanonicalBasis{3, Float64}()
        O = TensOrtho(
            250.0e3, 180.0e3, 150.0e3, 90.0e3, 80.0e3, 70.0e3,
            60.0e3, 50.0e3, 40.0e3, frame
        )
        # The frame field must be a concrete type: TensOrtho{T} decouples the
        # frame eltype from T (for ForwardDiff) via a *type parameter*, not by
        # erasing it to the abstract `OrthonormalBasis{3}` — an abstract field
        # boxes every access and defeats inference (v0.2.3 regression).
        @test isconcretetype(fieldtype(typeof(O), :frame))
        @test typeof(O) === TensOrtho{Float64, CanonicalBasis{3, Float64}}

        # get_array/getindex must not allocate wildly more than the necessary
        # 81-Float64 array (~650 B). The boxed-frame regression allocated
        # 283 264 B per call (≈440×) because every field access on the
        # abstract `frame` triggered dynamic dispatch.
        get_array(O)  # warm up
        @test (@allocated get_array(O)) < 2000
        O[1, 1, 1, 1]
        @test (@allocated O[1, 1, 1, 1]) < 2000

        # Closed-form inv (3×3 block adjugate + shear reciprocals) must match
        # the dense 6×6 Kelvin-Mandel inverse to machine precision.
        Km = TensND.KM_material(O)
        Km_inv_dense = inv(Km)
        Km_inv_closed = TensND.KM_material(inv(O))
        @test maximum(abs.(Km_inv_closed .- Km_inv_dense)) < 1.0e-10
        @test maximum(abs.(Km * Km_inv_closed - Matrix(I, 6, 6))) < 1.0e-12
    end

end  # "Walpole & Ortho tensors"
