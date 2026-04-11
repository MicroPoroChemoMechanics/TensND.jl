@testsection "Walpole & Ortho tensors" begin

    # ─── helpers shared across sub-sections ───────────────────────────────────
    n3 = 𝐞(Val(3), Val(3), Val(Float64))   # e₃ as Float64 Vec
    n3s = 𝐞(Val(3), Val(3), Val(Sym))       # e₃ as Sym Vec
    atol_num = 1e-12

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensWalpole — construction & traits" begin
        W1, W2, W3, W4, W5, W6 = Walpole(n3)
        @test W1 isa TensWalpole{Float64,6}
        @test W2 isa TensWalpole{Float64,6}
        @test size(W1) == (3, 3, 3, 3)
        @test getbasis(W1) isa CanonicalBasis{3,Float64}
        @test getvar(W1) == (:cont, :cont, :cont, :cont)

        # N=5 (symmetric) basis
        W1s, W2s, W3s, W4s, W5s = Walpole(n3; sym = true)
        @test W1s isa TensWalpole{Float64,5}
        @test size(W1s) == (3, 3, 3, 3)

        # individual constructors
        @test tensW1(n3) isa TensWalpole{Float64,6}
        @test tensW3(n3) isa TensWalpole{Float64,6}

        # axis accessor
        @test getaxis(W1) == (0.0, 0.0, 1.0)
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensWalpole — getarray & KM (n=e₃, Float64)" begin
        W1, W2, W3, W4, W5, W6 = Walpole(n3)
        sq2 = sqrt(2.0)

        # W₁ = nₙ⊗nₙ  →  only entry [3,3,3,3] = 1
        A1 = getarray(W1)
        @test A1[3,3,3,3] ≈ 1.0   atol=atol_num
        @test all(abs.(A1[i,j,k,l]) < atol_num
                  for i in 1:3, j in 1:3, k in 1:3, l in 1:3
                  if !(i==3 && j==3 && k==3 && l==3))

        # W₂ = (nT⊗nT)/2 with n=e₃ → nT = diag(1,1,0)
        # so W₂[i,j,k,l] = δᵢⱼ(1-δᵢ₃)δₖₗ(1-δₖ₃)/2
        A2 = getarray(W2)
        @test A2[1,1,1,1] ≈ 0.5   atol=atol_num
        @test A2[1,1,2,2] ≈ 0.5   atol=atol_num
        @test A2[2,2,2,2] ≈ 0.5   atol=atol_num
        @test A2[3,3,3,3] ≈ 0.0   atol=atol_num

        # W₃[3,3,1,1] = 1/√2
        A3 = getarray(W3)
        @test A3[3,3,1,1] ≈ 1.0/sq2   atol=atol_num
        @test A3[3,3,2,2] ≈ 1.0/sq2   atol=atol_num
        @test A3[3,3,3,3] ≈ 0.0       atol=atol_num

        # W₄[1,1,3,3] = 1/√2
        A4 = getarray(W4)
        @test A4[1,1,3,3] ≈ 1.0/sq2   atol=atol_num

        # W₅: shear in transverse plane (e₁,e₂)
        # W₅[1,2,1,2] = W₅[2,1,1,2] = W₅[1,2,2,1] = W₅[2,1,2,1] = 1/2
        A5 = getarray(W5)
        @test A5[1,2,1,2] ≈ 0.5   atol=atol_num
        @test A5[1,1,1,1] ≈ 0.5   atol=atol_num
        @test A5[3,3,3,3] ≈ 0.0   atol=atol_num

        # W₆: shear between transverse and axial
        A6 = getarray(W6)
        @test A6[1,3,1,3] ≈ 0.5   atol=atol_num
        @test A6[2,3,2,3] ≈ 0.5   atol=atol_num
        @test A6[1,1,1,1] ≈ 0.0   atol=atol_num
        @test A6[3,3,3,3] ≈ 0.0   atol=atol_num

        # KM structure for n=e₃: blocks should separate
        L = TensWalpole(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, n3)
        Km = KM(L)
        @test size(Km) == (6, 6)
        # Off-diagonal shear coupling should be zero for axis n=e₃
        @test abs(Km[1,4]) < atol_num
        @test abs(Km[1,5]) < atol_num
        @test abs(Km[4,6]) < atol_num
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensWalpole — Walpole product rule (Float64)" begin
        W1, W2, W3, W4, W5, W6 = Walpole(n3)

        # Idempotents: W₁⊡W₁=W₁, W₂⊡W₂=W₂, W₅⊡W₅=W₅, W₆⊡W₆=W₆
        @test opequal(getarray(W1 ⊡ W1), getarray(W1))
        @test opequal(getarray(W2 ⊡ W2), getarray(W2))
        @test opequal(getarray(W5 ⊡ W5), getarray(W5))
        @test opequal(getarray(W6 ⊡ W6), getarray(W6))

        # Cross products
        @test opequal(getarray(W3 ⊡ W4), getarray(W1))
        @test opequal(getarray(W4 ⊡ W3), getarray(W2))

        # Zero cross products between incompatible blocks
        zero4 = zeros(3, 3, 3, 3)
        @test opequal(getarray(W1 ⊡ W2), zero4)
        @test opequal(getarray(W1 ⊡ W5), zero4)
        @test opequal(getarray(W5 ⊡ W6), zero4)
        @test opequal(getarray(W6 ⊡ W5), zero4)

        # General product: Walpole vs direct array contraction
        L = TensWalpole(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, n3)
        M = TensWalpole(0.5, 1.5, 2.0, 0.3, 0.8, 1.2, n3)
        LM_walpole = getarray(L ⊡ M)
        LM_direct  = Tensor{4,3}(getarray(L)) ⊡ Tensor{4,3}(getarray(M))
        @test opequal(LM_walpole, Array(LM_direct))
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensWalpole — inverse (Float64)" begin
        𝕀 = tensId4(Val(3), Val(Float64))
        n = [1.0/sqrt(3.0), 1.0/sqrt(3.0), 1.0/sqrt(3.0)]

        # N=5 (symmetric)
        L5 = TensWalpole(2.0, 3.0, 1.0, 4.0, 5.0, n)
        Li5 = inv(L5)
        @test Li5 isa TensWalpole{Float64,5}
        LLi5 = getarray(L5 ⊡ Li5)
        𝕀arr = getarray(𝕀)
        @test opequal(LLi5, 𝕀arr)

        # N=6 (general)
        L6 = TensWalpole(2.0, 3.0, 1.0, 0.5, 4.0, 5.0, n)
        Li6 = inv(L6)
        @test Li6 isa TensWalpole{Float64,6}
        LLi6 = getarray(L6 ⊡ Li6)
        @test opequal(LLi6, 𝕀arr)
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensWalpole — fromISO" begin
        𝕀, 𝕁, 𝕂 = ISO(Val(3), Val(Float64))
        α, β = 3.0, 2.0
        ℂiso = α * 𝕁 + β * 𝕂

        # For any axis, fromISO should give the same array as the isotropic tensor
        for n ∈ ([0.0, 0.0, 1.0], [1.0, 0.0, 0.0],
                 [1.0/sqrt(3), 1.0/sqrt(3), 1.0/sqrt(3)])
            Wiso = fromISO(ℂiso, n)
            @test Wiso isa TensWalpole{Float64,5}
            @test isTI(Wiso)
            @test opequal(getarray(Wiso), getarray(ℂiso))
        end

        # Symbolic — use Sym ISO tensors to avoid Float64/Sym residuals
        αs, βs = symbols("α β", real = true)
        𝕁s, 𝕂s = tensJ4(Val(3), Val(Sym)), tensK4(Val(3), Val(Sym))
        ℂisos = αs * 𝕁s + βs * 𝕂s
        Wisos = fromISO(ℂisos, n3s)
        @test Wisos isa TensWalpole{<:Any,5}
        for i in 1:3, j in 1:3, k in 1:3, l in 1:3
            @test simplify(Wisos[i,j,k,l] - ℂisos[i,j,k,l]) == 0
        end
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensWalpole — isISO / isTI / isOrtho" begin
        W = tensW1(n3)
        @test !isISO(W)
        @test  isTI(W)
        @test !isOrtho(W)
        𝕀, 𝕁, 𝕂 = ISO(Val(3), Val(Float64))
        @test !isTI(𝕁)
        @test !isOrtho(𝕁)
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensOrtho — isISO / isTI / isOrtho" begin
        frame3 = CanonicalBasis{3,Float64}()
        t = TensOrtho(10., 8., 9., 3., 2., 4., 2.5, 3., 1.5, frame3)
        @test !isISO(t)
        @test !isTI(t)
        @test  isOrtho(t)
        # universal fallback
        @test !isOrtho(42)
        @test !isOrtho("string")
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensWalpole — show" begin
        L5 = TensWalpole(1.0, 2.0, 0.5, 3.0, 4.0, n3)   # N=5
        L6 = TensWalpole(1.0, 2.0, 0.5, 0.3, 3.0, 4.0, n3)   # N=6
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
        frame3 = CanonicalBasis{3,Float64}()
        t = TensOrtho(10., 8., 9., 3., 2., 4., 2.5, 3., 1.5, frame3)
        buf = IOBuffer()
        show(buf, t)
        s = String(take!(buf))
        @test contains(s, "P₁⊗P₁") && contains(s, "frame")
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensWalpole — tsimplify (symbolic, _rebuild)" begin
        ℓ₁, ℓ₂, ℓ₃ = symbols("ℓ₁ ℓ₂ ℓ₃", real = true)
        L = TensWalpole(ℓ₁, ℓ₂, ℓ₃, ℓ₁+ℓ₂, ℓ₂+ℓ₃, n3s)   # N=6
        Ls = tsimplify(L)
        @test Ls isa TensWalpole   # _rebuild preserves type
        @test get_ℓ(Ls)[1] == ℓ₁   # simplification is a no-op here
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensOrtho — tsimplify (_rebuild)" begin
        # TensOrtho frames are always numeric; test _rebuild with Float64 data.
        # tsimplify on a non-symbolic NTuple is a no-op, but _rebuild must still
        # reconstruct the TensOrtho, verifying the dispatch path.
        frame3 = CanonicalBasis{3,Float64}()
        t = TensOrtho(10., 8., 9., 3., 2., 4., 2.5, 3., 1.5, frame3)
        ts = tsimplify(t)
        @test ts isa TensOrtho   # _rebuild preserves type
        @test getdata(ts)[1] ≈ 10.0
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensWalpole — arithmetic" begin
        W1, W2, W3, W4, W5, W6 = Walpole(n3)
        L = TensWalpole(1.0, 2.0, 0.5, 0.5, 3.0, 4.0, n3)   # N=6
        M = TensWalpole(0.5, 1.0, 0.25, 0.25, 1.5, 2.0, n3)

        @test opequal(getarray(L + M), getarray(L) .+ getarray(M))
        @test opequal(getarray(L - M), getarray(L) .- getarray(M))
        @test opequal(getarray(2.0 * L), 2.0 .* getarray(L))
        @test opequal(getarray(-L), .-getarray(L))

        # Symmetric N=5
        Ls = TensWalpole(1.0, 2.0, 0.5, 3.0, 4.0, n3)
        Ms = TensWalpole(0.5, 1.0, 0.25, 1.5, 2.0, n3)
        @test (Ls + Ms) isa TensWalpole{Float64,5}
        @test opequal(getarray(Ls + Ms), getarray(Ls) .+ getarray(Ms))
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensWalpole — symbolic inverse" begin
        ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆ = symbols("ℓ₁ ℓ₂ ℓ₃ ℓ₅ ℓ₆", real = true)
        L = TensWalpole(ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆, n3s)
        Li = inv(L)
        @test Li isa TensWalpole{<:Any,5}
        # L⊡inv(L) should be identity
        𝕀 = tensId4(Val(3), Val(Sym))
        prod = L ⊡ Li
        for i in 1:3, j in 1:3, k in 1:3, l in 1:3
            @test simplify(prod[i,j,k,l] - 𝕀[i,j,k,l]) == 0
        end
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensOrtho — construction & traits" begin
        frame3 = CanonicalBasis{3,Float64}()
        C11, C22, C33 = 10.0, 8.0, 12.0
        C12, C13, C23 = 3.0, 4.0, 2.5
        C44, C55, C66 = 2.0, 3.0, 1.5
        t = TensOrtho(C11, C22, C33, C12, C13, C23, C44, C55, C66, frame3)
        @test t isa TensOrtho{Float64}
        @test size(t) == (3, 3, 3, 3)
        @test getbasis(t) isa CanonicalBasis{3,Float64}
        @test getvar(t) == (:cont, :cont, :cont, :cont)
        @test getframe(t) === frame3
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensOrtho — KM_material (canonical frame)" begin
        frame3 = CanonicalBasis{3,Float64}()
        C11, C22, C33 = 10.0, 8.0, 12.0
        C12, C13, C23 = 3.0, 4.0, 2.5
        C44, C55, C66 = 2.0, 3.0, 1.5
        t = TensOrtho(C11, C22, C33, C12, C13, C23, C44, C55, C66, frame3)
        Km = KM_material(t)
        @test size(Km) == (6, 6)
        # Diagonal blocks
        @test Km[1,1] ≈ C11  atol=atol_num
        @test Km[2,2] ≈ C22  atol=atol_num
        @test Km[3,3] ≈ C33  atol=atol_num
        @test Km[4,4] ≈ 2*C44 atol=atol_num
        @test Km[5,5] ≈ 2*C55 atol=atol_num
        @test Km[6,6] ≈ 2*C66 atol=atol_num
        # Off-diagonal within normal block
        @test Km[1,2] ≈ C12  atol=atol_num
        @test Km[1,3] ≈ C13  atol=atol_num
        @test Km[2,3] ≈ C23  atol=atol_num
        # Zeros between normal and shear blocks
        @test Km[1,4] ≈ 0.0  atol=atol_num
        @test Km[2,5] ≈ 0.0  atol=atol_num
        @test Km[3,6] ≈ 0.0  atol=atol_num
        @test Km[4,5] ≈ 0.0  atol=atol_num
        @test Km[4,6] ≈ 0.0  atol=atol_num
        @test Km[5,6] ≈ 0.0  atol=atol_num
        # Symmetry
        @test Km ≈ Km'  atol=atol_num
        # KM in canonical frame should match
        Km2 = KM(t)
        @test Km2 ≈ Km  atol=atol_num
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensOrtho — from KM matrix" begin
        frame3 = CanonicalBasis{3,Float64}()
        C11, C22, C33 = 10.0, 8.0, 12.0
        C12, C13, C23 = 3.0, 4.0, 2.5
        C44, C55, C66 = 2.0, 3.0, 1.5
        t1 = TensOrtho(C11, C22, C33, C12, C13, C23, C44, C55, C66, frame3)
        Km = KM_material(t1)
        t2 = TensOrtho(Km, frame3)
        @test opequal(getarray(t1), getarray(t2))
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensOrtho — inverse (canonical frame)" begin
        frame3 = CanonicalBasis{3,Float64}()
        C11, C22, C33 = 10.0, 8.0, 12.0
        C12, C13, C23 = 3.0, 4.0, 2.5
        C44, C55, C66 = 2.0, 3.0, 1.5
        t = TensOrtho(C11, C22, C33, C12, C13, C23, C44, C55, C66, frame3)
        ti = inv(t)
        @test ti isa TensOrtho{Float64}
        𝕀 = tensId4(Val(3), Val(Float64))
        A  = Tensor{4,3}(getarray(t))
        Ai = Tensor{4,3}(getarray(ti))
        prod = Array(A ⊡ Ai)
        𝕀arr = getarray(𝕀)
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
        @test Km[1,4] ≈ 0.0  atol=atol_num
        @test Km[4,5] ≈ 0.0  atol=atol_num
        @test Km[5,6] ≈ 0.0  atol=atol_num
        # change_tens to canonical basis should give a valid Tens
        tc = change_tens(t, CanonicalBasis{3,Float64}())
        @test tc isa AbstractTens
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensOrtho — TI consistency" begin
        # A TI tensor (C11=C22, C13=C23, C44=C55, C66=(C11-C12)/2)
        # should give same array as corresponding TensWalpole (with n=e₃)
        frame3 = CanonicalBasis{3,Float64}()
        C11 = 10.0; C33 = 12.0; C12 = 3.0; C13 = 2.5; C44 = 2.0
        C22 = C11; C23 = C13; C55 = C44; C66 = (C11 - C12) / 2
        to = TensOrtho(C11, C22, C33, C12, C13, C23, C44, C55, C66, frame3)

        # Build equivalent TensWalpole (n=e₃) via the engineering constants.
        # For a TI material with n=e₃:
        #   C₁₁=C₂₂=(ℓ₂+ℓ₅)/2, C₃₃=ℓ₁, C₁₂=(ℓ₂-ℓ₅)/2, C₁₃=C₂₃=ℓ₃/√2, C₄₄=ℓ₆/2, C₆₆=ℓ₅/2
        # Inverting: ℓ₂=C₁₁+C₁₂, ℓ₅=C₁₁−C₁₂=2C₆₆, ℓ₃=C₁₃√2, ℓ₆=2C₄₄
        sq2 = sqrt(2.0)
        ℓ₁ = C33
        ℓ₅ = C66 * 2       # = C11 - C12
        ℓ₂ = C11 + C12     # NOT (C11+C12)/2
        ℓ₃ = C13 * sq2
        ℓ₆ = C44 * 2
        tw = TensWalpole(ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆, n3)
        @test opequal(getarray(to), getarray(tw))
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensWalpole — dcontract with TensISO" begin
        𝕀, 𝕁, 𝕂 = ISO(Val(3), Val(Float64))
        α, β = 3.0, 2.0
        ℂiso = α * 𝕁 + β * 𝕂
        L = TensWalpole(2.0, 3.0, 1.0, 4.0, 5.0, n3)   # N=5

        # L⊡ℂiso via Walpole product rule (convert ISO first)
        res_w = L ⊡ ℂiso
        @test res_w isa TensWalpole
        # Compare with direct array contraction
        res_direct = Tensor{4,3}(getarray(L)) ⊡ Tensor{4,3}(getarray(ℂiso))
        @test opequal(getarray(res_w), Array(res_direct))
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "tensTI & argTI — numeric round-trip" begin
        C₁₁₁₁, C₁₁₂₂, C₁₁₃₃, C₃₃₃₃, C₂₃₂₃ = 10.0, 3.0, 2.5, 12.0, 2.0
        ℂ = tensTI(C₁₁₁₁, C₁₁₂₂, C₁₁₃₃, C₃₃₃₃, C₂₃₂₃, n3)
        @test ℂ isa TensWalpole{Float64,5}

        c₁₁₁₁, c₁₁₂₂, c₁₁₃₃, c₃₃₃₃, c₂₃₂₃ = argTI(ℂ)
        @test c₁₁₁₁ ≈ C₁₁₁₁  atol=atol_num
        @test c₁₁₂₂ ≈ C₁₁₂₂  atol=atol_num
        @test c₁₁₃₃ ≈ C₁₁₃₃  atol=atol_num
        @test c₃₃₃₃ ≈ C₃₃₃₃  atol=atol_num
        @test c₂₃₂₃ ≈ C₂₃₂₃  atol=atol_num
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "tensTI — consistency with TensOrtho" begin
        frame3 = CanonicalBasis{3,Float64}()
        C₁₁₁₁ = 10.0; C₃₃₃₃ = 12.0; C₁₁₂₂ = 3.0; C₁₁₃₃ = 2.5; C₂₃₂₃ = 2.0
        C₆₆ = (C₁₁₁₁ - C₁₁₂₂) / 2
        to = TensOrtho(C₁₁₁₁, C₁₁₁₁, C₃₃₃₃, C₁₁₂₂, C₁₁₃₃, C₁₁₃₃,
                       C₂₃₂₃, C₂₃₂₃, C₆₆, frame3)
        tw = tensTI(C₁₁₁₁, C₁₁₂₂, C₁₁₃₃, C₃₃₃₃, C₂₃₂₃, n3)
        @test opequal(getarray(to), getarray(tw))
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "tensTI_eng & argTI_eng — numeric round-trip" begin
        E₁, E₃, ν₁₂, ν₃₁, G₃₁ = 100.0, 200.0, 0.25, 0.15, 40.0
        𝕊 = tensTI_eng(E₁, E₃, ν₁₂, ν₃₁, G₃₁, n3)
        @test 𝕊 isa TensWalpole{Float64,5}

        e₁, e₃, n₁₂, n₃₁, g₃₁ = argTI_eng(𝕊)
        @test e₁ ≈ E₁   atol=atol_num
        @test e₃ ≈ E₃   atol=atol_num
        @test n₁₂ ≈ ν₁₂  atol=atol_num
        @test n₃₁ ≈ ν₃₁  atol=atol_num
        @test g₃₁ ≈ G₃₁  atol=atol_num
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "tensTI_eng — compliance components check" begin
        E₁, E₃, ν₁₂, ν₃₁, G₃₁ = 100.0, 200.0, 0.25, 0.15, 40.0
        𝕊 = tensTI_eng(E₁, E₃, ν₁₂, ν₃₁, G₃₁, n3)

        S₁₁₁₁, S₁₁₂₂, S₁₁₃₃, S₃₃₃₃, S₂₃₂₃ = argTI(𝕊)
        @test S₁₁₁₁ ≈ 1/E₁          atol=atol_num
        @test S₃₃₃₃ ≈ 1/E₃          atol=atol_num
        @test S₁₁₂₂ ≈ -ν₁₂/E₁      atol=atol_num
        @test S₁₁₃₃ ≈ -ν₃₁/E₃      atol=atol_num
        @test S₂₃₂₃ ≈ 1/(4*G₃₁)    atol=atol_num

        # inv(𝕊) ⊡ 𝕊 = 𝕀
        ℂ = inv(𝕊)
        𝕀 = tensId4(Val(3), Val(Float64))
        prod = ℂ ⊡ 𝕊
        @test opequal(getarray(prod), getarray(𝕀))
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "tensTI & argTI — symbolic round-trip" begin
        C₁, C₂, C₃, C₄, C₅ = symbols("C₁₁₁₁ C₁₁₂₂ C₁₁₃₃ C₃₃₃₃ C₂₃₂₃", real = true)
        ℂ = tensTI(C₁, C₂, C₃, C₄, C₅, n3s)
        @test ℂ isa TensWalpole{<:Any,5}
        c₁, c₂, c₃, c₄, c₅ = argTI(ℂ)
        @test simplify(c₁ - C₁) == 0
        @test simplify(c₂ - C₂) == 0
        @test simplify(c₃ - C₃) == 0
        @test simplify(c₄ - C₄) == 0
        @test simplify(c₅ - C₅) == 0
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "tensTI_Hoenig & argTI_Hoenig — numeric round-trip" begin
        E, ν₁, ν₂, H, Γ = 100.0, 0.25, 0.15, 2.0, 3.0
        𝕊 = tensTI_Hoenig(E, ν₁, ν₂, H, Γ, n3)
        @test 𝕊 isa TensWalpole{Float64,5}

        e, n1, n2, h, g = argTI_Hoenig(𝕊)
        @test e  ≈ E   atol=atol_num
        @test n1 ≈ ν₁  atol=atol_num
        @test n2 ≈ ν₂  atol=atol_num
        @test h  ≈ H   atol=atol_num
        @test g  ≈ Γ   atol=atol_num
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "tensTI_Hoenig — consistency with tensTI_eng" begin
        E, ν₁, ν₂, H, Γ = 100.0, 0.25, 0.15, 2.0, 3.0
        𝕊h = tensTI_Hoenig(E, ν₁, ν₂, H, Γ, n3)

        # Convert Hoenig → standard engineering constants
        E₁  = E
        E₃  = E * H
        ν₁₂ = ν₁
        ν₃₁ = H * ν₂
        G₃₁ = E * Γ / (2 * (1 + ν₁))
        𝕊e = tensTI_eng(E₁, E₃, ν₁₂, ν₃₁, G₃₁, n3)

        @test opequal(getarray(𝕊h), getarray(𝕊e))
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "tensTI_Hoenig — isotropic limit" begin
        # Setting ν₁=ν₂=ν, H=Γ=1 should give the isotropic compliance tensor
        E, ν = 100.0, 0.3
        𝕊 = tensTI_Hoenig(E, ν, ν, 1.0, 1.0, n3)
        k = E / (3 * (1 - 2ν))
        μ = E / (2 * (1 + ν))
        𝕁, 𝕂 = tensJ4(Val(3), Val(Float64)), tensK4(Val(3), Val(Float64))
        𝕊iso = (1 / (3k)) * 𝕁 + (1 / (2μ)) * 𝕂
        @test opequal(getarray(𝕊), getarray(𝕊iso))
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "otimes — TensTI{2} self-product" begin
        n = [0., 0., 1.]
        A = TensTI{2}(5.0, 8.0, n)
        R = otimes(A)
        @test R isa TensWalpole{Float64,5}
        sq2 = sqrt(2.0)
        @test getdata(R) == (64.0, 50.0, sq2*40.0, 0.0, 0.0)
        @test opequal(getarray(R), otimes(getarray(A), getarray(A)))
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "otimes — TensTI{2} × TensTI{2} (same axis)" begin
        n = [0., 0., 1.]
        A = TensTI{2}(5.0, 8.0, n)
        B = TensTI{2}(3.0, 2.0, n)
        R = A ⊗ B
        @test R isa TensWalpole{Float64,6}
        sq2 = sqrt(2.0)
        # ℓ₁=b₁b₂, ℓ₂=2a₁a₂, ℓ₃=√2·b₁a₂, ℓ₄=√2·a₁b₂
        @test R.data[1] ≈ 16.0     atol=atol_num   # 8*2
        @test R.data[2] ≈ 30.0     atol=atol_num   # 2*5*3
        @test R.data[3] ≈ sq2*24.0 atol=atol_num   # √2*8*3
        @test R.data[4] ≈ sq2*10.0 atol=atol_num   # √2*5*2
        @test R.data[5] ≈ 0.0      atol=atol_num
        @test R.data[6] ≈ 0.0      atol=atol_num
        @test opequal(getarray(R), otimes(getarray(A), getarray(B)))
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "otimes — TensISO{2,3} × TensTI{2} and reverse" begin
        n = [0., 0., 1.]
        I2 = tensId2(Val(3), Val(Float64))   # TensISO{2,3}(1.0)
        B = TensTI{2}(5.0, 8.0, n)
        R1 = I2 ⊗ B
        @test R1 isa TensWalpole{Float64,6}
        @test opequal(getarray(R1), otimes(getarray(I2), getarray(B)))
        R2 = B ⊗ I2
        @test R2 isa TensWalpole{Float64,6}
        @test opequal(getarray(R2), otimes(getarray(B), getarray(I2)))
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "otimes — TensTI{2} isotropic limit" begin
        n = [0., 0., 1.]
        λ = 3.0
        A = TensTI{2}(λ, λ, n)        # isotropic: a == b
        Rw = otimes(A)                  # → TensWalpole{T,5}
        Riso = otimes(TensISO{3}(λ))   # → TensISO{4,3}
        @test opequal(getarray(Rw), getarray(Riso))
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "otimes — TensTI{2} non-canonical axis" begin
        n = [1.0, 1.0, 1.0] / sqrt(3.0)
        A = TensTI{2}(4.0, 7.0, n)
        B = TensTI{2}(2.0, 5.0, n)
        R = A ⊗ B
        @test R isa TensWalpole{Float64,6}
        @test opequal(getarray(R), otimes(getarray(A), getarray(B)))
        Rs = otimes(A)
        @test Rs isa TensWalpole{Float64,5}
        @test opequal(getarray(Rs), otimes(getarray(A), getarray(A)))
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "otimes — TensTI{2} different axes (fallback)" begin
        n1 = [0., 0., 1.]
        n2 = [1., 0., 0.]
        A = TensTI{2}(5.0, 8.0, n1)
        B = TensTI{2}(3.0, 2.0, n2)
        R = A ⊗ B
        # Different axes → generic fallback → Tens (not TensWalpole)
        @test !(R isa TensWalpole)
        @test opequal(getarray(R), otimes(getarray(A), getarray(B)))
    end

end  # "Walpole & Ortho tensors"
