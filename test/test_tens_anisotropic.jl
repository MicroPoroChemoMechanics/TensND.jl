import Symbolics
# ─── helpers shared across sub-sections ───────────────────────────────────
# The sections are top-level blocks rather than the body of one enclosing
# `@testsection`: a top-level block holding loops is compiled as one function,
# and the whole file as one block took more than 3 GB to compile, against
# 0.93 GB section by section.
const _WO_N3 = 𝐞(Val(3), Val(3), Val(Float64))   # e₃ as Float64 Vec
const _WO_N3S = 𝐞(Val(3), Val(3), Val(Sym))       # e₃ as Sym Vec
const _WO_ATOL = 1.0e-12

# ═══════════════════════════════════════════════════════════════════════════
@testsection "Walpole & Ortho tensors: TensTI{4} — construction & traits" begin
    W1, W2, W3, W4, W5, W6 = walpole_basis(_WO_N3)
    @test W1 isa TensTI{4, Float64, 6}
    @test W2 isa TensTI{4, Float64, 6}
    @test size(W1) == (3, 3, 3, 3)
    @test get_basis(W1) isa CanonicalBasis{3, Float64}
    @test get_var(W1) == (:cont, :cont, :cont, :cont)

    # N=5 (symmetric) basis
    W1s, W2s, W3s, W4s, W5s = walpole_basis_sym(_WO_N3)
    @test W1s isa TensTI{4, Float64, 5}
    @test size(W1s) == (3, 3, 3, 3)

    # individual constructors
    @test tens_W1(_WO_N3) isa TensTI{4, Float64, 6}
    @test tens_W3(_WO_N3) isa TensTI{4, Float64, 6}

    # axis accessor
    @test axis(W1) == (0.0, 0.0, 1.0)
end

# ═══════════════════════════════════════════════════════════════════════════
@testsection "Walpole & Ortho tensors: TensTI{4} — get_array & KM (n=e₃, Float64)" begin
    W1, W2, W3, W4, W5, W6 = walpole_basis(_WO_N3)
    sq2 = sqrt(2.0)

    # W₁ = nₙ⊗nₙ  →  only entry [3,3,3,3] = 1
    A1 = get_array(W1)
    @test A1[3, 3, 3, 3] ≈ 1.0   atol = _WO_ATOL
    @test all(
        abs.(A1[i, j, k, l]) < _WO_ATOL
            for i in 1:3, j in 1:3, k in 1:3, l in 1:3
            if !(i == 3 && j == 3 && k == 3 && l == 3)
    )

    # W₂ = (nT⊗nT)/2 with n=e₃ → nT = diag(1,1,0)
    # so W₂[i,j,k,l] = δᵢⱼ(1-δᵢ₃)δₖₗ(1-δₖ₃)/2
    A2 = get_array(W2)
    @test A2[1, 1, 1, 1] ≈ 0.5   atol = _WO_ATOL
    @test A2[1, 1, 2, 2] ≈ 0.5   atol = _WO_ATOL
    @test A2[2, 2, 2, 2] ≈ 0.5   atol = _WO_ATOL
    @test A2[3, 3, 3, 3] ≈ 0.0   atol = _WO_ATOL

    # W₃[3,3,1,1] = 1/√2
    A3 = get_array(W3)
    @test A3[3, 3, 1, 1] ≈ 1.0 / sq2   atol = _WO_ATOL
    @test A3[3, 3, 2, 2] ≈ 1.0 / sq2   atol = _WO_ATOL
    @test A3[3, 3, 3, 3] ≈ 0.0       atol = _WO_ATOL

    # W₄[1,1,3,3] = 1/√2
    A4 = get_array(W4)
    @test A4[1, 1, 3, 3] ≈ 1.0 / sq2   atol = _WO_ATOL

    # W₅: shear in transverse plane (e₁,e₂)
    # W₅[1,2,1,2] = W₅[2,1,1,2] = W₅[1,2,2,1] = W₅[2,1,2,1] = 1/2
    A5 = get_array(W5)
    @test A5[1, 2, 1, 2] ≈ 0.5   atol = _WO_ATOL
    @test A5[1, 1, 1, 1] ≈ 0.5   atol = _WO_ATOL
    @test A5[3, 3, 3, 3] ≈ 0.0   atol = _WO_ATOL

    # W₆: shear between transverse and axial
    A6 = get_array(W6)
    @test A6[1, 3, 1, 3] ≈ 0.5   atol = _WO_ATOL
    @test A6[2, 3, 2, 3] ≈ 0.5   atol = _WO_ATOL
    @test A6[1, 1, 1, 1] ≈ 0.0   atol = _WO_ATOL
    @test A6[3, 3, 3, 3] ≈ 0.0   atol = _WO_ATOL

    # KM structure for n=e₃: blocks should separate
    L = TensTI{4}(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, _WO_N3)
    Km = mandel(L)
    @test size(Km) == (6, 6)
    # Off-diagonal shear coupling should be zero for axis n=e₃
    @test abs(Km[1, 4]) < _WO_ATOL
    @test abs(Km[1, 5]) < _WO_ATOL
    @test abs(Km[4, 6]) < _WO_ATOL
end

# ═══════════════════════════════════════════════════════════════════════════
@testsection "Walpole & Ortho tensors: TensTI{4} — Walpole product rule (Float64)" begin
    W1, W2, W3, W4, W5, W6 = walpole_basis(_WO_N3)

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
    L = TensTI{4}(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, _WO_N3)
    M = TensTI{4}(0.5, 1.5, 2.0, 0.3, 0.8, 1.2, _WO_N3)
    LM_walpole = get_array(L ⊡ M)
    LM_direct = Tensor{4, 3}(get_array(L)) ⊡ Tensor{4, 3}(get_array(M))
    @test opequal(LM_walpole, Array(LM_direct))
end

# ═══════════════════════════════════════════════════════════════════════════
@testsection "Walpole & Ortho tensors: TensTI{4} — inverse (Float64)" begin
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
@testsection "Walpole & Ortho tensors: TensTI{4} — iso_to_ti" begin
    𝕀, 𝕁, 𝕂 = iso_projectors(Val(3), Val(Float64))
    α, β = 3.0, 2.0
    ℂiso = α * 𝕁 + β * 𝕂

    # For any axis, iso_to_ti should give the same array as the isotropic tensor
    for n in (
            [0.0, 0.0, 1.0], [1.0, 0.0, 0.0],
            [1.0 / sqrt(3), 1.0 / sqrt(3), 1.0 / sqrt(3)],
        )
        Wiso = iso_to_ti(ℂiso, n)
        @test Wiso isa TensTI{4, Float64, 5}
        @test istransverselyisotropic(Wiso)
        @test opequal(get_array(Wiso), get_array(ℂiso))
    end

    # Symbolic — use Sym ISO tensors to avoid Float64/Sym residuals
    αs, βs = symbols("α β", real = true)
    𝕁s, 𝕂s = tens_J4(Val(3), Val(Sym)), tens_K4(Val(3), Val(Sym))
    ℂisos = αs * 𝕁s + βs * 𝕂s
    Wisos = iso_to_ti(ℂisos, _WO_N3S)
    @test Wisos isa TensTI{4, <:Any, 5}
    @test all(==(0), simplify.(get_array(Wisos) .- get_array(ℂisos)))
end

# ═══════════════════════════════════════════════════════════════════════════
@testsection "Walpole & Ortho tensors: TensTI{4} — isisotropic / istransverselyisotropic / isorthotropic" begin
    W = tens_W1(_WO_N3)
    @test !isisotropic(W)
    @test  istransverselyisotropic(W)
    @test !isorthotropic(W)
    𝕀, 𝕁, 𝕂 = iso_projectors(Val(3), Val(Float64))
    @test !istransverselyisotropic(𝕁)
    @test !isorthotropic(𝕁)
end

# ═══════════════════════════════════════════════════════════════════════════
@testsection "Walpole & Ortho tensors: TensOrtho — isisotropic / istransverselyisotropic / isorthotropic" begin
    frame3 = CanonicalBasis{3, Float64}()
    t = TensOrtho(10.0, 8.0, 9.0, 3.0, 2.0, 4.0, 2.5, 3.0, 1.5, frame3)
    @test !isisotropic(t)
    @test !istransverselyisotropic(t)
    @test  isorthotropic(t)
    # universal fallback
    @test !isorthotropic(42)
    @test !isorthotropic("string")
end

# ═══════════════════════════════════════════════════════════════════════════
@testsection "Walpole & Ortho tensors: TensTI{4} — show" begin
    L5 = TensTI{4}(1.0, 2.0, 0.5, 3.0, 4.0, _WO_N3)   # N=5
    L6 = TensTI{4}(1.0, 2.0, 0.5, 0.3, 3.0, 4.0, _WO_N3)   # N=6
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
@testsection "Walpole & Ortho tensors: TensOrtho — show" begin
    frame3 = CanonicalBasis{3, Float64}()
    t = TensOrtho(10.0, 8.0, 9.0, 3.0, 2.0, 4.0, 2.5, 3.0, 1.5, frame3)
    buf = IOBuffer()
    show(buf, t)
    s = String(take!(buf))
    @test contains(s, "P₁⊗P₁") && contains(s, "frame")
end

# ═══════════════════════════════════════════════════════════════════════════
@testsection "Walpole & Ortho tensors: TensTI{4} — tsimplify (symbolic, _rebuild)" begin
    ℓ₁, ℓ₂, ℓ₃ = symbols("ℓ₁ ℓ₂ ℓ₃", real = true)
    L = TensTI{4}(ℓ₁, ℓ₂, ℓ₃, ℓ₁ + ℓ₂, ℓ₂ + ℓ₃, _WO_N3S)   # N=6
    Ls = tsimplify(L)
    @test Ls isa TensTI{4}   # _rebuild preserves type
    @test walpole_coeffs(Ls)[1] == ℓ₁   # simplification is a no-op here
end

# ═══════════════════════════════════════════════════════════════════════════
@testsection "Walpole & Ortho tensors: TensOrtho — tsimplify (_rebuild)" begin
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
@testsection "Walpole & Ortho tensors: TensTI{4} — arithmetic" begin
    W1, W2, W3, W4, W5, W6 = walpole_basis(_WO_N3)
    L = TensTI{4}(1.0, 2.0, 0.5, 0.5, 3.0, 4.0, _WO_N3)   # N=6
    M = TensTI{4}(0.5, 1.0, 0.25, 0.25, 1.5, 2.0, _WO_N3)

    @test opequal(get_array(L + M), get_array(L) .+ get_array(M))
    @test opequal(get_array(L - M), get_array(L) .- get_array(M))
    @test opequal(get_array(2.0 * L), 2.0 .* get_array(L))
    @test opequal(get_array(-L), .-get_array(L))

    # Symmetric N=5
    Ls = TensTI{4}(1.0, 2.0, 0.5, 3.0, 4.0, _WO_N3)
    Ms = TensTI{4}(0.5, 1.0, 0.25, 1.5, 2.0, _WO_N3)
    @test (Ls + Ms) isa TensTI{4, Float64, 5}
    @test opequal(get_array(Ls + Ms), get_array(Ls) .+ get_array(Ms))
end

# ═══════════════════════════════════════════════════════════════════════════
@testsection "Walpole & Ortho tensors: TensTI{4} — symbolic inverse" begin
    ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆ = symbols("ℓ₁ ℓ₂ ℓ₃ ℓ₅ ℓ₆", real = true)
    L = TensTI{4}(ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆, _WO_N3S)
    Li = inv(L)
    @test Li isa TensTI{4, <:Any, 5}
    # L⊡inv(L) should be identity
    𝕀 = tens_Id4(Val(3), Val(Sym))
    prod = L ⊡ Li
    @test all(==(0), simplify.(get_array(prod) .- get_array(𝕀)))
end

# ═══════════════════════════════════════════════════════════════════════════
@testsection "Walpole & Ortho tensors: TensOrtho — construction & traits" begin
    frame3 = CanonicalBasis{3, Float64}()
    C11, C22, C33 = 10.0, 8.0, 12.0
    C12, C13, C23 = 3.0, 4.0, 2.5
    C44, C55, C66 = 2.0, 3.0, 1.5
    t = TensOrtho(C11, C22, C33, C12, C13, C23, C44, C55, C66, frame3)
    @test t isa TensOrtho{4, 3, Float64}
    @test size(t) == (3, 3, 3, 3)
    @test get_basis(t) isa CanonicalBasis{3, Float64}
    @test get_var(t) == (:cont, :cont, :cont, :cont)
    @test frame(t) === frame3
end

# ═══════════════════════════════════════════════════════════════════════════
@testsection "Walpole & Ortho tensors: TensOrtho — mandel_material (canonical frame)" begin
    frame3 = CanonicalBasis{3, Float64}()
    C11, C22, C33 = 10.0, 8.0, 12.0
    C12, C13, C23 = 3.0, 4.0, 2.5
    C44, C55, C66 = 2.0, 3.0, 1.5
    t = TensOrtho(C11, C22, C33, C12, C13, C23, C44, C55, C66, frame3)
    Km = mandel_material(t)
    @test size(Km) == (6, 6)
    # Diagonal blocks
    @test Km[1, 1] ≈ C11  atol = _WO_ATOL
    @test Km[2, 2] ≈ C22  atol = _WO_ATOL
    @test Km[3, 3] ≈ C33  atol = _WO_ATOL
    @test Km[4, 4] ≈ 2 * C44 atol = _WO_ATOL
    @test Km[5, 5] ≈ 2 * C55 atol = _WO_ATOL
    @test Km[6, 6] ≈ 2 * C66 atol = _WO_ATOL
    # Off-diagonal within normal block
    @test Km[1, 2] ≈ C12  atol = _WO_ATOL
    @test Km[1, 3] ≈ C13  atol = _WO_ATOL
    @test Km[2, 3] ≈ C23  atol = _WO_ATOL
    # Zeros between normal and shear blocks
    @test Km[1, 4] ≈ 0.0  atol = _WO_ATOL
    @test Km[2, 5] ≈ 0.0  atol = _WO_ATOL
    @test Km[3, 6] ≈ 0.0  atol = _WO_ATOL
    @test Km[4, 5] ≈ 0.0  atol = _WO_ATOL
    @test Km[4, 6] ≈ 0.0  atol = _WO_ATOL
    @test Km[5, 6] ≈ 0.0  atol = _WO_ATOL
    # Symmetry
    @test Km ≈ Km'  atol = _WO_ATOL
    # KM in canonical frame should match
    Km2 = mandel(t)
    @test Km2 ≈ Km  atol = _WO_ATOL
end

# ═══════════════════════════════════════════════════════════════════════════
@testsection "Walpole & Ortho tensors: TensOrtho — from KM matrix" begin
    frame3 = CanonicalBasis{3, Float64}()
    C11, C22, C33 = 10.0, 8.0, 12.0
    C12, C13, C23 = 3.0, 4.0, 2.5
    C44, C55, C66 = 2.0, 3.0, 1.5
    t1 = TensOrtho(C11, C22, C33, C12, C13, C23, C44, C55, C66, frame3)
    Km = mandel_material(t1)
    t2 = TensOrtho(Km, frame3)
    @test opequal(get_array(t1), get_array(t2))
end

# ═══════════════════════════════════════════════════════════════════════════
@testsection "Walpole & Ortho tensors: TensOrtho — inverse (canonical frame)" begin
    frame3 = CanonicalBasis{3, Float64}()
    C11, C22, C33 = 10.0, 8.0, 12.0
    C12, C13, C23 = 3.0, 4.0, 2.5
    C44, C55, C66 = 2.0, 3.0, 1.5
    t = TensOrtho(C11, C22, C33, C12, C13, C23, C44, C55, C66, frame3)
    ti = inv(t)
    @test ti isa TensOrtho{4, 3, Float64}
    𝕀 = tens_Id4(Val(3), Val(Float64))
    A = Tensor{4, 3}(get_array(t))
    Ai = Tensor{4, 3}(get_array(ti))
    prod = Array(A ⊡ Ai)
    𝕀arr = get_array(𝕀)
    @test opequal(prod, 𝕀arr)
end

# ═══════════════════════════════════════════════════════════════════════════
@testsection "Walpole & Ortho tensors: TensOrtho — rotated frame" begin
    # Rotate material frame: e₁→e₂, e₂→e₃, e₃→e₁ (cyclic permutation)
    R = Float64[0 0 1; 1 0 0; 0 1 0]
    frame_rot = RotatedBasis(R)
    C11, C22, C33 = 10.0, 10.0, 12.0
    C12, C13, C23 = 3.0, 3.0, 3.0
    C44, C55, C66 = 2.0, 2.0, 2.0
    t = TensOrtho(C11, C22, C33, C12, C13, C23, C44, C55, C66, frame_rot)
    # mandel_material should still have block-diagonal structure
    Km = mandel_material(t)
    @test Km[1, 4] ≈ 0.0  atol = _WO_ATOL
    @test Km[4, 5] ≈ 0.0  atol = _WO_ATOL
    @test Km[5, 6] ≈ 0.0  atol = _WO_ATOL
    # express to canonical basis should give a valid Tens
    tc = express(t, CanonicalBasis{3, Float64}())
    @test tc isa AbstractTens
end

# ═══════════════════════════════════════════════════════════════════════════
@testsection "Walpole & Ortho tensors: TensOrtho — TI consistency" begin
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
    tw = TensTI{4}(ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆, _WO_N3)
    @test opequal(get_array(to), get_array(tw))
end

# ═══════════════════════════════════════════════════════════════════════════
@testsection "Walpole & Ortho tensors: TensTI{4} — dcontract with TensISO" begin
    𝕀, 𝕁, 𝕂 = iso_projectors(Val(3), Val(Float64))
    α, β = 3.0, 2.0
    ℂiso = α * 𝕁 + β * 𝕂
    L = TensTI{4}(2.0, 3.0, 1.0, 4.0, 5.0, _WO_N3)   # N=5

    # L⊡ℂiso via Walpole product rule (convert ISO first)
    res_w = L ⊡ ℂiso
    @test res_w isa TensTI{4}
    # Compare with direct array contraction
    res_direct = Tensor{4, 3}(get_array(L)) ⊡ Tensor{4, 3}(get_array(ℂiso))
    @test opequal(get_array(res_w), Array(res_direct))
end

# ═══════════════════════════════════════════════════════════════════════════
@testsection "Walpole & Ortho tensors: ti_from_constants & ti_constants — numeric round-trip" begin
    C₁₁₁₁, C₁₁₂₂, C₁₁₃₃, C₃₃₃₃, C₂₃₂₃ = 10.0, 3.0, 2.5, 12.0, 2.0
    ℂ = ti_from_constants(C₁₁₁₁, C₁₁₂₂, C₁₁₃₃, C₃₃₃₃, C₂₃₂₃, _WO_N3)
    @test ℂ isa TensTI{4, Float64, 5}

    c₁₁₁₁, c₁₁₂₂, c₁₁₃₃, c₃₃₃₃, c₂₃₂₃ = ti_constants(ℂ)
    @test c₁₁₁₁ ≈ C₁₁₁₁  atol = _WO_ATOL
    @test c₁₁₂₂ ≈ C₁₁₂₂  atol = _WO_ATOL
    @test c₁₁₃₃ ≈ C₁₁₃₃  atol = _WO_ATOL
    @test c₃₃₃₃ ≈ C₃₃₃₃  atol = _WO_ATOL
    @test c₂₃₂₃ ≈ C₂₃₂₃  atol = _WO_ATOL
end

# ═══════════════════════════════════════════════════════════════════════════
@testsection "Walpole & Ortho tensors: ti_from_constants — consistency with TensOrtho" begin
    frame3 = CanonicalBasis{3, Float64}()
    C₁₁₁₁ = 10.0; C₃₃₃₃ = 12.0; C₁₁₂₂ = 3.0; C₁₁₃₃ = 2.5; C₂₃₂₃ = 2.0
    C₆₆ = (C₁₁₁₁ - C₁₁₂₂) / 2
    to = TensOrtho(
        C₁₁₁₁, C₁₁₁₁, C₃₃₃₃, C₁₁₂₂, C₁₁₃₃, C₁₁₃₃,
        C₂₃₂₃, C₂₃₂₃, C₆₆, frame3
    )
    tw = ti_from_constants(C₁₁₁₁, C₁₁₂₂, C₁₁₃₃, C₃₃₃₃, C₂₃₂₃, _WO_N3)
    @test opequal(get_array(to), get_array(tw))
end

# ═══════════════════════════════════════════════════════════════════════════
@testsection "Walpole & Ortho tensors: ti_from_constants & ti_constants — numeric round-trip" begin
    E₁, E₃, ν₁₂, ν₃₁, G₃₁ = 100.0, 200.0, 0.25, 0.15, 40.0
    𝕊 = ti_from_constants(E₁, E₃, ν₁₂, ν₃₁, G₃₁, _WO_N3; form = :engineering)
    @test 𝕊 isa TensTI{4, Float64, 5}

    e₁, e₃, n₁₂, n₃₁, g₃₁ = ti_constants(𝕊; form = :engineering)
    @test e₁ ≈ E₁   atol = _WO_ATOL
    @test e₃ ≈ E₃   atol = _WO_ATOL
    @test n₁₂ ≈ ν₁₂  atol = _WO_ATOL
    @test n₃₁ ≈ ν₃₁  atol = _WO_ATOL
    @test g₃₁ ≈ G₃₁  atol = _WO_ATOL
end

# ═══════════════════════════════════════════════════════════════════════════
@testsection "Walpole & Ortho tensors: ti_from_constants — compliance components check" begin
    E₁, E₃, ν₁₂, ν₃₁, G₃₁ = 100.0, 200.0, 0.25, 0.15, 40.0
    𝕊 = ti_from_constants(E₁, E₃, ν₁₂, ν₃₁, G₃₁, _WO_N3; form = :engineering)

    S₁₁₁₁, S₁₁₂₂, S₁₁₃₃, S₃₃₃₃, S₂₃₂₃ = ti_constants(𝕊)
    @test S₁₁₁₁ ≈ 1 / E₁          atol = _WO_ATOL
    @test S₃₃₃₃ ≈ 1 / E₃          atol = _WO_ATOL
    @test S₁₁₂₂ ≈ -ν₁₂ / E₁      atol = _WO_ATOL
    @test S₁₁₃₃ ≈ -ν₃₁ / E₃      atol = _WO_ATOL
    @test S₂₃₂₃ ≈ 1 / (4 * G₃₁)    atol = _WO_ATOL

    # inv(𝕊) ⊡ 𝕊 = 𝕀
    ℂ = inv(𝕊)
    𝕀 = tens_Id4(Val(3), Val(Float64))
    prod = ℂ ⊡ 𝕊
    @test opequal(get_array(prod), get_array(𝕀))
end

# ═══════════════════════════════════════════════════════════════════════════
@testsection "Walpole & Ortho tensors: ti_from_constants & ti_constants — symbolic round-trip" begin
    C₁, C₂, C₃, C₄, C₅ = symbols("C₁₁₁₁ C₁₁₂₂ C₁₁₃₃ C₃₃₃₃ C₂₃₂₃", real = true)
    ℂ = ti_from_constants(C₁, C₂, C₃, C₄, C₅, _WO_N3S)
    @test ℂ isa TensTI{4, <:Any, 5}
    c₁, c₂, c₃, c₄, c₅ = ti_constants(ℂ)
    @test simplify(c₁ - C₁) == 0
    @test simplify(c₂ - C₂) == 0
    @test simplify(c₃ - C₃) == 0
    @test simplify(c₄ - C₄) == 0
    @test simplify(c₅ - C₅) == 0
end

# ═══════════════════════════════════════════════════════════════════════════
@testsection "Walpole & Ortho tensors: ti_from_constants & ti_constants — numeric round-trip" begin
    E, ν₁, ν₂, H, Γ = 100.0, 0.25, 0.15, 2.0, 3.0
    𝕊 = ti_from_constants(E, ν₁, ν₂, H, Γ, _WO_N3; form = :hoenig)
    @test 𝕊 isa TensTI{4, Float64, 5}

    # A form that does not exist is refused, by name, in both directions.
    @test_throws "unknown form `:hoening`" ti_constants(𝕊; form = :hoening)
    @test_throws ArgumentError ti_from_constants(E, ν₁, ν₂, H, Γ, _WO_N3; form = :voigt)
    e, n1, n2, h, g = ti_constants(𝕊; form = :hoenig)
    @test e ≈ E   atol = _WO_ATOL
    @test n1 ≈ ν₁  atol = _WO_ATOL
    @test n2 ≈ ν₂  atol = _WO_ATOL
    @test h ≈ H   atol = _WO_ATOL
    @test g ≈ Γ   atol = _WO_ATOL
end

# ═══════════════════════════════════════════════════════════════════════════
@testsection "Walpole & Ortho tensors: ti_from_constants — consistency with ti_from_constants" begin
    E, ν₁, ν₂, H, Γ = 100.0, 0.25, 0.15, 2.0, 3.0
    𝕊h = ti_from_constants(E, ν₁, ν₂, H, Γ, _WO_N3; form = :hoenig)

    # Convert Hoenig → standard engineering constants
    E₁ = E
    E₃ = E * H
    ν₁₂ = ν₁
    ν₃₁ = H * ν₂
    G₃₁ = E * Γ / (2 * (1 + ν₁))
    𝕊e = ti_from_constants(E₁, E₃, ν₁₂, ν₃₁, G₃₁, _WO_N3; form = :engineering)

    @test opequal(get_array(𝕊h), get_array(𝕊e))
end

# ═══════════════════════════════════════════════════════════════════════════
@testsection "Walpole & Ortho tensors: ti_from_constants — isotropic limit" begin
    # Setting ν₁=ν₂=ν, H=Γ=1 should give the isotropic compliance tensor
    E, ν = 100.0, 0.3
    𝕊 = ti_from_constants(E, ν, ν, 1.0, 1.0, _WO_N3; form = :hoenig)
    k = E / (3 * (1 - 2ν))
    μ = E / (2 * (1 + ν))
    𝕁, 𝕂 = tens_J4(Val(3), Val(Float64)), tens_K4(Val(3), Val(Float64))
    𝕊iso = (1 / (3k)) * 𝕁 + (1 / (2μ)) * 𝕂
    @test opequal(get_array(𝕊), get_array(𝕊iso))
end

# ═══════════════════════════════════════════════════════════════════════════
@testsection "Walpole & Ortho tensors: otimes — TensTI{2} self-product" begin
    n = [0.0, 0.0, 1.0]
    A = TensTI{2}(5.0, 8.0, n)
    R = otimes(A)
    @test R isa TensTI{4, Float64, 5}
    sq2 = sqrt(2.0)
    @test get_data(R) == (64.0, 50.0, sq2 * 40.0, 0.0, 0.0)
    @test opequal(get_array(R), otimes(get_array(A), get_array(A)))
end

# ═══════════════════════════════════════════════════════════════════════════
@testsection "Walpole & Ortho tensors: otimes — TensTI{2} × TensTI{2} (same axis)" begin
    n = [0.0, 0.0, 1.0]
    A = TensTI{2}(5.0, 8.0, n)
    B = TensTI{2}(3.0, 2.0, n)
    R = A ⊗ B
    @test R isa TensTI{4, Float64, 6}
    sq2 = sqrt(2.0)
    # ℓ₁=b₁b₂, ℓ₂=2a₁a₂, ℓ₃=√2·b₁a₂, ℓ₄=√2·a₁b₂
    @test R.data[1] ≈ 16.0     atol = _WO_ATOL   # 8*2
    @test R.data[2] ≈ 30.0     atol = _WO_ATOL   # 2*5*3
    @test R.data[3] ≈ sq2 * 24.0 atol = _WO_ATOL   # √2*8*3
    @test R.data[4] ≈ sq2 * 10.0 atol = _WO_ATOL   # √2*5*2
    @test R.data[5] ≈ 0.0      atol = _WO_ATOL
    @test R.data[6] ≈ 0.0      atol = _WO_ATOL
    @test opequal(get_array(R), otimes(get_array(A), get_array(B)))
end

# ═══════════════════════════════════════════════════════════════════════════
@testsection "Walpole & Ortho tensors: otimes — TensISO{2,3} × TensTI{2} and reverse" begin
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
@testsection "Walpole & Ortho tensors: otimes — TensTI{2} isotropic limit" begin
    n = [0.0, 0.0, 1.0]
    λ = 3.0
    A = TensTI{2}(λ, λ, n)        # isotropic: a == b
    Rw = otimes(A)                  # → TensTI{4, T, 5}
    Riso = otimes(TensISO{3}(λ))   # → TensISO{4,3}
    @test opequal(get_array(Rw), get_array(Riso))
end

# ═══════════════════════════════════════════════════════════════════════════
@testsection "Walpole & Ortho tensors: otimes — TensTI{2} non-canonical axis" begin
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
@testsection "Walpole & Ortho tensors: otimes — TensTI{2} different axes (fallback)" begin
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
@testsection "Walpole & Ortho tensors: Cross-type dispatch — promotion output types" begin
    frame = CanonicalBasis{3, Float64}()
    n = (0.0, 0.0, 1.0)

    I4 = TensISO{3}(2.0, 3.0)   # α=2, β=3
    W = ti_from_constants(10.0, 3.0, 2.5, 12.0, 2.0, n)                 # TensTI{4, T, 5}
    O = TensOrtho(10.0, 8.0, 12.0, 3.0, 2.5, 1.5, 2.0, 3.0, 3.5, frame)

    # ── + / − between ISO and Ortho ────────────────────────────────────
    @test I4 + O isa TensOrtho{4, 3, Float64}
    @test O + I4 isa TensOrtho{4, 3, Float64}
    @test I4 - O isa TensOrtho{4, 3, Float64}
    @test O - I4 isa TensOrtho{4, 3, Float64}

    # Numerical agreement with generic array arithmetic
    @test opequal(get_array(I4 + O), get_array(I4) + get_array(O))
    @test opequal(get_array(O - I4), get_array(O) - get_array(I4))

    # ── + / − between Walpole(N=5) and Ortho (aligned axis) ─────────────
    @test W + O isa TensOrtho{4, 3, Float64}
    @test O - W isa TensOrtho{4, 3, Float64}
    @test opequal(get_array(W + O), get_array(W) + get_array(O))
    @test opequal(get_array(O - W), get_array(O) - get_array(W))

    # Reference mismatch → the generic route, the sum of the canonical arrays
    n_mis = (1.0, 1.0, 1.0) ./ sqrt(3.0)
    Wmis = ti_from_constants(10.0, 3.0, 2.5, 12.0, 2.0, n_mis)
    @test !((Wmis + O) isa TensOrtho) && opequal(canonical_components(Wmis + O), canonical_components(Wmis) + canonical_components(O))
    @test opequal(canonical_components(O - Wmis), canonical_components(O) - canonical_components(Wmis))

    # ── TensTI{4} N=5 + N=6 (same axis) → N=6 ─────────────────────────
    W5 = ti_from_constants(10.0, 3.0, 2.5, 12.0, 2.0, n)          # N=5
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
@testsection "Walpole & Ortho tensors: Cross-type dispatch — cross-order double contraction" begin
    n = (0.0, 0.0, 1.0)
    I4 = TensISO{3}(2.0, 3.0)
    W = ti_from_constants(10.0, 3.0, 2.5, 12.0, 2.0, n)
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
    Wm = ti_from_constants(10.0, 3.0, 2.5, 12.0, 2.0, n2)
    T2m = TensTI{2}(4.0, 7.0, n)
    Rm = Wm ⊡ T2m
    @test !(Rm isa TensTI)
    @test opequal(Array(get_array(Rm)), dcontract(get_array(Wm), get_array(T2m)))
end

# ═══════════════════════════════════════════════════════════════════════════
@testsection "Walpole & Ortho tensors: Cross-type dispatch — single contraction (·)" begin
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
@testsection "Walpole & Ortho tensors: Cross-type dispatch — inverse consistency" begin
    n = (0.0, 0.0, 1.0)
    frame = CanonicalBasis{3, Float64}()
    Id4 = tens_Id4(Val(3), Val(Float64))

    # TensTI{4}
    W = ti_from_constants(100.0, 30.0, 25.0, 120.0, 20.0, n)
    @test opequal(get_array(inv(W) ⊡ W), get_array(Id4))

    # TensOrtho
    O = TensOrtho(10.0, 8.0, 12.0, 3.0, 2.5, 1.5, 2.0, 3.0, 3.5, frame)
    @test opequal(get_array(inv(O)) ⊡ get_array(O), get_array(Id4))

    # TensISO{4}
    I4 = TensISO{3}(2.0, 3.0)
    @test opequal(get_array(inv(I4) ⊡ I4), get_array(Id4))
end

# ═══════════════════════════════════════════════════════════════════════════
@testsection "Walpole & Ortho tensors: Promotion helpers — iso_to_ortho / walpole_to_ortho" begin
    frame = CanonicalBasis{3, Float64}()

    # ── iso_to_ortho: ISO tensor rebuilt as TensOrtho yields identical array
    I4 = TensISO{3}(2.5, 4.0)
    O_from_iso = iso_to_ortho(I4, frame)
    @test O_from_iso isa TensOrtho{4, 3, Float64}
    @test opequal(get_array(I4), get_array(O_from_iso))

    # ── walpole_to_ortho (aligned axis)
    n = (0.0, 0.0, 1.0)
    W = ti_from_constants(10.0, 3.0, 2.5, 12.0, 2.0, n)
    k = TensND._axis_on_frame_index(n, frame)
    @test k == 3
    O_from_w = walpole_to_ortho(W, frame, k)
    @test O_from_w isa TensOrtho{4, 3, Float64}
    @test opequal(get_array(W), get_array(O_from_w))

    # Axis along e₁
    n1 = (1.0, 0.0, 0.0)
    W1 = ti_from_constants(10.0, 3.0, 2.5, 12.0, 2.0, n1)
    O_from_w1 = walpole_to_ortho(W1, frame, 1)
    @test opequal(get_array(W1), get_array(O_from_w1))

    # Axis along e₂, and an index that names no axis of the frame
    n2 = (0.0, 1.0, 0.0)
    W2 = ti_from_constants(10.0, 3.0, 2.5, 12.0, 2.0, n2)
    @test opequal(get_array(W2), get_array(walpole_to_ortho(W2, frame, 2)))
    @test_throws ErrorException walpole_to_ortho(W2, frame, 4)

    # A TI tensor is stored in the canonical frame, so its components there,
    # whatever the variance, are its array; an axis given as a 1×3 matrix
    # reads as the vector it holds.
    @test components(W, (:cont, :cont, :cont, :cont)) == get_array(W)
    @test TensND._extract_vec([0.0 0.0 1.0]) == (0.0, 0.0, 1.0)
end

# ═══════════════════════════════════════════════════════════════════════════
@testsection "Walpole & Ortho tensors: TensTI{4,T,8} — full axially-invariant algebra" begin
    ez = (0.0, 0.0, 1.0)
    nax = (0.6, 0.0, 0.8)   # unit, non-canonical axis

    # Mandel patterns of the antisymmetric generators (axis = e₃)
    M7 = mandel(tens_W7(ez))
    M8 = mandel(tens_W8(ez))
    s2 = 1 / sqrt(2.0)
    @test M7[4, 5] ≈ -1.0 atol = _WO_ATOL
    @test M7[5, 4] ≈ 1.0 atol = _WO_ATOL
    @test sum(abs, M7) ≈ 2.0 atol = _WO_ATOL
    @test M8[6, 1] ≈ s2 atol = _WO_ATOL
    @test M8[6, 2] ≈ -s2 atol = _WO_ATOL
    @test M8[1, 6] ≈ -s2 atol = _WO_ATOL
    @test M8[2, 6] ≈ s2 atol = _WO_ATOL
    @test sum(abs, M8) ≈ 4s2 atol = _WO_ATOL

    # constructor / accessors / lifts
    A8 = TensTI{4}(1.0, 2.0, 0.3, -0.4, 1.5, 2.5, 0.7, -0.2, nax)
    @test A8 isa TensTI{4, Float64, 8}
    @test walpole_coeffs(A8, Val(8)) == (1.0, 2.0, 0.3, -0.4, 1.5, 2.5, 0.7, -0.2)
    A5 = TensTI{4}(1.0, 2.0, 0.5, 3.0, 4.0, nax)
    @test walpole_coeffs(A5, Val(8)) == (1.0, 2.0, 0.5, 0.5, 3.0, 4.0, 0.0, 0.0)
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
    @test maximum(abs, mandel(Id) .- I(6)) < 1.0e-12

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
    f = x -> walpole_coeffs(inv(TensTI{4}(x, 2.0, 0.3, -0.4, 1.5, 2.5, 0.7, -0.2, ez)), Val(8))[1]
    g = FD.derivative(f, 1.0)
    h = 1.0e-6
    @test g ≈ (f(1.0 + h) - f(1.0 - h)) / 2h atol = 1.0e-6
end

# ═══════════════════════════════════════════════════════════════════════════
@testsection "Walpole & Ortho tensors: TensTI{2,T,3} — axially-invariant 2nd order with rotation part" begin
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
    @test r ≈ rsym atol = 1.0e-13

    # mixed-axis dot falls back generically
    t3c = TensTI{2}(1.5, -1.0, -0.3, ez)
    dm = t3 ⋅ t3c
    @test !(dm isa TensTI)
    @test maximum(abs, Array(get_array(dm)) .- (get_array(t3) * get_array(t3c))) < 1.0e-12
end

# ═══════════════════════════════════════════════════════════════════════════
@testsection "Walpole & Ortho tensors: TensOrtho — ForwardDiff compatibility (frame eltype decoupled)" begin
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
@testsection "Walpole & Ortho tensors: TensOrtho — concrete frame field, closed-form inv" begin
    frame = CanonicalBasis{3, Float64}()
    O = TensOrtho(
        250.0e3, 180.0e3, 150.0e3, 90.0e3, 80.0e3, 70.0e3,
        60.0e3, 50.0e3, 40.0e3, frame
    )
    # The frame field must be a concrete type: TensOrtho decouples the
    # frame eltype from T (for ForwardDiff) via a *type parameter*, not by
    # erasing it to the abstract `OrthonormalBasis{3}` — an abstract field
    # boxes every access and defeats inference.
    @test isconcretetype(fieldtype(typeof(O), :frame))
    @test typeof(O) === TensOrtho{4, 3, Float64, 9, CanonicalBasis{3, Float64}}

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
    Km = TensND.mandel_material(O)
    Km_inv_dense = inv(Km)
    Km_inv_closed = TensND.mandel_material(inv(O))
    @test maximum(abs.(Km_inv_closed .- Km_inv_dense)) < 1.0e-10
    @test maximum(abs.(Km * Km_inv_closed - Matrix(I, 6, 6))) < 1.0e-12
end

@testsection "Walpole & Ortho tensors: TensOrtho — Kelvin-Mandel congruence of the material frame" begin
    # `mandel(t)` is the canonical-frame matrix, `mandel_material(t)` the
    # block-diagonal one in the material frame.  The closed-form
    # `dcontract` rests entirely on the congruence relating them, so pin it
    # explicitly — and on a ROTATED frame, since in the canonical frame
    # `Q == I` and every convention error hides.
    for angles in ((0.0, 0.0, 0.0), (0.7, -0.4, 1.1), (-1.3, 2.2, 0.35))
        ℬ = RotatedBasis(angles...)
        O = TensOrtho(160.0, 120.0, 95.0, 55.0, 48.0, 42.0, 33.0, 29.0, 27.0, ℬ)
        E = TensND.vecbasis(TensND.frame(O), :cov)
        Q = TensND._km_congruence(E)
        @test maximum(abs, Q * transpose(Q) - I) < 1.0e-14      # orthogonal
        @test maximum(abs, mandel(O) - Q * TensND.mandel_material(O) * transpose(Q)) < 1.0e-12
    end
end

@testsection "Walpole & Ortho tensors: TensOrtho — closed-form double contraction" begin
    dense(t) = [t[i, j, k, l] for i in 1:3, j in 1:3, k in 1:3, l in 1:3]
    ref(A, B) = Tensors.dcontract(TensND.generic_tens(A), TensND.generic_tens(B))

    for angles in ((0.0, 0.0, 0.0), (0.7, -0.4, 1.1), (-1.3, 2.2, 0.35))
        ℬ = RotatedBasis(angles...)
        A = TensOrtho(160.0, 120.0, 95.0, 55.0, 48.0, 42.0, 33.0, 29.0, 27.0, ℬ)
        B = TensOrtho(11.0, 7.0, 5.0, 3.0, 2.0, 1.5, 4.0, 3.5, 2.5, ℬ)
        R = dense(ref(A, B))
        @test maximum(abs, dense(A ⊡ B) - R) / maximum(abs, R) < 1.0e-14

        # Isotropic operands are promoted onto the material frame.
        I4 = TensISO{3}(3.0, 2.0)
        for (l, r) in ((I4, A), (A, I4))
            Rlr = dense(ref(l, r))
            @test maximum(abs, dense(l ⊡ r) - Rlr) / maximum(abs, Rlr) < 1.0e-14
        end

        # `A ⊡ inv(A)` is the identity on symmetric 2nd-order tensors.
        @test maximum(abs, dense(A ⊡ inv(A)) - dense(TensND.generic_tens(TensISO{3}(1.0, 1.0)))) < 1.0e-12
    end

    # Different material frames: the product is generally fully anisotropic
    # and must fall back to the generic route, bit for bit.
    A = TensOrtho(160.0, 120.0, 95.0, 55.0, 48.0, 42.0, 33.0, 29.0, 27.0, RotatedBasis(0.7, -0.4, 1.1))
    B = TensOrtho(11.0, 7.0, 5.0, 3.0, 2.0, 1.5, 4.0, 3.5, 2.5, RotatedBasis(-1.3, 2.2, 0.35))
    @test dense(A ⊡ B) == dense(ref(A, B))

    # The result is orthotropic WITHOUT major symmetry: the product of two
    # symmetric 3×3 Kelvin-Mandel blocks is not symmetric unless they
    # commute.  This is why `dcontract` returns the twelve-constant form,
    # the same widening as `TensTI{4}` N=5 → N=6.
    ℬ = RotatedBasis(0.7, -0.4, 1.1)
    A = TensOrtho(160.0, 120.0, 95.0, 55.0, 48.0, 42.0, 33.0, 29.0, 27.0, ℬ)
    B = TensOrtho(11.0, 7.0, 5.0, 3.0, 2.0, 1.5, 4.0, 3.5, 2.5, ℬ)
    P = A ⊡ B
    @test P isa TensOrtho{4, 3, Float64, 12}
    p = get_array(P)
    @test maximum(abs, p - permutedims(p, (3, 4, 1, 2))) > 1.0
    @test Tensors.isminorsymmetric(P)
    @test p ≈ permutedims(p, (2, 1, 3, 4)) && p ≈ permutedims(p, (1, 2, 4, 3))
end

# A dual number compares by its value only: a tensor symmetric in value whose
# derivative is not minor-symmetric must keep its full storage, or the
# derivative is symmetrized. Reference: the derivative of x ↦ S₀ + x B, which
# is B, at first and second order.
@testsection "Storage of dual tensors symmetric in value only" begin
    rng = Random.Xoshiro(20261128)
    B = Tensor{4, 3}(randn(rng, 3, 3, 3, 3))
    S₀ = Tensor{4, 3}(Tensors.symmetric(Tensor{4, 3}(randn(rng, 3, 3, 3, 3))))
    f(x) = Array(canonical_components(Tens(S₀ + x * B)))
    @test ForwardDiff.derivative(f, 0.0) ≈ Array(B)
    g(x) = Array(canonical_components(Tens(S₀ + x^2 * B)))
    @test ForwardDiff.derivative(y -> ForwardDiff.derivative(g, y), 0.0) ≈ 2 * Array(B)
    # a dual tensor symmetric in value and derivatives keeps the symmetric storage
    h(x) = Tens(S₀ + x * Tensor{4, 3}(Tensors.symmetric(B)))
    @test get_array(h(ForwardDiff.Dual(0.0, 1.0))) isa SymmetricTensor
    @test get_array(Tens(S₀ + ForwardDiff.Dual(0.0, 1.0) * B)) isa Tensor
    # A derivative that vanishes in exact arithmetic but carries the round-off of
    # the others (as another processor leaves it) does not decide the storage:
    # its asymmetry is noise at the scale of the derivatives.
    A = Array(Tensor{4, 3}(Tensors.symmetric(B)))
    noise = zeros(3, 3, 3, 3)
    noise[1, 2, 1, 1] = eps(maximum(abs, A))
    duals = ForwardDiff.Dual.(Array(S₀), A, noise)
    @test get_array(Tens(Tensor{4, 3}(duals))) isa SymmetricTensor
    # nested duals: the same judgment at every order
    nested = ForwardDiff.Dual.(ForwardDiff.Dual.(Array(S₀), A), ForwardDiff.Dual.(noise, A))
    @test get_array(Tens(Tensor{4, 3}(nested))) isa SymmetricTensor
    nested_asym = ForwardDiff.Dual.(ForwardDiff.Dual.(Array(S₀), A), ForwardDiff.Dual.(noise, Array(B)))
    @test get_array(Tens(Tensor{4, 3}(nested_asym))) isa Tensor
    # exact derivatives are compared exactly
    Q = Array(Tensor{4, 3}(Tensors.symmetric(Tensor{4, 3}(rand(rng, -3:3, 3, 3, 3, 3) .// 2))))
    @test get_array(Tens(Tensor{4, 3}(ForwardDiff.Dual.(Q, Q)))) isa SymmetricTensor
    Qa = copy(Q); Qa[1, 2, 1, 1] += 1 // 10^20
    @test get_array(Tens(Tensor{4, 3}(ForwardDiff.Dual.(Q, Qa)))) isa Tensor
    # a dual number without partials
    @test get_array(Tens(Tensor{4, 3}(ForwardDiff.Dual.(Array(S₀))))) isa SymmetricTensor
end

@testsection "Walpole & Ortho tensors: minor-symmetry storage does not hinge on round-off" begin
    # Regression. `Tensors.issymmetric` compares components exactly, so the
    # choice between a 36-component `SymmetricTensor` and an 81-component
    # `Tensor` used to be decided by the last bit: writing a structured
    # tensor about a non-canonical axis leaves minor-antisymmetric residue
    # of order 1e-16, and 8.9e-16 on components of size 30 landed on the
    # full form where 4.4e-16 landed on the symmetric one.
    #
    # The two are not interchangeable. `inv` of the full fourth-order form
    # solves a 9x9 system whose three minor-antisymmetric directions are
    # null, so it threw `SingularException` on a well-posed problem — and
    # which of the two a model got depended on the axis it happened to use.
    for θ in range(0, π; length = 13), ϕ in (0.0, 0.7, 2.1)
        ℬ = RotatedBasis(θ, ϕ, 0.0)
        n = vecbasis(ℬ)[:, 3]
        for t in (
                TensTI{4}(31.0, 39.0, 11.313708498984763, 23.0, 23.0, n),
                TensTI{4}(2.0, 3.0, 1.0, 1.5, 4.0, 5.0, n),
                TensOrtho(160.0, 120.0, 95.0, 55.0, 48.0, 42.0, 33.0, 29.0, 27.0, ℬ),
            )
            @test get_array(TensND.generic_tens(t)) isa Tensors.SymmetricTensor
        end
    end

    # The failure this came from: a tensor whose VALUES are isotropic but
    # which carries a tilted transversely isotropic axis. It is bit-for-bit
    # the same tensor as the one built about e₃, so it must invert the same
    # way rather than raise.
    δ(n) = TensTI{4}(-31.0, -39.0, -11.313708498984763, -23.0, -23.0, n)
    e₃ = (0.0, 0.0, 1.0)
    tilted = vecbasis(RotatedBasis(0.4, 0.7, 0.0))[:, 3]
    @test get_array(δ(e₃)) ≈ get_array(δ(tilted)) atol = 1.0e-12

    𝕀 = TensISO{3}(1.0, 1.0)
    ℙ = TensND.generic_tens(TensISO{3}(0.05, 0.03))
    for n in (e₃, tilted)
        S = 𝕀 + ℙ ⊡ δ(n)
        @test get_array(TensND.generic_tens(S)) isa Tensors.SymmetricTensor
        A = inv(S)                      # this is what used to throw
        @test A ⊡ S ≈ 𝕀 atol = 1.0e-10
    end
end

# ═══════════════════════════════════════════════════════════════════════════
@testsection "Walpole & Ortho tensors: express into a rotated basis actually rotates" begin
    # REGRESSION. `express(::TensTI, ℬ)` and `express(::TensOrtho, ℬ)`
    # used to wrap the stored components — which live in the CANONICAL frame
    # — in the target basis without transforming them. For a canonical `ℬ`
    # that is right; for a rotated one it relabels the components and
    # silently returns a DIFFERENT physical tensor. Only `TensISO` may skip
    # the rotation.
    #
    # The oracle is the generic path: the very same components wrapped as a
    # plain `Tens` in the canonical basis, for which the rotation has always
    # been correct.
    ℬrot = RotatedBasis(0.4, 0.7, 0.3)

    structured = (
        "TensTI{4}" => TensTI{4}(20.0, 30.0, 4.0, 5.0, 8.0, (0.0, 0.0, 1.0)),
        "TensTI{2}" => TensTI{2}(3.0, 7.0, (0.0, 0.0, 1.0)),
        "TensOrtho" => TensOrtho(
            20.0, 22.0, 25.0, 8.0, 9.0, 10.0, 6.0, 7.0, 5.0,
            CanonicalBasis{3, Float64}()
        ),
    )

    for (name, t) in structured
        generic = Tens(get_array(t))                 # same tensor, generic storage
        @test get_array(generic) ≈ get_array(t) atol = _WO_ATOL

        rot_struct = express(t, ℬrot)
        rot_generic = express(generic, ℬrot)

        # The two storages must agree component by component.
        @test get_array(rot_struct) ≈ get_array(rot_generic) atol = _WO_ATOL

        # …and the rotation must be a real one: a non-trivial basis change
        # cannot leave the components untouched, which is exactly what the
        # bug did.
        @test !isapprox(get_array(rot_struct), get_array(t); atol = 1.0e-8)

        # Round trip back to the canonical basis restores the tensor.
        back = express(rot_struct, CanonicalBasis{3, Float64}())
        @test get_array(back) ≈ get_array(t) atol = 1.0e-10

        # A canonical target stays a no-op.
        @test get_array(express(t, CanonicalBasis{3, Float64}())) ≈
            get_array(t) atol = _WO_ATOL

        # `components(t, ℬ, var)` is the primitive behind `express` and
        # must tell the same story.
        v = ntuple(_ -> :cont, Val(length(size(t))))
        @test components(t, ℬrot, v) ≈ components(generic, ℬrot, v) atol = _WO_ATOL
    end

    # `TensISO` is genuinely invariant — its fast path is legitimate and must
    # be preserved.
    for iso in (TensISO{3}(2.0), TensISO{3}(3.0, 4.0))
        @test get_array(express(iso, ℬrot)) ≈ get_array(iso) atol = _WO_ATOL
    end
end

@testsection "ti_constants of the six- and eight-parameter forms" begin
    n = (0.36, 0.48, 0.8)
    ℓ = (3.0, 2.0, 0.7, 4.0, 1.5)
    five = ti_constants(TensTI{4, Float64, 5}(ℓ, n))
    # The same tensor in the longer forms: the same five constants.
    @test ti_constants(TensTI{4, Float64, 6}((3.0, 2.0, 0.7, 0.7, 4.0, 1.5), n)) == five
    @test ti_constants(TensTI{4, Float64, 8}((3.0, 2.0, 0.7, 0.7, 4.0, 1.5, 0.0, 0.0), n)) == five
    # Not major-symmetric: refused, not read as if ℓ₄, ℓ₇, ℓ₈ were not there.
    @test_throws ArgumentError ti_constants(TensTI{4, Float64, 6}((3.0, 2.0, 0.7, 0.9, 4.0, 1.5), n))
    @test_throws ArgumentError ti_constants(TensTI{4, Float64, 8}((3.0, 2.0, 0.7, 0.7, 4.0, 1.5, 0.2, 0.0), n))
    @test ti_constants(majorsymmetric(TensTI{4, Float64, 6}((3.0, 2.0, 0.7, 0.9, 4.0, 1.5), n)))[3] ≈ 0.8 / √2
    # A product of TI tensors about one axis is major-symmetric to rounding.
    A = TensTI{4, Float64, 5}(ℓ, n)
    P = A ⊡ A
    @test collect(ti_constants(P)) ≈ collect(ti_constants(majorsymmetric(P)))
    # Duals and symbols.
    @test ForwardDiff.derivative(x -> ti_constants(TensTI{4, typeof(x), 6}((x, 2.0, 0.7, 0.7, 4.0, 1.5), n))[4], 3.0) == 1
    TensND.Symbolics.@variables a b
    @test isequal(ti_constants(TensTI{4, TensND.Symbolics.Num, 6}((a, b, a, a, b, a), n))[4], a)
    @test_throws ArgumentError ti_constants(TensTI{4, TensND.Symbolics.Num, 6}((a, b, a, b, b, a), n))
end

@testsection "TI constants of a stiffness, directly, against the inverse" begin
    rng = Xoshiro(20261030)
    n = (0.0, 0.6, 0.8)
    for _ in 1:200
        # A positive definite TI stiffness, drawn on its Walpole blocks.
        A = randn(rng, 2, 2)
        ℓ = A * transpose(A) + 0.1I
        C = TensTI{4}(ℓ[2, 2], ℓ[1, 1], ℓ[1, 2], 0.1 + rand(rng), 0.1 + rand(rng), n)
        S = inv(C)
        for form in (:engineering, :hoenig)
            @test all(isapprox.(ti_constants(C; form, role = :stiffness), ti_constants(S; form); rtol = 1.0e-9))
            c = ti_constants(S; form)
            @test all(isapprox.(get_data(ti_from_constants(c..., n; form, role = :stiffness)), get_data(C); rtol = 1.0e-9))
            @test ti_constants(S; form, role = :compliance) === c
        end
        # `role` has no effect on the components.
        @test ti_constants(C; role = :stiffness) === ti_constants(C)
        @test ti_from_constants(ti_constants(C)..., n; role = :stiffness) == ti_from_constants(ti_constants(C)..., n)
    end
    @test_throws "unknown role" ti_constants(TensTI{4}(2.0, 3.0, 0.5, 1.0, 1.0, n); form = :hoenig, role = :rigidity)
    # Exactly, on symbolic components.
    c11, c12, c13, c33, c44 = symbols("C11 C12 C13 C33 C44", positive = true)
    Cs = ti_from_constants(c11, c12, c13, c33, c44, [0, 0, 1])
    e1, e3, n12, n31, g = ti_constants(Cs; form = :engineering, role = :stiffness)
    D = c11 * c33 - c13^2
    @test simplify(e3 - (c33 - 2 * c13^2 / (c11 + c12))) == 0
    @test simplify(e1 - (c11 - c12) * (c11 * c33 + c12 * c33 - 2 * c13^2) / D) == 0
    @test simplify(n31 - c13 / (c11 + c12)) == 0 && g == c44
    # Derivatives through the direct formulas.
    C0 = TensTI{4}(30.0, 20.0, 5.0, 8.0, 6.0, n)
    dE3 = ForwardDiff.derivative(x -> ti_constants(TensTI{4}(x, 20.0, 5.0, 8.0, 6.0, n); form = :engineering, role = :stiffness)[2], 30.0)
    @test dE3 ≈ 1    # E₃ = C₃₃₃₃ − …, and ℓ₁ = C₃₃₃₃
    stiff_eng(t) = ti_constants(t; form = :engineering, role = :stiffness)
    @test (@inferred stiff_eng(C0)) isa NTuple{5, Float64}
end

# ═══════════════════════════════════════════════════════════════════════════
# Hill's moduli (k, l, m, n, p). References independent of the code under test:
# the components read by contraction with a frame drawn about a random axis; the
# compliance obtained by `inv`, whose components give the engineering constants;
# the isotropic tensor.
@testsection "Walpole & Ortho tensors: Hill's moduli" begin
    rng = Xoshiro(20261126)
    contraction(C, a, b, c, d) = sum(C[i, j, k, l] * a[i] * b[j] * c[k] * d[l] for i in 1:3, j in 1:3, k in 1:3, l in 1:3)
    for _ in 1:20
        e₃ = normalize(randn(rng, 3))
        e₁ = normalize(cross(e₃, randn(rng, 3)))
        e₂ = cross(e₃, e₁)
        m, p = 0.1 + rand(rng), 0.1 + rand(rng)
        k, n = 1 + rand(rng), 1 + rand(rng)
        l = (2 * rand(rng) - 1) * sqrt(k * n) / 2          # k n > l²: positive definite
        ℂ = ti_from_constants(k, l, m, n, p, e₃; form = :hill)
        @test ℂ isa TensTI{4, Float64, 5}
        C = canonical_components(ℂ)
        @test (contraction(C, e₁, e₁, e₁, e₁) + contraction(C, e₁, e₁, e₂, e₂)) / 2 ≈ k
        @test contraction(C, e₁, e₁, e₃, e₃) ≈ l
        @test (contraction(C, e₁, e₁, e₁, e₁) - contraction(C, e₁, e₁, e₂, e₂)) / 2 ≈ m
        @test contraction(C, e₁, e₂, e₁, e₂) ≈ m
        @test contraction(C, e₃, e₃, e₃, e₃) ≈ n
        @test contraction(C, e₂, e₃, e₂, e₃) ≈ p
        # Read back: by name, as a tuple, to the last bit; `role` has no effect.
        h = hill_moduli(ℂ)
        @test keys(h) == (:k, :l, :m, :n, :p)
        @test values(h) === ti_constants(ℂ; form = :hill) === ti_constants(ℂ; form = :hill, role = :stiffness)
        @test h.k == k && h.m == m && h.n == n && h.p == p && h.l ≈ l
        @test ti_from_constants(k, l, m, n, p, e₃; form = :hill, role = :stiffness) == ℂ
        # The other forms of the same tensor.
        @test collect(ti_constants(ℂ)) ≈ [k + m, k - m, l, n, p]
        @test ti_from_constants(k + m, k - m, l, n, p, e₃) ≈ ℂ
        # The engineering constants of the material whose stiffness is ℂ, from
        # the components of its compliance.
        S = canonical_components(inv(ℂ))
        E₁, E₃, ν₁₂, ν₃₁, G₃₁ = ti_constants(ℂ; form = :engineering, role = :stiffness)
        @test E₃ ≈ 1 / contraction(S, e₃, e₃, e₃, e₃) ≈ n - l^2 / k
        @test E₁ ≈ 1 / contraction(S, e₁, e₁, e₁, e₁) ≈ 4 * m * (k * n - l^2) / ((k + m) * n - l^2)
        @test ν₃₁ ≈ -E₃ * contraction(S, e₁, e₁, e₃, e₃) ≈ l / (2k)
        @test ν₁₂ ≈ -E₁ * contraction(S, e₁, e₁, e₂, e₂) ≈ ((k - m) * n - l^2) / ((k + m) * n - l^2)
        @test G₃₁ ≈ p
        # … and back to Hill's moduli through each form and each role.
        for form in (:engineering, :hoenig)
            back = ti_from_constants(ti_constants(ℂ; form, role = :stiffness)..., e₃; form, role = :stiffness)
            @test collect(values(hill_moduli(back))) ≈ [k, l, m, n, p]
            𝕊 = ti_from_constants(ti_constants(ℂ; form, role = :stiffness)..., e₃; form)
            @test collect(values(hill_moduli(inv(𝕊)))) ≈ [k, l, m, n, p]
        end
        # A product of TI tensors is stored on six coefficients: read when it is
        # major-symmetric, refused otherwise.
        @test collect(values(hill_moduli(ℂ ⊡ one(ℂ)))) ≈ [k, l, m, n, p]
        @test_throws "not major-symmetric" hill_moduli(ℂ ⊡ ti_from_constants(n, l / 2, p, k, m, e₃; form = :hill))
    end

    # The isotropic tensor, about any axis.
    κ, μ = 10.0, 3.0
    iso = hill_moduli(iso_stiffness(κ, μ))
    @test collect(values(iso)) ≈ [κ + μ / 3, κ - 2μ / 3, μ, κ + 4μ / 3, μ]
    @test canonical_components(ti_from_constants(iso..., normalize([1.0, -2.0, 0.5]); form = :hill)) ≈ canonical_components(iso_stiffness(κ, μ))
    @test collect(values(hill_moduli(iso_to_ti(iso_stiffness(κ, μ), _WO_N3)))) ≈ collect(values(iso))

    # Inferred, for numbers and for dual numbers.
    from_hill(k) = ti_from_constants(k, 5.0, 11.0, 12.0, 4.0, (0.0, 0.0, 1.0); form = :hill)
    @test (@inferred from_hill(19.0)) isa TensTI{4, Float64, 5}
    @test (@inferred hill_moduli(from_hill(19.0))) === (k = 19.0, l = 5.0, m = 11.0, n = 12.0, p = 4.0)
    @test (@inferred hill_moduli(iso_stiffness(κ, μ))) isa NamedTuple{(:k, :l, :m, :n, :p), NTuple{5, Float64}}
    # Integers keep an exact tensor: only `l` meets √2.
    @test hill_moduli(ti_from_constants(19, 5, 11, 12, 4, (0, 0, 1); form = :hill)) == (k = 19.0, l = 5.0, m = 11.0, n = 12.0, p = 4.0)

    # Derivatives, against the closed forms above: first order, then nested.
    E₃(k, l, n) = ti_constants(ti_from_constants(k, l, 11.0, n, 4.0, (0.0, 0.0, 1.0); form = :hill); form = :engineering, role = :stiffness)[2]
    @test ForwardDiff.gradient(x -> E₃(x...), [19.0, 5.0, 12.0]) ≈ [5.0^2 / 19.0^2, -2 * 5.0 / 19.0, 1.0]
    @test ForwardDiff.derivative(l -> ForwardDiff.derivative(k -> E₃(k, l, 12.0), 19.0), 5.0) ≈ 2 * 5.0 / 19.0^2
    @test ForwardDiff.derivative(k -> hill_moduli(inv(inv(from_hill(k)))).k, 19.0) ≈ 1
    @test ForwardDiff.derivative(κ -> hill_moduli(iso_stiffness(κ, μ)).l, κ) ≈ 1

    # Exactly, on SymPy symbols.
    k, l, m, n, p = symbols("k l m n p", positive = true)
    ℂ = ti_from_constants(k, l, m, n, p, _WO_N3S; form = :hill)
    @test get_data(ℂ) == (n, 2k, sqrt(Sym(2)) * l, 2m, 2p)
    @test all(iszero, simplify.(collect(values(hill_moduli(ℂ))) .- [k, l, m, n, p]))
    @test all(iszero, simplify.(collect(ti_constants(ℂ)) .- [k + m, k - m, l, n, p]))
    E₁, E₃ˢ, ν₁₂, ν₃₁, G₃₁ = ti_constants(ℂ; form = :engineering, role = :stiffness)
    @test iszero(simplify(E₃ˢ - (n - l^2 / k))) && iszero(simplify(ν₃₁ - l / (2k))) && G₃₁ == p
    @test iszero(simplify(E₁ - 4 * m * (k * n - l^2) / ((k + m) * n - l^2)))
    @test iszero(simplify(ν₁₂ - ((k - m) * n - l^2) / ((k + m) * n - l^2)))
    # The compliance has Hill's "moduli" too: those of its own components.
    hS = hill_moduli(inv(ℂ))
    @test iszero(simplify(hS.n - k / (k * n - l^2))) && iszero(simplify(hS.k - n / (4 * (k * n - l^2))))
    @test iszero(simplify(hS.l + l / (2 * (k * n - l^2)))) && iszero(simplify(hS.m - 1 / (4m))) && iszero(simplify(hS.p - 1 / (4p)))
    κˢ, μˢ = symbols("kappa mu", positive = true)
    @test all(iszero, simplify.(collect(values(hill_moduli(iso_stiffness(κˢ, μˢ)))) .- [κˢ + μˢ / 3, κˢ - 2μˢ / 3, μˢ, κˢ + 4μˢ / 3, μˢ]))

    # And on Symbolics' ones.
    Symbolics.@variables kₛ lₛ mₛ nₛ pₛ
    ℂₛ = ti_from_constants(kₛ, lₛ, mₛ, nₛ, pₛ, (0, 0, 1); form = :hill)
    @test eltype(ℂₛ) === Symbolics.Num
    at = Dict(kₛ => 19.0, lₛ => 5.0, mₛ => 11.0, nₛ => 12.0, pₛ => 4.0)
    num(x) = Symbolics.value(Symbolics.substitute(x, at; fold = Val(true)))
    @test collect(num.(values(hill_moduli(ℂₛ)))) ≈ [19.0, 5.0, 11.0, 12.0, 4.0]
    @test collect(num.(ti_constants(ℂₛ; form = :hill))) ≈ [19.0, 5.0, 11.0, 12.0, 4.0]
    @test collect(num.(ti_constants(ℂₛ))) ≈ [30.0, 8.0, 5.0, 12.0, 4.0]
end

# ═══════════════════════════════════════════════════════════════════════════
# The Walpole coefficients as a form of the five constants. References: the
# coefficients the tensor was built with; its components by contraction; a
# product of TI tensors, stored on six, whose symmetric part is known.
@testsection "Walpole & Ortho tensors: the Walpole form" begin
    rng = Xoshiro(20261127)
    contraction(C, a, b, c, d) = sum(C[i, j, k, l] * a[i] * b[j] * c[k] * d[l] for i in 1:3, j in 1:3, k in 1:3, l in 1:3)
    for _ in 1:20
        e₃ = normalize(randn(rng, 3))
        e₁ = normalize(cross(e₃, randn(rng, 3)))
        ℓ = (1 + rand(rng), 2 + rand(rng), rand(rng) - 0.5, 0.1 + rand(rng), 0.1 + rand(rng))
        ℂ = ti_from_constants(ℓ..., e₃; form = :walpole)
        @test ℂ isa TensTI{4, Float64, 5} && ℂ == TensTI{4}(ℓ..., e₃)
        C = canonical_components(ℂ)
        @test contraction(C, e₃, e₃, e₃, e₃) ≈ ℓ[1]                 # ℓ₁ = C₃₃₃₃
        @test contraction(C, e₁, e₁, e₃, e₃) ≈ ℓ[3] / √2            # ℓ₃ = √2 C₁₁₃₃
        # Read back in the three lengths, to the last bit; `role` has no effect.
        @test ti_constants(ℂ; form = :walpole) === walpole_coeffs(ℂ, Val(5)) === ℓ
        @test ti_constants(ℂ; form = :walpole, role = :stiffness) === ℓ
        @test ti_from_constants(ℓ..., e₃; form = :walpole, role = :stiffness) == ℂ
        @test walpole_coeffs(ℂ) === (ℓ[1], ℓ[2], ℓ[3], ℓ[3], ℓ[4], ℓ[5])
        @test walpole_coeffs(ℂ, Val(8)) === (ℓ[1], ℓ[2], ℓ[3], ℓ[3], ℓ[4], ℓ[5], 0.0, 0.0)
        # The other forms of the same tensor.
        h = hill_moduli(ℂ)
        @test collect(ℓ) ≈ [h.n, 2h.k, √2 * h.l, 2h.m, 2h.p]
        @test ti_from_constants(ti_constants(ℂ)..., e₃) ≈ ℂ
        # A tensor stored on six or eight coefficients: five when it is
        # major-symmetric, refused otherwise.
        six = ℂ ⊡ one(ℂ)
        @test six isa TensTI{4, Float64, 6} && collect(walpole_coeffs(six, Val(5))) ≈ collect(ℓ)
        skew = ℂ ⊡ ti_from_constants(ℓ[2], ℓ[1], ℓ[3] / 2, ℓ[5], ℓ[4], e₃; form = :walpole)
        @test_throws "not major-symmetric" walpole_coeffs(skew, Val(5))
        @test_throws "not major-symmetric" ti_constants(skew; form = :walpole)
        @test length(walpole_coeffs(skew)) == 6 && !(walpole_coeffs(skew)[3] ≈ walpole_coeffs(skew)[4])
    end

    from_walpole(ℓ₁) = ti_from_constants(ℓ₁, 38.0, 7.0, 22.0, 8.0, (0.0, 0.0, 1.0); form = :walpole)
    @test (@inferred from_walpole(12.0)) isa TensTI{4, Float64, 5}
    five(t) = walpole_coeffs(t, Val(5))
    @test (@inferred five(from_walpole(12.0))) === (12.0, 38.0, 7.0, 22.0, 8.0)
    @test (@inferred five(from_walpole(12.0) ⊡ one(from_walpole(12.0)))) isa NTuple{5, Float64}

    # Derivatives: ℓ₁ is read back as given, and E₃ = ℓ₁ − ℓ₃²/ℓ₂ of a stiffness.
    @test ForwardDiff.derivative(x -> ti_constants(from_walpole(x); form = :walpole)[1], 12.0) == 1
    E₃(ℓ₁, ℓ₃) = ti_constants(ti_from_constants(ℓ₁, 38.0, ℓ₃, 22.0, 8.0, (0.0, 0.0, 1.0); form = :walpole); form = :engineering, role = :stiffness)[2]
    @test ForwardDiff.gradient(x -> E₃(x...), [12.0, 7.0]) ≈ [1, -2 * 7.0 / 38.0]
    @test ForwardDiff.derivative(a -> ForwardDiff.derivative(b -> E₃(a, b), 7.0), 12.0) == 0
    @test ForwardDiff.derivative(b -> ForwardDiff.derivative(c -> E₃(12.0, c), b), 7.0) ≈ -2 / 38.0

    # Exactly, on SymPy symbols and on Symbolics' ones.
    ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆ = symbols("l1 l2 l3 l5 l6", positive = true)
    ℂ = ti_from_constants(ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆, _WO_N3S; form = :walpole)
    @test ti_constants(ℂ; form = :walpole) == (ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆)
    h = hill_moduli(ℂ)
    @test all(iszero, simplify.([h.n - ℓ₁, 2h.k - ℓ₂, sqrt(Sym(2)) * h.l - ℓ₃, 2h.m - ℓ₅, 2h.p - ℓ₆]))
    @test all(iszero, simplify.(collect(walpole_coeffs(inv(inv(ℂ)), Val(5))) .- [ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆]))
    Symbolics.@variables a₁ a₂ a₃ a₅ a₆
    ℂₛ = ti_from_constants(a₁, a₂, a₃, a₅, a₆, (0, 0, 1); form = :walpole)
    @test eltype(ℂₛ) === Symbolics.Num
    @test all(isequal.(ti_constants(ℂₛ; form = :walpole), (a₁, a₂, a₃, a₅, a₆)))
end
