# Conventions asserted by the documentation.
#
# Every test here corresponds to a statement made on a Theory page. Their job is
# not to find bugs in new code but to make it impossible for the code to drift
# away from what the documentation says without something going red. See
# `docs/src/developer/testing_conventions.md`.

@testsection "Documented conventions" begin

    arr(x) = get_array(x)
    n = (0.0, 0.0, 1.0)
    W = walpole_basis(n)
    Ws = walpole_basis_sym(n)
    𝕀, 𝕁, 𝕂 = ISO(Val(3), Val(Float64))
    s2 = sqrt(2)

    # ── The Walpole basis ────────────────────────────────────────────────────
    # docs/src/theory/walpole.md

    @testset "Walpole: decomposition of 𝕀, 𝕁, 𝕂" begin
        @test arr(𝕀) ≈ arr(W[1]) + arr(W[2]) + arr(W[5]) + arr(W[6])
        @test arr(𝕁) ≈ (arr(W[1]) + 2arr(W[2]) + s2 * arr(W[3]) + s2 * arr(W[4])) / 3
        @test arr(𝕂) ≈ (2arr(W[1]) + arr(W[2]) - s2 * (arr(W[3]) + arr(W[4]))) / 3 +
            arr(W[5]) + arr(W[6])

        # And the two widely repeated identities that are FALSE. If either of
        # these ever passes, the basis convention has changed and the theory
        # page must be revisited.
        @test !(arr(𝕀) ≈ sum(arr(W[i]) for i in 1:6))
        @test !(arr(𝕁) ≈ arr(W[1]) + arr(W[2]))
    end

    @testset "Walpole: Gram matrix is diagonal (1,1,1,1,2,2)" begin
        gram = [sum(arr(W[i]) .* arr(W[j])) for i in 1:6, j in 1:6]
        @test gram ≈ Diagonal([1.0, 1.0, 1.0, 1.0, 2.0, 2.0])
    end

    @testset "Walpole: multiplication table" begin
        # (𝕎₁,𝕎₂,𝕎₃,𝕎₄) behave as the matrix units (E₁₁,E₂₂,E₁₂,E₂₁);
        # 𝕎₅ and 𝕎₆ are orthogonal idempotents.
        expected = Dict(
            (1, 1) => 1, (1, 3) => 3,
            (2, 2) => 2, (2, 4) => 4,
            (3, 2) => 3, (3, 4) => 1,
            (4, 1) => 4, (4, 3) => 2,
            (5, 5) => 5, (6, 6) => 6,
        )
        for i in 1:6, j in 1:6
            P = arr(W[i]) ⊡ arr(W[j])
            if haskey(expected, (i, j))
                @test P ≈ arr(W[expected[(i, j)]])
            else
                @test norm(P) < 1.0e-12
            end
        end
    end

    @testset "Walpole: symmetrized basis and its Gram" begin
        @test arr(Ws[1]) ≈ arr(W[1])
        @test arr(Ws[2]) ≈ arr(W[2])
        @test arr(Ws[3]) ≈ arr(W[3]) + arr(W[4])
        @test arr(Ws[4]) ≈ arr(W[5])
        @test arr(Ws[5]) ≈ arr(W[6])
        # These denominators are exactly the ones in `_project_TI_KM`.
        @test [sum(arr(Ws[i]) .* arr(Ws[i])) for i in 1:5] ≈ [1.0, 1.0, 2.0, 2.0, 2.0]
    end

    @testset "Walpole: N=5 ⊡ N=5 widens to N=6" begin
        A = TensTI{4}(2.0, 1.0, 0.5, 0.9, 1.3, n)
        B = TensTI{4}(3.0, 1.5, 0.2, 1.1, 0.4, n)
        P = A ⊡ B
        @test P isa TensTI{4, Float64, 6}
        @test arr(P) ≈ arr(A) ⊡ arr(B)
        ℓ = get_ℓ(P)
        @test !isapprox(ℓ[3], ℓ[4])          # major symmetry genuinely lost
    end

    # ── The extended Walpole algebra ─────────────────────────────────────────
    # docs/src/theory/walpole_extended.md

    @testset "𝕎₇, 𝕎₈ are major-antisymmetric" begin
        for Wx in (tens_W7(n), tens_W8(n))
            a = arr(Wx)
            @test norm(a + permutedims(a, (3, 4, 1, 2))) < 1.0e-12
        end
    end

    @testset "𝕎₇, 𝕎₈ annihilate axially invariant order-2 tensors — and only those" begin
        𝟏ₙ = [0.0 0 0; 0 0 0; 0 0 1]
        𝟏ₜ = [1.0 0 0; 0 1 0; 0 0 0]
        for Wx in (tens_W7(n), tens_W8(n))
            @test norm(arr(Wx) ⊡ 𝟏ₙ) < 1.0e-12
            @test norm(arr(Wx) ⊡ 𝟏ₜ) < 1.0e-12
            @test norm(arr(Wx) ⊡ (3.0𝟏ₜ + 5.0𝟏ₙ)) < 1.0e-12
        end
        # A *general* symmetric tensor is NOT annihilated. A source comment once
        # claimed otherwise; the documentation now says the opposite, so this is
        # the test that keeps it honest.
        g = [2.0 0.7 0.3; 0.7 1.5 0.4; 0.3 0.4 3.0]
        @test norm(arr(tens_W7(n)) ⊡ g) > 1.0e-3
        @test norm(arr(tens_W8(n)) ⊡ g) > 1.0e-3
    end

    @testset "the N=8 space is closed under ⊡ and inv" begin
        A = TensTI{4}(2.0, 1.0, 0.5, 0.7, 0.9, 1.3, 0.4, 0.6, n)
        B = TensTI{4}(3.0, 1.2, 0.3, 0.8, 1.1, 0.5, 0.2, 0.9, n)
        P = A ⊡ B
        @test P isa TensTI{4, Float64, 8}
        @test arr(P) ≈ arr(A) ⊡ arr(B)
        @test arr(A ⊡ inv(A)) ≈ arr(𝕀)

        # The product rule: one 2×2 block and two complex numbers.
        ℓA, ℓB, ℓP = get_ℓ8(A), get_ℓ8(B), get_ℓ8(P)
        LA = [ℓA[1] ℓA[3]; ℓA[4] ℓA[2]]
        LB = [ℓB[1] ℓB[3]; ℓB[4] ℓB[2]]
        @test LA * LB ≈ [ℓP[1] ℓP[3]; ℓP[4] ℓP[2]]
        @test (ℓA[6] + im * ℓA[7]) * (ℓB[6] + im * ℓB[7]) ≈ ℓP[6] + im * ℓP[7]
        @test (ℓA[5] + im * ℓA[8]) * (ℓB[5] + im * ℓB[8]) ≈ ℓP[5] + im * ℓP[8]
    end

    @testset "get_ℓ8 zero-pads N=5 and N=6" begin
        ℓ = get_ℓ8(TensTI{4}(2.0, 1.0, 0.5, 0.9, 1.3, n))
        @test ℓ[3] == ℓ[4]
        @test ℓ[7] == 0 && ℓ[8] == 0
    end

    # ── Kelvin–Mandel ────────────────────────────────────────────────────────
    # docs/src/theory/kelvin_mandel.md

    @testset "Kelvin–Mandel is an isometry" begin
        A = TensISO{3}(3.0, 2.0)
        B = TensISO{3}(5.0, 7.0)
        @test sqrt(sum(arr(A) .^ 2)) ≈ norm(KM(A))
        @test KM(A ⊡ B) ≈ KM(A) * KM(B)
        @test KM(inv(A)) ≈ inv(Matrix(KM(A)))
        @test KM(𝕀) ≈ Matrix(1.0I, 6, 6)
    end

    @testset "the KM rotation is orthogonal, and the congruence holds" begin
        θ, ϕ, ψ = 0.37, 0.83, 0.29
        Q = TensND._KM_rotation(θ, ϕ, ψ)
        @test Q'Q ≈ I
        @test Q ≈ KM(rot6(θ, ϕ, ψ))
        @test arr(rot6(θ, ϕ, ψ)) ≈ arr(Tens(Matrix(rot3(θ, ϕ, ψ))) ⊠ˢ Tens(Matrix(rot3(θ, ϕ, ψ))))

        # Mat(t) = Q · Mat_material(t) · Qᵀ, checked on a ROTATED frame: in the
        # canonical frame Q is the identity and every convention error hides.
        t = TensOrtho(10.0, 8.0, 9.0, 3.0, 2.0, 4.0, 2.5, 3.0, 1.5, Basis(θ, ϕ, ψ))
        Qc = TensND._km_congruence(frame(t))
        @test Qc'Qc ≈ I
        @test KM(t) ≈ Qc * KM_material(t) * transpose(Qc)
        @test !(Matrix(Qc) ≈ Matrix(1.0I, 6, 6))    # the frame really is rotated
    end

    # ── Order-4 identities ───────────────────────────────────────────────────
    # docs/src/theory/tensor_algebra.md

    @testset "𝟙 = 1⊠1 and 𝕀 = 1⊠ˢ1" begin
        𝟏 = Matrix(1.0I, 3, 3)
        m = [1.0 0.4 0.2; 0.1 2.0 0.3; 0.5 0.6 3.0]     # not symmetric
        @test (𝟏 ⊠ 𝟏) ⊡ m ≈ m                            # identity of ALL order-2
        @test (𝟏 ⊠ˢ 𝟏) ⊡ m ≈ (m + m') / 2                # identity of the SYMMETRIC ones
        @test arr(tens_Id4(Val(3), Val(Float64))) ≈ 𝟏 ⊠ˢ 𝟏
    end

    @testset "⊠ˢ does not invert termwise" begin
        𝟏 = Matrix(1.0I, 3, 3)
        Isym = 𝟏 ⊠ˢ 𝟏
        a = [1.0 0.4 0.2; 0.1 2.0 0.3; 0.5 0.6 3.0]
        chk(x, y) = norm((x ⊠ˢ y) ⊡ (inv(x) ⊠ˢ inv(y)) - Isym)
        @test chk(a, a) < 1.0e-10                     # proportional: holds
        @test chk(a, 3a) < 1.0e-10
        @test chk(a, 𝟏) > 1.0e-3                      # b = 1 does NOT help
        d1, d2 = Diagonal([1.0, 2.0, 3.0]), Diagonal([4.0, 5.0, 6.0])
        @test norm(d1 * d2 - d2 * d1) < 1.0e-14        # they commute
        @test chk(Matrix(d1), Matrix(d2)) > 1.0e-3     # yet it still fails
    end

    @testset "𝕁 ⊙ 𝕁 = 1, 𝕂 ⊙ 𝕂 = d(d+1)/2 − 1" begin
        @test 𝕁 ⊙ 𝕁 ≈ 1
        @test 𝕂 ⊙ 𝕂 ≈ 5
        @test 𝕁 ⊙ 𝕂 ≈ 0
        for d in 2:3
            @test tens_K4(Val(d), Val(Float64)) ⊙ tens_K4(Val(d), Val(Float64)) ≈
                d * (d + 1) / 2 - 1
        end
    end

    # ── Projection ───────────────────────────────────────────────────────────
    # docs/src/theory/projection.md

    @testset "the TI projection is the orthogonal projection on the Walpole basis" begin
        # The closed forms of `_project_TI_KM` must equal ⟨A,𝕎ˢᵢ⟩ / gᵢ.
        C = [
            10.0 3.0 2.5 0 0 0; 3.0 10.0 2.5 0 0 0; 2.5 2.5 12.0 0 0 0
            0 0 0 4.0 0 0; 0 0 0 0 4.0 0; 0 0 0 0 0 7.0
        ]
        A = arr(inv_KM(C))
        g = [1.0, 1.0, 2.0, 2.0, 2.0]
        by_hand = ntuple(i -> sum(A .* arr(Ws[i])) / g[i], 5)
        @test collect(TensND._project_TI_KM(C)) ≈ collect(by_hand)
    end

    @testset "proj_tens round trips and returns (B, d, drel)" begin
        nv = [0.0, 0.0, 1.0]
        C = tens_TI(10.0, 3.0, 2.5, 12.0, 2.0, nv)
        B, d, drel = proj_tens(:TI, arr(C), nv)
        @test B isa TensTI{4, Float64, 5}
        @test drel < 1.0e-12
        @test drel ≈ d / sqrt(sum(arr(C) .^ 2))
    end

    @testset "ForwardDiff passes through a projection" begin
        # Regression: `ForwardDiff.Dual` is not `<: AbstractFloat`, so before
        # `ApproxType` existed a few ulp of round-off made a Dual-valued tensor
        # look non-minor-symmetric, `_KM_of_array` built a 9×9 matrix and this
        # threw a DimensionMismatch.
        nf = [0.0, 0.0, 1.0]
        f(α) = get_data(
            proj_tens(:TI, arr(tens_TI(10.0, 3.0, 2.5, 12.0, 2.0, [sin(α), 0.0, cos(α)])), nf)[1]
        )[1]
        for α in (0.3, 0.7)
            fd = (f(α + 1.0e-6) - f(α - 1.0e-6)) / 2.0e-6
            @test ForwardDiff.derivative(f, α) ≈ fd rtol = 1.0e-5
        end

        # And through an ORTHO projection, where the material frame used to be
        # converted to the field's element type and produced NaN derivatives at
        # the gimbal-lock point of the canonical frame.
        for fr in (CanonicalBasis{3, Float64}(), Basis(0.3, 0.8, 0.2))
            g(α) = get_data(
                proj_tens(:ORTHO, arr(tens_TI(10.0, 3.0, 2.5, 12.0, 2.0, [sin(α), 0.0, cos(α)])), fr)[1]
            )[1]
            ad = ForwardDiff.derivative(g, 0.5)
            @test isfinite(ad)
            @test ad ≈ (g(0.5 + 1.0e-6) - g(0.5 - 1.0e-6)) / 2.0e-6 rtol = 1.0e-5
        end
    end

    # ── Curvilinear differential calculus ────────────────────────────────────
    # docs/src/theory/curvilinear.md

    @testset "Γᵏᵢⱼ = ∂ᵢ𝐚ⱼ ⋅ 𝐚ᵏ, with Γ[i,j,k] = Γᵏᵢⱼ" begin
        CS = coorsys_spherical()
        coords = getcoords(CS)
        Γ = Christoffel(CS)
        for i in 1:3, j in 1:3, k in 1:3
            direct = ∂(natvec(CS, j, :cov), coords[i], CS) ⋅ natvec(CS, k, :cont)
            @test iszero(tsimplify(direct - Γ[i, j, k]))
        end
        # Symmetric in the lower pair.
        @test all(iszero, tsimplify.(Γ - permutedims(Γ, (2, 1, 3))))
    end

    @testset "operator index placement" begin
        # GRAD appends the derivative index on the RIGHT: (∇v)ᵢⱼ = ∂ⱼvᵢ.
        # DIV contracts the LAST index:                  (div σ)ᵢ = ∂ⱼσᵢⱼ.
        # Both are invisible on symmetric fields, hence the deliberately
        # non-symmetric one below.
        Cartesian = coorsys_cartesian()
        X = getcoords(Cartesian)
        E = unitvec(Cartesian)
        v = sum(SymFunction("v$i", real = true)(X...) * E[i] for i in 1:3)
        Gv = get_array(GRAD(v, Cartesian))
        for i in 1:3, j in 1:3
            @test iszero(tsimplify(Gv[i, j] - diff(SymFunction("v$i", real = true)(X...), X[j])))
        end

        t = Tens(Tensor{2, 3}((i, j) -> SymFunction("t$i$j", real = true)(X...)))
        dv = get_array(DIV(t, Cartesian))
        for i in 1:3
            expected = sum(diff(SymFunction("t$i$j", real = true)(X...), X[j]) for j in 1:3)
            @test iszero(tsimplify(dv[i] - expected))
        end
    end

    @testset "non-orthogonal charts" begin
        # Every predefined system is orthogonal, and on an orthogonal chart the
        # covariant and contravariant natural vectors coincide — so a swap
        # between the two is invisible there. It was not invisible here: with
        # the variances exchanged, `natvec(CS, i, :cov)` returned the *dual*
        # vector instead of ∂OM/∂qⁱ and the Laplacian of a harmonic function
        # came out nonzero on any skew chart.
        u, v = symbols("u v", real = true)
        CS = CoorSystemSym(Tens([u + v^2, v]), (u, v))          # x = u+v², y = v

        @test !isorthonormal(normalized_basis(CS))               # genuinely skew

        # aᵢ = ∂OM/∂qⁱ
        @test all(iszero, tsimplify(components_canon(natvec(CS, 1, :cov)) - [Sym(1), Sym(0)]))
        @test all(iszero, tsimplify(components_canon(natvec(CS, 2, :cov)) - [2v, Sym(1)]))
        # duality aᵢ ⋅ aʲ = δ
        for i in 1:2, j in 1:2
            @test iszero(tsimplify(natvec(CS, i, :cov) ⋅ natvec(CS, j, :cont) - (i == j)))
        end
        # the induced metric
        @test all(iszero, tsimplify(metric(natural_basis(CS), :cov) - [1 2v; 2v 1 + 4v^2]))

        # x² − y² and x·y are harmonic; expressed in (u, v) they must stay so.
        @test iszero(simplify(LAPLACE((u + v^2)^2 - v^2, CS)))
        @test iszero(simplify(LAPLACE((u + v^2) * v, CS)))
    end

    @testset "Lamé coefficients of the predefined systems" begin
        @test all(isone, Lame(coorsys_cartesian()))
        r, θ = getcoords(coorsys_polar())
        @test collect(Lame(coorsys_polar())) == [1, r]
        θs, ϕs, rs = getcoords(coorsys_spherical())
        @test collect(Lame(coorsys_spherical())) == [rs, rs * sin(θs), 1]
    end

    @testset "classical operator identities" begin
        # The compact assertion form of the `operator_identities` tutorial. Each
        # identity holds trivially in Cartesian coordinates, tests the
        # Christoffel symbols on a curvilinear chart, and tests the
        # covariant/contravariant bookkeeping on the skew one — which is the
        # only chart where the natural basis and its dual differ.
        uu, vv = symbols("u v", real = true)
        for CS in (
                coorsys_cartesian(),
                coorsys_polar(),
                coorsys_spherical(),
                CoorSystemSym(Tens([uu + vv^2, vv]), (uu, vv)),   # non-orthogonal
            )
            d = get_dim(CS)
            q = getcoords(CS)
            f = SymFunction("f", real = true)(q...)
            g = SymFunction("g", real = true)(q...)
            𝐯 = sum(SymFunction("v$i", real = true)(q...) * unitvec(CS, i) for i in 1:d)
            𝟏 = tens_Id2(Val(d), Val(Sym))

            zero_res(x) = iszero(tsimplify(maximum(abs, tsimplify.(get_array(x)))))
            zero_sc(x) = iszero(tsimplify(x))

            @test zero_sc(LAPLACE(f, CS) - DIV(GRAD(f, CS), CS))
            @test zero_sc(LAPLACE(f, CS) - tr(HESS(f, CS)))
            @test zero_res(GRAD(f * g, CS) - (f * GRAD(g, CS) + g * GRAD(f, CS)))
            @test zero_sc(
                LAPLACE(f * g, CS) -
                    (f * LAPLACE(g, CS) + 2 * (GRAD(f, CS) ⋅ GRAD(g, CS)) + g * LAPLACE(f, CS))
            )
            @test zero_sc(DIV(f * 𝐯, CS) - (GRAD(f, CS) ⋅ 𝐯 + f * DIV(𝐯, CS)))
            @test zero_res(DIV(f * 𝟏, CS) - GRAD(f, CS))
            @test zero_sc(tr(GRAD(𝐯, CS)) - DIV(𝐯, CS))
            # `∂𝟏 = 0` is metric compatibility of the connection.
            @test zero_res(∂(𝟏, 1, CS))
            # Leibniz: `∂` is a derivation.
            𝐚 = sum(SymFunction("a$i", real = true)(q...) * unitvec(CS, i) for i in 1:d)
            @test zero_res(∂(𝐚 ⊗ 𝐯, 1, CS) - (∂(𝐚, 1, CS) ⊗ 𝐯 + 𝐚 ⊗ ∂(𝐯, 1, CS)))
        end
    end

    @testset "harmonic functions are annihilated" begin
        # The single strongest check on a coordinate system: it fails unless
        # every connection term is right.
        Polar = coorsys_polar()
        r, θ = getcoords(Polar)
        nn = symbols("n", integer = true)
        @test iszero(simplify(LAPLACE(r^nn * cos(nn * θ), Polar)))

        Spherical = coorsys_spherical()
        θs, ϕs, rs = getcoords(Spherical)
        @test iszero(tsimplify(LAPLACE(1 / rs, Spherical)))
        @test iszero(tsimplify(tr(HESS(1 / rs, Spherical)) - LAPLACE(1 / rs, Spherical)))
    end

    @testset "symbolic and numerical operators agree" begin
        CS = coorsys_spherical()
        θs, ϕs, rs = getcoords(CS)
        Γ = Christoffel(CS)
        x₀ = [0.7, 1.1, 2.3]
        Γnum = Christoffel(coorsys_spherical_num(), x₀)
        Γsym = [
            Float64(tsimplify(Γ[i, j, k]).subs(Dict(θs => x₀[1], ϕs => x₀[2], rs => x₀[3])))
                for i in 1:3, j in 1:3, k in 1:3
        ]
        @test Γnum ≈ Γsym
        @test collect(Lame(coorsys_spherical_num(), x₀)) ≈ [x₀[3], x₀[3] * sin(x₀[1]), 1.0]
    end

    # ── Submanifolds ─────────────────────────────────────────────────────────
    # docs/src/theory/submanifolds.md

    @testset "GRAD(n) = −b on a sphere, as tensors" begin
        R = symbols("R", positive = true)
        θ = symbols("θ", positive = true)
        ϕ = symbols("ϕ", real = true)
        Sphere = SubManifoldSym(
            Tens(R * [sin(θ)cos(ϕ), sin(θ)sin(ϕ), cos(θ)]), (θ, ϕ), (), (R,);
            rules = Dict(abs(sin(θ)) => sin(θ)),
        )
        b = curvature(Sphere)
        G = GRAD(normal(Sphere), Sphere)
        ℬ = natural_basis(Sphere)
        # The identity holds between TENSORS: it must be checked in a common
        # basis and variance, never by subtracting `get_array`s.
        for var in ((:cov, :cov), (:cont, :cov), (:cov, :cont), (:cont, :cont))
            @test all(iszero, tsimplify.(components(G, ℬ, var) + components(b, ℬ, var)))
        end
        # b = −a/R for a sphere with the outward normal.
        @test all(iszero, tsimplify.(get_array(b) + get_array(submetric(Sphere)) / R))
    end

end
