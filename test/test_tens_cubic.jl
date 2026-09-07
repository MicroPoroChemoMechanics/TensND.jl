@testsection "Cubic tensors" begin

    can = CanonicalBasis{3, Float64}()
    rot = RotatedBasis(0.3, 0.4, 0.5)
    atol_num = 1.0e-12

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensCubic — construction & traits" begin
        C = tens_cubic(10.0, 4.0, 2.0, can)
        @test C isa TensCubic{Float64}
        @test size(C) == (3, 3, 3, 3)
        @test length(C) == 81
        @test get_basis(C) isa CanonicalBasis{3, Float64}
        @test get_var(C) == (:cont, :cont, :cont, :cont)
        @test frame(C) === can
        @test reference(C) === can
        @test symmetry(C) == :CUBIC
        @test is_CUBIC(C)
        @test !is_ISO(C)
        @test !is_TI(C)
        @test !is_ORTHO(C)
        @test !is_CUBIC(TensISO{3}(1.0, 1.0))

        # (α, β, γ) = (C₁₁+2C₁₂, C₁₁-C₁₂, 2C₄₄)
        @test get_data(C) == (18.0, 6.0, 4.0)
        @test arg_cubic(C) == (10.0, 4.0, 2.0)
        @test all(arg_cubic(tens_cubic(arg_cubic(C)..., can)) .≈ arg_cubic(C))

        # From a Kelvin-Mandel matrix in the cube frame.
        @test arg_cubic(TensCubic(KM_material(C), can)) == arg_cubic(C)

        # Element types promote, and the frame's is independent of the data's.
        @test tens_cubic(10, 4, 2, can) isa TensCubic{Float64}
        d = ForwardDiff.Dual{:t}(10.0, 1.0)
        Cd = tens_cubic(d, 4.0, 2.0, can)
        @test eltype(Cd) <: ForwardDiff.Dual
        @test eltype(vecbasis(frame(Cd), :cov)) === Float64

        @test occursin("𝕁", sprint(show, C))
        @test occursin("cube frame", sprint(show, C))
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "Kelvin-Mandel and components" begin
        C = tens_cubic(10.0, 4.0, 2.0, can)
        M = Matrix(KM_material(C))
        @test M[1, 1] == M[2, 2] == M[3, 3] == 10.0
        @test M[1, 2] == M[1, 3] == M[2, 3] == 4.0
        @test M[4, 4] == M[5, 5] == M[6, 6] == 4.0      # 2C₄₄
        @test all(iszero, M[1:3, 4:6])
        # On the canonical frame the two Kelvin-Mandel forms coincide.
        @test Matrix(KM(C)) ≈ M atol = atol_num

        # On a rotated cube frame they are related by the Kelvin-Mandel
        # congruence of the frame — the same identity `TensOrtho` asserts.
        Cr = tens_cubic(10.0, 4.0, 2.0, rot)
        Q = Matrix(TensND._km_congruence(Matrix(vecbasis(rot, :cov))))
        @test Matrix(KM(Cr)) ≈ Q * Matrix(KM_material(Cr)) * transpose(Q) atol = atol_num
        # ...and they are genuinely different, or the test above says nothing.
        @test norm(Matrix(KM(Cr)) - Matrix(KM(C))) > 0.5

        # `getindex` is the closed form `get_array` uses, component by component.
        A = get_array(Cr)
        @test size(A) == (3, 3, 3, 3)
        @test all(A[i, j, k, l] ≈ Cr[i, j, k, l] for i in 1:3, j in 1:3, k in 1:3, l in 1:3)
        @test Array(Cr) ≈ A atol = atol_num
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "Symmetries — major symmetry is automatic" begin
        # Under the octahedral group the Kelvin-Mandel space splits as
        # A1g + Eg + T2g with multiplicity one, so the commutant is spanned by
        # three *symmetric* projectors and every cubic tensor is
        # major-symmetric — whatever the coefficients, and on any frame.
        for fr in (can, rot), (a, b, c) in ((18.0, 6.0, 4.0), (1.0, -3.0, 7.0))
            C = TensCubic(a, b, c, fr)
            @test issymmetric(C)
            @test isminorsymmetric(C)
            @test ismajorsymmetric(C)
            K = Matrix(KM(C))
            @test norm(K - transpose(K)) < atol_num * max(norm(K), 1)
        end
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "Algebra — componentwise on the three projectors" begin
        A = TensCubic(18.0, 6.0, 4.0, can)
        B = TensCubic(3.0, -1.0, 2.0, can)

        @test get_data(A + B) == (21.0, 5.0, 6.0)
        @test get_data(A - B) == (15.0, 7.0, 2.0)
        @test get_data(2A) == (36.0, 12.0, 8.0)
        @test get_data(A * 2) == (36.0, 12.0, 8.0)
        @test get_data(A / 2) == (9.0, 3.0, 2.0)
        @test get_data(-A) == (-18.0, -6.0, -4.0)
        @test (A + B) isa TensCubic

        # The inverse is componentwise, and it really is an inverse.
        @test get_data(inv(A)) == (1 / 18, 1 / 6, 1 / 4)
        @test Matrix(KM(A ⊡ inv(A))) ≈ Matrix(1.0I, 6, 6) atol = atol_num
        @test get_data(A^-1) == get_data(inv(A))
        @test get_data(A^2) == get_data(A) .^ 2
        @test get_data(A^1) == get_data(A)
        @test get_data(A^0) == (1.0, 1.0, 1.0)
        @test transpose(A) === A
        @test adjoint(A) === A

        # `one` is the symmetric identity; `zero` stays in the class and on the
        # frame rather than narrowing to a TensISO, which would silently change
        # the type of an accumulator.
        @test Matrix(KM(one(A))) ≈ Matrix(1.0I, 6, 6) atol = atol_num
        @test zero(A) isa TensCubic
        @test frame(zero(A)) === can
        @test all(iszero, get_data(zero(A)))
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "The product stays in the class" begin
        # Neither `TensTI{4}` (which widens from N=5 to N=6) nor `TensOrtho`
        # (which leaves its class) manages this.
        A = TensCubic(18.0, 6.0, 4.0, rot)
        B = TensCubic(3.0, 2.0, 5.0, rot)
        P = A ⊡ B
        @test P isa TensCubic
        @test get_data(P) == (54.0, 12.0, 20.0)
        @test frame(P) === rot
        @test ismajorsymmetric(P)
        # ...and it agrees with the dense route it replaces.
        @test Matrix(KM(P)) ≈ Matrix(KM(Tens(get_array(A)) ⊡ Tens(get_array(B)))) atol =
            1.0e-10
        # Cubic tensors on the same cube commute, which is what makes the
        # product stay in the class.
        @test Matrix(KM(A ⊡ B)) ≈ Matrix(KM(B ⊡ A)) atol = atol_num
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "Frames describing the same cube" begin
        C = tens_cubic(10.0, 4.0, 2.0, can)
        # A signed permutation of the axes is a symmetry of the cube, so it
        # describes the *same* tensor with the *same* coefficients. This is
        # specific to cubic symmetry: permuting an orthotropic frame permutes
        # C₁₁, C₂₂, C₃₃.
        perm = Basis([0.0 0.0 1.0; 1.0 0.0 0.0; 0.0 -1.0 0.0])
        Cp = tens_cubic(10.0, 4.0, 2.0, perm)
        @test TensND._same_cube_frame(C, Cp)
        @test Matrix(KM(Cp)) ≈ Matrix(KM(C)) atol = 1.0e-10
        @test (C + Cp) isa TensCubic
        @test get_data(C + Cp) == 2 .* get_data(C)

        # A general rotation is not, and the sum is refused rather than wrong.
        Cr = tens_cubic(10.0, 4.0, 2.0, rot)
        @test !TensND._same_cube_frame(C, Cr)
        @test_throws AssertionError C + Cr
        # The product falls back to the dense route instead of pretending.
        @test !((C ⊡ Cr) isa TensCubic)
        @test Matrix(KM(C ⊡ Cr)) ≈
            Matrix(KM(Tens(get_array(C)) ⊡ Tens(get_array(Cr)))) atol = 1.0e-10
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "Promotion — ISO ⊂ CUBIC ⊂ ORTHO" begin
        iso = TensISO{3}(18.0, 6.0)
        C = tens_cubic(10.0, 4.0, 2.0, can)

        # An isotropic tensor is cubic about *every* cube: α𝕁 + β𝕂 = α𝕁 + β𝔼 + β𝕋.
        Ci = iso_to_cubic(iso, rot)
        @test get_data(Ci) == (18.0, 6.0, 6.0)
        @test Matrix(KM(Ci)) ≈ Matrix(KM(iso)) atol = atol_num
        @test cubic_anisotropy(Ci) ≈ 0 atol = atol_num

        @test (iso + C) isa TensCubic
        @test (C + iso) isa TensCubic
        @test get_data(C + iso) == (36.0, 12.0, 10.0)
        @test Matrix(KM(C - iso)) ≈ Matrix(KM(C)) - Matrix(KM(iso)) atol = atol_num
        @test (iso ⊡ C) isa TensCubic
        @test get_data(iso ⊡ C) == (18.0 * 18.0, 36.0, 24.0)
        @test (C ⊡ iso) isa TensCubic

        # A cubic tensor is orthotropic with the three families equal.
        O = cubic_to_ortho(C)
        @test O isa TensOrtho
        @test get_data(O) == (10.0, 10.0, 10.0, 4.0, 4.0, 4.0, 2.0, 2.0, 2.0)
        @test Matrix(KM(O)) ≈ Matrix(KM(C)) atol = atol_num
        @test (C + O) isa TensOrtho
        @test (O - C) isa TensOrtho
        @test Matrix(KM(C + O)) ≈ 2 * Matrix(KM(C)) atol = atol_num

        # Different frames: refused, with the frame named rather than silently
        # producing something wrong.
        @test_throws ArgumentError C + cubic_to_ortho(tens_cubic(1.0, 0.0, 1.0, rot))

        # CUBIC and TI are incomparable — their intersection is the isotropic
        # class — so the sum falls through to the unstructured route on purpose.
        ti = tens_TI(10.0, 3.0, 2.5, 12.0, 2.0, [0.0, 0.0, 1.0])
        @test !((C + ti) isa TensCubic)
        @test Matrix(KM(C + ti)) ≈ Matrix(KM(C)) + Matrix(KM(ti)) atol = 1.0e-10
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "Projection onto the cubic class" begin
        C = tens_cubic(10.0, 4.0, 2.0, rot)
        B, d, drel = proj_tens(:CUBIC, get_array(C), rot)
        @test B isa TensCubic
        @test drel < 1.0e-14
        @test all(isapprox.(get_data(B), get_data(C); atol = 1.0e-10))

        # An isotropic tensor is exactly cubic; a TI one is not.
        _, _, dri = proj_tens(:CUBIC, get_array(TensISO{3}(18.0, 6.0)), can)
        @test dri < 1.0e-14
        ti = tens_TI(10.0, 3.0, 2.5, 12.0, 2.0, [0.0, 0.0, 1.0])
        _, _, drt = proj_tens(:CUBIC, get_array(ti), can)
        @test drt > 1.0e-2

        # The projection is orthogonal: the residual is orthogonal to the class.
        Araw = get_array(ti)
        Bti, _, _ = proj_tens(:CUBIC, Araw, can)
        R = Araw .- get_array(Bti)
        @test abs(sum(R .* get_array(Bti))) < 1.0e-10 * sum(abs2, Araw)
        # ...and idempotent.
        _, d2, _ = proj_tens(:CUBIC, get_array(Bti), can)
        @test d2 < 1.0e-12

        # `best_fit_cubic` is the projection component, and at order two the
        # cubic class *is* the isotropic class.
        @test best_fit_cubic(Tens(get_array(C)), rot) isa TensCubic
        A2 = Tens([2.0 0.3 0.0; 0.3 1.0 0.0; 0.0 0.0 0.5])
        @test best_fit_cubic(A2, can) isa TensISO
        # `KM` of a 2nd-order tensor is the 6-vector, not a matrix.
        @test collect(KM(best_fit_cubic(A2, can))) ≈ collect(KM(best_fit_iso(A2))) atol =
            atol_num

        # The Kelvin-Mandel aliases round-trip.
        @test cubic_params_from_KM(Matrix(KM_material(C))) == (18.0, 6.0, 4.0)
        @test Matrix(KM_from_cubic_params(18.0, 6.0, 4.0)) ≈ Matrix(KM_material(C)) atol =
            atol_num

        # A zero tensor projects to zero without dividing by it.
        Bz, dz, drz = proj_tens(:CUBIC, zeros(3, 3, 3, 3), can)
        @test all(iszero, get_data(Bz))
        @test dz == 0
        @test drz == 0
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "Rotational average" begin
        C = tens_cubic(10.0, 4.0, 2.0, can)
        a = isotropify(C)
        @test a isa TensISO{4, 3}
        # α𝕁 + ((2β + 3γ)/5)𝕂, the weighted mean over the dimensions of Eg and T2g.
        @test get_data(a) == (18.0, (2 * 6.0 + 3 * 4.0) / 5)
        # Exactly the generic route, which goes through all 81 components.
        @test Matrix(KM(a)) ≈ Matrix(KM(isotropify(Tens(get_array(C))))) atol = atol_num
        # The average cannot depend on where the cube points.
        Cr = tens_cubic(10.0, 4.0, 2.0, rot)
        @test Matrix(KM(isotropify(Cr))) ≈ Matrix(KM(a)) atol = 1.0e-10
        # A cubic tensor that is already isotropic is its own average.
        Ci = tens_cubic(10.0, 4.0, 3.0, can)          # C₁₁ = C₁₂ + 2C₄₄
        @test cubic_anisotropy(Ci) ≈ 0 atol = atol_num
        @test Matrix(KM(isotropify(Ci))) ≈ Matrix(KM(Ci)) atol = 1.0e-10
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "cubic_anisotropy" begin
        @test cubic_anisotropy(tens_cubic(10.0, 4.0, 3.0, can)) ≈ 0 atol = atol_num
        @test cubic_anisotropy(tens_cubic(10.0, 4.0, 2.0, can)) ≈ (10 - 4 - 4) / 10
        @test cubic_anisotropy(tens_cubic(10.0, 4.0, 4.0, can)) ≈ (10 - 4 - 8) / 10
        # It is what the isotropic average discards: zero anisotropy ⟺ the
        # average changes nothing.
        for (C11, C12, C44) in ((10.0, 4.0, 3.0), (10.0, 4.0, 2.0), (5.0, 1.0, 4.0))
            C = tens_cubic(C11, C12, C44, can)
            same = norm(Matrix(KM(isotropify(C))) - Matrix(KM(C))) < 1.0e-10
            @test same == (abs(cubic_anisotropy(C)) < 1.0e-12)
        end
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "ForwardDiff through the class" begin
        # d/dC₄₄ of the inverse's γ coefficient: γ = 2C₄₄, so 1/γ has slope
        # -2/γ² = -1/(2C₄₄²).
        f = C44 -> get_data(inv(tens_cubic(10.0, 4.0, C44, can)))[3]
        @test ForwardDiff.derivative(f, 2.0) ≈ -1 / (2 * 2.0^2)
        # Through a product and the Kelvin-Mandel matrix.
        g = C11 -> Matrix(
            KM(
                tens_cubic(C11, 4.0, 2.0, can) ⊡
                    tens_cubic(1.0, 0.5, 0.25, can)
            )
        )[1, 1]
        @test isfinite(ForwardDiff.derivative(g, 10.0))
        @test ForwardDiff.derivative(g, 10.0) != 0
        # Through the projection, on a plain Float64 frame — the case that used
        # to return NaN for ORTHO by forcing the frame to the data's type.
        h = C11 -> get_data(proj_tens(Val(:CUBIC), get_array(tens_cubic(C11, 4.0, 2.0, can)), can)[1])[1]
        @test ForwardDiff.derivative(h, 10.0) ≈ 1        # α = C₁₁ + 2C₁₂
        # And through the rotational average.
        k = C44 -> get_data(isotropify(tens_cubic(10.0, 4.0, C44, can)))[2]
        @test ForwardDiff.derivative(k, 2.0) ≈ 6 / 5     # (2β + 3γ)/5, γ = 2C₄₄
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "Symbolic element type" begin
        @syms C11::positive C12::real C44::positive
        Cs = tens_cubic(C11, C12, C44, can)
        # `Sym` is parametric (`Sym{PyObject}`), so match on the eltype.
        @test Cs isa TensCubic
        @test eltype(Cs) <: Sym
        @test eltype(vecbasis(frame(Cs), :cov)) === Float64
        # The shared symbolic operations reach the data through `_rebuild`.
        Cn = tsubs(Cs, C11 => 10, C12 => 4, C44 => 2)
        @test all(isapprox.(Float64.(arg_cubic(Cn)), (10.0, 4.0, 2.0); atol = atol_num))
        # The algebra is exact in the symbols.
        @test tsimplify(get_data(Cs ⊡ inv(Cs))[1] - 1) == 0
        @test tsimplify(get_data(Cs ⊡ inv(Cs))[3] - 1) == 0
        # The isotropic average, symbolically.
        @test tsimplify(get_data(isotropify(Cs))[1] - (C11 + 2C12)) == 0
    end
end
