@testsection "Cubic tensors" begin

    can = CanonicalBasis{3, Float64}()
    rot = RotatedBasis(0.3, 0.4, 0.5)
    atol_num = 1.0e-12

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "TensCubic — construction & traits" begin
        C = cubic_from_constants(10.0, 4.0, 2.0, can)
        @test C isa TensCubic{Float64}
        @test size(C) == (3, 3, 3, 3)
        @test length(C) == 81
        @test get_basis(C) isa CanonicalBasis{3, Float64}
        @test get_var(C) == (:cont, :cont, :cont, :cont)
        @test frame(C) === can
        @test reference(C) === can
        @test symmetry(C) == :CUBIC
        @test iscubic(C)
        @test !isisotropic(C)
        @test !istransverselyisotropic(C)
        @test !isorthotropic(C)
        # `iscubic` on a non-`TensCubic` is the *value*-level question, exactly
        # as `isorthotropic` is, so an isotropic tensor answers `true`: it really is
        # cubic, about every cube.
        @test iscubic(TensISO{3}(1.0, 1.0), can)
        @test isorthotropic(TensISO{3}(1.0, 1.0), can)

        # (α, β, γ) = (C₁₁+2C₁₂, C₁₁-C₁₂, 2C₄₄)
        @test get_data(C) == (18.0, 6.0, 4.0)
        @test cubic_constants(C) == (10.0, 4.0, 2.0)
        @test all(cubic_constants(cubic_from_constants(cubic_constants(C)..., can)) .≈ cubic_constants(C))

        # From a Kelvin-Mandel matrix in the cube frame.
        @test cubic_constants(TensCubic(mandel_material(C), can)) == cubic_constants(C)

        # Element types promote, and the frame's is independent of the data's.
        @test cubic_from_constants(10, 4, 2, can) isa TensCubic{Float64}
        d = ForwardDiff.Dual{:t}(10.0, 1.0)
        Cd = cubic_from_constants(d, 4.0, 2.0, can)
        @test eltype(Cd) <: ForwardDiff.Dual
        @test eltype(vecbasis(frame(Cd), :cov)) === Float64

        @test occursin("𝕁", sprint(show, C))
        @test occursin("cube frame", sprint(show, C))
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "Kelvin-Mandel and components" begin
        C = cubic_from_constants(10.0, 4.0, 2.0, can)
        M = Matrix(mandel_material(C))
        @test M[1, 1] == M[2, 2] == M[3, 3] == 10.0
        @test M[1, 2] == M[1, 3] == M[2, 3] == 4.0
        @test M[4, 4] == M[5, 5] == M[6, 6] == 4.0      # 2C₄₄
        @test all(iszero, M[1:3, 4:6])
        # On the canonical frame the two Kelvin-Mandel forms coincide.
        @test Matrix(mandel(C)) ≈ M atol = atol_num

        # On a rotated cube frame they are related by the Kelvin-Mandel
        # congruence of the frame — the same identity `TensOrtho` asserts.
        Cr = cubic_from_constants(10.0, 4.0, 2.0, rot)
        Q = Matrix(TensND._km_congruence(Matrix(vecbasis(rot, :cov))))
        @test Matrix(mandel(Cr)) ≈ Q * Matrix(mandel_material(Cr)) * transpose(Q) atol = atol_num
        # ...and they are genuinely different, or the test above says nothing.
        @test norm(Matrix(mandel(Cr)) - Matrix(mandel(C))) > 0.5

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
            K = Matrix(mandel(C))
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
        @test Matrix(mandel(A ⊡ inv(A))) ≈ Matrix(1.0I, 6, 6) atol = atol_num
        @test get_data(A^-1) == get_data(inv(A))
        @test get_data(A^2) == get_data(A) .^ 2
        @test get_data(A^1) == get_data(A)
        @test get_data(A^0) == (1.0, 1.0, 1.0)
        @test transpose(A) === A
        @test adjoint(A) === A

        # `one` is the symmetric identity; `zero` stays in the class and on the
        # frame rather than narrowing to a TensISO, which would silently change
        # the type of an accumulator.
        @test Matrix(mandel(one(A))) ≈ Matrix(1.0I, 6, 6) atol = atol_num
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
        @test Matrix(mandel(P)) ≈ Matrix(mandel(Tens(get_array(A)) ⊡ Tens(get_array(B)))) atol =
            1.0e-10
        # Cubic tensors on the same cube commute, which is what makes the
        # product stay in the class.
        @test Matrix(mandel(A ⊡ B)) ≈ Matrix(mandel(B ⊡ A)) atol = atol_num
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "Frames describing the same cube" begin
        C = cubic_from_constants(10.0, 4.0, 2.0, can)
        # A signed permutation of the axes is a symmetry of the cube, so it
        # describes the *same* tensor with the *same* coefficients. This is
        # specific to cubic symmetry: permuting an orthotropic frame permutes
        # C₁₁, C₂₂, C₃₃.
        perm = Basis([0.0 0.0 1.0; 1.0 0.0 0.0; 0.0 -1.0 0.0])
        Cp = cubic_from_constants(10.0, 4.0, 2.0, perm)
        @test TensND._same_cube_frame(C, Cp)
        @test Matrix(mandel(Cp)) ≈ Matrix(mandel(C)) atol = 1.0e-10
        @test (C + Cp) isa TensCubic
        @test get_data(C + Cp) == 2 .* get_data(C)

        # A general rotation is not, and the sum takes the generic route: the
        # sum of the canonical arrays, which no cubic tensor holds.
        Cr = cubic_from_constants(10.0, 4.0, 2.0, rot)
        @test !TensND._same_cube_frame(C, Cr)
        @test !((C + Cr) isa TensCubic) && canonical_components(C + Cr) ≈ canonical_components(C) + canonical_components(Cr)
        @test canonical_components(C - Cr) ≈ canonical_components(C) - canonical_components(Cr)
        # The product falls back to the dense route instead of pretending.
        @test !((C ⊡ Cr) isa TensCubic)
        @test Matrix(mandel(C ⊡ Cr)) ≈
            Matrix(mandel(Tens(get_array(C)) ⊡ Tens(get_array(Cr)))) atol = 1.0e-10

        # The hard case for the test, and the reason it is written with two
        # thresholds rather than one: a frame that is *almost* a signed
        # permutation. Rotate by an angle small enough that every large entry
        # of FᴬᵀFᴮ passes the `≈ 1` test -- cos(10⁻⁶) differs from one by
        # 5·10⁻¹³, well inside the 10⁻¹⁰ tolerance -- and what is left to
        # notice is the off-diagonal residue alone. The check must reject it:
        # two cube frames are the same cube or they are not, and there is no
        # meaningful sense in which a cube tilted by 10⁻⁶ rad is the first one.
        θ = 1.0e-6
        near = Basis([cos(θ) -sin(θ) 0.0; sin(θ) cos(θ) 0.0; 0.0 0.0 1.0])
        Cn = cubic_from_constants(10.0, 4.0, 2.0, near)
        @test !TensND._same_cube_frame(C, Cn)
        @test !((C + Cn) isa TensCubic) && canonical_components(C + Cn) ≈ canonical_components(C) + canonical_components(Cn)
        @test !((C ⊡ Cn) isa TensCubic)
        # And the fallback is not an approximation of the structured route: it
        # is the dense answer, which differs from `2 .* get_data(C)` here.
        @test Matrix(mandel(C ⊡ Cn)) ≈
            Matrix(mandel(Tens(get_array(C)) ⊡ Tens(get_array(Cn)))) atol = 1.0e-10
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "Promotion — ISO ⊂ CUBIC ⊂ ORTHO" begin
        iso = TensISO{3}(18.0, 6.0)
        C = cubic_from_constants(10.0, 4.0, 2.0, can)

        # An isotropic tensor is cubic about *every* cube: α𝕁 + β𝕂 = α𝕁 + β𝔼 + β𝕋.
        Ci = iso_to_cubic(iso, rot)
        @test get_data(Ci) == (18.0, 6.0, 6.0)
        @test Matrix(mandel(Ci)) ≈ Matrix(mandel(iso)) atol = atol_num
        @test cubic_anisotropy(Ci) ≈ 0 atol = atol_num

        @test (iso + C) isa TensCubic
        @test (C + iso) isa TensCubic
        @test get_data(C + iso) == (36.0, 12.0, 10.0)
        @test Matrix(mandel(C - iso)) ≈ Matrix(mandel(C)) - Matrix(mandel(iso)) atol = atol_num
        @test (iso ⊡ C) isa TensCubic
        @test get_data(iso ⊡ C) == (18.0 * 18.0, 36.0, 24.0)
        @test (C ⊡ iso) isa TensCubic

        # A cubic tensor is orthotropic with the three families equal.
        O = cubic_to_ortho(C)
        @test O isa TensOrtho
        @test get_data(O) == (10.0, 10.0, 10.0, 4.0, 4.0, 4.0, 2.0, 2.0, 2.0)
        @test Matrix(mandel(O)) ≈ Matrix(mandel(C)) atol = atol_num
        @test (C + O) isa TensOrtho
        @test (O - C) isa TensOrtho
        @test Matrix(mandel(C + O)) ≈ 2 * Matrix(mandel(C)) atol = atol_num

        # Different frames: the generic route, the sum of the canonical arrays.
        Or = cubic_to_ortho(cubic_from_constants(1.0, 0.0, 1.0, rot))
        @test !((C + Or) isa TensOrtho) && canonical_components(C + Or) ≈ canonical_components(C) + canonical_components(Or)
        @test canonical_components(Or - C) ≈ canonical_components(Or) - canonical_components(C)

        # CUBIC and TI are incomparable — their intersection is the isotropic
        # class — so the sum falls through to the unstructured route on purpose.
        ti = ti_from_constants(10.0, 3.0, 2.5, 12.0, 2.0, [0.0, 0.0, 1.0])
        @test !((C + ti) isa TensCubic)
        @test Matrix(mandel(C + ti)) ≈ Matrix(mandel(C)) + Matrix(mandel(ti)) atol = 1.0e-10
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "Projection onto the cubic class" begin
        C = cubic_from_constants(10.0, 4.0, 2.0, rot)
        B, d, drel = project(get_array(C), :CUBIC, rot)
        @test B isa TensCubic
        @test drel < 1.0e-14
        @test all(isapprox.(get_data(B), get_data(C); atol = 1.0e-10))

        # An isotropic tensor is exactly cubic; a TI one is not.
        _, _, dri = project(get_array(TensISO{3}(18.0, 6.0)), :CUBIC, can)
        @test dri < 1.0e-14
        ti = ti_from_constants(10.0, 3.0, 2.5, 12.0, 2.0, [0.0, 0.0, 1.0])
        _, _, drt = project(get_array(ti), :CUBIC, can)
        @test drt > 1.0e-2

        # The projection is orthogonal: the residual is orthogonal to the class.
        Araw = get_array(ti)
        Bti, _, _ = project(Araw, :CUBIC, can)
        R = Araw .- get_array(Bti)
        @test abs(sum(R .* get_array(Bti))) < 1.0e-10 * sum(abs2, Araw)
        # ...and idempotent.
        _, d2, _ = project(get_array(Bti), :CUBIC, can)
        @test d2 < 1.0e-12

        # `best_fit_cubic` is the projection component, and at order two the
        # cubic class *is* the isotropic class.
        @test best_fit_cubic(Tens(get_array(C)), rot) isa TensCubic
        A2 = Tens([2.0 0.3 0.0; 0.3 1.0 0.0; 0.0 0.0 0.5])
        @test best_fit_cubic(A2, can) isa TensISO
        # `mandel` of a 2nd-order tensor is the 6-vector, not a matrix.
        @test collect(mandel(best_fit_cubic(A2, can))) ≈ collect(mandel(best_fit_iso(A2))) atol =
            atol_num

        # The Kelvin-Mandel aliases round-trip.
        @test cubic_params_from_mandel(Matrix(mandel_material(C))) == (18.0, 6.0, 4.0)
        @test Matrix(mandel_from_cubic_params(18.0, 6.0, 4.0)) ≈ Matrix(mandel_material(C)) atol =
            atol_num

        # A zero tensor projects to zero without dividing by it.
        Bz, dz, drz = project(zeros(3, 3, 3, 3), :CUBIC, can)
        @test all(iszero, get_data(Bz))
        @test dz == 0
        @test drz == 0
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "Value-level predicates, and why cubic has no eigencandidate" begin
        C = cubic_from_constants(10.0, 4.0, 2.0, rot)
        A = get_array(C)

        # With the frame given: exact, and no optimizer involved.
        @test iscubic(A, rot)
        @test iscubic(C, rot)
        @test !iscubic(A, can)                    # not cubic about *that* cube
        @test iscubic(TensISO{3}(18.0, 6.0), can) # isotropic ⊂ cubic
        @test !iscubic(get_array(ti_from_constants(10.0, 3.0, 2.5, 12.0, 2.0, [0.0, 0.0, 1.0])), can)

        # The reason there is no cheap free-frame route, pinned rather than
        # asserted: **every second-order contraction of a cubic tensor is
        # isotropic**, so no eigenframe carries any information about where the
        # cube points.
        for D in (
                [sum(A[i, i, k, l] for i in 1:3) for k in 1:3, l in 1:3],
                [sum(A[i, k, i, l] for i in 1:3) for k in 1:3, l in 1:3],
            )
            @test norm(D - (tr(D) / 3) * I) < 1.0e-10 * abs(tr(D))
        end
        # ...and the consequence: handing the orthotropic eigencandidate a cubic
        # tensor returns an arbitrary frame. The free-orientation search
        # therefore starts from the screened orientations alone; it is tested
        # in `test_orientation_search.jl`.
        @test project(A, Val(:CUBIC), TensND._candidate_ORTHO_frame(A))[3] > 1.0e-2
        @test project(A, Val(:CUBIC), rot)[3] < 1.0e-12

        @test !istransverselyisotropic(A)
        @test isorthotropic(A, rot)          # cubic ⊂ orthotropic, frame given
        @test !isisotropic(A)
        @test isisotropic(get_array(TensISO{3}(18.0, 6.0)))
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "best_sym_tens reports the class" begin
        # The cascade is ordered by number of constants -- 2, 3, 5, 9 -- and not
        # by inclusion, which it cannot be: CUBIC and TI are incomparable.
        _, _, _, sym = best_sym_tens(TensISO{3}(18.0, 6.0))
        @test sym == :ISO
        _, _, drel, sym = best_sym_tens(Tens(get_array(cubic_from_constants(10.0, 4.0, 2.0, can))))
        @test sym == :CUBIC
        @test drel < 1.0e-8
        _, _, _, sym = best_sym_tens(ti_from_constants(10.0, 3.0, 2.5, 12.0, 2.0, [0.0, 0.0, 1.0]))
        @test sym == :TI
        # A structured `TensCubic` is reported from its **own** frame, with no
        # search: `best_sym_tens` reads `frame(t)` when the container has one.
        # A rotated cubic tensor handed over as a plain array cannot be detected
        # on the cheap path -- see the eigencandidate section above -- and needs
        # `optimize_angles = true`, which is exercised in `test_nlopt_ext.jl`.
        _, _, drel, sym = best_sym_tens(cubic_from_constants(10.0, 4.0, 2.0, rot))
        @test sym == :CUBIC
        @test drel < 1.0e-8

        # With a fixed orientation, the classes the argument cannot serve are
        # dropped rather than raising a `MethodError` from inside `project`:
        # an axis leaves ISO and TI, a frame leaves ISO, CUBIC and ORTHO.
        _, _, _, sym = best_sym_tens(
            Tens(get_array(ti_from_constants(10.0, 3.0, 2.5, 12.0, 2.0, [0.0, 0.0, 1.0]))),
            [0.0, 0.0, 1.0],
        )
        @test sym == :TI
        _, _, _, sym = best_sym_tens(Tens(get_array(cubic_from_constants(10.0, 4.0, 2.0, can))), can)
        @test sym == :CUBIC
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "express / components" begin
        # Stored with `get_basis` canonical and the cube frame kept separately,
        # so re-expressing in another basis is a genuine rotation of the
        # components -- relabeling them would return a different tensor.
        C = cubic_from_constants(10.0, 4.0, 2.0, rot)
        @test components(C) == get_array(C)
        @test components(C, (:cont, :cont, :cont, :cont)) == get_array(C)
        # Into a rotated basis: the components change, and they change the way
        # the generic path changes them.
        b = RotatedBasis(0.7, -0.2, 0.4)
        moved = express(C, b)
        ref = express(Tens(get_array(C)), b)
        @test Array(get_array(moved)) ≈ Array(get_array(ref)) atol = 1.0e-12
        @test Array(components(C, b, (:cont, :cont, :cont, :cont))) ≈
            Array(components(Tens(get_array(C)), b, (:cont, :cont, :cont, :cont))) atol =
            1.0e-12
        # The physical tensor is unchanged, which is the point of rotating
        # rather than relabeling — compared in a common basis, since `mandel` of a
        # tensor stored in a rotated basis returns the components *in that
        # basis*.
        @test Matrix(mandel(canonical(moved))) ≈ Matrix(mandel(C)) atol = 1.0e-12
        # And the canonical target is a no-op.
        @test Array(get_array(express(C, can))) ≈ Array(get_array(C)) atol = 1.0e-12
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "Rotational average" begin
        C = cubic_from_constants(10.0, 4.0, 2.0, can)
        a = isotropify(C)
        @test a isa TensISO{4, 3}
        # α𝕁 + ((2β + 3γ)/5)𝕂, the weighted mean over the dimensions of Eg and T2g.
        @test get_data(a) == (18.0, (2 * 6.0 + 3 * 4.0) / 5)
        # Exactly the generic route, which goes through all 81 components.
        @test Matrix(mandel(a)) ≈ Matrix(mandel(isotropify(Tens(get_array(C))))) atol = atol_num
        # The average cannot depend on where the cube points.
        Cr = cubic_from_constants(10.0, 4.0, 2.0, rot)
        @test Matrix(mandel(isotropify(Cr))) ≈ Matrix(mandel(a)) atol = 1.0e-10
        # A cubic tensor that is already isotropic is its own average.
        Ci = cubic_from_constants(10.0, 4.0, 3.0, can)          # C₁₁ = C₁₂ + 2C₄₄
        @test cubic_anisotropy(Ci) ≈ 0 atol = atol_num
        @test Matrix(mandel(isotropify(Ci))) ≈ Matrix(mandel(Ci)) atol = 1.0e-10
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "cubic_anisotropy" begin
        @test cubic_anisotropy(cubic_from_constants(10.0, 4.0, 3.0, can)) ≈ 0 atol = atol_num
        @test cubic_anisotropy(cubic_from_constants(10.0, 4.0, 2.0, can)) ≈ (10 - 4 - 4) / 10
        @test cubic_anisotropy(cubic_from_constants(10.0, 4.0, 4.0, can)) ≈ (10 - 4 - 8) / 10
        # It is what the isotropic average discards: zero anisotropy ⟺ the
        # average changes nothing.
        for (C11, C12, C44) in ((10.0, 4.0, 3.0), (10.0, 4.0, 2.0), (5.0, 1.0, 4.0))
            C = cubic_from_constants(C11, C12, C44, can)
            same = norm(Matrix(mandel(isotropify(C))) - Matrix(mandel(C))) < 1.0e-10
            @test same == (abs(cubic_anisotropy(C)) < 1.0e-12)
        end
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "ForwardDiff through the class" begin
        # d/dC₄₄ of the inverse's γ coefficient: γ = 2C₄₄, so 1/γ has slope
        # -2/γ² = -1/(2C₄₄²).
        f = C44 -> get_data(inv(cubic_from_constants(10.0, 4.0, C44, can)))[3]
        @test ForwardDiff.derivative(f, 2.0) ≈ -1 / (2 * 2.0^2)
        # Through a product and the Kelvin-Mandel matrix.
        g = C11 -> Matrix(
            mandel(
                cubic_from_constants(C11, 4.0, 2.0, can) ⊡
                    cubic_from_constants(1.0, 0.5, 0.25, can)
            )
        )[1, 1]
        @test isfinite(ForwardDiff.derivative(g, 10.0))
        @test ForwardDiff.derivative(g, 10.0) != 0
        # Through the projection, on a plain Float64 frame — the case that used
        # to return NaN for ORTHO by forcing the frame to the data's type.
        h = C11 -> get_data(project(get_array(cubic_from_constants(C11, 4.0, 2.0, can)), Val(:CUBIC), can)[1])[1]
        @test ForwardDiff.derivative(h, 10.0) ≈ 1        # α = C₁₁ + 2C₁₂
        # And through the rotational average.
        k = C44 -> get_data(isotropify(cubic_from_constants(10.0, 4.0, C44, can)))[2]
        @test ForwardDiff.derivative(k, 2.0) ≈ 6 / 5     # (2β + 3γ)/5, γ = 2C₄₄
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testsection "Symbolic element type" begin
        @syms C11::positive C12::real C44::positive
        Cs = cubic_from_constants(C11, C12, C44, can)
        # `Sym` is parametric (`Sym{Py}`), so match on the eltype.
        @test Cs isa TensCubic
        @test eltype(Cs) <: Sym
        @test eltype(vecbasis(frame(Cs), :cov)) === Float64
        # The shared symbolic operations reach the data through `_rebuild`.
        Cn = tsubs(Cs, C11 => 10, C12 => 4, C44 => 2)
        @test all(isapprox.(Float64.(cubic_constants(Cn)), (10.0, 4.0, 2.0); atol = atol_num))
        # The algebra is exact in the symbols.
        @test tsimplify(get_data(Cs ⊡ inv(Cs))[1] - 1) == 0
        @test tsimplify(get_data(Cs ⊡ inv(Cs))[3] - 1) == 0
        # The isotropic average, symbolically.
        @test tsimplify(get_data(isotropify(Cs))[1] - (C11 + 2C12)) == 0
    end
end
