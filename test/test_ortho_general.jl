# The orthotropic shapes beyond `TensOrtho{4, T, 9}`: `TensOrtho{4, T, 12}`, an
# orthotropic tensor without major symmetry (a normal Kelvin-Mandel block of any
# 3×3 and three shear moduli), and `TensOrtho{2, T, 3}`, diagonal in its frame;
# with them the transpose of the TI shapes. Every operation is checked against an independent
# reference: the material-frame Kelvin-Mandel matrix turned into components by
# `from_mandel` and rotated to the canonical frame index by index, or the generic
# route on the canonical components.

# The canonical components of the tensor whose material-frame Kelvin-Mandel
# matrix is `K6`, the frame vectors being the columns of `E`.
function _og_ref(K6, E)
    Am = get_array(from_mandel(Matrix(K6)))
    T = promote_type(eltype(Am), eltype(E))
    R = zeros(T, 3, 3, 3, 3)
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3, a in 1:3, b in 1:3, c in 1:3, d in 1:3
        R[i, j, k, l] += E[i, a] * E[j, b] * E[k, c] * E[l, d] * Am[a, b, c, d]
    end
    return R
end
_og_km(N, s) = [N zeros(3, 3); zeros(3, 3) Diagonal(2 .* collect(s))]
_og_dense(t) = [t[i, j, k, l] for i in 1:3, j in 1:3, k in 1:3, l in 1:3]
_og_rel(a, b) = maximum(abs, a - b) / maximum(abs, b)
_og_frame(ℬ) = Matrix(TensND.vecbasis(ℬ, :cov))

const _OG_N = [160.0 55.0 48.0; 31.0 120.0 42.0; 22.0 17.0 95.0]
const _OG_S = (33.0, 29.0, 27.0)
_og_tens(N, s, ℬ) = TensOrtho(
    N[1, 1], N[2, 2], N[3, 3], N[1, 2], N[1, 3], N[2, 3], s..., N[2, 1], N[3, 1], N[3, 2], ℬ
)

@testsection "Orthotropic shapes" begin
    @testsection "construction and components" begin
        for angles in ((0.0, 0.0, 0.0), (0.7, -0.4, 1.1), (-1.3, 2.2, 0.35))
            ℬ = RotatedBasis(angles...)
            A = _og_tens(_OG_N, _OG_S, ℬ)
            @test A isa TensOrtho{4, Float64, 12, typeof(ℬ)}
            R = _og_ref(_og_km(_OG_N, _OG_S), _og_frame(ℬ))
            @test _og_rel(_og_dense(A), R) < 1.0e-14
            @test _og_rel(get_array(A), R) < 1.0e-14
            @test TensND.mandel_material(A) == _og_km(_OG_N, _OG_S)
            @test _og_rel(mandel(A), mandel(TensND.generic_tens(A))) < 1.0e-14
            # From the Kelvin-Mandel matrix: the whole block with `Val(12)`.
            @test TensOrtho(_og_km(_OG_N, _OG_S), ℬ, Val(12)) == A
            @test TensOrtho(_og_km(_OG_N, _OG_S), ℬ, Val(9)) == TensOrtho(_og_km(_OG_N, _OG_S), ℬ)
            # Change of basis in closed form, against the generic route.
            ℬ2 = RotatedBasis(0.3, 1.2, -0.8)
            @test _og_rel(
                get_array(change_tens(A, ℬ2)),
                get_array(change_tens(TensND.generic_tens(A), ℬ2)),
            ) < 1.0e-13
        end
        # The nine-constant form keeps its type, now with its length.
        O = TensOrtho(1.0, 2.0, 3.0, 0.1, 0.2, 0.3, 4.0, 5.0, 6.0, CanonicalBasis{3, Float64}())
        @test typeof(O) === TensOrtho{4, Float64, 9, CanonicalBasis{3, Float64}}
        @test TensOrtho{4, Float64, 9}(get_data(O), CanonicalBasis{3, Float64}()) === O
        @test_throws "a TensOrtho of order 2 stores 3 coefficients and one of order 4 stores 9 or 12; got order 4 with 10" TensOrtho{4, Float64}(
            ntuple(Float64, 10), CanonicalBasis{3, Float64}()
        )
        # A twelve-constant tensor with a symmetric block is the nine-constant one.
        Os = TensND._ortho_general(O)
        @test Os isa TensOrtho{4, Float64, 12}
        @test _og_dense(Os) == _og_dense(O)
        @test issymmetric(Os) && Tensors.ismajorsymmetric(Os)
        @test !issymmetric(_og_tens(_OG_N, _OG_S, CanonicalBasis{3, Float64}()))
        @test !Tensors.ismajorsymmetric(_og_tens(_OG_N, _OG_S, CanonicalBasis{3, Float64}()))
        @test Tensors.isminorsymmetric(_og_tens(_OG_N, _OG_S, CanonicalBasis{3, Float64}()))
        # `rebuild` with twelve coefficients gives the general form.
        @test rebuild(O, (get_data(O)..., 0.4, 0.5, 0.6)) isa TensOrtho{4, Float64, 12}
        @test length(param_norms(Os)) == 12
        @test norm(Os) ≈ norm(param_norms(Os) .* get_data(Os))
    end

    @testsection "arithmetic" begin
        ℬ = RotatedBasis(0.7, -0.4, 1.1)
        A = _og_tens(_OG_N, _OG_S, ℬ)
        N2 = [11.0 3.0 2.0; 1.0 7.0 1.5; 0.5 0.25 5.0]
        B = _og_tens(N2, (4.0, 3.5, 2.5), ℬ)
        O = TensOrtho(11.0, 7.0, 5.0, 3.0, 2.0, 1.5, 4.0, 3.5, 2.5, ℬ)
        I4 = TensISO{3}(3.0, 2.0)
        Cℬ = TensCubic{Float64}((30.0, 10.0, 8.0), ℬ)
        for (l, r) in ((A, B), (A, O), (O, A), (I4, A), (A, I4), (Cℬ, A), (A, Cℬ))
            @test l + r isa TensOrtho{4, Float64, 12}
            @test _og_rel(_og_dense(l + r), _og_dense(l) + _og_dense(r)) < 1.0e-14
            @test _og_rel(_og_dense(l - r), _og_dense(l) - _og_dense(r)) < 1.0e-14
        end
        @test _og_dense(2 * A) == 2 * _og_dense(A)
        @test _og_dense(A / 2) ≈ _og_dense(A) / 2
        @test _og_dense(-A) == -_og_dense(A)
        @test (A + B) isa TensOrtho{4, Float64, 12}
        @test (O + O) isa TensOrtho{4, Float64, 9}
        @test_throws ArgumentError A + _og_tens(N2, (4.0, 3.5, 2.5), RotatedBasis(0.1, 0.2, 0.3))
        @test (@inferred A + O) isa TensOrtho
        @test (@inferred A - B) isa TensOrtho
    end

    @testsection "double contraction" begin
        ref(A, B) = get_array(Tensors.dcontract(TensND.generic_tens(A), TensND.generic_tens(B)))
        for angles in ((0.0, 0.0, 0.0), (0.7, -0.4, 1.1), (-1.3, 2.2, 0.35))
            ℬ = RotatedBasis(angles...)
            A = _og_tens(_OG_N, _OG_S, ℬ)
            O = TensOrtho(11.0, 7.0, 5.0, 3.0, 2.0, 1.5, 4.0, 3.5, 2.5, ℬ)
            I4 = TensISO{3}(3.0, 2.0)
            Cℬ = TensCubic{Float64}((30.0, 10.0, 8.0), ℬ)
            for (l, r) in ((A, A), (A, O), (O, A), (O, O), (I4, A), (A, I4), (I4, O), (Cℬ, O), (A, Cℬ))
                P = l ⊡ r
                # The product is orthotropic on the same frame, twelve constants:
                # a `TensCanonical` before.
                @test P isa TensOrtho{4, Float64, 12}
                @test TensND.frame(P) == ℬ
                @test _og_rel(_og_dense(P), ref(l, r)) < 1.0e-13
            end
            @test TensND.mandel_material(A ⊡ O) ≈ TensND.mandel_material(A) * TensND.mandel_material(O)
            # Two symmetric blocks that do not commute: no major symmetry.
            P = O ⊡ TensOrtho(160.0, 120.0, 95.0, 55.0, 48.0, 42.0, 33.0, 29.0, 27.0, ℬ)
            p = get_array(P)
            @test maximum(abs, p - permutedims(p, (3, 4, 1, 2))) > 1.0
            @test !Tensors.ismajorsymmetric(P)
            # Same frame or not is a run-time question: a two-member union.
            @test (@inferred Union{TensOrtho, TensND.TensCanonical} A ⊡ O) isa TensOrtho
        end
        # Different frames: the generic route.
        A = _og_tens(_OG_N, _OG_S, RotatedBasis(0.7, -0.4, 1.1))
        B = _og_tens(_OG_N', _OG_S, RotatedBasis(-1.3, 2.2, 0.35))
        @test !(A ⊡ B isa TensOrtho)
        @test _og_dense(A ⊡ B) == ref(A, B)
    end

    @testsection "inverse, transpose, symmetric part" begin
        for angles in ((0.0, 0.0, 0.0), (0.7, -0.4, 1.1))
            ℬ = RotatedBasis(angles...)
            A = _og_tens(_OG_N, _OG_S, ℬ)
            Ai = inv(A)
            @test Ai isa TensOrtho{4, Float64, 12}
            @test TensND.mandel_material(Ai) ≈ inv(_og_km(_OG_N, _OG_S)) rtol = 1.0e-14
            @test _og_dense(A ⊡ Ai) ≈ _og_dense(TensISO{3}(1.0, 1.0)) atol = 1.0e-14
            @test A^-1 == Ai
            @test (@inferred inv(A)) isa TensOrtho
            At = transpose(A)
            @test At isa TensOrtho{4, Float64, 12}
            @test _og_dense(At) ≈ permutedims(_og_dense(A), (3, 4, 1, 2)) rtol = 1.0e-14
            @test transpose(At) == A
            S = Tensors.majorsymmetric(A)
            @test S isa TensOrtho{4, Float64, 9}
            @test _og_dense(S) ≈ (_og_dense(A) + permutedims(_og_dense(A), (3, 4, 1, 2))) / 2 rtol = 1.0e-14
            O = S
            @test transpose(O) === O
            @test Tensors.majorsymmetric(O) === O
            # The orthogonal projection on the nine-constant class is the
            # symmetric part: the skew part is orthogonal to every symmetric tensor.
            @test _og_rel(_og_dense(best_fit_ortho(A, ℬ)), _og_dense(S)) < 1.0e-12
            @test _og_rel(get_array(isotropify(A)), get_array(isotropify(TensND.generic_tens(A)))) < 1.0e-13
        end
        # The adjoint conjugates, then transposes.
        Az = _og_tens(_OG_N .+ im .* _OG_N', _OG_S .+ (0.5im, 0.0im, 1.0im), CanonicalBasis{3, Float64}())
        @test _og_dense(Az') ≈ conj.(permutedims(_og_dense(Az), (3, 4, 1, 2)))
    end

    @testsection "transversely isotropic operands about a frame axis" begin
        ℬ = RotatedBasis(0.7, -0.4, 1.1)
        E = _og_frame(ℬ)
        ref(A, B) = get_array(Tensors.dcontract(TensND.generic_tens(A), TensND.generic_tens(B)))
        A = _og_tens(_OG_N, _OG_S, ℬ)
        for k in 1:3
            n = E[:, k]
            for W in (TensTI{4}(31.0, 39.0, 11.0, 23.0, 21.0, n), TensTI{4}(31.0, 39.0, 11.0, 6.0, 23.0, 21.0, n))
                Wo = TensND.walpole_to_ortho(W, ℬ, k)
                @test Wo isa TensOrtho{4, Float64, length(get_data(W)) == 5 ? 9 : 12}
                @test _og_rel(_og_dense(Wo), _og_dense(W)) < 1.0e-13
                for (l, r) in ((W, A), (A, W))
                    @test l + r isa TensOrtho{4, Float64, 12}
                    @test _og_rel(_og_dense(l + r), _og_dense(l) + _og_dense(r)) < 1.0e-13
                    @test _og_rel(_og_dense(l - r), _og_dense(l) - _og_dense(r)) < 1.0e-13
                    @test l ⊡ r isa TensOrtho{4, Float64, 12}
                    @test _og_rel(_og_dense(l ⊡ r), ref(l, r)) < 1.0e-13
                end
            end
        end
        # An axis off the frame: sums are refused, products take the generic route.
        n = normalize([1.0, 2.0, 2.0])
        for W in (TensTI{4}(31.0, 39.0, 11.0, 23.0, 21.0, n), TensTI{4}(31.0, 39.0, 11.0, 6.0, 23.0, 21.0, n))
            @test_throws "requires the Walpole axis" W + A
            @test_throws "requires the Walpole axis" A - W
            for (l, r) in ((W, A), (A, W))
                @test !(l ⊡ r isa TensOrtho)
                @test _og_rel(_og_dense(l ⊡ r), ref(l, r)) < 1.0e-13
            end
        end
        @test_throws "axis_idx must be 1, 2 or 3" TensND.walpole_to_ortho(TensTI{4}(31.0, 39.0, 11.0, 6.0, 23.0, 21.0, n), ℬ, 4)
    end

    @testsection "display" begin
        A = _og_tens(_OG_N, _OG_S, CanonicalBasis{3, Float64}())
        s = sprint(show, A)
        @test occursin("(31.0) P₂⊗P₁", s) && occursin("(55.0) P₁⊗P₂", s) && occursin("frame:", s)
        @test occursin("(17.0) P₃⊗P₂", s)
        @test occursin("KM (material frame)", sprint(show, MIME"text/plain"(), A))
        printed = mktemp() do path, io
            redirect_stdout(() -> pprint(A), io)
            flush(io)
            read(path, String)
        end
        @test occursin("(31.0) P₂⊗P₁", printed) && occursin("frame:", printed)
    end

    @testsection "derivatives through the twelve-constant algebra" begin
        ℬ = RotatedBasis(0.7, -0.4, 1.1)
        N1 = [1.0 2.0 0.5; -1.0 0.5 0.25; 0.75 1.5 -0.5]
        s1 = (0.5, 1.0, 1.5)
        A(x) = _og_tens(_OG_N .+ x .* N1, _OG_S .+ x .* s1, ℬ)
        O = TensOrtho(11.0, 7.0, 5.0, 3.0, 2.0, 1.5, 4.0, 3.5, 2.5, ℬ)
        K1 = _og_km(N1, s1)
        KO = TensND.mandel_material(O)
        # The product is linear in `x`: its derivative is the product of the
        # derivative, in the material frame.
        dP = ForwardDiff.derivative(x -> Matrix(TensND.mandel_material(A(x) ⊡ O)), 0.3)
        @test dP ≈ K1 * KO rtol = 1.0e-13
        # d²(A⁻¹)/dx² = 2 A⁻¹ A' A⁻¹ A' A⁻¹, nested duals under distinct tags.
        K(x) = _og_km(_OG_N .+ x .* N1, _OG_S .+ x .* s1)
        Ki = inv(K(0.3))
        d2 = ForwardDiff.derivative(
            y -> ForwardDiff.derivative(x -> Matrix(TensND.mandel_material(inv(A(x)))), y), 0.3
        )
        @test d2 ≈ 2 * Ki * K1 * Ki * K1 * Ki rtol = 1.0e-10
        @test (@inferred inv(A(ForwardDiff.Dual(0.3, 1.0)))) isa TensOrtho
    end

    @testsection "symbolic coefficients" begin
        ℬ = CanonicalBasis{3, Float64}()
        x = SymPy.symbols("x", real = true)
        As = _og_tens(_OG_N .+ x .* [1 0 0; 2 0 0; 0 0 1], _OG_S, ℬ)
        O = TensOrtho(11.0, 7.0, 5.0, 3.0, 2.0, 1.5, 4.0, 3.5, 2.5, ℬ)
        An = _og_tens(_OG_N .+ 0.5 .* [1 0 0; 2 0 0; 0 0 1], _OG_S, ℬ)
        sub(t) = [Float64(c.subs(x, 0.5)) for c in get_data(t)]
        @test sub(As ⊡ O) ≈ collect(get_data(An ⊡ O)) rtol = 1.0e-13
        @test sub(inv(As)) ≈ collect(get_data(inv(An))) rtol = 1.0e-12
        @test sub(transpose(As)) ≈ collect(get_data(transpose(An)))
        @test sub(Tensors.majorsymmetric(As)) ≈ collect(get_data(Tensors.majorsymmetric(An)))
        TensND.Symbolics.@variables y
        Ay = _og_tens(_OG_N .+ y .* [1 0 0; 2 0 0; 0 0 1], _OG_S, ℬ)
        val(t) = [TensND.Symbolics.value(TensND.Symbolics.substitute(c, Dict(y => 0.5); fold = Val(true))) for c in get_data(t)]
        @test val(Ay ⊡ O) ≈ collect(get_data(An ⊡ O)) rtol = 1.0e-13
        @test val(inv(Ay)) ≈ collect(get_data(inv(An))) rtol = 1.0e-12
    end

    @testsection "averages over the orthotropic group: D₂ and D₂ₕ agree" begin
        # The rotations of D₂ (π about each frame axis) and D₂ₕ = D₂ × {±𝟏},
        # written in the canonical frame.
        ℬ = RotatedBasis(0.7, -0.4, 1.1)
        E = _og_frame(ℬ)
        D2 = [E * Diagonal(collect(d)) * transpose(E) for d in ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))]
        D2h = vcat(D2, [-Q for Q in D2])
        rot2(A, Q) = Q * A * transpose(Q)
        function rot4(A, Q)
            R = zeros(eltype(A), 3, 3, 3, 3)
            for i in 1:3, j in 1:3, k in 1:3, l in 1:3, a in 1:3, b in 1:3, c in 1:3, d in 1:3
                R[i, j, k, l] += Q[i, a] * Q[j, b] * Q[k, c] * Q[l, d] * A[a, b, c, d]
            end
            return R
        end
        avg(f, A, G) = sum(f(A, Q) for Q in G) / length(G)
        # Order 2: the inversion acts trivially, and the average is diagonal in
        # the frame, hence symmetric even for a non-symmetric tensor.
        A2 = [1.0 2.0 3.0; -4.0 5.0 6.0; 0.5 -7.0 8.0]
        a2 = avg(rot2, A2, D2)
        @test a2 ≈ avg(rot2, A2, D2h)
        M = transpose(E) * a2 * E
        @test M ≈ Diagonal(diag(M)) atol = 1.0e-14
        @test a2 ≈ transpose(a2)
        @test get_array(TensOrtho{2}(diag(M)..., ℬ)) ≈ a2
        # Order 4, a minor-symmetric tensor without major symmetry: the average
        # has a 3×3 normal block and three shear moduli, the twelve constants.
        G = get_array(TensND.generic_tens(Tens(SymmetricTensor{4, 3}((i, j, k, l) -> sin(i + 2j) * cos(3k - l) + i * j / (k + l)))))
        a4 = avg(rot4, Array(G), D2)
        @test a4 ≈ avg(rot4, Array(G), D2h)
        Km = mandel(change_tens(Tens(a4), ℬ))
        @test Km[1:3, 4:6] ≈ zeros(3, 3) atol = 1.0e-13
        @test Km[4:6, 1:3] ≈ zeros(3, 3) atol = 1.0e-13
        @test Km[4:6, 4:6] ≈ Diagonal(diag(Km[4:6, 4:6])) atol = 1.0e-13
        @test maximum(abs, Km[1:3, 1:3] - transpose(Km[1:3, 1:3])) > 0.1
        @test _og_rel(_og_dense(TensOrtho(Km, ℬ, Val(12))), a4) < 1.0e-13
    end

    @testsection "order 2" begin
        for angles in ((0.0, 0.0, 0.0), (0.7, -0.4, 1.1))
            ℬ = RotatedBasis(angles...)
            E = _og_frame(ℬ)
            A = TensOrtho{2}(3.0, 5.0, 7.0, ℬ)
            B = TensOrtho{2}(2.0, -1.0, 0.5, ℬ)
            ref(a) = E * Diagonal(collect(a)) * transpose(E)
            @test A isa TensOrtho{2, Float64, 3}
            @test get_array(A) ≈ ref((3.0, 5.0, 7.0)) rtol = 1.0e-14
            @test _og_rel([A[i, j] for i in 1:3, j in 1:3], ref((3.0, 5.0, 7.0))) < 1.0e-14
            @test mandel(A) ≈ mandel(TensND.generic_tens(A))
            @test tr(A) == 15.0 && issymmetric(A) && transpose(A) === A && Tensors.symmetric(A) === A
            @test !isisotropic(A) && isisotropic(TensOrtho{2}(2.0, 2.0, 2.0, ℬ)) && isorthotropic(A) && symmetry(A) == :ORTHO
            # Algebra on a shared frame.
            @test (A + B) isa TensOrtho{2, Float64, 3} && get_array(A + B) ≈ get_array(A) + get_array(B)
            @test get_array(A - B) ≈ get_array(A) - get_array(B)
            @test get_array(2A) ≈ 2 * get_array(A) && get_array(-A) == -get_array(A)
            @test inv(A) isa TensOrtho{2} && get_array(inv(A)) ≈ inv(get_array(A))
            @test A^-1 == inv(A)
            @test (A ⋅ B) isa TensOrtho{2} && get_array(A ⋅ B) ≈ get_array(A) * get_array(B)
            @test Tensors.dcontract(A, B) ≈ sum(get_array(A) .* get_array(B))
            v = [1.0, -2.0, 0.5]
            @test collect(A ⋅ v) ≈ get_array(A) * v
            @test collect(get_array(A ⋅ Tens(v))) ≈ get_array(A) * v
            # Isotropic and coaxial TI operands.
            I2 = TensISO{3}(4.0)
            for (l, r) in ((I2, A), (A, I2))
                @test l + r isa TensOrtho{2} && get_array(l + r) ≈ get_array(l) + get_array(r)
                @test get_array(l - r) ≈ get_array(l) - get_array(r)
                @test l ⋅ r isa TensOrtho{2} && get_array(l ⋅ r) ≈ get_array(l) * get_array(r)
            end
            for k in 1:3
                W = TensTI{2}(1.5, 9.0, E[:, k])
                @test get_array(TensND.walpole_to_ortho(W, ℬ, k)) ≈ get_array(W)
                for (l, r) in ((W, A), (A, W))
                    @test l + r isa TensOrtho{2} && get_array(l + r) ≈ get_array(l) + get_array(r)
                    @test get_array(l - r) ≈ get_array(l) - get_array(r)
                    @test l ⋅ r isa TensOrtho{2} && get_array(l ⋅ r) ≈ get_array(l) * get_array(r)
                end
            end
            Wo = TensTI{2}(1.5, 9.0, normalize([1.0, 2.0, 2.0]))
            @test !(Wo ⋅ A isa TensOrtho) && get_array(Wo ⋅ A) ≈ get_array(Wo) * get_array(A)
            @test !(A ⋅ Wo isa TensOrtho) && get_array(A ⋅ Wo) ≈ get_array(A) * get_array(Wo)
            @test_throws "requires the TI axis" TensTI{2}(1.5, 9.0, normalize([1.0, 2.0, 2.0])) + A
            @test_throws "requires the TI axis" A - TensTI{2}(1.5, 9.0, normalize([1.0, 2.0, 2.0]))
            @test_throws "axis_idx must be 1, 2 or 3" TensND.walpole_to_ortho(TensTI{2}(1.0, 2.0, E[:, 1]), ℬ, 0)
            # Double contraction with the order-4 classes, against the generic route.
            C = _og_tens(_OG_N, _OG_S, ℬ)
            refd(X, Y) = get_array(Tensors.dcontract(TensND.generic_tens(X), TensND.generic_tens(Y)))
            I4 = TensISO{3}(3.0, 2.0)
            Cu = TensCubic{Float64}((30.0, 10.0, 8.0), ℬ)
            W4 = TensTI{4}(31.0, 39.0, 11.0, 6.0, 23.0, 21.0, E[:, 2])
            for X in (C, Tensors.majorsymmetric(C), I4, Cu, W4)
                @test X ⊡ A isa TensOrtho{2} && _og_rel(get_array(X ⊡ A), refd(X, A)) < 1.0e-13
                @test A ⊡ X isa TensOrtho{2} && _og_rel(get_array(A ⊡ X), refd(A, X)) < 1.0e-13
            end
            @test I2 ⊡ C isa TensOrtho{2} && _og_rel(get_array(I2 ⊡ C), refd(I2, C)) < 1.0e-13
            @test C ⊡ I2 isa TensOrtho{2} && _og_rel(get_array(C ⊡ I2), refd(C, I2)) < 1.0e-13
            # Tensor products: self product major-symmetric, general otherwise.
            @test Tensors.otimes(A) isa TensOrtho{4, Float64, 9}
            @test _og_dense(Tensors.otimes(A)) ≈ _og_dense(A ⊗ A)
            @test A ⊗ B isa TensOrtho{4, Float64, 12}
            @test I2 ⊗ A isa TensOrtho{4, Float64, 12} && A ⊗ I2 isa TensOrtho{4, Float64, 12}
            for (l, r) in ((A, A), (A, B), (I2, A), (A, I2))
                @test _og_rel(_og_dense(l ⊗ r), get_array(Tensors.otimes(TensND.generic_tens(l), TensND.generic_tens(r)))) < 1.0e-14
            end
            # Change of basis in closed form.
            ℬ2 = RotatedBasis(0.3, 1.2, -0.8)
            @test get_array(change_tens(A, ℬ2)) ≈ get_array(change_tens(TensND.generic_tens(A), ℬ2)) rtol = 1.0e-13
            @test components(A, ℬ2, (:cov, :cov)) ≈ get_array(change_tens(TensND.generic_tens(A), ℬ2)) rtol = 1.0e-13
            # The order-2 projection returns the class.
            P, d, _ = proj_tens(Val(:ORTHO), get_array(A) + [0 1 0; 1 0 0; 0 0 0], ℬ)
            @test P isa TensOrtho{2}
            @test best_fit_ortho(A, ℬ) isa TensOrtho{2} && get_array(best_fit_ortho(A, ℬ)) ≈ get_array(A)
        end
        # Different frames: the generic route.
        A = TensOrtho{2}(3.0, 5.0, 7.0, RotatedBasis(0.7, -0.4, 1.1))
        B = TensOrtho{2}(2.0, -1.0, 0.5, RotatedBasis(0.1, 0.2, 0.3))
        @test !(A ⋅ B isa TensOrtho) && get_array(A ⋅ B) ≈ get_array(A) * get_array(B)
        @test !(A ⊗ B isa TensOrtho)
        @test_throws ArgumentError A + B
        C = _og_tens(_OG_N, _OG_S, RotatedBasis(0.1, 0.2, 0.3))
        @test !(C ⊡ A isa TensOrtho) && !(A ⊡ C isa TensOrtho)
        @test occursin("(5.0) P₂", sprint(show, A))
        @test_throws "a TensOrtho of order 2 stores 3 coefficients" TensOrtho{2, Float64}((1.0, 2.0), CanonicalBasis{3, Float64}())
        # Derivatives and symbols.
        ℬ = RotatedBasis(0.7, -0.4, 1.1)
        f(x) = get_data(inv(TensOrtho{2}(x, 2x, 3.0, ℬ) ⋅ TensOrtho{2}(1.0, x^2, 2.0, ℬ)))
        @test ForwardDiff.derivative(x -> f(x)[1], 1.5) ≈ -1 / 1.5^2
        @test ForwardDiff.derivative(x -> f(x)[2], 1.5) ≈ -3 / (2 * 1.5^4)
        @test (@inferred inv(TensOrtho{2}(ForwardDiff.Dual(1.0, 1.0), 2.0, 3.0, ℬ))) isa TensOrtho{2}
        x = SymPy.symbols("x", real = true)
        As = TensOrtho{2}(x, 2.0, 3.0, CanonicalBasis{3, Float64}())
        @test Float64(get_data(inv(As))[1].subs(x, 0.5)) ≈ 2.0
        @test Float64(get_data(As ⊡ _og_tens(_OG_N, _OG_S, CanonicalBasis{3, Float64}()))[1].subs(x, 0.5)) ≈ 0.5 * 160 + 2 * 31 + 3 * 22
    end

    @testsection "major transpose of the transversely isotropic shapes" begin
        n = normalize([0.36, 0.48, 0.8])
        for W in (
                TensTI{4}(31.0, 39.0, 11.0, 23.0, 21.0, n), TensTI{4}(31.0, 39.0, 11.0, 6.0, 23.0, 21.0, n),
                TensTI{4}(31.0, 39.0, 11.0, 6.0, 23.0, 21.0, 4.0, -3.0, n),
            )
            Wt = transpose(W)
            @test typeof(Wt) == typeof(W)
            @test get_array(Wt) ≈ permutedims(get_array(W), (3, 4, 1, 2)) atol = 1.0e-13
            @test transpose(Wt) == W
        end
        for w in (TensTI{2}(2.0, 5.0, n), TensTI{2}(2.0, 5.0, 1.5, n))
            @test typeof(transpose(w)) == typeof(w)
            @test get_array(transpose(w)) ≈ transpose(get_array(w))
        end
        Wz = TensTI{4}(31.0 + 1im, 39.0, 11.0, 6.0im, 23.0, 21.0, n)
        @test get_array(Wz') ≈ conj.(permutedims(get_array(Wz), (3, 4, 1, 2)))
    end
end
