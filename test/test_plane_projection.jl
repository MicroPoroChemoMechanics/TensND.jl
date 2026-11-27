# =============================================================================
#  test_plane_projection.jl — projections of a tensor of the plane onto the
#  isotropic, square and orthotropic classes, with the frame given or free, the
#  symmetry search and the average over the rotations of the plane.
#
#  References, independent of the code under test:
#  - a projection is orthogonal: its residual is orthogonal to every tensor of
#    the class in the frame, whose basis is written here from Kelvin-Mandel
#    unit matrices rotated by explicit index sums;
#  - with the frame free, the distance is compared with the least distance over
#    a fine grid of frames, and a tensor of a class is recovered whatever the
#    frame it was built on;
#  - the average over the rotations is compared with a quadrature over angles;
#  - the transformation of the harmonics the closed forms rest on is checked
#    with symbols.
# =============================================================================

const PLP_SEED = 20261128
const PLP_PAIRS = ((1, 1), (2, 2), (1, 2))

plp_rot(θ) = [cos(θ) -sin(θ); sin(θ) cos(θ)]

"The array of the plane tensor of Kelvin-Mandel matrix `M` in the frame of rotation `R`."
function plp_array(M, R)
    w(I) = I == 3 ? sqrt(2) : 1.0
    A = zeros(2, 2, 2, 2)
    for (I, (i, j)) in enumerate(PLP_PAIRS), (J, (k, l)) in enumerate(PLP_PAIRS)
        for (a, b) in ((i, j), (j, i)), (c, d) in ((k, l), (l, k))
            A[a, b, c, d] = M[I, J] / (w(I) * w(J))
        end
    end
    return ein"ia,jb,kc,ld,abcd->ijkl"(R, R, R, R, A)
end

"A basis of the major-symmetric orthotropic tensors of the frame, then the two more of the general ones (as Kelvin-Mandel unit matrices)."
plp_unit(I, J) = (M = zeros(3, 3); M[I, J] = 1.0; M[J, I] = 1.0; M)
const PLP_ORTHO = (plp_unit(1, 1), plp_unit(2, 2), plp_unit(1, 2), plp_unit(3, 3))
const PLP_SQUARE = (plp_unit(1, 1) + plp_unit(2, 2), plp_unit(1, 2), plp_unit(3, 3))

"A random array of order 4 with the minor symmetries only."
function plp_minor(rng)
    A = randn(rng, 2, 2, 2, 2)
    A = (A + permutedims(A, (2, 1, 3, 4))) / 2
    return (A + permutedims(A, (1, 2, 4, 3))) / 2
end
plp_major(rng) = (M = randn(rng, 3, 3); plp_array(M + M', plp_rot(0.0)))

plp_inferred_free(A, v) = @inferred project(A, v)
plp_inferred_fixed(A, v, ℬ) = @inferred project(A, v, ℬ)

@testsection "Projections in the plane" begin
    rng = Xoshiro(PLP_SEED)

    @testsection "frame given: the residual is orthogonal to the class" begin
        for _ in 1:4
            θ = 2π * rand(rng)
            ℬ, R = Basis(θ), plp_rot(θ)
            for A in (plp_minor(rng), plp_major(rng), randn(rng, 2, 2, 2, 2))
                for (v, class, T) in ((Val(:ORTHO), PLP_ORTHO, TensOrtho{4, 2, Float64, 4}), (Val(:CUBIC), PLP_SQUARE, TensCubic{2, Float64}))
                    B, d, drel = project(A, v, ℬ)
                    @test B isa T
                    @test frame(B) == ℬ
                    res = A - get_array(B)
                    @test d ≈ norm(res) && drel ≈ d / norm(A)
                    @test all(E -> abs(sum(res .* plp_array(E, R))) < 1.0e-12 * norm(A), class)
                    # A projection: applied twice, it changes nothing.
                    @test project(get_array(B), v, ℬ)[2] < 1.0e-12 * norm(A)
                end
            end
            # Order 2: the diagonal in the frame, and the isotropic part for the square class.
            a = randn(rng, 2, 2)
            Ba, da, _ = project(a, Val(:ORTHO), ℬ)
            @test Ba isa TensOrtho{2, 2, Float64, 2}
            @test collect(get_data(Ba)) ≈ diag(R' * a * R)
            @test da ≈ norm(a - R * Diagonal(diag(R' * a * R)) * R')
            @test project(a, Val(:CUBIC), ℬ)[1] ≈ TensISO{2}(tr(a) / 2)
        end
        # Symbols by their name, and through a tensor.
        A = plp_major(rng)
        @test project(A, :ORTHO, Basis(0.3))[1] == project(A, Val(:ORTHO), Basis(0.3))[1]
        @test best_fit_ortho(Tens(A), Basis(0.3)) == project(A, Val(:ORTHO), Basis(0.3))[1]
        @test best_fit_cubic(Tens(A), Basis(0.3)) == project(A, Val(:CUBIC), Basis(0.3))[1]
        @test best_fit_ortho(Tens(randn(rng, 2, 2)), Basis(0.3)) isa TensOrtho{2, 2}
        @test best_fit_cubic(Tens(randn(rng, 2, 2)), Basis(0.3)) isa TensISO{2, 2}
        # A tensor written in a rotated basis: the same tensor, the same projection.
        t = express(Tens(A), Basis(1.1))
        @test mandel(project(t, :ORTHO, Basis(0.3))[1]) ≈ mandel(project(A, :ORTHO, Basis(0.3))[1])
        @test mandel(project(t, :CUBIC, Basis(0.3))[1]) ≈ mandel(project(A, :CUBIC, Basis(0.3))[1])
    end

    @testsection "frame free: the least distance over the frames" begin
        grid = range(0, π / 2; length = 2001)
        for _ in 1:6
            A = plp_major(rng)
            for v in (Val(:ORTHO), Val(:CUBIC))
                B, d, drel = project(A, v)
                @test d ≈ norm(A - get_array(B))
                dgrid = minimum(θ -> project(A, v, Basis(θ))[2], grid)
                @test d ≤ dgrid + 1.0e-12
                @test d ≥ dgrid - 1.0e-4 * norm(A)        # the grid is within its step of the least
            end
            # Order 2: all of the symmetric part is kept.
            a = randn(rng, 2, 2)
            Ba, da, _ = project(a, Val(:ORTHO))
            @test da ≈ norm(a - a') / 2
            @test get_array(Ba) ≈ (a + a') / 2
        end
        # A tensor of a class is recovered, whatever its frame.
        for _ in 1:5
            θ = 2π * rand(rng)
            O = TensOrtho(10.0 + rand(rng), 6.0 + rand(rng), 3.0 * rand(rng), 2.0 + rand(rng), Basis(θ))
            Q = cubic_from_constants(10.0 + rand(rng), 4.0 * rand(rng), 2.0 + rand(rng), Basis(θ))
            L = iso_stiffness(Val(2), 3.0 + rand(rng), 1.0 + rand(rng))
            for (t, v) in ((O, :ORTHO), (Q, :CUBIC), (Q, :ORTHO), (L, :ORTHO), (L, :CUBIC))
                B, d, _ = project(Tens(get_array(t)), v)
                @test d < 1.0e-12 * norm(get_array(t))
                @test mandel(B) ≈ mandel(t)
            end
            # The frame found is that of the tensor, up to a quarter turn (an eighth for the square).
            E = vecbasis(frame(project(Tens(get_array(O)), :ORTHO)[1]), :cov)
            @test abs(abs(E[1, 1] * cos(θ) + E[2, 1] * sin(θ)) - 0.5) ≈ 0.5 atol = 1.0e-9
        end
        # The harmonic of order 4 alone, or of order 2 alone: the quartic is singular there.
        only4 = plp_array([1.0 -1 0; -1 1 0; 0 0 0.5], plp_rot(0.4))            # square, no harmonic of order 2
        @test project(only4, Val(:ORTHO))[2] < 1.0e-12
        only2 = plp_array([3.0 0 0; 0 1 0; 0 0 2.0], plp_rot(0.9))               # b = (3 + 1 − 4)/4 = 0
        @test project(only2, Val(:ORTHO))[2] < 1.0e-12
        # A null tensor is of every class.
        for v in (Val(:ORTHO), Val(:CUBIC))
            B, d, drel = project(zeros(2, 2, 2, 2), v)
            @test iszero(d) && iszero(drel) && all(iszero, get_data(B))
            @test all(iszero, get_data(project(zeros(2, 2, 2, 2), v, Basis(0.3))[1]))
        end
        @test all(iszero, get_data(project(zeros(2, 2), Val(:ORTHO))[1]))
        @test all(iszero, get_data(project(zeros(2, 2), Val(:ORTHO), Basis(0.3))[1]))
        # Space is not affected.
        @test project(get_array(rand(rng, Tens{4, 3})), :ORTHO)[1] isa TensOrtho{4, 3}
    end

    @testsection "symmetry search and predicates" begin
        θ = 0.7
        O = Tens(get_array(TensOrtho(10.0, 6.0, 3.0, 2.0, Basis(θ))))
        Q = Tens(get_array(cubic_from_constants(10.0, 4.0, 2.0, Basis(θ))))
        L = Tens(get_array(iso_stiffness(Val(2), 3.0, 1.0)))
        G = Tens(plp_major(rng))
        @test last.((best_sym_tens(L), best_sym_tens(Q), best_sym_tens(O), best_sym_tens(G))) == (:ISO, :CUBIC, :ORTHO, :ANISO)
        @test best_sym_tens(Q)[1] isa TensCubic{2} && best_sym_tens(O)[1] isa TensOrtho{4, 2}
        @test best_sym_tens(O; optimize_angles = true)[4] === :ORTHO
        @test best_sym_tens(O; proj = (:ISO, :CUBIC))[4] === :ANISO
        @test best_sym_tens(O, Basis(θ))[4] === :ORTHO
        @test best_sym_tens(O, Basis(θ + 0.3))[4] === :ANISO
        @test best_sym_tens(express(O, Basis(1.9)))[4] === :ORTHO
        @test mandel(best_sym_tens(express(O, Basis(1.9)))[1]) ≈ mandel(O)
        a = Tens(get_array(TensOrtho{2}(2.0, 5.0, Basis(θ))))
        @test best_sym_tens(a)[4] === :ORTHO
        @test best_sym_tens(a, Basis(θ))[4] === :ORTHO
        @test best_sym_tens(Tens([2.0 0; 0 2.0]))[4] === :ISO
        @test best_sym_tens(Tens([2.0 1; -1 2.0]))[4] === :ANISO
        @test (isorthotropic(O), isorthotropic(Q), isorthotropic(G)) == (true, true, false)
        @test (iscubic(O), iscubic(Q), isisotropic(Q), isisotropic(L)) == (false, true, false, true)
        @test (isorthotropic(O, Basis(θ)), isorthotropic(O, Basis(θ + 0.3)), iscubic(Q, Basis(θ)), iscubic(O, Basis(θ))) == (true, false, true, false)
        @test isorthotropic(get_array(O)) && isorthotropic(get_array(O), Basis(θ))
        # The plane has no transversely isotropic class.
        @test_throws ArgumentError istransverselyisotropic(O)
        @test_throws ArgumentError istransverselyisotropic(O, (1.0, 0.0))
        @test_throws ArgumentError project(O, :TI, (1.0, 0.0))
        @test_throws ArgumentError project(O, Val(:TI))
        @test_throws ArgumentError project(get_array(O), Val(:TI))
        @test_throws ArgumentError project(get_array(O), Val(:TI), (1.0, 0.0))
        @test_throws ArgumentError project(get_array(a), Val(:TI))
        @test_throws ArgumentError project(get_array(a), Val(:TI), (1.0, 0.0))
        @test_throws ArgumentError project(a, :TI)
        @test_throws ArgumentError best_sym_tens(O; proj = (:ISO, :TI))
        @test_throws ArgumentError best_sym_tens(O, Basis(θ); proj = (:TI,))
        @test_throws "no transversely isotropic class" istransverselyisotropic(get_array(O))
    end

    @testsection "average over the rotations of the plane" begin
        n = 256
        θs = (0:(n - 1)) .* (2π / n)
        for _ in 1:3
            # Order 4, minor symmetries only: the part on 𝕎 is kept.
            G = plp_minor(rng)
            Mavg = sum(θ -> mandel_rotation(θ)' * mandel(Tens(G); symmetrize = true) * mandel_rotation(θ), θs) / n
            H = transverse_isotropify(Tens(G))
            @test H isa AbstractTens{4, 2}
            @test mandel(H; symmetrize = true) ≈ Mavg
            @test norm(Mavg - Mavg') > 1.0e-3                                       # it has no major symmetry
            @test mandel(isotropify(Tens(G))) ≈ (Mavg + Mavg') / 2                  # the reflections remove it
            @test mandel(transverse_isotropify(H); symmetrize = true) ≈ Mavg        # an average is a projection
            # With the major symmetry the two averages coincide.
            S = plp_major(rng)
            @test mandel(transverse_isotropify(Tens(S)); symmetrize = true) ≈ mandel(isotropify(Tens(S)))
            # Order 2: the trace and the skew part.
            g = randn(rng, 2, 2)
            gavg = sum(θ -> plp_rot(θ)' * g * plp_rot(θ), θs) / n
            @test get_array(transverse_isotropify(Tens(g))) ≈ gavg
            @test gavg ≈ [tr(g) / 2 (g[1, 2] - g[2, 1]) / 2; (g[2, 1] - g[1, 2]) / 2 tr(g) / 2]
            # In a rotated basis the average is the same tensor.
            @test canonical_components(transverse_isotropify(express(Tens(g), Basis(0.8)))) ≈ gavg
        end
    end

    @testsection "inference, derivatives, symbols" begin
        A = SymmetricTensor{4, 2}(plp_major(rng))
        ℬ = Basis(0.6)
        @test plp_inferred_free(A, Val(:CUBIC)) isa Tuple{TensCubic{2, Float64, RotatedBasis{2, Float64}}, Float64, Float64}
        @test plp_inferred_free(A, Val(:ORTHO)) isa Tuple{TensOrtho{4, 2, Float64, 4, RotatedBasis{2, Float64}}, Float64, Float64}
        @test plp_inferred_fixed(Array(A), Val(:ORTHO), ℬ)[1] isa TensOrtho{4, 2, Float64, 4}
        @test plp_inferred_fixed(Array(A), Val(:CUBIC), ℬ)[1] isa TensCubic{2, Float64}
        # The projection is linear in the tensor: its derivative along a direction is the projection of the direction.
        A₀, dA = plp_major(rng), plp_major(rng)
        for v in (Val(:ORTHO), Val(:CUBIC))
            D = ForwardDiff.derivative(x -> collect(get_data(project(A₀ + x * dA, v, ℬ)[1])), 0.0)
            @test D ≈ collect(get_data(project(dA, v, ℬ)[1]))
        end
        # With respect to the angle of the frame: C₁₁′(θ) of the orthotropic projection is the
        # entry (1,1) of the Kelvin-Mandel matrix in the frame, first and second derivatives.
        M₀ = mandel(Tens(A₀))
        c11(θ) = get_data(project(A₀, Val(:ORTHO), Basis(θ))[1])[1]
        c11_ref(θ) = (mandel_rotation(θ)' * M₀ * mandel_rotation(θ))[1, 1]
        @test ForwardDiff.derivative(c11, 0.4) ≈ ForwardDiff.derivative(c11_ref, 0.4)
        @test ForwardDiff.derivative(θ -> ForwardDiff.derivative(c11, θ), 0.4) ≈
            ForwardDiff.derivative(θ -> ForwardDiff.derivative(c11_ref, θ), 0.4)

        # The harmonics under a rotation of the frame by θ (header of the source file).
        θ = symbols("θ", real = true)
        C1111, C2222, C1122, C1212, C1112, C2212 = symbols("C1111 C2222 C1122 C1212 C1112 C2212", real = true)
        r2 = sqrt(Sym(2))
        M = [C1111 C1122 r2 * C1112; C1122 C2222 r2 * C2212; r2 * C1112 r2 * C2212 2C1212]
        d₁₂, δ, a, b = TensND._plane_harmonics(M)
        @test simplify.([d₁₂ - (C1112 + C2212), δ - (C1111 - C2222) / 2, a - (C1112 - C2212), b - (C1111 + C2222 - 2C1122 - 4C1212) / 4]) == zeros(4)
        Q = mandel_rotation(θ)
        d₁₂′, δ′, a′, b′ = TensND._plane_harmonics(transpose(Q) * M * Q)
        @test simplify(d₁₂′ - (d₁₂ * cos(2θ) - δ * sin(2θ))) == 0
        @test simplify(δ′ - (δ * cos(2θ) + d₁₂ * sin(2θ))) == 0
        @test simplify(a′ - (a * cos(4θ) - b * sin(4θ))) == 0
        @test simplify(b′ - (b * cos(4θ) + a * sin(4θ))) == 0
        # A projection with symbols, on the canonical frame.
        c = symbols("c11 c22 c12 c66", positive = true)
        ℬˢ = CanonicalBasis{2, Sym}()
        Aˢ = get_array(TensOrtho(c..., ℬˢ))
        α, β, γ = get_data(project(Aˢ, Val(:CUBIC), ℬˢ)[1])
        @test simplify.([α - ((c[1] + c[2]) / 2 + c[3]), β - ((c[1] + c[2]) / 2 - c[3]), γ - 2c[4]]) == zeros(3)
        @test simplify.(collect(get_data(project(Aˢ, Val(:ORTHO), ℬˢ)[1])) .- collect(c)) == zeros(4)
    end
end
