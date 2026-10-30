# =============================================================================
#  test_tens_rotational_average.jl — exact SO(3) and azimuthal averages.
#
#  Oracles:
#  1. the generic array method, which builds the 3×3×3×3 components and
#     contracts them, is the reference for every structured fast path;
#  2. a discrete azimuthal quadrature: the average of a 4th-order tensor over
#     N ≥ 5 equally spaced angles about the axis equals the continuous average
#     exactly, the integrand being a trigonometric polynomial of degree ≤ 4;
#  3. the two tensors whose isotropic average is known by inspection, 𝕀 and 𝕁.
# =============================================================================

_rand_minor4(rng) = begin
    a = randn(rng, 3, 3, 3, 3)
    (a + permutedims(a, (2, 1, 3, 4)) + permutedims(a, (1, 2, 4, 3)) + permutedims(a, (2, 1, 4, 3))) / 4
end

_rot_axis(n, φ) = begin
    nv = collect(n) ./ norm(collect(n))
    K = [0 -nv[3] nv[2]; nv[3] 0 -nv[1]; -nv[2] nv[1] 0]
    # `Matrix(1.0I, 3, 3) + …`, NOT `I(3) .+ …`: on Julia 1.14-DEV a broadcast
    # against `I(3)::Diagonal{Bool}` preserves the `Diagonal` structure and
    # silently DROPS every off-diagonal entry, so the "rotation" came back
    # diagonal (det 0.596) and the quadrature oracle was wrong while the code
    # under test was right. Ordinary `+` on a dense matrix is version-proof.
    Matrix(1.0I, 3, 3) + sin(φ) * K + (1 - cos(φ)) * (K * K)
end

_rotate4_ref(arr, R) = begin
    out = zeros(3, 3, 3, 3)
    for a in 1:3, b in 1:3, c in 1:3, d in 1:3
        s = 0.0
        for i in 1:3, j in 1:3, k in 1:3, l in 1:3
            s += R[i, a] * R[j, b] * R[k, c] * R[l, d] * arr[i, j, k, l]
        end
        out[a, b, c, d] = s
    end
    out
end

@testsection "rotational averages" begin

    rng = MersenneTwister(4711)
    nax = (0.36, -0.48, 0.8)

    @testset "isotropify — structured fast paths == generic array path" begin
        # The generic method is reached through `Tens`, whose array carries no
        # structure to dispatch on.
        cases = (
            TensTI{4}(1.3, 2.1, 0.7, 1.7, 0.4, nax),                        # N = 5
            TensTI{4}(1.3, 2.1, 0.7, -0.4, 1.7, 0.9, nax),                  # N = 6
            # N = 8: ℓ₇ and ℓ₈ must contribute to neither invariant.
            TensTI{4}(1.3, 2.1, 0.7, -0.4, 1.7, 0.9, 0.55, -0.22, nax),
        )
        for t in cases
            fast = get_data(isotropify(t))
            gen = get_data(isotropify(Tens(get_array(t))))
            @test fast[1] ≈ gen[1] atol = 1.0e-12
            @test fast[2] ≈ gen[2] atol = 1.0e-12
        end

        t2 = TensTI{2}(1.4, -0.6, nax)
        @test get_data(isotropify(t2))[1] ≈ get_data(isotropify(Tens(get_array(t2))))[1]
        t2c = TensTI{2}(1.4, -0.6, 0.9, nax)      # antisymmetric part: traceless
        @test get_data(isotropify(t2c))[1] ≈ get_data(isotropify(t2))[1]

        # An isotropic tensor is its own average, at either order.
        @test get_data(isotropify(TensISO{3}(2.0, 5.0))) == (2.0, 5.0)
        @test get_data(isotropify(TensISO{3}(3.0))) == (3.0,)
    end

    @testset "isotropify — the two averages known by inspection" begin
        for (t, want) in (
                (tens_Id4(Val(3), Val(Float64)), (1.0, 1.0)),
                (tens_J4(Val(3), Val(Float64)), (1.0, 0.0)),
            )
            d = get_data(isotropify(iso_to_ti(t, nax)))
            @test d[1] ≈ want[1] atol = 1.0e-14
            @test d[2] ≈ want[2] atol = 1.0e-14
        end
    end

    @testset "transverse_isotropify vs discrete quadrature" begin
        arr = _rand_minor4(rng)
        N = 9
        avg = zeros(3, 3, 3, 3)
        for m in 0:(N - 1)
            avg .+= _rotate4_ref(arr, permutedims(_rot_axis(nax, 2π * m / N)))
        end
        avg ./= N
        got = transverse_isotropify(Tens(arr), nax)
        @test maximum(abs, get_array(got) .- avg) < 1.0e-11
    end

    @testset "ti8_params_from_mandel — exact read-off in the e₃ frame" begin
        t = TensTI{4}(1.3, 2.1, 0.7, -0.4, 1.7, 0.9, 0.55, -0.22, (0.0, 0.0, 1.0))
        p = ti8_params_from_mandel(mandel(get_array(t); symmetrize = true))
        @test all(p .≈ TensND.walpole_coeffs(t, Val(8)))
        # Round trip through the Kelvin-Mandel form.
        @test maximum(abs, mandel_from_ti8_params(p) .- mandel(get_array(t); symmetrize = true)) < 1.0e-13
        # It is a projection, not just a read-off: a matrix that is NOT axially
        # invariant comes back as its azimuthal average.
        arr = _rand_minor4(rng)
        q = ti8_params_from_mandel(mandel(arr; symmetrize = true))
        @test all(
            q .≈ ti8_params_from_mandel(
                mandel(get_array(transverse_isotropify(Tens(arr), (0.0, 0.0, 1.0))); symmetrize = true)
            )
        )
    end

    @testset "mandel(·; symmetrize = true) / from_mandel round trip" begin
        arr = _rand_minor4(rng)
        @test maximum(abs, get_array(from_mandel(mandel(arr; symmetrize = true))) .- arr) < 1.0e-13
    end

    @testset "the Mandel block forms agree with the tensor ones" begin
        arr = _rand_minor4(rng)
        M = mandel(arr; symmetrize = true)
        α, β = iso_average_mandel(M)
        d = get_data(isotropify(Tens(arr)))
        @test α ≈ d[1]
        @test β ≈ d[2]
        Mti = ti_average_mandel(M, nax)
        @test maximum(
            abs, Mti .- mandel(get_array(transverse_isotropify(Tens(arr), nax)); symmetrize = true)
        ) < 1.0e-11
    end

    @testset "a symbolic axis reaches the azimuthal average" begin
        # `_axis_frame` used to pick its in-plane reference with `argmin(abs)`,
        # which a symbolic axis cannot answer.
        θ = symbols("theta", real = true)
        n = (sin(θ), Sym(0), cos(θ))
        t = TensISO{3}(Sym(3), Sym(2))
        got = transverse_isotropify(t, n)
        @test got isa TensTI{4}
        @test all(iszero, tsimplify.(collect(axis(got)) .- collect(n)))
    end

    @testset "ForwardDiff through the averages" begin
        arr = _rand_minor4(rng)
        f = x -> get_data(transverse_isotropify(Tens(arr .* x), nax))[7]
        g = ForwardDiff.derivative(f, 2.0)
        h = 1.0e-6
        @test g ≈ (f(2.0 + h) - f(2.0 - h)) / 2h atol = 1.0e-6
    end

    @testset "best-fit projections are the projection component of project" begin
        arr = _rand_minor4(rng)
        t = Tens(arr)
        @test get_data(best_fit_iso(t)) == get_data(project(t, Val(:ISO))[1])
        @test get_data(best_fit_ti(t, nax)) == get_data(project(t, Val(:TI), nax)[1])
        # On a minor-symmetric tensor the isotropic fit IS the SO(3) average.
        @test all(get_data(best_fit_iso(t)) .≈ get_data(isotropify(t)))
        # The TI fit is NOT the azimuthal average: it forces major symmetry.
        ti_avg = transverse_isotropify(t, nax)
        ℓ = TensND.walpole_coeffs(ti_avg)
        @test !(ℓ[3] ≈ ℓ[4])                      # the average keeps ℓ₃ ≠ ℓ₄
        @test length(get_data(best_fit_ti(t, nax))) == 5   # the fit cannot
    end

    @testset "an integer axis" begin
        # The normalized axis was converted back to the type of `n`, which
        # threw `InexactError` for an integer axis such as (1, 1, 1).
        t = TensISO{3}(3.0, 2.0)
        @test get_data(transverse_isotropify(t, (1, 1, 1))) ==
            get_data(transverse_isotropify(t, (1.0, 1.0, 1.0)))
        M = Matrix(mandel(t))
        @test TensND.ti_average_mandel(M, (1, 1, 1)) == TensND.ti_average_mandel(M, (1.0, 1.0, 1.0))
    end
end

@testsection "Mandel matrices of the minor-symmetric part, in 3-D and 2-D, and the TI average of a 3×3" begin
    rng = Xoshiro(20261027)
    w = (1, 1, 1, sqrt(2), sqrt(2), sqrt(2))
    pairs = ((1, 1), (2, 2), (3, 3), (2, 3), (1, 3), (1, 2))
    # On components that have the minor symmetries, the entries are read exactly.
    S = Array(SymmetricTensor{4, 3}((i, j, k, l) -> randn(rng)))
    @test mandel(S; symmetrize = true) == [w[I] * w[J] * S[pairs[I]..., pairs[J]...] for I in 1:6, J in 1:6]
    # On any array, the minor-symmetric part, against Tensors.jl's.
    A = randn(rng, 3, 3, 3, 3)
    # Inferred where the keyword is written at the call, as it is in practice.
    symmetrized(A) = mandel(A; symmetrize = true)
    M = @inferred symmetrized(A)
    @test M ≈ tomandel(minorsymmetric(Tensor{4, 3}(A)))
    @test mandel(get_array(from_mandel(M)); symmetrize = true) ≈ M
    # A tensor: in its own basis, or in the one given, as `mandel`.
    C = Tens(SymmetricTensor{4, 3}((i, j, k, l) -> randn(rng)), RotatedBasis(0.3, 0.2, 0.1))
    @test mandel(C; symmetrize = true) ≈ mandel(C)
    can = CanonicalBasis{3, Float64}()
    @test mandel(C, can; symmetrize = true) ≈ mandel(C, can)
    @test mandel(C, can, (:cont, :cont, :cont, :cont); symmetrize = true) ≈ mandel(C, can)
    T = TensTI{4}(1.3, 2.1, 0.7, 1.7, 0.9, (0.0, 0.6, 0.8))
    @test mandel(T; symmetrize = true) ≈ mandel(T)
    @test mandel(Tensor{4, 3}(A); symmetrize = true) ≈ M
    # Order 2: the symmetric part.
    G = randn(rng, 3, 3)
    @test mandel(G; symmetrize = true) ≈ tomandel(symmetric(Tensor{2, 3}(G)))
    G2 = randn(rng, 2, 2)
    @test mandel(G2; symmetrize = true) ≈ tomandel(symmetric(Tensor{2, 2}(G2)))
    # 2-D: the Kelvin-Mandel matrix of an isotropic 2-D tensor, and the round trip.
    C2 = TensISO{2}(3.0, 2.0)
    @test mandel(Array(get_array(C2)); symmetrize = true) ≈ mandel(C2)
    B = randn(rng, 2, 2, 2, 2)
    @test mandel(B; symmetrize = true) ≈ tomandel(minorsymmetric(Tensor{4, 2}(B)))
    @test get_array(from_mandel(mandel(C2))) ≈ Array(get_array(C2))
    # No Kelvin-Mandel form in another dimension.
    @test_throws DimensionMismatch mandel(randn(rng, 4, 4, 4, 4); symmetrize = true)
    @test_throws DimensionMismatch mandel(randn(rng, 4, 4); symmetrize = true)
    # TI average of a 3×3 block, against the tensor-level average.
    n = random_axis(rng)
    @test ti_average_matrix(G, n) ≈ Array(get_array(transverse_isotropify(Tens(G), n)))
    @test ti_average_matrix(ti_average_matrix(G, n), n) ≈ ti_average_matrix(G, n)
    # Duals and symbols.
    a = (A[2, 3, 2, 3] + A[3, 2, 2, 3] + A[2, 3, 3, 2] + A[3, 2, 3, 2]) / 2
    @test ForwardDiff.derivative(x -> mandel(A .* x; symmetrize = true)[4, 4], 2.0) ≈ a
    @syms s::real
    @test mandel(Array(get_array(TensISO{2}(s, 2s))); symmetrize = true)[3, 3] == 2s
end

@testsection "The average of a TI tensor about another axis" begin
    # Not invariant about the new axis, so averaged in full: the same as the
    # average of its generic form.
    rng = Xoshiro(20261028)
    n1, n2 = random_axis(rng), random_axis(rng)
    A = TensTI{4}(10.0, 7.0, 3.0, 4.0, 5.0, n1)
    @test transverse_isotropify(A, n2) ≈ transverse_isotropify(generic_tens(A), n2)
    @test !(transverse_isotropify(A, n2) ≈ generic_tens(A))
end
