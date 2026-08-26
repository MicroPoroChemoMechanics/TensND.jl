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
    b = zeros(3, 3, 3, 3)
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        b[i, j, k, l] = (a[i, j, k, l] + a[j, i, k, l] + a[i, j, l, k] + a[j, i, l, k]) / 4
    end
    b
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
            d = get_data(isotropify(fromISO(t, nax)))
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

    @testset "ti8_params_from_KM — exact read-off in the e₃ frame" begin
        t = TensTI{4}(1.3, 2.1, 0.7, -0.4, 1.7, 0.9, 0.55, -0.22, (0.0, 0.0, 1.0))
        p = ti8_params_from_KM(mandel66_minor(get_array(t)))
        @test all(p .≈ TensND.get_ℓ8(t))
        # Round trip through the Kelvin-Mandel form.
        @test maximum(abs, KM_from_ti8_params(p) .- mandel66_minor(get_array(t))) < 1.0e-13
        # It is a projection, not just a read-off: a matrix that is NOT axially
        # invariant comes back as its azimuthal average.
        arr = _rand_minor4(rng)
        q = ti8_params_from_KM(mandel66_minor(arr))
        @test all(
            q .≈ ti8_params_from_KM(
                mandel66_minor(get_array(transverse_isotropify(Tens(arr), (0.0, 0.0, 1.0))))
            )
        )
    end

    @testset "mandel66_minor / array_from_mandel66 round trip" begin
        arr = _rand_minor4(rng)
        @test maximum(abs, array_from_mandel66(mandel66_minor(arr)) .- arr) < 1.0e-13
    end

    @testset "the Mandel block forms agree with the tensor ones" begin
        arr = _rand_minor4(rng)
        M = mandel66_minor(arr)
        α, β = iso_average_mandel66(M)
        d = get_data(isotropify(Tens(arr)))
        @test α ≈ d[1]
        @test β ≈ d[2]
        Mti = ti_average_mandel66(M, nax)
        @test maximum(
            abs, Mti .- mandel66_minor(get_array(transverse_isotropify(Tens(arr), nax)))
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

    @testset "best-fit projections are the projection component of proj_tens" begin
        arr = _rand_minor4(rng)
        t = Tens(arr)
        @test get_data(best_fit_iso(t)) == get_data(proj_tens(Val(:ISO), t)[1])
        @test get_data(best_fit_ti(t, nax)) == get_data(proj_tens(Val(:TI), t, nax)[1])
        # On a minor-symmetric tensor the isotropic fit IS the SO(3) average.
        @test all(get_data(best_fit_iso(t)) .≈ get_data(isotropify(t)))
        # The TI fit is NOT the azimuthal average: it forces major symmetry.
        ti_avg = transverse_isotropify(t, nax)
        ℓ = TensND.get_ℓ(ti_avg)
        @test !(ℓ[3] ≈ ℓ[4])                      # the average keeps ℓ₃ ≠ ℓ₄
        @test length(get_data(best_fit_ti(t, nax))) == 5   # the fit cannot
    end
end
