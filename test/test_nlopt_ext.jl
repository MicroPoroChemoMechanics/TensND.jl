# Loading NLopt activates TensNDNLoptExt.  This file MUST be included last:
# `test_tens_projection.jl` asserts the no-NLopt behaviour, and a package
# stays loaded for the rest of the session once `using` has run.
using NLopt

@testsection "NLopt extension — rotation-optimized projections" begin

    # The extension overrides the error-throwing fallbacks of
    # tens_projection.jl.  `test_tens_projection.jl` deliberately checks the
    # *absence* of NLopt; here NLopt is loaded, so the real optimizers run.

    @testsection "proj_tens(:TI, A) — order 4, exact TI recovered" begin
        # A tensor that IS transversely isotropic about a tilted axis must be
        # recovered exactly: the optimized projection has zero residual and
        # finds (up to sign) the right axis.
        for n_true in ([0.0, 0.0, 1.0], [1 / √2, 1 / √2, 0.0], [0.3, -0.4, √(1 - 0.25)])
            n_true = n_true ./ norm(n_true)
            A = TensTI{4}(3.0, 5.0, 1.5, 2.0, 2.5, Tuple(n_true))
            arr = get_array(A)

            B, d, drel = proj_tens(:TI, arr)
            @test B isa TensTI{4}
            @test drel < 1.0e-6
            @test d < 1.0e-6 * norm(arr)

            # The recovered axis spans the same line (sign is arbitrary).
            n_found = collect(axis(B))
            @test abs(dot(n_found, n_true)) ≈ 1.0 atol = 1.0e-5

            # And the rebuilt tensor matches the original.
            @test maximum(abs, get_array(B) - arr) < 1.0e-6 * maximum(abs, arr)
        end
    end

    @testsection "proj_tens(:TI, A) — order 4, isotropic input" begin
        # An isotropic tensor is TI about *every* axis: residual must vanish
        # whatever axis the optimizer settles on.
        A = TensISO{3}(30.0, 12.0)
        arr = Array(get_array(A))
        B, d, drel = proj_tens(:TI, arr)
        @test drel < 1.0e-6
        @test maximum(abs, get_array(B) - arr) < 1.0e-6 * maximum(abs, arr)
    end

    @testsection "proj_tens(:TI, A) — order 4, generic input" begin
        # A genuinely non-TI tensor: the optimizer must beat (or match) the
        # canonical-axis projection, never do worse.
        arr = Array(get_array(TensOrtho(
            20.0, 8.0, 6.0, 30.0, 7.0, 40.0, 5.0, 6.0, 7.0,
            CanonicalBasis{3, Float64}()
        )))
        B_opt, _, drel_opt = proj_tens(:TI, arr)
        @test B_opt isa TensTI{4}
        @test 0.0 ≤ drel_opt ≤ 1.0

        for n in ([0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0])
            _, _, drel_fixed = proj_tens(Val(:TI), arr, Tuple(n))
            @test drel_opt ≤ drel_fixed + 1.0e-8
        end
    end

    @testsection "proj_tens(:TI, A) — order 2" begin
        n_true = [1 / √2, 1 / √2, 0.0]
        A = TensTI{2}(5.0, 8.0, Tuple(n_true))
        arr = Array(get_array(A))

        B, d, drel = proj_tens(:TI, arr)
        @test B isa TensTI{2}
        @test drel < 1.0e-6
        @test abs(dot(collect(axis(B)), n_true)) ≈ 1.0 atol = 1.0e-5
        @test maximum(abs, get_array(B) - arr) < 1.0e-6 * maximum(abs, arr)
    end

    @testsection "proj_tens(:ORTHO, A) — order 4, exact ORTHO recovered" begin
        # Orthotropic in the canonical frame, then rotated: the optimizer must
        # find a frame giving a vanishing residual.
        C = TensOrtho(
            20.0, 8.0, 6.0, 30.0, 7.0, 40.0, 5.0, 6.0, 7.0,
            CanonicalBasis{3, Float64}()
        )
        arr = Array(get_array(C))

        B, d, drel = proj_tens(:ORTHO, arr)
        @test B isa TensOrtho
        @test drel < 1.0e-6
        @test maximum(abs, get_array(B) - arr) < 1.0e-6 * maximum(abs, arr)
    end

    @testsection "proj_tens(:ORTHO, A) — order 4, rotated frame" begin
        ℬ = Basis(0.3, 0.7, 0.2)
        C = TensOrtho(20.0, 8.0, 6.0, 30.0, 7.0, 40.0, 5.0, 6.0, 7.0, ℬ)
        arr = Array(get_array(C))

        B, d, drel = proj_tens(:ORTHO, arr)
        @test B isa TensOrtho
        @test drel < 1.0e-5
        @test maximum(abs, get_array(B) - arr) < 1.0e-5 * maximum(abs, arr)
    end

    @testsection "proj_tens(:ORTHO, A) — order 2" begin
        A = [4.0 0.0 0.0; 0.0 7.0 0.0; 0.0 0.0 11.0]
        B, d, drel = proj_tens(:ORTHO, A)
        # A diagonal 2nd-order tensor is orthotropic in the canonical frame:
        # the projection is exact.
        @test drel < 1.0e-6
        @test maximum(abs, B - A) < 1.0e-6 * maximum(abs, A)
    end

    @testsection "proj_tens — zero tensor degenerate branches" begin
        # `sqnorm ≈ 0` short-circuits before the optimizer in all four methods.
        Z4 = zeros(3, 3, 3, 3)
        B, d, drel = proj_tens(:TI, Z4)
        @test B isa TensTI{4}
        @test d == 0.0 && drel == 0.0
        @test all(iszero, get_array(B))

        B, d, drel = proj_tens(:ORTHO, Z4)
        @test d == 0.0 && drel == 0.0
        @test all(iszero, get_array(B))

        Z2 = zeros(3, 3)
        B, d, drel = proj_tens(:TI, Z2)
        @test B isa TensTI{2}
        @test d == 0.0 && drel == 0.0

        B, d, drel = proj_tens(:ORTHO, Z2)
        @test d == 0.0 && drel == 0.0
        @test all(iszero, B)
    end

    @testsection "best_sym_tens — optimize_angles path" begin
        # With NLopt loaded, `optimize_angles = true` goes through the
        # extension instead of throwing.
        n_true = [0.0, 1 / √2, 1 / √2]
        t = Tens(get_array(TensTI{4}(3.0, 5.0, 1.5, 2.0, 2.5, Tuple(n_true))))

        B, d, drel, sym = best_sym_tens(t; optimize_angles = true)
        @test sym === :TI
        @test drel < 1.0e-6
        @test maximum(abs, get_array(B) - get_array(t)) < 1.0e-5 * maximum(abs, get_array(t))
    end

end
