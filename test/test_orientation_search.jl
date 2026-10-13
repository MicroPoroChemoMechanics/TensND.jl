# The projections with a free orientation (`src/orientation_search.jl`): the
# axis or frame recovered exactly for tensors of the class, never worse than a
# fixed orientation, global on random tensors against a dense search, the
# closed forms of order 2 against an exhaustive one, and determinism.

@testsection "Orientation search — projections with a free orientation" begin

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
        arr = Array(
            get_array(
                TensOrtho(
                    20.0, 8.0, 6.0, 30.0, 7.0, 40.0, 5.0, 6.0, 7.0,
                    CanonicalBasis{3, Float64}()
                )
            )
        )
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

    @testsection "proj_tens(:CUBIC, A) — order 4, rotated cube recovered" begin
        # The one class with no eigenstructure candidate to start from: every
        # second-order contraction of a cubic tensor is isotropic, so the search
        # begins from the screened orientations alone (`_cubic_starts`). It
        # still recovers a rotated cube exactly, which is the point.
        for fr in (
                CanonicalBasis{3, Float64}(),
                RotatedBasis(0.3, 0.4, 0.5),
                RotatedBasis(1.1, -0.7, 0.2),
            )
            C = tens_cubic(10.0, 4.0, 2.0, fr)
            arr = get_array(C)
            B, d, drel = proj_tens(:CUBIC, arr)
            @test B isa TensCubic
            @test drel < 1.0e-6
            # The coefficients are frame-independent under the octahedral group,
            # so they come back whichever of the 24 equivalent frames is found.
            @test all(isapprox.(get_data(B), get_data(C); rtol = 1.0e-5))
            @test is_CUBIC(arr)
        end

        # An isotropic tensor is cubic about every cube; a transversely
        # isotropic one about none.
        @test is_CUBIC(get_array(TensISO{3}(18.0, 6.0)))
        @test !is_CUBIC(get_array(tens_TI(10.0, 3.0, 2.5, 12.0, 2.0, [0.0, 0.0, 1.0])))

        # Optimizing can only improve on a fixed frame, the canonical one being
        # refined among the others.
        arr = get_array(tens_cubic(10.0, 4.0, 2.0, RotatedBasis(0.3, 0.4, 0.5)))
        @test proj_tens(:CUBIC, arr)[3] ≤
            proj_tens(:CUBIC, arr, CanonicalBasis{3, Float64}())[3] + 1.0e-12

        # And `best_sym_tens` finds it once it is allowed to search, where the
        # cheap path cannot.
        _, _, _, sym = best_sym_tens(Tens(arr); optimize_angles = true)
        @test sym == :CUBIC

        # Zero tensor: the degenerate branch returns zero without dividing by it.
        B, d, drel = proj_tens(:CUBIC, zeros(3, 3, 3, 3))
        @test all(iszero, get_data(B))
        @test d == 0 && drel == 0
    end

    @testsection "proj_tens — determinism and monotonicity" begin
        # No random search is involved, so a call is repeatable to the bit, and
        # the result is never worse than a fixed orientation the search refines.

        arr_ortho = Array(
            get_array(
                TensOrtho(
                    20.0, 8.0, 6.0, 30.0, 7.0, 40.0, 5.0, 6.0, 7.0,
                    Basis(0.3, 0.7, 0.2)
                )
            )
        )
        # A tensor with no symmetry at all, so the optimum is a genuine
        # compromise rather than an exactly recoverable frame.
        M = Float64[
            31.0 6.0 5.0 1.2 0.8 2.1
            6.0 27.0 4.5 0.7 1.9 1.1
            5.0 4.5 35.0 2.3 1.4 0.6
            1.2 0.7 2.3 12.0 1.1 0.9
            0.8 1.9 1.4 1.1 14.0 1.7
            2.1 1.1 0.6 0.9 1.7 16.0
        ]
        arr_aniso = Array(get_array(inv_KM(M)))

        for (name, arr) in (("exactly orthotropic", arr_ortho), ("fully anisotropic", arr_aniso))
            for sym in (:TI, :ORTHO)
                # 1. Repeatable: same input, bit-identical output.
                results = [proj_tens(sym, arr)[3] for _ in 1:4]
                @test all(==(results[1]), results)

                # 2. Never worse than the projection about the canonical
                #    orientation, which is always refined: a guarantee by
                #    construction, not a hope about the optimizer.
                fixed = sym === :TI ?
                    proj_tens(Val(:TI), arr, (0.0, 0.0, 1.0))[3] :
                    proj_tens(Val(:ORTHO), arr, CanonicalBasis{3, Float64}())[3]
                @test results[1] ≤ fixed + 1.0e-8
            end
        end

        # Order 2 obeys the same two properties.
        A2 = Float64[5.0 1.0 2.0; 1.0 8.0 3.0; 2.0 3.0 12.0]
        for sym in (:TI, :ORTHO)
            results = [proj_tens(sym, A2)[3] for _ in 1:4]
            @test all(==(results[1]), results)
        end
        # Diagonalizing a symmetric 3×3 *is* its orthotropic projection.
        @test proj_tens(:ORTHO, A2)[3] < 1.0e-10
    end

    @testsection "Order 2 — closed forms against an exhaustive search" begin
        # The best orthotropic frame is the eigenframe of the symmetric part,
        # and the best TI axis an eigenvector of an extreme eigenvalue: no
        # orientation of a dense set of axes does better.
        dense = TensND._hemisphere_spiral(4000)
        # A symmetric tensor, a non-symmetric one, and one with a double
        # eigenvalue.
        for A in (
                Float64[5.0 1.0 2.0; 1.0 8.0 3.0; 2.0 3.0 12.0],
                Float64[3.0 -2.0 0.5; 1.0 -1.0 4.0; 0.2 3.0 7.0],
                Float64[4.0 0.0 0.0; 0.0 4.0 0.0; 0.0 0.0 9.0],
            )
            drel = proj_tens(:TI, A)[3]
            best = minimum(x -> proj_tens(Val(:TI), A, TensND._n_from_angles(x[1], x[2]))[3], dense)
            @test drel ≤ best + 1.0e-12
            # The orthotropic projection keeps all of the symmetric part.
            Aₛ = (A + A') / 2
            @test proj_tens(:ORTHO, A)[3] ≈ norm(A - Aₛ) / norm(A) atol = 1.0e-12
        end
    end

    @testsection "Order 4 — the search is global on random tensors" begin
        # Against a dense search (3,000 axes; 600 axes × 12 turns for the
        # frames, each of the 12 best refined), on random anisotropic
        # stiffnesses: the minimum found is the dense one, to rounding. A
        # multi-start from a dozen points missed it for one TI tensor in twenty.
        rng = Random.Xoshiro(20261011)
        axes = TensND._hemisphere_spiral(3000)
        frames = [TensND.SVector(a[1], a[2], ψ) for a in TensND._hemisphere_spiral(600) for ψ in range(0, π / 2; length = 13)[1:12]]
        objective = Dict(:TI => TensND._obj_TI4, :ORTHO => TensND._obj_ORTHO4, :CUBIC => TensND._obj_CUBIC4)
        for _ in 1:8
            M = randn(rng, 6, 6)
            A = Array(frommandel(SymmetricTensor{4, 3}, M * M' + 6I))
            C = TensND._KM_of_array(A)
            s = sum(abs2, C)
            for S in (:TI, :ORTHO, :CUBIC)
                j(x) = objective[S](x, C, s)
                grid = S === :TI ? axes : frames
                vals = j.(grid)
                dense = minimum(j(TensND._refine_angles(j, x)) for x in grid[partialsortperm(vals, 1:12)])
                @test proj_tens(S, A)[3]^2 ≤ dense + 1.0e-12
            end
        end
    end

    @testsection "Newton refinement — a smooth minimum to machine precision" begin
        # A function of two angles with a minimum at (0.4, -1.1) and a saddle
        # elsewhere: the refinement reaches the minimum from nearby, and never
        # returns a point worse than its start.
        f(x) = (1 - cos(x[1] - 0.4)) + 2 * (1 - cos(x[2] + 1.1)) + 0.1 * sin(x[1]) * sin(x[2])^2
        x = TensND._refine_angles(f, TensND.SVector(1.0, -0.2))
        @test norm(ForwardDiff.gradient(f, x)) < 1.0e-10
        for x0 in (TensND.SVector(3.0, 2.0), TensND.SVector(0.0, 0.0), TensND.SVector(-2.5, 1.5))
            @test f(TensND._refine_angles(f, x0)) ≤ f(x0)
        end
    end

    @testsection "best_sym_tens — optimize_angles path" begin
        # `optimize_angles = true` searches the orientation.
        n_true = [0.0, 1 / √2, 1 / √2]
        t = Tens(get_array(TensTI{4}(3.0, 5.0, 1.5, 2.0, 2.5, Tuple(n_true))))

        B, d, drel, sym = best_sym_tens(t; optimize_angles = true)
        @test sym === :TI
        @test drel < 1.0e-6
        @test B ≈ t rtol = 1.0e-5
    end

end
