# =============================================================================
#  test_repatriated.jl — small tensor utilities that client packages used to
#  write for themselves: typed constructors, the exact membership of a
#  symmetry class, and the constants of a transversely isotropic tensor
#  without the major symmetry.
#
#  References: the components read on the full array, explicit index sums, and
#  the major-symmetric forms already tested elsewhere.
# =============================================================================

using StaticArrays: SMatrix

rep_inferred_static(a, ℬ) = @inferred static_tens(a, ℬ)
rep_inferred_symmetric(a, ℬ) = @inferred symmetric_tens(a, ℬ)
rep_inferred_hill(t) = @inferred hill_moduli(t, Val(6))

@testsection "Tensor utilities shared with client packages" begin
    rng = Xoshiro(20261129)

    @testsection "typed constructors" begin
        for d in (2, 3)
            ℬ = random_basis(rng, Val(d))
            a = randn(rng, d, d)
            A = randn(rng, d, d, d, d)
            As = (A + permutedims(A, (2, 1, 3, 4))) / 2
            As = (As + permutedims(As, (1, 2, 4, 3))) / 2
            # The components of `Tens(a, ℬ)`, with a type known from the arguments.
            @test rep_inferred_symmetric((a + a') / 2, ℬ) == Tens((a + a') / 2, ℬ)
            @test get_array(rep_inferred_symmetric(a, ℬ)) ≈ (a + a') / 2          # the mean over the symmetry
            @test get_array(rep_inferred_symmetric(a, ℬ)) isa SymmetricTensor{2, d, Float64}
            @test get_array(rep_inferred_symmetric(As, ℬ)) isa SymmetricTensor{4, d, Float64}
            @test get_array(rep_inferred_symmetric(A, ℬ)) ≈ As
            @test rep_inferred_symmetric(As, ℬ) ≈ Tens(As, ℬ)
            @test static_tens(a, ℬ) == Tens(a, ℬ) && static_tens(As, ℬ) == Tens(As, ℬ)
            @test get_basis(static_tens(a, ℬ)) == ℬ
            @test static_tens(SMatrix{d, d}(a), CanonicalBasis{d, Float64}()) == Tens(a)
            @test Base.return_types(static_tens, (Matrix{Float64}, typeof(ℬ))) |> only |> isconcretetype |> !   # two storages
        end
        # Dual numbers and symbols keep their type.
        s = [symbols("s$(min(i, j))$(max(i, j))", real = true) for i in 1:2, j in 1:2]
        @test eltype(symmetric_tens(s, CanonicalBasis{2, Sym}())) <: Sym
        dA = ForwardDiff.derivative(x -> get_array(symmetric_tens([x 2x; 0 x^2], CanonicalBasis{2, Float64}())), 3.0)
        @test dA ≈ [1.0 1.0; 1.0 6.0]
    end

    @testsection "exact membership of a class" begin
        ℬ = random_basis(rng)
        O = TensOrtho(10.0, 6.0, 8.0, 3.0, 2.0, 1.0, 4.0, 5.0, 6.0, ℬ)
        t = Tens(get_array(O))
        @test project_exact(t, :ORTHO, ℬ) ≈ O
        @test project_exact(t, :ORTHO, ℬ) isa TensOrtho{4, 3}
        @test project_exact(t, Val(:ISO)) === nothing
        @test project_exact(t, :CUBIC, ℬ) === nothing
        @test project_exact(t, :ISO; rtol = 10.0) isa TensISO{4, 3}                  # the tolerance is the caller's
        L = iso_stiffness(3.0, 1.0)
        @test project_exact(Tens(get_array(L)), :ISO) ≈ L
        n = random_axis(rng)
        W = ti_from_constants(19.0, 5.0, 11.0, 12.0, 4.0, n; form = :hill)
        @test project_exact(Tens(get_array(W)), :TI, n) ≈ W
        @test project_exact(Tens(get_array(W)), :TI, random_axis(rng)) === nothing
        @test all(iszero, project_exact(Tens(zeros(3, 3, 3, 3)), :ISO))             # the null tensor is of every class
        # The plane.
        O₂ = TensOrtho(10.0, 6.0, 3.0, 2.0, Basis(0.4))
        @test project_exact(Tens(get_array(O₂)), :ORTHO, Basis(0.4)) ≈ O₂
        @test project_exact(Tens(get_array(O₂)), :ISO) === nothing
        # Symbols: a distance cannot be compared; the isotropic tensor is of every class.
        k, μ = symbols("k μ", positive = true)
        Lˢ = iso_stiffness(k, μ)
        @test project_exact(Lˢ, :ISO) == Lˢ
        @test project_exact(Tens(get_array(Lˢ)), :ISO) === nothing
        # Dual numbers are compared by their values.
        d = ForwardDiff.derivative(x -> bulk_modulus(project_exact(Tens(get_array(iso_stiffness(x, 1.0))), :ISO)), 3.0)
        @test d ≈ 1
    end

    @testsection "transversely isotropic constants without the major symmetry" begin
        k, l, m, n, p, lt = 19.0, 5.0, 11.0, 12.0, 4.0, 7.0
        e₃ = (0.0, 0.0, 1.0)
        t = ti_from_constants(k, l, m, n, p, lt, e₃; form = :hill)
        @test t isa TensTI{4, Float64, 6}
        C = get_array(t)
        # The definitions of the moduli, read on the components.
        @test (C[1, 1, 1, 1] + C[1, 1, 2, 2]) / 2 ≈ k && (C[1, 1, 1, 1] - C[1, 1, 2, 2]) / 2 ≈ m
        @test C[1, 1, 3, 3] ≈ l && C[3, 3, 1, 1] ≈ lt
        @test C[3, 3, 3, 3] ≈ n && C[2, 3, 2, 3] ≈ p && C[1, 2, 1, 2] ≈ m
        h = rep_inferred_hill(t)
        @test keys(h) == (:k, :l, :m, :n, :p, :lt)
        @test collect(h) ≈ [k, l, m, n, p, lt]
        @test_throws ArgumentError hill_moduli(t)                                    # five moduli: major symmetry needed
        # The other two forms build the same tensor.
        @test ti_from_constants(k + m, k - m, l, lt, n, p, e₃; form = :components) ≈ t
        s = sqrt(2.0)
        @test ti_from_constants(n, 2k, s * lt, s * l, 2m, 2p, e₃; form = :walpole) ≈ t
        @test_throws ArgumentError ti_from_constants(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, e₃; form = :engineering)
        @test_throws ArgumentError ti_from_constants(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, e₃; form = :voigt)
        # With the major symmetry the sixth modulus is the second, on five or eight coefficients.
        t₅ = ti_from_constants(k, l, m, n, p, e₃; form = :hill)
        @test collect(hill_moduli(t₅, Val(6))) ≈ [k, l, m, n, p, l]
        @test ti_from_constants(k, l, m, n, p, l, e₃; form = :hill) ≈ t₅
        L = iso_stiffness(10.0, 3.0)
        @test collect(hill_moduli(L, Val(6))) ≈ [collect(hill_moduli(L)); hill_moduli(L).l]
        # A product of two stiffnesses has no major symmetry: its moduli from the components.
        ax = random_axis(rng)
        P = ti_from_constants(19.0, 5.0, 11.0, 12.0, 4.0, ax; form = :hill) ⊡ ti_from_constants(8.0, 2.0, 3.0, 15.0, 6.0, ax; form = :hill)
        hP = hill_moduli(P, Val(6))
        @test !(hP.l ≈ hP.lt)
        @test ti_from_constants(hP..., ax; form = :hill) ≈ P
        # The couplings of the rotations alone have no Hill moduli.
        t₈ = TensTI{4}(1.0, 2.0, 0.3, 0.4, 5.0, 6.0, 0.1, 0.2, e₃)
        @test_throws ArgumentError hill_moduli(t₈, Val(6))
        @test collect(hill_moduli(TensTI{4}(1.0, 2.0, 0.3, 0.4, 5.0, 6.0, 0.0, 0.0, e₃), Val(6))) ≈ [1.0, 0.4 / s, 2.5, 1.0, 3.0, 0.3 / s]
        # Derivatives and symbols.
        g = ForwardDiff.gradient(x -> hill_moduli(ti_from_constants(x..., e₃; form = :hill), Val(6)).lt, [k, l, m, n, p, lt])
        @test g ≈ [0, 0, 0, 0, 0, 1]
        ks = symbols("k l m n p l_t", positive = true)
        hs = hill_moduli(ti_from_constants(ks..., e₃; form = :hill), Val(6))
        @test all(iszero, simplify.(collect(hs) .- collect(ks)))
    end
end
