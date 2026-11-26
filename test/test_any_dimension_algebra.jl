# =============================================================================
#  test_any_dimension_algebra.jl — the tensor algebra beyond dimension 3:
#  transposes, products with an isotropic tensor, the isotropic symmetry
#  basis, random generic tensors, the vectors of the canonical basis and the
#  Cartesian coordinate systems.
#
#  Tensors.jl stops at dimension 3: beyond, the components are held by plain
#  arrays. The references are explicit sums over the indices (OMEinsum),
#  Tensors.jl itself in dimensions 2 and 3, where the array kernels must give
#  what it gives, and derivatives taken by SymPy.
# =============================================================================

using StaticArrays: SMatrix, SVector

const ADA_SEED = 20261127

"The full array of `α 𝕁 + β 𝕂` in dimension `d`, from its definition."
function ada_iso_array(α, β, d)
    δ(i, j) = i == j ? 1 : 0
    J(i, j, k, l) = δ(i, j) * δ(k, l) / d
    Is(i, j, k, l) = (δ(i, k) * δ(j, l) + δ(i, l) * δ(j, k)) / 2
    return [α * J(i, j, k, l) + β * (Is(i, j, k, l) - J(i, j, k, l)) for i in 1:d, j in 1:d, k in 1:d, l in 1:d]
end

ada_inferred_basis(L) = @inferred symmetry_basis(L)
ada_inferred_layout(L) = @inferred symmetry_layout(L)

@testsection "Tensor algebra beyond dimension 3" begin
    rng = Xoshiro(ADA_SEED)

    @testsection "transposes on array storage" begin
        for d in (4, 5)
            a, A = randn(rng, d, d), randn(rng, d, d, d, d)
            ta, tA = Tens(a), Tens(A)
            @test transpose(ta) isa AbstractTens{2, d}
            @test get_array(transpose(ta)) == permutedims(a)
            @test get_array(transpose(tA)) == permutedims(A, (2, 1, 4, 3))
            @test get_array(minortranspose(tA)) == permutedims(A, (2, 1, 4, 3))
            @test get_array(majortranspose(tA)) == permutedims(A, (3, 4, 1, 2))
            @test transpose(transpose(tA)) == tA
            @test majortranspose(majortranspose(tA)) == tA
            # In a rotated basis the components are permuted in that basis.
            ℬ = random_basis(rng, Val(d))
            @test transpose(Tens(a, ℬ)) == Tens(permutedims(a), ℬ)
            @test majortranspose(Tens(A, ℬ)) == Tens(permutedims(A, (3, 4, 1, 2)), ℬ)
            @test minortranspose(Tens(A, ℬ)) == Tens(permutedims(A, (2, 1, 4, 3)), ℬ)
            @test transpose(Tens(A, ℬ)) == Tens(permutedims(A, (2, 1, 4, 3)), ℬ)
            # In a basis that is not orthonormal the variances follow the indices:
            # the canonical components of the transpose are the transposed ones.
            𝒢 = Basis(randn(rng, d, d) + d * I)
            for var in ((:cov, :cont), (:cont, :cont), (:cov, :cov))
                t = Tens(a, 𝒢, var)
                @test get_var(transpose(t)) == reverse(var)
                @test canonical_components(transpose(t)) ≈ permutedims(canonical_components(t))
            end
            T4 = Tens(A, 𝒢, (:cov, :cont, :cont, :cov))
            @test get_var(majortranspose(T4)) == (:cont, :cov, :cov, :cont)
            @test canonical_components(majortranspose(T4)) ≈ permutedims(canonical_components(T4), (3, 4, 1, 2))
            @test get_var(minortranspose(T4)) == (:cont, :cov, :cov, :cont)
            @test canonical_components(minortranspose(T4)) ≈ permutedims(canonical_components(T4), (2, 1, 4, 3))
            @test canonical_components(transpose(T4)) ≈ permutedims(canonical_components(T4), (2, 1, 4, 3))
        end
        # The array kernels give what Tensors.jl gives in dimensions 2 and 3.
        for d in (2, 3)
            x2, x4 = Tensor{2, d}(randn(rng, d, d)), Tensor{4, d}(randn(rng, d, d, d, d))
            @test TensND._transpose_array(Array(x2)) == Array(transpose(x2))
            @test TensND._transpose_array(Array(x4)) == Array(transpose(x4))
            @test TensND._minortranspose_array(Array(x4)) == Array(minortranspose(x4))
            @test TensND._majortranspose_array(Array(x4)) == Array(majortranspose(x4))
        end
        # Dual numbers and symbols go through: a transpose only moves components.
        A = randn(rng, 4, 4, 4, 4)
        @test ForwardDiff.derivative(x -> get_array(majortranspose(Tens(x .* A))), 0.7) == permutedims(A, (3, 4, 1, 2))
        s = [symbols("s$i$j", real = true) for i in 1:4, j in 1:4]
        @test get_array(transpose(Tens(s))) == permutedims(s)
    end

    @testsection "products with an isotropic tensor" begin
        for d in (4, 5)
            α, β = 3.0 + rand(rng), 2.0 + rand(rng)
            L, Lfull = TensISO{d}(α, β), ada_iso_array(α, β, d)
            a, A = randn(rng, d, d), randn(rng, d, d, d, d)
            @test Array(L) ≈ Lfull
            # Order 2: a tensor comes back, not a bare matrix.
            @test L ⊡ Tens(a) isa AbstractTens{2, d}
            @test Tens(a) ⊡ L isa AbstractTens{2, d}
            @test get_array(L ⊡ Tens(a)) ≈ ein"ijkl,kl->ij"(Lfull, a)
            @test get_array(Tens(a) ⊡ L) ≈ ein"ij,ijkl->kl"(a, Lfull)
            ℬ = random_basis(rng, Val(d))
            @test L ⊡ Tens(a, ℬ) isa AbstractTens{2, d}
            @test get_array(L ⊡ Tens(a, ℬ)) ≈ ein"ijkl,kl->ij"(Lfull, a)
            # Order 4.
            @test get_array(L ⊡ Tens(A)) ≈ ein"ijkl,klmn->ijmn"(Lfull, A)
            @test get_array(Tens(A) ⊡ L) ≈ ein"ijkl,klmn->ijmn"(A, Lfull)
            @test (L ⊙ Tens(A)) ≈ ein"ijkl,ijkl->"(Lfull, A)[]
            @test (Tens(A) ⊙ L) ≈ ein"ijkl,ijkl->"(A, Lfull)[]
            # Sum, inverse, Kelvin-Mandel matrix.
            @test get_array(L + Tens(A)) ≈ Lfull + A
            @test Array(inv(L)) ≈ ada_iso_array(1 / α, 1 / β, d)
            @test mandel(L) ≈ mandel(Tens(Lfull))
        end
        # The full contraction of two isotropic tensors, in closed form.
        for d in 2:5
            A, B = TensISO{d}(1.3, 0.7), TensISO{d}(-0.4, 2.1)
            @test (@inferred A ⊙ B) ≈ ein"ijkl,ijkl->"(ada_iso_array(1.3, 0.7, d), ada_iso_array(-0.4, 2.1, d))[]
            @test ForwardDiff.gradient(x -> TensISO{d}(x[1], x[2]) ⊙ B, [1.3, 0.7]) ≈ [-0.4, (d * (d + 1) ÷ 2 - 1) * 2.1]
        end
        # Derivatives with respect to the moduli, first and second order.
        a = randn(rng, 4, 4)
        energy(x) = Tens(a) ⊡ (TensISO{4}(x[1], x[2]) ⊡ Tens(a))
        sa = (a + a') / 2
        g = ForwardDiff.gradient(energy, [3.0, 2.0])
        @test g ≈ [tr(a)^2 / 4, sum(abs2, sa) - tr(a)^2 / 4]
        @test ForwardDiff.hessian(x -> energy(x)^2, [3.0, 2.0]) ≈ 2 * g * g'
        # With symbols.
        α, β = symbols("α β", positive = true)
        Lˢ = TensISO{4}(α, β)
        e = Tens([symbols("ε$(min(i, j))$(max(i, j))", real = true) for i in 1:4, j in 1:4])
        σ = Lˢ ⊡ e
        @test σ isa AbstractTens{2, 4}
        @test simplify(tr(σ) - α * tr(e)) == 0
        @test simplify(σ[1, 2] - β * e[1, 2]) == 0
        @test simplify(Lˢ ⊙ Lˢ - (α^2 + 9β^2)) == 0
    end

    @testsection "isotropic symmetry basis in dimension d" begin
        for d in (4, 5)
            n = d * (d + 1) ÷ 2
            L = TensISO{d}(3.0, 2.0)
            Q = ada_inferred_basis(L)
            @test Q isa SMatrix{n, n, Float64}
            @test Q' * Q ≈ I
            @test Q[:, 1] ≈ mandel(tens_Id2(Val(d), Val(Float64))) / sqrt(d)
            @test Q' * mandel(L) * Q ≈ Diagonal([3.0; fill(2.0, n - 1)])
            layout = ada_inferred_layout(L)
            @test layout == ((SVector(1),), ntuple(i -> SVector(i + 1), n - 1))
            @test symmetry_matrix(L, symmetry_blocks(L)) ≈ mandel(L)
            # The blocks of any matrix: the projections on 𝕁 and on 𝕂, read on the first copy.
            M = mandel(Tens(ada_iso_array(1.7, 0.3, d)))
            b = symmetry_blocks(L, M)
            @test b[1][1] ≈ 1.7 && b[2][1] ≈ 0.3
            @test ForwardDiff.derivative(x -> symmetry_blocks(L, mandel(TensISO{d}(x, 2.0)))[1][1], 3.0) ≈ 1
        end
        # In the plane and in space the general basis is the one written out.
        @test TensND._iso_mandel_basis(Float64, Val(2)) ≈ symmetry_basis(TensISO{2}(3.0, 2.0))
        @test TensND._iso_mandel_basis(Float64, Val(3)) ≈ symmetry_basis(TensISO{3}(3.0, 2.0))
        Qb = TensND._iso_mandel_basis(BigFloat, Val(4))
        @test eltype(Qb) === BigFloat
        @test norm(Qb' * Qb - I) < 1.0e-70
    end

    @testsection "random generic tensors" begin
        for d in (4, 5)
            n = d * (d + 1) ÷ 2
            a = rand(rng, Tens{2, d})
            @test a isa AbstractTens{2, d, Float64}
            @test get_array(a) == permutedims(get_array(a))
            @test rand(rng, Tens{2, d, Float32}) isa AbstractTens{2, d, Float32}
            g = rand(rng, Tens{2, d}; symmetry = Val(:none))
            @test get_array(g) != permutedims(get_array(g))
            p = rand(rng, Tens{2, d}; definite = true)
            @test all(>(0), eigvals(Symmetric(get_array(p))))
            A = rand(rng, Tens{4, d})
            @test A isa AbstractTens{4, d, Float64}
            @test size(mandel(A)) == (n, n)
            @test mandel(A) ≈ mandel(A)'
            @test all(>(0), eigvals(Symmetric(mandel(rand(rng, Tens{4, d}; definite = true)))))
            B = rand(rng, Tens{4, d}; symmetry = Val(:minor))
            @test get_array(B) ≈ permutedims(get_array(B), (2, 1, 4, 3))
            @test !(mandel(B) ≈ mandel(B)')
            C = rand(rng, Tens{4, d}; symmetry = Val(:none))
            @test size(get_array(C)) == (d, d, d, d)
            @test !(get_array(C) ≈ permutedims(get_array(C), (2, 1, 3, 4)))
            @test_throws ArgumentError rand(rng, Tens{2, d}; symmetry = Val(:none), definite = true)
        end
        # Up to dimension 3 the components stay in Tensors.jl.
        @test get_array(rand(rng, Tens{2, 3})) isa SymmetricTensor{2, 3, Float64}
        @test get_array(rand(rng, Tens{2, 2}; symmetry = Val(:none))) isa Tensor{2, 2, Float64}
        @test get_array(rand(rng, Tens{4, 3}; symmetry = Val(:none))) isa Tensor{4, 3, Float64}
    end

    @testsection "vectors of the canonical basis" begin
        for d in (4, 5)
            e = ntuple(i -> 𝐞(i, d, Float64), d)
            @test all(v -> v isa AbstractTens{1, d, Float64}, e)
            @test [e[i] ⋅ e[j] for i in 1:d, j in 1:d] == I
            @test get_array(𝐞(Val(2), Val(d), Val(Float64))) == [j == 2 for j in 1:d]
            @test eltype(𝐞(1, d)) <: Sym                       # the default number type of the session
            @test sum(i * e[i] for i in 1:d) == Tens(collect(1.0:d))
            coords, vectors, ℬ = init_cartesian(d)
            @test coords == ntuple(i -> symbols("x$i", real = true), d)
            @test length(vectors) == d && ℬ isa CanonicalBasis{d}
            @test init_cartesian(Val(d)) == init_cartesian(d)
        end
        @test get_array(𝐞(1, 3, Float64)) isa Vec{3, Float64}
        @test get_array(𝐞(2, 2, Float64)) isa Vec{2, Float64}
        # Without a symbolic package the default coordinates are refused.
        @test_throws ArgumentError TensND._default_cartesian_coords(Val(4), Float64)
    end

    @testsection "Cartesian coordinate system in dimension 4, with symbols" begin
        CS = coorsys_cartesian(Val(4))
        x = get_coords(CS)
        @test x == ntuple(i -> symbols("x$i", real = true), 4)
        @test all(isone, lame_coefficients(CS))
        f = x[1]^2 * x[4] + x[2] * x[3]^3 + sin(x[1] * x[2])
        ∂(g, i) = diff(g, x[i])
        @test all(i -> simplify(grad(f, CS)[i] - ∂(f, i)) == 0, 1:4)
        @test simplify(laplacian(f, CS) - sum(∂(∂(f, i), i) for i in 1:4)) == 0
        H = hess(f, CS)
        @test all(simplify(H[i, j] - ∂(∂(f, i), j)) == 0 for i in 1:4, j in 1:4)
        𝐞ᵢ = unitvec(CS)
        u = ntuple(i -> x[i]^2 * x[mod1(i + 1, 4)], 4)
        𝐮 = sum(u[i] * 𝐞ᵢ[i] for i in 1:4)
        G = grad(𝐮, CS)
        @test all(simplify(G[i, j] - ∂(u[i], j)) == 0 for i in 1:4, j in 1:4)
        E = symgrad(𝐮, CS)
        @test all(simplify(E[i, j] - (∂(u[i], j) + ∂(u[j], i)) / 2) == 0 for i in 1:4, j in 1:4)
        @test simplify(divergence(𝐮, CS) - sum(∂(u[i], i) for i in 1:4)) == 0
        # The plane, by its dimension or by its coordinates.
        @test get_coords(coorsys_cartesian(Val(2))) == symbols("x y", real = true)
        @test get_coords(coorsys_cartesian(Val(3))) == get_coords(coorsys_cartesian())
    end

    @testsection "Pointwise Cartesian coordinate system in any dimension" begin
        for d in (2, 4)
            CS = coorsys_cartesian(PointwiseGeometry, Val(d))
            @test CS isa CoordinateSystem{d, Float64}
            x₀ = collect(1.0:d) ./ 3
            f(x) = x[1]^2 * x[d] + sin(x[1] * x[2])
            ∇f = ForwardDiff.gradient(f, x₀)
            @test get_array(grad(f, CS)(x₀)) ≈ ∇f
            @test laplacian(f, CS)(x₀) ≈ tr(ForwardDiff.hessian(f, x₀))
            @test get_array(hess(f, CS)(x₀)) ≈ ForwardDiff.hessian(f, x₀)
            u(x) = [x[i]^2 * x[mod1(i + 1, d)] for i in 1:d]
            J = ForwardDiff.jacobian(u, x₀)
            @test get_array(grad(u, CS)(x₀)) ≈ J
            @test get_array(symgrad(u, CS)(x₀)) ≈ (J + J') / 2
            @test divergence(u, CS)(x₀) ≈ tr(J)
            # A derivative with respect to a parameter of the field, through the operator.
            @test ForwardDiff.derivative(p -> laplacian(x -> p * x[1]^2 * x[d], CS)(x₀), 1.5) ≈ 2 * x₀[d]
        end
        @test coorsys_cartesian(PointwiseGeometry, Val(2), Float32) isa CoordinateSystem{2, Float32}
        # Without a dimension the system is three-dimensional.
        CS₃ = coorsys_cartesian(PointwiseGeometry)
        @test CS₃ isa CoordinateSystem{3, Float64}
        g(x) = x[1] * x[2]^2 * x[3]^3
        @test get_array(grad(g, CS₃)([1.0, 2.0, 3.0])) == get_array(grad(g, coorsys_cartesian(PointwiseGeometry, Val(3)))([1.0, 2.0, 3.0]))
        @test get_array(grad(g, CS₃)([1.0, 2.0, 3.0])) ≈ [108.0, 108.0, 108.0]
    end
end
