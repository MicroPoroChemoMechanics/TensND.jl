# =============================================================================
#  test_plane_space.jl — a tensor of space in one of its planes, a tensor of
#  the plane in space: `plane_block`, `plane_stress_block`, `space_embedding`,
#  and the kernels of the plane next to those of space (`acoustic_components`
#  in the plane, `sym2_inv`, `sym3_adj`).
#
#  References, independent of the code under test:
#  - the sub-arrays of a component array on the plane indices, taken with index
#    vectors, and the canonical components of a tensor written in a rotated
#    basis by explicit index sums (`reference_contractions.jl`);
#  - the plane-stress stiffness as the inverse of the in-plane block of the
#    inverse of the 6×6 Kelvin-Mandel matrix, the slots written out below;
#  - the structured reductions `plane_strain_stiffness` and
#    `plane_stress_stiffness`;
#  - the out-of-plane fillers written in MeanFieldHomogenization, as matrices;
#  - the acoustic tensor by an explicit contraction, inverses and adjugates by
#    `inv` and `det`;
#  - the plane-strain and plane-stress matrices of an isotropic material in
#    closed form, with symbols.
# =============================================================================

using StaticArrays: SArray, SMatrix, SVector

const PLS_SEED = 20261126
# The axes of the plane normal to e_k in the cyclic order, and the Kelvin-Mandel
# slots of space (11, 22, 33, √2·23, √2·13, √2·12): those of the plane in the
# order (11, 22, √2·12) of the plane, then e_k⊗e_k and the two out-of-plane shears.
const PLS_AXES = ([2, 3], [3, 1], [1, 2])
const PLS_IP = ([2, 3, 4], [3, 1, 5], [1, 2, 6])
const PLS_OP = ([1, 5, 6], [2, 4, 6], [3, 4, 5])

pls_slot3(i, j) = i == j ? i : 9 - i - j
pls_slot2(i, j) = i == j ? i : 3
pls_w(i, j) = i == j ? 1.0 : sqrt(2.0)

"Component array of space of the 6×6 Kelvin-Mandel matrix `M`, by the definition of its components."
pls_array(M) = [M[pls_slot3(i, j), pls_slot3(k, l)] / (pls_w(i, j) * pls_w(k, l)) for i in 1:3, j in 1:3, k in 1:3, l in 1:3]
"Component array of the plane of the 3×3 Kelvin-Mandel matrix `M`."
pls_array2(M) = [M[pls_slot2(i, j), pls_slot2(k, l)] / (pls_w(i, j) * pls_w(k, l)) for i in 1:2, j in 1:2, k in 1:2, l in 1:2]
"6×6 Kelvin-Mandel matrix of a minor-symmetric array of space, by definition."
pls_km(A) = [
    pls_w(p...) * pls_w(q...) * A[p..., q...]
        for p in ((1, 1), (2, 2), (3, 3), (2, 3), (1, 3), (1, 2)), q in ((1, 1), (2, 2), (3, 3), (2, 3), (1, 3), (1, 2))
]

"The minor-symmetric part of an array of order 4, summed by pairs so that it is exactly minor-symmetric."
pls_minor(A) = [((A[i, j, k, l] + A[j, i, k, l]) + (A[i, j, l, k] + A[j, i, l, k])) / 4 for i in 1:3, j in 1:3, k in 1:3, l in 1:3]

"A symmetric positive definite matrix of size `n`."
pls_spd(rng, n) = (G = randn(rng, n, n); G * transpose(G) + n * I)

"Canonical components of the tensor whose components in the basis of columns `R` are `A`."
pls_canonical(A::AbstractArray{<:Any, 4}, R) = ref_rotate(R, A)
pls_canonical(A::AbstractMatrix, R) = R * A * transpose(R)

pls_iszero(M) = all(x -> iszero(simplify(x)), M)

# The out-of-plane fillers of MeanFieldHomogenization, written as there: a
# stiffness of the plane (e₂, e₃) of Kelvin-Mandel matrix `M`, a conductivity
# `K` of that plane, a gradient concentration `S` of that plane.
function pls_mfh_space_stiffness(M)
    s = (M[1, 1] + M[2, 2]) / 2
    o = zero(s)
    return SMatrix{6, 6}(
        s, o, o, o, o, o,
        o, M[1, 1], M[2, 1], M[3, 1], o, o,
        o, M[1, 2], M[2, 2], M[3, 2], o, o,
        o, M[1, 3], M[2, 3], M[3, 3], o, o,
        o, o, o, o, s, o,
        o, o, o, o, o, s,
    )
end
function pls_mfh_space_conductivity(K)
    s = (K[1, 1] + K[2, 2]) / 2
    o = zero(s)
    return SMatrix{3, 3}(s, o, o, o, K[1, 1], K[2, 1], o, K[1, 2], K[2, 2])
end
pls_mfh_hole(S) = SVector{6}(1.0, S[1, 1], S[2, 2], sqrt(2.0) * S[2, 1], 0.0, 0.0)

pls_inferred_block(t) = @inferred plane_block(t, 1)
pls_inferred_stress(t) = @inferred plane_stress_block(t, 2)
pls_inferred_embed2(t) = @inferred space_embedding(t, 1; axial = 1.0)
pls_inferred_embed4(t) = @inferred space_embedding(t, 3; axial = 1.0, shear = 2.0)
pls_inferred_acoustic(C, ζ) = @inferred acoustic_components(C, ζ)

@testsection "Tensors of space in a plane, tensors of the plane in space" begin
    rng = Xoshiro(PLS_SEED)

    @testsection "plane_block: the restriction to a plane" begin
        ℬ = random_basis(rng)
        R = collect(vecbasis(ℬ))
        for k in 1:3
            ix = PLS_AXES[k]
            A = randn(rng, 3, 3, 3, 3)
            As = pls_minor(A)
            for X in (A, As)
                @test Array(get_array(plane_block(Tens(X), k))) == X[ix, ix, ix, ix]
            end
            @test get_array(plane_block(Tens(As), k)) isa SymmetricTensor{4, 2, Float64}
            a = randn(rng, 3, 3)
            as = (a + transpose(a)) / 2
            for x in (a, as)
                @test Array(get_array(plane_block(Tens(x), k))) == x[ix, ix]
            end
            @test get_array(plane_block(Tens(as), k)) isa SymmetricTensor{2, 2, Float64}
            # Written in a rotated basis: the block of its canonical components.
            @test Array(get_array(plane_block(Tens(As, ℬ), k))) ≈ pls_canonical(As, R)[ix, ix, ix, ix]
            @test Array(get_array(plane_block(Tens(a, ℬ), k))) ≈ pls_canonical(a, R)[ix, ix]
            # Through the Kelvin-Mandel matrix: the in-plane block (plane strain).
            M = pls_spd(rng, 6)
            @test mandel(plane_block(Tens(pls_array(M)), k)) ≈ M[PLS_IP[k], PLS_IP[k]]
        end
        @test plane_block(Tens(pls_minor(randn(rng, 3, 3, 3, 3)))) isa AbstractTens{4, 2}
        @test_throws ArgumentError plane_block(TensISO{3}(30.0, 8.0), 4)
        @test_throws ArgumentError plane_block(Tens(randn(rng, 3, 3)), 0)
    end

    @testsection "structured classes: plane_strain_stiffness and plane_stress_stiffness" begin
        CB3 = CanonicalBasis{3, Float64}()
        N = [160.0 55.0 48.0; 31.0 120.0 42.0; 22.0 17.0 95.0]
        Ns = (N + transpose(N)) / 2
        S = (33.0, 29.0, 27.0)
        ℂ9 = TensOrtho(Ns[1, 1], Ns[2, 2], Ns[3, 3], Ns[1, 2], Ns[1, 3], Ns[2, 3], S..., CB3)
        ℂ12 = TensOrtho(N[1, 1], N[2, 2], N[3, 3], N[1, 2], N[1, 3], N[2, 3], S..., N[2, 1], N[3, 1], N[3, 2], CB3)
        Cc = cubic_from_constants(170.0, 120.0, 75.0, CB3)
        Ci = TensISO{3}(30.0, 8.0)
        for k in 1:3
            for ℂ in (ℂ9, ℂ12)
                @test mandel(plane_block(ℂ, k)) ≈ mandel(plane_strain_stiffness(ℂ, k))
                @test mandel(plane_stress_block(ℂ, k)) ≈ mandel(plane_stress_stiffness(ℂ, k))
            end
            @test mandel(plane_block(Cc, k)) ≈ mandel(plane_strain_stiffness(Cc))
            @test mandel(plane_stress_block(Cc, k)) ≈ mandel(plane_stress_stiffness(Cc))
            @test mandel(plane_block(Ci, k)) ≈ mandel(plane_strain_stiffness(Ci))
            @test mandel(plane_stress_block(Ci, k)) ≈ mandel(plane_stress_stiffness(Ci))
        end
        # Transversely isotropic about e₃ in the plane normal to it, and about e₁
        # in the plane (e₁, e₂) that contains the axis.
        ℂti3 = ti_from_constants(10.0, 4.0, 3.0, 20.0, 2.5, [0.0, 0.0, 1.0]; form = :hill)
        ℂti1 = ti_from_constants(10.0, 4.0, 3.0, 20.0, 2.5, [1.0, 0.0, 0.0]; form = :hill)
        @test mandel(plane_block(ℂti3, 3)) ≈ mandel(plane_strain_stiffness(ℂti3))
        @test mandel(plane_block(ℂti1, 3)) ≈ mandel(plane_strain_stiffness(ℂti1, Val(:axial)))
        @test mandel(plane_stress_block(ℂti1, 3)) ≈ mandel(plane_stress_stiffness(ℂti1, Val(:axial)))
        @test get_array(plane_block(ℂti1, 3)) isa SymmetricTensor{4, 2, Float64}
    end

    @testsection "plane_stress_block: the inverse of the in-plane compliance" begin
        ℬ = random_basis(rng)
        R = collect(vecbasis(ℬ))
        for k in 1:3
            p = PLS_IP[k]
            M = pls_spd(rng, 6)
            H = randn(rng, 6, 6)
            Mn = M + (H - transpose(H)) / 2               # no major symmetry
            for X in (M, Mn)
                @test mandel(plane_stress_block(Tens(pls_array(X)), k)) ≈ inv(inv(X)[p, p])
            end
            # Written in a rotated basis.
            Cc = pls_canonical(pls_array(M), R)
            @test mandel(plane_stress_block(Tens(pls_array(M), ℬ), k)) ≈ inv(inv(pls_km(Cc))[p, p])
        end
    end

    @testsection "space_embedding" begin
        for k in 1:3
            ix, p, (n, s₁, s₂) = PLS_AXES[k], PLS_IP[k], PLS_OP[k]
            @test n == k
            # Order 2, with and without symmetry.
            a = randn(rng, 2, 2)
            for x in (a, (a + transpose(a)) / 2), α in (0.0, 1.0, 2.7)
                t = Tens(x)
                X = zeros(3, 3)
                X[ix, ix] = x
                X[k, k] = α
                E = space_embedding(t, k; axial = α)
                @test Array(get_array(E)) == X
                @test Array(get_array(plane_block(E, k))) == x
            end
            X = zeros(3, 3)
            X[ix, ix] = a
            @test Array(get_array(space_embedding(Tens(a), k))) == X
            # Order 4: the Kelvin-Mandel matrix assembled slot by slot.
            M3 = pls_spd(rng, 3)
            t = Tens(pls_array2(M3))
            for (α, sh) in ((0.0, 0.0), (1.0, 2.0), (3.5, 0.7))
                X = zeros(6, 6)
                X[p, p] = M3
                X[n, n] = α
                X[s₁, s₁] = X[s₂, s₂] = sh
                E = space_embedding(t, k; axial = α, shear = sh)
                @test get_array(E) isa SymmetricTensor{4, 3, Float64}
                @test mandel(E) ≈ X
                @test get_array(plane_block(E, k)) == get_array(t)
            end
            X = zeros(6, 6)
            X[p, p] = M3
            @test mandel(space_embedding(t, k)) ≈ X
        end
        # The fillers of MeanFieldHomogenization, in the plane (e₂, e₃).
        M3 = pls_spd(rng, 3)
        s = (M3[1, 1] + M3[2, 2]) / 2
        @test mandel(space_embedding(Tens(pls_array2(M3)), 1; axial = s, shear = s)) ≈ pls_mfh_space_stiffness(M3)
        K = pls_spd(rng, 2)
        @test Array(get_array(space_embedding(Tens(K), 1; axial = (K[1, 1] + K[2, 2]) / 2))) ≈ pls_mfh_space_conductivity(K)
        Sg = pls_spd(rng, 2)
        @test mandel(space_embedding(Tens(Sg), 1; axial = 1.0)) ≈ pls_mfh_hole(Sg)
        P2 = pls_spd(rng, 2)
        @test Array(get_array(space_embedding(Tens(P2), 1))) ≈ cat(0.0, P2; dims = (1, 2))
        @test_throws ArgumentError space_embedding(Tens(P2), 4)
    end

    @testsection "acoustic tensor of the plane, sym2_inv, sym3_adj" begin
        ζ = randn(rng, 2)
        C = pls_array2(pls_spd(rng, 3))
        A = randn(rng, 2, 2, 2, 2)                         # no symmetry
        for X in (C, A)
            K = ref_acoustic(ζ, X, ζ)
            ref = (K[1, 1], K[2, 2], K[1, 2])
            @test all(acoustic_components(X, ζ, Val(2)) .≈ ref)
            @test all(acoustic_components(SArray{Tuple{2, 2, 2, 2}}(X), ζ) .≈ ref)
            @test all(acoustic_components(Tensor{4, 2}(X), Tuple(ζ)) .≈ ref)
        end
        @test all(acoustic_components(SymmetricTensor{4, 2}(C), ζ) .≈ acoustic_components(C, ζ, Val(2)))
        # Space keeps its method.
        ζ3 = randn(rng, 3)
        C3 = pls_array(pls_spd(rng, 6))
        K3 = ref_acoustic(ζ3, C3, ζ3)
        @test all(acoustic_components(C3, ζ3) .≈ (K3[1, 1], K3[2, 2], K3[3, 3], K3[1, 2], K3[1, 3], K3[2, 3]))
        # Inverses and adjugate.
        K2 = pls_spd(rng, 2)
        @test collect(sym2_inv(K2[1, 1], K2[2, 2], K2[1, 2])) ≈ vec(inv(K2))
        K = pls_spd(rng, 3)
        e = (K[1, 1], K[2, 2], K[3, 3], K[1, 2], K[1, 3], K[2, 3])
        adj = det(K) * inv(K)
        @test collect(sym3_adj(e...)) ≈ [adj[1, 1], adj[2, 2], adj[3, 3], adj[1, 2], adj[1, 3], adj[2, 3]]
        @test collect(sym3_inv(e...)) ≈ vec(inv(K))
        @test pls_inferred_acoustic(SymmetricTensor{4, 2}(C), Tuple(ζ)) isa NTuple{3, Float64}
        @test pls_inferred_acoustic(SymmetricTensor{4, 2}(C), (ForwardDiff.Dual(1.0, 1.0), 0.5)) isa NTuple{3, ForwardDiff.Dual}
    end

    @testsection "inference" begin
        ℬ = random_basis(rng)
        As = pls_minor(randn(rng, 3, 3, 3, 3))
        for t in (Tens(As), Tens(As, ℬ), TensISO{3}(30.0, 8.0))
            @test pls_inferred_block(t) isa AbstractTens{4, 2, Float64}
            @test pls_inferred_stress(t) isa AbstractTens{4, 2, Float64}
        end
        a = randn(rng, 3, 3)
        @test pls_inferred_block(Tens((a + transpose(a)) / 2)) isa AbstractTens{2, 2, Float64}
        @test pls_inferred_embed2(Tens(pls_spd(rng, 2))) isa AbstractTens{2, 3, Float64}
        @test pls_inferred_embed4(Tens(pls_array2(pls_spd(rng, 3)))) isa AbstractTens{4, 3, Float64}
        @test pls_inferred_stress(Tens(ForwardDiff.Dual.(As, 1.0))) isa AbstractTens{4, 2, <:ForwardDiff.Dual}
    end

    @testsection "automatic differentiation" begin
        M, G = pls_spd(rng, 6), randn(rng, 6, 6)
        dM = (G + transpose(G)) / 2
        for k in 1:3
            p = PLS_IP[k]
            f = x -> Matrix(mandel(plane_stress_block(Tens(pls_array(M + x * dM)), k)))
            g = x -> inv(inv(M + x * dM)[p, p])
            @test ForwardDiff.derivative(f, 0.0) ≈ ForwardDiff.derivative(g, 0.0)
            @test ForwardDiff.derivative(y -> ForwardDiff.derivative(f, y), 0.0) ≈
                ForwardDiff.derivative(y -> ForwardDiff.derivative(g, y), 0.0)
        end
        # The keywords of the embedding carry a derivative too.
        t = Tens(pls_array2(pls_spd(rng, 3)))
        @test ForwardDiff.derivative(x -> mandel(space_embedding(t, 2; axial = x, shear = 3x))[2, 2], 1.5) ≈ 1
        @test ForwardDiff.derivative(x -> mandel(space_embedding(t, 2; axial = x, shear = 3x))[4, 4], 1.5) ≈ 3
    end

    @testsection "symbolic: isotropic plane strain and plane stress" begin
        E, ν = symbols("E ν", positive = true)
        C = iso_stiffness_young_poisson(E, ν)
        Pε = E / ((1 + ν) * (1 - 2ν)) * [(1 - ν) ν 0; ν (1 - ν) 0; 0 0 (1 - 2ν)]
        Pσ = E / (1 - ν^2) * [1 ν 0; ν 1 0; 0 0 (1 - ν)]
        for k in 1:3
            @test pls_iszero(mandel(plane_block(C, k)) - Pε)
            @test pls_iszero(mandel(plane_stress_block(C, k)) - Pσ)
        end
    end
end
