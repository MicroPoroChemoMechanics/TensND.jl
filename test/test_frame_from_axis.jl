# frame_from_axis: a right-handed orthonormal frame about an axis, by
# Gram–Schmidt of a reference vector.

@testsection "frame_from_axis" begin
    rng = Xoshiro(20261024)
    frame(ℬ) = Matrix(vecbasis(ℬ))
    for _ in 1:5
        n = randn(rng, 3)
        a = randn(rng, 3)
        for ref in (nothing, a)
            R = frame(frame_from_axis(n; in_plane = ref))
            @test transpose(R) * R ≈ I
            @test det(R) ≈ 1
            @test R[:, 3] ≈ n / norm(n)
            # ℓ is in the plane of the reference and the axis, on its side.
            r = ref === nothing ? R[:, 1] : a
            @test abs(dot(cross(R[:, 1], r), R[:, 3])) < 1.0e-12 * norm(r)
            @test dot(R[:, 1], r) > 0
        end
    end
    # The default reference is the canonical axis least aligned with n.
    @test frame(frame_from_axis([0.1, 0.9, 0.3]))[:, 1] ≈ normalize([1, 0, 0] - 0.1 / norm([0.1, 0.9, 0.3])^2 * [0.1, 0.9, 0.3])
    # The canonical axis gives the canonical basis, an integer axis a float frame.
    @test frame_from_axis((0, 0, 1)) isa CanonicalBasis
    @test frame_from_axis(Tens([0.0, 0.0, 2.0])) isa CanonicalBasis
    @test eltype(frame_from_axis([1, 2, 2])) == Float64
    @test_throws ArgumentError frame_from_axis([0.0, 0.0, 0.0])
    @test_throws ArgumentError frame_from_axis([0.0, 0.0, 1.0]; in_plane = [0.0, 0.0, 3.0])

    # Derivatives: the axis vector moves with its argument.
    v(x) = frame(frame_from_axis([cos(x), sin(x), 0.5]))[:, 3]
    x₀ = 0.7
    @test ForwardDiff.derivative(v, x₀) ≈ ForwardDiff.derivative(x -> [cos(x), sin(x), 0.5] / norm([cos(x), sin(x), 0.5]), x₀)

    # Symbolic: an exact frame; an axis identically along e₁ takes e₂ as reference.
    TensND.Symbolics.@variables θ
    Rs = frame(frame_from_axis((sin(θ), 0, cos(θ))))
    @test isequal(TensND.Symbolics.simplify(Rs[1, 3] - sin(θ)), 0)
    @syms s::positive
    Re = frame(frame_from_axis((s, 0, 0)))
    @test Re[:, 3] == Sym[1, 0, 0]
    @test Re[:, 1] == Sym[0, 1, 0]
end

# permute_basis: the vectors of a basis permuted, its orientation kept by
# reversing one vector after an odd permutation.
@testsection "permute_basis" begin
    rng = Xoshiro(20261023)
    cols(ℬ) = Matrix(vecbasis(ℬ))
    perms = ((1, 2, 3), (2, 3, 1), (3, 1, 2), (2, 1, 3), (1, 3, 2), (3, 2, 1))
    for _ in 1:3
        Q = Matrix(qr(randn(rng, 3, 3)).Q)
        R = det(Q) < 0 ? Q * Diagonal([1.0, 1.0, -1.0]) : Q
        L = R * Diagonal([1.0, 1.0, -1.0])            # a left-handed basis too
        for M in (R, L), σ in perms, neg in (3, 2)
            P = cols(permute_basis(RotatedBasis(M), σ; negate = neg))
            odd = σ in ((2, 1, 3), (1, 3, 2), (3, 2, 1))
            # Exactly the columns asked, the `negate`-th reversed after an odd
            # permutation, and the orientation of the original basis.
            @test all(P[:, k] == (odd && k == neg ? -M[:, σ[k]] : M[:, σ[k]]) for k in 1:3)
            @test det(P) ≈ det(M)
        end
    end
    # The identity of the canonical basis is the canonical basis.
    @test permute_basis(CanonicalBasis{3, Float64}(), (1, 2, 3)) isa CanonicalBasis{3, Float64}
    @test cols(permute_basis(CanonicalBasis{3, Float64}(), (2, 1, 3))) == [0.0 1 0; 1 0 0; 0 0 -1]
    # Dimension 2: the swap, the second vector reversed.
    @test cols(permute_basis(RotatedBasis(0.3), (2, 1))) == hcat(cols(RotatedBasis(0.3))[:, 2], -cols(RotatedBasis(0.3))[:, 1])
    @test_throws ArgumentError permute_basis(RotatedBasis(0.3, 0.2, 0.1), (1, 1, 3))
    @test_throws ArgumentError permute_basis(RotatedBasis(0.3, 0.2, 0.1), (2, 1, 3); negate = 4)
    # A dual frame: the derivative of the permuted basis is the permuted
    # derivative.
    D(f) = ForwardDiff.derivative(f, 0.7)
    dR = D(θ -> cols(RotatedBasis(θ, 0.4, 1.1)))
    @test D(θ -> cols(permute_basis(RotatedBasis(θ, 0.4, 1.1), (3, 2, 1)))) == hcat(dR[:, 3], dR[:, 2], -dR[:, 1])
    # Symbolic frames, exactly.
    @syms θ::real ϕ::real
    Rs = cols(RotatedBasis(θ, ϕ, 0))
    @test cols(permute_basis(RotatedBasis(θ, ϕ, 0), (2, 3, 1))) == hcat(Rs[:, 2], Rs[:, 3], Rs[:, 1])
    TensND.Symbolics.@variables α
    Rn = cols(RotatedBasis(α))
    @test isequal(cols(permute_basis(RotatedBasis(α), (2, 1))), hcat(Rn[:, 2], -Rn[:, 1]))
end
