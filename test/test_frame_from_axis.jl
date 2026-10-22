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
