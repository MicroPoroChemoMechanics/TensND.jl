# The sum and the difference of two structured tensors, for every pair of
# classes, on frames and axes drawn independently: the canonical components of
# the result are the sums and differences of those of the operands, whatever
# route the pair takes (its own class when the references agree, the generic
# one otherwise), and no pair refuses. The same with dual coefficients, the
# derivative of the result being the operand the coefficients scale. Pairs on
# shared references keep their class, and so do pairs whose references share a
# symmetry of a class: frames equal up to a signed permutation of the axes,
# opposite or orthogonal TI axes, a TI axis along an axis of a cube.

using TensND, Test, LinearAlgebra, ForwardDiff, Random

_sb_rng = MersenneTwister(20261029)
_sb_frame() = random_basis(_sb_rng)       # uniform orientations
_sb_axis() = random_axis(_sb_rng)

# One tensor of each class of order 4, its coefficients scaled by x, on its own
# random frame or axis; then of order 2.
_sb_order4(x, ℬ, n) = (
    TensISO{3}(3.0x, 2.0x),
    tens_ti(10.0x, 3.0x, 2.5x, 12.0x, 2.0x, n),
    TensTI{4, typeof(x), 6}((3.0x, 2.0x, 0.7x, 0.9x, 4.0x, 1.5x), n),
    TensOrtho(9.0x, 8.0x, 7.0x, 1.0x, 2.0x, 3.0x, 4.0x, 5.0x, 6.0x, ℬ),
    TensOrtho(160.0x, 120.0x, 95.0x, 55.0x, 48.0x, 42.0x, 33.0x, 29.0x, 27.0x, 31.0x, 22.0x, 17.0x, ℬ),
    tens_cubic(10.0x, 4.0x, 2.0x, ℬ),
)
_sb_order2(x, ℬ, n) = (TensISO{3}(3.0x), TensTI{2}(1.5x, 9.0x, n), TensOrtho{2}(3.0x, 5.0x, 7.0x, ℬ))

_sb_close(a, b) = maximum(abs, a - b) ≤ 1.0e-12 * max(1, maximum(abs, b))

@testset "structured ± structured on any references" begin
    for builder in (_sb_order4, _sb_order2), _ in 1:3
        As = builder(1.0, _sb_frame(), _sb_axis())
        Bs = builder(1.0, _sb_frame(), _sb_axis())
        for A in As, B in Bs
            a, b = Array(components_canon(A)), Array(components_canon(B))
            @test _sb_close(Array(components_canon(A + B)), a + b)
            @test _sb_close(Array(components_canon(A - B)), a - b)
        end
    end
end

@testset "shared references keep the class" begin
    ℬ, n = _sb_frame(), _sb_axis()
    for (A, B) in zip(_sb_order4(1.0, ℬ, n), _sb_order4(2.0, ℬ, n))
        @test typeof(A + B).name === typeof(A).name
    end
end

@testset "derivatives through the generic route" begin
    ℬ₁, ℬ₂, n₁, n₂ = _sb_frame(), _sb_frame(), _sb_axis(), _sb_axis()
    for (k, l) in Iterators.product(1:6, 1:6)
        B = _sb_order4(1.0, ℬ₂, n₂)[l]
        # d/dx (x A + B) = A: the derivative carried by every component.
        d = ForwardDiff.derivative(x -> Array(components_canon(_sb_order4(x, ℬ₁, n₁)[k] + B)), 1.0)
        @test _sb_close(d, Array(components_canon(_sb_order4(1.0, ℬ₁, n₁)[k])))
    end
end

@testset "shared symmetries keep a class" begin
    ℬ = _sb_frame()
    ℬp = permute_basis(ℬ, (2, 3, 1))              # the same axes, permuted, one negated
    E = Matrix(TensND.vecbasis(ℬ, :cov))
    same(A, B, op) = _sb_close(Array(components_canon(op(A, B))), op(Array(components_canon(A)), Array(components_canon(B))))
    O9(ℱ) = TensOrtho(9.0, 8.0, 7.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, ℱ)
    O12(ℱ) = TensOrtho(160.0, 120.0, 95.0, 55.0, 48.0, 42.0, 33.0, 29.0, 27.0, 31.0, 22.0, 17.0, ℱ)
    for op in (+, -)
        # Orthotropic tensors on frames that are the same up to a permutation.
        for (A, B) in ((O9(ℬ), O9(ℬp)), (O12(ℬ), O12(ℬp)), (O9(ℬ), O12(ℬp)), (TensOrtho{2}(3.0, 5.0, 7.0, ℬ), TensOrtho{2}(1.0, 2.0, 4.0, ℬp)))
            @test op(A, B) isa TensOrtho && same(A, B, op)
        end
        # A cubic tensor and an orthotropic one on the same cube, permuted.
        C = tens_cubic(10.0, 4.0, 2.0, ℬ)
        @test op(C, O9(ℬp)) isa TensOrtho && same(C, O9(ℬp), op)
        @test op(O12(ℬp), C) isa TensOrtho && same(O12(ℬp), C, op)
        # TI tensors about opposite axes, the shapes with terms odd in the axis
        # included, and about orthogonal axes.
        n, m = Tuple(E[:, 1]), Tuple(E[:, 2])
        for (A, B) in (
                (tens_ti(10.0, 3.0, 2.5, 12.0, 2.0, n), tens_ti(7.0, 1.0, 2.0, 9.0, 3.0, .-n)),
                (TensTI{2, Float64, 3}((1.0, 2.0, 0.5), n), TensTI{2, Float64, 3}((3.0, 1.0, 0.7), .-n)),
                (TensTI{4, Float64, 8}((3.0, 2.0, 0.7, 0.9, 4.0, 1.5, 0.3, 0.2), n), TensTI{4, Float64, 8}((1.0, 5.0, 0.4, 0.6, 2.0, 1.0, 0.5, 0.1), .-n)),
            )
            @test op(A, B) isa TensTI && same(A, B, op)
        end
        for (A, B) in (
                (tens_ti(10.0, 3.0, 2.5, 12.0, 2.0, n), tens_ti(7.0, 1.0, 2.0, 9.0, 3.0, m)),
                (TensTI{4, Float64, 6}((3.0, 2.0, 0.7, 0.9, 4.0, 1.5), n), TensTI{4, Float64, 6}((1.0, 5.0, 0.4, 0.6, 2.0, 1.0), m)),
                (TensTI{2}(1.5, 9.0, n), TensTI{2}(2.0, 4.0, m)),
            )
            @test op(A, B) isa TensOrtho && same(A, B, op)
        end
        # Orthogonal axes leave the shapes odd in the axis no common symmetry.
        A3, B3 = TensTI{2, Float64, 3}((1.0, 2.0, 0.5), n), TensTI{2, Float64, 3}((3.0, 1.0, 0.7), m)
        @test !(op(A3, B3) isa TensTI) && !(op(A3, B3) isa TensOrtho) && same(A3, B3, op)
        # A TI axis along an axis of the cube.
        for T in (tens_ti(10.0, 3.0, 2.5, 12.0, 2.0, m), TensTI{4, Float64, 6}((3.0, 2.0, 0.7, 0.9, 4.0, 1.5), m))
            @test op(C, T) isa TensOrtho && same(C, T, op)
            @test op(T, C) isa TensOrtho && same(T, C, op)
        end
    end
end

@testset "the generic route computed two ways" begin
    # The structured tensors sharing no class sum by their canonical Kelvin–Mandel
    # matrices; the route by their 81 canonical components gives the same tensor.
    by_arrays(op, A, B) = invoke(TensND._generic_binary, Tuple{Any, TensND.AbstractTens{4, 3}, TensND.AbstractTens{4, 3}}, op, A, B)
    As, Bs = _sb_order4(1.0, _sb_frame(), _sb_axis()), _sb_order4(1.0, _sb_frame(), _sb_axis())
    for A in As[2:end], B in Bs[2:end], op in (+, -)
        @test _sb_close(Array(components_canon(TensND._structured_sum(op, A, B))), Array(components_canon(by_arrays(op, A, B))))
    end
    # Symbolic coefficients on numeric frames that share nothing: the same
    # route, and the values it gives once substituted.
    Symbolics = TensND.Symbolics
    TensND.Symbolics.@variables a b c
    ℬ₁, ℬ₂ = _sb_frame(), _sb_frame()
    S = tens_cubic(a, b, c, ℬ₁) + tens_cubic(2a, b, 3c, ℬ₂)
    N = tens_cubic(10.0, 4.0, 2.0, ℬ₁) + tens_cubic(20.0, 4.0, 6.0, ℬ₂)
    values = Dict(a => 10.0, b => 4.0, c => 2.0)
    Sv = [Symbolics.value(Symbolics.substitute(x, values)) for x in Array(components_canon(S))]
    @test _sb_close(Float64.(Sv), Array(components_canon(N)))
end
