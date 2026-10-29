# The sum and the difference of two structured tensors, for every pair of
# classes, on frames and axes drawn independently: the canonical components of
# the result are the sums and differences of those of the operands, whatever
# route the pair takes (its own class when the references agree, the generic
# one otherwise), and no pair refuses. The same with dual coefficients, the
# derivative of the result being the operand the coefficients scale. Pairs on
# shared references keep their class.

using TensND, Test, LinearAlgebra, ForwardDiff, Random

_sb_rng = MersenneTwister(20261029)
_sb_frame() = RotatedBasis(π * rand(_sb_rng), 2π * rand(_sb_rng), 2π * rand(_sb_rng))
_sb_axis() = Tuple(normalize(randn(_sb_rng, 3)))

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
