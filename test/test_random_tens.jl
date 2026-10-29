# Random tensors. Each class comes back in its own type, inferred, with the
# shape asked for; the positive definite draws are positive definite and the
# shapes without major symmetry refuse them; a fixed axis or frame is kept.
# The orientations are uniform: the moments of a random rotation are those of
# the Haar measure, its angle follows the law of that measure, a random axis is
# uniform on the sphere (Archimedes: its height is uniform on [-1, 1]), and the
# mean of a crystal over random frames is its isotropic projection, an
# independent closed form. Uniform Euler angles fail the same moments, so the
# checks discriminate. The laws of the coefficients are invariant: the
# Kelvin–Mandel moments of a generic tensor are unchanged by a rotation.

using TensND, Test, LinearAlgebra, Random

_rt_rng = Xoshiro(20261029)
_rt_mean(x) = sum(x) / length(x)
_rt_mean(f, x) = sum(f, x) / length(x)

# The largest gap between the empirical distribution of `x` and the continuous
# distribution function `F`: the Kolmogorov–Smirnov statistic.
function _rt_ks(x, F)
    s = sort(x)
    n = length(s)
    return maximum(max(abs(F(v) - (i - 1) / n), abs(F(v) - i / n)) for (i, v) in enumerate(s))
end
# Its 99.9 % quantile for n draws, asymptotically.
_rt_ks_bound(n) = 1.95 / sqrt(n)

@testset "classes, shapes, element types" begin
    for (class, expected) in (
            (TensISO{4, 3}, TensISO{4, 3, Float64, 2}), (TensISO{2, 3}, TensISO{2, 3, Float64, 1}),
            (TensISO{4, 2, Float32}, TensISO{4, 2, Float32, 2}),
            (TensTI{4}, TensTI{4, Float64, 5}), (TensTI{4, Float64, 6}, TensTI{4, Float64, 6}),
            (TensTI{4, Float64, 8}, TensTI{4, Float64, 8}), (TensTI{4, Float32, 5}, TensTI{4, Float32, 5}),
            (TensTI{2}, TensTI{2, Float64, 2}), (TensTI{2, Float64, 3}, TensTI{2, Float64, 3}),
            (TensOrtho{4}, TensOrtho{4, Float64, 9}), (TensOrtho{4, Float64, 12}, TensOrtho{4, Float64, 12}),
            (TensOrtho{2}, TensOrtho{2, Float64, 3}), (TensCubic, TensCubic{Float64}),
            (TensCubic{Float32}, TensCubic{Float32}),
            (Tens{4, 3}, TensND.AbstractTens{4, 3, Float64}), (Tens{4, 2}, TensND.AbstractTens{4, 2, Float64}),
            (Tens{2, 3}, TensND.AbstractTens{2, 3, Float64}), (Tens{2, 3, Float32}, TensND.AbstractTens{2, 3, Float32}),
        )
        t = @inferred rand(_rt_rng, class)
        @test t isa expected
    end
    @test @inferred(rand(_rt_rng, Tens{4, 3}; symmetry = Val(:minor))) isa TensND.AbstractTens{4, 3, Float64}
    # Without symmetry, `Tens` chooses its storage from the values, as everywhere
    # in the package: the union of the two.
    @test @inferred(TensND.AbstractTens{4, 3, Float64}, rand(_rt_rng, Tens{4, 3}; symmetry = Val(:none))) isa TensND.AbstractTens{4, 3, Float64}
    @test @inferred(TensND.AbstractTens{2, 3, Float64}, rand(_rt_rng, Tens{2, 3}; symmetry = Val(:none))) isa TensND.AbstractTens{2, 3, Float64}
    @test rand(TensCubic) isa TensCubic{Float64}      # the default generator
    @test rand(TensCubic; definite = true) isa TensCubic{Float64}
    # The array forms of `rand`, with the default law: the draws of the scalar
    # form, one per entry, from the same stream.
    v = rand(Xoshiro(3), TensTI{4}, 4)
    w = let rng = Xoshiro(3)
        [rand(rng, TensTI{4}) for _ in 1:4]
    end
    @test length(v) == 4 && all(t -> t isa TensTI{4, Float64, 5}, v)
    @test all(components_canon(a) == components_canon(b) for (a, b) in zip(v, w))
    m = Matrix{TensISO{4, 3, Float64}}(undef, 2, 3)
    @test rand!(Xoshiro(4), m) === m && all(t -> t isa TensISO{4, 3, Float64}, m)
    @test !isdefined(TensND, :random_tens)
    @test @inferred(random_rotation(_rt_rng, Float32)) isa TensND.StaticArrays.SMatrix{3, 3, Float32}
    @test @inferred(random_axis(_rt_rng, Float32)) isa NTuple{3, Float32}
    @test @inferred(random_basis(_rt_rng)) isa RotatedBasis{3, Float64}
    @test random_rotation() isa TensND.StaticArrays.SMatrix && random_axis() isa NTuple{3} && random_basis() isa RotatedBasis
    # The same seed, the same tensor.
    @test components_canon(rand(Xoshiro(7), TensOrtho{4})) == components_canon(rand(Xoshiro(7), TensOrtho{4}))
end

@testset "symmetries of the generic draws" begin
    A = Array(components_canon(rand(_rt_rng, Tens{4, 3})))
    @test A ≈ permutedims(A, (2, 1, 3, 4)) && A ≈ permutedims(A, (3, 4, 1, 2))
    A = Array(components_canon(rand(_rt_rng, Tens{4, 3}; symmetry = Val(:minor))))
    @test A ≈ permutedims(A, (2, 1, 3, 4)) && A ≈ permutedims(A, (1, 2, 4, 3)) && !(A ≈ permutedims(A, (3, 4, 1, 2)))
    A = Array(components_canon(rand(_rt_rng, Tens{4, 3}; symmetry = Val(:none))))
    @test !(A ≈ permutedims(A, (2, 1, 3, 4)))
    A = Array(components_canon(rand(_rt_rng, Tens{2, 3})))
    @test A ≈ transpose(A)
end

@testset "positive definite draws, and the refusals" begin
    spectrum(t) = eigvals(Symmetric(Matrix(mandel(t))))
    for class in (TensISO{4, 3}, TensTI{4}, TensOrtho{4}, TensCubic, Tens{4, 3}), _ in 1:50
        @test minimum(spectrum(rand(_rt_rng, class; definite = true))) > 0
    end
    for class in (TensISO{2, 3}, TensTI{2}, TensOrtho{2}, Tens{2, 3}), _ in 1:50
        @test minimum(eigvals(Symmetric(Array(components_canon(rand(_rt_rng, class; definite = true)))))) > 0
    end
    # Without the request, indefinite tensors come too.
    @test any(minimum(spectrum(rand(_rt_rng, TensOrtho{4}))) < 0 for _ in 1:50)
    for class in (TensTI{4, Float64, 6}, TensTI{4, Float64, 8}, TensTI{2, Float64, 3}, TensOrtho{4, Float64, 12})
        @test_throws ArgumentError rand(_rt_rng, class; definite = true)
    end
    @test_throws ArgumentError rand(_rt_rng, Tens{4, 3}; symmetry = Val(:minor), definite = true)
    @test_throws ArgumentError rand(_rt_rng, Tens{4, 3}; symmetry = Val(:none), definite = true)
    @test_throws ArgumentError rand(_rt_rng, Tens{2, 3}; symmetry = Val(:none), definite = true)
    @test_throws ArgumentError rand(_rt_rng, Tens{4, 3}; symmetry = Val(:symmetric))
    @test_throws ArgumentError rand(_rt_rng, Tens{2, 3}; symmetry = Val(:major))
end

@testset "a fixed direction is kept" begin
    n, ℬ = (0.0, 0.6, 0.8), RotatedBasis(0.3, 1.2, 2.0)
    for class in (TensTI{4}, TensTI{4, Float64, 6}, TensTI{4, Float64, 8}, TensTI{2}, TensTI{2, Float64, 3})
        @test axis(rand(_rt_rng, class; axis = n)) == n
        # A frame gives its third vector as the axis, the convention of `frame_from_axis`.
        @test collect(axis(rand(_rt_rng, class; frame = ℬ))) ≈ TensND.vecbasis(ℬ, :cov)[:, 3]
        @test_throws ArgumentError rand(_rt_rng, class; axis = n, frame = ℬ)
    end
    for class in (TensOrtho{4}, TensOrtho{4, Float64, 12}, TensOrtho{2}, TensCubic)
        @test frame(rand(_rt_rng, class; frame = ℬ)) == ℬ
    end
    @test frame(rand(_rt_rng, TensCubic; frame = CanonicalBasis{3, Float64}())) isa CanonicalBasis
end

@testset "uniform orientations" begin
    n = 20_000
    Rs = [random_rotation(_rt_rng) for _ in 1:n]
    @test all(R -> transpose(R) * R ≈ I && det(R) ≈ 1, Rs)
    # Haar moments: every entry of mean 0 and of mean square 1/3; the trace,
    # the character of an irreducible representation, of mean 0 and mean square 1.
    σ = sqrt(1 / (3n))
    @test maximum(abs, _rt_mean(Rs)) < 5σ
    @test maximum(abs, _rt_mean(R -> R .^ 2, Rs) .- 1 / 3) < 5 * sqrt(4 / 45 / n)
    @test abs(_rt_mean(tr, Rs)) < 5 / sqrt(n) && abs(_rt_mean(R -> tr(R)^2, Rs) - 1) < 5 * sqrt(2 / n)
    # The angle ω of the rotation has the distribution function (ω − sin ω)/π.
    ω = [acos(clamp((tr(R) - 1) / 2, -1, 1)) for R in Rs]
    @test _rt_ks(ω, w -> (w - sin(w)) / π) < _rt_ks_bound(n)
    # Uniform Euler angles are not uniform rotations: the axis e₃ of the frame
    # has the height cos θ, of mean square 1/2 instead of 1/3.
    euler = [Matrix(TensND.vecbasis(RotatedBasis(π * rand(_rt_rng), 2π * rand(_rt_rng), 2π * rand(_rt_rng)))) for _ in 1:n]
    @test abs(_rt_mean(R -> R[3, 3]^2, euler) - 1 / 3) > 0.1
    @test abs(_rt_mean(R -> R[3, 3]^2, Rs) - 1 / 3) < 0.01
    # A random axis: uniform on the sphere, its height uniform on [-1, 1].
    ns = [collect(random_axis(_rt_rng)) for _ in 1:n]
    @test all(v -> norm(v) ≈ 1, ns)
    @test norm(_rt_mean(ns)) < 5 * sqrt(1 / (3n)) && norm(_rt_mean(v -> v * transpose(v), ns) - I / 3) < 0.02
    @test _rt_ks([v[3] for v in ns], z -> (z + 1) / 2) < _rt_ks_bound(n)
    @test _rt_ks([atan(v[2], v[1]) for v in ns], φ -> (φ + π) / (2π)) < _rt_ks_bound(n)
    # The frame of a random basis is a random rotation.
    ℬs = [Matrix(TensND.vecbasis(random_basis(_rt_rng))) for _ in 1:n]
    @test abs(_rt_mean(R -> R[3, 3]^2, ℬs) - 1 / 3) < 0.01 && abs(_rt_mean(R -> tr(R)^2, ℬs) - 1) < 5 * sqrt(2 / n)
end

@testset "the mean over random orientations is the isotropic projection" begin
    # A crystal of fixed constants on random frames, and a transversely
    # isotropic tensor about random axes: their mean is the isotropic part of
    # either, computed in closed form by `isotropify`.
    n = 4000
    C = tens_cubic(10.0, 4.0, 2.0, CanonicalBasis{3, Float64}())
    T = tens_ti(10.0, 3.0, 2.5, 12.0, 2.0, (0.0, 0.0, 1.0))
    for (draw, ref) in (
            (() -> tens_cubic(10.0, 4.0, 2.0, random_basis(_rt_rng)), C),
            (() -> tens_ti(10.0, 3.0, 2.5, 12.0, 2.0, random_axis(_rt_rng)), T),
        )
        m = sum(_ -> Array(components_canon(draw())), 1:n) / n
        iso = Array(components_canon(isotropify(ref)))
        @test maximum(abs, m - iso) < 0.05 * maximum(abs, iso)
    end
end

@testset "the laws of the coefficients are invariant under rotations" begin
    # The Kelvin–Mandel matrix of a generic tensor is drawn from the Gaussian
    # orthogonal ensemble: diagonal entries of mean square 1, off-diagonal 1/2,
    # whatever the basis it is read in.
    n = 20_000
    ℛ = RotatedBasis(0.7, 2.1, -0.4)
    canonical = [Matrix(mandel(rand(_rt_rng, Tens{4, 3}))) for _ in 1:n]
    rotated = [Matrix(mandel(change_tens(rand(_rt_rng, Tens{4, 3}), ℛ))) for _ in 1:n]
    for sample in (canonical, rotated)
        s = _rt_mean(K -> K .^ 2, sample)
        tol = 5 * sqrt(2 / length(sample))
        @test maximum(abs, diag(s) .- 1) < tol
        @test maximum(abs, (s - Diagonal(diag(s))) .- (1 .- I(6)) / 2) < tol
    end
end
