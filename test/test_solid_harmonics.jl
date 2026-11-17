# Solid harmonics, vector harmonics, solutions of the Lamé equation and their
# translation. References independent of the Cartesian recurrences under test:
# closed forms of the low degrees, the meridian harmonics of the spherical chart,
# the differential relations that define the families (gradient, divergence,
# curl) and the Lamé equation through ForwardDiff, and the direct evaluation of
# a translated harmonic or solution.

import Symbolics

# Real and imaginary parts differentiated apart: the fields are complex.
_so_jac(f, x) = ForwardDiff.jacobian(v -> real.(f(v)), x) + im * ForwardDiff.jacobian(v -> imag.(f(v)), x)
_so_div(f, x) = tr(_so_jac(f, x))
function _so_curl(f, x)
    J = _so_jac(f, x)
    return [J[3, 2] - J[2, 3], J[1, 3] - J[3, 1], J[2, 1] - J[1, 2]]
end
# The Lamé operator 2(1-ν)/(1-2ν) ∇(∇·u) - ∇×∇×u.
function _so_lame(f, x, ν)
    graddiv = ForwardDiff.gradient(v -> real(_so_div(f, v)), x) + im * ForwardDiff.gradient(v -> imag(_so_div(f, v)), x)
    return 2 * (1 - ν) / (1 - 2ν) * graddiv - _so_curl(v -> _so_curl(f, v), x)
end

@testsection "Solid harmonics" begin
    x = [0.3, -0.7, 0.45]

    @testsection "Closed forms of the low degrees" begin
        w, z, r = complex(x[1], x[2]), x[3], norm(x)
        y = regular_solid_harmonics(x, 3)
        Y = irregular_solid_harmonics(x, 3)
        @test y[0, 0] ≈ 1 && y[1, 0] ≈ z && y[1, 1] ≈ w / 2 && y[1, -1] ≈ -conj(w) / 2
        @test y[2, 0] ≈ (3z^2 - r^2) / 4 && y[2, 1] ≈ z * w / 2 && y[2, 2] ≈ w^2 / 8 && y[2, -2] ≈ conj(w)^2 / 8
        @test Y[0, 0] ≈ 1 / r && Y[1, 0] ≈ z / r^3 && Y[1, 1] ≈ w / r^3 && Y[1, -1] ≈ -conj(w) / r^3
        @test Y[2, 0] ≈ (3z^2 - r^2) / r^5 && Y[2, 1] ≈ 3z * w / r^5
        # an order beyond the degree reads as zero, a degree beyond the table is refused
        @test y[2, 3] == 0 && Y[1, -2] == 0
        @test_throws BoundsError y[4, 0]
        @test_throws BoundsError y[-1, 0]
        @test regular_solid_harmonics(zeros(3), 2)[0, 0] == 1
        # a table from its values
        h = SolidHarmonicTable((n, m) -> y[n, m], 3)
        @test h.values == y.values && h[2, -1] == y[2, -1]
        @test y.values[solid_harmonic_index(2, 1)] == y[2, 1] && solid_harmonic_index(3, 3) == length(y.values)
        @test_throws DomainError solid_harmonic_index(2, 3)
    end

    @testsection "The meridian harmonics of the spherical chart" begin
        ξ = cartesian_to_coords(SphericalChart(), x)       # (θ, ϕ, r)
        N = 8
        y = regular_solid_harmonics(x, N)
        Y = irregular_solid_harmonics(x, N)
        for m in 0:N
            I = interior_harmonics(SphericalChart(), N, m, ξ)
            E = exterior_harmonics(SphericalChart(), N, m, ξ)
            for n in m:N
                @test y[n, m] ≈ (-1)^m * I[n + 1] * cis(m * ξ[2]) / factorial(n + m) rtol = 1.0e-12
                @test Y[n, m] ≈ (-1)^m * factorial(n - m) * E[n + 1] * cis(m * ξ[2]) rtol = 1.0e-12
            end
        end
        # Y_n^m = (n-m)! (n+m)! y_n^m / r^{2n+1}, and the norm
        for n in 0:N, m in (-n):n
            @test Y[n, m] ≈ solid_harmonic_norm(n, m)^2 * y[n, m] / norm(x)^(2n + 1) rtol = 1.0e-12
        end
        @test solid_harmonic_norm(2, 1) ≈ sqrt(6) && solid_harmonic_norm(3, -3) ≈ sqrt(720)
    end

    @testsection "Derivatives and harmonicity" begin
        N = 5
        y = regular_solid_harmonics(x, N + 1)
        Y = irregular_solid_harmonics(x, N + 1)
        dy = ForwardDiff.derivative(h -> real.(regular_solid_harmonics(x .+ h .* [0, 0, 1], N).values), 0.0)
        dY = ForwardDiff.derivative(h -> imag.(irregular_solid_harmonics(x .+ h .* [0, 0, 1], N).values), 0.0)
        k = 0
        for t in 0:N, s in 0:t
            k += 1
            # ∂₃ y_t^s = y_{t-1}^s and ∂₃ Y_t^s = -Y_{t+1}^s (Hobson)
            @test dy[k] ≈ (t > 0 ? real(y[t - 1, s]) : 0.0) atol = 1.0e-13
            @test dY[k] ≈ -imag(Y[t + 1, s]) rtol = 1.0e-12 atol = 1.0e-13
        end
        lap(f) = tr(ForwardDiff.hessian(f, x))
        @test abs(lap(v -> real(regular_solid_harmonics(v, 4)[4, 2]))) < 1.0e-12
        @test abs(lap(v -> imag(irregular_solid_harmonics(v, 4)[4, 3]))) < 1.0e-9
    end

    @testsection "Translation of a decaying harmonic" begin
        # Y_k^l(r + R) = Σ (-1)^{t+s} Y_{k+t}^{l-s}(R) y_t^s(r), ‖r‖ < ‖R‖
        R, r = [1.1, -0.6, 0.9], [0.15, 0.2, -0.1]
        N = 30
        YR = irregular_solid_harmonics(R, N + 3)
        yr = regular_solid_harmonics(r, N)
        direct = irregular_solid_harmonics(R .+ r, 3)
        for (k, l) in ((1, 0), (2, -1), (3, 3))
            @test sum((-1)^(t + s) * YR[k + t, l - s] * yr[t, s] for t in 0:N for s in (-t):t) ≈ direct[k, l] rtol = 1.0e-12
        end
    end

    @testsection "Types, symbols and refusals" begin
        @test @inferred(regular_solid_harmonics([1.0, 2.0, 3.0], 3)) isa SolidHarmonicTable{ComplexF64}
        @test @inferred(irregular_solid_harmonics([1, 2, 3], 3)) isa SolidHarmonicTable{ComplexF64}
        @test regular_solid_harmonics(Float32[1, 2, 3], 2) isa SolidHarmonicTable{ComplexF32}
        Symbolics.@variables a b c
        ys = regular_solid_harmonics([a, b, c], 3)
        at = Dict(a => 0.3, b => -0.7, c => 0.45)
        val(e) = complex(
            Float64(Symbolics.value(Symbolics.substitute(real(e), at; fold = Val(true)))),
            Float64(Symbolics.value(Symbolics.substitute(imag(e), at; fold = Val(true))))
        )
        @test val(ys[3, 1]) ≈ regular_solid_harmonics(x, 3)[3, 1] rtol = 1.0e-12
        @test_throws DimensionMismatch regular_solid_harmonics([1.0, 2.0], 2)
        @test_throws DomainError regular_solid_harmonics([1.0, 2.0, 3.0], -1)
        @test_throws DomainError irregular_solid_harmonics(zeros(3), 2)
    end
end

@testsection "Vector harmonics and Lamé solutions" begin
    x, ν = [0.3, -0.7, 0.45], 0.3
    Yt(v, n) = irregular_solid_harmonics(v, n)
    yt(v, n) = regular_solid_harmonics(v, n)
    grad(f) = ForwardDiff.gradient(v -> real(f(v)), x) + im * ForwardDiff.gradient(v -> imag(f(v)), x)

    @testsection "The relations that define the families" begin
        for (t, s) in ((1, 0), (2, 1), (3, -2), (4, 4))
            Y1(v) = irregular_vector_harmonic(1, t, s, Yt(v, t + 1))
            Y2(v) = irregular_vector_harmonic(2, t, s, Yt(v, t + 1))
            Y3(v) = irregular_vector_harmonic(3, t, s, Yt(v, t + 1))
            y1(v) = regular_vector_harmonic(1, t, s, yt(v, t + 1))
            y2(v) = regular_vector_harmonic(2, t, s, yt(v, t + 1))
            y3(v) = regular_vector_harmonic(3, t, s, yt(v, t + 1))
            @test Y1(x) ≈ grad(v -> Yt(v, t)[t, s]) rtol = 1.0e-11
            @test y1(x) ≈ grad(v -> yt(v, t)[t, s]) rtol = 1.0e-11
            @test Y2(x) ≈ cross(Y1(x), x) / t rtol = 1.0e-12          # rotation fields
            @test y2(x) ≈ -cross(y1(x), x) / (t + 1) rtol = 1.0e-12
            @test abs(_so_div(Y2, x)) < 1.0e-9 && abs(_so_div(y2, x)) < 1.0e-11
            @test _so_div(Y3, x) ≈ -Yt(x, t)[t, s] rtol = 1.0e-10
            @test _so_div(y3, x) ≈ -yt(x, t)[t, s] rtol = 1.0e-10
            @test _so_curl(Y2, x) ≈ -Y1(x) rtol = 1.0e-10
            @test _so_curl(Y3, x) ≈ Y2(x) rtol = 1.0e-10
            @test _so_curl(y3, x) ≈ y2(x) rtol = 1.0e-10
            # x·Y⁽³⁾ = r² Y_t^s/(2t-1) and x·y⁽³⁾ = -r² y_t^s/(2t+3)
            @test sum(x .* Y3(x)) ≈ sum(abs2, x) * Yt(x, t)[t, s] / (2t - 1) rtol = 1.0e-12
            @test sum(x .* y3(x)) ≈ -sum(abs2, x) * yt(x, t)[t, s] / (2t + 3) rtol = 1.0e-12
        end
        @test_throws DomainError regular_vector_harmonic(4, 1, 0, yt(x, 2))
        @test_throws DomainError irregular_vector_harmonic(0, 1, 0, Yt(x, 2))
    end

    @testsection "The Lamé equation and the low degrees" begin
        for (t, s) in ((1, 1), (2, 0), (3, -1)), i in 1:3
            U(v) = irregular_lame_solution(i, t, s, v, ν)
            u(v) = regular_lame_solution(i, t, s, v, ν)
            @test norm(_so_lame(U, x, ν)) < 1.0e-7 * norm(U(x)) / norm(x)^2
            @test norm(_so_lame(u, x, ν)) < 1.0e-9
        end
        # the dilatation, the translations, a rotation, a uniform deviatoric strain
        @test regular_lame_solution(3, 0, 0, x, ν) ≈ (4ν - 2) / 3 * x
        @test regular_lame_solution(1, 1, 0, x, ν) ≈ [0, 0, 1]
        @test regular_lame_solution(1, 1, 1, x, ν) ≈ [1, im, 0] / 2
        @test regular_lame_solution(2, 1, 0, x, ν) ≈ [x[2], -x[1], 0] / 2
        @test regular_lame_solution(1, 2, 0, x, ν) ≈ [-x[1] / 2, -x[2] / 2, x[3]]
        @test irregular_lame_solution(1, 0, 0, x, ν) ≈ -x / norm(x)^3     # a center of dilatation
        # divergence and curl of the third family (kushch2013, Eq. (4.35))
        for (t, s) in ((2, 1), (3, 0))
            @test _so_div(v -> irregular_lame_solution(3, t, s, v, ν), x) ≈ 2 * (2ν - 1) * Yt(x, t)[t, s] rtol = 1.0e-9
            @test _so_curl(v -> regular_lame_solution(3, t, s, v, ν), x) ≈ 4 * (1 - ν) * regular_lame_solution(2, t, s, x, ν) rtol = 1.0e-9
        end
        @test @inferred(regular_lame_solution(3, 2, 1, x, ν)) isa AbstractVector{ComplexF64}
        @test @inferred(irregular_lame_solution(2, 2, 1, x, ν)) isa AbstractVector{ComplexF64}
        # the Poisson ratio may carry a derivative: dU⁽³⁾/dν = -4 Y⁽³⁾
        @test ForwardDiff.derivative(n -> real(irregular_lame_solution(3, 2, 1, x, n)), ν) ≈
            -4 * real(irregular_vector_harmonic(3, 2, 1, Yt(x, 3))) rtol = 1.0e-12
    end

    @testsection "Translation of the Lamé solutions" begin
        R, r = [1.1, -0.6, 0.9], [0.12, 0.2, -0.1]
        kmax = 26
        Y = irregular_solid_harmonics(R, kmax + 4)
        for (i, t, s) in ((1, 0, 0), (1, 2, 1), (2, 1, 0), (2, 3, 2), (3, 1, 0), (3, 1, 1), (3, 2, 0), (3, 2, 2), (3, 3, -1))
            series = zeros(ComplexF64, 3)
            for k in 0:kmax, l in (-k):k, j in 1:3
                a, b = lame_solution_translation(i, j, t, k, s, l, Y, R)
                series += (-1)^(k + l) * (a + 4 * (1 - ν) * b) * regular_lame_solution(j, k, l, r, ν)
            end
            @test series ≈ irregular_lame_solution(i, t, s, R + r, ν) rtol = 1.0e-12
        end
        # the families above the one translated do not appear
        @test lame_solution_translation(1, 2, 2, 3, 0, 1, Y, R) == (0, 0)
        @test lame_solution_translation(2, 3, 2, 3, 0, 1, Y, R) == (0, 0)
        # the derivative along the position of the second center, against the
        # derivative of the direct evaluation of the translated solution
        f(h) = (
            Rh = R + h * [0, 0, 1]; Yh = irregular_solid_harmonics(Rh, 34);
            sum((-1)^(k + l) * sum(lame_solution_translation(3, j, 2, k, 1, l, Yh, Rh) .* (1, 4 * (1 - ν))) * regular_lame_solution(j, k, l, r, ν) for k in 0:30 for l in (-k):k for j in 1:3)
        )
        @test ForwardDiff.derivative(h -> real(f(h)), 0.0) ≈
            ForwardDiff.derivative(h -> real(irregular_lame_solution(3, 2, 1, R + h * [0, 0, 1] + r, ν)), 0.0) rtol = 1.0e-9
    end
end
