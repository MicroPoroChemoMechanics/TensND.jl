# The pointwise submanifold: against the exact one on the sphere, against closed
# forms of the curvatures of a torus, a cylinder and a circle, the
# Laplace–Beltrami operator of the sphere, and derivatives with respect to a
# parameter of the surface.

@testsection "Pointwise submanifolds" begin
    sphere_OM(R) = u -> R * [sin(u[1])cos(u[2]), sin(u[1])sin(u[2]), cos(u[1])]
    R = 2.0
    S = Submanifold(sphere_OM(R), 3)
    u = [0.7, 0.3]
    @test S isa Submanifold{3, Float64} && nderiv(S) == 2

    @testset "the sphere, against the exact submanifold" begin
        θ, ϕ = symbols("θ ϕ", real = true)
        Rs = symbols("R", positive = true)
        E = Submanifold(Tens(Rs * [sin(θ)cos(ϕ), sin(θ)sin(ϕ), cos(θ)]), (θ, ϕ), (), (Rs,))
        at(x) = Float64(x.subs(θ, u[1]).subs(ϕ, u[2]).subs(Rs, R))
        @test christoffel(S, u) ≈ at.(christoffel(E)) atol = 1.0e-14
        @test connection(S, u) ≈ at.(connection(E)) atol = 1.0e-14
        @test Array(curvature(S, u)) ≈ at.(Array(curvature(E))) atol = 1.0e-14
        @test Array(submetric(S, u)) ≈ at.(Array(submetric(E))) atol = 1.0e-14
        @test lame_coefficients(S, u) ≈ at.(collect(lame_coefficients(E)))
        @test canonical_components(normal(S, u)) ≈ [sin(u[1])cos(u[2]), sin(u[1])sin(u[2]), cos(u[1])]
        # The surface gradient of f = cos θ sin ϕ: (∂θ f / R) 𝐞θ + (∂ϕ f / (R sin θ)) 𝐞ϕ.
        f(v) = cos(v[1]) * sin(v[2])
        eθ, eϕ = (canonical_components(unitvec(S, u, i)) for i in 1:2)
        expected = (-sin(u[1]) * sin(u[2]) / R) * eθ + (cos(u[1]) * cos(u[2]) / (R * sin(u[1]))) * eϕ
        @test canonical_components(grad(f, S)(u)) ≈ expected
    end

    @testset "Laplace–Beltrami and divergence on the sphere" begin
        # cos θ is a spherical harmonic of degree 1: Δ cos θ = -2 cos θ / R².
        @test laplacian(v -> cos(v[1]), S)(u) ≈ -2cos(u[1]) / R^2
        # The divergence of the unit normal over the surface: −(principal curvatures) = 2/R.
        @test divergence(v -> [0.0, 0.0, 1.0], S)(u) ≈ 2 / R
    end

    @testset "curvatures of a torus, a cylinder and a circle" begin
        # Torus of radii a > r, (θ, ϕ) around the tube and about the axis:
        # Gaussian curvature cos θ / (r (a + r cos θ)).
        a, r = 3.0, 1.0
        torus = Submanifold(v -> [(a + r * cos(v[1]))cos(v[2]), (a + r * cos(v[1]))sin(v[2]), r * sin(v[1])], 3)
        w = [0.4, 1.1]
        χ = lame_coefficients(torus, w)
        B = Array(curvature(torus, w))
        shape = [B[α, β] / (χ[α] * χ[β]) for α in 1:2, β in 1:2]       # in the orthonormal frame
        @test det(shape) ≈ cos(w[1]) / (r * (a + r * cos(w[1])))
        # A cylinder of radius r: one principal curvature 1/r, one zero.
        cylinder = Submanifold(v -> [r * cos(v[1]), r * sin(v[1]), v[2]], 3)
        C = Array(curvature(cylinder, [0.3, 2.0]))
        χc = lame_coefficients(cylinder, [0.3, 2.0])
        @test sort(abs.([C[1, 1] / χc[1]^2, C[2, 2] / χc[2]^2])) ≈ [0.0, 1 / r] atol = 1.0e-14
        # A circle in the plane, a hypersurface of dimension 1: curvature 1/ρ.
        ρ = 1.5
        circle = Submanifold(v -> [ρ * cos(v[1]), ρ * sin(v[1])], 2)
        @test nderiv(circle) == 1
        κ = Array(curvature(circle, [0.8]))[1, 1] / lame_coefficients(circle, [0.8])[1]^2
        @test abs(κ) ≈ 1 / ρ
    end

    @testset "derivatives with respect to a parameter" begin
        # d/dR (Δ cos θ) = 4 cos θ / R³.
        g(r) = laplacian(v -> cos(v[1]), Submanifold(sphere_OM(r), 3))(u)
        @test ForwardDiff.derivative(g, R) ≈ 4cos(u[1]) / R^3
        @test ForwardDiff.derivative(r -> Array(curvature(Submanifold(sphere_OM(r), 3), u))[1, 1], R) ≈ -1
    end

    @test_throws ArgumentError Submanifold(v -> [v[1]], 1)
end
