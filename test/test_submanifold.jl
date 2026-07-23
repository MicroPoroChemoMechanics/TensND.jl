@testsection "Submanifolds" begin

    # `SubManifoldSym` builds the induced geometry of a hypersurface given by a
    # symbolic parametrization `OM(coords)`: unit normal, first fundamental
    # form (`submetric`), second fundamental form (`curvature`) and the
    # connection coefficients (`Riemann`).
    #
    # Assertions are stated as geometric invariants (unit normal, orthogonality
    # to the tangent plane, Gaussian curvature) so they hold regardless of the
    # `Abs(sin θ)` factors SymPy keeps for want of a sign assumption.

    # Substitute a concrete interior point to discharge those `Abs`.
    _at(x, θ, ϕ, R; θ₀ = Sym(1) / 3, ϕ₀ = Sym(2) / 5, R₀ = 2) =
        tsimplify(tsubs(x, θ => θ₀, ϕ => ϕ₀, R => R₀))

    @testsection "Sphere — normal, fundamental forms, curvature" begin
        R = symbols("R", positive = true)
        θ, ϕ = symbols("θ ϕ", real = true)
        OM = R * Tens(Vec{3}([sin(θ)cos(ϕ), sin(θ)sin(ϕ), cos(θ)]))
        SM = SubManifoldSym(OM, (θ, ϕ))

        @test SM isa SubManifoldSym{3}

        n = components_canon(normal(SM))
        # Unit normal.
        @test _at(sum(n .^ 2), θ, ϕ, R) == 1

        # Orthogonal to both tangent vectors ∂OM/∂θ and ∂OM/∂ϕ.
        OMc = components_canon(OM)
        for c in (θ, ϕ)
            tang = [tdiff(OMc[i], c) for i in 1:3]
            @test _at(sum(n .* tang), θ, ϕ, R) == 0
        end

        # First fundamental form of the sphere: diag(R², R² sin²θ), padded
        # with a zero in the normal direction.
        a = Array(components(submetric(SM)))
        @test tsimplify(a[1, 1] - R^2) == 0
        @test tsimplify(a[2, 2] - R^2 * sin(θ)^2) == 0
        @test tsimplify(a[1, 2]) == 0
        @test tsimplify(a[3, 3]) == 0

        # Gaussian curvature K = det(b)/det(a) = 1/R² — the `Abs` cancels since
        # it enters squared.
        b = Array(components(curvature(SM)))
        K = tsimplify(
            (b[1, 1] * b[2, 2] - b[1, 2] * b[2, 1]) /
                (a[1, 1] * a[2, 2] - a[1, 2] * a[2, 1])
        )
        @test tsimplify(K - 1 / R^2) == 0

        # Mean curvature H = -1/R with the outward normal convention used here
        # (b is negative definite for the outward normal).
        @test _at(b[1, 1], θ, ϕ, R) < 0
    end

    @testsection "Sphere — connection coefficients" begin
        R = symbols("R", positive = true)
        θ, ϕ = symbols("θ ϕ", real = true)
        OM = R * Tens(Vec{3}([sin(θ)cos(ϕ), sin(θ)sin(ϕ), cos(θ)]))
        SM = SubManifoldSym(OM, (θ, ϕ))

        Γ = Riemann(SM)
        @test size(Γ) == (2, 2, 2)

        # Christoffel symbols of the round sphere in (θ, ϕ):
        #   Γ^θ_ϕϕ = -sin θ cos θ,   Γ^ϕ_θϕ = Γ^ϕ_ϕθ = cot θ,  others 0.
        @test tsimplify(Γ[2, 2, 1] + sin(θ)cos(θ)) == 0
        @test tsimplify(Γ[1, 2, 2] - cos(θ) / sin(θ)) == 0
        @test tsimplify(Γ[2, 1, 2] - cos(θ) / sin(θ)) == 0
        @test tsimplify(Γ[1, 1, 1]) == 0
        @test tsimplify(Γ[1, 1, 2]) == 0
        @test tsimplify(Γ[2, 2, 2]) == 0
    end

    @testsection "Cylinder — vanishing Gaussian curvature" begin
        R = symbols("R", positive = true)
        ϕ, z = symbols("ϕ z", real = true)
        OM = Tens(Vec{3}([R * cos(ϕ), R * sin(ϕ), z]))
        SM = SubManifoldSym(OM, (ϕ, z))

        a = Array(components(submetric(SM)))
        b = Array(components(curvature(SM)))

        # First fundamental form diag(R², 1).
        @test tsimplify(a[1, 1] - R^2) == 0
        @test tsimplify(a[2, 2] - 1) == 0
        @test tsimplify(a[1, 2]) == 0

        # A cylinder is developable: det(b) = 0, hence K = 0.
        @test tsimplify(b[1, 1] * b[2, 2] - b[1, 2] * b[2, 1]) == 0
        # …but it is curved: b ≠ 0 along the hoop direction.
        @test tsimplify(b[1, 1]) != 0
        @test tsimplify(b[2, 2]) == 0

        # Unit normal, radial.
        n = components_canon(normal(SM))
        @test tsimplify(sum(n .^ 2) - 1) == 0
    end

    @testsection "Plane — flat submanifold" begin
        x, y = symbols("x y", real = true)
        OM = Tens(Vec{3}([x, y, Sym(0)]))
        SM = SubManifoldSym(OM, (x, y))

        a = Array(components(submetric(SM)))
        b = Array(components(curvature(SM)))

        # Euclidean induced metric.
        @test tsimplify(a[1, 1] - 1) == 0
        @test tsimplify(a[2, 2] - 1) == 0
        @test tsimplify(a[1, 2]) == 0

        # A plane has no curvature at all.
        @test all(tsimplify.(b) .== 0)

        # …and no connection coefficients.
        @test all(tsimplify.(Riemann(SM)) .== 0)

        # Normal is the constant e₃.
        n = components_canon(normal(SM))
        @test tsimplify(n[1]) == 0
        @test tsimplify(n[2]) == 0
        @test tsimplify(abs(n[3]) - 1) == 0
    end

    @testsection "Differential operators on a submanifold" begin
        # On the plane the submanifold operators must collapse to the plain
        # Cartesian ones — the cheapest non-trivial check of ∂/GRAD/DIV/LAPLACE.
        x, y = symbols("x y", real = true)
        OM = Tens(Vec{3}([x, y, Sym(0)]))
        SM = SubManifoldSym(OM, (x, y))

        f = x^2 * y + 3y
        @test tsimplify(∂(f, 1, SM) - 2x * y) == 0
        @test tsimplify(∂(f, 2, SM) - (x^2 + 3)) == 0

        # Same through the coordinate-symbol form of ∂.
        @test tsimplify(∂(f, x, SM) - 2x * y) == 0
        @test tsimplify(∂(f, y, SM) - (x^2 + 3)) == 0
        # A symbol that is not a coordinate gives zero.
        @test tsimplify(∂(f, symbols("w", real = true), SM)) == 0

        # GRAD of a scalar, expressed in the canonical basis.
        g = components_canon(GRAD(f, SM))
        @test tsimplify(g[1] - 2x * y) == 0
        @test tsimplify(g[2] - (x^2 + 3)) == 0

        # LAPLACE f = ∂²f/∂x² + ∂²f/∂y² = 2y.
        @test tsimplify(LAPLACE(f, SM) - 2y) == 0

        # HESS is the gradient of the gradient; its trace is the Laplacian.
        H = components_canon(HESS(f, SM))
        @test tsimplify(H[1, 1] + H[2, 2] - 2y) == 0

        # DIV of a vector field on the plane.
        V = Tens(Vec{3}([x^2, y^3, Sym(0)]))
        @test tsimplify(DIV(V, SM) - (2x + 3y^2)) == 0

        # SYMGRAD of a vector field is the symmetric part of GRAD.
        S = components_canon(SYMGRAD(V, SM))
        @test tsimplify(S[1, 2] - S[2, 1]) == 0
    end

end
