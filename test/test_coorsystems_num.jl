@testsection "Numerical coordinate systems (AD)" begin

    tol = 1.0e-10

    # ─────────────────────────────────────────────────────────────
    # Cartesian 3D
    # ─────────────────────────────────────────────────────────────
    @testsection "Cartesian (numeric)" begin
        CS = coorsys_cartesian_num()
        x₀ = [1.0, 2.0, 3.0]

        # Christoffel symbols vanish in Cartesian
        Γ = CS.Γ_func(x₀)
        @test all(abs.(Γ) .< tol)

        # normalized_basis returns a CanonicalBasis (identity rotation)
        @test normalized_basis(CS, x₀) isa CanonicalBasis

        # natvec and unitvec return AbstractTens
        @test natvec(CS, x₀, 1, :cov) isa AbstractTens
        @test natvec(CS, x₀, 2, :cont) isa AbstractTens
        @test unitvec(CS, x₀, 3) isa AbstractTens

        # Gradient of scalar: ∇f = [∂f/∂x, ∂f/∂y, ∂f/∂z]
        f = x -> x[1]^2 + x[2] * x[3]
        g = GRAD(f, CS)(x₀)
        @test g isa AbstractTens
        @test Array(g) ≈ [2 * x₀[1], x₀[3], x₀[2]]

        # Laplacian of scalar: ∇²f = 2
        @test LAPLACE(f, CS)(x₀) ≈ 2.0

        # Divergence of vector field v = [x, y, z] → div = 3
        v = x -> x
        @test DIV(v, CS)(x₀) ≈ 3.0

        # SYMGRAD of linear vector field v=[x,y,z]: sym(I) = I
        sg = SYMGRAD(v, CS)(x₀)
        @test sg isa AbstractTens
        @test Array(sg) ≈ Matrix(I, 3, 3)

        # HESS of quadratic scalar
        fq = x -> x[1]^2 + 3 * x[2]^2
        H = HESS(fq, CS)(x₀)
        @test H isa AbstractTens
        @test Array(H)[1, 1] ≈ 2.0  atol = 1.0e-8
        @test Array(H)[2, 2] ≈ 6.0  atol = 1.0e-8
        @test Array(H)[3, 3] ≈ 0.0  atol = 1.0e-8
    end

    # ─────────────────────────────────────────────────────────────
    # Polar 2D
    # ─────────────────────────────────────────────────────────────
    @testsection "Polar (numeric)" begin
        CS = coorsys_polar_num()
        r₀, θ₀ = 2.0, π / 4
        x₀ = [r₀, θ₀]

        # Christoffel symbols for polar: Γ²₁₂ = Γ²₂₁ = 1/r, Γ¹₂₂ = -r
        # Convention: Γ[i,j,k] = Γᵏᵢⱼ
        Γ = CS.Γ_func(x₀)
        @test Γ[2, 2, 1] ≈ -r₀   atol = tol   # Γ¹₂₂ = -r
        @test Γ[1, 2, 2] ≈ 1 / r₀  atol = tol   # Γ²₁₂ = 1/r
        @test Γ[2, 1, 2] ≈ 1 / r₀  atol = tol   # Γ²₂₁ = 1/r
        @test abs(Γ[1, 1, 1]) < tol
        @test abs(Γ[1, 1, 2]) < tol

        # normalized_basis at non-trivial angle returns a RotatedBasis
        ℬ = normalized_basis(CS, x₀)
        @test ℬ isa RotatedBasis
        @test size(ℬ) == (2, 2)

        # natvec: contravariant natural vector aʳ = eʳ/χʳ = eʳ (χ₁=1)
        a1_cont = natvec(CS, x₀, 1, :cont)
        @test a1_cont isa AbstractTens
        @test abs(Array(a1_cont)[1] - 1.0) < tol   # χ₁=1 → 1/χ₁=1

        # Laplacian of scalar f(r,θ) = r²: ∇²(r²) = 4
        f = x -> x[1]^2
        @test LAPLACE(f, CS)(x₀) ≈ 4.0 atol = 1.0e-8

        # Laplacian of f = r²cos²θ: ∇²f = 2 (standard result in polar)
        f2 = x -> x[1]^2 * cos(x[2])^2
        lap = LAPLACE(f2, CS)(x₀)
        # ∇²(r²) = 4, ∇²(r²cos2θ) = 0  → ∇²f2 = 2
        @test lap ≈ 2.0 atol = 1.0e-8

        # GRAD of scalar: ∇(r²) = [2r, 0] in (𝐞ʳ, 𝐞ᶿ) components
        g = GRAD(x -> x[1]^2, CS)(x₀)
        @test g isa AbstractTens
        @test Array(g)[1] ≈ 2 * r₀  atol = 1.0e-8
        @test abs(Array(g)[2]) < 1.0e-8
    end

    # ─────────────────────────────────────────────────────────────
    # Spherical 3D — consistency with symbolic results
    # ─────────────────────────────────────────────────────────────
    @testsection "Spherical (numeric)" begin
        CS = coorsys_spherical_num()
        θ₀, ϕ₀, r₀ = π / 3, π / 4, 2.0
        x₀ = [θ₀, ϕ₀, r₀]

        # Christoffel symbols (non-zero for spherical, coords = θ,ϕ,r)
        Γ = CS.Γ_func(x₀)

        # Known non-zero Christoffel symbols for spherical (θ,ϕ,r) ordering,
        # Lamé = (r, r·sinθ, 1):
        @test Γ[1, 1, 3] ≈ -r₀                    atol = 1.0e-8
        @test Γ[2, 2, 3] ≈ -r₀ * sin(θ₀)^2         atol = 1.0e-8
        @test Γ[2, 2, 1] ≈ -sin(θ₀) * cos(θ₀)      atol = 1.0e-8
        @test Γ[1, 3, 1] ≈ 1 / r₀                   atol = 1.0e-8
        @test Γ[3, 1, 1] ≈ 1 / r₀                   atol = 1.0e-8
        @test Γ[2, 3, 2] ≈ 1 / r₀                   atol = 1.0e-8
        @test Γ[3, 2, 2] ≈ 1 / r₀                   atol = 1.0e-8
        @test Γ[1, 2, 2] ≈ cos(θ₀) / sin(θ₀)       atol = 1.0e-8
        @test Γ[2, 1, 2] ≈ cos(θ₀) / sin(θ₀)       atol = 1.0e-8

        # normalized_basis returns a RotatedBasis
        @test normalized_basis(CS, x₀) isa RotatedBasis

        # natvec: covariant aᵣ = χᵣ * 𝐞ʳ  (χ₃ = 1 for r-coord)
        ar_cov = natvec(CS, x₀, 3, :cov)
        @test ar_cov isa AbstractTens
        @test Array(ar_cov)[3] ≈ 1.0  atol = tol   # 3rd component = χ₃ = 1

        # natvec: contravariant aᶿ  (χ₁ = r)
        at_cont = natvec(CS, x₀, 1, :cont)
        @test at_cont isa AbstractTens
        @test Array(at_cont)[1] ≈ 1 / r₀  atol = tol

        # Lamé coefficients accessor
        χ = Lame(CS, x₀)
        @test χ[1] ≈ r₀         atol = tol
        @test χ[2] ≈ r₀ * sin(θ₀) atol = tol
        @test χ[3] ≈ 1.0        atol = tol

        # Laplacian of scalar f = r²: ∇²(r²) = 6  (in 3D)
        f = x -> x[3]^2
        @test LAPLACE(f, CS)(x₀) ≈ 6.0 atol = 1.0e-7

        # Laplacian of f = r: ∇²(r) = 2/r
        fr = x -> x[3]
        @test LAPLACE(fr, CS)(x₀) ≈ 2 / r₀ atol = 1.0e-7

        # GRAD of scalar f = r: ∇r = 𝐞ʳ → only 3rd component (r-component) = 1
        grad_r = GRAD(fr, CS)(x₀)
        @test grad_r isa AbstractTens
        @test abs(Array(grad_r)[1]) < 1.0e-7   # θ component = 0
        @test abs(Array(grad_r)[2]) < 1.0e-7   # ϕ component = 0
        @test Array(grad_r)[3] ≈ 1.0  atol = 1.0e-7  # r component = 1

        # HESS of f = r²: H = diag(0, 0, 2) in normalized spherical basis
        # Actually ∇∇(r²) = 2𝐞ʳ⊗𝐞ʳ + (2/r)𝐞ᶿ⊗𝐞ᶿ + (2/r)𝐞ᵠ⊗𝐞ᵠ … let's just check symmetry
        H = HESS(f, CS)(x₀)
        @test H isa AbstractTens
        Harr = Array(H)
        @test Harr ≈ Harr'  atol = 1.0e-7   # symmetry

        # Divergence of radial vector field v = (0,0,vᵣ(r)) with vᵣ = r²
        # div(r²𝐞ʳ) = (1/r²) d(r²·r²)/dr = 4r
        vr_field = x -> [zero(x[1]), zero(x[1]), x[3]^2]
        @test DIV(vr_field, CS)(x₀) ≈ 4 * r₀ atol = 1.0e-7

        # SYMGRAD of radial displacement field u = [0, 0, u_r(r)]
        # ε_rr = du_r/dr = 2r at r₀=2 → ε_rr = 4 (no: u_r = r^2 → du_r/dr = 2r)
        u_func = x -> [zero(x[1]), zero(x[1]), x[3]^2]
        ε = SYMGRAD(u_func, CS)(x₀)
        @test ε isa AbstractTens
        @test Array(ε)[3, 3] ≈ 2 * r₀  atol = 1.0e-7  # ε_rr = du_r/dr = 2r
        @test Array(ε)[1, 2] ≈ 0.0   atol = 1.0e-7  # no off-diagonal

        # DIV of matrix field: stress equilibrium check (Lamé problem)
        a_lame = 1.0; b_lame = 3.0; p_i = 1.0
        κ_val = 2.0; μ_val = 1.0
        λ_val = κ_val - 2μ_val / 3
        A_lame = p_i * a_lame^3 / (3κ_val * (b_lame^3 - a_lame^3))
        B_lame = p_i * a_lame^3 * b_lame^3 / (4μ_val * (b_lame^3 - a_lame^3))
        u_lame = x -> [zero(x[1]), zero(x[1]), A_lame * x[3] + B_lame / x[3]^2]
        ε_func = x -> Array(SYMGRAD(u_lame, CS)(x))
        σ_func = x -> begin
            ε = ε_func(x)
            tr_ε = sum(ε[i, i] for i in 1:3)
            [λ_val * tr_ε * (i == j ? 1.0 : 0.0) + 2μ_val * ε[i, j] for i in 1:3, j in 1:3]
        end
        div_σ = DIV(σ_func, CS)([π / 3, π / 4, 1.5])
        @test div_σ isa AbstractTens
        @test norm(Array(div_σ)) < 1.0e-13   # equilibrium
    end

    # ─────────────────────────────────────────────────────────────
    # Generic constructor from OM function
    # ─────────────────────────────────────────────────────────────
    @testsection "Generic OM constructor — cylindrical" begin
        # Build cylindrical from OM(r,θ,z) = (r·cosθ, r·sinθ, z)
        OM_cyl = x -> [x[1] * cos(x[2]), x[1] * sin(x[2]), x[3]]
        CS_gen = CoorSystemNum(OM_cyl, 3)
        CS_ref = coorsys_cylindrical_num()

        r₀, θ₀, z₀ = 2.0, π / 6, 1.0
        x₀ = [r₀, θ₀, z₀]

        Γ_gen = CS_gen.Γ_func(x₀)
        Γ_ref = CS_ref.Γ_func(x₀)
        @test Γ_gen ≈ Γ_ref atol = 1.0e-8

        χ_gen = CS_gen.χ_func(x₀)
        χ_ref = CS_ref.χ_func(x₀)
        @test χ_gen ≈ χ_ref atol = 1.0e-10

        # Laplacian of r² in cylindrical = 4
        f = x -> x[1]^2
        @test LAPLACE(f, CS_gen)(x₀) ≈ 4.0 atol = 1.0e-7
    end

end
