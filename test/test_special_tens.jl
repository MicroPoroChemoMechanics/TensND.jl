@testsection "Special tensors" begin

    # ── Rotations ────────────────────────────────────────────────────────────

    @testset "rot3 / rot6 consistency" begin
        θ, ϕ, ψ = symbols("θ ϕ ψ", real = true)
        cθ, cϕ, cψ, sθ, sϕ, sψ = symbols("cθ cϕ cψ sθ sϕ sψ", real = true)
        d = Dict(cos(θ) => cθ, cos(ϕ) => cϕ, cos(ψ) => cψ, sin(θ) => sθ, sin(ϕ) => sϕ, sin(ψ) => sψ)
        R = Tens(tsubs(rot3(θ, ϕ, ψ), d...))
        R6 = inv_KM(tsubs(KM(rot6(θ, ϕ, ψ)), d...))
        @test R6 == R ⊠ˢ R
    end

    @testset "rot3 is Z-Y-Z and its third column is the axis" begin
        θ, ϕ, ψ = 0.37, 0.83, 0.29
        R = Matrix(rot3(θ, ϕ, ψ))
        @test R ≈ RotZ(ϕ) * RotY(θ) * RotZ(ψ)
        @test R[:, 3] ≈ [sin(θ)cos(ϕ), sin(θ)sin(ϕ), cos(θ)]
        @test R'R ≈ I
        @test det(R) ≈ 1
        # The columns of R are the vectors of the corresponding RotatedBasis.
        @test vecbasis(Basis(θ, ϕ, ψ), :cov) ≈ R
    end

    @testset "rot3 parametrization degeneracies" begin
        θ, ϕ, ψ = 0.41, 1.13, 0.77
        # (θ, ϕ, ψ) and (-θ, ϕ+π, ψ+π) name the same rotation.
        @test Matrix(rot3(θ, ϕ, ψ)) ≈ Matrix(rot3(-θ, ϕ + π, ψ + π))
        # At θ = 0 only ϕ + ψ matters; at θ = π only ϕ - ψ.
        @test Matrix(rot3(0.0, 0.3, 0.9)) ≈ Matrix(rot3(0.0, 0.7, 0.5))
        @test Matrix(rot3(π, 0.9, 0.3)) ≈ Matrix(rot3(π, 1.1, 0.5))
    end

    @testset "rot2" begin
        θ = 0.6
        @test get_array(Tens(rot2(θ))) ≈ [cos(θ) -sin(θ); sin(θ) cos(θ)]
        @test Matrix(get_array(Tens(rot2(θ)))) * Matrix(get_array(Tens(rot2(-θ)))) ≈ I
    end

    # ── Levi-Civita ──────────────────────────────────────────────────────────

    @testset "LeviCivita" begin
        a = LeviCivita(Float64)
        @test a[1, 2, 3] == 1
        @test a[2, 3, 1] == 1 && a[3, 1, 2] == 1          # even permutations
        @test a[2, 1, 3] == -1 && a[1, 3, 2] == -1        # odd permutations
        @test a[1, 1, 2] == 0 && a[3, 3, 3] == 0          # repeated indices
        # Total antisymmetry.
        @test a ≈ -permutedims(a, (2, 1, 3))
        @test a ≈ -permutedims(a, (1, 3, 2))
        # The ε–δ identity: εᵢⱼₖ εₗₘₖ = δᵢₗδⱼₘ − δᵢₘδⱼₗ.
        δ = Matrix(1.0I, 3, 3)
        lhs = [sum(a[i, j, k] * a[l, m, k] for k in 1:3) for i in 1:3, j in 1:3, l in 1:3, m in 1:3]
        rhs = [δ[i, l] * δ[j, m] - δ[i, m] * δ[j, l] for i in 1:3, j in 1:3, l in 1:3, m in 1:3]
        @test lhs ≈ rhs
        # The cross product it generates.
        u, v = [1.0, 2.0, 3.0], [4.0, 5.0, 6.0]
        @test [sum(a[i, j, k] * u[j] * v[k] for j in 1:3, k in 1:3) for i in 1:3] ≈ cross(u, v)
    end

    # ── Predefined unit vectors ──────────────────────────────────────────────

    @testset "𝐞 (canonical)" begin
        for i in 1:3
            e = 𝐞(Val(i), Val(3), Val(Float64))
            @test get_array(e) ≈ [j == i ? 1.0 : 0.0 for j in 1:3]
        end
    end

    @testset "𝐞ᵖ, 𝐞ᶜ, 𝐞ˢ are orthonormal frames" begin
        θs = symbols("θ", real = true)
        # Polar
        p1, p2 = 𝐞ᵖ(Val(1), θs), 𝐞ᵖ(Val(2), θs)
        @test tsimplify(p1 ⋅ p1) == 1
        @test tsimplify(p2 ⋅ p2) == 1
        @test tsimplify(p1 ⋅ p2) == 0
        # Cylindrical
        c1, c2, c3 = 𝐞ᶜ(Val(1), θs), 𝐞ᶜ(Val(2), θs), 𝐞ᶜ(Val(3), θs)
        for (x, y) in ((c1, c1), (c2, c2), (c3, c3))
            @test tsimplify(x ⋅ y) == 1
        end
        for (x, y) in ((c1, c2), (c1, c3), (c2, c3))
            @test tsimplify(x ⋅ y) == 0
        end
    end

    @testset "𝐞ˢ spherical frame" begin
        θs, ϕs = symbols("θ ϕ", real = true)
        s1, s2, s3 = 𝐞ˢ(Val(1), θs, ϕs), 𝐞ˢ(Val(2), θs, ϕs), 𝐞ˢ(Val(3), θs, ϕs)
        for x in (s1, s2, s3)
            @test tsimplify(x ⋅ x) == 1
        end
        for (x, y) in ((s1, s2), (s1, s3), (s2, s3))
            @test tsimplify(x ⋅ y) == 0
        end
        # The third vector is the radial one, and θ = ϕ = 0 gives the canonical
        # basis *in the canonical order* — the reason for the (θ, ϕ, r) ordering.
        @test tsimplify(tsubs(components_canon(s3), θs => Sym(0), ϕs => Sym(0))) == [0, 0, 1]
        @test tsimplify(tsubs(components_canon(s1), θs => Sym(0), ϕs => Sym(0))) == [1, 0, 0]
        @test tsimplify(tsubs(components_canon(s2), θs => Sym(0), ϕs => Sym(0))) == [0, 1, 0]
    end

    # ── init_* constructors ──────────────────────────────────────────────────

    @testset "init_* return (coords, vectors, basis)" begin
        for (f, dim) in (
                (init_polar, 2),
                (init_cylindrical, 3),
                (init_spherical, 3),
                (init_rotated, 3),
            )
            coords, vecs, ℬ = f()
            @test length(coords) == dim
            @test length(vecs) == dim
            @test ℬ isa TensND.AbstractBasis
            @test get_dim(ℬ) == dim
            @test isorthonormal(ℬ)
        end

        coords, vecs, ℬ = init_cartesian()
        @test length(coords) == 3 && length(vecs) == 3
        @test ℬ isa CanonicalBasis
    end

    @testset "init_* `canonical` keyword" begin
        # canonical = false (default): components are given in the local frame,
        # so the i-th vector is the i-th column of the identity.
        _, (𝐞ʳ, 𝐞ᶿ), _ = init_polar()
        @test get_array(𝐞ʳ) == [1, 0]
        # canonical = true: components are given in the canonical basis.
        (r, θ), (𝐞ʳᶜ, _), _ = init_polar(; canonical = true)
        @test tsimplify(get_array(𝐞ʳᶜ)) == [cos(θ), sin(θ)]
    end

    @testset "init_spherical ordering is (θ, ϕ, r)" begin
        (θ, ϕ, r), (𝐞ᶿ, 𝐞ᵠ, 𝐞ʳ), ℬˢ = init_spherical()
        @test string(θ) == "θ" && string(ϕ) == "ϕ" && string(r) == "r"
        # 𝐞ʳ is the radial direction.
        @test tsimplify(components_canon(𝐞ʳ)) ==
            [sin(θ)cos(ϕ), sin(θ)sin(ϕ), cos(θ)]
    end

end
