@testsection "Coordinate systems" begin
    s∂ = tsimplify ∘ ∂
    (x, y, z), (𝐞₁, 𝐞₂, 𝐞₃), ℬ = init_cartesian()
    (θ, ϕ, r), (𝐞ᶿ, 𝐞ᵠ, 𝐞ʳ), ℬˢ = init_spherical()

    @testsection "Usual coordinate systems" begin
        @test components(𝐞ʳ ⊗ 𝐞ᵠ, ℬˢ) == components(𝐞₃ ⊗ 𝐞₂, ℬ) == canonical_components(𝐞₃ ⊗ 𝐞₂)
    end

    @testsection "Partial derivatives" begin
        @test s∂(𝐞ʳ, θ) == 𝐞ᶿ
        @test s∂(𝐞ʳ, ϕ) == sin(θ) * 𝐞ᵠ
        @test s∂(𝐞ᵠ ⊗ 𝐞ᶿ, ϕ) == s∂(𝐞ᵠ, ϕ) ⊗ 𝐞ᶿ + 𝐞ᵠ ⊗ s∂(𝐞ᶿ, ϕ)
        @test s∂(𝐞ʳ ⊗ˢ 𝐞ᵠ, ϕ) == s∂(𝐞ʳ, ϕ) ⊗ˢ 𝐞ᵠ + 𝐞ʳ ⊗ˢ s∂(𝐞ᵠ, ϕ)
    end

    @testsection "Coordinate systems" begin
        # Cartesian
        Cartesian = coorsys_cartesian()
        𝐗 = get_coords(Cartesian)
        𝐄 = unitvec(Cartesian)
        ℬ = normalized_basis(Cartesian)
        𝛔 = Tens(SymmetricTensor{2, 3}((i, j) -> SymFunction("σ$i$j", real = true)(𝐗...)))
        @test divergence(𝛔, Cartesian) ==
            sum([sum([∂(𝛔[i, j], 𝐗[j]) for j in 1:3]) * 𝐄[i] for i in 1:3])

        # Polar
        Polar = coorsys_polar()
        r, θ = get_coords(Polar)
        𝐞ʳ, 𝐞ᶿ = unitvec(Polar)
        ℬᵖ = normalized_basis(Polar)
        f = SymFunction("f", real = true)(r, θ)
        @test tsimplify(laplacian(f, Polar)) ==
            tsimplify(∂(r * ∂(f, r), r) / r + ∂(f, θ, θ) / r^2)

        # Cylindrical
        Cylindrical = coorsys_cylindrical()
        rθz = get_coords(Cylindrical)
        𝐞ʳ, 𝐞ᶿ, 𝐞ᶻ = unitvec(Cylindrical)
        ℬᶜ = normalized_basis(Cylindrical)
        r, θ, z = rθz
        𝐯 = Tens(Vec{3}(i -> SymFunction("v$(rθz[i])", real = true)(rθz...)), ℬᶜ)
        vʳ, vᶿ, vᶻ = get_array(𝐯)
        @test tsimplify(divergence(𝐯, Cylindrical)) ==
            tsimplify(∂(vʳ, r) + vʳ / r + ∂(vᶿ, θ) / r + ∂(vᶻ, z))

        # Spherical
        Spherical = coorsys_spherical()
        θ, ϕ, r = get_coords(Spherical)
        𝐞ᶿ, 𝐞ᵠ, 𝐞ʳ = unitvec(Spherical)
        ℬˢ = normalized_basis(Spherical)
        for σⁱʲ in ("σʳʳ", "σᶿᶿ", "σᵠᵠ")
            @eval $(Symbol(σⁱʲ)) = SymFunction($σⁱʲ, real = true)($r)
        end
        𝛔 = σʳʳ * 𝐞ʳ ⊗ 𝐞ʳ + σᶿᶿ * 𝐞ᶿ ⊗ 𝐞ᶿ + σᵠᵠ * 𝐞ᵠ ⊗ 𝐞ᵠ
        div𝛔 = tsimplify(divergence(𝛔, Spherical))
        @test tsimplify(div𝛔 ⋅ 𝐞ʳ) == tsimplify(∂(σʳʳ, r) + (2σʳʳ - σᶿᶿ - σᵠᵠ) / r)

        # Concentric sphere - hydrostatic part
        θ, ϕ, r = get_coords(Spherical)
        𝐞ᶿ, 𝐞ᵠ, 𝐞ʳ = unitvec(Spherical)
        ℬˢ = normalized_basis(Spherical)
        𝕀, 𝕁, 𝕂 = iso_projectors(Val(3), Val(Sym))
        𝟏 = tens_Id2(Val(3), Val(Sym))
        k, μ = symbols("k μ", positive = true)
        λ = k - 2μ / 3
        ℂ = 3k * 𝕁 + 2μ * 𝕂
        u = SymFunction("u", real = true)(r)
        𝐮 = u * 𝐞ʳ
        𝛆 = tsimplify(symgrad(𝐮, Spherical))
        𝛔 = tsimplify(ℂ ⊡ 𝛆)
        # 𝛔 = tsimplify(λ * tr(𝛆) * 𝟏 + 2μ * 𝛆)
        @test dsolve(tfactor(tsimplify(divergence(𝛔, Spherical) ⋅ 𝐞ʳ)), u) ==
            Eq(u, symbols("C1") / r^2 + symbols("C2") * r)

        # Spheroidal
        Spheroidal = coorsys_spheroidal()
        OM = position_vector(Spheroidal)
        @test tsimplify(laplacian(OM[1]^2, Spheroidal)) == 2


    end


    @testsection "The default chart" begin
        # The chart, its letter and its coordinate names were held in three
        # references, set one after the other, and `unset_coorsys!` forgot only
        # two of them: the letter of a forgotten chart stayed in `pprint`.
        printed(x) = mktemp() do path, io
            redirect_stdout(() -> pprint(x), io)
            flush(io)
            read(path, String)
        end
        unset_coorsys!()
        @test_throws ErrorException default_coorsys()
        Polar = coorsys_polar()
        r, θ = get_coords(Polar)
        f = SymFunction("f", real = true)(r, θ)
        @test set_coorsys!(Polar; vec = '𝐮') === Polar
        @test default_coorsys() === Polar
        @test grad(f) == grad(f, Polar)
        v = Tens(Sym[1, 0])
        @test occursin('𝐮', printed(v))
        unset_coorsys!()
        @test_throws ErrorException default_coorsys()
        @test !occursin('𝐮', printed(v))
    end
end
