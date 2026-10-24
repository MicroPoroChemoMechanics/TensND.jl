@testsection "Isotropic tensors" begin

    # ── Type predicates ────────────────────────────────────────────────────────
    @testsection "TensISO — type predicates" begin
        𝟏 = tens_Id2(Val(3), Val(Float64))
        𝕀, 𝕁, 𝕂 = iso_projectors(Val(3), Val(Float64))
        @test  isisotropic(𝕀)
        @test  isisotropic(𝟏)
        @test !istransverselyisotropic(𝕀)
        @test !istransverselyisotropic(𝟏)
        @test !isorthotropic(𝕀)
        @test !isorthotropic(𝟏)
    end

    # ── Display (show) ─────────────────────────────────────────────────────────
    @testsection "TensISO — show" begin
        𝟏 = tens_Id2(Val(3), Val(Float64))
        𝕀, 𝕁, 𝕂 = iso_projectors(Val(3), Val(Float64))
        # show should write to the provided IO, not to stdout
        buf4 = IOBuffer()
        show(buf4, 𝕁 + 𝕂)          # 4th-order
        s4 = String(take!(buf4))
        @test contains(s4, "𝕁") || contains(s4, "𝕂")

        buf2 = IOBuffer()
        show(buf2, 𝟏)               # 2nd-order
        s2 = String(take!(buf2))
        @test contains(s2, "𝟏")
    end

    for T in (Sym, Float64), dim in (2, 3)
        @testsection "type $T, dim $dim" begin
            𝟏 = tens_Id2(Val(dim), Val(T))
            @test opequal(𝟏, 1I)
            @test opequal(tr(𝟏), dim)

            if T == Sym
                α = symbols("α", real = true)
            else
                α = rand()
                while opequal(α, one(T)) || opequal(α, zero(T))
                    α = rand()
                end
            end
            @test α * 𝟏 isa TensISO{2}
            t = α * 𝟏 + (1 - α) * 𝟏
            @test opequal(t, 𝟏)
            @test t isa TensISO{2}
            @test opequal(inv(α * 𝟏), inv(α) * 𝟏)
            @test isisotropic(t)

            𝕀, 𝕁, 𝕂 = iso_projectors(Val(dim), Val(T))
            @test opequal(𝕁 + 𝕂, 𝕀)
            @test opequal((𝟏 ⊗ 𝟏) / dim, 𝕁)
            if T == Sym
                α, β = symbols("α β", real = true)
            else
                α = rand()
                β = rand()
                while opequal(α + β, one(T)) || opequal(α, zero(T)) || opequal(β, zero(T))
                    α = rand()
                    β = rand()
                end
            end
            𝕋 = α * 𝕁 + β * 𝕂
            @test 𝕋 isa TensISO{4}
            @test opequal(inv(𝕋), inv(α) * 𝕁 + inv(β) * 𝕂)
            @test isisotropic(𝕋)

            if T == Sym && dim == 3
                E, ν = symbols("E ν", real = true)
                k = E / 3(1 - 2ν)
                μ = E / 2(1 + ν)
                λ = E * ν / ((1 + ν) * (1 - 2ν))
                ℂ = tsimplify(3k * 𝕁 + 2μ * 𝕂)
                @test ℂ == tsimplify(TensISO{dim}(3k, 2μ))
                𝕊 = tsimplify(inv(ℂ))
                @test tsimplify.(mandel(𝕊)) == [
                    1 / E -ν / E -ν / E 0 0 0
                    -ν / E 1 / E -ν / E 0 0 0
                    -ν / E -ν / E 1 / E 0 0 0
                    0 0 0 (1 + ν) / E 0 0
                    0 0 0 0 (1 + ν) / E 0
                    0 0 0 0 0 (1 + ν) / E
                ]
                @test tsimplify(ℂ ⊡ 𝕊) == 𝕀

                n = 𝐞(3)
                Eᵒᵉᵈᵒ = E * (1 - ν) / ((1 + ν) * (1 - 2ν))
                Kref = tsimplify.([μ 0 0; 0 μ 0; 0 0 Eᵒᵉᵈᵒ])
                @test tfactor(n ⋅ ℂ ⋅ n) == tfactor(dotdot(n, ℂ, n)) == Kref
                # `dotdot` is TensND's own export, not only Tensors': the suite
                # also loads Tensors, which used to hide that it was missing.
                @test :dotdot in names(TensND)
                @test TensND.dotdot === Tensors.dotdot
                # Hooke law
                for i in 1:3, j in 1:3
                    @eval $(Symbol("ε$i$j")) = symbols($"ε$i$j", real = true)
                end
                𝛆 = Tens(SymmetricTensor{2, 3}((i, j) -> eval(Symbol("ε$i$j"))))
                𝛔 = ℂ ⊡ 𝛆
                @test tfactor(𝛔) == tfactor(λ * tr(𝛆) * 𝟏 + 2μ * 𝛆)
                @test tfactor(tsimplify(𝛔 ⊡ 𝛆)) == tfactor(tsimplify(λ * tr(𝛆)^2 + 2μ * 𝛆 ⊡ 𝛆))

                @test 𝕀 == 𝟏 ⊠ˢ 𝟏
                @test 3𝕁 == 𝟏 ⊗ 𝟏
                @test 𝕀 ⊙ 𝕀 == 6
                @test 𝕁 ⊙ 𝕀 == 𝕁 ⊙ 𝕁 == 1
                @test 𝕂 ⊙ 𝕀 == 𝕂 ⊙ 𝕂 == 5
                @test 𝕂 ⊙ 𝕁 == 𝕁 ⊙ 𝕂 == 0
                @test tsimplify(ℂ ⊙ 𝕁) == tsimplify(3k)
                @test tsimplify(ℂ ⊙ 𝕂) == tsimplify(10μ)


            end

        end
    end

    # ══════════════════════════════════════════════════════════════════════════
    @testsection "TensISO products return a Tens, in the operand's basis" begin
        # REGRESSION. `dcontract(::TensISO{4}, ::TensOrthonormal{2})` and
        # `dotdot(::AbstractTens{1}, ::TensISO{4}, ::AbstractTens{1})` used to
        # add a `Tens` to `LinearAlgebra.I`. A `UniformScaling` is not a tensor,
        # so the sum fell back to plain array arithmetic and returned a bare
        # `Array`: the wrapper was lost, and with it the BASIS — a result
        # computed in a rotated basis came back as an unlabeled array that the
        # caller would read as canonical components.
        #
        # The identity is now built in the operand's own basis (`_id2_like`), so
        # nothing is re-expressed in the canonical frame either.
        ℂ = TensISO{3}(3 * 30.0, 2 * 18.0)
        ℬ = RotatedBasis(0.4, 0.7, 0.3)
        mk(b) = Tens(SymmetricTensor{2, 3}((i, j) -> 1.0e-3 * (i == j ? i : 0.5)), b)

        for b in (CanonicalBasis{3, Float64}(), ℬ)
            𝛆 = mk(b)
            for r in (ℂ ⊡ 𝛆, 𝛆 ⊡ ℂ)
                @test r isa AbstractTens{2, 3}
                @test get_basis(r) == b                     # basis preserved
            end
            # …and the value matches the dense route.
            dense = TensND.generic_tens(ℂ) ⊡ 𝛆
            @test get_array(ℂ ⊡ 𝛆) ≈ get_array(dense) atol = 1.0e-14
        end

        # Structured × structured keeps the structured type.
        @test ℂ ⊡ TensISO{3}(1.0e-3) isa TensISO

        # `dotdot` with two vectors goes through the same identity.
        𝐯₁ = Tens(Vec{3}((1.0, 2.0, 3.0)))
        𝐯₂ = Tens(Vec{3}((0.5, -1.0, 2.0)))
        d = dotdot(𝐯₁, ℂ, 𝐯₂)
        @test d isa AbstractTens{2, 3}
        @test get_array(d) ≈ get_array(dotdot(𝐯₁, TensND.generic_tens(ℂ), 𝐯₂)) atol = 1.0e-12

        # Symbolic and ForwardDiff element types must survive the same path —
        # the `UniformScaling` route had its own `SymType` patches, so this is
        # where a regression would show.
        k, μ = symbols("k μ", real = true)
        ℂs = TensISO{3}(3k, 2μ)
        𝛆s = Tens(SymmetricTensor{2, 3}((i, j) -> symbols("e$(i)$(j)", real = true)))
        rs = ℂs ⊡ 𝛆s
        @test rs isa AbstractTens{2, 3}
        @test eltype(rs) <: Sym

        gd = ForwardDiff.derivative(
            x -> get_array(TensISO{3}(3x, 2 * 18.0) ⊡ mk(CanonicalBasis{3, Float64}()))[1, 1],
            30.0,
        )
        @test isfinite(gd)
    end

end


@testsection "Isotropic tensors — Kelvin-Mandel in closed form, symbolically too" begin
    # REGRESSION. `mandel(::TensISO)` used to go through
    # `tomandel(SymmetricTensor{order,dim}(A))`, and building a
    # `SymmetricTensor` from a full component array makes `Tensors` validate the
    # symmetry. That check cannot be answered on a `Symbolics.Num`, concluded
    # "not symmetric", and `mandel` — hence `show`, which goes through it — threw an
    # `InexactError` on any symbolic isotropic tensor. It broke a Documenter
    # `@example` block that merely displayed one.
    old_KM(A::TensISO{order, dim}) where {order, dim} =
        tomandel(SymmetricTensor{order, dim}(A))

    for dim in (2, 3), (a, b) in ((3.0, 2.0), (7.5, 0.25), (1.0, 1.0))
        A4 = TensISO{dim}(a, b)
        A2 = TensISO{dim}(a)
        @test mandel(A4) ≈ old_KM(A4)
        @test mandel(A2) ≈ old_KM(A2)
        @test eltype(mandel(A4)) === eltype(old_KM(A4))
        @test size(mandel(A4)) == size(old_KM(A4))
    end

    # ForwardDiff still reaches through, and the round trip closes.
    @test ForwardDiff.derivative(x -> mandel(TensISO{3}(3.0, 2x))[4, 4], 1.5) ≈ 2.0
    A = TensISO{3}(3.3, 1.7)
    @test maximum(abs, get_array(from_mandel(mandel(A))) .- get_array(A)) < 1.0e-14

    # SymPy and Symbolics both, at both orders, plus `show` — the actual trigger.
    @syms ka::positive mua::positive
    @test mandel(TensISO{3}(3ka, 2mua)) isa Matrix
    @test mandel(TensISO{3}(ka)) isa Vector
    @test repr(MIME"text/latex"(), TensISO{3}(3ka, 2mua)) isa String

    TensND.Symbolics.@variables kn mun
    An = TensISO{3}(3kn, 2mun)
    @test mandel(An) isa Matrix
    @test mandel(TensISO{3}(kn)) isa Vector
    @test repr(MIME"text/latex"(), An) isa String
    # TensND's symmetry predicates answer on a `Num` (Tensors.jl's compare with `==`).
    @test TensND._ismajorsym(Tensor{4, 3}(get_array(An)))
    @test TensND._isminorsym(Tensor{4, 3}(get_array(An)))
end

@testsection "Isotropic tensors — exact and integer element types" begin
    # `mandel` sized its matrix to the type of the moduli and added the spherical
    # part `(α - β) / dim` into it: integer moduli threw `InexactError`.
    @test mandel(TensISO{3}(3, 2)) == mandel(TensISO{3}(3.0, 2.0))
    @test mandel(TensISO{3}(3 // 1, 2 // 1)) isa Matrix{Rational{Int}}

    # `isotropify` divided the deviatoric part by `5.0`, a float: rational
    # moduli came back as floats and symbolic ones as `2.0*μ`.
    @test get_data(isotropify(get_array(TensISO{3}(3 // 1, 2 // 1)))) === (3 // 1, 2 // 1)
    @syms k::positive μ::positive
    β = get_data(isotropify(get_array(TensISO{3}(3k, 2μ))))[2]
    @test string(β) == "2*μ"
end

@testsection "isotropify in dimensions 2 and 4" begin
    # Exact on an isotropic array of any dimension: the moduli come back.
    for d in (2, 4)
        A = get_array(TensISO{d}(3 // 1, 2 // 1))
        @test get_data(isotropify(A)) === (3 // 1, 2 // 1)
        @test isotropify(A) isa TensISO{4, d}
    end
end

@testsection "TensISO with a UniformScaling, and the double contraction of two" begin
    A = TensISO{3}(2.0, 5.0)
    𝕀 = Array(one(SymmetricTensor{4, 3}))
    @test Array(get_array(A + 3I)) ≈ Array(get_array(A)) .+ 3 .* 𝕀
    @test Array(get_array(A - 3I)) ≈ Array(get_array(A)) .- 3 .* 𝕀
    @test Array(get_array(A * 3I)) ≈ 3 .* Array(get_array(A))
    a, b = TensISO{3}(4.0), TensISO{3}(5.0)
    @test dcontract(a, b) ≈ dcontract(get_array(a), get_array(b)) ≈ 60
end

@testsection "Isotropic elastic moduli" begin
    k, μ = 12.0, 5.0
    ℂ = iso_stiffness(k, μ)
    @test ℂ isa TensISO{4, 3}
    @test all(bulk_shear(ℂ) .≈ (k, μ))
    # Independent reference: the Lamé form, λ = k − 2μ/3,
    # E = μ(3λ + 2μ)/(λ + μ), ν = λ/(2(λ + μ)).
    λ = k - 2μ / 3
    E, ν = young_poisson(ℂ)
    @test E ≈ μ * (3λ + 2μ) / (λ + μ)
    @test ν ≈ λ / (2 * (λ + μ))
    @test iso_stiffness_young_poisson(E, ν) ≈ ℂ
    # The stiffness is 3k 𝕁 + 2μ 𝕂, component by component.
    @test get_array(ℂ) ≈ 3k * get_array(tens_J4(Val(3), Val(Float64))) + 2μ * get_array(tens_K4(Val(3), Val(Float64)))
    # A compliance: the moduli of its inverse.
    @test all(bulk_shear(inv(inv(ℂ))) .≈ (k, μ))

    # Derivatives: dE/dμ = 27k²/(3k + μ)², d²E/dμ² = −54k²/(3k + μ)³.
    Eμ(m) = young_poisson(iso_stiffness(k, m))[1]
    @test ForwardDiff.derivative(Eμ, μ) ≈ 27k^2 / (3k + μ)^2
    @test ForwardDiff.derivative(m -> ForwardDiff.derivative(Eμ, m), μ) ≈ -54k^2 / (3k + μ)^3
    @test ForwardDiff.derivative(e -> young_poisson(iso_stiffness_young_poisson(e, 0.3))[1], 7.0) ≈ 1

    # Symbolic: the closed forms, exactly.
    @syms ks::positive μs::positive
    Es, νs = young_poisson(iso_stiffness(ks, μs))
    @test simplify(Es - 9ks * μs / (3ks + μs)) == 0
    @test simplify(νs - (3ks - 2μs) / (2 * (3ks + μs))) == 0
    @test all(simplify.(collect(bulk_shear(iso_stiffness(ks, μs))) .- [ks, μs]) .== 0)
    @test simplify(poisson_ratio(ks, μs) - νs) == 0
    @test all(shear_poisson(ℂ) .≈ (μ, ν))
end

@testsection "Isotropic moduli in 2-D: the plane-strain restriction" begin
    # Independent reference: the in-plane components of a 3-D isotropic
    # stiffness, read back as a 2-D isotropic tensor α 𝕁 + β 𝕂 by
    # α = C₁₁₁₁ + C₁₁₂₂, β = C₁₁₁₁ − C₁₁₂₂. Its moduli are those of the 3-D
    # material: ν, E and μ, and the in-plane bulk modulus λ + μ.
    for (k, μ) in ((12.0, 5.0), (3.0, 2.5), (40.0, 1.0))
        C3 = get_array(iso_stiffness(k, μ))
        C2 = TensISO{2}(C3[1, 1, 1, 1] + C3[1, 1, 2, 2], C3[1, 1, 1, 1] - C3[1, 1, 2, 2])
        @test C2 isa TensISO{4, 2}
        @test all(young_poisson(C2) .≈ young_poisson(iso_stiffness(k, μ)))
        @test all(shear_poisson(C2) .≈ shear_poisson(iso_stiffness(k, μ)))
        @test all(bulk_shear(C2) .≈ (k - 2μ / 3 + μ, μ))
        # The 2-D tensor is the restriction, component by component.
        @test get_array(C2) ≈ C3[1:2, 1:2, 1:2, 1:2]
    end
    # Derivatives along k, through the restriction; symbols, exactly.
    ν2(k) = shear_poisson(TensISO{2}(2 * (k + 5.0 / 3), 10.0))[2]
    @test ForwardDiff.derivative(ν2, 12.0) ≈ ForwardDiff.derivative(k -> poisson_ratio(k, 5.0), 12.0)
    @syms ks::positive μs::positive
    a, b = 2 * (ks + μs / 3), 2μs
    @test simplify(shear_poisson(TensISO{2}(a, b))[2] - poisson_ratio(ks, μs)) == 0
    @test simplify(young_poisson(TensISO{2}(a, b))[1] - young_poisson(iso_stiffness(ks, μs))[1]) == 0
end

@testsection "Components of an isotropic tensor with variances" begin
    # Canonical basis: covariant and contravariant components coincide, so the
    # variance pattern changes nothing. The reference is the same tensor stored
    # as a generic `Tens`, which goes through the general change of variance.
    for (A, var) in (
            (TensISO{3}(3 * 7.0, 2 * 3.0), (:cov, :cont, :cov, :cont)),
            (TensISO{3}(2.0), (:cov, :cont)),
        )
        @test components(A, var) == get_array(A)
        @test components(A, var) ≈ components(Tens(get_array(A)), var)
    end
end
