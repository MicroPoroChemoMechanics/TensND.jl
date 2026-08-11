@testsection "Isotropic tensors" begin

    # ── Type predicates ────────────────────────────────────────────────────────
    @testsection "TensISO — type predicates" begin
        𝟏 = tens_Id2(Val(3), Val(Float64))
        𝕀, 𝕁, 𝕂 = ISO(Val(3), Val(Float64))
        @test  is_ISO(𝕀)
        @test  is_ISO(𝟏)
        @test !is_TI(𝕀)
        @test !is_TI(𝟏)
        @test !is_ORTHO(𝕀)
        @test !is_ORTHO(𝟏)
    end

    # ── Display (show) ─────────────────────────────────────────────────────────
    @testsection "TensISO — show" begin
        𝟏 = tens_Id2(Val(3), Val(Float64))
        𝕀, 𝕁, 𝕂 = ISO(Val(3), Val(Float64))
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
            @test is_ISO(t)

            𝕀, 𝕁, 𝕂 = ISO(Val(dim), Val(T))
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
            @test is_ISO(𝕋)

            if T == Sym && dim == 3
                E, ν = symbols("E ν", real = true)
                k = E / 3(1 - 2ν)
                μ = E / 2(1 + ν)
                λ = E * ν / ((1 + ν) * (1 - 2ν))
                ℂ = tsimplify(3k * 𝕁 + 2μ * 𝕂)
                @test ℂ == tsimplify(TensISO{dim}(3k, 2μ))
                𝕊 = tsimplify(inv(ℂ))
                @test tsimplify.(KM(𝕊)) == [
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
            dense = TensND._generic_tens(ℂ) ⊡ 𝛆
            @test get_array(ℂ ⊡ 𝛆) ≈ get_array(dense) atol = 1.0e-14
        end

        # Structured × structured keeps the structured type.
        @test ℂ ⊡ TensISO{3}(1.0e-3) isa TensISO

        # `dotdot` with two vectors goes through the same identity.
        𝐯₁ = Tens(Vec{3}((1.0, 2.0, 3.0)))
        𝐯₂ = Tens(Vec{3}((0.5, -1.0, 2.0)))
        d = dotdot(𝐯₁, ℂ, 𝐯₂)
        @test d isa AbstractTens{2, 3}
        @test get_array(d) ≈ get_array(dotdot(𝐯₁, TensND._generic_tens(ℂ), 𝐯₂)) atol = 1.0e-12

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
