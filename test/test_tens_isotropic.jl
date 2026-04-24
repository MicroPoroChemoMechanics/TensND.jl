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
end
