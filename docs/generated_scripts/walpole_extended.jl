import Pkg
Pkg.activate(joinpath(@__DIR__, "..", "docs"); io = devnull)
# The `docs` environment declares every dependency the tutorials use
# (SymPy, Symbolics, NLopt, BenchmarkTools, ...) and resolves TensND
# itself through `[sources]`. These four lines are stripped from the
# generated page and notebook, which already run inside that project.

using TensND
using LinearAlgebra
using Printf
using Random

Random.seed!(20260804)
n = (0.0, 0.0, 1.0)
arr(x) = get_array(x)

W7, W8 = tens_W7(n), tens_W8(n)

function show_KM(label, t)
    M = KM(t)
    println(label)
    for r in 1:6
        println("   ", join([abs(M[r, c]) < 1.0e-12 ? "   .   " : @sprintf("%7.4f", M[r, c]) for c in 1:6], " "))
    end
    return println()
end

show_KM("𝕎7 =", W7)
show_KM("𝕎8 =", W8)

(
    norm(arr(W7) + permutedims(arr(W7), (3, 4, 1, 2))),
    norm(arr(W8) + permutedims(arr(W8), (3, 4, 1, 2))),
)

𝟏ₙ = [0.0 0 0; 0 0 0; 0 0 1]
𝟏ₜ = [1.0 0 0; 0 1 0; 0 0 0]
gen = (x -> (x + x') / 2)(rand(3, 3))          # a general symmetric tensor
ti = 3.0 * 𝟏ₜ + 5.0 * 𝟏ₙ                       # axially invariant

for (name, M) in (("𝟏ₙ", 𝟏ₙ), ("𝟏ₜ", 𝟏ₜ), ("3·𝟏ₜ + 5·𝟏ₙ", ti), ("generic symmetric", gen))
    @printf "%-20s  ‖𝕎₇:M‖ = %.3e   ‖𝕎₈:M‖ = %.3e\n" name norm(arr(W7) ⊡ M) norm(arr(W8) ⊡ M)
end

A = TensTI{4}(2.0, 1.0, 0.5, 0.7, 0.9, 1.3, 0.4, 0.6, n)
B = TensTI{4}(3.0, 1.2, 0.3, 0.8, 1.1, 0.5, 0.2, 0.9, n)

println("typeof(A)   = ", typeof(A))
println("typeof(A⊡B) = ", typeof(A ⊡ B))
println("‖A⊡B − dense‖ = ", norm(arr(A ⊡ B) - arr(A) ⊡ arr(B)))
println("‖A⊡inv(A) − 𝕀‖ = ", norm(arr(A ⊡ inv(A)) - arr(tens_Id4(Val(3), Val(Float64)))))

ℓA, ℓB, ℓP = get_ℓ8(A), get_ℓ8(B), get_ℓ8(A ⊡ B)
LA = [ℓA[1] ℓA[3]; ℓA[4] ℓA[2]]
LB = [ℓB[1] ℓB[3]; ℓB[4] ℓB[2]]

println("2×2 block : ", round.(LA * LB, digits = 10), "  vs  ", round.([ℓP[1] ℓP[3]; ℓP[4] ℓP[2]], digits = 10))
println("z₁ = ℓ₆+iℓ₇ : ", round((ℓA[6] + im * ℓA[7]) * (ℓB[6] + im * ℓB[7]), digits = 10), "  vs  ", round(ℓP[6] + im * ℓP[7], digits = 10))
println("z₂ = ℓ₅+iℓ₈ : ", round((ℓA[5] + im * ℓA[8]) * (ℓB[5] + im * ℓB[8]), digits = 10), "  vs  ", round(ℓP[5] + im * ℓP[8], digits = 10))

ℓI = get_ℓ8(inv(A))
LI = [ℓI[1] ℓI[3]; ℓI[4] ℓI[2]]
println("L⁻¹  : ", round.(inv(LA) - LI, digits = 12))
println("1/z₁ : ", round(1 / (ℓA[6] + im * ℓA[7]) - (ℓI[6] + im * ℓI[7]), digits = 12))
println("1/z₂ : ", round(1 / (ℓA[5] + im * ℓA[8]) - (ℓI[5] + im * ℓI[8]), digits = 12))

println("get_ℓ(A)  = ", round.(get_ℓ(A), digits = 4))
println("get_ℓ8(A) = ", round.(get_ℓ8(A), digits = 4))

A5 = TensTI{4}(2.0, 1.0, 0.5, 0.9, 1.3, n)
println("get_ℓ8 of an N=5 tensor = ", round.(get_ℓ8(A5), digits = 6))

Bfit, d, drel = proj_tens(:TI, arr(A), [0.0, 0.0, 1.0])
@printf "best fit → (ℓ₁,ℓ₂,ℓ₃,ℓ₅,ℓ₆) = %s\n" string(round.(get_data(Bfit), digits = 6))
@printf "relative distance discarded  = %.4f\n" drel
println("ℓ₃−ℓ₄ in A : ", round(ℓA[3] - ℓA[4], digits = 6))
println("ℓ₇, ℓ₈ in A : ", round(ℓA[7], digits = 6), ", ", round(ℓA[8], digits = 6))

A2 = TensTI{2}(3.0, 5.0, 1.7, n)
println("typeof = ", typeof(A2))
arr(A2)

A2s = TensTI{2}(3.0, 5.0, n)
println("typeof = ", typeof(A2s))
arr(A2s)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
