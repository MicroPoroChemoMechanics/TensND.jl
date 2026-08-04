import Pkg
Pkg.activate(joinpath(@__DIR__, "..", "docs"); io = devnull)
# The `docs` environment declares every dependency the tutorials use
# (SymPy, Symbolics, NLopt, BenchmarkTools, ...) and resolves TensND
# itself through `[sources]`. These four lines are stripped from the
# generated page and notebook, which already run inside that project.

using TensND
using LinearAlgebra
using SymPy

Spherical = coorsys_spherical()
θ, ϕ, r = getcoords(Spherical)
𝐞ᶿ, 𝐞ᵠ, 𝐞ʳ = unitvec(Spherical)
@set_coorsys Spherical

𝐞₁, 𝐞₂, 𝐞₃ = unitvec(coorsys_cartesian())
𝕀, 𝕁, 𝕂 = iso_projectors(Val(3), Val(Sym))
𝟏 = tens_Id2(Val(3), Val(Sym))

k, μ = symbols("k μ", positive = true)
ℂ = 3k * 𝕁 + 2μ * 𝕂

u = SymFunction("u", real = true)
𝐮 = u(r) * 𝐞ʳ

𝛆 = SYMGRAD(𝐮)
𝛔 = ℂ ⊡ 𝛆
𝐓 = 𝛔 ⋅ 𝐞ʳ
eq = factor(simplify(DIV(𝛔) ⋅ 𝐞ʳ))

sol = dsolve(eq, u(r))
û = sol.rhs()

T̂ = tsimplify(tsimplify(subs(𝐓 ⋅ 𝐞ʳ, u(r) => û)))

remote_angle_functions(𝐄) = let fʳ = simplify(𝐞ʳ ⋅ 𝐄 ⋅ 𝐞ʳ)
    (diff(fʳ, θ) / 2, diff(fʳ, ϕ) / (2sin(θ)), fʳ)
end

uᶿ = SymFunction("uᶿ", real = true)
uᵠ = SymFunction("uᵠ", real = true)
uʳ = SymFunction("uʳ", real = true)
α, Λ = symbols("α Λ", real = true)

fᶿ, _, fʳ = remote_angle_functions(𝟏 - 3𝐞₃ ⊗ 𝐞₃)
(fᶿ, fʳ)

𝐮ᵈ = uᶿ(r) * fᶿ * 𝐞ᶿ + uʳ(r) * fʳ * 𝐞ʳ
𝛔ᵈ = ℂ ⊡ SYMGRAD(𝐮ᵈ)
𝐓ᵈ = 𝛔ᵈ ⋅ 𝐞ʳ
div𝛔ᵈ = DIV(𝛔ᵈ)

eqᶿ = tsimplify(div𝛔ᵈ ⋅ 𝐞ᶿ / fᶿ)
eqʳ = tsimplify(div𝛔ᵈ ⋅ 𝐞ʳ / fʳ)

eqs = tsimplify.(subs.([eqᶿ, eqʳ], uᶿ(r) => r^α, uʳ(r) => Λ * r^α))
αΛ = solve([e.doit() for e in eqs], [α, Λ])

[(pair[1], simplify(pair[2])) for pair in αΛ]

ûᶿ = sum(Sym("C$(i + 2)") * r^αΛ[i][1] for i in 1:length(αΛ))
ûʳ = sum(Sym("C$(i + 2)") * αΛ[i][2] * r^αΛ[i][1] for i in 1:length(αΛ))
(ûᶿ, ûʳ)

T̂ᶿ = tsimplify(tsimplify(subs(simplify(𝐓ᵈ ⋅ 𝐞ᶿ / fᶿ), uᶿ(r) => ûᶿ, uʳ(r) => ûʳ)))

T̂ʳ = tsimplify(tsimplify(subs(simplify(𝐓ᵈ ⋅ 𝐞ʳ / fʳ), uᶿ(r) => ûᶿ, uʳ(r) => ûʳ)))

fᶿ₂, fᵠ₂, fʳ₂ = remote_angle_functions(𝐞₁ ⊗ 𝐞₁ - 𝐞₂ ⊗ 𝐞₂)
(fᶿ₂, fᵠ₂, fʳ₂)

𝐮ˢ = uᶿ(r) * fᶿ₂ * 𝐞ᶿ + uᵠ(r) * fᵠ₂ * 𝐞ᵠ + uʳ(r) * fʳ₂ * 𝐞ʳ
div𝛔ˢ = DIV(ℂ ⊡ SYMGRAD(𝐮ˢ))

eqᶿˢ = tsimplify(div𝛔ˢ ⋅ 𝐞ᶿ / fᶿ₂)
eqᵠˢ = tsimplify(div𝛔ˢ ⋅ 𝐞ᵠ / fᵠ₂)
eqʳˢ = tsimplify(div𝛔ˢ ⋅ 𝐞ʳ / fʳ₂)

X = symbols("X", real = true)
uᵠsol = solve(tsimplify(diff(subs(eqᵠˢ, sin(θ)^2 => 1 / X), X)), uᵠ(r))[1]

eqs₂ = tsimplify.(subs.([eqᶿˢ, eqʳˢ], uᵠ(r) => r^α, uᶿ(r) => r^α, uʳ(r) => Λ * r^α))
αΛ₂ = solve([e.doit() for e in eqs₂], [α, Λ])

sort(string.(first.(αΛ))) == sort(string.(first.(αΛ₂)))

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
