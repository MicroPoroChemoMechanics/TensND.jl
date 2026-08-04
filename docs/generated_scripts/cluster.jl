import Pkg
Pkg.activate(joinpath(@__DIR__, "..", "docs"); io = devnull)
# The `docs` environment declares every dependency the tutorials use
# (SymPy, Symbolics, NLopt, BenchmarkTools, ...) and resolves TensND
# itself through `[sources]`. These four lines are stripped from the
# generated page and notebook, which already run inside that project.

using TensND
using LinearAlgebra
using SymPy

𝕀, 𝕁, 𝕂 = iso_projectors(Val(3), Val(Sym))
𝟏 = tens_Id2(Val(3), Val(Sym))

E = symbols("E", positive = true)
ν = symbols("ν", real = true)
k = E / (3(1 - 2ν))
μ = E / (2(1 + ν))

θ, ϕ, r = symbols("θ", real = true), symbols("ϕ", real = true), symbols("r", positive = true)
S = coorsys_spherical((θ, ϕ, r))
𝐱 = getOM(S)

θ′, ϕ′, r′ = symbols("θ′", real = true), symbols("ϕ′", real = true), symbols("r′", positive = true)
S′ = coorsys_spherical((θ′, ϕ′, r′))

Cartesian = coorsys_cartesian(symbols("x y z", real = true))
𝐞₁, 𝐞₂, 𝐞₃ = unitvec(Cartesian)
R = symbols("R", positive = true)

𝐱′ = getOM(S′) + R * 𝐞₃

Δ𝐱 = 𝐱 - 𝐱′
ρ = tsimplify(norm(Δ𝐱))

ρ

𝐍 = Δ𝐱 / ρ
𝐧 = tsimplify(change_tens_canon(𝐍))

ℾ = 1 / (16PI * μ * (1 - ν) * ρ^3) * (
    -3𝕁 + 2(1 - 2ν) * 𝕀 + 3(𝟏 ⊗ 𝐧 ⊗ 𝐧 + 𝐧 ⊗ 𝐧 ⊗ 𝟏)
        + 12ν * 𝐧 ⊗ˢ 𝟏 ⊗ˢ 𝐧 - 15𝐧 ⊗ 𝐧 ⊗ 𝐧 ⊗ 𝐧
)

tsimplify(𝐧 ⋅ 𝐧)

tsimplify(subs(ρ, θ => Sym(0), θ′ => Sym(0)))

tsimplify(subs(get_array(𝐧)[3], θ => Sym(0), θ′ => Sym(0)))

ρ_same = tsimplify(subs(ρ, R => Sym(0), θ′ => θ, ϕ′ => ϕ))
tsimplify(ρ_same - abs(r - r′))

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
