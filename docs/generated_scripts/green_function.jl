import Pkg
Pkg.activate(joinpath(@__DIR__, "..", "docs"); io = devnull)
# The `docs` environment declares every dependency the tutorials use
# (SymPy, Symbolics, NLopt, BenchmarkTools, ...) and resolves TensND
# itself through `[sources]`. These four lines are stripped from the
# generated page and notebook, which already run inside that project.

using TensND
using LinearAlgebra
using SymPy
using Tensors

Polar = coorsys_polar()
r, θ = getcoords(Polar)
𝐞ʳ, 𝐞ᶿ = unitvec(Polar)
@set_coorsys Polar
ℬᵖ = normalized_basis(Polar)

𝕀₂, 𝕁₂, 𝕂₂ = iso_projectors(Val(2), Val(Sym))
𝟏₂ = tens_Id2(Val(2), Val(Sym))

E = symbols("E", positive = true)
ν = symbols("ν", real = true)
k = E / (3(1 - 2ν))
μ = E / (2(1 + ν))
λ = k - 2μ / 3

𝐆 = tsimplify(1 / (8 * PI * μ * (1 - ν)) * (𝐞ʳ ⊗ 𝐞ʳ - (3 - 4ν) * log(r) * 𝟏₂))

HG = -tsimplify(HESS(𝐆))
aHG = get_array(HG)
𝕄 = SymmetricTensor{4, 2}((i, j, k, l) -> (aHG[i, k, j, l] + aHG[j, k, i, l] + aHG[i, l, j, k] + aHG[j, l, i, k]) / 4)
ℾ = tsimplify(Tens(𝕄, ℬᵖ))

ℾ₂ = tsimplify(
    1 / (8PI * μ * (1 - ν) * r^2) * (
        -2𝕁₂ + 2(1 - 2ν) * 𝕀₂ + 2(𝟏₂ ⊗ 𝐞ʳ ⊗ 𝐞ʳ + 𝐞ʳ ⊗ 𝐞ʳ ⊗ 𝟏₂)
            + 8ν * 𝐞ʳ ⊗ˢ 𝟏₂ ⊗ˢ 𝐞ʳ - 8𝐞ʳ ⊗ 𝐞ʳ ⊗ 𝐞ʳ ⊗ 𝐞ʳ
    )
)

tsimplify(ℾ - ℾ₂)

ℂ₂ = 2λ * 𝕁₂ + 2μ * 𝕀₂
𝕜 = tsimplify(ℾ ⊡ ℂ₂)
get_array(𝕜)[1, 1, 1, 1]

Spherical = coorsys_spherical()
θs, ϕs, rs = getcoords(Spherical)
𝐞ᶿ, 𝐞ᵠ, 𝐞ʳˢ = unitvec(Spherical)
ℬˢ = normalized_basis(Spherical)
@set_coorsys Spherical

𝕀, 𝕁, 𝕂 = iso_projectors(Val(3), Val(Sym))
𝟏 = tens_Id2(Val(3), Val(Sym))

𝐆₃ = 1 / (8PI * μ * (3k + 4μ) * rs) * ((3k + 7μ) * 𝟏 + (3k + μ) * 𝐞ʳˢ ⊗ 𝐞ʳˢ)
𝐆₃ᵥ = 1 / (16PI * μ * (1 - ν) * rs) * ((3 - 4ν) * 𝟏 + 𝐞ʳˢ ⊗ 𝐞ʳˢ)

tsimplify(𝐆₃ - 𝐆₃ᵥ)

HG₃ = -tsimplify(HESS(𝐆₃))
aHG₃ = get_array(HG₃)
𝕄₃ = SymmetricTensor{4, 3}((i, j, k, l) -> (aHG₃[i, k, j, l] + aHG₃[j, k, i, l] + aHG₃[i, l, j, k] + aHG₃[j, l, i, k]) / 4)
ℾ₃ = tsimplify(Tens(𝕄₃, ℬˢ))

ℾ₃ᶜ = tsimplify(
    1 / (16PI * μ * (1 - ν) * rs^3) * (
        -3𝕁 + 2(1 - 2ν) * 𝕀 + 3(𝟏 ⊗ 𝐞ʳˢ ⊗ 𝐞ʳˢ + 𝐞ʳˢ ⊗ 𝐞ʳˢ ⊗ 𝟏)
            + 12ν * 𝐞ʳˢ ⊗ˢ 𝟏 ⊗ˢ 𝐞ʳˢ - 15𝐞ʳˢ ⊗ 𝐞ʳˢ ⊗ 𝐞ʳˢ ⊗ 𝐞ʳˢ
    )
)

tsimplify(ℾ₃ - ℾ₃ᶜ)

Cartesian = coorsys_cartesian(symbols("x y z", real = true))
𝐞₁, 𝐞₂, 𝐞₃ = unitvec(Cartesian)
F = symbols("F", real = true)

J = tsimplify(det(𝟏 + F * GRAD(𝐆₃ ⋅ 𝐞₁)))
factor(tsimplify(subs(J, θs => PI / 2, ϕs => 0)))

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
