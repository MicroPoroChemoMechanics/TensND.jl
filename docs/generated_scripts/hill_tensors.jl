import Pkg
Pkg.activate(joinpath(@__DIR__, "..", "docs"); io = devnull)
# The `docs` environment declares every dependency the tutorials use
# (SymPy, Symbolics, NLopt, BenchmarkTools, ...) and resolves TensND
# itself through `[sources]`. These four lines are stripped from the
# generated page and notebook, which already run inside that project.

using TensND
using LinearAlgebra
using SymPy

Cartesian = coorsys_cartesian(symbols("x y z", real = true))
𝐞₁, 𝐞₂, 𝐞₃ = unitvec(Cartesian)

Spherical = coorsys_spherical((symbols("θ ϕ", real = true)..., symbols("ξ", positive = true)))
θ, ϕ, ξ = getcoords(Spherical)
𝐞ᶿ, 𝐞ᵠ, 𝐞ʳ = unitvec(Spherical)

𝕀, 𝕁, 𝕂 = iso_projectors(Val(3), Val(Sym))
𝟏 = tens_Id2(Val(3), Val(Sym))
𝛏 = getOM(Spherical)

λ = symbols("λ", real = true)
μ = symbols("μ", positive = true)
ℂ = 3λ * 𝕁 + 2μ * 𝕀

𝐊 = 𝛏 ⋅ ℂ ⋅ 𝛏

tsimplify(𝐊)

ℾ = 𝛏 ⊗ˢ 𝐊^(-1) ⊗ˢ 𝛏
𝚲 = tsimplify(ℂ ⊡ ℾ ⊡ ℂ)

𝚲₂ = tsimplify(
    λ^2 / (λ + 2μ) * 𝟏 ⊗ 𝟏
        + 2λ * μ / (λ + 2μ) * (𝟏 ⊗ 𝐞ʳ ⊗ 𝐞ʳ + 𝐞ʳ ⊗ 𝐞ʳ ⊗ 𝟏)
        + 4μ * (𝐞ʳ ⊗ˢ 𝟏 ⊗ˢ 𝐞ʳ - (λ + μ) / (λ + 2μ) * 𝐞ʳ ⊗ 𝐞ʳ ⊗ 𝐞ʳ ⊗ 𝐞ʳ)
)

intrinsic(tsimplify(𝚲 - 𝚲₂), Spherical)

C₁₁₁₁, C₁₁₂₂, C₁₁₃₃, C₃₃₃₃, C₂₃₂₃ = symbols("C₁₁₁₁ C₁₁₂₂ C₁₁₃₃ C₃₃₃₃ C₂₃₂₃", positive = true)
n = 𝐞₃
ℂᵗⁱ = tens_TI(C₁₁₁₁, C₁₁₂₂, C₁₁₃₃, C₃₃₃₃, C₂₃₂₃, [Sym(0), Sym(0), Sym(1)])

typeof(ℂᵗⁱ), get_ℓ(ℂᵗⁱ)

𝐊ᵃˣ = tsimplify(𝐞₃ ⋅ ℂᵗⁱ ⋅ 𝐞₃)
get_array(𝐊ᵃˣ)

𝐊ᵗ = tsimplify(𝐞₁ ⋅ ℂᵗⁱ ⋅ 𝐞₁)
get_array(𝐊ᵗ)

iso_subs = Dict(
    C₁₁₁₁ => λ + 2μ, C₃₃₃₃ => λ + 2μ,
    C₁₁₂₂ => λ, C₁₁₃₃ => λ, C₂₃₂₃ => μ,
)
(tsimplify(subs.(get_array(𝐊ᵃˣ), iso_subs...)), tsimplify(subs.(get_array(𝐊ᵗ), iso_subs...)))

ℾᵗⁱ = 𝐞₃ ⊗ˢ inv(𝐊ᵃˣ) ⊗ˢ 𝐞₃
𝚲ᵗⁱ = tsimplify(ℂᵗⁱ ⊡ ℾᵗⁱ ⊡ ℂᵗⁱ)

get_array(𝚲ᵗⁱ)[3, 3, 3, 3]

Λᵗⁱ_iso = tsimplify(subs.(get_array(𝚲ᵗⁱ), iso_subs...))
Λ_iso_axis = tsimplify(subs.(get_array(𝚲₂), θ => Sym(0), ϕ => Sym(0)))

tsimplify(Λᵗⁱ_iso - Λ_iso_axis)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
