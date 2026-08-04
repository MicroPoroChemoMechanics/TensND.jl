import Pkg
Pkg.activate(joinpath(@__DIR__, "..", "docs"); io = devnull)
# The `docs` environment declares every dependency the tutorials use
# (SymPy, Symbolics, NLopt, BenchmarkTools, ...) and resolves TensND
# itself through `[sources]`. These four lines are stripped from the
# generated page and notebook, which already run inside that project.

using TensND
using LinearAlgebra
using SymPy

Polar = coorsys_polar()
r, θ = getcoords(Polar)
𝐞ʳ, 𝐞ᶿ = unitvec(Polar)
@set_coorsys Polar

f = SymFunction("f", real = true)(r, θ)
LAPLACE(f)

n = symbols("n", integer = true)
simplify(LAPLACE(r^n * cos(n * θ)))

m = symbols("m", integer = true)
simplify(LAPLACE(r^n * sin(m * θ)) / (r^(n - 2) * sin(m * θ)))

H = simplify(HESS(r^n))

simplify(tr(H) - LAPLACE(r^n))

Cylindrical = coorsys_cylindrical()
rc, θc, zc = getcoords(Cylindrical)
𝐞ʳᶜ, 𝐞ᶿᶜ, 𝐞ᶻᶜ = unitvec(Cylindrical)
@set_coorsys Cylindrical

ξʳ = SymFunction("ξʳ", real = true)(rc, zc)
ξᶻ = SymFunction("ξᶻ", real = true)(rc, zc)
𝛏 = ξʳ * 𝐞ʳᶜ + ξᶻ * 𝐞ᶻᶜ

𝛆 = tsimplify(SYMGRAD(𝛏))

get_array(𝛆)[2, 2]

Spherical = coorsys_spherical()
θs, ϕs, rs = getcoords(Spherical)
𝐞ᶿ, 𝐞ᵠ, 𝐞ʳˢ = unitvec(Spherical)
@set_coorsys Spherical

getcoords(Spherical), Lame(Spherical)

tsimplify(GRAD(rs)), tsimplify(DIV(𝐞ʳˢ))

tsimplify(GRAD(𝐞ʳˢ) - (𝐞ᶿ ⊗ 𝐞ᶿ + 𝐞ᵠ ⊗ 𝐞ᵠ) / rs)

tsimplify(LAPLACE(1 / rs))

σʳʳ = SymFunction("σʳʳ", real = true)(rs)
σᶿᶿ = SymFunction("σᶿᶿ", real = true)(rs)
σᵠᵠ = SymFunction("σᵠᵠ", real = true)(rs)
𝛔 = σʳʳ * 𝐞ʳˢ ⊗ 𝐞ʳˢ + σᶿᶿ * 𝐞ᶿ ⊗ 𝐞ᶿ + σᵠᵠ * 𝐞ᵠ ⊗ 𝐞ᵠ

simplify(DIV(𝛔))

Cartesian = coorsys_cartesian()
X = getcoords(Cartesian)
E = unitvec(Cartesian)
@set_coorsys Cartesian

v = sum(SymFunction("v$i", real = true)(X...) * E[i] for i in 1:3)
Gv = get_array(GRAD(v))

(Gv[1, 2], diff(SymFunction("v1", real = true)(X...), X[2]))

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
