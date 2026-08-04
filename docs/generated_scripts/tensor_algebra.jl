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
using Random

Random.seed!(20260804)

a = Tens(Vec{3}((i,) -> symbols("a$i", real = true)))
b = Tens(Vec{3}((i,) -> symbols("b$i", real = true)))
a ⊗ b

a ⊗ˢ b

tsimplify(get_array(a ⊗ b) - get_array(a ⊗ˢ b))

𝟏 = Matrix(1.0I, 3, 3)
𝕀 = get_array(tens_Id4(Val(3), Val(Float64)))

norm(𝕀 - 𝟏 ⊠ˢ 𝟏)

m = rand(3, 3)
norm((𝟏 ⊠ 𝟏) ⊡ m - m), norm((𝟏 ⊠ˢ 𝟏) ⊡ m - m), norm((𝟏 ⊠ˢ 𝟏) ⊡ m - (m + m') / 2)

A, B, C, D = rand(3, 3), rand(3, 3), rand(3, 3), rand(3, 3)

ids = [
    "(a⊠b):(c⊠d) = (a·c)⊠(b·d)" => norm((A ⊠ B) ⊡ (C ⊠ D) - (A * C) ⊠ (B * D)),
    "(a⊠b):c     = a·c·ᵗb" => norm((A ⊠ B) ⊡ C - A * C * B'),
    "(a⊠ˢb):c    = (a·c·ᵗb + a·ᵗc·ᵗb)/2" => norm((A ⊠ˢ B) ⊡ C - (A * C * B' + A * C' * B') / 2),
    "(a⊗b):(c⊗d) = (b:c) a⊗d" => norm((A ⊗ B) ⊡ (C ⊗ D) - sum(B .* C) * (A ⊗ D)),
    "(a⊠b):(a⁻¹⊠b⁻¹) = 1⊠1" => norm((A ⊠ B) ⊡ (inv(A) ⊠ inv(B)) - 𝟏 ⊠ 𝟏),
]
for (name, residual) in ids
    println(rpad(name, 38), " residual = ", residual)
end

chk(x, y) = norm((x ⊠ˢ y) ⊡ (inv(x) ⊠ˢ inv(y)) - 𝟏 ⊠ˢ 𝟏)
D1, D2 = diagm(rand(3)), diagm(rand(3))     # diagonal ⟹ they commute

for (name, r) in [
        "b = a          " => chk(A, A),
        "b = 3a         " => chk(A, 3A),
        "a, b commuting " => chk(D1, D2),
        "b = 1          " => chk(A, 𝟏),
    ]
    println(name, " residual = ", round(r, sigdigits = 4))
end

𝕁 = tens_J4(Val(3), Val(Float64))
𝕂 = tens_K4(Val(3), Val(Float64))

(𝕁 ⊙ 𝕁, 𝕂 ⊙ 𝕂, 𝕁 ⊙ 𝕂)

[(d, tens_K4(Val(d), Val(Float64)) ⊙ tens_K4(Val(d), Val(Float64)), d * (d + 1) ÷ 2 - 1) for d in 2:3]

k, μ = symbols("k μ", positive = true)
𝕀s, 𝕁s, 𝕂s = iso_projectors(Val(3), Val(Sym))
ℂ = 3k * 𝕁s + 2μ * 𝕂s

(tsimplify(ℂ ⊙ 𝕁s / 3), tsimplify(ℂ ⊙ 𝕂s / 10))

t = TensOrtho(10.0, 8.0, 9.0, 3.0, 2.0, 4.0, 2.5, 3.0, 1.5, CanonicalBasis{3, Float64}())
Ct = get_array(t)
kiso = (Ct ⊙ get_array(𝕁)) / 3
μiso = (Ct ⊙ get_array(𝕂)) / 10
(kiso, μiso)

Biso, d, drel = proj_tens(:ISO, Ct)
println("closest isotropic tensor : ", Biso)
println("relative distance        : ", round(drel, sigdigits = 4))

norm(get_array(Biso) - (3kiso * get_array(𝕁) + 2μiso * get_array(𝕂)))

Siso, _, _ = proj_tens(:ISO, get_array(inv(t)))
lhs = get_array(inv(Biso))          # ISO(ℂ) then invert
rhs = get_array(Siso)               # invert then ISO

println("‖inv(ISO(ℂ)) − ISO(inv(ℂ))‖ = ", round(norm(lhs - rhs), sigdigits = 4))
println("relative                    = ", round(norm(lhs - rhs) / norm(rhs), sigdigits = 4))

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
