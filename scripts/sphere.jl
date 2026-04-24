# ============================================================================
#  Isotropic sphere under classical remote loadings (symbolic)
#
#  Constructs the displacement/stress fields inside a hollow isotropic
#  sphere for three canonical remote strain states:
#    1. Spherical:              𝔼 = 𝟏
#    2. Deviatoric axisymmetric: 𝔼 = 𝟏 − 3 𝐞₃⊗𝐞₃
#    3. Pure shear:              𝔼 = 𝐞₁⊗𝐞₁ − 𝐞₂⊗𝐞₂
#
#  For each case, solves the Navier equilibrium equations symbolically
#  (via SymPy `dsolve`/`solve`) and recovers the radial traction 𝐓⋅𝐞ʳ
#  plus the Lamé exponents (α, Λ) of the displacement ansatz uₐ(r) = C·rᵅ.
# ============================================================================

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)

using TensND, LinearAlgebra, SymPy, Tensors, OMEinsum, Rotations

Spherical = coorsys_spherical()
θ, ϕ, r = getcoords(Spherical)
𝐞ᶿ, 𝐞ᵠ, 𝐞ʳ = unitvec(Spherical)
ℬˢ = normalized_basis(Spherical)
𝐱 = getOM(Spherical)
@set_coorsys Spherical

𝐞₁, 𝐞₂, 𝐞₃ = unitvec(coorsys_cartesian())
𝕀, 𝕁, 𝕂 = iso_projectors(Val(3), Val(Sym))
𝟏 = tens_Id2(Val(3), Val(Sym))
k, μ = symbols("k μ", positive = true)
ℂ = 3k*𝕁 + 2μ*𝕂
remote_angle_functions(𝐄) = (fʳ=simplify(𝐞ʳ ⋅ 𝐄 ⋅ 𝐞ʳ) ; (diff(fʳ, θ)/2, diff(fʳ, ϕ)/(2sin(θ)), fʳ)) # return fᶿ, fᵠ, fʳ
u = SymFunction("u", real = true)
uᶿ = SymFunction("uᶿ", real = true)
uᵠ = SymFunction("uᵠ", real = true)
uʳ = SymFunction("uʳ", real = true)
α, Λ = symbols("α Λ", real = true)

# Spherical loading 𝔼 = 𝟏
𝐮 = u(r) * 𝐞ʳ
𝛆 = SYMGRAD(𝐮, Spherical)
𝛔 = ℂ ⊡ 𝛆
𝐓 = 𝛔 ⋅ 𝐞ʳ
div𝛔 = DIV(𝛔, Spherical)
eq = factor(simplify(div𝛔 ⋅ 𝐞ʳ))
sol = dsolve(eq, u(r))
û = sol.rhs()
T̂ = tsimplify(tsimplify(subs(𝐓 ⋅ 𝐞ʳ, u(r) => û)))

# Deviatoric axisymmetric loading 𝔼 = 𝟏 - 3𝐞₃⊗𝐞₃
fᶿ, _, fʳ = remote_angle_functions(𝟏 - 3𝐞₃⊗𝐞₃)
𝐮 = uᶿ(r) * fᶿ * 𝐞ᶿ + uʳ(r) * fʳ * 𝐞ʳ
𝛆 = SYMGRAD(𝐮, Spherical)
𝛔 = ℂ ⊡ 𝛆
𝐓 = 𝛔 ⋅ 𝐞ʳ
div𝛔 = DIV(𝛔, Spherical)
eqᶿ = tsimplify(div𝛔 ⋅ 𝐞ᶿ / fᶿ)
eqʳ = tsimplify(div𝛔 ⋅ 𝐞ʳ / fʳ)
eq = tsimplify.(subs.([eqᶿ,eqʳ], uᶿ(r) => r^α, uʳ(r) => Λ*r^α))
αΛ = solve([eqi.doit() for eqi ∈ eq], [α, Λ])
ûᶿ = sum([Sym("C$(i+2)") * r^αΛ[i][1] for i ∈ 1:length(αΛ)])
ûʳ = sum([Sym("C$(i+2)") * αΛ[i][2] * r^αΛ[i][1] for i ∈ 1:length(αΛ)])
T̂ᶿ = tsimplify(tsimplify(subs(simplify(𝐓 ⋅ 𝐞ᶿ / fᶿ), uᶿ(r) => ûᶿ, uʳ(r) => ûʳ)))
T̂ʳ = tsimplify(tsimplify(subs(simplify(𝐓 ⋅ 𝐞ʳ / fʳ), uᶿ(r) => ûᶿ, uʳ(r) => ûʳ)))

# Deviatoric pure shear loading 𝔼 = 𝐞₁⊗𝐞₁ - 𝐞₂⊗𝐞₂
fᶿ, fᵠ, fʳ = remote_angle_functions(𝐞₁⊗𝐞₁ - 𝐞₂⊗𝐞₂)
𝐮 = uᶿ(r) * fᶿ * 𝐞ᶿ + uᵠ(r) * fᵠ * 𝐞ᵠ + uʳ(r) * fʳ * 𝐞ʳ
𝛆 = SYMGRAD(𝐮, Spherical)
𝛔 = ℂ ⊡ 𝛆
𝐓 = 𝛔 ⋅ 𝐞ʳ
div𝛔 = DIV(𝛔, Spherical)
eqᶿ = tsimplify(div𝛔 ⋅ 𝐞ᶿ / fᶿ)
eqᵠ = tsimplify(div𝛔 ⋅ 𝐞ᵠ / fᵠ)
eqʳ = tsimplify(div𝛔 ⋅ 𝐞ʳ / fʳ)
X = symbols("X", real = true)
uᵠsol = solve(tsimplify(diff(subs(eqᵠ, sin(θ)^2 => 1/X), X)), uᵠ(r))[1]  # shows that uᵠ(r)=uᶿ(r)
eq = tsimplify.(subs.([eqᶿ,eqʳ], uᵠ(r) => r^α, uᶿ(r) => r^α, uʳ(r) => Λ*r^α))
αΛ = solve([eqi.doit() for eqi ∈ eq], [α, Λ])
ûᶿ = sum([Sym("C$(i+2)") * r^αΛ[i][1] for i ∈ 1:length(αΛ)])
ûʳ = sum([Sym("C$(i+2)") * αΛ[i][2] * r^αΛ[i][1] for i ∈ 1:length(αΛ)])
T̂ᶿ = tsimplify(tsimplify(subs(simplify(𝐓 ⋅ 𝐞ᶿ / fᶿ), uᶿ(r) => ûᶿ, uᵠ(r) => ûᶿ, uʳ(r) => ûʳ)))
T̂ʳ = tsimplify(tsimplify(subs(simplify(𝐓 ⋅ 𝐞ʳ / fʳ), uᶿ(r) => ûᶿ, uᵠ(r) => ûᶿ, uʳ(r) => ûʳ)))


E₁, E₂, E₃ = symbols("E₁, E₂, E₃", real = true)
𝐄 = E₁ * 𝐞₁⊗𝐞₁ + E₂ * 𝐞₂⊗𝐞₂ + E₃ * 𝐞₃⊗𝐞₃
