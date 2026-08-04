```@meta
EditURL = "../../../../scripts/22_sphere_problems.jl"
```

# The elastic sphere, solved symbolically

The classical building block of micromechanics: an isotropic sphere — or a
stack of concentric isotropic layers — under a uniform remote strain. Because
the material is isotropic, the problem splits by the symmetry of the remote
loading, and each part reduces to an ordinary differential equation in the
radius that SymPy solves in closed form.

Everything here is produced by the differential operators of
[Curvilinear differential calculus](@ref th-curvilinear): `SYMGRAD` gives the
strain, `DIV` the equilibrium equation. Nothing is transcribed from a
textbook. Background: [mura1987](@cite).

````@example sphere_problems
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
````

## The hydrostatic problem

For a remote strain ``\mathbb{E}^\infty\propto\boldsymbol{1}`` the
displacement is purely radial, ``\underline{u}=u(r)\,\underline{e}^r``.

````@example sphere_problems
u = SymFunction("u", real = true)
𝐮 = u(r) * 𝐞ʳ
````

The strain follows from `SYMGRAD`, the stress from the constitutive law, and
the equilibrium equation from `DIV`:

````@example sphere_problems
𝛆 = SYMGRAD(𝐮)
𝛔 = ℂ ⊡ 𝛆
𝐓 = 𝛔 ⋅ 𝐞ʳ
eq = factor(simplify(DIV(𝛔) ⋅ 𝐞ʳ))
````

A second-order Euler equation in ``r``, solved directly:

````@example sphere_problems
sol = dsolve(eq, u(r))
û = sol.rhs()
````

The two exponents are ``r`` and ``r^{-2}``, i.e. the familiar
``u(r)=C_1r+C_2/r^2``. The radial traction carries the interface and boundary
conditions of a layered assemblage:

````@example sphere_problems
T̂ = tsimplify(tsimplify(subs(𝐓 ⋅ 𝐞ʳ, u(r) => û)))
````

A solid sphere is regular at the origin, so the ``r^{-2}`` constant vanishes
and the remaining term is ``3k`` times the uniform strain — a state of uniform
hydrostatic stress, as it must be.

## The deviatoric axisymmetric problem

For ``\mathbb{E}^\infty=\boldsymbol{1}-3\,\underline{e}_3\otimes\underline{e}_3``
the angular dependence is fixed by the loading and only the radial profiles
remain unknown. The angular functions are generated from the remote strain
itself:

````@example sphere_problems
remote_angle_functions(𝐄) = let fʳ = simplify(𝐞ʳ ⋅ 𝐄 ⋅ 𝐞ʳ)
    (diff(fʳ, θ) / 2, diff(fʳ, ϕ) / (2sin(θ)), fʳ)
end

uᶿ = SymFunction("uᶿ", real = true)
uᵠ = SymFunction("uᵠ", real = true)
uʳ = SymFunction("uʳ", real = true)
α, Λ = symbols("α Λ", real = true)

fᶿ, _, fʳ = remote_angle_functions(𝟏 - 3𝐞₃ ⊗ 𝐞₃)
(fᶿ, fʳ)
````

````@example sphere_problems
𝐮ᵈ = uᶿ(r) * fᶿ * 𝐞ᶿ + uʳ(r) * fʳ * 𝐞ʳ
𝛔ᵈ = ℂ ⊡ SYMGRAD(𝐮ᵈ)
𝐓ᵈ = 𝛔ᵈ ⋅ 𝐞ʳ
div𝛔ᵈ = DIV(𝛔ᵈ)
````

The angular dependence factors out exactly: dividing by ``f^\theta`` and
``f^r`` leaves two coupled radial equations.

````@example sphere_problems
eqᶿ = tsimplify(div𝛔ᵈ ⋅ 𝐞ᶿ / fᶿ)
eqʳ = tsimplify(div𝛔ᵈ ⋅ 𝐞ʳ / fʳ)
````

## The Lamé exponents

Substituting the power ansatz ``u^\theta=r^\alpha``, ``u^r=\Lambda r^\alpha``
turns the differential system into an algebraic one for ``(\alpha,\Lambda)``:

````@example sphere_problems
eqs = tsimplify.(subs.([eqᶿ, eqʳ], uᶿ(r) => r^α, uʳ(r) => Λ * r^α))
αΛ = solve([e.doit() for e in eqs], [α, Λ])
````

Four solutions — the four exponents of the deviatoric problem. Two are regular
at the origin and two at infinity, which is exactly what a layered assemblage
needs: two constants per layer, fixed by continuity at each interface.

````@example sphere_problems
[(pair[1], simplify(pair[2])) for pair in αΛ]
````

The general solution is their combination:

````@example sphere_problems
ûᶿ = sum(Sym("C$(i + 2)") * r^αΛ[i][1] for i in 1:length(αΛ))
ûʳ = sum(Sym("C$(i + 2)") * αΛ[i][2] * r^αΛ[i][1] for i in 1:length(αΛ))
(ûᶿ, ûʳ)
````

and the tractions that carry the interface conditions:

````@example sphere_problems
T̂ᶿ = tsimplify(tsimplify(subs(simplify(𝐓ᵈ ⋅ 𝐞ᶿ / fᶿ), uᶿ(r) => ûᶿ, uʳ(r) => ûʳ)))
````

````@example sphere_problems
T̂ʳ = tsimplify(tsimplify(subs(simplify(𝐓ᵈ ⋅ 𝐞ʳ / fʳ), uᶿ(r) => ûᶿ, uʳ(r) => ûʳ)))
````

## Pure shear gives the same exponents

A different deviatoric loading,
``\mathbb{E}^\infty=\underline{e}_1\otimes\underline{e}_1
-\underline{e}_2\otimes\underline{e}_2``, now with an azimuthal component.
Isotropy demands that it produce the *same* radial exponents; only the angular
functions differ.

````@example sphere_problems
fᶿ₂, fᵠ₂, fʳ₂ = remote_angle_functions(𝐞₁ ⊗ 𝐞₁ - 𝐞₂ ⊗ 𝐞₂)
(fᶿ₂, fᵠ₂, fʳ₂)
````

````@example sphere_problems
𝐮ˢ = uᶿ(r) * fᶿ₂ * 𝐞ᶿ + uᵠ(r) * fᵠ₂ * 𝐞ᵠ + uʳ(r) * fʳ₂ * 𝐞ʳ
div𝛔ˢ = DIV(ℂ ⊡ SYMGRAD(𝐮ˢ))

eqᶿˢ = tsimplify(div𝛔ˢ ⋅ 𝐞ᶿ / fᶿ₂)
eqᵠˢ = tsimplify(div𝛔ˢ ⋅ 𝐞ᵠ / fᵠ₂)
eqʳˢ = tsimplify(div𝛔ˢ ⋅ 𝐞ʳ / fʳ₂)
````

The azimuthal equation forces ``u^\varphi=u^\theta``: the two transverse
profiles are not independent.

````@example sphere_problems
X = symbols("X", real = true)
uᵠsol = solve(tsimplify(diff(subs(eqᵠˢ, sin(θ)^2 => 1 / X), X)), uᵠ(r))[1]
````

With that identification the remaining system reproduces the same exponents as
the axisymmetric case:

````@example sphere_problems
eqs₂ = tsimplify.(subs.([eqᶿˢ, eqʳˢ], uᵠ(r) => r^α, uᶿ(r) => r^α, uʳ(r) => Λ * r^α))
αΛ₂ = solve([e.doit() for e in eqs₂], [α, Λ])

sort(string.(first.(αΛ))) == sort(string.(first.(αΛ₂)))
````

The deviatoric response of an isotropic sphere therefore depends on the
*symmetry class* of the remote loading, not on its particular orientation —
which is what lets an ``N``-layer assemblage be solved once and reused for any
deviatoric loading.

## Assembling ``N`` layers

In a stack of concentric layers, layer ``i`` occupies
``R_{i-1}\le r\le R_i`` with its own moduli ``(k_i,\mu_i)`` and carries the
constants found above. Continuity of the displacement and of the radial
traction at each interface gives two scalar equations per interface in the
hydrostatic problem and four in the deviatoric one; regularity at the center
and the remote condition close the system.

The expressions `T̂`, `T̂ᶿ` and `T̂ʳ` derived above are exactly the quantities
to be matched, so the assembly is ordinary linear algebra on the constants —
no further tensor calculus is required.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

