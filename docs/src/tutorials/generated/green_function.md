```@meta
EditURL = "../../../../scripts/20_green_function.jl"
```

# Elastic Green's functions

The Kelvin fundamental solution — the displacement produced by a point force
in an infinite isotropic medium — in plane strain and in three dimensions,
built symbolically and checked against its closed form. The second-gradient
object ``\mathbb{\Gamma}=-\mathrm{HESS}(\boldsymbol{G})``, symmetrized, is the
kernel every micromechanical Green operator is built from.

This exercises the differential operators of
[Curvilinear differential calculus](@ref th-curvilinear) on a genuinely
non-trivial field: a closed-form identity that only holds if every Christoffel
term is right. Background on the elastic Green operator: [mura1987](@cite).

````@example green_function
using TensND
using LinearAlgebra
using SymPy
using Tensors
````

## Plane strain (2-D)

In polar coordinates, with ``\nu`` the Poisson ratio and ``\mu`` the shear
modulus:

```math
\boldsymbol{G}=\frac{1}{8\pi\mu(1-\nu)}
\Bigl(\underline{e}^r\otimes\underline{e}^r-(3-4\nu)\log r\;\boldsymbol{1}\Bigr)
```

````@example green_function
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
````

## The Green operator ``\mathbb{\Gamma}``

``\mathbb{\Gamma}=-\mathrm{HESS}(\boldsymbol{G})``, symmetrized over both index
pairs so that it acts on symmetric strain tensors:

````@example green_function
HG = -tsimplify(HESS(𝐆))
aHG = get_array(HG)
𝕄 = SymmetricTensor{4, 2}((i, j, k, l) -> (aHG[i, k, j, l] + aHG[j, k, i, l] + aHG[i, l, j, k] + aHG[j, l, i, k]) / 4)
ℾ = tsimplify(Tens(𝕄, ℬᵖ))
````

The closed form it must reproduce:

```math
\mathbb{\Gamma}=\frac{1}{8\pi\mu(1-\nu)r^{2}}
\Bigl(-2\mathbb{J}+2(1-2\nu)\mathbb{I}
+2(\boldsymbol{1}\otimes\underline{e}^r\otimes\underline{e}^r
  +\underline{e}^r\otimes\underline{e}^r\otimes\boldsymbol{1})
+8\nu\,\underline{e}^r\stackrel{s}{\otimes}\boldsymbol{1}\stackrel{s}{\otimes}\underline{e}^r
-8\,\underline{e}^r{}^{\otimes4}\Bigr)
```

````@example green_function
ℾ₂ = tsimplify(
    1 / (8PI * μ * (1 - ν) * r^2) * (
        -2𝕁₂ + 2(1 - 2ν) * 𝕀₂ + 2(𝟏₂ ⊗ 𝐞ʳ ⊗ 𝐞ʳ + 𝐞ʳ ⊗ 𝐞ʳ ⊗ 𝟏₂)
            + 8ν * 𝐞ʳ ⊗ˢ 𝟏₂ ⊗ˢ 𝐞ʳ - 8𝐞ʳ ⊗ 𝐞ʳ ⊗ 𝐞ʳ ⊗ 𝐞ʳ
    )
)

tsimplify(ℾ - ℾ₂)
````

Identically zero: the operator route and the closed form agree.

## Contraction with the stiffness

``\mathbb{k}=\mathbb{\Gamma}:\mathbb{C}`` is the object entering the
Lippmann–Schwinger equation of micromechanics.

````@example green_function
ℂ₂ = 2λ * 𝕁₂ + 2μ * 𝕀₂
𝕜 = tsimplify(ℾ ⊡ ℂ₂)
get_array(𝕜)[1, 1, 1, 1]
````

## Three dimensions

```math
\boldsymbol{G}=\frac{1}{16\pi\mu(1-\nu)r}
\Bigl((3-4\nu)\boldsymbol{1}+\underline{e}^r\otimes\underline{e}^r\Bigr)
```

equivalently written with the bulk modulus, and the two forms must agree:

````@example green_function
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
````

The same construction in 3-D:

````@example green_function
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
````

Note the pattern: the 2-D coefficients ``(-2,\,2,\,2,\,8,\,-8)`` become
``(-3,\,2,\,3,\,12,\,-15)`` in 3-D, and ``r^{-2}`` becomes ``r^{-3}``.

## The Jacobian of the induced deformation

For a unit point force along ``\underline{e}_1``, the deformation gradient is
``\boldsymbol{1}+F\,\nabla(\boldsymbol{G}\cdot\underline{e}_1)`` and its
determinant measures the local volume change. Evaluated near the equator
``\theta=\pi/2``, ``\varphi=0``:

````@example green_function
Cartesian = coorsys_cartesian(symbols("x y z", real = true))
𝐞₁, 𝐞₂, 𝐞₃ = unitvec(Cartesian)
F = symbols("F", real = true)

J = tsimplify(det(𝟏 + F * GRAD(𝐆₃ ⋅ 𝐞₁)))
factor(tsimplify(subs(J, θs => PI / 2, ϕs => 0)))
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

