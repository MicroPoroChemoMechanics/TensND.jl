```@meta
EditURL = "../../../../scripts/20_green_function.jl"
```

# [Elastic Green's functions](@id tut-green-function)

!!! tip "See also"
    **Theory**: [Tensors of the plane, of space and of any dimension](@ref th-dimensions) · [Isotropic tensors](@ref th-isotropic) · [The Walpole basis](@ref th-walpole)

    **Manual**: [Two dimensions, and any dimension](@ref man-dimensions) · [Structured tensors](@ref man-structured) · [Elastic moduli](@ref man-elastic-moduli)

    **Tutorials**: [The symmetry classes of a plane section](@ref tut-plane-classes) · [The elastic sphere, solved symbolically](@ref tut-sphere-problems) · [A two-point Green tensor between offset frames](@ref tut-cluster) · [Acoustic and Hill tensors](@ref tut-hill-tensors)

The Kelvin fundamental solution — the displacement produced by a point force
in an infinite isotropic medium — in plane strain and in three dimensions,
built symbolically and checked against its closed form. The second-gradient
object ``\mathbb{\Gamma}=-\mathrm{Hess}(\boldsymbol{G})``, symmetrized, is the
kernel every micromechanical Green operator is built from.

This exercises the differential operators of
[Curvilinear differential calculus](@ref th-curvilinear) on a genuinely
non-trivial field: a closed-form identity that only holds if every Christoffel
term is right. Background on the elastic Green operator: [mura1987](@citet).

````@example green_function
using TensND
using LinearAlgebra
using SymPyPythonCall
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
r, θ = get_coords(Polar)
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

``\mathbb{\Gamma}=-\mathrm{Hess}(\boldsymbol{G})``, symmetrized over both index
pairs so that it acts on symmetric strain tensors:

````@example green_function
HG = -tsimplify(hess(𝐆))
aHG = get_array(HG)
𝕄 = SymmetricTensor{4, 2}((i, j, k, l) -> (aHG[i, k, j, l] + aHG[j, k, i, l] + aHG[i, l, j, k] + aHG[j, l, i, k]) / 4)
ℾ = tsimplify(Tens(𝕄, ℬᵖ))
````

The closed form it must reproduce:

```math
\begin{aligned}
\mathbb{\Gamma}=\frac{1}{8\pi\mu(1-\nu)r^{2}}
\Bigl(&-2\mathbb{J}+2(1-2\nu)\mathbb{I}
+2(\boldsymbol{1}\otimes\underline{e}^r\otimes\underline{e}^r
  +\underline{e}^r\otimes\underline{e}^r\otimes\boldsymbol{1})\\
&{}+8\nu\,\underline{e}^r\stackrel{s}{\otimes}\boldsymbol{1}\stackrel{s}{\otimes}\underline{e}^r
-8\,\underline{e}^r{}^{\otimes4}\Bigr)
\end{aligned}
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

Identically zero: the operator method and the closed form agree.

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
θs, ϕs, rs = get_coords(Spherical)
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
HG₃ = -tsimplify(hess(𝐆₃))
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

J = tsimplify(det(𝟏 + F * grad(𝐆₃ ⋅ 𝐞₁)))
factor(tsimplify(subs(J, θs => PI / 2, ϕs => 0)))
````

## The point force is recovered

Away from the source the displacement ``\boldsymbol{G}\cdot\underline{a}`` of
a unit force ``\underline{a}`` is in equilibrium, and the tractions on any
closed line (in the plane) or surface (in space) around the source balance the
force: their resultant is ``-\underline{a}``. In the plane, with
``\underline{a}=\underline{e}_1`` written in the polar basis and the stress of
plane strain ``\boldsymbol{\sigma}=\lambda\operatorname{tr}\boldsymbol{\varepsilon}\,\boldsymbol{1}+2\mu\,\boldsymbol{\varepsilon}``:

````@example green_function
𝐩ʳ, 𝐩ᶿ = unitvec(Polar)
𝐚 = cos(θ) * 𝐩ʳ - sin(θ) * 𝐩ᶿ
𝐛 = sin(θ) * 𝐩ʳ + cos(θ) * 𝐩ᶿ
𝛆 = tsimplify(symgrad(tsimplify(𝐆 ⋅ 𝐚), Polar))
𝛔 = tsimplify(λ * tr(𝛆) * 𝟏₂ + 2μ * 𝛆)
get_array(tsimplify(divergence(𝛔, Polar)))
````

The resultant of the tractions ``\boldsymbol{\sigma}\cdot\underline{e}^r`` on
the circle of radius ``r``, along the force and normal to it:

````@example green_function
𝐭 = 𝛔 ⋅ 𝐩ʳ
simplify(integrate(simplify((𝐭 ⋅ 𝐚) * r), (θ, 0, 2PI))), simplify(integrate(simplify((𝐭 ⋅ 𝐛) * r), (θ, 0, 2PI)))
````

In space, with ``\underline{a}=\underline{e}_3`` in the spherical basis and
the resultant on the sphere of radius ``r``:

````@example green_function
𝐚₃ = cos(θs) * 𝐞ʳˢ - sin(θs) * 𝐞ᶿ
𝛆₃ = tsimplify(symgrad(tsimplify(𝐆₃ ⋅ 𝐚₃), Spherical))
𝛔₃ = tsimplify(λ * tr(𝛆₃) * 𝟏 + 2μ * 𝛆₃)
get_array(tsimplify(divergence(𝛔₃, Spherical)))
````

````@example green_function
𝐭₃ = 𝛔₃ ⋅ 𝐞ʳˢ
simplify(integrate(simplify((𝐭₃ ⋅ 𝐚₃) * rs^2 * sin(θs)), (θs, 0, PI), (ϕs, 0, 2PI)))
````

## Conduction, in the plane and in space

The Green function of conduction is the temperature of a unit point source
in a medium of conductivity ``\kappa``: ``-\log r/(2\pi\kappa)`` in the
plane, ``1/(4\pi\kappa r)`` in space. Each is harmonic away from the source,
and the flux ``-\kappa\,\mathrm{grad}\,g`` through a circle or a sphere is the
unit source:

````@example green_function
κ = symbols("κ", positive = true)
g₂ = -log(r) / (2PI * κ)
g₃ = 1 / (4PI * κ * rs)
simplify(laplacian(g₂, Polar)), simplify(laplacian(g₃, Spherical))
````

````@example green_function
flux₂ = integrate(simplify(-κ * (grad(g₂, Polar) ⋅ 𝐩ʳ) * r), (θ, 0, 2PI))
flux₃ = integrate(simplify(-κ * (grad(g₃, Spherical) ⋅ 𝐞ʳˢ) * rs^2 * sin(θs)), (θs, 0, PI), (ϕs, 0, 2PI))
flux₂, flux₃
````

## An orthotropic plane

With an orthotropic conductivity
``\boldsymbol{K}=\kappa_1\,\underline{e}_1\otimes\underline{e}_1+\kappa_2\,\underline{e}_2\otimes\underline{e}_2``,
a tensor of the plane, the change of variables ``x/\sqrt{\kappa_1}``,
``y/\sqrt{\kappa_2}`` brings the equation back to the isotropic one, and

```math
g=-\frac{1}{4\pi\sqrt{\kappa_1\kappa_2}}\,\log\Bigl(\frac{x^2}{\kappa_1}+\frac{y^2}{\kappa_2}\Bigr) .
```

It satisfies ``\mathrm{div}(\boldsymbol{K}\cdot\mathrm{grad}\,g)=0`` away from
the source, in the Cartesian system of the plane:

````@example green_function
Plane = coorsys_cartesian(Val(2))
x, y = get_coords(Plane)
κ₁, κ₂ = symbols("κ₁ κ₂", positive = true)
𝐊 = TensOrtho{2}(κ₁, κ₂, CanonicalBasis{2, Sym}())
gₒ = -log(x^2 / κ₁ + y^2 / κ₂) / (4PI * sqrt(κ₁ * κ₂))
𝐪 = -(𝐊 ⋅ grad(gₒ, Plane))
simplify(divergence(𝐪, Plane))
````

and the flux through the unit circle is the unit source, here for
``\kappa_1=3`` and ``\kappa_2=1`` (the integrand is periodic: the trapezoidal
rule converges faster than any power of the step):

````@example green_function
τ = symbols("τ", real = true)
qₙ = simplify(subs(get_array(𝐪)[1] * x + get_array(𝐪)[2] * y, x => cos(τ), y => sin(τ), κ₁ => 3, κ₂ => 1))
f = lambdify(qₙ, (τ,))
sum(f, range(0, 2π; length = 201)[1:200]) * 2π / 200
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

