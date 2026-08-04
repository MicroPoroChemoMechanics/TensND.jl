```@meta
EditURL = "../../../../scripts/10_symbolic_operators.jl"
```

# Differential operators, symbolically

The five operators of [Curvilinear differential calculus](@ref th-curvilinear)
— `GRAD`, `SYMGRAD`, `DIV`, `LAPLACE`, `HESS` — applied in polar, cylindrical
and spherical coordinates, with exact SymPy derivatives. The results are
closed-form expressions, valid everywhere, not values at a point.

The definitions, once:

```math
\mathrm{GRAD}(t)=\sum_i\partial_i t\otimes\underline{a}^i,
\qquad
\mathrm{SYMGRAD}(t)=\sum_i\partial_i t\stackrel{s}{\otimes}\underline{a}^i,
\qquad
\mathrm{DIV}(t)=\sum_i\partial_i t\cdot\underline{a}^i .
```

````@example symbolic_operators
using TensND
using LinearAlgebra
using SymPy
````

## Setting a default coordinate system

[`@set_coorsys`](@ref) makes a system the default, so the operators can be
called with one argument.

````@example symbolic_operators
Polar = coorsys_polar()
r, θ = getcoords(Polar)
𝐞ʳ, 𝐞ᶿ = unitvec(Polar)
@set_coorsys Polar
````

## The polar Laplacian

Applied to an arbitrary function, `LAPLACE` reproduces the textbook operator
``\partial_{rr}f+\tfrac{1}{r}\partial_r f+\tfrac{1}{r^2}\partial_{\theta\theta}f``:

````@example symbolic_operators
f = SymFunction("f", real = true)(r, θ)
LAPLACE(f)
````

Harmonic functions in the plane are ``r^n\cos n\theta`` and ``r^n\sin n\theta``;
their Laplacian vanishes identically for **symbolic** ``n``:

````@example symbolic_operators
n = symbols("n", integer = true)
simplify(LAPLACE(r^n * cos(n * θ)))
````

while a mismatched pair of exponents does not vanish, and the residual is the
familiar ``(n^2-m^2)r^{n-2}``:

````@example symbolic_operators
m = symbols("m", integer = true)
simplify(LAPLACE(r^n * sin(m * θ)) / (r^(n - 2) * sin(m * θ)))
````

## The Hessian

``\mathrm{HESS}=\mathrm{GRAD}\circ\mathrm{GRAD}``, and its trace is the
Laplacian:

````@example symbolic_operators
H = simplify(HESS(r^n))
````

````@example symbolic_operators
simplify(tr(H) - LAPLACE(r^n))
````

## Cylindrical coordinates: an axisymmetric strain tensor

`SYMGRAD` applied to a displacement field is directly the linearized strain
tensor. For an axisymmetric field ``\underline{\xi}(r,z)=\xi^r\underline{e}^r
+\xi^z\underline{e}^z``, the hoop strain ``\varepsilon_{\theta\theta}=\xi^r/r``
appears from the Christoffel terms alone — there is no ``\theta`` derivative
anywhere in the field.

````@example symbolic_operators
Cylindrical = coorsys_cylindrical()
rc, θc, zc = getcoords(Cylindrical)
𝐞ʳᶜ, 𝐞ᶿᶜ, 𝐞ᶻᶜ = unitvec(Cylindrical)
@set_coorsys Cylindrical

ξʳ = SymFunction("ξʳ", real = true)(rc, zc)
ξᶻ = SymFunction("ξᶻ", real = true)(rc, zc)
𝛏 = ξʳ * 𝐞ʳᶜ + ξᶻ * 𝐞ᶻᶜ

𝛆 = tsimplify(SYMGRAD(𝛏))
````

Read off the hoop component:

````@example symbolic_operators
get_array(𝛆)[2, 2]
````

## Spherical coordinates

Recall the unusual ordering ``(\theta,\varphi,r)``, chosen so that
``\theta=\varphi=0`` gives the canonical basis in the canonical order.

````@example symbolic_operators
Spherical = coorsys_spherical()
θs, ϕs, rs = getcoords(Spherical)
𝐞ᶿ, 𝐞ᵠ, 𝐞ʳˢ = unitvec(Spherical)
@set_coorsys Spherical

getcoords(Spherical), Lame(Spherical)
````

The gradient of the radius is the radial unit vector, and the divergence of
that vector is ``2/r``:

````@example symbolic_operators
tsimplify(GRAD(rs)), tsimplify(DIV(𝐞ʳˢ))
````

The gradient of ``\underline{e}^r`` is the transverse projector divided by
``r`` — the relation that makes a sphere's curvature ``-1/r``, see
[Submanifolds](@ref th-submanifolds):

````@example symbolic_operators
tsimplify(GRAD(𝐞ʳˢ) - (𝐞ᶿ ⊗ 𝐞ᶿ + 𝐞ᵠ ⊗ 𝐞ᵠ) / rs)
````

``1/r`` is the kernel of the three-dimensional Laplace equation:

````@example symbolic_operators
tsimplify(LAPLACE(1 / rs))
````

## Radial equilibrium of a spherical stress state

A diagonal stress field in the spherical frame. Its divergence gives the
equilibrium equation of a spherically symmetric problem — the ``2(\sigma^{rr}
-\sigma^{\theta\theta})/r`` term coming entirely from the connection.

````@example symbolic_operators
σʳʳ = SymFunction("σʳʳ", real = true)(rs)
σᶿᶿ = SymFunction("σᶿᶿ", real = true)(rs)
σᵠᵠ = SymFunction("σᵠᵠ", real = true)(rs)
𝛔 = σʳʳ * 𝐞ʳˢ ⊗ 𝐞ʳˢ + σᶿᶿ * 𝐞ᶿ ⊗ 𝐞ᶿ + σᵠᵠ * 𝐞ᵠ ⊗ 𝐞ᵠ

simplify(DIV(𝛔))
````

## Index placement is a convention

`GRAD` appends the derivative index **on the right**, so
``(\nabla\underline{v})_{ij}=\partial_j v_i``, and `DIV` contracts the **last**
index, ``(\mathrm{DIV}\,\boldsymbol{\sigma})_i=\partial_j\sigma_{ij}``. In
Cartesian coordinates, where all Christoffel symbols vanish, this is easy to
read off:

````@example symbolic_operators
Cartesian = coorsys_cartesian()
X = getcoords(Cartesian)
E = unitvec(Cartesian)
@set_coorsys Cartesian

v = sum(SymFunction("v$i", real = true)(X...) * E[i] for i in 1:3)
Gv = get_array(GRAD(v))

(Gv[1, 2], diff(SymFunction("v1", real = true)(X...), X[2]))
````

The ``(1,2)`` entry is ``\partial_y v_1``, not ``\partial_x v_2``. A library
using the opposite convention differs by a transpose.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

