# Kelvin fundamental solution

The **Kelvin fundamental solution** (Green's function) of linear isotropic elasticity is the
second-order tensor field $\mathbf{G}(\mathbf{x})$ such that the displacement produced at $\mathbf{x}$
by a unit point force $\mathbf{f} = \mathbf{e}_i$ applied at the origin is $u_j = G_{ji}(\mathbf{x})$.
It satisfies the distributional equation

```math
\operatorname{div}(\mathbb{C}:\nabla^s\mathbf{G}) + \delta(\mathbf{x})\,\mathbf{1} = \mathbf{0},
```

where $\mathbb{C} = \lambda\,\mathbf{1}\otimes\mathbf{1} + 2\mu\,\mathbb{I}$ is the elastic stiffness,
$\delta$ the Dirac distribution, and $\mathbb{I}$ the symmetric fourth-order identity tensor.

The **fourth-order Green operator** $\boldsymbol{\Gamma}$ (Mindlin–Somigliana tensor) is the
opposite of the Hessian of $\mathbf{G}$:

```math
\boldsymbol{\Gamma} = -\nabla\nabla\mathbf{G},
\qquad
\Gamma_{ijkl} = -\partial_k\partial_l G_{ij}.
```

It plays a central role in micromechanics (Eshelby theory, polarization tensors) and in the
iterative solution of heterogeneous elasticity problems.

## 2D — Plane strain (polar coordinates)

For **plane strain** the Kelvin solution reads

```math
\mathbf{G} = \frac{1}{8\pi\mu(1-\nu)}
\Bigl(\mathbf{e}_r\otimes\mathbf{e}_r - (3-4\nu)\ln r\,\mathbf{1}\Bigr),
```

where $r = \|\mathbf{x}\|$, $\mathbf{e}_r = \mathbf{x}/r$ is the radial unit vector, $\mu$ the shear
modulus, and $\nu$ Poisson's ratio.

### Setup

```@example gf2d
using TensND, LinearAlgebra, SymPy, Tensors, OMEinsum, Rotations

Polar = coorsys_polar()
r, θ = getcoords(Polar)
𝐞ʳ, 𝐞ᶿ = unitvec(Polar)
@set_coorsys Polar
ℬˢ = normalized_basis(Polar)
𝕀, 𝕁, 𝕂 = ISO(Val(2), Val(Sym))
𝟏 = tens_Id2(Val(2), Val(Sym))

E, k, μ = symbols("E k μ", positive = true)
ν, κ = symbols("ν κ", real = true)
k = E / (3(1-2ν))
μ = E / (2(1+ν))
λ = k - 2μ/3
nothing  # hide
```

### Green's function

The second-order tensor $\mathbf{G}$ is assembled in polar coordinates using the normalized basis
and simplified symbolically:

```@example gf2d
𝐆 = tsimplify(1/(8 * PI * μ * (1-ν)) * (𝐞ʳ ⊗ 𝐞ʳ - (3-4ν) * log(r) * 𝟏))
```

### Fourth-order Green operator

The Green operator $\boldsymbol{\Gamma}$ is obtained from the Hessian of $\mathbf{G}$ (see
[`HESS`](@ref)):

```@example gf2d
HG = -tsimplify(HESS(𝐆))
aHG = get_array(HG)
𝕄 = SymmetricTensor{4,2}((i,j,k,l) -> (aHG[i,k,j,l] + aHG[j,k,i,l] + aHG[i,l,j,k] + aHG[j,l,i,k]) / 4)
ℾ = tsimplify(Tens(𝕄, ℬˢ))
```

The known closed-form expression is

```math
\boldsymbol{\Gamma} = \frac{1}{8\pi\mu(1-\nu)r^2}
\Bigl(
  -2\mathbb{J}
  + 2(1-2\nu)\mathbb{I}
  + 2\bigl(\mathbf{1}\otimes\mathbf{e}_r\otimes\mathbf{e}_r
           + \mathbf{e}_r\otimes\mathbf{e}_r\otimes\mathbf{1}\bigr)
  + 8\nu\,\mathbf{e}_r\overset{s}{\otimes}\mathbf{1}\overset{s}{\otimes}\mathbf{e}_r
  - 8\,\mathbf{e}_r\otimes\mathbf{e}_r\otimes\mathbf{e}_r\otimes\mathbf{e}_r
\Bigr),
```

where $\mathbb{J} = \tfrac{1}{2}\mathbf{1}\otimes\mathbf{1}$ (in 2D) is the spherical projector and
$\mathbf{a}\overset{s}{\otimes}\mathbf{B}\overset{s}{\otimes}\mathbf{c}$ denotes the symmetrized outer
product of a vector, a second-order tensor, and a vector. Symbolic verification gives zero:

```@example gf2d
ℾ₂ = tsimplify(1/(8PI * μ * (1-ν) * r^2) *
    (-2𝕁 + 2(1-2ν)*𝕀 + 2(𝟏⊗𝐞ʳ⊗𝐞ʳ + 𝐞ʳ⊗𝐞ʳ⊗𝟏) + 8ν*𝐞ʳ⊗ˢ𝟏⊗ˢ𝐞ʳ - 8𝐞ʳ⊗𝐞ʳ⊗𝐞ʳ⊗𝐞ʳ))
tsimplify(ℾ - ℾ₂)   # → zero tensor
```

### Contraction with the elastic stiffness

The double contraction $\boldsymbol{\Gamma}:\mathbb{C}$ (where
$\mathbb{C} = 2\lambda\mathbb{J} + 2\mu\mathbb{I}$) appears in the definition of the
polarization tensor in micromechanics:

```@example gf2d
ℂ = 2λ * 𝕁 + 2μ * 𝕀
𝕜 = tsimplify(ℾ ⊡ ℂ)
```

The result can be expressed in Cartesian coordinates via the substitution dictionary
(using the Kolosov constant $\kappa = 3 - 4\nu$ for plane strain):

```@example gf2d
Cartesian = coorsys_cartesian(symbols("x y", real = true))
x₁, x₂ = getcoords(Cartesian)
d = Dict(r      => sqrt(x₁^2 + x₂^2),
         sin(θ) => x₂ / sqrt(x₁^2 + x₂^2),
         cos(θ) => x₁ / sqrt(x₁^2 + x₂^2),
         ν      => (3 - κ) / 4)
nothing  # hide
```

## 3D — Full space (spherical coordinates)

In three dimensions the Kelvin solution admits two equivalent forms. In terms of bulk modulus
$\kappa = k$ and shear modulus $\mu$:

```math
\mathbf{G} = \frac{1}{8\pi\mu(3\kappa+4\mu)\,r}
\Bigl((3\kappa+7\mu)\,\mathbf{1} + (3\kappa+\mu)\,\mathbf{e}_r\otimes\mathbf{e}_r\Bigr),
```

and equivalently in terms of Poisson's ratio $\nu$:

```math
\mathbf{G} = \frac{1}{16\pi\mu(1-\nu)\,r}
\Bigl((3-4\nu)\,\mathbf{1} + \mathbf{e}_r\otimes\mathbf{e}_r\Bigr).
```

### Setup

```@example gf3d
using TensND, LinearAlgebra, SymPy, Tensors, OMEinsum, Rotations

Spherical = coorsys_spherical()
θ, ϕ, r = getcoords(Spherical)
𝐞ᶿ, 𝐞ᵠ, 𝐞ʳ = unitvec(Spherical)
ℬˢ = normalized_basis(Spherical)
@set_coorsys Spherical
𝕀, 𝕁, 𝕂 = ISO(Val(3), Val(Sym))
𝟏 = tens_Id2(Val(3), Val(Sym))

E, k, μ = symbols("E k μ", positive = true)
ν = symbols("ν", real = true)
k = E / (3(1-2ν))
μ = E / (2(1+ν))
λ = k - 2μ/3
nothing  # hide
```

### Green's function and equivalence of forms

Both forms are constructed and their equivalence checked:

```@example gf3d
𝐆  = 1/(8PI * μ * (3k+4μ) * r) * ((3k+7μ) * 𝟏 + (3k+μ) * 𝐞ʳ⊗𝐞ʳ)
𝐆₂ = 1/(16PI * μ * (1-ν) * r) * ((3-4ν) * 𝟏 + 𝐞ʳ⊗𝐞ʳ)
tsimplify(𝐆 - 𝐆₂)   # → zero tensor
```

### Fourth-order Green operator

The Green operator is computed from the Hessian of $\mathbf{G}$:

```@example gf3d
HG = -tsimplify(HESS(𝐆))
aHG = get_array(HG)
𝕄 = SymmetricTensor{4,3}((i,j,k,l) -> (aHG[i,k,j,l] + aHG[j,k,i,l] + aHG[i,l,j,k] + aHG[j,l,i,k]) / 4)
ℾ = tsimplify(Tens(𝕄, ℬˢ))
```

The known closed-form is

```math
\boldsymbol{\Gamma} = \frac{1}{16\pi\mu(1-\nu)r^3}
\Bigl(
  -3\mathbb{J}
  + 2(1-2\nu)\mathbb{I}
  + 3\bigl(\mathbf{1}\otimes\mathbf{e}_r\otimes\mathbf{e}_r
           + \mathbf{e}_r\otimes\mathbf{e}_r\otimes\mathbf{1}\bigr)
  + 12\nu\,\mathbf{e}_r\overset{s}{\otimes}\mathbf{1}\overset{s}{\otimes}\mathbf{e}_r
  - 15\,\mathbf{e}_r\otimes\mathbf{e}_r\otimes\mathbf{e}_r\otimes\mathbf{e}_r
\Bigr),
```

where $\mathbb{J} = \tfrac{1}{3}\mathbf{1}\otimes\mathbf{1}$ (in 3D). Symbolic verification:

```@example gf3d
ℾ₂ = tsimplify(1/(16PI * μ * (1-ν) * r^3) *
    (-3𝕁 + 2(1-2ν)*𝕀 + 3(𝟏⊗𝐞ʳ⊗𝐞ʳ + 𝐞ʳ⊗𝐞ʳ⊗𝟏) + 12ν*𝐞ʳ⊗ˢ𝟏⊗ˢ𝐞ʳ - 15𝐞ʳ⊗𝐞ʳ⊗𝐞ʳ⊗𝐞ʳ))
tsimplify(ℾ - ℾ₂)   # → zero tensor
```

### Jacobian of the displacement gradient

For a concentrated force $F\,\mathbf{e}_1$ applied at the origin, the displacement field is
$\mathbf{u} = F\,\mathbf{G}\cdot\mathbf{e}_1$. The determinant of the deformation gradient
$\mathbf{1} + \nabla\mathbf{u}$ (Jacobian) quantifies local volume change. It can be computed
symbolically in spherical coordinates and then evaluated at a specific direction:

```@example gf3d
Cartesian = coorsys_cartesian(symbols("x y z", real = true))
𝐞₁, 𝐞₂, 𝐞₃ = unitvec(Cartesian)
F = symbols("F", real = true)
J = tsimplify(det(𝟏 + F * GRAD(𝐆 ⋅ 𝐞₁)))
factor(tsimplify(subs(J, θ => PI/2, ϕ => 0)))   # on the x-axis (θ=π/2, ϕ=0)
```
