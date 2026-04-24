# N-layer sphere

The **N-layer sphere** (or Composite Sphere Assemblage) is a classical micromechanics model used to
estimate effective elastic properties of composites. The sphere consists of $N$ concentric isotropic
layers: layer $i$ occupies $R_{i-1} \le r \le R_i$ ($R_0 = 0$) with Lamé constants $(\lambda_i, \mu_i)$, or equivalently bulk modulus $\kappa_i$ and shear modulus $\mu_i$.

In each layer the equilibrium equation reads

```math
\operatorname{div}\boldsymbol{\sigma} = \mathbf{0},
\qquad
\boldsymbol{\sigma} = \lambda\operatorname{tr}(\boldsymbol{\varepsilon})\,\mathbf{1} + 2\mu\,\boldsymbol{\varepsilon},
\qquad \lambda = \kappa - \tfrac{2\mu}{3}.
```

Under a uniform far-field strain $\mathbf{E}^\infty$, the displacement field decomposes into two independent problems according to the isotropic symmetry of the layers:

- **Spherical** (hydrostatic) part — $\mathbf{E}^\infty \propto \mathbf{1}$
- **Deviatoric** part — $\mathbf{E}^\infty = \mathbf{1} - 3\mathbf{e}_3\otimes\mathbf{e}_3$

## Setup

```@example nlayer
using TensND, LinearAlgebra, SymPy, Tensors, OMEinsum, Rotations

Spherical = coorsys_spherical()
θ, ϕ, r = getcoords(Spherical)
𝐞ᶿ, 𝐞ᵠ, 𝐞ʳ = unitvec(Spherical)
ℬˢ = normalized_basis(Spherical)
𝐱 = getOM(Spherical)
𝐞₁, 𝐞₂, 𝐞₃ = unitvec(coorsys_cartesian())
𝕀, 𝕁, 𝕂 = ISO(Val(3),Val(Sym))
𝟏 = tens_Id2(Val(3),Val(Sym))
k, μ = symbols("k μ", positive = true)
λ = k - 2μ/3
nothing  # hide
```

## Spherical (hydrostatic) problem

For a hydrostatic far-field strain $\mathbf{E}^\infty \propto \mathbf{1}$, by symmetry the displacement is purely radial:

```math
\mathbf{u}^{sph} = u(r)\,\mathbf{e}_r.
```

The equilibrium equation $\operatorname{div}\boldsymbol{\sigma} = \mathbf{0}$ reduces to a 2nd-order ODE for $u(r)$
whose general solution is

```math
u(r) = C_1\,r + \frac{C_2}{r^2}.
```

The radial traction $\hat{T}^{sph}(r) = (\boldsymbol{\sigma}\cdot\mathbf{e}_r)\cdot\mathbf{e}_r$ provides the
interface and boundary conditions.

```@example nlayer
u = SymFunction("u", real = true)
𝐮ˢᵖʰ = u(r) * 𝐞ʳ
𝛆ˢᵖʰ = SYMGRAD(𝐮ˢᵖʰ, Spherical)
𝛔ˢᵖʰ = λ * tr(𝛆ˢᵖʰ) * 𝟏 + 2μ * 𝛆ˢᵖʰ
𝐓ˢᵖʰ = 𝛔ˢᵖʰ ⋅ 𝐞ʳ
div𝛔ˢᵖʰ = DIV(𝛔ˢᵖʰ, Spherical) ;
eqˢᵖʰ = factor(simplify(div𝛔ˢᵖʰ ⋅ 𝐞ʳ))
```

```@example nlayer
solˢᵖʰ = dsolve(eqˢᵖʰ, u(r))
ûˢᵖʰ = solˢᵖʰ.rhs()
T̂ˢᵖʰ = factor(simplify(subs(𝐓ˢᵖʰ ⋅ 𝐞ʳ, u(r) => ûˢᵖʰ)))
```

## Deviatoric problem

For the deviatoric loading defined by

```math
\mathbf{E} = \mathbf{1} - 3\,\mathbf{e}_3\otimes\mathbf{e}_3,
```

the angular dependence factorizes and the displacement takes the form

```math
\mathbf{u}^{dev} = u^\theta(r)\,f^\theta\,\mathbf{e}_\theta + u^r(r)\,f^r\,\mathbf{e}_r,
```

where the angular factors are scalar projections of $\mathbf{E}$:

```math
f^\theta = \mathbf{e}_\theta \cdot \mathbf{E} \cdot \mathbf{e}_r,
\qquad
f^r = \mathbf{e}_r \cdot \mathbf{E} \cdot \mathbf{e}_r.
```

The equilibrium equations reduce to a $2\times 2$ ODE system for $(u^\theta(r),u^r(r))$.
The power law ansatz $u^\theta = r^\alpha$, $u^r = \Lambda r^\alpha$ yields four solutions $(\alpha_i,\Lambda_i)$:

```math
u^\theta(r) = \sum_{i=1}^{4} C_{i+2}\,r^{\alpha_i},
\qquad
u^r(r) = \sum_{i=1}^{4} C_{i+2}\,\Lambda_i\,r^{\alpha_i}.
```

The tangential and radial tractions $\hat{T}^\theta$, $\hat{T}^r$ (divided by $f^\theta$, $f^r$ respectively)
are then expressed in terms of the four constants $C_3,\ldots,C_6$.

```@example nlayer
𝐄 = 𝟏 - 3𝐞₃⊗𝐞₃
fᶿ = simplify(𝐞ᶿ ⋅ 𝐄 ⋅ 𝐞ʳ)
fʳ = simplify(𝐞ʳ ⋅ 𝐄 ⋅ 𝐞ʳ)
uᶿ = SymFunction("uᶿ", real = true)
uʳ = SymFunction("uʳ", real = true)
𝐮ᵈᵉᵛ = uᶿ(r) * fᶿ * 𝐞ᶿ + uʳ(r) * fʳ * 𝐞ʳ
𝛆ᵈᵉᵛ = SYMGRAD(𝐮ᵈᵉᵛ, Spherical)
𝛔ᵈᵉᵛ = λ * tr(𝛆ᵈᵉᵛ) * 𝟏 + 2μ * 𝛆ᵈᵉᵛ
𝐓ᵈᵉᵛ = 𝛔ᵈᵉᵛ ⋅ 𝐞ʳ
div𝛔ᵈᵉᵛ = simplify(DIV(𝛔ᵈᵉᵛ, Spherical))
eqᶿᵈᵉᵛ = factor(simplify(div𝛔ᵈᵉᵛ ⋅ 𝐞ᶿ / fᶿ))
eqʳᵈᵉᵛ = factor(simplify(div𝛔ᵈᵉᵛ ⋅ 𝐞ʳ / fʳ))
α, Λ = symbols("α Λ", real = true)
eqᵈᵉᵛ = factor.(simplify.(subs.([eqᶿᵈᵉᵛ,eqʳᵈᵉᵛ], uᶿ(r) => r^α, uʳ(r) => Λ*r^α)))
αΛ = solve([eq.doit() for eq ∈ eqᵈᵉᵛ], [α, Λ])
ûᶿᵈᵉᵛ = sum([Sym("C$(i+2)") * r^αΛ[i][1] for i ∈ 1:length(αΛ)])
ûʳᵈᵉᵛ = sum([Sym("C$(i+2)") * αΛ[i][2] * r^αΛ[i][1] for i ∈ 1:length(αΛ)])
T̂ᶿᵈᵉᵛ = factor(simplify(subs(simplify(𝐓ᵈᵉᵛ ⋅ 𝐞ᶿ / fᶿ), uᶿ(r) => ûᶿᵈᵉᵛ, uʳ(r) => ûʳᵈᵉᵛ)))
T̂ʳᵈᵉᵛ = factor(simplify(subs(simplify(𝐓ᵈᵉᵛ ⋅ 𝐞ʳ / fʳ), uᶿ(r) => ûᶿᵈᵉᵛ, uʳ(r) => ûʳᵈᵉᵛ)))
```
