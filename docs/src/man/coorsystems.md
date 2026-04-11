# Coordinate systems and differential operators

## Symbolic coordinate systems (`CoorSystemSym`)

Symbolic coordinate systems support exact derivation of differential operators using SymPy.
A `CoorSystemSym` stores the position vector `OM`, coordinate symbols, the natural basis
``\mathbf{a}_i = \partial_i \mathbf{OM}``, the normalized (unit) basis ``\mathbf{e}_i = \mathbf{a}_i / \|\mathbf{a}_i\|``,
Lame coefficients ``\chi_i = \|\mathbf{a}_i\|``, and Christoffel symbols
``\Gamma_{ij}^k = \partial_i \mathbf{a}_j \cdot \mathbf{a}^k``.

### Predefined symbolic systems

| Constructor | Coordinates | Description |
| ----------- | ----------- | ----------- |
| `coorsys_cartesian()` | ``(x, y)`` or ``(x, y, z)`` | Cartesian |
| `coorsys_polar()` | ``(r, \theta)`` | Polar |
| `coorsys_cylindrical()` | ``(r, \theta, z)`` | Cylindrical |
| `coorsys_spherical()` | ``(\theta, \varphi, r)`` | Spherical |
| `coorsys_spheroidal(c)` | ``(\varphi, p, q)`` | Prolate spheroidal |

### Example: polar Laplacian

```@repl coorsys
using TensND, SymPy
Polar = coorsys_polar() ; r, θ = getcoords(Polar) ; 𝐞ʳ, 𝐞ᶿ = unitvec(Polar) ;
@set_coorsys Polar
LAPLACE(SymFunction("f", real = true)(r, θ))
n = symbols("n", integer = true)
simplify(HESS(r^n))
```

### Example: spherical divergence

```@repl coorsys2
using TensND, SymPy
Spherical = coorsys_spherical() ; θ, ϕ, r = getcoords(Spherical) ; 𝐞ᶿ, 𝐞ᵠ, 𝐞ʳ = unitvec(Spherical) ;
@set_coorsys Spherical
Christoffel(Spherical)
ℬˢ = normalized_basis(Spherical)
σʳʳ = SymFunction("σʳʳ", real = true)(r) ;
σᶿᶿ = SymFunction("σᶿᶿ", real = true)(r) ;
σᵠᵠ = SymFunction("σᵠᵠ", real = true)(r) ;
𝛔 = σʳʳ * 𝐞ʳ ⊗ 𝐞ʳ + σᶿᶿ * 𝐞ᶿ ⊗ 𝐞ᶿ + σᵠᵠ * 𝐞ᵠ ⊗ 𝐞ᵠ
div𝛔 = simplify(DIV(𝛔))
```

### Differential operators

The following operators are available after calling `@set_coorsys CS`:

| Operator | Signature | Description |
| -------- | --------- | ----------- |
| `GRAD(f)` | scalar or tensor field | Gradient |
| `SYMGRAD(v)` | vector field | Symmetric gradient ``(\nabla \mathbf{v} + \nabla \mathbf{v}^T)/2`` |
| `DIV(t)` | tensor field (order >= 1) | Divergence |
| `LAPLACE(f)` | scalar or tensor field | Laplacian ``\nabla^2 f`` |
| `HESS(f)` | scalar field | Hessian ``\nabla \nabla f`` |

### Accessors

- `getcoords(CS)` / `getcoords(CS, i)`: coordinate symbols
- `getOM(CS)`: position vector
- `normalized_basis(CS)`: orthonormal basis ``(\mathbf{e}_i)``
- `natural_basis(CS)`: natural basis (Basis type including ``\mathbf{a}_i`` and ``\mathbf{a}^i``)
- `unitvec(CS)` / `unitvec(CS, i)`: unit basis vectors as tensors
- `natvec(CS, i, :cov)` / `natvec(CS, i, :cont)`: natural basis vectors
- `Lame(CS)` / `Lame(CS, i)`: Lame coefficients ``\chi_i``
- `Christoffel(CS)`: Christoffel symbols as a 3D array

## Numerical coordinate systems (`CoorSystemNum`)

Numerical coordinate systems evaluate differential operators pointwise using automatic
differentiation via `ForwardDiff`. They do not require a symbolic setup and work with
any numeric type (including `ForwardDiff.Dual` for nested differentiation).

### Predefined numerical systems

| Constructor | Coordinates | Description |
| ----------- | ----------- | ----------- |
| `coorsys_cartesian_num(dim)` | ``(x_1, \ldots, x_d)`` | Cartesian |
| `coorsys_polar_num()` | ``(r, \theta)`` | Polar |
| `coorsys_cylindrical_num()` | ``(r, \theta, z)`` | Cylindrical |
| `coorsys_spherical_num()` | ``(\theta, \varphi, r)`` | Spherical |

### Custom numerical systems

A `CoorSystemNum` can be built from any position vector function `OM(x)`:

```julia
CS = CoorSystemNum(x -> [x[1]*cos(x[2]), x[1]*sin(x[2])], 2)  # polar coordinates
```

### Pointwise evaluation

All operators are functions that return **functions** (lazy evaluation). They must be
called with the coordinate point as argument:

```julia
CS = coorsys_polar_num()
f(x) = x[1]^2              # r^2
lap_f = LAPLACE(CS, f)      # returns a function
lap_f([2.0, π/4])           # evaluate at (r=2, theta=pi/4) => 4.0
```

### Operators (same as symbolic, pointwise)

| Operator | Signature | Returns |
| -------- | --------- | ------- |
| `GRAD(CS, f)` | scalar/tensor field function | function `x -> gradient` |
| `SYMGRAD(CS, v)` | vector field function | function `x -> sym. gradient` |
| `DIV(CS, t)` | tensor field function | function `x -> divergence` |
| `LAPLACE(CS, f)` | scalar/tensor field function | function `x -> Laplacian` |
| `HESS(CS, f)` | scalar field function | function `x -> Hessian` |

### Pointwise accessors

- `normalized_basis(CS, x0)`: orthonormal basis at point `x0`
- `unitvec(CS, x0, i)`: unit vector ``\mathbf{e}_i`` at point `x0`
- `natvec(CS, x0, i, :cov)`: natural covariant vector at point `x0`
- `Lame(CS, x0)` / `Lame(CS, x0, i)`: Lame coefficients at point `x0`

For a detailed tutorial with validation examples (polar, spherical, cylindrical, Lame
problem on a hollow sphere), see [Numerical differential operators](@ref).
