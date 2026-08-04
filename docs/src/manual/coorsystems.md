# [Coordinate systems (symbolic)](@id man-coorsystems)

Differential operators on tensor fields, with exact symbolic derivatives. The
mathematics — natural basis, Lamé coefficients, Christoffel symbols, the index
placement of each operator — is on
[Curvilinear differential calculus](@ref th-curvilinear).

## Predefined systems

| Constructor | Coordinates | Lamé ``\chi_i`` |
| :--- | :--- | :--- |
| [`coorsys_cartesian`](@ref) | ``(x,y)`` or ``(x,y,z)`` | ``(1,1,1)`` |
| [`coorsys_polar`](@ref) | ``(r,\theta)`` | ``(1,r)`` |
| [`coorsys_cylindrical`](@ref) | ``(r,\theta,z)`` | ``(1,r,1)`` |
| [`coorsys_spherical`](@ref) | ``(\theta,\varphi,r)`` | ``(r,\;r\sin\theta,\;1)`` |
| [`coorsys_spheroidal`](@ref) | ``(\varphi,p,q)`` | prolate spheroidal |

!!! note "Spherical coordinates are ordered ``(\theta,\varphi,r)``"
    Not ``(r,\theta,\varphi)``, so that ``\theta=\varphi=0`` reproduces the
    canonical basis in the canonical order.

## The operators

```@example cs
using TensND, SymPy

Polar = coorsys_polar()
r, θ = getcoords(Polar)
𝐞ʳ, 𝐞ᶿ = unitvec(Polar)
@set_coorsys Polar

LAPLACE(SymFunction("f", real = true)(r, θ))
```

[`@set_coorsys`](@ref) installs a default system so the operators take one
argument; without it, pass the system explicitly, `LAPLACE(f, Polar)`.

| Call | Does |
| :--- | :--- |
| [`@set_coorsys`](@ref)` CS` / [`set_coorsys!`](@ref)`(CS)` | make `CS` the default |
| [`default_coorsys`](@ref)`()` | the current default (throws if none) |
| [`unset_coorsys!`](@ref)`()` | forget it |

!!! note "`∂` is deliberately not affected by the default"
    `∂(t, x)` always means the plain derivative of `t` with respect to the
    symbol `x`; for the covariant derivative pass the chart explicitly,
    `∂(t, x, CS)`.

    Earlier versions had `@set_coorsys` **define methods** in the `TensND`
    module, including one for `∂(t, x)`. That method was more specific than the
    plain-derivative fallback, so the same call silently changed meaning — and
    since both `CoorSystemSym` and `SubManifoldSym` differentiate a position
    vector while being built, constructing *any* new chart after a
    `@set_coorsys` failed. The default is now simply stored, and the
    single-argument operator methods are defined once.

| Operator | Order | Meaning |
| :--- | :--- | :--- |
| [`GRAD`](@ref) | ``+1`` | gradient, derivative index appended **on the right** |
| [`SYMGRAD`](@ref) | ``+1`` | symmetrized gradient; for a displacement, the strain |
| [`DIV`](@ref) | ``-1`` | divergence, contracting the **last** index |
| [`LAPLACE`](@ref) | ``0`` | ``\mathrm{DIV}\circ\mathrm{GRAD}`` |
| [`HESS`](@ref) | ``+2`` | ``\mathrm{GRAD}\circ\mathrm{GRAD}`` |

```@example cs
n = symbols("n", integer = true)
simplify(LAPLACE(r^n * cos(n * θ)))
```

## Accessors

| Function | Returns |
| :--- | :--- |
| [`getcoords`](@ref) | the coordinate symbols |
| [`getOM`](@ref) | the position vector |
| [`unitvec`](@ref) | vectors of the normalized basis |
| [`natvec`](@ref) | vectors of the natural basis (`:cov` or `:cont`) |
| [`normalized_basis`](@ref), [`natural_basis`](@ref) | the two bases |
| [`Lame`](@ref) | the ``\chi_i`` |
| [`Christoffel`](@ref) | the array ``\Gamma[i,j,k]=\Gamma^k_{ij}`` |
| [`∂`](@ref) | the covariant derivative along one coordinate |

```@example cs2
using TensND, SymPy
Spherical = coorsys_spherical()
Lame(Spherical)
```

```@example cs2
Christoffel(Spherical)[1, 1, 3]     # Γʳ_θθ = −r
```

## Custom systems

```julia
CoorSystemSym(OM, coords, tmp_coords = (), params = ();
              rules = Dict(), tmp_var = Dict(), to_coords = Dict())
```

The position vector and the coordinate symbols suffice; the natural basis, Lamé
coefficients and Christoffel symbols are derived. A variant takes the normalized
basis and the ``\chi_i`` directly, when they are known in closed form.

```@example cs3
using TensND, SymPy

x, y = symbols("x y", real = true)
a = symbols("a", positive = true)
Parabolic = CoorSystemSym(Tens([(x^2 - y^2) / 2, x * y]), (x, y))
Lame(Parabolic)
```

## Taming symbolic expressions

Non-trivial charts produce metric expressions SymPy will not reduce unaided, and
the Christoffel symbols become unusable nested radicals. Four optional arguments
control this:

| Argument | Role |
| :--- | :--- |
| `tmp_coords` | auxiliary symbols standing for compound expressions |
| `params` | constants appearing in ``\underline{OM}`` |
| `tmp_var` | substitutions replacing expressions by those symbols |
| `to_coords` | how to eliminate them again before differentiating |
| `rules` | rewrite rules applied after each simplification |

The prolate spheroidal system is the standard illustration, and
[`coorsys_spheroidal`](@ref) is built exactly this way:

```@example cs3
Spheroidal = coorsys_spheroidal()
Lame(Spheroidal)
```

```@example cs3
simplify(LAPLACE(getOM(Spheroidal)[1]^2, Spheroidal))
```

`rules` is also what resolves sign ambiguities such as ``|\sin\theta|``, which
SymPy cannot reduce without knowing that ``\theta\in(0,\pi)`` — see
[Surfaces: fundamental forms and curvature](@ref).

## Symbolic or numerical?

For pointwise evaluation, sweeps and gradients, use
[`CoorSystemNum`](@ref) instead — same operator names, derivatives by automatic
differentiation. See
[Numerical coordinate systems](@ref man-coorsystems-num).
