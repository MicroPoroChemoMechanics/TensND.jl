# [Coordinate systems (symbolic)](@id man-coorsystems)

!!! tip "See also"
    **Theory**: [Curvilinear differential calculus](@ref th-curvilinear) · [Classical coordinate systems](@ref th-predefined-systems) · [Spheroidal coordinates](@ref th-spheroidal) · [Ellipsoidal coordinates](@ref th-ellipsoidal)

    **Manual**: [Coordinate systems (numerical)](@ref man-coorsystems-num) · [Charts and coordinate conversions](@ref man-charts)

    **Tutorials**: [Differential operators, symbolically](@ref tut-symbolic-operators) · [Christoffel symbols](@ref tut-christoffel) · [Differential operators, numerically](@ref tut-numerical-operators) · [Classical identities of the differential operators](@ref tut-operator-identities) · [The operators in the classical systems](@ref tut-operators-catalog) · [The spheroid in spheroidal coordinates](@ref tut-spheroidal-appendix) · [Ellipsoidal coordinates](@ref tut-ellipsoidal)

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
| [`coorsys_spheroidal`](@ref) | ``(\varphi,p,q)`` | prolate spheroidal, [Spheroidal coordinates](@ref th-spheroidal) |
| [`coorsys_spheroidal_oblate`](@ref) | ``(\varphi,p,\tau)`` | oblate spheroidal |
| [`coorsys_ellipsoidal`](@ref) | ``(\rho,\mu,\nu)`` | ellipsoidal, [Ellipsoidal coordinates](@ref th-ellipsoidal) |

Each curvilinear system is drawn, with its normalized basis and its surfaces of
constant coordinate, in an interactive figure of its Theory page:
[Classical coordinate systems](@ref th-predefined-systems),
[Spheroidal coordinates](@ref th-spheroidal),
[Ellipsoidal coordinates](@ref th-ellipsoidal). A system that is not in the
table is built from its position vector, see [Custom systems](@ref) below.

!!! note "Load SymPyPythonCall for the default coordinates"
    Called without coordinates (`coorsys_spherical()`), a constructor names its
    coordinates with SymPy symbols, which needs `using TensND, SymPyPythonCall`.
    With Symbolics, pass the variables, `coorsys_spherical(Tuple(@variables θ ϕ r))`,
    or the type, `coorsys_spherical(Val(Num))`.

!!! note "Ordering of the spherical coordinates"
    The spherical coordinates are ordered ``(\theta,\varphi,r)``, not
    ``(r,\theta,\varphi)``, so that ``\theta=\varphi=0`` reproduces the
    canonical basis in the canonical order.

## The operators

```@example cs
using TensND, SymPyPythonCall

Polar = coorsys_polar()
r, θ = get_coords(Polar)
𝐞ʳ, 𝐞ᶿ = unitvec(Polar)
@set_coorsys Polar

laplacian(SymFunction("f", real = true)(r, θ))
```

[`@set_coorsys`](@ref) installs a default system so the operators take one
argument; without it, pass the system explicitly, `laplacian(f, Polar)`.

| Call | Does |
| :--- | :--- |
| [`@set_coorsys`](@ref)` CS` / [`set_coorsys!`](@ref)`(CS)` | make `CS` the default |
| [`default_coorsys`](@ref)`()` | the current default (throws if none) |
| [`unset_coorsys!`](@ref)`()` | forget it |

!!! note "`∂` is deliberately not affected by the default"
    `∂(t, x)` always means the plain derivative of `t` with respect to the
    symbol `x`; for the covariant derivative pass the chart explicitly,
    `∂(t, x, CS)`.

    Defining a method for `∂(t, x)` in the `TensND` module would be more specific
    than the plain-derivative fallback, so that the same call would silently
    change meaning, and since both `CoordinateSystem` and `Submanifold`
    differentiate a position vector while being built, constructing any new
    chart after a `@set_coorsys` would fail. The default is therefore simply
    stored, and the single-argument operator methods are defined once.

| Operator | Order | Meaning |
| :--- | :--- | :--- |
| [`grad`](@ref) | ``+1`` | gradient, derivative index appended **on the right** |
| [`symgrad`](@ref) | ``+1`` | symmetrized gradient; for a displacement, the strain |
| [`divergence`](@ref) | ``-1`` | divergence, contracting the **last** index |
| [`laplacian`](@ref) | ``0`` | ``\mathrm{div}\circ\mathrm{grad}`` |
| [`hess`](@ref) | ``+2`` | ``\mathrm{grad}\circ\mathrm{grad}`` |

```@example cs
n = symbols("n", integer = true)
simplify(laplacian(r^n * cos(n * θ)))
```

## Accessors

| Function | Returns |
| :--- | :--- |
| [`get_coords`](@ref) | the coordinate symbols |
| [`position_vector`](@ref) | the position vector |
| [`unitvec`](@ref) | vectors of the normalized basis |
| [`natvec`](@ref) | vectors of the natural basis (`:cov` or `:cont`) |
| [`normalized_basis`](@ref), [`natural_basis`](@ref) | the two bases |
| [`lame_coefficients`](@ref) | the ``\chi_i`` |
| [`christoffel`](@ref) | the array ``\Gamma[i,j,k]=\Gamma^k_{ij}`` |
| [`∂`](@ref) | the covariant derivative along one coordinate |

```@example cs2
using TensND, SymPyPythonCall
Spherical = coorsys_spherical()
lame_coefficients(Spherical)
```

```@example cs2
christoffel(Spherical)[1, 1, 3]     # Γʳ_θθ = −r
```

## Custom systems

A system that is not predefined is built from its position vector, by one of
the routes of [Curvilinear differential calculus](@ref th-curvilinear).

### From the position vector alone

```julia
CoordinateSystem(OM, coords, tmp_coords = (), params = ();
              rules = Dict(), tmp_var = Dict(), to_coords = Dict())
```

The position vector and the coordinate symbols suffice: the natural basis, the
Lamé coefficients, the normalized basis and the Christoffel symbols are
derived. The parabolic coordinates of the plane,
``\underline{OM}=\bigl((x^2-y^2)/2,\;x\,y\bigr)``:

```@example cs3
using TensND, SymPyPythonCall

x, y = symbols("x y", real = true)
OM = Tens([(x^2 - y^2) / 2, x * y])
Parabolic = CoordinateSystem(OM, (x, y))
lame_coefficients(Parabolic)
```

### From the position vector, the basis and the Lamé coefficients

```julia
CoordinateSystem(OM, coords, ℬ, χᵢ, tmp_coords = (), params = (); kwargs...)
```

When the unit vectors and the Lamé coefficients are known in closed form,
passing them spares the computer algebra system the work of finding them, and
keeps them in the form they were written in. Only the Christoffel symbols are
computed. Here ``\underline{a}_x=(x,y)`` and ``\underline{a}_y=(-y,x)``, both
of norm ``\chi=\sqrt{x^2+y^2}``:

```@example cs3
χ = sqrt(x^2 + y^2)
ℬ = RotatedBasis([x/χ -y/χ; y/χ x/χ])
Parabolic = CoordinateSystem(OM, (x, y), ℬ, (χ, χ))
f = SymFunction("f", real = true)(x, y)
tsimplify(laplacian(f, Parabolic))
```

Every predefined system is built this way.

### From a function

When no closed form is wanted, the position vector is given as a function of
the coordinates and the geometry is evaluated at each point by automatic
differentiation: see [Numerical coordinate systems](@ref man-coorsystems-num).

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

The prolate spheroidal system is the standard illustration:
[`coorsys_spheroidal`](@ref) names ``\bar{p}=\sqrt{1-p^2}`` and
``\bar{q}=\sqrt{q^2-1}``, and writes the position vector with them.

```@example cs3
Spheroidal = coorsys_spheroidal()
position_vector(Spheroidal)
```

```@example cs3
lame_coefficients(Spheroidal)
```

```@example cs3
simplify(laplacian(position_vector(Spheroidal)[1]^2, Spheroidal))
```

The ellipsoidal system names its six square roots in the same way:

```@example cs3
Ellipsoidal = coorsys_ellipsoidal()
position_vector(Ellipsoidal)
```

`rules` is also what resolves sign ambiguities such as ``|\sin\theta|``, which
SymPy cannot reduce without knowing that ``\theta\in(0,\pi)`` — see
[Surfaces: fundamental forms and curvature](@ref tut-submanifolds).

## Symbolic or numerical?

For pointwise evaluation, sweeps and gradients, use
the pointwise [`CoordinateSystem`](@ref) instead — same operator names, derivatives by automatic
differentiation. See
[Numerical coordinate systems](@ref man-coorsystems-num).
