# [Charts and coordinate conversions](@id man-charts)

!!! tip "See also"
    **Theory**: [Spheroidal coordinates](@ref th-spheroidal) · [Ellipsoidal coordinates](@ref th-ellipsoidal) · [Curvilinear differential calculus](@ref th-curvilinear) · [Changes of coordinates](@ref th-coordinate-changes) · [Spheroidal harmonics and their re-expansion](@ref th-spheroidal-harmonics)

    **Manual**: [Coordinate systems (symbolic)](@ref man-coorsystems) · [Spheroidal harmonics](@ref man-spheroidal-harmonics)

    **Tutorials**: [The spheroid in spheroidal coordinates](@ref tut-spheroidal-appendix) · [Ellipsoidal coordinates](@ref tut-ellipsoidal) · [Changing the focal distance of spheroidal harmonics](@ref tut-spheroidal-reexpansion)

How to pass from the coordinates of a point to its Cartesian components and
back, and from one coordinate system to another. The formulas are on
[Changes of coordinates](@ref th-coordinate-changes).

## A chart

A chart holds the parameters of a coordinate system and maps numbers to
numbers. The predefined ones:

| Chart | Coordinates | Built with |
| :---- | :---------- | :--------- |
| [`PolarChart`](@ref) | ``(r,\theta)`` | `PolarChart()` |
| [`CylindricalChart`](@ref) | ``(r,\theta,z)`` | `CylindricalChart()` |
| [`SphericalChart`](@ref) | ``(\theta,\varphi,r)`` | `SphericalChart()` |
| [`ProlateSpheroidalChart`](@ref) | ``(\varphi,p,q)`` | `ProlateSpheroidalChart(c)` |
| [`OblateSpheroidalChart`](@ref) | ``(\varphi,p,\tau)`` | `OblateSpheroidalChart(c)` |
| [`EllipsoidalChart`](@ref) | ``(\rho,\mu,\nu)`` | `EllipsoidalChart(h₂, h₃)` or `EllipsoidalChart(a₁, a₂, a₃)` |

The coordinates are ordered as in the system of the same name
([Coordinate systems (symbolic)](@ref man-coorsystems)).

## From coordinates to a point and back

[`cartesian_to_coords`](@ref) returns the coordinates of a point, and
[`coords_to_cartesian`](@ref) the point of given coordinates:

```@example charts
using TensND

chart = ProlateSpheroidalChart(1.5)      # semi-focal distance c = 1.5
x = [0.4, -0.3, 1.2]
ξ = cartesian_to_coords(chart, x)        # (ϕ, p, q)
```

```@example charts
coords_to_cartesian(chart, ξ)
```

The coordinate ``q`` is half the sum of the distances from the point to the
foci ``\pm c\,\underline{e}_3``, divided by ``c``:

```@example charts
using LinearAlgebra
(norm(x - [0, 0, 1.5]) + norm(x + [0, 0, 1.5])) / (2 * 1.5)
```

## From one system to another

[`change_coords`](@ref) composes the two maps. The same point in an oblate
chart of semi-focal distance `0.8`, then in the spherical chart:

```@example charts
change_coords(OblateSpheroidalChart(0.8), chart, ξ)
```

```@example charts
change_coords(SphericalChart(), chart, ξ)       # (θ, ϕ, r)
```

## The ellipsoidal chart

An ellipsoidal chart is defined by its two semi-focal distances ``h_2>h_3``, or
by the semi-axes ``a_1>a_2>a_3`` of one ellipsoid of the family, on which
``\rho=a_1``:

```@example charts
ell = EllipsoidalChart(3.0, 2.0, 1.0)
x = [1.8, 1.2, 0.5]
x[1]^2 / 9 + x[2]^2 / 4 + x[3]^2          # inside the ellipsoid: < 1
```

```@example charts
ρ, μ, ν = cartesian_to_coords(ell, x)
```

The chart covers the quarter space ``x_2\ge0``, ``x_3\ge0``; a point of
another quarter gets the coordinates of its mirror image, ``\nu`` keeping the
sign of ``x_1``:

```@example charts
coords_to_cartesian(ell, cartesian_to_coords(ell, [-1.8, -1.2, 0.5]))
```

## Derivatives

Both maps accept `ForwardDiff` duals, in the point and in the parameters of
the chart. The Jacobian of the inverse map has the gradients of the
coordinates as rows, ``\underline{e}_i/\chi_i``:

```@example charts
using ForwardDiff
J = ForwardDiff.jacobian(x -> cartesian_to_coords(chart, x), [0.4, -0.3, 1.2])
```

and the derivative with respect to the semi-focal distance tells how the
coordinates of a fixed point move with the foci:

```@example charts
ForwardDiff.derivative(c -> cartesian_to_coords(ProlateSpheroidalChart(c), [0.4, -0.3, 1.2]), 1.5)
```

## Symbols

The direct map has no branch and accepts symbols; it returns the position
vector of the system:

```@example charts
using SymPyPythonCall
c, ϕ, p, q = symbols("c ϕ p q", real = true)
coords_to_cartesian(ProlateSpheroidalChart(c), (ϕ, p, q))
```

The inverse maps select a root or a sign from values: they take numbers.

## A pointwise system from a chart

The spheroidal and ellipsoidal systems exist in pointwise form, with the
operators of [Numerical coordinate systems](@ref man-coorsystems-num):

```@example charts
CS = coorsys_ellipsoidal(PointwiseGeometry, 2.0, 1.0)     # h₂ = 2, h₃ = 1
lame_coefficients(CS, [2.5, 1.5, 0.4])
```

The harmonics that separate in the spheroidal charts, and their change from
one chart to another, are on
[Spheroidal harmonics](@ref man-spheroidal-harmonics).
