# [Charts and coordinate conversions](@id man-charts)

!!! tip "See also"
    **Theory**: [Curvilinear differential calculus](@ref th-curvilinear) · [Classical coordinate systems](@ref th-predefined-systems) · [Spheroidal coordinates](@ref th-spheroidal) · [Ellipsoidal coordinates](@ref th-ellipsoidal) · [Toroidal coordinates](@ref th-toroidal) · [Toroidal functions](@ref th-toroidal-harmonics) · [Changes of coordinates](@ref th-coordinate-changes) · [Spheroidal harmonics and their re-expansion](@ref th-spheroidal-harmonics)

    **Manual**: [Charts, coordinate systems and surfaces](@ref man-coordinates-overview) · [Coordinate systems (symbolic)](@ref man-coorsystems) · [Coordinate systems (numerical)](@ref man-coorsystems-num) · [Toroidal functions](@ref man-toroidal-harmonics) · [Spheroidal harmonics](@ref man-spheroidal-harmonics)

    **Tutorials**: [Differential operators, symbolically](@ref tut-symbolic-operators) · [Christoffel symbols](@ref tut-christoffel) · [Differential operators, numerically](@ref tut-numerical-operators) · [Classical identities of the differential operators](@ref tut-operator-identities) · [The operators in the classical systems](@ref tut-operators-catalog) · [The spheroid in spheroidal coordinates](@ref tut-spheroidal-appendix) · [Ellipsoidal coordinates](@ref tut-ellipsoidal) · [An insulating torus in a uniform flux](@ref tut-insulating-torus) · [Changing the focal distance of spheroidal harmonics](@ref tut-spheroidal-reexpansion)

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
| [`ToroidalChart`](@ref) | ``(\varphi,\alpha,\beta)`` | `ToroidalChart(a)` or `ToroidalChart(R, r)` |

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

## Unit vectors and Lamé coefficients

At a point of given coordinates, [`unit_vectors`](@ref) returns the unit
vectors ``\underline{e}_i = (\partial\underline{x}/\partial\xi_i)/\chi_i`` as the
columns of a matrix of Cartesian components, in the order of the coordinates,
and [`lame_coefficients`](@ref lame_coefficients(::AbstractChart, ::Any)) the
coefficients ``\chi_i = \lVert\partial\underline{x}/\partial\xi_i\rVert``. Both
are closed forms in static arrays, which allocate nothing; on the sphere of
radius 2, at the colatitude ``\pi/3`` and the azimuth ``\pi/4``, the columns
are ``(\underline{e}_\theta, \underline{e}_\varphi, \underline{e}_r)``:

```@example charts
ξ = (π / 3, π / 4, 2.0)
unit_vectors(SphericalChart(), ξ)
```

```@example charts
lame_coefficients(SphericalChart(), ξ)                    # (r, r sin θ, 1)
```

## [A tensor that turns with the point](@id man-chart-fields)

A material deposited layer by layer on a sphere, or wound around an axis, has
the same tensor of properties at every point in the frame of the unit vectors
above, which turns with the point. [`ChartTens`](@ref) holds such a field: a
*pattern*, whose components in the canonical basis are those of the field in
the local frame, and the chart. With
``\boldsymbol{Q}(\underline{x})`` the matrix of the unit vectors at the point
``\underline{x}``, the field of a pattern ``\boldsymbol{A}`` of order 2 is
``\boldsymbol{Q}\cdot\boldsymbol{A}\cdot\boldsymbol{Q}^{\!T}``, and likewise
at any order.

A conductivity of value 5 along the radius of a sphere and 1 across it has a
transversely isotropic pattern, whose axis is the radial vector of the frame
(the third one for the spherical chart, [`radial_axis`](@ref)):

```@example charts
F = ChartTens(TensTI{2}(1.0, 5.0, radial_axis(SphericalChart())), SphericalChart())
x = [0.3, -0.5, 0.8]
K = tensor_at(F, x)                       # the tensor at the point x, of axis x/‖x‖
n = x / norm(x)
canonical_components(K) ≈ 5 * n * n' + (I - n * n')
```

[`radial_ti`](@ref) writes the fields that are transversely isotropic about the
radial direction, of a sphere or of a cylinder, from the two coefficients of an
order-2 tensor or from the five moduli of a stiffness: with ``r`` the radial
direction and ``t``, ``t'`` two orthogonal directions across it, Hill's moduli
are ``n = C_{rrrr}``, ``l = C_{rrtt}``,
``k = \tfrac12(C_{tttt} + C_{ttt't'})``,
``m = \tfrac12(C_{tttt} - C_{ttt't'})`` and ``p = C_{rtrt}``.

```@example charts
C = radial_ti(CylindricalChart(); n = 12.0, l = 5.0, k = 19.0, m = 11.0, p = 4.0)
eᵣ = Tens([x[1], x[2], 0.0] / hypot(x[1], x[2]))
eᵣ ⋅ (eᵣ ⋅ tensor_at(C, x) ⋅ eᵣ) ⋅ eᵣ      # n, the stiffness along the radius of the cylinder
```

Fields of the same chart add, and a field is multiplied by a number and
inverted, as their patterns are. In an exact coordinate system the field is a
tensor of constant components in the normalized basis, on which the
differential operators act
([Coordinate systems](@ref man-coorsystems)):

```@example charts
Sph = coorsys_spherical()
θ, ϕ, r = get_coords(Sph)
@syms kₜ::positive kᵣ::positive s::real
𝐊 = tensor_at(radial_ti(SphericalChart(), kₜ, kᵣ), Sph)
T = r^s * cos(θ)
factor(simplify(divergence(𝐊 ⋅ grad(T, Sph), Sph) / (r^(s - 2) * cos(θ))))
```

The temperature ``r^s\cos\theta`` is harmonic for that conductivity when
``k_r\,s(s + 1) = 2k_t``: the exponent is no longer an integer.

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
