# [Charts, coordinate systems and surfaces](@id man-coordinates-overview)

!!! tip "See also"
    **Theory**: [Curvilinear differential calculus](@ref th-curvilinear) · [Classical coordinate systems](@ref th-predefined-systems) · [Submanifolds](@ref th-submanifolds)

    **Manual**: [Coordinate systems (symbolic)](@ref man-coorsystems) · [Coordinate systems (numerical)](@ref man-coorsystems-num) · [Charts and coordinate conversions](@ref man-charts) · [Submanifolds](@ref man-submanifolds)

    **Tutorials**: [Differential operators, symbolically](@ref tut-symbolic-operators) · [Christoffel symbols](@ref tut-christoffel) · [Differential operators, numerically](@ref tut-numerical-operators) · [Classical identities of the differential operators](@ref tut-operator-identities) · [The operators in the classical systems](@ref tut-operators-catalog) · [Surfaces: fundamental forms and curvature](@ref tut-submanifolds) · [The spheroid in spheroidal coordinates](@ref tut-spheroidal-appendix)

The package has three kinds of objects for curvilinear coordinates. They answer
different questions, and each of the last two comes in an exact and a
pointwise form. This page says which to take; the pages that follow give the
details.

| to … | take | built with | it holds |
| :--- | :--- | :--- | :--- |
| go from coordinates to a point, and back | a **chart** | `SphericalChart()`, `ToroidalChart(a)`, … | the parameters of the system; two maps between numbers |
| differentiate a field in the coordinates of the space | a **coordinate system** | `coorsys_spherical()` (exact), `coorsys_spherical(PointwiseGeometry)` (pointwise) | the bases, the Lamé coefficients, the Christoffel symbols and the operators `grad`, `symgrad`, `divergence`, `laplacian`, `hess` |
| work on a surface embedded in the space | a **submanifold** | `Submanifold(OM, coords)` (exact), `Submanifold(OM_func, dim)` (pointwise) | the tangent frame and the normal, the two fundamental forms, the same operators along the surface |

## A chart: coordinates and points

A chart maps the coordinates of a point to its Cartesian components and back.
It knows nothing of bases or derivatives, works on numbers and duals, and is
what places a point, or locates one in an inclusion
([Charts and coordinate conversions](@ref man-charts)):

```@example overview
using TensND
chart = SphericalChart()                       # (θ, φ, r)
x = coords_to_cartesian(chart, (π / 3, π / 4, 2.0))
cartesian_to_coords(chart, x)
```

## A coordinate system: fields and their derivatives

A coordinate system adds the differential geometry. Its **exact** form works on
expressions of the coordinate symbols and gives closed forms
([Coordinate systems (symbolic)](@ref man-coorsystems)); its **pointwise**
form takes a field as a function of the coordinates and returns the derivative
as a function, evaluated by automatic differentiation
([Coordinate systems (numerical)](@ref man-coorsystems-num)). The operators
have the same names in both:

```@example overview
using SymPyPythonCall
Spherical = coorsys_spherical()
θ, φ, r = get_coords(Spherical)
laplacian(r^3, Spherical)                      # exact: an expression
```

```@example overview
pointwise = coorsys_spherical(PointwiseGeometry)
laplacian(x -> x[3]^3, pointwise)([π / 3, π / 4, 2.0])     # pointwise: a number, 12 r
```

## A submanifold: a surface in the space

A submanifold is parametrized by its own coordinates, one fewer than the space,
and completes its tangent frame with the unit normal. It gives the first and
the second fundamental forms, hence the curvatures, and the operators along the
surface ([Submanifolds](@ref man-submanifolds)). The exact form is built from
the position vector as an expression, the pointwise one from a function:

```@example overview
R = 2.0
sphere = Submanifold(u -> R * [sin(u[1])cos(u[2]), sin(u[1])sin(u[2]), cos(u[1])], 3)
u = [π / 3, π / 4]
laplacian(v -> cos(v[1]), sphere)(u), -2cos(u[1]) / R^2      # Laplace–Beltrami of cos θ
```

## Exact or pointwise

| | exact | pointwise |
| :--- | :--- | :--- |
| fields | expressions of `Sym` or `Num` symbols | Julia functions of the coordinates |
| results | closed forms, valid everywhere | numbers at a point |
| cost | symbolic simplification, once | automatic differentiation, at each point |
| numbers | symbols | `Float64`, `ForwardDiff.Dual` (derivatives with respect to a parameter) |
| use | derive a solution, a Hill tensor, an identity | evaluate a field, a curvature, a residual, in a loop or under differentiation |
