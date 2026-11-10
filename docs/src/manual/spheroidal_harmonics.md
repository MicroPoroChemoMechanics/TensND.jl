# [Spheroidal harmonics](@id man-spheroidal-harmonics)

!!! tip "See also"
    **Theory**: [Spheroidal harmonics and their re-expansion](@ref th-spheroidal-harmonics) · [Spheroidal coordinates](@ref th-spheroidal)

    **Manual**: [Charts and coordinate conversions](@ref man-charts)

    **Tutorials**: [Changing the focal distance of spheroidal harmonics](@ref tut-spheroidal-reexpansion) · [The spheroid in spheroidal coordinates](@ref tut-spheroidal-appendix)

How to evaluate the harmonic functions that separate in a spheroidal chart,
and how to express a harmonic of one chart in the harmonics of another chart
of the same axis. The definitions and the proofs are on
[Spheroidal harmonics and their re-expansion](@ref th-spheroidal-harmonics).

## Legendre functions

Every function returns the values for all the degrees up to `N` at once, the
degree `n` at the index `n + 1`; the entries of degree lower than the order
`m` are zero.

```@example sh
using TensND

legendre_P(3, 1, 0.5)          # P₀¹, P₁¹, P₂¹, P₃¹ at x = 0.5
```

The radial factors depend on the chart: ``P_n^m(q)`` and ``Q_n^m(q)`` in a
prolate chart, their real counterparts ``\mathcal{P}_n^m(\tau)`` and
``\mathcal{Q}_n^m(\tau)`` in an oblate chart, ``r^n`` and ``r^{-n-1}`` in the
spherical chart.

```@example sh
chart = ProlateSpheroidalChart(1.0)
radial_P(chart, 3, 1, 2.0)
```

```@example sh
radial_Q(chart, 3, 1, 2.0)
```

[`radial_Q`](@ref) stays accurate at high degree, where the upward recurrence
would lose the function to rounding:

```@example sh
radial_Q(chart, 60, 0, 2.0)[[1, 21, 41, 61]]
```

## Harmonics at a point

[`interior_harmonics`](@ref) and [`exterior_harmonics`](@ref) take the
coordinates of the point in the chart and return the products
``P_n^m(p)\,\mathsf{P}_n^m(\zeta)`` and ``P_n^m(p)\,\mathsf{Q}_n^m(\zeta)``.
Multiplied by ``\cos m\varphi`` or ``\sin m\varphi``, they are harmonic.

```@example sh
x = [0.7, -0.4, 2.1]
ξ = cartesian_to_coords(chart, x)              # (ϕ, p, q)
interior_harmonics(chart, 3, 1, ξ)
```

```@example sh
exterior_harmonics(chart, 3, 1, ξ)
```

## Interior harmonics in another chart

An interior harmonic of degree ``n`` of a chart is a finite sum of the
interior harmonics of degrees ``n, n-2, \ldots`` of any coaxial chart.
[`interior_reexpansion`](@ref) returns the coefficients, the target chart
first:

```@example sh
from, to = ProlateSpheroidalChart(1.0), ProlateSpheroidalChart(1.5)
L = interior_reexpansion(to, from, 4, 0)       # degrees 4, 2, 0
```

```@example sh
ξ₂ = change_coords(to, from, ξ)
interior_harmonics(from, 4, 0, ξ)[5], sum(L .* interior_harmonics(to, 4, 0, ξ₂)[5:-2:1])
```

The two charts may be prolate, oblate or spherical in any combination. In the
spherical chart the harmonics are the solid harmonics
``r^nP_n^m(\cos\theta)``:

```@example sh
interior_reexpansion(SphericalChart(), from, 2, 0)    # in r² P₂(cos θ) and 1
```

## Exterior harmonics in another chart

An exterior harmonic of degree ``n`` is an infinite series of the exterior
harmonics of degrees ``n, n+2, \ldots`` of the other chart.
[`exterior_reexpansion`](@ref) returns its first `imax + 1` coefficients:

```@example sh
K = exterior_reexpansion(to, from, 2, 1, 30)
K[1:4]
```

```@example sh
H = exterior_harmonics(to, 2 + 2 * 30, 1, ξ₂)[3:2:end]     # degrees 2, 4, …, 62
exact = exterior_harmonics(from, 2, 1, ξ)[3]
[abs(sum(K[1:N+1] .* H[1:N+1]) / exact - 1) for N in (2, 5, 10, 20)]
```

The series converges beyond the coordinate surface of the target chart that
passes through the edge of the focal set of the source chart, which
[`reexpansion_bound`](@ref) returns as a value of the radial coordinate:

```@example sh
reexpansion_bound(to, from), ξ₂[3]       # the point is beyond the bound
```

```@example sh
reexpansion_bound(from, to)              # the other way: q₁ must exceed c₂/c₁
```

## Derivatives and symbols

The Legendre functions and the coefficients accept `ForwardDiff` duals. The
functions of the first kind and the coefficients, which involve no iteration,
also accept symbols:

```@example sh
using SymPyPythonCall
c₁, c₂ = symbols("c₁ c₂", positive = true)
interior_reexpansion(ProlateSpheroidalChart(c₂), ProlateSpheroidalChart(c₁), 2, 0)
```

[`radial_Q`](@ref) sizes an iteration on the value of its argument: it takes
numbers.
