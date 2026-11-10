# [Adding a coordinate system](@id dev-adding-coorsystem)

The most common extension. There are five approaches, in increasing order of
effort and of control.

## 1. Numerical, from a position map

The cheapest. Supply a differentiable `OM(x)` and everything else is derived by
automatic differentiation:

```julia
CS = CoordinateSystem(x -> [a * cosh(x[1]) * cos(x[2]), a * sinh(x[1]) * sin(x[2])], 2)
```

Nothing else is required — no symbolic setup, no simplification rules. Use this
whenever pointwise values are enough. See
[Numerical coordinate systems](@ref man-coorsystems-num).

## 2. Symbolic, from a position vector

```julia
CoordinateSystem(OM, coords)
```

The natural basis ``\underline{a}_i=\partial_i\underline{OM}``, the Lamé
coefficients, the dual basis and the Christoffel symbols are all computed. This
works out of the box for charts whose metric SymPy can simplify — polar,
cylindrical, parabolic.

## 3. Symbolic, with simplification machinery

For anything harder, the intermediate expressions must be tamed or the result is
correct and unusable. Four optional arguments do it:

| Argument | Role |
| :--- | :--- |
| `tmp_coords` | auxiliary symbols standing for compound expressions |
| `params` | constants appearing in ``\underline{OM}`` |
| `tmp_var` | substitutions replacing expressions by those symbols |
| `to_coords` | how to eliminate them again before differentiating |
| `rules` | rewrite rules applied after each simplification |

The prolate spheroidal chart is the reference example — see the source of
[`coorsys_spheroidal`](@ref):

```julia
ϕ, p = symbols("ϕ p", real = true)
p̄, q, q̄, c = symbols("p̄ q q̄ c", positive = true)
OM = Tens(c * [p̄ * q̄ * cos(ϕ), p̄ * q̄ * sin(ϕ), p * q])

Spheroidal = CoordinateSystem(
    OM, (ϕ, p, q), (p̄, q̄), (c,);
    tmp_var  = Dict(1 - p^2 => p̄^2, q^2 - 1 => q̄^2),
    to_coords = Dict(p̄ => √(1 - p^2), q̄ => √(q^2 - 1)),
)
```

The pattern: introduce a symbol for each radical, declare it `positive` so SymPy
does not carry an absolute value, tell `tmp_var` how to recognize it, and
`to_coords` how to undo it before differentiating.

`rules` is the escape hatch for sign ambiguities that assumptions cannot
express, such as ``|\sin\theta|=\sin\theta`` on ``(0,\pi)``.

## 4. Supplying the basis directly

When the normalized basis and the Lamé coefficients are known in closed form,
give them instead of letting them be derived:

```julia
CoordinateSystem(OM, coords, bnorm, χᵢ)
```

This is what the predefined systems do: only the Christoffel symbols remain to
be computed, and the basis keeps the form it was written in.

## 5. A predefined system, exact and pointwise at once

The spheroidal and ellipsoidal systems are written once, as a **chart**
(`src/charts.jl`): a small type that holds the parameters and three functions
of the coordinates, generic in the number type.

| Function | Returns |
| :--- | :--- |
| [`coords_to_cartesian`](@ref)`(chart, ξ)` | the Cartesian components of the point |
| `TensND._lame(chart, ξ)` | the Lamé coefficients, as a tuple |
| `TensND._frame(chart, ξ)` | the unit vectors, as the columns of a static matrix |

Called with symbols, they give the exact system (`coorsys_ellipsoidal` passes
them to the constructor of approach 4, with the position vector written in
auxiliary symbols as in approach 3); called with numbers, the pointwise one
(`TensND._pointwise_system(chart, T)`). A fourth method,
[`cartesian_to_coords`](@ref), inverts the chart. The ellipsoidal system
([Ellipsoidal coordinates](@ref th-ellipsoidal)) is the reference example: its
six square roots are six auxiliary symbols, and its matrix of unit vectors is
written entry by entry, since Symbolics takes an `hcat` of its own numbers as a
symbolic array.

## Checklist for a new system

1. Does [`lame_coefficients`](@ref) return something readable? If not, go to approach 3.
2. Do the [`christoffel`](@ref) symbols match
   ``\partial_i\underline{a}_j\cdot\underline{a}^k``? The check is one
   comprehension — see [Christoffel symbols](@ref tut-christoffel).
3. Is a known harmonic function annihilated by [`laplacian`](@ref)? This is the
   strongest single test: it fails unless every connection term is right.
4. If a numerical counterpart exists, do the two agree at a point?

5. If the system has a chart, do `_lame` and `_frame` match the Jacobian of
   `coords_to_cartesian`, and does `cartesian_to_coords` invert it, degenerate
   points included?

Steps 2–4 are exactly what `test/test_coorsystems.jl`,
`test/test_coorsystems_num.jl` and `test/test_coorsystems_confocal.jl` do for
the predefined systems, and step 5 what `test/test_charts.jl` does.
