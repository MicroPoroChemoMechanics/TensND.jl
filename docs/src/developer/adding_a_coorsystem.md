# [Adding a coordinate system](@id dev-adding-coorsystem)

The most common extension. There are three routes, in increasing order of
effort and of control.

## 1. Numerical, from a position map

The cheapest. Supply a differentiable `OM(x)` and everything else is derived by
automatic differentiation:

```julia
CS = CoorSystemNum(x -> [a * cosh(x[1]) * cos(x[2]), a * sinh(x[1]) * sin(x[2])], 2)
```

Nothing else is required — no symbolic setup, no simplification rules. Use this
whenever pointwise values are enough. See
[Numerical coordinate systems](@ref man-coorsystems-num).

## 2. Symbolic, from a position vector

```julia
CoorSystemSym(OM, coords)
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

Spheroidal = CoorSystemSym(
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
CoorSystemSym(OM, coords, bnorm, χᵢ)
```

This is what the predefined systems do, and it is much faster than
differentiating the position vector.

## Checklist for a new system

1. Does [`Lame`](@ref) return something readable? If not, go to route 3.
2. Do the [`Christoffel`](@ref) symbols match
   ``\partial_i\underline{a}_j\cdot\underline{a}^k``? The check is one
   comprehension — see [Christoffel symbols](@ref tut-christoffel).
3. Is a known harmonic function annihilated by [`LAPLACE`](@ref)? This is the
   strongest single test: it fails unless every connection term is right.
4. If a numerical counterpart exists, do the two agree at a point?

Steps 2–4 are exactly what `test/test_coorsystems.jl` and
`test/test_coorsystems_num.jl` do for the predefined systems.
