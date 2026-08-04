# [Symbolic and numeric](@id man-symbolic-numeric)

`TensND` is generic in its element type. The same code paths run on `Float64`,
on `ForwardDiff.Dual`, on SymPy's `Sym` and on Symbolics' `Num`, which is what
lets a derivation be done symbolically and the result evaluated — or
differentiated — numerically without rewriting anything.

## The four scalar worlds

| Element type | Comes from | Use for |
| :--- | :--- | :--- |
| `Float64` (and any `Real`) | Base | numerical work |
| `ForwardDiff.Dual` | `ForwardDiff` | derivatives and sensitivities |
| `Sym` | `SymPy` | symbolic derivation, `dsolve`, `solve` |
| `Num` | `Symbolics` | native Julia CAS, code generation |

```@example sn
using TensND, SymPy

k, μ = symbols("k μ", positive = true)
ℂ = TensISO{3}(3k, 2μ)
get_data(inv(ℂ))
```

The same with Symbolics:

```@example sn2
using TensND, Symbolics

Symbolics.@variables kn μn
get_data(inv(TensISO{3}(3kn, 2μn)))
```

## What is generic, and what is not

| | tensor algebra | structured types | projection, fixed axis | projection, optimized | symbolic operators | numerical operators |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| `Float64` | ✓ | ✓ | ✓ | ✓ | — | ✓ |
| `ForwardDiff.Dual` | ✓ | ✓ | ✓ | — | — | ✓ |
| `Sym` | ✓ | ✓ | ✓ | — | ✓ | — |
| `Num` | ✓ | ✓ | ✓ | partial | ✓ | — |

The one genuine restriction is the **orientation search**: it is a numerical
optimization, so it needs a numeric element type. Projection onto a class with a
*given* axis or frame is closed-form and works symbolically.

## Symbolic helpers

Six functions apply a symbolic operation elementwise to a tensor, an array or a
scalar, and are **no-ops on numeric types** — so generic code can call them
unconditionally:

| Function | Does | SymPy | Symbolics |
| :--- | :--- | :--- | :--- |
| [`tsimplify`](@ref) | simplify | ✓ | ✓ |
| [`tfactor`](@ref) | factorize | ✓ | — |
| [`tsubs`](@ref) | substitute | ✓ | ✓ |
| [`tdiff`](@ref) | differentiate | ✓ | ✓ |
| [`ttrigsimp`](@ref) | trigonometric simplification | ✓ | — |
| [`texpand_trig`](@ref) | trigonometric expansion | ✓ | — |

```@example sn
tsimplify(3.0), tsimplify(k + k)
```

Being no-ops on `Float64` is the point: a function written with `tsimplify`
sprinkled through it runs unchanged, and at full speed, on numeric input.

## Differentiating

`ForwardDiff` composes with the tensor algebra, the structured types and the
projections:

```@example sn3
using TensND, ForwardDiff, LinearAlgebra

n = [0.0, 0.0, 1.0]
f(α) = get_data(proj_tens(:TI, get_array(tens_TI(10.0, 3.0, 2.5, 12.0, 2.0, [sin(α), 0.0, cos(α)])), n)[1])[1]

ad = ForwardDiff.derivative(f, 0.3)
fd = (f(0.3 + 1.0e-6) - f(0.3 - 1.0e-6)) / 2.0e-6
(ad, fd)
```

and with the numerical differential operators, in both directions — with respect
to the evaluation point and with respect to a parameter carried by the field.
See [Numerical coordinate systems](@ref man-coorsystems-num).

## Interoperability with `Tensors.jl`

`TensND` stores its data in `Tensors.jl` containers whenever the shape allows,
so the two compose directly and index symmetries are preserved:

```@example sn4
using TensND, Tensors, LinearAlgebra

st = SymmetricTensor{2, 3}((i, j) -> Float64(i + j))
t = Tens(st)
typeof(get_array(t))
```

```@example sn4
C4 = SymmetricTensor{4, 3}((i, j, k, l) -> Float64(i + j + k + l))
norm(get_array(inv_KM(KM(Tens(C4)))) - get_array(Tens(C4)))
```

## Performance notes

Symbolic element types are orders of magnitude slower than numeric ones, and
expression size — not operation count — dominates. Three practical consequences:

- **simplify early and often.** An unsimplified intermediate propagates into
  every later expression. This is what the `rules` / `tmp_var` / `to_coords`
  machinery of [`CoorSystemSym`](@ref) exists for.
- **derive symbolically, evaluate numerically.** Obtain the closed form once,
  then substitute with [`tsubs`](@ref) or generate a numeric function from it.
- **prefer a structured type over a dense array** whenever the physics allows —
  the measured gains are one to three orders of magnitude, see
  [Performance of the structured types](@ref).

For pointwise numerical work, [`CoorSystemNum`](@ref) avoids the symbolic layer
entirely and costs a constant amount per point.
