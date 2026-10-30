# [Symbolic and numeric](@id man-symbolic-numeric)

<!-- crosslinks: written by .github/scripts/crosslinks.py from docs/crosslinks.toml; do not edit -->
!!! tip "See also"
    **Manual**: [Elastic moduli](@ref man-elastic-moduli)

    **Tutorials**: [Differentiation and interoperability](@ref tut-ad-interop) · [Performance of the structured types](@ref tut-performance) · [Differential operators, symbolically](@ref tut-symbolic-operators) · [The elastic sphere, solved symbolically](@ref tut-sphere-problems)
<!-- /crosslinks -->

`TensND` is generic in its element type. The same code paths run on `Float64`,
on `ForwardDiff.Dual`, on SymPy's `Sym` and on Symbolics' `Num`, which is what
lets a derivation be done symbolically and the result evaluated — or
differentiated — numerically without rewriting anything.

## The four scalar worlds

| Element type | Load | Use for |
| :--- | :--- | :--- |
| `Float64` (and any `Real`) | `using TensND` | numerical work |
| `ForwardDiff.Dual` | `using TensND, ForwardDiff` | derivatives and sensitivities |
| `Sym` | `using TensND, SymPyPythonCall` | symbolic derivation, `dsolve`, `solve` |
| `Num` | `using TensND, Symbolics` | native Julia CAS, code generation |

`SymPyPythonCall` must be loaded explicitly: it is an optional dependency of
TensND, whose SymPy-specific methods (the symbolic simplifications of a `Sym`,
the default coordinates of the coordinate systems) it activates.

### Integer and rational coefficients

A structured tensor stores its coefficients in the type its components need.
The isotropic and cubic projectors and the Walpole tensors ``\mathbb{W}_2`` and
``\mathbb{W}_5`` carry ``1/2`` or ``1/3`` in their components, and ``\mathbb{W}_3``,
``\mathbb{W}_4`` carry ``1/\sqrt{2}``. A fourth-order `TensTI` built from
integers or rationals is therefore stored in floats, and a fourth-order
`TensISO` or a `TensCubic` built from integers too, while rationals stay exact
there; `BigInt` becomes `BigFloat`, and every other type is kept as it is:

```@example sn
using TensND
(eltype(TensISO{3}(3, 2)), eltype(TensISO{3}(3 // 1, 2 // 1)),
    eltype(TensTI{4}(1 // 1, 2 // 1, 1 // 1, 1 // 1, 1 // 1, (0, 0, 1))),
    eltype(walpole_basis((0, 0, 1))[3]))
```

A tensor whose components are integers when its coefficients are, a `TensTI`
or a `TensISO` of order 2 or a `TensOrtho`, keeps an integer type, and its
inverse is in the type of the reciprocal, as `inv` of an integer matrix is.

```@example sn
using TensND, SymPyPythonCall

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
f(α) = get_data(project(get_array(ti_from_constants(10.0, 3.0, 2.5, 12.0, 2.0, [sin(α), 0.0, cos(α)])), :TI, n)[1])[1]

ad = ForwardDiff.derivative(f, 0.3)
fd = (f(0.3 + 1.0e-6) - f(0.3 - 1.0e-6)) / 2.0e-6
(ad, fd)
```

and with the numerical differential operators, in both directions — with respect
to the evaluation point and with respect to a parameter carried by the field.
See [Numerical coordinate systems](@ref man-coorsystems-num).

## Interoperability with `Tensors.jl`

`TensND` stores its data in `Tensors.jl` containers whenever the shape allows
(order 1, 2 or 4, dimension 1 to 3; any other shape keeps a plain array), so
the two compose directly and index symmetries are preserved:

```@example sn4
using TensND, Tensors, LinearAlgebra

st = SymmetricTensor{2, 3}((i, j) -> Float64(i + j))
t = Tens(st)
typeof(get_array(t))
```

```@example sn4
C4 = SymmetricTensor{4, 3}((i, j, k, l) -> Float64(i + j + k + l))
norm(get_array(from_mandel(mandel(Tens(C4)))) - get_array(Tens(C4)))
```

## Performance notes

Symbolic element types are orders of magnitude slower than numeric ones, and
expression size — not operation count — dominates. Three practical consequences:

- **simplify early and often.** An unsimplified intermediate propagates into
  every later expression. This is what the `rules` / `tmp_var` / `to_coords`
  machinery of [`CoordinateSystem`](@ref) exists for.
- **derive symbolically, evaluate numerically.** Obtain the closed form once,
  then substitute with [`tsubs`](@ref) or generate a numeric function from it.
- **prefer a structured type over a dense array** whenever the physics allows —
  the measured gains are one to three orders of magnitude, see
  [Performance of the structured types](@ref tut-performance).

For pointwise numerical work, the pointwise [`CoordinateSystem`](@ref) avoids the symbolic layer
entirely and costs a constant amount per point.
