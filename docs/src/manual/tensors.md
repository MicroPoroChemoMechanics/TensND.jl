# [Tensors](@id man-tensors)

Building tensors, moving their components between bases and variances, and
reading them back. The underlying notions are on
[Bases and variance](@ref th-bases-variance) and
[Tensor algebra](@ref th-tensor-algebra).

## Construction

A tensor is an array, a basis and a variance tuple. Supplying fewer arguments
fills in the defaults — the canonical basis and fully contravariant variance.

```@example tens
using TensND, SymPy, Tensors

V = Tens(Tensor{1, 3}(i -> symbols("v$i", real = true)))
T = Tens(Tensor{2, 3}((i, j) -> symbols("t$i$j", real = true)))
typeof(V), typeof(T)
```

The type follows the **basis**: a tensor on a canonical basis is a
`TensCanonical`, on a rotated one a `TensRotated`, on an orthogonal one a
`TensOrthogonal`, and only the general case is a `Tens`. Only the last two carry
a variance tuple, because on an orthonormal basis the distinction collapses.

```@example tens
ℬ = Basis(Sym[1 0 0; 0 1 0; 0 1 1])
typeof(Tens(get_array(T), ℬ, (:cov, :cov)))
```

Data may come from a plain `Array`, or from a `Tensors.jl` container
(`Tensor`, `SymmetricTensor`, `Vec`), in which case the symmetry is preserved.

## Components

| Call | Returns |
| :--- | :--- |
| [`components`](@ref)`(t, ℬ, var)` | components on `ℬ` with variance `var` |
| [`components`](@ref)`(t, var)` | components on the tensor's own basis |
| [`components`](@ref)`(t, ℬ)` | keeps the tensor's variance |
| [`components_canon`](@ref)`(t)` | components in the canonical basis |
| [`get_array`](@ref)`(t)` | the **stored** array, as is |

```@example tens
components(V, ℬ, (:cont,))
```

```@example tens
components(V, ℬ, (:cov,))
```

!!! warning "`get_array` is not `components`"
    [`get_array`](@ref) returns the array **as stored**, in the tensor's own
    basis and variance. Two tensors that are mathematically equal can have
    completely different `get_array`s. Comparing raw arrays is therefore
    meaningless unless both tensors are known to share a basis *and* a variance —
    bring them to a common one with [`components`](@ref) or
    [`change_tens`](@ref) first. This is a real source of bugs; a worked
    example is on [Submanifolds](@ref th-submanifolds).

## Changing basis

[`change_tens`](@ref) returns a *tensor* re-expressed on another basis, where
[`components`](@ref) returns a bare array:

```@example tens
𝐞ᶿ, 𝐞ᵠ, 𝐞ʳ = init_spherical()[2]
ℬˢ = init_spherical()[3]
a = Tens(Vec{3}((i,) -> symbols("a$i", real = true)))
A = Tens(rot3(symbols("θ", real = true), symbols("ϕ", real = true)) * get_array(a))
simplify(change_tens(A, ℬˢ))
```

[`change_tens_canon`](@ref) is the shorthand for the canonical basis.

## Accessors

| Function | Returns |
| :--- | :--- |
| [`get_order`](@ref) | the order |
| [`get_dim`](@ref) | the dimension |
| [`get_basis`](@ref) | the basis |
| [`get_var`](@ref) | the variance tuple |
| [`get_data`](@ref) | the stored data (coefficients, for structured types) |
| [`get_array`](@ref) | the full component array |

## Products

The full table, with index formulas, is on
[Tensor algebra](@ref th-tensor-algebra). In practice:

```@example tens2
using TensND, SymPy, Tensors, LinearAlgebra

𝟏 = tens_Id2(Val(3), Val(Sym))
𝕀, 𝕁, 𝕂 = ISO(Val(3), Val(Sym))
𝕀 == 𝟏 ⊠ˢ 𝟏, 𝕁 == (𝟏 ⊗ 𝟏) / 3
```

```@example tens2
a = Tens(Vec{3}((i,) -> symbols("a$i", real = true)))
b = Tens(Vec{3}((i,) -> symbols("b$i", real = true)))
get_array(a ⊗ˢ b)
```

| Julia | Function | Contracted indices |
| :--- | :--- | :--- |
| `⊗` | `otimes` | none |
| `⊗ˢ` | `sotimes` | none, symmetrized |
| `⊠` | `otimesu` | none |
| `⊠ˢ` | `otimesul` | none, symmetrized |
| `⋅` | `dot` | 1 |
| `⊡` | `dcontract` | 2 |
| `⊙` | `qcontract` | 4 |

## Display

[`print_tensor`](@ref) prints a tensor **expanded on its basis** — as a sum of
basis products rather than as a component array, which is far more readable for
symbolic results:

```@example tens2
Spherical = coorsys_spherical()
𝐞ᶿ, 𝐞ᵠ, 𝐞ʳ = unitvec(Spherical)
print_tensor(𝐞ʳ ⊗ 𝐞ʳ + 2 * 𝐞ᶿ ⊗ 𝐞ᶿ, Spherical)
```

Applied to a scalar it falls back to the rich `text/plain` display, so a
symbolic expression comes out in SymPy's two-dimensional form:

```@example tens2
LAPLACE(1 / getcoords(Spherical)[3], Spherical) |> print_tensor
```

Structured types print in their compact algebraic form instead — see
[Structured tensors](@ref man-structured).

## Kelvin–Mandel

[`KM`](@ref) maps a symmetric order-2 tensor to a 6-vector and a minor-symmetric
order-4 tensor to a ``6\times6`` matrix; [`inv_KM`](@ref) inverts it. The
convention and what it buys are on
[Kelvin–Mandel representation](@ref th-kelvin-mandel).

```@example tens2
C = tens_TI(10.0, 3.0, 2.5, 12.0, 2.0, [0.0, 0.0, 1.0])
KM(C)
```

## Symmetry predicates

`issymmetric`, `isminorsymmetric` and `ismajorsymmetric` test the index
symmetries; [`is_ISO`](@ref), [`is_TI`](@ref), [`is_ORTHO`](@ref) test the
**material** symmetry class and are documented under
[Projection](@ref man-projection).
