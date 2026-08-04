# [Structured tensors](@id man-structured)

The compact symmetry-class types. Each stores a handful of scalars instead of
``d^4`` components and computes products and inverses in closed form. The
algebra behind each is on [Isotropic tensors](@ref th-isotropic),
[The Walpole basis](@ref th-walpole),
[The extended Walpole algebra](@ref th-walpole-extended) and
[Orthotropy](@ref th-orthotropy).

## Storage

| Type | Stored | Orientation | Class |
| :--- | :--- | :--- | :--- |
| `TensISO{2,dim}` | 1 scalar ``\lambda`` | none | ``\lambda\boldsymbol{1}`` |
| `TensISO{4,dim}` | 2 scalars ``(\alpha,\beta)`` | none | ``\alpha\mathbb{J}+\beta\mathbb{K}`` |
| `TensTI{2,T,2}` | ``(a,b)`` | axis ``\underline{n}`` | ``a\boldsymbol{1}_T+b\boldsymbol{1}_n`` |
| `TensTI{2,T,3}` | ``(a,b,c)`` | axis ``\underline{n}`` | ``+\,c\,\boldsymbol{w}``, the antisymmetric part |
| `TensTI{4,T,5}` | ``(\ell_1,\ell_2,\ell_3,\ell_5,\ell_6)`` | axis ``\underline{n}`` | major-symmetric TI |
| `TensTI{4,T,6}` | ``(\ell_1,\ldots,\ell_6)`` | axis ``\underline{n}`` | general TI |
| `TensTI{4,T,8}` | ``(\ell_1,\ldots,\ell_8)`` | axis ``\underline{n}`` | full axially invariant |
| `TensOrtho{T,B}` | 9 constants | frame `B` | orthotropic |

```@example struct
using TensND, LinearAlgebra

iso = TensISO{3}(3 * 20.0, 2 * 8.0)      # 3k 𝕁 + 2μ 𝕂
```

The `show` method prints the algebraic form, not 81 numbers:

```@example struct
n = (0.0, 0.0, 1.0)
ti = TensTI{4}(12.0, 13.0, 3.5, 7.0, 4.0, n)
```

```@example struct
ℬ = CanonicalBasis{3, Float64}()
ort = TensOrtho(10.0, 8.0, 9.0, 3.0, 2.0, 4.0, 2.5, 3.0, 1.5, ℬ)
```

## Constructors

| Class | Constructors |
| :--- | :--- |
| ISO | [`TensISO`](@ref), [`tens_Id2`](@ref), [`tens_Id4`](@ref), [`tens_J4`](@ref), [`tens_K4`](@ref), [`ISO`](@ref), [`iso_projectors`](@ref) |
| TI order 4 | [`TensTI`](@ref), [`tens_TI`](@ref), [`tens_TI_eng`](@ref), [`tens_TI_Hoenig`](@ref), [`fromISO`](@ref), [`tens_W1`](@ref)…[`tens_W8`](@ref), [`walpole_basis`](@ref), [`walpole_basis_sym`](@ref) |
| TI order 2 | [`TensTI`](@ref)`{2}(a, b, n)` and `(a, b, c, n)` |
| ORTHO | [`TensOrtho`](@ref), [`iso_to_ortho`](@ref), [`walpole_to_ortho`](@ref) |

## Accessors

| Function | Applies to | Returns |
| :--- | :--- | :--- |
| [`get_data`](@ref) | all | the stored coefficient tuple |
| [`get_ℓ`](@ref) | `TensTI{4}` | the six classical Walpole coefficients |
| [`get_ℓ8`](@ref) | `TensTI{4}` | all eight, zero-padded |
| [`axis`](@ref) | `TensTI` | the symmetry axis |
| [`frame`](@ref) | `TensOrtho` | the material frame |
| [`reference`](@ref) | `TensTI`, `TensOrtho` | axis or frame, whichever applies |
| [`symmetry`](@ref) | all | the class as a `Symbol` |
| [`KM_material`](@ref) | `TensOrtho` | the block-diagonal matrix in the material frame |

```@example struct
get_ℓ(ti), axis(ti)
```

```@example struct
KM_material(ort)
```

## Closure: what an operation returns

This is the table to know. `⊡` between structured types stays structured
**only while the result stays in the class**; otherwise it falls back to the
generic route, exactly and without warning.

| `A ⊡ B` | ISO | TI{4,5} | TI{4,6} | TI{4,8} | ORTHO |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **ISO** | ISO | TI{4,6} | TI{4,6} | TI{4,8} | generic |
| **TI{4,5}** | TI{4,6} | TI{4,6} | TI{4,6} | TI{4,8} | generic |
| **TI{4,6}** | TI{4,6} | TI{4,6} | TI{4,6} | TI{4,8} | generic |
| **TI{4,8}** | TI{4,8} | TI{4,8} | TI{4,8} | TI{4,8} | generic |
| **ORTHO** | generic | generic | generic | generic | generic |

Two widenings deserve an explanation, and both have the same cause: **the
major-symmetric tensors of a class do not form a subalgebra**.

- `TI{4,5} ⊡ TI{4,5} → TI{4,6}`. The Walpole product is a ``2\times2`` matrix
  product, and ``L\,M\neq{}^{t}(L\,M)`` unless the blocks commute, so the
  product of two major-symmetric TI tensors generally has ``\ell_3\neq\ell_4``.
- `ORTHO ⊡ ORTHO → generic`. The product keeps the block structure but loses
  major symmetry: twelve independent constants where `TensOrtho` stores nine.
  Rather than add a twelve-parameter container, the generic tensor is returned.

Inversion, by contrast, **always stays in the class**:

```@example struct
typeof(inv(iso)), typeof(inv(ti)), typeof(inv(ort))
```

## Exact promotions

Moving *up* the lattice is exact — an isotropic tensor really is transversely
isotropic about every axis:

```@example struct
typeof(fromISO(iso, [0.0, 0.0, 1.0]))
```

```@example struct
typeof(iso_to_ortho(iso, ℬ)), typeof(walpole_to_ortho(ti, ℬ, 3))
```

Moving *down* is approximation, and is [Projection](@ref man-projection).

## Same axis, same frame

Addition and subtraction require the two operands to share their orientation —
adding two orthotropic tensors expressed in **different** material frames throws:

```julia
julia> ort + TensOrtho(6., 5., 7., 2., 1., 3., 1.2, 2.2, 0.8, Basis(0.3, 0.7, 0.2))
ERROR: AssertionError: TensOrtho operation requires the same material frame
```

The assertion is deliberate: silently promoting to a generic tensor would hide a
modeling error. Convert explicitly if that is what you mean.

## When to use which

| Situation | Type |
| :--- | :--- |
| the physics guarantees isotropy | `TensISO` |
| a stiffness or compliance with an axis of symmetry | `TensTI{4,T,5}` |
| an object with an axis but no major symmetry (a concentration tensor) | `TensTI{4,T,6}` or `{4,T,8}` |
| nine constants and a material frame | `TensOrtho` |
| anything else | a plain `Tens` |

The measured cost of each is in [Performance of the structured types](@ref).
