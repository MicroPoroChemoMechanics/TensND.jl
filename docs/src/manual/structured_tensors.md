# [Structured tensors](@id man-structured)

The compact symmetry-class types. Each stores a handful of scalars instead of
``d^4`` components and computes products and inverses in closed form. The
algebra behind each is on [Isotropic tensors](@ref th-isotropic),
[The Walpole basis](@ref th-walpole),
[The extended Walpole algebra](@ref th-walpole-extended),
[Orthotropy](@ref th-orthotropy) and [Cubic symmetry](@ref th-cubic).

## Storage

Every class is parametrized the same way, by its order and the number `N` of
stored coefficients; the invariance group decides which tensors a shape holds.
The last three columns say what a product (`⊡` at order 4, `⋅` at order 2) of
two tensors of the shape on the same reference, an inverse and a major
transpose return; a sum keeps the shape.

| Type | Order | `N` | Stored | Reference | Group | Product | Inverse | Transpose |
| :--- | :-: | :-: | :--- | :--- | :--- | :--- | :--- | :--- |
| `TensISO{2,dim}` | 2 | 1 | ``\lambda`` | none | ``O(3)`` | same | same | itself |
| `TensISO{4,dim}` | 4 | 2 | ``(\alpha,\beta)`` | none | ``O(3)`` | same | same | itself |
| `TensCubic{T,B}` | 4 | 3 | ``(\alpha,\beta,\gamma)`` | cube frame | ``O_h`` | same | same | itself |
| `TensTI{2,T,2}` | 2 | 2 | ``(a,b)`` | axis ``\underline{n}`` | ``D_{\infty h}`` | same | same | itself |
| `TensTI{2,T,3}` | 2 | 3 | ``(a,b,c)`` | axis ``\underline{n}`` | ``SO(2)`` | same | same | ``c\to-c`` |
| `TensTI{4,T,5}` | 4 | 5 | ``(\ell_1,\ell_2,\ell_3,\ell_5,\ell_6)`` | axis ``\underline{n}`` | ``D_{\infty h}`` | `N=6` | same | itself |
| `TensTI{4,T,6}` | 4 | 6 | ``(\ell_1,\ldots,\ell_6)`` | axis ``\underline{n}`` | ``D_{\infty h}`` | same | same | ``\ell_3\leftrightarrow\ell_4`` |
| `TensTI{4,T,8}` | 4 | 8 | ``(\ell_1,\ldots,\ell_8)`` | axis ``\underline{n}`` | ``SO(2)`` | same | same | ``\ell_3\leftrightarrow\ell_4``, ``\ell_{7,8}\to-\ell_{7,8}`` |
| `TensOrtho{2,T,3}` | 2 | 3 | ``(a_1,a_2,a_3)`` | frame | ``D_{2h}`` | same | same | itself |
| `TensOrtho{4,T,9}` | 4 | 9 | ``(C_{11},\ldots,C_{66})`` | frame | ``D_{2h}`` | `N=12` | same | itself |
| `TensOrtho{4,T,12}` | 4 | 12 | ``+\,(C_{21},C_{31},C_{32})`` | frame | ``D_{2h}`` | same | same | normal block transposed |

The major-symmetric shapes are the ones whose product widens (`N=5 → 6`,
`N=9 → 12`): their product is a product of two symmetric blocks.
`Tensors.majorsymmetric` goes back (`N=6, 8 → 5`, `N=12 → 9`). The ``SO(2)``
shapes hold the terms odd in ``\underline{n}`` that the mirrors of
``D_{\infty h}`` forbid; orthotropy has no such shapes, the rotations of
``D_2`` and the full ``D_{2h}`` having the same invariant tensors at even order.

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

```@example struct
cub = tens_cubic(10.0, 4.0, 2.0, ℬ)      # (C₁₁, C₁₂, C₄₄)
```

## Constructors

| Class | Constructors |
| :--- | :--- |
| ISO | [`TensISO`](@ref), [`tens_Id2`](@ref), [`tens_Id4`](@ref), [`tens_J4`](@ref), [`tens_K4`](@ref), [`iso_projectors`](@ref), [`iso_projectors`](@ref) |
| TI order 4 | [`TensTI`](@ref), [`tens_ti`](@ref), [`tens_ti_eng`](@ref), [`tens_ti_hoenig`](@ref), [`iso_to_ti`](@ref), [`tens_W1`](@ref)…[`tens_W8`](@ref), [`walpole_basis`](@ref), [`walpole_basis_sym`](@ref) |
| TI order 2 | [`TensTI`](@ref)`{2}(a, b, n)` and `(a, b, c, n)` |
| CUBIC | [`TensCubic`](@ref), [`tens_cubic`](@ref), [`iso_to_cubic`](@ref) |
| ORTHO order 4 | [`TensOrtho`](@ref)`{4}(C₁₁, …, C₆₆, frame)` (nine constants) and `(…, C₂₁, C₃₁, C₃₂, frame)` (twelve), `TensOrtho(KM, frame[, Val(12)])`, [`iso_to_ortho`](@ref), [`walpole_to_ortho`](@ref), [`cubic_to_ortho`](@ref) |
| ORTHO order 2 | [`TensOrtho`](@ref)`{2}(a₁, a₂, a₃, frame)`, [`iso_to_ortho`](@ref), [`walpole_to_ortho`](@ref) |

## Accessors

| Function | Applies to | Returns |
| :--- | :--- | :--- |
| [`get_data`](@ref) | all | the stored coefficient tuple |
| [`get_ℓ`](@ref) | `TensTI{4}` | the six classical Walpole coefficients |
| [`get_ℓ8`](@ref) | `TensTI{4}` | all eight, zero-padded |
| [`axis`](@ref) | `TensTI` | the symmetry axis |
| [`frame`](@ref) | `TensOrtho`, `TensCubic` | the material or cube frame |
| [`reference`](@ref) | `TensTI`, `TensOrtho`, `TensCubic` | axis or frame, whichever applies |
| [`symmetry`](@ref) | all | the class as a `Symbol` |
| [`mandel_material`](@ref) | `TensOrtho`, `TensCubic` | the sparse matrix in the class's own frame |
| [`arg_cubic`](@ref) | `TensCubic` | the constants ``(C_{11},C_{12},C_{44})`` |
| [`cubic_anisotropy`](@ref) | `TensCubic` | the Zener-type departure from isotropy |

```@example struct
get_ℓ(ti), axis(ti)
```

```@example struct
mandel_material(ort)
```

## Closure: what an operation returns

This is the table to know. `⊡` between structured types stays structured
**only while the result stays in the class**; otherwise it falls back to the
generic route, exactly and without warning.

| `A ⊡ B` | ISO | CUBIC | TI{4,5} | TI{4,6} | TI{4,8} | ORTHO |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **ISO** | ISO | **CUBIC** | TI{4,6} | TI{4,6} | TI{4,8} | ORTHO{12} |
| **CUBIC** | **CUBIC** | **CUBIC** | generic | generic | generic | ORTHO{12} |
| **TI{4,5}** | TI{4,6} | generic | TI{4,6} | TI{4,6} | TI{4,8} | ORTHO{12} |
| **TI{4,6}** | TI{4,6} | generic | TI{4,6} | TI{4,6} | TI{4,8} | ORTHO{12} |
| **TI{4,8}** | TI{4,8} | generic | TI{4,8} | TI{4,8} | TI{4,8} | generic |
| **ORTHO** | ORTHO{12} | ORTHO{12} | ORTHO{12} | ORTHO{12} | generic | ORTHO{12} |

The `CUBIC` row and column need the two frames to describe the same cube, the
`iso_projectors`/`TI`/`ORTHO` ones the same axis or frame, and a `TI ⊡ ORTHO` product a
TI axis along one of the frame axes; otherwise the generic route is taken.
`ORTHO{12}` stands for `TensOrtho{4, T, 12}`, whichever form the operands have.

Two widenings deserve an explanation, and both have the same cause: **the
major-symmetric tensors of a class do not form a subalgebra**.

- `TI{4,5} ⊡ TI{4,5} → TI{4,6}`. The Walpole product is a ``2\times2`` matrix
  product, and ``L\,M\neq{}^{t}(L\,M)`` unless the blocks commute, so the
  product of two major-symmetric TI tensors generally has ``\ell_3\neq\ell_4``.
- `ORTHO ⊡ ORTHO → ORTHO{12}`. In the material frame the Kelvin-Mandel matrix
  is a ``3\times3`` normal block plus three shear moduli, and the product of two
  symmetric blocks is not symmetric unless they commute: twelve independent
  constants where the major-symmetric form stores nine. The product is one
  ``3\times3`` matrix product and three scalar ones, and inversion inverts the
  block. `Tensors.majorsymmetric` returns the nine-constant symmetric part.

**`CUBIC ⊡ CUBIC → CUBIC` is the one entry that does not widen**, and it is
worth knowing why: each irreducible representation of the octahedral group
appears with multiplicity one, so the three projectors are orthogonal, any two
cubic tensors about the same cube commute, and the product is cubic *and* still
major-symmetric. Cubic symmetry is the only class here that is a genuine
commutative subalgebra. See [Cubic symmetry](@ref th-cubic).

`CUBIC ⊡ TI` is generic for the opposite reason: the two classes are
incomparable, their intersection being the isotropic class.

At order 2, and between orders, the same rules hold on a shared reference:

| Operation | Result |
| :--- | :--- |
| `ORTHO{2} ⋅ ORTHO{2}`, with an `ISO{2}` or a coaxial `TI{2,2}` | `ORTHO{2}` |
| `ORTHO{4} ⊡ ORTHO{2}`, `ORTHO{2} ⊡ ORTHO{4}`, the order-4 operand `iso_projectors`, `CUBIC` or a coaxial `TI{4,5}`/`TI{4,6}` | `ORTHO{2}`: the normal block acting on ``(a_1,a_2,a_3)`` |
| `ORTHO{2} ⊗ ORTHO{2}` | `ORTHO{4,12}`; `ORTHO{4,9}` for ``\boldsymbol{a}\otimes\boldsymbol{a}`` |

as `TI{4} ⊡ TI{2,2} → TI{2,2}` and `TI{2,2} ⊗ TI{2,2} → TI{4,6}` (`TI{4,5}` for a
self product).

Inversion, by contrast, **always stays in the class**:

```@example struct
typeof(inv(iso)), typeof(inv(ti)), typeof(inv(ort)), typeof(inv(cub))
```

## Exact promotions

Moving *up* the lattice is exact — an isotropic tensor really is transversely
isotropic about every axis:

```@example struct
typeof(iso_to_ti(iso, [0.0, 0.0, 1.0]))
```

```@example struct
typeof(iso_to_ortho(iso, ℬ)), typeof(walpole_to_ortho(ti, ℬ, 3))
```

```@example struct
typeof(iso_to_cubic(iso, ℬ)), typeof(cubic_to_ortho(cub))
```

Moving *down* is approximation, and is [Projection](@ref man-projection).

## Same axis, same frame

Addition and subtraction require the two operands to share their orientation —
adding two orthotropic tensors expressed in **different** material frames throws:

```julia
julia> ort + TensOrtho(6., 5., 7., 2., 1., 3., 1.2, 2.2, 0.8, Basis(0.3, 0.7, 0.2))
ERROR: ArgumentError: TensOrtho operation requires the same material frame
```

The assertion is deliberate: silently promoting to a generic tensor would hide a
modeling error. Convert explicitly if that is what you mean.

`TensCubic` is the one class where the test is looser than equality, and
legitimately so: the three projectors are invariant under the whole octahedral
group, so two frames related by a **signed permutation** of the axes describe
the same tensor with the same coefficients and are accepted. Permuting the axes
of an orthotropic tensor permutes ``C_{11},C_{22},C_{33}``; permuting those of a
cubic one changes nothing.

## When to use which

| Situation | Type |
| :--- | :--- |
| the physics guarantees isotropy | `TensISO` |
| a stiffness or compliance with an axis of symmetry | `TensTI{4,T,5}` |
| an object with an axis but no major symmetry (a concentration tensor) | `TensTI{4,T,6}` or `{4,T,8}` |
| a cubic crystal, a cubic array, a cube-symmetric shape | `TensCubic` |
| nine constants and a material frame | `TensOrtho{4,T,9}` |
| an orthotropic object without major symmetry (a product, a concentration tensor) | `TensOrtho{4,T,12}` |
| an order-2 tensor diagonal in a material frame | `TensOrtho{2,T,3}` |
| anything else | a plain `Tens` |

The measured cost of each is in [Performance of the structured types](@ref).
