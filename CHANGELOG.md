# Changelog

## v0.2.0 — API unification & TI type fusion (breaking)

### Breaking changes

#### Type fusion: `TensWalpole` removed, merged into `TensTI{4, T, N}`

The historical `TensWalpole{T, N}` has been absorbed into the parametric
`TensTI{order, T, N}` family.  One single type now covers all TI tensors:

| Was                           | Now                          |
|-------------------------------|------------------------------|
| `TensTI{2, T, 2}`             | unchanged                    |
| `TensWalpole{T, 5}`           | `TensTI{4, T, 5}`            |
| `TensWalpole{T, 6}`           | `TensTI{4, T, 6}`            |
| `TensWalpole(ℓ₁,…,ℓ₅,  n)`    | `tens_TI(…)` or `TensTI{4}(ℓ₁,…,ℓ₅, n)` |
| `TensWalpole(ℓ₁,…,ℓ₆, n)`     | `TensTI{4}(ℓ₁,…,ℓ₆, n)`      |

All methods previously dispatched on `::TensWalpole` now dispatch on
`::TensTI{4}` (or `::TensTI{4, <:Any, N}` when the N matters).

#### Naming policy: hybrid `snake_case + UPPERCASE` acronyms

| Was            | Now             |
|----------------|------------------|
| `isISO`        | `is_ISO`        |
| `isTI`         | `is_TI`         |
| `isOrtho`      | `is_ORTHO`      |
| `getaxis`      | `axis`          |
| `getframe`     | `frame`         |
| `getdata`      | `get_data`      |
| `getarray`     | `get_array`     |
| `getbasis`     | `get_basis`     |
| `getvar`       | `get_var`       |
| `getdim`       | `get_dim`       |
| `getorder`     | `get_order`     |
| `tensId2`      | `tens_Id2`      |
| `tensId4`      | `tens_Id4`      |
| `tensJ4`       | `tens_J4`       |
| `tensK4`       | `tens_K4`       |
| `tensTI`       | `tens_TI`       |
| `argTI`        | `arg_TI`        |
| `tensTI_eng`   | `tens_TI_eng`   |
| `argTI_eng`    | `arg_TI_eng`    |
| `tensTI_Hoenig`| `tens_TI_Hoenig`|
| `argTI_Hoenig` | `arg_TI_Hoenig` |
| `tensW1`…`tensW6` | `tens_W1`…`tens_W6` |
| `tensbasis`    | `tens_basis`    |
| `invKM`        | `inv_KM`        |

Type names (`TensISO`, `TensTI`, `TensOrtho`, `Tens`, etc.) follow Julia's
standard `PascalCase` convention and are unchanged.

#### `Walpole(n; sym = true)` split into two functions

To make the return arity predictable from the name, `Walpole(n)` has been
split into:

- `walpole_basis(n)` → `(W₁, W₂, W₃, W₄, W₅, W₆)` (6-tuple, general)
- `walpole_basis_sym(n)` → `(W₁ˢ, W₂ˢ, W₃ˢ, W₄ˢ, W₅ˢ)` (5-tuple, major-sym)

The old `Walpole(n; sym::Bool = false)` is kept as a dispatching alias for
backward compatibility.

#### `ISO()` → `iso_projectors()`

`iso_projectors(Val(dim), Val(T))` returns the `(𝕀, 𝕁, 𝕂)` triple.
`ISO(args...)` is kept as a legacy alias.

### Additions

- **`symmetry(t) :: Symbol`** — single-call query of the material symmetry
  class imposed by the container type (`:ISO`, `:TI`, `:ORTHO`, `:ANISO`).
- **`reference(t)`** — unified accessor that returns `axis(t)` for TI-family,
  `frame(t)` for Ortho, or `nothing` for ISO / unstructured tensors.
- Value-level predicates `is_ISO(A::AbstractArray; ε)`, `is_TI(A, n; ε)`,
  `is_ORTHO(A, frame; ε)` with optional `optimize_angles` kwarg for
  `is_TI` / `is_ORTHO` without reference argument.
- `best_sym_tens(t; …, optimize_angles = false)` — no longer requires NLopt
  by default; the cheap path derives axis/frame candidates from the
  Kelvin-Mandel eigenstructure of the trace tensor.
- Cross-type dispatch extensions (`TensISO ⊡ TensTI{2}`, `TensTI{4} ⊡ TensTI{2}`,
  `dot(TensTI{2}, TensTI{2})`, `TensISO + TensTI{2}`, `TensWalpole{N=5} ±
  TensOrtho` with aligned axis, etc.).

### Migration guide

1. Run `sed -i -E 's/\b(isISO|getaxis|getframe|tensId2|tensId4|tensJ4|tensK4|tensTI|argTI|tensbasis|invKM)\b/.../'`
   on user code with the table above.
2. Replace any explicit `TensWalpole{T, N}` type annotations with
   `TensTI{4, T, N}`; bare `TensWalpole` usages with `TensTI{4}`.
3. `Walpole(n)` → `walpole_basis(n)`; `Walpole(n; sym = true)` →
   `walpole_basis_sym(n)` (or keep the alias — the old signature still works).
4. `ISO(...)` still works; prefer `iso_projectors(...)` in new code.
