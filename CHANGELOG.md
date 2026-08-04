# Changelog

## v0.3.0 — non-orthogonal charts fixed, documentation overhaul, two renames

Minor bump rather than patch: two public names change (both keep deprecated
aliases) and `@set_coorsys` no longer defines methods, which is observable.

### Fixed

- **Non-orthogonal coordinate systems were silently wrong.**
  `_build_basis_vectors` exchanged the `:cov` and `:cont` variances, so
  `natvec(CS, i, :cov)` returned the *dual* vector `𝐚ⁱ` instead of `∂OM/∂qⁱ`,
  and `unitvec` likewise. On an orthogonal chart the two coincide, which is why
  all five predefined systems and the whole test suite hid it; on `x = u+v²,
  y = v` the Laplacian of a harmonic function came out nonzero. Charts really
  are arbitrary now. Pinned by the `non-orthogonal charts` testset.

- **`@set_coorsys` broke every subsequent coordinate-system construction.** It
  `@eval`ed single-argument methods into the module, including one for
  `∂(t, x)` that was more specific than the plain-derivative fallback. Since
  `CoorSystemSym` and `SubManifoldSym` both differentiate a position vector
  while being built, that method captured the call, failed to find the
  coordinate among *its own*, returned `zero(t)`, and left the frame matrix
  singular — `NonInvertibleMatrixError`. The default chart is now simply stored
  (`set_coorsys!` / `default_coorsys` / `unset_coorsys!`) and the
  single-argument operator methods are defined once. `∂(t, x)` keeps its single
  meaning: the plain derivative.

- **Automatic differentiation through `proj_tens` was impossible.**
  `ForwardDiff.Dual <: Real` but not `<: AbstractFloat`, so the tolerant
  symmetry predicates never applied to it; a few ulp of round-off made a
  `Dual`-valued tensor look non-minor-symmetric, `_KM_of_array` built a `9×9`
  matrix instead of a `6×6` one, and every call died with a
  `DimensionMismatch`. Predicates now dispatch on `ApproxType`, which includes
  `Dual`.

- **`proj_tens(:ORTHO, A, frame)` returned `NaN` derivatives.** The material
  frame was converted to the *field's* element type, so a `Float64` frame became
  a `Dual` one with zero partials and `angles` evaluated its inverse
  trigonometry at the gimbal-lock point of the canonical frame. The frame keeps
  its own type.

- **`CoorSystemNum` could not be differentiated with respect to a field
  parameter.** The internal buffers took their element type from the Lamé
  coefficients alone; they now promote with the field value.

- **`show(::MIME"text/latex", ::AbstractTens{4})` threw for numeric tensors**
  in any front-end requesting `text/latex` (Documenter, Jupyter). Restricted to
  symbolic element types, which is what it was written for.

- **`DIV(𝟏, CS)` was a `MethodError`**: `tens_Id2(Val(3), Val(Sym))` carries the
  abstract `Sym` while a chart carries `Sym{PyCall.PyObject}`, and the operator
  signatures forced the two to be identical. The field's element type and the
  chart's are now independent.

- `SubManifoldSym` carried a verbatim copy of the `CoorSystemSym` docstring, and
  the `Tens` and `coorsys_cartesian` examples referenced undefined names.

### Renamed

- `Riemann(SM)` → **`connection(SM)`**. It returns the connection coefficients
  (Christoffel symbols) of the induced metric — never a Riemann curvature
  tensor.
- `intrinsic(t)` → **`print_tensor(t)`**. What it prints is explicitly *basis
  dependent*, which is the opposite of intrinsic.

Both old names still work and forward, with a deprecation warning.

### Removed

- `arraytype` was exported but never defined anywhere; `using TensND;
  arraytype` has always been an `UndefVarError`, so removing the export cannot
  break working code.

### Changed

- The symbolic `∂` and the five differential operators are defined **once**,
  over `AbstractCoorSystem`, parametrized by the new `nderiv` trait (`dim` for a
  chart, `dim-1` for a submanifold). `submanifold.jl` used to duplicate all six
  definitions verbatim.
- `CoorSystemSym` and `SubManifoldSym` now share one internal `ChartCore` block
  instead of declaring the same fourteen fields twice.
- Surface indices (`1 … dim-1`) are written `α, β, γ` throughout the
  documentation, ambient indices `i, j, k`.

### Added

- `set_coorsys!`, `default_coorsys`, `unset_coorsys!`, `nderiv`.
- Documentation rebuilt on the MeanFieldHom model: 14 Theory pages, 11 Manual
  pages, 16 Literate tutorials (each also a notebook and a standalone script),
  3 Developer pages and 12 curated API pages, with `checkdocs = :exports` as a
  guard. Bibliography expanded from 4 to 14 entries, every one verified against
  Crossref — the previous `hoenig1979` DOI did not resolve.
- Three chapters with no counterpart in the Echoes manual: bases and variance,
  curvilinear differential calculus, submanifolds. A fourth, the
  eight-dimensional axially invariant space behind `TensTI{4,T,8}`, appears in
  no reference we know of.
- `test/test_conventions.jl`: pins every convention the documentation asserts —
  the Walpole identities (the widely repeated `𝕀 = Σᵢ𝕎ᵢ` and `𝕁 = 𝕎₁+𝕎₂` are
  **false**), the Kelvin–Mandel congruence on a *rotated* frame, the operator
  index placement, `GRAD(𝐧) = −𝐛`, and the classical operator identities on four
  charts including a deliberately skew one.
- `test/test_special_tens.jl` filled out: Levi-Civita, the `𝐞ᵖ/𝐞ᶜ/𝐞ˢ` frames,
  every `init_*`, `rot2`/`rot3`/`rot6`.

## v0.2.6 — TensOrtho closed forms, basis-comparison fix, generic-conversion fix

### Added

- **Closed-form `dcontract(::TensOrtho, ::TensOrtho)`** (and
  `TensISO{4,3} ⊡ TensOrtho`, promoted onto the material frame). Double
  contraction *is* the matrix product in Kelvin-Mandel, and in a shared material
  frame both operands are block-diagonal — `[3×3 sym] ⊕ diag(2C₄₄,2C₅₅,2C₆₆)` —
  so the product is one 3×3 product plus three scalar products instead of two
  dense 81-component expansions. Different material frames still fall back to
  the generic route, bit for bit, since the product is then generally fully
  anisotropic.

  **The result is deliberately not a `TensOrtho`.** The product of two symmetric
  3×3 Kelvin-Mandel blocks is not symmetric unless they commute, so `A ⊡ B` is
  orthotropic *without major symmetry*: 12 independent constants where
  `TensOrtho` stores 9. This is the same widening that already makes
  `dcontract(::TensTI{4}, ::TensTI{4})` return N=6 rather than N=5. Rather than
  introduce a 12-parameter container, the method returns the very
  `TensCanonical` the generic route produced, agreeing with it to 6.2e-16.

- **`_km_congruence`**, the 6×6 Kelvin-Mandel congruence of a material frame,
  built straight from the frame vectors (no angle recovery). It is what relates
  the two Kelvin-Mandel views of a `TensOrtho`:

      KM(t) == Q * KM_material(t) * transpose(Q)

  `Q` is the Kelvin-Mandel representation of `R ⊠ˢ R`, i.e. `KM(rot6(θ,ϕ,ψ))`,
  and it is **orthogonal** — that is precisely what Kelvin-Mandel buys over
  Voigt, where the analogous matrix is not. The identity is pinned in
  `test/test_tens_walpole.jl` on rotated frames, because in the canonical frame
  `Q == I` and every convention error hides.

  `KM_material` now returns an `SMatrix`.

### Performance

- **`getindex(::TensOrtho, i,j,k,l)` was `get_array(t)[i,j,k,l]`** — an
  `Array{T,4}` allocation plus an 81-iteration loop (~40 multiplications each)
  **per scalar access**. Since `TensOrtho <: AbstractArray`, every generic
  algorithm (`Array(t)`, `norm`, broadcasting, `show`, `tomandel`, …) walks the
  tensor element by element, so an O(81) traversal cost O(81²) and 81
  allocations. `TensOrtho` was the last structured type still on the dense
  path; `TensISO` and `TensTI` already had closed forms.

  The per-entry expression is now factored into `_ortho_entry`, shared by
  `get_array` (loop) and `getindex` (single entry), so there is one source of
  truth for the formula. Measured, values bitwise identical:

  | | before | after |
  |---|---|---|
  | `getindex` | 2792 ns / 2304 B | **5 ns / 96 B** |
  | `Array(t)` (81 accesses) | 48 810 ns / 60 448 B | **2240 ns / 832 B** |
  | `get_array` | 575 ns | **311 ns** |

- **`tensor_or_array` built a non-concrete type.** `dim = size(tab, 1)` is a
  *runtime* value, so `Tensor{order, dim}(tab)` was constructed dynamically:
  3147 ns and 3120 bytes for an 81-element array, against 311 ns to produce
  that array in the first place. Every structured tensor reaches the generic
  route through this function (`change_tens` → `same_basis` → every binary
  operation), so the cost was paid twice per `⊡` between structured operands.
  The runtime `dim` is now funnelled through `Val` so each branch constructs a
  concrete type (2 and 3 specialized; anything else keeps the old path).

  `tensor_or_array` 3147 → **393 ns**; `TensOrtho ⊡ TensOrtho`
  10 240 → **4541 ns** (−56 %), 8976 → 4816 B.

- **Comparing two bases cost more than the tensor algebra it guarded.**
  `AbstractBasis <: AbstractMatrix`, and its `getindex` reached the stored
  matrix through `vecbasis(ℬ, :cov)` — the `Symbol` overload, which builds
  `Val(var)` from a *runtime* value, so reading one entry went through a
  dynamic dispatch (67 ns). The `AbstractArray` fallback for `==` then did
  `2·dim²` of them: **1117 ns to compare two 3×3 bases**, which dominated every
  `_check_same_reference` and therefore every binary operation on `TensOrtho`.

  `getindex` now names `Val(:cov)` literally, and `==` compares the primal
  bases in one shot. Semantics are unchanged — still an elementwise comparison,
  so a `RotatedBasis` holding the identity still equals a `CanonicalBasis`.

  Basis `==` 1117 → **162 ns**. Combined with the closed form above,
  `TensOrtho ⊡ TensOrtho` 4541 → **164 ns / 304 B**, i.e. **−98.8 %** against
  the 13 860 ns / 9136 B this contraction cost before v0.2.6.

- **`inv_KM` inferred as `Any`.** The shape-to-type mapping went through
  `select_type_KM[size(v)]`; indexing a `Dict{…, UnionAll}` yields a type as a
  runtime *value*, so `frommandel` was reached by dynamic dispatch. Each shape
  now names its tensor type literally at its own call site, so every branch is
  concrete. (Returning the type from a helper would not work — the union would
  then be over `UnionAll` values, which carries no type information.)
  `select_type_KM` is kept as the documented table.

  `inv_KM` on a 6×6: 228 → **37.5 ns**; inferred return type `Any` →
  `TensCanonical{4,dim,Float64}`. Round-trips for all eight shapes (orders 2
  and 4, symmetric and not, dims 2 and 3) are covered by tests.

### Fixed

- **`proj_tens(:TI, A)` and `proj_tens(:ORTHO, A)` were nondeterministic, and
  intermittently wrong.** The NLopt extension opened with `GD_MLSL`, a
  *stochastic* global optimizer whose generator NLopt seeds from the clock. Two
  calls on the same array returned different angles, and roughly 1 call in 200
  returned a spurious local minimum: on a tensor that is exactly orthotropic
  about the frame `Basis(0.3, 0.7, 0.2)`, `drel ≈ 0.4389` instead of `≈ 1e-13`.
  This surfaced as a CI failure on a single matrix entry, but had nothing to do
  with the platform or the Julia version — it would eventually hit any of them.

  `GD_MLSL` is replaced by a deterministic multi-start. The starting points are
  the eigenstructure candidate (`_candidate_TI_axis` / `_candidate_ORTHO_frame`,
  which are *exact* whenever the tensor genuinely has the symmetry being
  projected onto — the failing case above is one, and they return it to 1.5e-15
  with no optimization at all) followed by a fixed angular grid containing the
  canonical axes and the canonical frame. Each start is refined with
  `LD_TNEWTON`, and the objective is evaluated at every start *and* every refined
  start, so the answer is never worse than the best start.

  A challenger has to beat the incumbent by more than the objective's own
  evaluation noise (`1e-14`) to displace it, and the exact candidate is scanned
  first. Without that floor the fix would have been half a fix: `j = 1 − ‖B‖²/‖C‖²`
  is a difference of O(1) quantities, so near an exact symmetry it lands anywhere
  in ±1e-16 and a refined point could displace the exact candidate by "improving"
  purely in rounding. Visible on a symmetric 3×3, whose eigenframe *is* its
  orthotropic projection: the candidate gives `drel = 4e-16`, the noise-level
  winner `drel ≈ 1e-8` — the angular resolution of the refinement.

  Two properties now hold that did not before, and both are pinned in
  `test/test_nlopt_ext.jl`:

  - repeated calls on the same array return bit-identical results;
  - the optimized projection is never worse than the fixed-axis / fixed-frame
    projection along any grid point — a guarantee by construction, not a
    property of the optimizer's luck.

  Measured against the old strategy over 150 random anisotropic tensors, the
  multi-start never returned a worse objective (0/150) and ran 3-8× faster
  (ORTHO 6.5 → 2.3 ms, TI 8.3 → 1.0 ms), the global pass having been the bulk of
  the cost.

- **Corrupted `rot6` docstring** — stray text (`cde Liv Lehn ϕ`) had replaced
  `cϕ` in one entry of the displayed matrix.

- **`TensTI` could not be built with `Dual` coefficients and a real axis.** The
  inner constructor `TensTI{order,T,N}(::NTuple{N,T}, ::NTuple{3,T})` demanded a
  single element type for both the Walpole coefficients and the symmetry axis.
  But the axis is *geometry* — routinely a literal `(0.0, 0.0, 1.0)` — while the
  coefficients get promoted, most importantly to `ForwardDiff.Dual` whenever a
  homogenization scheme is differentiated with respect to a modulus. Every such
  call failed with

  ```
  MethodError: no method matching TensTI{4, Dual{…}, 5}(::NTuple{5, Dual{…}}, ::Tuple{Float64, Float64, Float64})
  ```

  which made `ForwardDiff` unusable through any scheme carrying a `TensTI` phase
  property. Fixed by an outer constructor that converts both tuples to `T`.


## v0.2.4 — TensOrtho: concrete frame field, closed-form inverse

### Fixed

- **Regression from v0.2.3**: decoupling the material frame's element type from
  the data element type was done by erasing the field to the abstract
  `frame::OrthonormalBasis{3}` (a `Union{CanonicalBasis{3,T},RotatedBasis{3,T}}
  where T`), which made the field non-concrete. Every access boxed and
  dispatched dynamically: `get_array(::TensOrtho)` and a single
  `getindex(::TensOrtho, i,j,k,l)` each allocated **283 264 bytes** (≈440× the
  81-`Float64` array actually needed), even behind a function barrier. Fixed by
  making the frame a genuine type parameter — `TensOrtho{T, B<:OrthonormalBasis{3}}`
  — which keeps `B` independent of `T` (preserving the v0.2.3 ForwardDiff fix)
  while restoring concreteness: `get_array`/`getindex` now allocate ~736 bytes
  (~385× reduction), matching the necessary 81-`Float64` array.

### Changed

- `inv(::TensOrtho)` now uses a closed-form inverse of the documented block
  structure — the upper 3×3 symmetric block via scalar adjugate/determinant,
  and each shear term as `Cₘₘ' = 1/(4Cₘₘ)` — instead of a dense 6×6
  `Matrix`-allocating LU factorization. Verified against the previous dense
  inverse to machine precision (`< 1e-10` on the Kelvin-Mandel matrix), and
  still `ForwardDiff`-compatible (verified against finite differences).

## v0.2.3 — TensOrtho ForwardDiff compatibility

### Fixed

- `TensOrtho` is now `ForwardDiff`-compatible: the material frame's element
  type is decoupled from the data element type (`frame::OrthonormalBasis{3}`
  instead of `{3, T}`), so differentiating w.r.t. the nine elastic constants
  (data `T = ForwardDiff.Dual`) no longer requires — and no longer fails to
  build — a Dual-typed geometric frame.

## v0.2.2 — Full axially-invariant TI algebra (additive)

### Added

- `TensTI{4, T, 8}` — the FULL 8-dimensional space of minor-symmetric
  4th-order tensors invariant under rotations about an axis (the commutant of
  the SO(2) action on Kelvin-Mandel space): the six Walpole coefficients plus
  two antisymmetric azimuthal generators `W₇` (m=1) and `W₈` (m=2). Closed
  under double contraction and inversion via a 2×2 block product and two
  complex products; `get_ℓ8`, `tens_W7`, `tens_W8` accessors; lifts from
  `N=5`/`N=6`. This is what an EXACT azimuthal average of a (generally
  non-major-symmetric) concentration tensor lives in — `ℓ₃ ≠ ℓ₄` and the
  antisymmetric couplings are no longer forced to zero.
- `TensTI{2, T, 3}` — 2nd-order axially-invariant tensor `a·nT + b·nₙ + c·w`
  (`w` the in-plane rotation generator), preserving the antisymmetric in-plane
  part; closed `dot`/`inv` (complex-number algebra in the plane ⊕ scalar on
  the axis).

### Changed

- Binary `±` and `dcontract`/`dot` between two structured TI tensors with
  DIFFERENT axes now fall back to a generic `Tens` result instead of throwing
  an axis-mismatch assertion. This enables accumulation of differently-axed TI
  contributions (e.g. multi-orientation self-consistent estimates). Same-axis
  behaviour is unchanged.

## v0.2.1 — Maintenance

- `[compat]` upper bound for `TimerOutputs` raised to `"0.5, 1"`.
- CI badge restored; Runic badge; DOI badge switched to shields.io with the
  concept DOI (was pointing to a stale per-version DOI).
- Installation instructions updated for registration in Julia's General
  registry (no registry to add beforehand).
- Confirmed each GitHub Release keeps archiving automatically to Zenodo's
  existing concept DOI `10.5281/zenodo.17985768` via the native
  GitHub↔Zenodo integration (no workflow or token needed).
- Retired the Codeberg return path: removed `.forgejo/` workflows and
  `docs/deploy_docs.jl`; GitHub is now the sole home.

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
