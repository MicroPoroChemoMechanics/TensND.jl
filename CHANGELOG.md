# Changelog

## v0.5.0 — cubic symmetry as a first-class class

`TensND` could store an isotropic tensor in 2 scalars, a transversely isotropic
one in 5, an orthotropic one in 9 — and had nothing for the class every
cube-symmetric object lands in. A cubic crystal, a cubic array of inclusions, a
supersphere: all have **three** independent constants, and until now went
through the 21-component generic route.

`TensCubic` fills that gap with the same surface as the other three: storage,
constructors, accessors, closed-form algebra, mixed-class promotions,
projection, exact rotation-group average, symbolic and `ForwardDiff` support,
documentation and tests.

### The class is the one that is *closed*

It is worth stating what is different about it, because it is not merely "one
more class".

Under the octahedral group the Kelvin-Mandel space splits as
`A1g + Eg + T2g`, dimensions `1 + 2 + 3`, three inequivalent irreducible
representations each of **multiplicity one**. The commutant is therefore
spanned by three mutually orthogonal projectors, and the algebra is
componentwise:

    C = α J + β E + γ T ,   C⁻¹ = α⁻¹J + β⁻¹E + γ⁻¹T ,
    A : B = α_A α_B J + β_A β_B E + γ_A γ_B T .

Two consequences no other class here enjoys:

- **the product of two cubic tensors is cubic**, where `TensTI{4,T,5} ⊡
  TensTI{4,T,5}` widens to `N=6` and `TensOrtho ⊡ TensOrtho` leaves its class
  entirely. Cubic symmetry is a genuine commutative subalgebra;
- **major symmetry is automatic**, the three projectors being symmetric. A
  strain localization tensor, which in general has none, recovers it as soon as
  the morphology is cube-symmetric.

The geometry, by contrast, is orthotropy's: unlike an isotropic tensor, a cubic
one is only cubic relative to its cube axes, so the type carries a frame. With
one difference that is specific to the class and is implemented rather than
merely noted: the three projectors are invariant under the *whole* octahedral
group, so two frames related by a **signed permutation** of the axes describe
the same tensor with the same coefficients, and binary operations accept them.
Permuting the axes of an orthotropic tensor permutes `C₁₁, C₂₂, C₃₃`; permuting
those of a cubic one changes nothing. The test is "`ᵗF_A F_B` is a signed
permutation matrix", nine comparisons rather than enumerating 24 rotations.

### At order two, cubic *is* isotropic

The octahedral group leaves no second-order tensor invariant but a multiple of
the identity, so there is no `TensCubic` at order 2 and
`proj_tens(Val(:CUBIC), ::AbstractArray{T,2}, frame)` returns the isotropic
projection.

Not a technicality: it is why a cube-symmetric pore has a single **scalar**
resistivity contribution while its compliance contribution needs three
constants — and therefore why a conduction computation on such a morphology
carries no anisotropy signal at all, whatever the shape does in elasticity.

### Added

- `TensCubic`, `tens_cubic(C₁₁, C₁₂, C₄₄, frame)`, `arg_cubic`,
  `cubic_anisotropy` (the Zener-type departure from isotropy), `is_CUBIC`.
- `KM` and `KM_material`, `get_array` and a closed-form `getindex`, `inv`,
  `one`, `zero`, `literal_pow`, `transpose`, the symmetry predicates,
  `symmetry` (`:CUBIC`), `reference`, `show` and `pprint`.
- Exact promotions `iso_to_cubic` and `cubic_to_ortho`, completing the lattice
  `TensISO ⊂ TensCubic ⊂ TensOrtho`, with the mixed arithmetic that keeps
  `TensISO ± TensCubic` and `TensISO ⊡ TensCubic` in the cubic class and
  `TensCubic ± TensOrtho` in the orthotropic one. `TensCubic` and `TensTI` are
  **incomparable** — their intersection is the isotropic class — so their sum
  falls through to the unstructured route on purpose.
- `proj_tens(Val(:CUBIC), A, frame)` at a fixed cube frame and
  `proj_tens(Val(:CUBIC), A)` optimizing the orientation, the second through
  `TensNDNLoptExt` exactly as the TI and orthotropic searches; `best_fit_cubic`;
  the value-level predicates `is_CUBIC(A, frame)` and `is_CUBIC(A)` alongside
  the type-level one; and the Kelvin-Mandel block helpers
  `cubic_params_from_KM` / `KM_from_cubic_params`.
- `best_sym_tens` now considers the class: the default cascade is
  `(:ISO, :CUBIC, :TI, :ORTHO)`, ordered by **number of constants** — 2, 3, 5,
  9 — and not by inclusion, which it cannot be, `:CUBIC` and `:TI` being
  incomparable. A tensor satisfying both is reported cubic, having the fewer
  constants. Its fixed-orientation variant takes one argument for two kinds of
  orientation, so the classes that argument cannot serve are now **dropped**
  rather than reaching `proj_tens` and raising a `MethodError` — which is what
  already happened whenever `:ORTHO` was reached with an axis.
- A closed-form `isotropify(::TensCubic) = α J + ((2β + 3γ)/5) K`, the weights
  being the dimensions of `Eg` and `T2g`. Exact, and it avoids expanding 81
  components — which on a symbolic element type is the difference between an
  answer and an expression swell.

### The one place cubic is *less* well served, and why

The transversely isotropic and orthotropic searches each start from an
eigenstructure candidate, exact whenever the tensor really has the symmetry
sought. Cubic symmetry admits none, and the reason is not numerical:
**every second-order contraction of a cubic tensor is isotropic.** Measured on
`tens_cubic(10, 4, 2, frame)`, in any frame, `C_iikl = 18 δ_kl` and
`C_ikil = 14 δ_kl`, so no eigenframe carries any information about where the
cube points — the cube axes are a genuinely fourth-order feature. Handing the
orthotropic candidate a cubic tensor returns an arbitrary frame, at a relative
residual of 0.099 against 3.5e-16 for the true one.

So the free-orientation search starts from the angular grid alone, which
contains the canonical frame; and a *rotated* cubic tensor handed over as a bare
array is not recognized on `best_sym_tens`'s cheap path — it needs
`optimize_angles = true`, or the frame. Both facts are stated in the theory page
and pinned by tests rather than left to be discovered.

### Documentation

- A theory page, `theory/cubic.md`, at the same depth as `theory/orthotropy.md`:
  where the three constants come from, why the algebra is closed, the frame
  subtlety, order two, the rotational average, and the class lattice.
- An API page, `api/cubic.md`; the projection family is listed with its
  siblings on `api/projection.md`.
- `manual/structured_tensors.md` extended throughout — storage, constructors,
  accessors, the closure table (with the `CUBIC ⊡ CUBIC → CUBIC` entry
  explained), the promotions, the frame rule and the "when to use which" table.
- `theory/orthotropy.md` now says why *it* loses closure where cubic keeps it,
  which is the same fact seen from the other side, and its lattice diagram
  includes the new class.

### Tests

161 new assertions in `test_tens_cubic.jl` and 27 more in
`test_nlopt_ext.jl`, in the shape of the existing ones: construction and traits,
the two Kelvin-Mandel forms and their congruence, automatic major symmetry on
several frames and coefficient sets, the componentwise algebra with `C ⊡ C⁻¹`
against the identity, closure of the product and its agreement with the dense
route, frames describing the same cube accepted and a general rotation refused,
every promotion, an orthogonal and idempotent projection with its residual
orthogonal to the class, the rotational average against the generic route and
its invariance under rotating the cube, `ForwardDiff` through the inverse, the
product, the projection and the average, and a symbolic element type.


## v0.4.1 — Aqua.jl, and the dispatch defects it found

TensND is now checked by [Aqua.jl](https://github.com/JuliaTesting/Aqua.jl) on
every run of the test suite. Aqua audits properties no functional test looks
at: method ambiguities, type parameters bound by nothing, undefined exports,
dead dependencies, missing `[compat]` bounds, and type piracy. The audit found
real defects, and this release fixes them. Nothing in the exported API changes
name or meaning.

### Bug fixes

- **Seventeen method ambiguities.** Each was a pair of signatures where neither
  is more specific than the other, so the call raised an ambiguity
  `MethodError` instead of running. They came from `AbstractTens` being an
  `AbstractMatrix` (`Basis(ℬ, :cov)`, `is_ORTHO(t, frame)`), from
  zero-dimensional arrays reaching the `Tens` constructors, and from the
  isotropic contraction pair. Every one of them now dispatches.

- **Contracting two isotropic tensors of different dimension** produced an
  ambiguity error the caller could do nothing with. It now throws a
  `DimensionMismatch` naming both dimensions.

- **Five signatures had a type parameter bound by nothing** — `NTuple{N, T}`
  also matches the empty tuple, which leaves `T` undetermined, and
  `Union{T, AbstractTens{order, dim}} where {T <: SymType}` leaves `T`
  undetermined for a tensor argument and `order`/`dim` undetermined for a
  scalar one. `GRAD`, `SYMGRAD`, `LAPLACE` and `HESS` are now two methods each
  instead of one union; `CoorSystemSym`, `coorsys_cartesian`, `TensISO` and
  `eltype_of` now require at least one element, which is what they always
  meant. Both argument types every one of these accepted before are still
  accepted.

### Dependencies

- **`TimerOutputs` moved from `[deps]` to `[extras]`.** It was used only by the
  test suite, never by `src/` or `ext/`, so installing TensND was pulling in a
  package it never loaded.

- **`[compat]` bounds added for `LinearAlgebra`, `Random` and `Test`.** Standard
  libraries need a bound like anything else.


## v0.4.0 — the symmetry-average layer moves in, and `tlimit`

Three groups of functions that had grown up inside
`MeanFieldHomogenization.jl` are pure tensor algebra with nothing
homogenization-specific about them, and they belong here: the exact
rotation-group averages, the Kelvin-Mandel helpers they are built on, and the
best-fit projections. Moving them removes an awkward seam — `MeanFieldHomogenization`
had to `import TensND: isotropify` to avoid two rival bindings of the same
name reaching its top level — and puts the *average* and the *projection* side
by side, where the difference between them can be stated once instead of being
re-explained at every call site.

### Breaking changes

Nothing in the public API was removed or changed in meaning. The bump is a
minor one because the exported surface grows substantially, and because under
1.0 Julia's resolver treats a minor bump as breaking: a downstream package
pinned to `TensND = "0.3"` must widen its bound to `"0.4"`.

One behavior *does* change, and it is the bug fix below: code that called
`tsimplify` / `tsubs` / `tdiff` on a `TensISO`, a `TensTI` or a `TensOrtho`
with symbolic entries used to get its input back unchanged and now gets the
operation actually applied.

### Bug fixes

- **`tsimplify`, `tfactor`, `tsubs`, `tdiff`, `ttrigsimp` and `texpand_trig`
  were silent no-ops on every structured tensor type.** `structured_tens_ops.jl`
  implements the whole family as `_rebuild(A, OP(get_data(A)))`, and `get_data`
  returns an `NTuple`. A `Tuple` is not an `AbstractArray`, so the elementwise
  methods never matched and the generic identity fallback `OP(x, args...) = x`
  caught the call: `tsimplify(::TensISO{4,3,Sym})` returned the tensor
  untouched and reported success. Only the generic `Tens`, whose data is an
  `Array`, ever simplified.

  ```julia
  t = TensISO{3}(sin(a)^2 + cos(a)^2, 2b)
  get_data(tsimplify(t))     # before: (sin(a)^2 + cos(a)^2, 2b)   after: (1, 2b)
  ```

  `Tuple` methods are now defined for the whole family.

- **Every symmetry predicate answered `false` on a `Symbolics.Num` tensor**, and
  `ismajorsymmetric` threw outright. `iszero(a - b)` is not an equality test on
  a `Num`: Symbolics folds `x - x` to zero only when `x` is an *atom*, so two
  components built from the same expression — `(3k − 2μ)/3` against
  `(3k − 2μ)/3` — left an unsimplified sum and `iszero` said `false`.
  `ismajorsymmetric` was worse still, testing `a - b == zero(Num)`, which builds
  a symbolic *equation* and made `!(…)` throw `TypeError: non-boolean (Num) used
  in boolean context`.

  The consequence was not local, and it is the exact Symbolics analog of the
  `ForwardDiff.Dual` breakage v0.3.6 cured: `_store_symmetric` kept the full 9×9
  form for a `Num` fourth-order tensor where a `Float64` one gets 6×6, so every
  consumer expecting the Kelvin-Mandel form of a minor-symmetric tensor was
  handed the wrong shape.

  The three predicates now go through `_num_equal`, which tries SymbolicUtils'
  *structural* equality first — exact, cheap, and decisive for the case that
  matters — and falls back to `iszero`. Neither calls `simplify`: this runs on
  every tensor construction.

- **`KM` of an isotropic tensor threw on symbolic entries**, and so did
  displaying one. `KM(::TensISO)` went through
  `tomandel(SymmetricTensor{order,dim}(A))`, and building a `SymmetricTensor`
  from a full component array makes `Tensors` validate the symmetry — a check
  that concluded "not symmetric" for the reason above and raised an
  `InexactError`. `show` goes through `KM`, so a Documenter `@example` block
  that merely *displayed* a symbolic isotropic tensor failed the whole
  documentation build.

  `KM(::TensISO)` is now closed form at both orders — `λ𝟏` is the vector whose
  first `dim` entries are `λ`, and `α𝕁 + β𝕂` is
  `β I_m + ((α−β)/dim) E` with `E` the leading `dim × dim` block of ones — so
  nothing is materialized at order `dim^4` and no symmetry has to be re-proved.
  Verified identical to the previous path on `Float64` at both orders and both
  dimensions, and still `ForwardDiff`-traversable.

- The new test suite's rotation oracle built its rotation matrix as
  `I(3) .+ sin(φ) .* K .+ …`. On Julia 1.14-DEV a broadcast against
  `I(3)::Diagonal{Bool}` **preserves the `Diagonal` structure and silently
  drops every off-diagonal entry**, so the "rotation" came back diagonal
  (`det = 0.596`) and the oracle disagreed with a correct answer by 0.6. The
  library itself was never affected — every intermediate of
  `transverse_isotropify` is bit-identical across 1.12, 1.13 and 1.14-DEV —
  but the trap is worth recording: use `Matrix(1.0I, n, n) + …`, not a
  broadcast against `I(n)`, whenever the result must be dense.

### Exact rotation-group averages (new here, moved from `MeanFieldHomogenization`)

`isotropify(t)` and `transverse_isotropify(t, n)` are the *exact* averages of a
minor-symmetric tensor over SO(3) and over the rotations about `n`. They assume
no major symmetry, which is what makes them the right operators on a
concentration or contribution tensor: the azimuthal average returns the full
eight-coefficient `TensTI{4,T,8}`, preserving `ℓ₃ ≠ ℓ₄` and the antisymmetric
couplings `ℓ₇`, `ℓ₈` that a symmetric TI parametrization drops silently. At
order 2 the antisymmetric in-plane part is likewise kept (`TensTI{2,T,3}`).

`isotropify` already existed here as a Frobenius projection on a raw array; the
new `AbstractTens` method coincides with it exactly on minor-symmetric input.

Two closed-form fast paths come with them, both new:

- `isotropify(::TensTI{4})` reads the two invariants off the Walpole
  coefficients — `T_iijj = ℓ₁ + 2ℓ₂ + √2(ℓ₃+ℓ₄)`, `T_ijij = ℓ₁ + ℓ₂ + 2ℓ₅ + 2ℓ₆`
  — instead of materializing 81 components. On a symbolic element type that is
  the difference between an expression one can read and one that fills a
  screen; on a numeric polycrystal self-consistent iteration it removes an
  array contraction per phase per iteration.
- `isotropify(::TensTI{2})` and `isotropify(::TensISO)`, same idea.

Also moved, and now public rather than internal:

| Name | Role |
| :--- | :--- |
| `mandel66_minor(arr)` / `array_from_mandel66(M)` | 6×6 Kelvin-Mandel ↔ 3×3×3×3, minor-symmetrizing on read, no major symmetry assumed |
| `ti8_params_from_KM(M)` / `KM_from_ti8_params(p)` | the eight coefficients of an axially-invariant tensor about `e₃` — the non-major-symmetric counterpart of `ti_params_from_KM`, which *projects* onto the five-coefficient span and would discard `ℓ₃ ≠ ℓ₄` |
| `ti_average_mandel66(M, n)` / `iso_average_mandel66(M)` | the same averages on a 6×6 block rather than on a tensor |

`ti8_params_from_KM` doubles as an exact **read-off** whenever the matrix is
already axially invariant, which is how a laminate recovers the symmetry class
of its localization tensors.

The in-plane reference of `_axis_frame` is now chosen structurally when the
axis is not comparable, so an azimuthal average about a symbolic axis no longer
runs into `argmin`.

### Best-fit projections (new here, moved from `MeanFieldHomogenization`)

`best_fit_iso(t)`, `best_fit_ti(t, axis)` and `best_fit_ortho(t, frame)` return
the projection alone, where `proj_tens` returns `(projection, distance,
relative distance)`. They sit next to `proj_tens` and the `*_params_from_KM`
conversions, and their docstrings state the one thing not to get wrong: a
projection is a least-squares fit and drops whatever does not fit, an average
is exact and lossless on the invariant subspace. Use the averages inside a
computation and the fits only to report parameters.

### `tlimit` — the limit passage, keeping the symmetry class

```julia
tlimit(x, s, v)          # limit as the symbol `s` tends to `v`
tlimit(x, s, v, "+")     # one-sided
```

The pass that was missing from the `tsimplify` / `tsubs` / `tdiff` family. Like
its siblings it is the identity on numeric types — a number has no free symbol
to send anywhere — and, applied to a structured tensor, it takes the limit of
the few canonical coefficients and rebuilds the same type rather than going
through `dim^order` components.

That is what makes a degenerate limit of a closed form practical: a vanishing
regularization, an incompressible phase (`k => oo`), a flat inclusion. There is
deliberately no `Symbolics.Num` method — Symbolics has no limit, and returning
the input unchanged would be a silent wrong answer, so a `Num` raises instead.

### `is_hard_numeric` (new here, moved from `MeanFieldHomogenization.Elliptic`)

Whether comparisons on a scalar type yield an honest `Bool`, so a value may
drive an `if`, an `&&` or a tolerance test. `T <: Real` is **not** that
predicate: `Symbolics.Num <: Real` yet answers no comparison, and
`ForwardDiff.Dual <: Real` and does. It is the companion of `ApproxType`, which
this package already owned — `ApproxType` says *use a tolerance rather than
exact equality*, `is_hard_numeric` says *a comparison is allowed at all* — so
the two now live in the same file.

## v0.3.6 — the round-off tolerance now covers `ForwardDiff.Dual`

### Bug fixes

- **Automatic differentiation was impossible through any structured tensor
  written about a non-canonical axis.** `KM` of a `Dual`-valued `TensISO` in a
  `RotatedBasis` returned a **9×9** matrix where the very same tensor in
  `Float64` returned a 6×6 one, so every downstream consumer that expects the
  Kelvin-Mandel form of a minor-symmetric tensor failed with a
  `DimensionMismatch` — and `inv` of the full form solves a 9×9 system whose
  three minor-antisymmetric directions are null.

  This is the case v0.3.2 left behind. That release made the storage decision
  relative to the tensor's own scale instead of exact, precisely because writing
  a structured tensor about a non-canonical axis leaves minor-antisymmetric
  residue of order `1e-16` — but it keyed the tolerant method on
  `AbstractFloat`, and its note justified the exact branch by "there is no
  round-off to absorb there". That is true of a symbolic or rational element
  type. It is **false of a `ForwardDiff.Dual`**, which subtypes `Real` without
  subtyping `AbstractFloat` and carries exactly the round-off of the float it is
  built on.

  `_store_symmetric` now keys on `ApproxType` — the union this package already
  defined for this very distinction (`array_utils.jl`), whose own comment
  describes the symptom — so a `Dual` is treated like the float underneath it.
  `Complex{<:AbstractFloat}`, also in that union, joins the tolerant side at the
  same time; its tolerance goes through the new internal `_approx_eps`, since
  `eps` is defined on a `Dual` type but not on a `Complex` one.

  Concretely, this is what makes `ForwardDiff` reach through a rotation:

  ```julia
  ℬʳ = RotatedBasis(0.3, 0.4, 0.1)
  f(x) = Matrix(KM(TensISO{3}(6.0, 2x), ℬʳ))[4, 4]
  ForwardDiff.derivative(f, 0.8)        # used to throw; now returns 2.0
  ```

  A tensor whose minor antisymmetry is genuine — above the threshold — is stored
  as before. Symbolic and rational element types are untouched: they keep the
  exact comparison, where they belong.

## v0.3.3 — the outer products no longer go through a contraction engine

### Performance

- **`otimes`, `otimesu`, `otimesl` and `sotimes` were 20–55× slower than the
  products they compute.** Each built an `EinCode` and called `einsum` — a
  contraction engine — for operations that contract nothing: their output
  indices are the union of the input indices, one multiplication per element
  and no summation anywhere. `otimesu` and `otimesl` merely *interleave* the
  two operands' indices; `sotimes` is the average of `otimes` and `otimesu`.

  The time was machinery — `EinCode` construction, code selection, dispatch —
  not arithmetic.

  All four are now a single broadcast. Both index lists are increasing, so
  neither operand needs its axes permuted: giving each one singleton axes where
  the other's indices sit places every element correctly. `sotimes` fuses its
  two terms into one broadcast, turning three allocations into one.

  | case | before | after |
  | :--- | ---: | ---: |
  | `otimes(3×3, 3×3)` | 3517 ns / 4144 B | **163 ns / 864 B** |
  | `otimes(3, 3)` | 3240 ns / 3200 B | **50.7 ns / 240 B** |
  | `otimesu(3×3, 3×3)` | 8337 ns / 7680 B | **158 ns / 864 B** |
  | `otimesl(3×3, 3×3)` | 8180 ns / 7632 B | **159 ns / 864 B** |
  | `sotimes(3×3, 3×3)` | 8810 ns / 10016 B | **307 ns / 992 B** |

  Verified **bit-for-bit against the einsum implementations they replace**, on
  nine shape combinations per function — including non-square, order-3 and
  mixed-order operands, and a first-order second operand — with `===` element
  equality, the `ForwardDiff.Dual` element type preserved, and the full test
  suite green. No behavior changes.

  Broadcasting also keeps these generic over `Dual` and symbolic element types,
  where BLAS could not have been used. `sotimes` deliberately keeps its
  division by 2 rather than a multiplication by `0.5`: the two are bit-identical
  in binary floating point but render differently for symbolic elements.

  Note for anyone reading an intermediate commit: a first draft of `otimes`
  used `vec(t1) .* transpose(vec(t2))`, which is the same layout and faster on
  plain arrays, but **raises on a first-order `Tensors` operand** — `transpose`
  of a `Tensors.Vec` is deliberately discontinued upstream. Singleton axes
  transpose nothing and have no such restriction; the slightly slower figure
  above is that correction.

- **`best_sym_tens` derived candidates it never read.** The transversely
  isotropic axis and the orthotropic frame were both computed on every call —
  each from its own fresh `Array(get_array(t))`, so the array was materialized
  three times — even when `proj` asked for neither. Both are now derived only
  when `proj` names the symmetry that reads them, and from the array already
  materialized.

  | case | before | after |
  | :--- | ---: | ---: |
  | `best_sym_tens(C; proj = (:ISO,))` | 11.3 µs / 14016 B | **7.1 µs / 7568 B** |
  | `best_sym_tens(C)` (all three) | 17.6 µs / 22640 B | **17.3 µs / 21168 B** |

  The full call gains only the duplicate materializations, since it genuinely
  needs both candidates; a restricted `proj` gains the eigen-decomposition it
  was never going to use.

### Tests

- New `test/test_tensor_products.jl` pins the four outer products against three
  independent references — the definition as explicit loops, the `einsum` form
  they replaced, and the algebraic identities they must satisfy — across
  shapes, element types (`Rational`, `Int`, `Dual`, symbolic) and, in
  particular, the `Tensors` types the callers actually pass.

  Because this is an N-dimensional package, the generic case is covered rather
  than assumed. **There is no order limit in the code** — `ntuple(…, Val(n))`
  and singleton axes are generic, and the products were checked by hand up to
  9 ⊗ 9 (order 18, 262144 elements). What the tests bound is their own run
  time, not the implementation: every order pair up to 5 ⊗ 3 at dim 3, a few
  lopsided splits up to 9 ⊗ 1 at dim 2, and — since it costs nothing and is a
  statement about index arithmetic alone — the property the broadcast rests on
  checked up to order 20 on each side: that both index lists are increasing, so
  neither operand's axes need reordering, and that together they cover the
  output positions exactly once. **Variance** is
  checked on a deliberately non-orthonormal basis, where co- and contravariant
  components differ: the variance tuple must be carried in index order, and
  `⊠` must interleave it exactly as it interleaves the indices. The invariant
  used is basis-independent — the canonical components of a product are the
  product of the canonical components.

  The contraction paths (`contract`, `dcontract`, `dotdot`, `qcontract`) still
  go through `einsum` and are deliberately unchanged; they are exercised at
  arbitrary orders too, so a future attempt to give them the same treatment
  cannot land unnoticed.

  Mutation-checked, three ways: swapping two index lists, reinstating the
  `transpose` draft, and hard-coding the output order to 4 (correct at 2 ⊗ 2,
  wrong above) each fail the suite. Also pins that restricting `proj` in
  `best_sym_tens` does not change the numbers reported for a given symmetry.

## v0.3.2 — minor symmetry is decided by the type, not by round-off

### Bug fixes

- **`inv` could throw `SingularException` on a well-posed problem.** Whether a
  fourth-order tensor was stored as a 36-component `SymmetricTensor` or an
  81-component `Tensor` was decided by `Tensors.issymmetric`, an exact
  comparison. Writing a structured tensor about a non-canonical axis leaves
  minor-antisymmetric residue of order `1e-16`, so a tensor whose components
  are of size 30 landed on the full form where the same tensor written about
  `e₃` landed on the symmetric one — and the full form inverts through a 9×9
  system whose three minor-antisymmetric directions are null.

  Two stiffnesses that were **equal to the last bit** therefore behaved
  differently: `inv(𝕀 + ℙ : δℂ)` returned a localization tensor for one and
  raised for the other. The symmetry test is now relative to the tensor's own
  scale (`16 eps` of its largest component), and the symmetric branch stores
  `Tensors.symmetric(t)` rather than a bare `convert`, so the residue is
  cleaned instead of carried. Non-float element types (symbolic, rational)
  keep the exact comparison: there is no round-off to absorb there.

  A tensor whose minor antisymmetry is genuine — above that threshold — is
  stored as before. Code that dispatched on the *storage* of a tensor
  antisymmetric at round-off level may now see `SymmetricTensor` where it saw
  `Tensor`; both remain `<: AbstractTens`, and the values agree.

### Internal

- `src/tens_walpole.jl` → **`src/tens_anisotropic.jl`** (and its test file with
  it). The name had outgrown its contents: the file holds `TensTI` *and*
  `TensOrtho`, and only the former is expressed in a Walpole basis. It is now
  the counterpart of `tens_isotropic.jl` — the structured anisotropic types.
  No code moved; nothing exported changed.

## v0.3.1 — `print_tensor` renamed to `pprint`

### Deprecations

Nothing breaks on upgrade: `print_tensor` and `intrinsic` are still exported and
still forward, so this is a patch release, not a minor one.

- `print_tensor(t)` → **`pprint(t)`** (v0.3.0's own rename of `intrinsic`, so a
  short-lived name). Two reasons. Its **fallback method accepts anything**, and
  a scalar field is the common case — `LAPLACE(f(r)) |> print_tensor` printed a
  scalar under a name promising a tensor. And `pprint` is the name already used
  across MPCM for a generic pretty-printer: `ChemistryLab.pprint` has methods on
  `Formula`, `Species`, `Reaction` and on a plain `AbstractMatrix`, so it is not
  a domain-specific name.

  Note that `pprint` is therefore exported by two MPCM packages with no
  dependency between them: in a session doing `using MeanFieldHomogenization, ChemistryLab`
  the unqualified `pprint` is ambiguous and must be written `TensND.pprint` or
  `ChemistryLab.pprint`.

### Bug fixes

- The docstring example printed its terms in the wrong order: spherical
  coordinates are ordered `(θ, ϕ, r)`, so `𝐞ʳ⊗𝐞ʳ + 2 𝐞ᶿ⊗𝐞ᶿ` displays as
  `(2)𝐞ᶿ⊗𝐞ᶿ + 𝐞ʳ⊗𝐞ʳ`. The block is ```` ```julia ````, hence not executed by
  Documenter, which is why nothing caught it.

### Additions

- `test/test_tens.jl` gains a `pprint` testset: the expanded layout, the
  omission of zero components, the `0` of a vanishing tensor, the `vec`/`coords`
  keywords, the variance deciding sub- or superscript, and both deprecated
  aliases resolving.

### Internal

- The formatting half of `pprint` is split into `_pprint_string`, so the layout
  is testable against a string instead of requiring `stdout` capture.

## v0.3.0 — non-orthogonal charts fixed, documentation overhaul, two renames

Minor bump rather than patch: two public names change (both keep deprecated
aliases) and `@set_coorsys` no longer defines methods, which is observable.

### Breaking changes

- `Riemann(SM)` → **`connection(SM)`**. It returns the connection coefficients
  (Christoffel symbols) of the induced metric — never a Riemann curvature
  tensor.
- `intrinsic(t)` → **`print_tensor(t)`**. What it prints is explicitly *basis
  dependent*, which is the opposite of intrinsic.

  Both old names still work and forward, with a deprecation warning.

- **`@set_coorsys` no longer `@eval`s methods into the module**; it stores the
  default chart, alongside the new `set_coorsys!` / `default_coorsys` /
  `unset_coorsys!`. The single-argument operator methods are defined once,
  unconditionally. The observable consequence: `∂(t, x)` now always means the
  plain derivative with respect to the symbol `x`, and the covariant derivative
  is spelled `∂(t, x, CS)`.

  This was also a bug fix. The `@eval`ed `∂(t, x)` method was more specific than
  the plain-derivative fallback, and since `CoorSystemSym` and `SubManifoldSym`
  both differentiate a position vector *while being built*, that method captured
  the call, failed to find the coordinate among *its own*, returned `zero(t)`,
  and left the frame matrix singular — `NonInvertibleMatrixError` on every chart
  constructed after a `@set_coorsys`.

- `arraytype` is no longer exported. It was exported but never defined anywhere,
  so `using TensND; arraytype` has always been an `UndefVarError`; removing the
  export cannot break working code.

### Bug fixes

- **Non-orthogonal coordinate systems were silently wrong.**
  `_build_basis_vectors` exchanged the `:cov` and `:cont` variances, so
  `natvec(CS, i, :cov)` returned the *dual* vector `𝐚ⁱ` instead of `∂OM/∂qⁱ`,
  and `unitvec` likewise. On an orthogonal chart the two coincide, which is why
  all five predefined systems and the whole test suite hid it; on `x = u+v²,
  y = v` the Laplacian of a harmonic function came out nonzero. Charts really
  are arbitrary now. Pinned by the `non-orthogonal charts` testset.

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

### Additions

- `set_coorsys!`, `default_coorsys`, `unset_coorsys!`, `nderiv`.
- Documentation rebuilt on the MeanFieldHomogenization model: 14 Theory pages, 11 Manual
  pages, 16 Literate tutorials (each also a notebook and a standalone script),
  3 Developer pages and 12 curated API pages, with `checkdocs = :exports` as a
  guard. Bibliography expanded from 4 to 14 entries, every one verified against
  Crossref — the previous `hoenig1979` DOI did not resolve. Surface indices
  (`1 … dim-1`) are written `α, β, γ` throughout, ambient indices `i, j, k`.
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

### Internal

- The symbolic `∂` and the five differential operators are defined **once**,
  over `AbstractCoorSystem`, parametrized by the new `nderiv` trait (`dim` for a
  chart, `dim-1` for a submanifold). `submanifold.jl` used to duplicate all six
  definitions verbatim.
- `CoorSystemSym` and `SubManifoldSym` now share one internal `ChartCore` block
  instead of declaring the same fourteen fields twice.

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
  behavior is unchanged.

## v0.2.1 — Maintenance

- `[compat]` upper bound for `TimerOutputs` raised to `"0.5, 1"`.
- CI badge restored; Runic badge; DOI badge switched to shields.io with the
  concept DOI (was pointing to a stale per-version DOI).
- Installation instructions updated for registration in Julia's General
  registry (no registry to add beforehand).
- Confirmed each GitHub Release keeps archiving automatically to Zenodo's
  existing concept DOI `10.5281/zenodo.17985768` via the native
  GitHub↔Zenodo integration (no workflow or token needed).

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
