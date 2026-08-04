# [Testing and conventions](@id dev-testing)

The test suite does two jobs: it checks that the code is correct, and it **pins
the conventions the documentation asserts**, so that the two cannot drift apart.

## Running the tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

`test/runtests.jl` seeds the generator (`Random.seed!`) so failures are
reproducible, wraps each file in a `@testsection` with `TimerOutputs`, and
includes the files in a **deliberate order**.

!!! warning "`test_nlopt_ext.jl` must stay last"
    Loading `NLopt` permanently activates `TensNDNLoptExt` for the session, and
    `test/test_tens_projection.jl` asserts the behavior of the *cheap*,
    no-NLopt path. Moving the extension test earlier makes the projection tests
    exercise a different code path than intended.

## Conventions pinned by tests

Each of the following is a statement made somewhere in the documentation that a
test would catch if the code changed under it.

| Convention | Where it is stated | Why it needs a test |
| :--- | :--- | :--- |
| ``\mathbb{I}=\mathbb{W}_1+\mathbb{W}_2+\mathbb{W}_5+\mathbb{W}_6`` | [The Walpole basis](@ref th-walpole) | the widely repeated ``\mathbb{I}=\sum_i\mathbb{W}_i`` is **false** |
| ``\mathbb{J}=\tfrac13(\mathbb{W}_1+2\mathbb{W}_2+\sqrt2\mathbb{W}_3+\sqrt2\mathbb{W}_4)`` | [The Walpole basis](@ref th-walpole) | likewise, ``\mathbb{J}\neq\mathbb{W}_1+\mathbb{W}_2`` |
| the Walpole Gram matrix is ``\mathrm{diag}(1,1,1,1,2,2)`` | [The Walpole basis](@ref th-walpole) | it is what makes every TI projection closed-form |
| ``\mathrm{Mat}(t)=Q\,\mathrm{Mat}_{\text{mat}}(t)\,{}^{t}Q`` | [Kelvin–Mandel](@ref th-kelvin-mandel) | must be checked on a **rotated** frame: in the canonical frame ``Q=\boldsymbol{1}`` and every convention error hides |
| ``\mathbb{W}_7,\mathbb{W}_8`` annihilate axially invariant order-2 tensors — and *only* those | [The extended Walpole algebra](@ref th-walpole-extended) | a source comment once claimed they annihilate every symmetric tensor, which is false |
| the ``N=8`` space is closed under `⊡` and `inv` | [The extended Walpole algebra](@ref th-walpole-extended) | it is the justification for the storage type |
| ``\Gamma^k_{ij}=\partial_i\underline{a}_j\cdot\underline{a}^k`` | [Curvilinear calculus](@ref th-curvilinear) | the storage order ``\Gamma[i,j,k]`` is easy to transpose |
| `GRAD` appends on the right, `DIV` contracts the last index | [Curvilinear calculus](@ref th-curvilinear) | invisible on symmetric fields, wrong by a transpose otherwise |
| ``\mathrm{GRAD}(\underline{n})=-\boldsymbol{b}`` on a sphere | [Submanifolds](@ref th-submanifolds) | ties the normal orientation to the curvature sign |

## What a good test looks like here

Prefer an **identity that fails unless everything is right** over a check
against a transcribed number.

- The Laplacian of a harmonic function vanishes only if every Christoffel term
  is correct. `LAPLACE(rⁿcos nθ) == 0` in polar coordinates, and the spheroidal
  harmonics ``P_n^m(p)P_n^m(q)\cos m\varphi``, are the strongest single checks
  on a coordinate system.
- A closed form matched against the operator route — the Green's function pages
  do this — exercises `HESS`, the symmetrization and the basis handling at once.
- A round trip (`arg_TI(tens_TI(x...)) == x`) catches parametrization errors
  that a one-way check cannot.
- Equilibrium residuals: `DIV(σ) ≈ 0` for an exact solution tests the full
  order-2 divergence including connection terms.

## Non-obvious traps

**Comparing `get_array`s.** Two mathematically equal tensors stored on different
bases have different component arrays. Always bring them to a common basis and
variance with [`components`](@ref) before subtracting — see
[Tensors](@ref man-tensors).

**Element types that are not `AbstractFloat`.** `ForwardDiff.Dual <: Real` but
**not** `<: AbstractFloat`. The tolerant symmetry predicates are therefore
declared on `ApproxType` (`src/array_utils.jl`), which includes `Dual`; before
that union existed, a few ulp of round-off made a `Dual`-valued tensor look
non-minor-symmetric, `_KM_of_array` built a ``9\times9`` matrix instead of a
``6\times6`` one, and every `proj_tens` call on a `Dual` died with a
`DimensionMismatch`. Any new predicate should use `ApproxType`, not
`AbstractFloat`.

**Symbolic assumptions.** `symbols("θ", real = true)` is not enough for SymPy to
reduce ``|\sin\theta|``; a test on a sphere will then compare against expressions
carrying ``\sin\theta/|\sin\theta|``. Declare the tightest assumptions available
and use the `rules` mechanism for the rest.

## Doctests

Docstring examples that print **deterministic** output — types, booleans, small
integers, exact tuples — are written as `jldoctest` and run by

```bash
julia --project=docs -e '
  using TensND, Documenter, LinearAlgebra, SymPy, Tensors, OMEinsum, Rotations
  DocMeta.setdocmeta!(TensND, :DocTestSetup,
      :(using TensND, LinearAlgebra, SymPy, Tensors, OMEinsum, Rotations); recursive = true)
  doctest(TensND)'
```

Examples printing `Float64` component dumps or symbolic expressions stay plain
` ```julia ` blocks: both are fragile across Julia and dependency versions, and
a doctest that has to be re-blessed on every upgrade trains people to re-bless
it without reading. Those examples are still executed — as `@example` blocks —
on the manual and tutorial pages.
