# [Tutorials](@id tut-index)

Every page in this section is generated from a runnable script in `scripts/` by
[Literate.jl](https://github.com/fredrikekre/Literate.jl), so each one can be
executed, edited and re-run locally:

```bash
julia scripts/03_walpole.jl
```

The generated pages come with a pre-run Jupyter notebook and a cleaned
standalone script under `docs/generated_notebooks/` and
`docs/generated_scripts/`.

Simplest first within each group.

## Tensor algebra

| Page | What it shows |
| :--- | :--- |
| [Bases, variance and the metric](generated/bases_variance.md) | an oblique basis, dual basis and metric, covariant vs contravariant components, and where variance stops mattering |
| [Tensor products, contractions and their identities](generated/tensor_algebra.md) | the ⊠/⊠ˢ algebra checked numerically, the identity that **fails**, isotropic projection, and why it does not commute with inversion |
| [The Walpole basis end to end](generated/walpole.md) | the six ``\mathbb{W}_i`` as ``6\times6`` matrices, the multiplication table, the correct decomposition of ``\mathbb{I},\mathbb{J},\mathbb{K}``, and three parametrizations of one carbon/epoxy ply |
| [The eight-dimensional axially invariant space](generated/walpole_extended.md) | ``\mathbb{W}_7,\mathbb{W}_8``, the complex product rule, and what a five-parameter best fit silently discards |
| [Projection and symmetry detection](generated/projection.md) | fixed axis vs optimized, determinism, the detection cascade, and the cost of each path |

## Differential calculus

| Page | What it shows |
| :--- | :--- |
| [Differential operators, symbolically](generated/symbolic_operators.md) | the five operators in polar, cylindrical and spherical coordinates; the hoop strain that comes from the connection alone; the index-placement convention |
| [Christoffel symbols](generated/christoffel.md) | the definition checked term by term, every predefined system tabulated, symbolic and numerical routes compared, spheroidal harmonics |
| [Differential operators, numerically](generated/numerical_operators.md) | pointwise evaluation by AD, the Lamé hollow-sphere problem with its equilibrium residual, and sensitivities through the operators |
| [Surfaces: fundamental forms and curvature](generated/submanifolds.md) | sphere, cylinder, plane and paraboloid; Gauss's *Theorema Egregium*; the array-comparison trap |
| [Classical identities of the differential operators](generated/operator_identities.md) | every identity **checked**, not quoted, on four charts — including a deliberately non-orthogonal one, the only case that tests the covariant/contravariant bookkeeping |

## Applications in mechanics

| Page | What it shows |
| :--- | :--- |
| [Elastic Green's functions](generated/green_function.md) | the Kelvin solution in 2-D and 3-D, ``\mathbb{\Gamma}=-\mathrm{HESS}(\boldsymbol{G})`` matched against its closed form |
| [A two-point Green tensor between offset frames](generated/cluster.md) | combining tensors that live in two different spherical charts |
| [The elastic sphere, solved symbolically](generated/sphere_problems.md) | hydrostatic and deviatoric problems, the Lamé exponents, and how they assemble into an ``N``-layer sphere |
| [Acoustic and Hill tensors](generated/hill_tensors.md) | ``\boldsymbol{K}=\underline{\xi}\cdot\mathbb{C}\cdot\underline{\xi}`` and ``\mathbb{\Lambda}=\mathbb{C}:\mathbb{\Gamma}:\mathbb{C}``, written once and run on an isotropic *and* a transversely isotropic stiffness |

## Types, differentiation and performance

| Page | What it shows |
| :--- | :--- |
| [Differentiation and interoperability](generated/ad_interop.md) | `ForwardDiff` through the algebra and through a projection, symbolic vs automatic derivatives, `Tensors.jl` interop |
| [Performance of the structured types](generated/performance.md) | measured speed-ups of the closed forms over the dense route, allocations, and where the structure stops paying |

## Where to start

If you are new to the library, read
[Getting started](@ref man-getting-started) first, then
[Bases, variance and the metric](generated/bases_variance.md) — variance is the
one notion a Cartesian-only tensor library does not have, and everything else
rests on it.
