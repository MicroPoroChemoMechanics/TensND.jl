```@meta
CurrentModule = TensND
```

# Documentation for [TensND](https://codeberg.org/MicroPoroChemoMechanics/TensND.jl)

![TensND.jl](assets/logo.svg)

*Symbolic and numerical tensor calculations in arbitrary coordinate systems.*

## Introduction

TensND.jl is a Julia package for tensor calculations of any order and dimension in arbitrary
coordinate systems. It supports both **symbolic computation** (via SymPy and Symbolics.jl) and
**numerical evaluation** (via ForwardDiff automatic differentiation).

### Key features

- **Basis types**: canonical, rotated, orthogonal, and fully general non-orthogonal bases
- **Tensor algebra**: products, contractions, change of basis, variance (covariant/contravariant) management
- **Structured tensor types**: isotropic (`TensISO`), transversely isotropic (`TensTI{4}`, `TensTI`), orthotropic (`TensOrtho`) with compact parametric storage and efficient algebra
- **Symmetry projection**: find the closest ISO, TI, or ORTHO tensor via Frobenius distance minimization; rotation-optimized via NLopt
- **Differential operators**: gradient, symmetric gradient, divergence, Laplacian, Hessian in curvilinear coordinate systems (symbolic and numerical)
- **Generic type system**: works with `Float64`, symbolic types (`Sym`, `Num`), and `ForwardDiff.Dual` numbers for automatic differentiation

The implementation is inspired by the Maple library [Tens3d](http://jean.garrigues.perso.centrale-marseille.fr/tens3d.html) developed by Jean Garrigues.

## Installation

The package can be installed with the Julia package manager. From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```julia
pkg> add TensND
```

Or, equivalently, via the `Pkg` API:

```julia
julia> import Pkg; Pkg.add("TensND")
```

For rotation-optimized tensor projections, also install the optional dependency:

```julia
pkg> add NLopt
```

## Manual outline

```@contents
Pages = [
    "man/getting_started.md",
    "man/bases.md",
    "man/tensors.md",
    "man/coorsystems.md",
]
Depth = 1
```

## Tutorials

```@contents
Pages = [
    "tuto/nlayersphere.md",
    "tuto/coorsystems_num.md",
    "tuto/green_function.md",
    "tuto/projection.md",
]
Depth = 1
```

## Citing TensND.jl

```latex
@misc{TensND.jl,
  author  = {Jean-François Barthélémy},
  title   = {TensND.jl},
  url     = {https://codeberg.org/MicroPoroChemoMechanics/TensND.jl},
  version = {v0.1.1},
  year    = {2021},
  month   = {8}
}
```

## Related packages

- [SymPy.jl](https://github.com/JuliaPy/SymPy.jl) — symbolic mathematics via Python/PyCall
- [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl) — native Julia CAS (also supported)
- [Tensors.jl](https://github.com/Ferrite-FEM/Tensors.jl) — low-level tensor storage
- [OMEinsum.jl](https://github.com/under-Peter/OMEinsum.jl) — Einstein summation
- [Rotations.jl](https://github.com/JuliaGeometry/Rotations.jl) — rotation matrices
- [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) — automatic differentiation
- [NLopt.jl](https://github.com/JuliaOpt/NLopt.jl) — rotation-optimized projections (optional)

## References

1. [Tens3d](http://jean.garrigues.perso.centrale-marseille.fr/tens3d.html) — Maple tensor library by Jean Garrigues
1. [Walpole basis decomposition](https://sbrisard.github.io/posts/20140226-decomposition_of_transverse_isotropic_fourth-rank_tensors.html) — blog post by S. Brisard
1. Hoenig, A. (1978). *The behavior of a flat elliptical crack in an anisotropic elastic body*. International Journal of Solids and Structures, 14, 925-934.
