```@meta
CurrentModule = TensND
```

# Documentation for [TensND](https://github.com/jfbarthelemy/TensND.jl)

*Package allowing tensor calculations in arbitrary coordinate systems.*

## Introduction

This Julia package provides tools to perform tensor calculations of any order and any dimension in arbitrary coordinate systems (cartesian, polar, cylindrical, spherical, spheroidal or any user defined coordinate systems...). In particular differential operators are available: gradient, symmetrized gradient, divergence, Laplace, Hessian. The implementation of this library is much inspired by the Maple library [Tens3d](http://jean.garrigues.perso.centrale-marseille.fr/tens3d.html) developed by Jean Garrigues.

## Installation

The package can be installed with the Julia package manager. From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```julia
pkg> add TensND
```

Or, equivalently, via the `Pkg` API:

```julia
julia> import Pkg; Pkg.add("TensND")
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
  url     = {https://github.com/jfbarthelemy/TensND.jl},
  version = {v0.1.1},
  year    = {2021},
  month   = {8}
}
```

## Related packages

- [SymPy.jl](https://github.com/JuliaPy/SymPy.jl)
- [Tensors.jl](https://github.com/Ferrite-FEM/Tensors.jl)
- [OMEinsum.jl](https://github.com/under-Peter/OMEinsum.jl)
- [Rotations.jl](https://github.com/JuliaGeometry/Rotations.jl)

## References

1. [Tens3d](http://jean.garrigues.perso.centrale-marseille.fr/tens3d.html)
1. [Sébastien Brisard's blog](https://sbrisard.github.io/posts/20140226-decomposition_of_transverse_isotropic_fourth-rank_tensors.html)