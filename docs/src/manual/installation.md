# [Installation](@id man-installation)

`TensND` is registered in the Julia General registry.

```julia
julia> import Pkg; Pkg.add("TensND")
```

or, from the Pkg REPL mode (`]`):

```julia
pkg> add TensND
```

## Symbolic backends

`SymPy` and `Symbolics` are ordinary dependencies, so both symbolic element
types are available without further installation. `SymPy` calls into Python
through `PyCall`; if its build fails, the usual remedy is

```julia
ENV["PYTHON"] = ""
import Pkg; Pkg.build("PyCall")
```

which installs a private Conda Python. See
[Symbolic and numeric](@ref man-symbolic-numeric) for how the four scalar
worlds — `Float64`, `ForwardDiff.Dual`, `Sym` and `Num` — relate.

## Checking the installation

```julia
using TensND

Spherical = coorsys_spherical()
θ, ϕ, r = getcoords(Spherical)
𝐞ᶿ, 𝐞ᵠ, 𝐞ʳ = unitvec(Spherical)
@set_coorsys Spherical

divergence(𝐞ʳ)          # 2/r
laplacian(1 / r)   # 0
```

If those two return `2/r` and `0`, the symbolic stack is working.

## Citing

```bibtex
@misc{TensND.jl,
  author  = {Jean-François Barthélémy},
  title   = {TensND.jl: symbolic and numerical tensor calculations
             in arbitrary coordinate systems},
  url     = {https://github.com/MicroPoroChemoMechanics/TensND.jl},
  doi     = {10.5281/zenodo.17985768},
  year    = {2026}
}
```

`CITATION.cff` in the repository root carries the same metadata in a
machine-readable form. Works cited by this documentation are collected on the
[References](@ref references) page.

## Related packages

| Package | Role here |
| :--- | :--- |
| [Tensors.jl](https://github.com/Ferrite-FEM/Tensors.jl) | low-level storage for small tensors |
| [SymPy.jl](https://github.com/JuliaPy/SymPy.jl) | symbolic backend |
| [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl) | native Julia CAS backend |
| [Rotations.jl](https://github.com/JuliaGeometry/Rotations.jl) | rotation representations |
| [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) | automatic differentiation |
