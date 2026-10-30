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

!!! tip "What to load for symbolic computations"
    | to work with | add | load |
    | :--- | :--- | :--- |
    | SymPy symbols (`Sym`: `symbols`, `@syms`, `SymFunction`), and the coordinate systems called without coordinates (`coorsys_spherical()`, `init_cartesian()`, …) | `pkg> add SymPyPythonCall` | `using TensND, SymPyPythonCall` |
    | Symbolics variables (`Num`: `@variables`) | `pkg> add Symbolics` | `using TensND, Symbolics` |
    | numbers and derivatives only | nothing | `using TensND` |

`Symbolics` is a dependency of TensND; `SymPyPythonCall` is not, since
version 0.9.0: what calls Python lives in an extension that loading
SymPyPythonCall activates, so that `using TensND` alone no longer loads Python.
Called without coordinates while SymPyPythonCall is not loaded, a coordinate
system says so in its error message.

`SymPyPythonCall` reaches Python through `PythonCall`, and `CondaPkg` installs
the Python module `sympy` on first use, in an environment of its own: no Python
needs to be set up beforehand. To use a Python of your own instead (one that has
`sympy`), set

```julia
ENV["JULIA_CONDAPKG_BACKEND"] = "Null"
ENV["JULIA_PYTHONCALL_EXE"] = "/path/to/python"
```

before loading `SymPyPythonCall`.

Code written for SymPy.jl changes `using SymPy` into `using SymPyPythonCall`; the
symbols of the two packages are different types.

See [Symbolic and numeric](@ref man-symbolic-numeric) for how the four scalar
worlds — `Float64`, `ForwardDiff.Dual`, `Sym` and `Num` — relate.

## Checking the installation

```julia
using TensND

Spherical = coorsys_spherical()
θ, ϕ, r = get_coords(Spherical)
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
| [SymPyPythonCall.jl](https://github.com/jverzani/SymPyPythonCall.jl) | symbolic backend (SymPy through PythonCall) |
| [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl) | native Julia CAS backend |
| [Rotations.jl](https://github.com/JuliaGeometry/Rotations.jl) | rotation representations |
| [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) | automatic differentiation |
