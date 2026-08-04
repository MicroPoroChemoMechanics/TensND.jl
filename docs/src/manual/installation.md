# [Installation](@id man-installation)

`TensND` is registered in the Julia General registry.

```julia
julia> import Pkg; Pkg.add("TensND")
```

or, from the Pkg REPL mode (`]`):

```julia
pkg> add TensND
```

## Optional dependency

The **orientation search** of [`proj_tens`](@ref) and [`best_sym_tens`](@ref) —
finding the transversely isotropic axis or the orthotropic frame when it is not
known in advance — is provided by a package extension:

```julia
pkg> add NLopt
```

Loading `NLopt` activates `TensNDNLoptExt` and enables the methods of
[`proj_tens`](@ref) that take no axis or frame:

```julia
using TensND, NLopt

B, d, drel = proj_tens(:TI, A)        # axis optimized
B, d, drel = proj_tens(:ORTHO, A)     # frame optimized
```

Without `NLopt`, everything else works, including projection onto a class with a
**given** axis or frame and the `optimize_angles = false` path of
[`best_sym_tens`](@ref), which infers the orientation from the Kelvin–Mandel
eigenstructure. See [Projection](@ref man-projection).

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

DIV(𝐞ʳ)          # 2/r
LAPLACE(1 / r)   # 0
```

If those two return `2/r` and `0`, the symbolic stack is working.
