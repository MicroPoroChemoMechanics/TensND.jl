```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: "TensND.jl"
  text: "Tensors on any coordinate system"
  tagline: Bases and variance, structured tensor types, symmetry projection, and differential operators in curvilinear coordinates — symbolic or numeric, one implementation.
  image:
    src: /logo.png
    alt: TensND
  actions:
    - theme: brand
      text: Get started
      link: /manual/installation
    - theme: alt
      text: Theory
      link: /theory/
    - theme: alt
      text: Tutorials
      link: /tutorials/
    - theme: alt
      text: View on GitHub
      link: https://github.com/MicroPoroChemoMechanics/TensND.jl

features:
  - icon: 📐
    title: Theory

    details: Bases and variance, the tensor products, Kelvin-Mandel storage, rotations, and the symmetry classes with their Walpole algebra.
    link: /theory/
  - icon: 🧰
    title: Manual
    details: Every type and operator with the call that produces it — bases, tensors, structured types, projection, coordinate systems, submanifolds.
    link: /manual/getting_started
  - icon: 🎓
    title: Tutorials
    details: Worked sequences, from tensor algebra to differential calculus and applications in mechanics, each runnable end to end.
    link: /tutorials/
  - icon: ⚡
    title: Structured types
    details: TensISO, TensTI, TensCubic and TensOrtho store 2, 5, 3 and 9 scalars and compute products and inverses in closed form, orders of magnitude faster than the dense route.
    link: /manual/structured_tensors
  - icon: 📖
    title: API reference
    details: Every exported type and function, grouped by topic.
    link: /api/
  - icon: 🛠️
    title: Developer guide
    details: How the package is put together, and how to add a coordinate system without touching the rest.
    link: /developer/architecture
---
```

## Installation

TensND.jl is registered in Julia's General registry. In the Julia REPL, press
`]` for the Pkg mode and add it:

```julia
pkg> add TensND
```

then load it:

```julia
using TensND
```

Nothing else is installed by hand. The symbolic backend, SymPy, is reached
through [SymPyPythonCall.jl](https://github.com/jverzani/SymPyPythonCall.jl):
on the first `using TensND`, [CondaPkg.jl](https://github.com/JuliaPy/CondaPkg.jl)
installs a private Python with the `sympy` module into the `.CondaPkg` folder of
the active environment. This needs an internet connection once; the next loads
are offline.

- **Offline afterwards, or behind a proxy that blocks it:** run the first load
  once where the connection works, or set `JULIA_CONDAPKG_OFFLINE=true` once the
  environment exists.
- **An existing Python of your own** (it must have `sympy`): set
  `JULIA_CONDAPKG_BACKEND=Null` and `JULIA_PYTHONCALL_EXE=/path/to/python`
  before loading TensND, or store them as preferences
  (`pkg> preference add CondaPkg backend=Null`,
  `pkg> preference add PythonCall exe=/path/to/python`). See the
  [CondaPkg](https://github.com/JuliaPy/CondaPkg.jl#preferences) and
  [PythonCall](https://juliapy.github.io/PythonCall.jl/stable/pythoncall/#pythoncall-config)
  documentation.

**Coming from SymPy.jl:** write `using SymPyPythonCall` instead of
`using SymPy`. The symbols of SymPy.jl (`Sym{PyObject}`) are another type than
those of SymPyPythonCall (`Sym{Py}`), which TensND uses.

## What it does

`TensND` handles tensors of any order in any coordinate system. A basis carries
its metric, so covariant and contravariant components are related rather than
assumed; the differential operators follow from the chart, symbolically or by
automatic differentiation.

Three structured types — [`TensISO`](@ref), [`TensTI`](@ref) and
[`TensCubic`](@ref), [`TensOrtho`](@ref) — store the 2, 5, 3 and 9 scalars a
symmetry class really has,
and compute products and inverses in closed form. The same code runs on
`Float64`, `ForwardDiff.Dual`, `SymPyPythonCall.Sym` and `Symbolics.Num`.

The design is inspired by the Maple library
[Tens3d](http://jean.garrigues.perso.centrale-marseille.fr/tens3d.html) of Jean
Garrigues.

```@example home
using TensND, SymPyPythonCall

Spherical = coorsys_spherical()
θ, ϕ, r = getcoords(Spherical)
𝐞ᶿ, 𝐞ᵠ, 𝐞ʳ = unitvec(Spherical)
@set_coorsys Spherical

σʳʳ = SymFunction("σʳʳ", real = true)(r)
σᶿᶿ = SymFunction("σᶿᶿ", real = true)(r)
𝛔 = σʳʳ * 𝐞ʳ ⊗ 𝐞ʳ + σᶿᶿ * (𝐞ᶿ ⊗ 𝐞ᶿ + 𝐞ᵠ ⊗ 𝐞ᵠ)

pprint(divergence(𝛔))
```

The equilibrium equation of a spherically symmetric stress state, derived rather
than transcribed.
