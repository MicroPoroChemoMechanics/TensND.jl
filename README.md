# TensND.jl

*Symbolic and numerical tensor calculations in arbitrary coordinate systems.*

[![Docs - Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://MicMacTools.github.io/TensND.jl/stable/)
[![Docs - Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://MicMacTools.github.io/TensND.jl/dev)
[![Build Status](https://github.com/MicMacTools/TensND.jl/workflows/CI/badge.svg?branch=main)](https://github.com/MicMacTools/TensND.jl/actions/workflows/CI.yml?query=branch%3Amain)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/MicMacTools/TensND.jl/blob/main/LICENSE)
[![code style: runic](https://img.shields.io/badge/code_style-%E1%9A%B1%E1%9A%A2%E1%9A%BE%E1%9B%81%E1%9A%B2-black)](https://github.com/fredrikekre/Runic.jl)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.17985768-blue)](https://doi.org/10.5281/zenodo.17985768)

## Introduction

TensND.jl is a Julia package for tensor calculations of any order and dimension in arbitrary coordinate systems (cartesian, polar, cylindrical, spherical, spheroidal, or user-defined). It supports both **symbolic computation** (via [SymPy.jl](https://github.com/JuliaPy/SymPy.jl) and [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl)) and **numerical evaluation** (via [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) automatic differentiation).

### Key features

- **Basis types**: canonical, rotated, orthogonal, and fully general (non-orthogonal, symbolic)
- **Tensor algebra**: products (`⊗`, `⊗ˢ`, `⊠`, `⊠ˢ`, `⋅`, `⊡`, `⊙`), change of basis, variance management
- **Structured tensors**: isotropic (`TensISO`), transversely isotropic (`TensTI{4}`, `TensTI`), orthotropic (`TensOrtho`) with compact storage and efficient algebra
- **Symmetry projection**: find the closest ISO, TI, or ORTHO tensor; rotation-optimized via [NLopt.jl](https://github.com/JuliaOpt/NLopt.jl)
- **Differential operators**: gradient, symmetric gradient, divergence, Laplacian, Hessian in curvilinear coordinates (symbolic and numerical)
- **Generic type system**: works with `Float64`, symbolic types (`Sym`, `Num`), and `ForwardDiff.Dual` for automatic differentiation

The implementation is inspired by the Maple library [Tens3d](http://jean.garrigues.perso.centrale-marseille.fr/tens3d.html) developed by Jean Garrigues.

The following example is provided to illustrate the purpose of the library

```julia
julia> using SymPy, TensND

julia> Spherical = coorsys_spherical() ; θ, ϕ, r = getcoords(Spherical) ; 𝐞ᶿ, 𝐞ᵠ, 𝐞ʳ = unitvec(Spherical) ;

julia> @set_coorsys Spherical

julia> GRAD(𝐞ʳ) |> intrinsic
(1/r)𝐞ᶿ⊗𝐞ᶿ + (1/r)𝐞ᵠ⊗𝐞ᵠ

julia> DIV(𝐞ʳ ⊗ 𝐞ʳ) |> intrinsic
(2/r)𝐞ʳ

julia> LAPLACE(1/r) |> intrinsic
0

julia> f = SymFunction("f", real = true)
f

julia> DIV(f(r) * 𝐞ʳ ⊗ 𝐞ʳ) |> intrinsic
(Derivative(f(r), r) + 2*f(r)/r)𝐞ʳ

julia> LAPLACE(f(r)) |> intrinsic
              d       
  2         2⋅──(f(r))
 d            dr
───(f(r)) + ──────────
  2             r
dr

julia> for σⁱʲ ∈ ("σʳʳ", "σᶿᶿ", "σᵠᵠ") @eval $(Symbol(σⁱʲ)) = SymFunction($σⁱʲ, real = true)($r) end

julia> 𝛔 = σʳʳ * 𝐞ʳ ⊗ 𝐞ʳ + σᶿᶿ * 𝐞ᶿ ⊗ 𝐞ᶿ + σᵠᵠ * 𝐞ᵠ ⊗ 𝐞ᵠ ; intrinsic(𝛔)
(σᶿᶿ(r))𝐞ᶿ⊗𝐞ᶿ + (σᵠᵠ(r))𝐞ᵠ⊗𝐞ᵠ + (σʳʳ(r))𝐞ʳ⊗𝐞ʳ

julia> div𝛔 = tsimplify(DIV(𝛔)) ; intrinsic(div𝛔)
((-σᵠᵠ(r) + σᶿᶿ(r))/(r*tan(θ)))𝐞ᶿ + ((r*Derivative(σʳʳ(r), r) + 2*σʳʳ(r) - σᵠᵠ(r) - σᶿᶿ(r))/r)𝐞ʳ
```

## Installation

The package can be installed with the Julia package manager. From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```julia
pkg> add TensND
```

Or, equivalently, via the `Pkg` API:

```julia
julia> import Pkg; Pkg.add("TensND")
```

## Documentation

- [**STABLE**](https://MicMacTools.github.io/TensND.jl/stable/) &mdash; **most recently tagged version of the documentation.**
- [**DEV**](https://MicMacTools.github.io/TensND.jl/dev) &mdash; **development version of the documentation.**

## Citation

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.17985768-blue)](https://doi.org/10.5281/zenodo.17985768)

See [CITATION.cff](CITATION.cff) for citation details.

## Acknowledgements

Parts of this codebase were developed with the assistance of Anthropic's
*Claude Code*, under the author's review and numerical validation.
