# [Architecture](@id dev-architecture)

How the source is organized, and why. Files are listed in the order
`src/TensND.jl` includes them, which is also their dependency order.

## The include chain

```mermaid
%%{init: {"flowchart": {"useMaxWidth": false}}}%%
flowchart TB
    AU["<b>array_utils.jl</b><br/>products and contractions<br/>on raw arrays (OMEinsum)"]
    BA["<b>bases.jl</b><br/>Basis hierarchy,<br/>metric, variance"]
    TE["<b>tens.jl</b><br/>AbstractTens, Tens,<br/>components, KM"]
    ISO["<b>tens_isotropic.jl</b><br/>TensISO"]
    ANI["<b>tens_anisotropic.jl</b><br/>TensTI, TensOrtho"]
    OPS["<b>structured_tens_ops.jl</b><br/>shared metaprogramming"]
    PRO["<b>structured_tens_promotion.jl</b><br/>cross-class promotion"]
    PRJ["<b>tens_projection.jl</b><br/>projection kernels"]
    SPE["<b>special_tens.jl</b><br/>LeviCivita, rot2/3/6, init_*"]
    CSS["<b>coorsystems.jl</b><br/>CoorSystemSym"]
    CSN["<b>coorsystems_num.jl</b><br/>CoorSystemNum"]
    SUB["<b>submanifold.jl</b><br/>SubManifoldSym"]
    AU --> BA --> TE --> ISO --> ANI --> OPS --> PRO --> PRJ --> SPE --> CSS --> CSN --> SUB
```

## Layers

**Raw arrays.** `array_utils.jl` implements every product as an OMEinsum
contraction code on `AbstractArray`s. It knows nothing about bases, and is
therefore generic in the element type by construction. Everything above it
inherits that genericity.

**Bases and tensors.** `bases.jl` and `tens.jl` add the basis and the variance.
The tensor type follows the basis type, which is what makes the cost of a
component conversion a property of the data rather than of the call site.

**Structured types.** `tens_isotropic.jl` and `tens_anisotropic.jl` store
coefficients instead of components and override the generic operations with
closed forms. `structured_tens_ops.jl` holds the metaprogramming they share —
every structured type exposes `get_data` and a `_rebuild`, and the shared code
generates the scalar arithmetic and the symbolic passes from those two.
`structured_tens_promotion.jl` handles mixed-class operations.

**Projection.** `tens_projection.jl` works entirely in the ``6\times6``
Kelvin–Mandel picture. The orientation search lives in the weak-dependency
extension `ext/TensNDNLoptExt.jl`, so `NLopt` is optional.

**Differential calculus.** The symbolic side is written **once**. `∂` and the
five operators are defined in `coorsystems.jl` over `AbstractCoorSystem`, and
the only thing a submanifold changes is how many directions it differentiates
along — the [`nderiv`](@ref) trait, `dim` for a chart and `dim-1` for a surface.
`submanifold.jl` therefore contains geometry and nothing else; it used to carry
a verbatim copy of all six definitions, which is how one bug fix had to be
applied twice and one variance bug survived in the copy.

Both types hold the same internal `ChartCore` — position vector, coordinates, the
two bases, the natural and dual frames, Lamé coefficients, connection array and
simplification settings — and a submanifold adds only the two fundamental forms.

`coorsystems_num.jl` keeps its own implementation, deliberately: its operators
take a *field* and return a *function of the point*, a different evaluation
model that a shared body would obscure rather than simplify.

## Two design decisions worth knowing

**Geometry and data have independent element types.** A `TensOrtho` carries its
material frame as a separate type parameter, so `Dual`-valued constants can sit
on a `Float64` frame. Forcing them to agree broke every attempt to differentiate
a scheme with respect to a modulus. The same reasoning applies to the symmetry
axis of a `TensTI` and to the frame passed to `proj_tens`, where converting the
frame to the field's element type once produced `NaN` derivatives — see the
comment in `proj_tens(::Val{:ORTHO}, …)`.

**Closed forms only while the result stays in the class.** The structured
methods return the generic tensor whenever the exact result leaves the class —
a product of two orthotropic tensors needs twelve constants, not nine. The
fallback is bit-for-bit the generic computation, never an approximation. The
full table is on [Structured tensors](@ref man-structured).

## The extension

`ext/TensNDNLoptExt.jl` overrides `TensND._proj_tens_opt`, the hook the
no-orientation methods of [`proj_tens`](@ref) dispatch to. Without `NLopt` those
methods throw a clear error; with it, they run the deterministic multi-start
described on [Projection onto a symmetry class](@ref th-projection).

Because loading `NLopt` **permanently** activates the extension for the session,
`test/test_nlopt_ext.jl` must run last: `test/test_tens_projection.jl` asserts
the no-NLopt behavior and would fail after it.
