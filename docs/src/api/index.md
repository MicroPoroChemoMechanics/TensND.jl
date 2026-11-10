# [API overview](@id api-index)

The public interface, grouped by theme. Every exported name appears on exactly
one of the pages below; the last one is a catch-all that lists anything not yet
curated, so a newly exported symbol can never go undocumented.

| Page | Covers |
| :--- | :--- |
| [Bases](bases.md) | basis types, metric, accessors, `init_*` |
| [Tensors](tensors.md) | `AbstractTens`, `Tens`, components, Kelvin-Mandel |
| [Structured tensors](structured.md) | `TensISO` and the isotropic projectors |
| [Walpole and orthotropy](walpole.md) | `TensTI`, `TensOrtho`, parametrizations |
| [Cubic symmetry](cubic.md) | `TensCubic`, three constants and a cube frame |
| [Projection](projection.md) | `project`, `best_sym_tens`, predicates |
| [Special tensors](special.md) | Levi-Civita, rotations, unit vectors, random tensors |
| [Coordinate systems](coorsystems.md) | `CoordinateSystem`, its exact geometry and the operators |
| [Numerical coordinate systems](coorsystems_num.md) | `PointwiseGeometry` |
| [Charts and coordinate conversions](charts.md) | the charts, `coords_to_cartesian`, `cartesian_to_coords`, `change_coords` |
| [Spheroidal harmonics](spheroidal_harmonics.md) | Legendre functions, harmonics, re-expansion coefficients |
| [Submanifolds](submanifold.md) | `Submanifold` and its accessors |
| [Symbolic helpers](symbolic.md) | `tsimplify`, `tsubs`, and the array algebra |
| [Full index](full_index.md) | everything, alphabetically |

`divergence` is Tensors.jl's function, extended to the coordinate systems, so
that `using TensND, Tensors` leaves one `divergence`. The Hessian is `hess`, not
`hessian`: SymPy, which TensND loads, exports a `hessian` of its own, and the
same name would have made `using TensND, SymPyPythonCall` ambiguous.
