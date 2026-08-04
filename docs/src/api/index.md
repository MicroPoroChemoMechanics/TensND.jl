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
| [Projection](projection.md) | `proj_tens`, `best_sym_tens`, predicates |
| [Special tensors](special.md) | Levi-Civita, rotations, unit vectors |
| [Coordinate systems](coorsystems.md) | `CoorSystemSym` and the operators |
| [Numerical coordinate systems](coorsystems_num.md) | `CoorSystemNum` |
| [Submanifolds](submanifold.md) | `SubManifoldSym` and its accessors |
| [Symbolic helpers](symbolic.md) | `tsimplify`, `tsubs`, and the array algebra |
| [Full index](full_index.md) | everything, alphabetically |
