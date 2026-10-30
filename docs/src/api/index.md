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
| [Submanifolds](submanifold.md) | `Submanifold` and its accessors |
| [Symbolic helpers](symbolic.md) | `tsimplify`, `tsubs`, and the array algebra |
| [Full index](full_index.md) | everything, alphabetically |

## [Renamed in 0.11](@id api-renamed-0-11)

TensND 0.11 names each operation after what it does. The former names are
removed, without a deprecation period, as are those 0.8 had
deprecated (`KM`, `tens_TI`, `Walpole`, `coorsys_*_num`, …: the
[changelog](https://github.com/MicroPoroChemoMechanics/TensND.jl/blob/main/CHANGELOG.md)
lists them with their successors). The moduli keep their symbols as exact
aliases: `E_nu`, `k_mu`, `mu_nu` and `iso_stiffness_E_nu` are
[`young_poisson`](@ref), [`bulk_shear`](@ref), [`shear_poisson`](@ref) and
[`iso_stiffness_young_poisson`](@ref).

| Name in 0.10 | Name in 0.11 |
| :--- | :--- |
| `change_tens(t, ℬ[, var])`, `change_tens(t, var)` | [`express`](@ref)`(t, ℬ[, var])`, `express(t, var)` |
| `change_tens_canon(t)`, `components_canon(t)` | [`canonical`](@ref)`(t)`, [`canonical_components`](@ref)`(t)` |
| `mandel66_minor(A)`, `mandel33_minor(A)`, `mandel66(A)` | [`mandel`](@ref)`(A; symmetrize = true)` |
| `array_from_mandel66(M)`, `array_from_mandel33(M)` | [`from_mandel`](@ref)`(M)` (its components: `get_array(from_mandel(M))`) |
| `tens_ti(C…, n)`, `ti_args(t)` | [`ti_from_constants`](@ref)`(C…, n)`, [`ti_constants`](@ref)`(t)` |
| `tens_ti_eng(…)`, `ti_args_eng(t)` | `ti_from_constants(…; form = :engineering)`, `ti_constants(t; form = :engineering)` |
| `tens_ti_hoenig(…)`, `ti_args_hoenig(t)` | `ti_from_constants(…; form = :hoenig)`, `ti_constants(t; form = :hoenig)` |
| `tens_cubic(C₁₁, C₁₂, C₄₄, frame)`, `arg_cubic(t)` | [`cubic_from_constants`](@ref), [`cubic_constants`](@ref) |
| `get_ℓ(t)`, `get_ℓ8(t)` | [`walpole_coeffs`](@ref)`(t)`, `walpole_coeffs(t, Val(8))` |
| `getcoords` | [`get_coords`](@ref) |
| `rot2(θ)`, `rot3(θ, ϕ, ψ)`, `rot6(θ, ϕ, ψ)` | [`rotation_matrix`](@ref)`(θ)`, `rotation_matrix(θ, ϕ, ψ)`, [`rotation_tensor4`](@ref)`(θ, ϕ, ψ)` |
| `proj_tens(:TI, A, n)` | [`project`](@ref)`(A, :TI, n)` |
| `ti_average_mandel66`, `iso_average_mandel66`, `ti_average_mat33` | [`ti_average_mandel`](@ref), [`iso_average_mandel`](@ref), [`ti_average_matrix`](@ref) |

`divergence` is Tensors.jl's function, extended to the coordinate systems, so
that `using TensND, Tensors` leaves one `divergence`. The Hessian is `hess`, not
`hessian`: SymPy, which TensND loads, exports a `hessian` of its own, and the
same name would have made `using TensND, SymPyPythonCall` ambiguous.
