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
| [Projection](projection.md) | `proj_tens`, `best_sym_tens`, predicates |
| [Special tensors](special.md) | Levi-Civita, rotations, unit vectors |
| [Coordinate systems](coorsystems.md) | `CoordinateSystem`, its exact geometry and the operators |
| [Numerical coordinate systems](coorsystems_num.md) | `PointwiseGeometry` |
| [Submanifolds](submanifold.md) | `Submanifold` and its accessors |
| [Symbolic helpers](symbolic.md) | `tsimplify`, `tsubs`, and the array algebra |
| [Full index](full_index.md) | everything, alphabetically |

## [Renamed in 0.8](@id api-renamed-0-8)

TensND 0.8 names its functions as Julia does: lowercase, a proper name included,
as `cholesky` or `schur` in LinearAlgebra; types keep their capitals. The former
names still work and warn once, naming their successor; they will be removed in
a later minor release. The moduli keep their symbols as exact aliases, without a
warning.

| Former name | Name in 0.8 |
| :--- | :--- |
| `KM`, `inv_KM` | [`mandel`](@ref), [`from_mandel`](@ref) |
| `KM_material`, `KM_rotation` | [`mandel_material`](@ref), [`mandel_rotation`](@ref) |
| `KM_from_ti_params`, `KM_from_ti8_params`, `KM_from_ortho_params`, `KM_from_cubic_params` | `mandel_from_ti_params`, `mandel_from_ti8_params`, `mandel_from_ortho_params`, `mandel_from_cubic_params` |
| `ti_params_from_KM`, `ti8_params_from_KM`, `ortho_params_from_KM`, `cubic_params_from_KM` | `ti_params_from_mandel`, `ti8_params_from_mandel`, `ortho_params_from_mandel`, `cubic_params_from_mandel` |
| `KM_IP`, `KM_OP` (constants) | [`MANDEL_IP`](@ref), [`MANDEL_OP`](@ref) |
| `km_ip_block`, `km_op_block`, `km_ip_embed`, `km_op_embed`, `km_from_blocks`, `km_inv_schur` | `mandel_ip_block`, `mandel_op_block`, `mandel_ip_embed`, `mandel_op_embed`, `mandel_from_blocks`, `mandel_inv_schur` |
| `is_ISO`, `is_TI`, `is_ORTHO`, `is_CUBIC` | [`isisotropic`](@ref), [`istransverselyisotropic`](@ref), [`isorthotropic`](@ref), [`iscubic`](@ref) |
| `tens_TI`, `tens_TI_eng`, `tens_TI_Hoenig` | [`tens_ti`](@ref), [`tens_ti_eng`](@ref), [`tens_ti_hoenig`](@ref) |
| `arg_TI`, `arg_TI_eng`, `arg_TI_Hoenig` | [`ti_args`](@ref), [`ti_args_eng`](@ref), [`ti_args_hoenig`](@ref) |
| `fromISO` | [`iso_to_ti`](@ref) |
| `ISO`, `Walpole` | [`iso_projectors`](@ref), [`walpole_basis`](@ref) (or [`walpole_basis_sym`](@ref) for `sym = true`) |
| `LeviCivita` | [`levi_civita`](@ref) |
| `Christoffel`, `Lame`, `getOM` | [`christoffel`](@ref), [`lame_coefficients`](@ref), [`position_vector`](@ref) |
| `GRAD`, `SYMGRAD`, `DIV`, `LAPLACE`, `HESS` | [`grad`](@ref), [`symgrad`](@ref), [`divergence`](@ref), [`laplacian`](@ref), [`hess`](@ref) |
| `CoorSystemSym`, `CoorSystemNum` (types and constructors) | [`CoordinateSystem`](@ref), of an [`ExactGeometry`](@ref) or a [`PointwiseGeometry`](@ref) |
| `coorsys_cartesian_num`, `coorsys_polar_num`, `coorsys_cylindrical_num`, `coorsys_spherical_num` | `coorsys_cartesian(PointwiseGeometry)`, … (with an optional element type, `coorsys_polar(PointwiseGeometry, Float32)`) |
| `SubManifoldSym` | [`Submanifold`](@ref) |
| `E_nu`, `k_mu`, `mu_nu`, `iso_stiffness_E_nu` | [`young_poisson`](@ref), [`bulk_shear`](@ref), [`shear_poisson`](@ref), [`iso_stiffness_young_poisson`](@ref); the former names stay as aliases |

`divergence` is Tensors.jl's function, extended to the coordinate systems, so
that `using TensND, Tensors` leaves one `divergence`. The Hessian is `hess`, not
`hessian`: SymPy, which TensND loads, exports a `hessian` of its own, and the
same name would have made `using TensND, SymPyPythonCall` ambiguous.
