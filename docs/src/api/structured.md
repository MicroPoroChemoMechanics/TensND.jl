# [Structured tensors — isotropic](@id api-structured)

Compact storage for isotropic tensors, in arbitrary dimension. Theory:
[Isotropic tensors](@ref th-isotropic); usage:
[Structured tensors](@ref man-structured).

```@docs
TensISO
tens_Id2
tens_Id4
tens_J4
tens_K4
iso_projectors
```

## [Blocks of a symmetry class](@id api-symmetry-blocks)

The small matrices whose products and inverses are those of the tensors of a
class; usage on [Blocks: a class as a sum of small matrices](@ref man-symmetry-blocks).

```@docs
symmetry_blocks
from_symmetry_blocks
symmetry_basis
symmetry_layout
symmetry_matrix
```

## [Isotropic elastic moduli](@id api-iso-moduli)

The stiffness and the compliance from any pair of the six moduli, every modulus
read from either (`role = :stiffness` or `:compliance`), and plane strain and
plane stress are explained on [Elastic moduli](@ref man-elastic-moduli).

```@docs
iso_stiffness
iso_stiffness_bulk_young
iso_stiffness_bulk_lame
iso_stiffness_bulk_poisson
iso_stiffness_bulk_pwave
iso_stiffness_young_lame
iso_stiffness_young_shear
iso_stiffness_young_poisson
iso_stiffness_young_pwave
iso_stiffness_lame_shear
iso_stiffness_lame_poisson
iso_stiffness_lame_pwave
iso_stiffness_shear_poisson
iso_stiffness_shear_pwave
iso_stiffness_poisson_pwave
iso_compliance
iso_compliance_bulk_young
iso_compliance_bulk_lame
iso_compliance_bulk_poisson
iso_compliance_bulk_pwave
iso_compliance_young_lame
iso_compliance_young_shear
iso_compliance_young_poisson
iso_compliance_young_pwave
iso_compliance_lame_shear
iso_compliance_lame_poisson
iso_compliance_lame_pwave
iso_compliance_shear_poisson
iso_compliance_shear_pwave
iso_compliance_poisson_pwave
lame_lambda
bulk_modulus
young_modulus
shear_modulus
poisson_ratio
pwave_modulus
bulk_young
bulk_lame
bulk_shear
bulk_poisson
bulk_pwave
young_lame
young_shear
young_poisson
young_pwave
lame_shear
lame_poisson
lame_pwave
shear_poisson
shear_pwave
poisson_pwave
plane_strain_stiffness
plane_stress_stiffness
plane_stress_to_plane_strain
plane_strain_to_plane_stress
```

`isotropify` is the exact rotation-group average, not a form of compact
storage, so it is documented with the other averages under
[Exact rotation-group averages](@ref api-projection).
