# [Elastic waves](@id api-waves)

Plane waves of a homogeneous elastic medium, their first-order forms in a
weakly anisotropic one, waves at a plane interface, and the small waves of a
deformed solid. Theory: [Plane waves](@ref th-plane-waves),
[Thomsen's parameters](@ref th-thomsen),
[Weak anisotropy of any symmetry](@ref th-weak-anisotropy),
[Shear-wave splitting](@ref th-splitting),
[Waves at an interface](@ref th-interface-waves),
[Characteristic roots](@ref th-characteristic-roots),
[Third-order elasticity](@ref th-third-order); usage:
[Plane waves](@ref man-plane-waves),
[Waves at an interface](@ref man-interface-waves),
[Third-order elasticity](@ref man-third-order).

## Plane waves

```@docs
PlaneWaves
plane_waves
phase_velocities
slowness
wave_modes
christoffel_invariants
```

## Thomsen's parameters and weak anisotropy

The form `:thomsen` of [`ti_from_constants`](@ref) and [`ti_constants`](@ref)
builds a stiffness from these parameters and reads them back.

```@docs
thomsen_parameters
anellipticity
weak_wave_modes
squared_moveout_velocities
weak_plane_waves
generalized_thomsen_parameters
shear_splitting
```

## Waves at an interface

```@docs
InterfaceWave
interface_waves
reflection_transmission
free_surface_reflection
weak_contrast_reflection
```

## Characteristic roots of a stiffness of the plane

```@docs
engineering_compliances
characteristic_pair
characteristic_roots
```

## Third-order elasticity

```@docs
ThirdOrderIsotropic
third_order_constants
third_order_from_constants
strain_energy
cauchy_stress
tangent_stiffness
acoustoelastic_tensor
prestressed_waves
acoustoelastic_coefficients
strain_stiffening
```
