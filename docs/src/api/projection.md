# [Projection](@id api-projection)

Closest tensor of a prescribed material symmetry, and symmetry detection.
Theory: [Projection onto a symmetry class](@ref th-projection); usage:
[Projection](@ref man-projection).

The orientation-optimizing methods of `proj_tens` require `NLopt`
(package extension `TensNDNLoptExt`).

```@docs
proj_tens
best_sym_tens
is_ISO
is_TI
is_ORTHO
ti_params_from_KM
KM_from_ti_params
ortho_params_from_KM
KM_from_ortho_params
```
