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
best_fit_iso
best_fit_ti
best_fit_ortho
```

## Exact rotation-group averages

A projection is a least-squares **fit** and drops whatever does not fit; a
rotation-group average is **exact** and lossless on the invariant subspace. The
two coincide only when the tensor already belongs to the class. Use the
averages inside a computation — on a concentration or contribution tensor,
which carries no major symmetry — and the fits only to report parameters.

```@docs
isotropify
transverse_isotropify
mandel66_minor
array_from_mandel66
ti8_params_from_KM
KM_from_ti8_params
ti_average_mandel66
iso_average_mandel66
```
