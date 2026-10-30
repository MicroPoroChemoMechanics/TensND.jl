# [Projection](@id api-projection)

Closest tensor of a prescribed material symmetry, and symmetry detection.
Theory: [Projection onto a symmetry class](@ref th-projection); usage:
[Projection](@ref man-projection).

```@docs
project
best_sym_tens
isisotropic
istransverselyisotropic
isorthotropic
iscubic
ti_params_from_mandel
mandel_from_ti_params
ti2_params_from_matrix
ortho_params_from_mandel
mandel_from_ortho_params
cubic_params_from_mandel
mandel_from_cubic_params
best_fit_iso
best_fit_ti
best_fit_ortho
best_fit_cubic
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
ti8_params_from_mandel
mandel_from_ti8_params
ti_average_mandel
ti_average_matrix
iso_average_mandel
```
