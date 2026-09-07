# [Cubic symmetry](@id api-cubic)

Compact storage for cubic tensors: three constants and a cube frame. Theory:
[Cubic symmetry](@ref th-cubic); usage:
[Structured tensors](@ref man-structured).

The class is the one place in this package where the algebra is **closed**: the
double contraction of two cubic tensors about the same cube is cubic, and every
cubic tensor is automatically major-symmetric. Both follow from the three
irreducible representations of the octahedral group appearing with multiplicity
one.

```@docs
TensCubic
tens_cubic
arg_cubic
cubic_anisotropy
iso_to_cubic
cubic_to_ortho
```

`best_fit_cubic`, `cubic_params_from_KM` and `KM_from_cubic_params` are the
projection family and live with their siblings under
[Projection](@ref api-projection). `isotropify` of a cubic tensor is an exact
rotation-group average, not a fit, and is documented there too.
