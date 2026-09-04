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
ISO
iso_projectors
```

`isotropify` is the exact rotation-group average, not a form of compact
storage, so it is documented with the other averages under
[Exact rotation-group averages](@ref api-projection).
