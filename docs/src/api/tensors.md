# [Tensors](@id api-tensors)

The tensor hierarchy, component conversion and the Kelvin-Mandel map. Theory:
[Tensor algebra](@ref th-tensor-algebra); usage: [Tensors](@ref man-tensors).

```@docs
AbstractTens
Tens
components
components_canon
change_tens
change_tens_canon
generic_tens
get_order
get_basis
get_var
get_data
get_array
sqnorm
pprint
mandel
from_mandel
mandel_rotation
diff_with_basis
tens_basis
```

## Tensors.jl and the Voigt forms

Components in the canonical basis, whatever basis a tensor is written in: the
form a finite-element code (Ferrite among others) and an Abaqus-style interface
expect.

```@docs
to_tensors
from_tensors
voigt_stress
voigt_strain
stress_from_voigt
strain_from_voigt
```

## Symmetric parts

Methods of the `Tensors.jl` functions, which keep the symmetry class of a
structured tensor where the class allows it.

```@docs
TensND.Tensors.majorsymmetric(::TensND.AbstractTens{4})
TensND.Tensors.symmetric(::TensND.AbstractTens{2})
```

## Real and imaginary parts

```@docs
Base.real(::TensND.AbstractTens)
Base.conj(::TensND.AbstractTens)
```

## Coefficients and flat vectors

```@docs
rebuild
param_norms
to_vector
from_vector
```

## Kelvin-Mandel blocks about a normal

In a frame whose third axis is a normal `n`, the in-plane and out-of-plane
symmetric tensors split the Kelvin-Mandel slots; the maps below work on that
partition, with no pivoting and no SVD (duals and symbols go through).

```@docs
MANDEL_IP
MANDEL_OP
mandel_ip_block
mandel_op_block
mandel_ip_embed
mandel_op_embed
mandel_from_blocks
mandel_inv_schur
plane_pinv
plane_pinv2
deviatoric_pinv
acoustic_tensor
compliance_op_block
```

