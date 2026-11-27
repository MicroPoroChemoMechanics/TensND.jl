# [Tensors](@id api-tensors)

The tensor hierarchy, component conversion and the Kelvin-Mandel map. Theory:
[Tensor algebra](@ref th-tensor-algebra); usage: [Tensors](@ref man-tensors).

```@docs
AbstractTens
Tens
static_tens
symmetric_tens
components
canonical_components
express
canonical
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
mandel_pairs
mandel_length
mandel_dimension
mandel_weights
diff_with_basis
tens_basis
```

## The matrix of a linear map

A tensor of order 4 with the minor symmetries, or of order 2, as the matrix of
the linear map it is, and what that matrix tells: the sign of the quadratic
form and the departure from the major symmetry.

```@docs
operator_matrix
from_operator_matrix
definiteness
major_asymmetry
```

## The two-index (Voigt) matrix

The ``6\times6`` matrices of a stiffness and of a compliance in the two-index
notation, pairs in the order ``11, 22, 33, 23, 13, 12``; usage:
[Tensors](@ref man-tensors).

```@docs
voigt_stiffness
voigt_compliance
stiffness_from_voigt
compliance_from_voigt
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
TensND.Tensors.minorsymmetric(::TensND.AbstractTens{4})
TensND.SymmetricFirstPair
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

## Symmetric tensors as Kelvin-Mandel vectors

A symmetric second-order tensor of dimension 3 as the static vector of its six
Kelvin-Mandel components, with the operations a constitutive law needs on it.

```@docs
Mandel
identity_mandel
trace_mandel
deviator
sqnorm_deviator
von_mises
equivalent_strain
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

## Small matrices by cofactors

Closed-form inverses of 2×2 and 3×3 matrices, the adjugate of a symmetric 3×3
one, and the acoustic tensor as the numbers an integral over the directions
needs: no pivot, no allocation, duals and symbols through.

```@docs
inv3
sym2_inv
sym3_inv
sym3_adj
acoustic_components
```

