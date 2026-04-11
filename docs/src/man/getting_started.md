# Getting started

## Brief description of the package

The package relies on the definition of

- **bases** which can be of the following types (`T` denotes the scalar type, subtype of `Number`)

  - `CanonicalBasis{dim,T}`: fundamental canonical basis in `Rdim` in which the metric tensor is the second-order identity
  - `RotatedBasis{dim,T}`: orthonormal basis in `Rdim` obtained by rotation of the canonical basis by means of one angle if `dim=2` or three Euler angles if `dim=3`, the metric tensor is again the second-order identity
  - `OrthogonalBasis{dim,T}`: orthogonal basis in `Rdim` obtained from a given orthonormal rotated basis by applying a scaling factor along each unit vector, the metric tensor is then diagonal
  - `Basis{dim,T}`: arbitrary basis not entering the previous cases

- **tensors**

  - a tensor is determined by a set of data (array or synthetic parameters) corresponding to its `order`, a basis and a tuple of variances
  - depending on the type of basis, the type of tensor can be `TensCanonical{order,dim,T,A}`, `TensRotated{order,dim,T,A}`, `TensOrthogonal{order,dim,T,A}` or `Tens{order,dim,T,A}` if the data are stored under the form of an array or a `Tensor` object (see [Tensors.jl](https://github.com/Ferrite-FEM/Tensors.jl))
  - `TensISO{order,dim,T,N}`: isotropic tensor stored as N scalar parameters (1 for order 2, 2 for order 4)
  - `TensWalpole{T,N}`: transversely isotropic 4th-order tensor in the Walpole basis (N=5 for major-symmetric, N=6 for general case)
  - `TensTI{order,T,N}`: transversely isotropic tensor (order 2 or 4)
  - `TensOrtho{T}`: orthotropic 4th-order tensor with 9 elastic constants in a material frame
  - Projection functions `proj_tens` and `best_sym_tens` allow finding the closest tensor of a given symmetry class (ISO, TI, ORTHO) with or without rotation optimization (the latter requires `NLopt`)

- **coordinate systems**

  - **Symbolic** (`CoorSystemSym`): a coordinate system contains all information required to perform differential operations on tensor fields: position vector `OM` expressed in the canonical basis, coordinate names, natural basis, normalized basis, Christoffel coefficients
  - **Numerical** (`CoorSystemNum`): pointwise evaluation of differential operators via automatic differentiation (`ForwardDiff`), without requiring symbolic setup
  - Predefined systems are available for both: cartesian, polar, cylindrical, spherical and spheroidal. The user can also define custom systems

All types support generic scalar types (`Sym`, `Num`, `Float64`, `ForwardDiff.Dual`, ...) enabling seamless interoperability between symbolic computation and automatic differentiation.

## Detailed manual

Before detailing explanations about the main features of `TensND`, it is worth recalling that the use of the libraries `TensND` and `SymPy` requires starting scripts by

```julia
julia> using TensND, SymPy
```

For numerical computations (without symbolic overhead), `SymPy` is not required:

```julia
julia> using TensND
```

The detailed manual is decomposed into the following chapters

```@contents
Pages = [
    "bases.md",
    "tensors.md",
    "coorsystems.md",
]
Depth = 1
```

## Tutorials

Several tutorials demonstrate advanced usage of TensND in various contexts:

```@contents
Pages = [
    "../tuto/nlayersphere.md",
    "../tuto/coorsystems_num.md",
    "../tuto/green_function.md",
    "../tuto/projection.md",
]
Depth = 1
```
