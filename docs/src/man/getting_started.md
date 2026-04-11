# Getting started

## Brief description of the package

The package relies on the definition of

- **bases** which can be of the following types (`T` denotes the scalar type, subtype of `Number`)

  - `CanonicalBasis{dim,T}`: fundamental canonical basis in `ℝᵈⁱᵐ` in which the metric tensor is the second-order identity
  - `RotatedBasis{dim,T}`: orthonormal basis in `ℝᵈⁱᵐ` obtained by rotation of the canonical basis by means of one angle if `dim=2` or three Euler angles if `dim=3`, the metric tensor is again the second-order identity
  - `OrthogonalBasis{dim,T}`: orthogonal basis in `ℝᵈⁱᵐ` obtained from a given orthonormal rotated basis by applying a scaling factor along each unit vector, the metric tensor is then diagonal
  - `Basis{dim,T}`: arbitrary basis not entering the previous cases

- **tensors**

  - a tensor is determined by a set of data (array or synthetic parameters) corresponding to its `order`, a basis and a tuple of variances
  - depending on the type of basis, the type of tensor can be `TensCanonical{order,dim,T,A}`, `TensRotated{order,dim,T,A}`, `TensOrthogonal{order,dim,T,A}` or `Tens{order,dim,T,A}` if the data are stored under the form of an array or a `Tensor` object (see [Tensors.jl](https://github.com/Ferrite-FEM/Tensors.jl)), or the type can be `TensISO{order,dim,T,N}` if the tensor is isotropic and data are stored under the form of a set of parameters (one for order 2 and two for order 4).
  - `TensWalpole{T,N}` if the tensor is transversely isotropic of order 4 (Walpole basis, N=5 for major-symmetric or N=6 for general case).
  - `TensTI{order,T,N}` if the tensor is transversely isotropic of order 2 (N=2, data = (a, b) with symmetry axis n).
  - `TensOrtho{T}` if the tensor is orthotropic (9 elastic constants in a material frame).
  - Projection functions `proj_tens` and `best_sym_tens` allow finding the closest tensor of a given symmetry class (ISO, TI, ORTHO) with or without rotation optimization (the latter requires `NLopt`).

- **coordinate systems**

  - a coordinate system contains all information required to perform differential operations on tensor fields: position vector `OM` expressed in the canonical basis, coordinate names `xⁱ`, natural basis `aᵢ=∂ᵢOM`, normalized basis `eᵢ=aᵢ/||aᵢ||`, Christoffel coefficients `Γᵢⱼᵏ=∂ᵢaⱼ⋅aᵏ` where `(aⁱ)(1≤i≤dim)` form the dual basis associated to the natural one
  - predefined coordinate systems are available: cartesian, polar, cylindrical, spherical and spheroidal but the user can define new systems

  **_NOTE:_**
      Note that for the moment the coordinate systems and differential operators have been implemented only for symbolic calculations (using [SymPy.jl](https://github.com/JuliaPy/SymPy.jl)). Numerical coordinate systems and differential operators based on automatic differentiation will be implemented in the future.

## Detailed manual

Before detailing explanations about the main features of `TensND`, it is worth recalling that the use of the libraries `TensND` and `SymPy` requires starting scripts by

```julia
julia> using TensND, SymPy
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
