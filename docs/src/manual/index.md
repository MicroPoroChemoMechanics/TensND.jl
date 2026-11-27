# [Manual](@id man-index)

The Manual says how to call the package, one kind of object per page, on the smallest example that shows it. Each page states what the call computes and refers to the [Theory](@ref th-index) for the mathematics.

## First steps

[Installation](@ref man-installation) installs the package, states what to load for symbolic computations and how to check the installation. [Getting started](@ref man-getting-started) is a one-page tour of bases, tensors, structured tensors and symmetry detection, each step leading to the chapter that develops it.

## Tensors and bases

[Bases](@ref man-bases) constructs canonical, rotated, orthogonal and general bases and reads their dual basis and metric. [Tensors](@ref man-tensors) creates tensors of any order on a basis, converts their components between covariant and contravariant forms, and applies the products and contractions. [Two dimensions, and any dimension](@ref man-dimensions) says what is built in the plane and in dimension 4 and above, call by call, and what exists in dimension 3 only. [Random tensors](@ref man-random-tensors) draws random tensors of any symmetry class, with uniformly distributed orientations, in the plane and in space, for tests and sampling.

## Material symmetries

[Structured tensors](@ref man-structured) stores isotropic, transversely isotropic, cubic and orthotropic tensors, and the square and orthotropic tensors of the plane, with their few constants and computes their products and inverses. [Parametrizations and conversions](@ref man-parametrizations) moves between components, engineering constants and ratios for one symmetry class. [Elastic moduli](@ref man-elastic-moduli) builds an isotropic stiffness or compliance from any pair of moduli and reads any modulus from either tensor, in any dimension, with the two conventions of the plane. [Projection and symmetry detection](@ref man-projection) finds the closest tensor of a given symmetry class, with the orientation given or optimized, and detects the class of a tensor. [Symbolic tensors of a material symmetry](@ref man-symbolic-materials) builds a tensor of each class with symbols, in each notation, and shows it as its matrix: what a constructor takes, in which order, and whether it builds a stiffness or a compliance.

## Elastic waves

[Plane waves](@ref man-plane-waves) computes the velocities, polarizations and group velocities of the waves of a stiffness, exactly or to first order in the anisotropy, with Thomsen's parameters and their generalization. [Waves at an interface](@ref man-interface-waves) finds the waves that share a tangential slowness and the reflection and transmission coefficients of a welded or slipping interface, and the characteristic roots of a stiffness of the plane. [Third-order elasticity](@ref man-third-order) gives the energy, stress and tangent of a solid isotropic at orders two and three and the waves it carries under stress.

## Coordinate systems

[Charts, coordinate systems and surfaces](@ref man-coordinates-overview) says which of the three objects to take, exact or pointwise. [Coordinate systems (symbolic)](@ref man-coorsystems) defines a curvilinear chart and computes its basis vectors, Christoffel symbols and differential operators in closed form. [Coordinate systems (numerical)](@ref man-coorsystems-num) evaluates the same operators at a point by automatic differentiation. [Charts and coordinate conversions](@ref man-charts) passes from the coordinates of a point to its Cartesian components and back, and from one system to another. [Spheroidal harmonics](@ref man-spheroidal-harmonics) evaluates the harmonic functions of the spheroidal systems and re-expands them from one system to another.

## Surfaces

[Submanifolds](@ref man-submanifolds) describes a surface embedded in space and computes its fundamental forms, curvatures and differential operators, in closed form or at a point.

## Symbolic and numeric

[Symbolic and numeric](@ref man-symbolic-numeric) describes the four element types (`Float64`, `ForwardDiff.Dual`, `Sym`, `Num`), what is generic and what is not, differentiation, and interoperability with `Tensors.jl`.
