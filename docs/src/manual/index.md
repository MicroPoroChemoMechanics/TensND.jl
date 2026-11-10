# [Manual](@id man-index)

The Manual says how to call the package, one kind of object per page, on the smallest example that shows it. Each page states what the call computes and refers to the [Theory](@ref th-index) for the mathematics.

## First steps

[Installation](@ref man-installation) installs the package, states what to load for symbolic computations and how to check the installation. [Getting started](@ref man-getting-started) is a one-page tour of bases, tensors, structured tensors and symmetry detection, each step leading to the chapter that develops it.

## Tensors and bases

[Bases](@ref man-bases) constructs canonical, rotated, orthogonal and general bases and reads their dual basis and metric. [Tensors](@ref man-tensors) creates tensors of any order on a basis, converts their components between covariant and contravariant forms, and applies the products and contractions. [Random tensors](@ref man-random-tensors) draws random tensors of any symmetry class, with uniformly distributed orientations, for tests and sampling.

## Material symmetries

[Structured tensors](@ref man-structured) stores isotropic, transversely isotropic, cubic and orthotropic tensors with their few constants and computes their products and inverses. [Parametrizations and conversions](@ref man-parametrizations) moves between components, engineering constants and ratios for one symmetry class. [Elastic moduli](@ref man-elastic-moduli) builds an isotropic stiffness or compliance from any pair of moduli and reads any modulus from either tensor, in three and two dimensions. [Projection and symmetry detection](@ref man-projection) finds the closest tensor of a given symmetry class, with the orientation given or optimized, and detects the class of a tensor.

## Coordinate systems

[Coordinate systems (symbolic)](@ref man-coorsystems) defines a curvilinear chart and computes its basis vectors, Christoffel symbols and differential operators in closed form. [Coordinate systems (numerical)](@ref man-coorsystems-num) evaluates the same operators at a point by automatic differentiation. [Charts and coordinate conversions](@ref man-charts) passes from the coordinates of a point to its Cartesian components and back, and from one system to another. [Spheroidal harmonics](@ref man-spheroidal-harmonics) evaluates the harmonic functions of the spheroidal systems and re-expands them from one system to another.

## Surfaces

[Submanifolds](@ref man-submanifolds) describes a surface embedded in space and computes its fundamental forms and curvatures.

## Symbolic and numeric

[Symbolic and numeric](@ref man-symbolic-numeric) describes the four element types (`Float64`, `ForwardDiff.Dual`, `Sym`, `Num`), what is generic and what is not, differentiation, and interoperability with `Tensors.jl`.
