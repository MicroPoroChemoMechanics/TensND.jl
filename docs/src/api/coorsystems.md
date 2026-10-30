# [Coordinate systems (symbolic)](@id api-coorsystems)

Symbolic differential calculus on a curvilinear chart. Theory:
[Curvilinear differential calculus](@ref th-curvilinear); usage:
[Coordinate systems](@ref man-coorsystems).

```@docs
CoordinateSystem
AbstractGeometry
ExactGeometry
@set_coorsys
set_coorsys!
default_coorsys
unset_coorsys!
grad
symgrad
divergence
laplacian
hess
∂
get_coords
position_vector
unitvec
natvec
normalized_basis
natural_basis
lame_coefficients
christoffel
nderiv
coorsys_cartesian
coorsys_polar
coorsys_cylindrical
coorsys_spherical
coorsys_spheroidal
```
