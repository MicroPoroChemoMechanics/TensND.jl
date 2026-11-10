# [Changes of coordinates](@id th-coordinate-changes)

!!! tip "See also"
    **Manual**: [Charts and coordinate conversions](@ref man-charts)

    **Tutorials**: [Ellipsoidal coordinates](@ref tut-ellipsoidal) · [Changing the focal distance of spheroidal harmonics](@ref tut-spheroidal-reexpansion)

A coordinate system is used in two directions: from the coordinates to the
point, which is the position vector, and from the point to its coordinates,
which a field given in another system requires. This page writes the second
direction in closed form for each predefined system, then composes the two to
pass from one system to another. A **chart** denotes here the pair of maps,
between numbers, that a system defines.

## From coordinates to a point and back

For a chart of coordinates ``\xi=(q^1,\ldots,q^d)``,
[`coords_to_cartesian`](@ref) returns the Cartesian components
``\underline{x}(\xi)`` of the position vector and
[`cartesian_to_coords`](@ref) the inverse map ``\xi(\underline{x})``. Their
Jacobian matrices are inverse of each other, and that of the first has the
natural vectors as columns: in an orthogonal system,

```math
\frac{\partial\underline{x}}{\partial\xi}=R\,\mathrm{diag}(\chi_i),
\qquad
\frac{\partial\xi}{\partial\underline{x}}=\mathrm{diag}\Bigl(\frac{1}{\chi_i}\Bigr)\,R^{\mathsf T},
```

where the columns of ``R`` are the unit vectors ``\underline{e}_i``. The row
``i`` of the second matrix is the gradient of the coordinate ``q^i``,
``\underline{a}^i=\underline{e}_i/\chi_i``. Both maps are differentiable by
automatic differentiation, with respect to the point and to the parameters of
the chart.

| Chart | Coordinates | Parameters | Singular set |
| :---- | :---------- | :--------- | :----------- |
| [`PolarChart`](@ref) | ``(r,\theta)`` | — | the origin |
| [`CylindricalChart`](@ref) | ``(r,\theta,z)`` | — | the axis |
| [`SphericalChart`](@ref) | ``(\theta,\varphi,r)`` | — | the axis |
| [`ProlateSpheroidalChart`](@ref) | ``(\varphi,p,q)`` | ``c`` | the axis |
| [`OblateSpheroidalChart`](@ref) | ``(\varphi,p,\tau)`` | ``c`` | the axis and the focal disk |
| [`EllipsoidalChart`](@ref) | ``(\rho,\mu,\nu)`` | ``h_2,h_3`` | the coordinate planes |

On its singular set a chart still returns coordinates, but its Jacobian is not
invertible there.

## Classical systems

With ``\varrho=\sqrt{x_1^2+x_2^2}``, the polar and cylindrical coordinates of a
point are ``r=\varrho``, ``\theta=\mathrm{atan2}(x_2,x_1)`` and ``z=x_3``; its
spherical coordinates, in the order of
[Classical coordinate systems](@ref th-predefined-systems), are

```math
\theta=\mathrm{atan2}(\varrho,x_3),
\qquad
\varphi=\mathrm{atan2}(x_2,x_1),
\qquad
r=\sqrt{\varrho^2+x_3^2} .
```

The two-argument arctangent keeps its accuracy near the axis, where an
arccosine of ``x_3/r`` would lose half of the digits.

## Prolate spheroidal system

By the construction of [Spheroidal coordinates](@ref th-spheroidal), ``2\,c\,q``
is the sum of the distances from the point to the two foci:

```math
q=\frac{r_++r_-}{2c},
\qquad
r_\pm=\sqrt{\varrho^2+(x_3\mp c)^2},
\qquad
p=\frac{x_3}{c\,q},
\qquad
\varphi=\mathrm{atan2}(x_2,x_1) .
```

The coordinate ``p`` is also the half-difference ``(r_--r_+)/(2c)``; it is
taken from ``x_3=c\,p\,q`` instead, since the difference of two close distances
loses accuracy far from the foci.

## Oblate spheroidal system

Eliminating ``p`` between ``\varrho^2=c^2\,(1-p^2)\,(\tau^2+1)`` and
``x_3^2=c^2\,p^2\,\tau^2`` shows that ``\tau^2`` and ``-p^2`` are the two roots
of

```math
c^2\,u^2-d\,u-x_3^2=0,
\qquad
d=\varrho^2+x_3^2-c^2 .
```

Their product is ``-x_3^2/c^2``. The root of larger modulus adds two terms of
the same sign,

```math
u_>=\frac{|d|+\sqrt{d^2+4\,c^2\,x_3^2}}{2\,c^2},
\qquad
u_<=\frac{x_3^2}{c^2\,u_>},
```

and the other follows from the product, so that no difference of close numbers
is formed. Outside the sphere of radius ``c`` (``d\ge0``), ``\tau^2=u_>`` and
``p^2=u_<``; inside, ``p^2=u_>`` and ``\tau^2=u_<``. The sign of ``p`` is that
of ``x_3``. In the focal disk ``\tau=0`` and ``p=\sqrt{1-\varrho^2/c^2}``: the
two faces of the disk, ``p`` and ``-p``, are the same points.

## Ellipsoidal system

The three coordinates are the roots of a cubic, solved in closed form in
[Ellipsoidal coordinates](@ref th-ellipsoidal). The chart covers the quarter
space ``x_2\ge0``, ``x_3\ge0``, and the coordinates returned for a point of
another quarter are those of its mirror image.

## From one system to another

The coordinates of a point in a chart ``\mathcal{C}_2``, given its coordinates
in a chart ``\mathcal{C}_1``, are the composition
``\xi_2=\xi_{\mathcal{C}_2}\bigl(\underline{x}_{\mathcal{C}_1}(\xi_1)\bigr)``,
which [`change_coords`](@ref) evaluates. Between two prolate systems of the
same axis and of semi-focal distances ``c_1`` and ``c_2``, the meridian
coordinates ``(p_2,q_2)`` only depend on ``(p_1,q_1)`` and on the ratio
``c_1/c_2``, and the azimuth is unchanged: with
``\varrho=c_1\sqrt{1-p_1^2}\sqrt{q_1^2-1}`` and ``x_3=c_1\,p_1\,q_1``,

```math
q_2=\frac{\sqrt{\varrho^2+(x_3-c_2)^2}+\sqrt{\varrho^2+(x_3+c_2)^2}}{2\,c_2},
\qquad
p_2=\frac{c_1\,p_1\,q_1}{c_2\,q_2} .
```

A spheroid ``q_1=\text{const}`` of the first system is not a coordinate surface
of the second: the two families are confocal with different focal sets. A field
written as a series of harmonics of one system is therefore re-expanded in the
harmonics of the other, which is the subject of
[Spheroidal harmonics and their re-expansion](@ref th-spheroidal-harmonics).
