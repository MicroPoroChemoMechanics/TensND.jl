# [Classical coordinate systems](@id th-predefined-systems)

!!! tip "See also"
    **Theory**: [Curvilinear differential calculus](@ref th-curvilinear)

    **Manual**: [Coordinate systems (symbolic)](@ref man-coorsystems) · [Coordinate systems (numerical)](@ref man-coorsystems-num)

    **Tutorials**: [Differential operators, symbolically](@ref tut-symbolic-operators) · [Christoffel symbols](@ref tut-christoffel) · [Differential operators, numerically](@ref tut-numerical-operators) · [Classical identities of the differential operators](@ref tut-operator-identities) · [The operators in the classical systems](@ref tut-operators-catalog)

The Cartesian, polar, cylindrical and spherical systems, each with its chart,
its unit vectors and the operators that
[Curvilinear differential calculus](@ref th-orthogonal-operators) gives for an
orthogonal system.

| System | Coordinates | ``\underline{OM}`` | Lamé ``\chi_i`` |
| :----- | :---------- | :----------------- | :-------------- |
| [`coorsys_cartesian`](@ref) | ``(x,y,z)`` | ``(x,\,y,\,z)`` | ``(1,1,1)`` |
| [`coorsys_polar`](@ref) | ``(r,\theta)`` | ``(r\cos\theta,\;r\sin\theta)`` | ``(1,\,r)`` |
| [`coorsys_cylindrical`](@ref) | ``(r,\theta,z)`` | ``(r\cos\theta,\;r\sin\theta,\;z)`` | ``(1,\,r,\,1)`` |
| [`coorsys_spherical`](@ref) | ``(\theta,\varphi,r)`` | ``(r\sin\theta\cos\varphi,\;r\sin\theta\sin\varphi,\;r\cos\theta)`` | ``(r,\;r\sin\theta,\;1)`` |

All four are orthogonal, so each normalized basis is orthonormal; its vectors
are named after their coordinate, ``\underline{e}^r``, ``\underline{e}^\theta``,
… — for an orthonormal basis the dual basis is the basis itself, and the
position of the label carries no meaning.

!!! warning "Ordering of the spherical coordinates"
    The spherical coordinates are ordered ``(\theta,\varphi,r)``, not
    ``(r,\theta,\varphi)``. The ordering is chosen so that
    ``\theta=\varphi=0`` reproduces the canonical basis **in the canonical
    order**, which makes the spherical frame a genuine
    [`RotatedBasis`](@ref) — see [Rotations](@ref th-rotations) — and lets
    spherical results be compared with Euler-angle results without a
    permutation. The same ordering applies to
    [`init_spherical`](@ref) and to ``\underline{e}^\theta,\underline{e}^\varphi,
    \underline{e}^r``.

The figures below draw each curvilinear system at a point ``M`` set by the
sliders: the canonical basis at the origin ``O``, the normalized basis at
``M``, the angles and lengths that locate ``M``, and, when the box is ticked,
the lines (in 2-D) or surfaces (in 3-D) on which one coordinate is constant.
Each coordinate keeps one color throughout. Since
``\underline{a}^i=\underline{a}_i/\chi_i^2`` is the gradient of ``q^i``, the
unit vector of a coordinate is normal to the surface on which that coordinate
is constant, and both are drawn in its color. A 3-D figure turns with the
mouse.

## Polar coordinates

With ``r>0`` and ``\theta`` taken on any interval of length ``2\pi`` (the
slider of the figure runs from ``-\pi`` to ``\pi``):

```math
\underline{e}^r=\cos\theta\,\underline{e}_1+\sin\theta\,\underline{e}_2,
\qquad
\underline{e}^\theta=-\sin\theta\,\underline{e}_1+\cos\theta\,\underline{e}_2 .
```

```@raw html
<ClientOnly><FrameFigure kind="polar" /></ClientOnly>
```

The Lamé coefficients ``(1,r)`` give, for a scalar field ``f(r,\theta)`` and a
vector field ``\underline{v}=v_r\,\underline{e}^r+v_\theta\,\underline{e}^\theta``,

```math
\underline{\nabla}f=\frac{\partial f}{\partial r}\,\underline{e}^r
+\frac{1}{r}\,\frac{\partial f}{\partial\theta}\,\underline{e}^\theta,
\qquad
\mathrm{div}\,\underline{v}=\frac{1}{r}\,\frac{\partial(r\,v_r)}{\partial r}
+\frac{1}{r}\,\frac{\partial v_\theta}{\partial\theta},
\qquad
\Delta f=\frac{1}{r}\,\frac{\partial}{\partial r}\Bigl(r\,\frac{\partial f}{\partial r}\Bigr)
+\frac{1}{r^2}\,\frac{\partial^2f}{\partial\theta^2},
```

and the non-vanishing Christoffel symbols are

```math
\Gamma^r_{\theta\theta}=-r,
\qquad
\Gamma^\theta_{r\theta}=\Gamma^\theta_{\theta r}=\frac{1}{r} .
```

## Cylindrical coordinates

The polar pair in the plane ``z=0``, completed by the axial direction:
``\underline{e}^r`` and ``\underline{e}^\theta`` as above, and
``\underline{e}^z=\underline{e}_3``. The surfaces of constant coordinate are the
cylinders of axis ``\underline{e}_3`` (``r``), the half-planes bounded by that
axis (``\theta``) and the planes normal to it (``z``).

```@raw html
<ClientOnly><FrameFigure kind="cylindrical" /></ClientOnly>
```

The operators are the polar ones completed by the axial terms,

```math
\mathrm{div}\,\underline{v}=\frac{1}{r}\,\frac{\partial(r\,v_r)}{\partial r}
+\frac{1}{r}\,\frac{\partial v_\theta}{\partial\theta}+\frac{\partial v_z}{\partial z},
\qquad
\Delta f=\frac{1}{r}\,\frac{\partial}{\partial r}\Bigl(r\,\frac{\partial f}{\partial r}\Bigr)
+\frac{1}{r^2}\,\frac{\partial^2f}{\partial\theta^2}+\frac{\partial^2f}{\partial z^2},
```

and the Christoffel symbols are the polar ones.

## Spherical coordinates

With ``r>0``, the colatitude ``\theta\in[0,\pi]`` measured from
``\underline{e}_3`` and the azimuth ``\varphi`` measured from ``\underline{e}_1``:

```math
\begin{aligned}
\underline{e}^\theta&=\cos\theta\cos\varphi\,\underline{e}_1
+\cos\theta\sin\varphi\,\underline{e}_2-\sin\theta\,\underline{e}_3,\\
\underline{e}^\varphi&=-\sin\varphi\,\underline{e}_1+\cos\varphi\,\underline{e}_2,\\
\underline{e}^r&=\sin\theta\cos\varphi\,\underline{e}_1
+\sin\theta\sin\varphi\,\underline{e}_2+\cos\theta\,\underline{e}_3 .
\end{aligned}
```

In this order the triad is the image of the canonical basis by the rotation
``R(\theta,\varphi,0)`` of [Rotations and Euler angles](@ref th-rotations). The
surfaces of constant coordinate are the cones of axis ``\underline{e}_3`` and
apex ``O`` (``\theta``), the half-planes bounded by that axis (``\varphi``) and
the spheres of center ``O`` (``r``).

```@raw html
<ClientOnly><FrameFigure kind="spherical" /></ClientOnly>
```

With the Lamé coefficients ``(r,\,r\sin\theta,\,1)`` and ``J=r^2\sin\theta``,

```math
\underline{\nabla}f=\frac{1}{r}\,\frac{\partial f}{\partial\theta}\,\underline{e}^\theta
+\frac{1}{r\sin\theta}\,\frac{\partial f}{\partial\varphi}\,\underline{e}^\varphi
+\frac{\partial f}{\partial r}\,\underline{e}^r ,
```

```math
\mathrm{div}\,\underline{v}=\frac{1}{r\sin\theta}\,\frac{\partial(\sin\theta\,v_\theta)}{\partial\theta}
+\frac{1}{r\sin\theta}\,\frac{\partial v_\varphi}{\partial\varphi}
+\frac{1}{r^2}\,\frac{\partial(r^2\,v_r)}{\partial r},
```

```math
\Delta f=\frac{1}{r^2\sin\theta}\,\frac{\partial}{\partial\theta}\Bigl(\sin\theta\,\frac{\partial f}{\partial\theta}\Bigr)
+\frac{1}{r^2\sin^2\theta}\,\frac{\partial^2f}{\partial\varphi^2}
+\frac{1}{r^2}\,\frac{\partial}{\partial r}\Bigl(r^2\,\frac{\partial f}{\partial r}\Bigr),
```

and the non-vanishing Christoffel symbols are

```math
\Gamma^r_{\theta\theta}=-r,\quad
\Gamma^r_{\varphi\varphi}=-r\sin^2\theta,\quad
\Gamma^\theta_{\varphi\varphi}=-\sin\theta\cos\theta,
\quad
\Gamma^\theta_{r\theta}=\Gamma^\theta_{\theta r}=\Gamma^\varphi_{r\varphi}=\Gamma^\varphi_{\varphi r}=\frac{1}{r},\quad
\Gamma^\varphi_{\theta\varphi}=\Gamma^\varphi_{\varphi\theta}=\cot\theta .
```

The same construction applied to confocal quadrics instead of spheres and
cones gives the systems of [Spheroidal coordinates](@ref th-spheroidal) and
[Ellipsoidal coordinates](@ref th-ellipsoidal).
