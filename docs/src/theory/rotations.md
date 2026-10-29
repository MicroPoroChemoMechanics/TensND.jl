# [Rotations and Euler angles](@id th-rotations)

Every symmetry class in this library is defined up to an orientation — an axis
for transverse isotropy, a full frame for orthotropy — so the angular convention
is part of the specification, not a detail. This page states the one `TensND`
implements, read off `_rot3_raw` in `src/tens_projection.jl` and `rot2`/`rot3`/
`rot6` in `src/special_tens.jl`.

```@raw html
<ClientOnly><FrameFigure kind="euler" /></ClientOnly>
```

The rotation ``R(\theta,\varphi,\psi)`` of this page: the canonical basis
``\underline{e}_i`` at the center of the unit sphere, the rotated basis
``\underline{e}'_i=R\cdot\underline{e}_i`` at the tip of
``\underline{n}=\underline{e}'_3``. The sliders set the three angles and the
scene turns with the mouse; how the three rotations compose is read in
[The three rotations](@ref th-rotations-composed).

## The convention: **Z–Y–Z**, in the order ``(\theta,\varphi,\psi)``

```math
R(\theta,\varphi,\psi)=R_z(\varphi)\,R_y(\theta)\,R_z(\psi)
```

Note the argument order: the function takes ``(\theta,\varphi,\psi)`` but the
*first* rotation applied is ``R_z(\psi)``. In `Rotations.jl` terms this is
exactly `RotZYZ(ϕ, θ, ψ)`.

Explicitly, with ``c_\bullet=\cos\bullet`` and ``s_\bullet=\sin\bullet``:

```math
R=
\begin{pmatrix}
c_\theta c_\psi c_\varphi-s_\psi s_\varphi &
-c_\theta c_\varphi s_\psi-c_\psi s_\varphi &
c_\varphi s_\theta\\[2pt]
c_\theta c_\psi s_\varphi+c_\varphi s_\psi &
-c_\theta s_\psi s_\varphi+c_\psi c_\varphi &
s_\theta s_\varphi\\[2pt]
-c_\psi s_\theta & s_\theta s_\psi & c_\theta
\end{pmatrix}
```

The columns of ``R`` are the new basis vectors expressed in the old frame, so
that [`Basis`](@ref)`(θ, ϕ, ψ)` builds the [`RotatedBasis`](@ref) they span.

### The third column is the axis

```math
\underline{n}=R\cdot\underline{e}_3
=\bigl(\sin\theta\cos\varphi,\;\sin\theta\sin\varphi,\;\cos\theta\bigr)
```

``\theta`` is therefore the **polar** angle from ``\underline{e}_3`` and
``\varphi`` the **azimuth**. This is what makes the convention the natural one
here: a transversely isotropic tensor is oriented by ``\underline{n}`` alone, so
its projection is parametrized by ``(\theta,\varphi)`` with ``\psi`` irrelevant,
while orthotropy needs all three. The degeneracy at ``\theta=0``, where
``\varphi`` no longer affects ``\underline{n}``, is why the orientation search
screens axes spread over a hemisphere rather than a grid in
``(\theta,\varphi)``, which would crowd them at the pole.

!!! note "Spherical coordinates use the same angles in a different slot order"
    The predefined spherical coordinate system returns its coordinates as
    ``(\theta,\varphi,r)``, not ``(r,\theta,\varphi)``, precisely so that
    ``\theta=\varphi=0`` reproduces the canonical basis in the canonical order.
    See [Curvilinear differential calculus](@ref th-curvilinear).

### [The three rotations](@id th-rotations-composed)

Read from the left, ``R=R_z(\varphi)\,R_y(\theta)\,R_z(\psi)`` is a chain of
rotations, each about an axis that the previous ones have moved:

1. ``\varphi`` about ``\underline{e}_3`` turns ``\underline{e}_2`` into the
   *line of nodes* ``\underline{e}^\varphi=-\sin\varphi\,\underline{e}_1+\cos\varphi\,\underline{e}_2``;
2. ``\theta`` about the line of nodes tilts ``\underline{e}_3`` into
   ``\underline{n}``, within the meridian plane of azimuth ``\varphi``, and
   carries the image of ``\underline{e}_1`` to
   ``\underline{e}^\theta=\cos\theta\cos\varphi\,\underline{e}_1+\cos\theta\sin\varphi\,\underline{e}_2-\sin\theta\,\underline{e}_3``;
3. ``\psi`` about ``\underline{n}`` turns the pair
   ``(\underline{e}^\theta,\underline{e}^\varphi)`` within the plane tangent to
   the unit sphere at ``\underline{n}``.

The rotated basis ``\underline{e}'_i=R\cdot\underline{e}_i`` is therefore

```math
\begin{aligned}
\underline{e}'_1&=\cos\psi\,\underline{e}^\theta+\sin\psi\,\underline{e}^\varphi,\\
\underline{e}'_2&=-\sin\psi\,\underline{e}^\theta+\cos\psi\,\underline{e}^\varphi,\\
\underline{e}'_3&=\underline{n} ,
\end{aligned}
```

and ``(\underline{e}^\theta,\underline{e}^\varphi,\underline{n})``, the frame
``R(\theta,\varphi,0)`` makes of the canonical basis, is the normalized basis of
the spherical coordinates at the point ``\underline{n}``
([Predefined systems](@ref th-predefined-systems)); the figure at the top of
this page draws ``\underline{e}^\theta`` and ``\underline{e}^\varphi`` dashed,
in the plane tangent to the sphere. At
``\theta=0`` the line of nodes is still drawn, but ``\varphi`` and ``\psi``
then turn about the same axis and only their sum matters: the gimbal lock met
again under *Recovering angles from a frame*.

## Rotating tensors of each order

A rotation acts on an order-``p`` tensor by rotating every index. `TensND`
provides the three cases that occur in practice:

| Function | Returns | Acts on |
| :------- | :------ | :------ |
| [`rot2`](@ref)`(θ)` | ``2\times2`` rotation | 2-D |
| [`rot3`](@ref)`(θ, ϕ, ψ)` | ``3\times3`` rotation ``R`` | vectors, order-2 tensors |
| [`rot6`](@ref)`(θ, ϕ, ψ)` | order-4 tensor | minor-symmetric order-4 tensors |

with

```math
\mathbb{R}=R\stackrel{s}{\boxtimes}R ,
\qquad
\mathbb{A}'=\mathbb{R}:\mathbb{A}:{}^{t}\mathbb{R} ,
```

the symmetrized box product being what preserves minor symmetry. In the
Kelvin–Mandel picture this order-4 object is the ``6\times6`` matrix

```math
Q=\mathrm{Mat}(\mathbb{R})=\mathrm{Mat}\bigl(R\stackrel{s}{\boxtimes}R\bigr),
\qquad {}^{t}Q\,Q=\boldsymbol{1},
```

sometimes called the Bond matrix. Its orthogonality — the property Voigt
notation lacks — is what
[Kelvin–Mandel representation](@ref th-kelvin-mandel) is about.

## [Uniform random rotations](@id th-random-rotations)

A rotation drawn "at random" should favor no orientation: its law must be
unchanged when every rotation is composed with a fixed one. On ``SO(3)`` there
is exactly one such probability, the Haar measure. In the angles of this page
it reads

```math
\mathrm{d}\mu = \frac{1}{8\pi^2}\,\sin\theta\,\mathrm{d}\theta\,\mathrm{d}\varphi\,\mathrm{d}\psi ,
\qquad \theta\in[0,\pi],\ \varphi,\psi\in[0,2\pi),
```

so drawing ``\theta`` uniformly is wrong: the axis
``\underline{n}=R\cdot\underline{e}_3``, of height ``\cos\theta``, would fall near
the poles as often as near the equator, where the area to fill is largest.
The height must be uniform on ``[-1,1]`` instead (the sphere and its
circumscribed cylinder have the same area between two horizontal planes), so
that ``\mathbb{E}[n_3^2]=1/3``, where uniform angles give ``1/2``.

[`random_rotation`](@ref) does without angles. Unit quaternions ``q`` form the
sphere ``S^3`` of ``\mathbb{R}^4``, and ``q\mapsto R_q`` maps them two to one
onto ``SO(3)``, products to products. Multiplying by a fixed unit quaternion,
on either side, is an orthogonal map of ``\mathbb{R}^4``, which preserves the
uniform law on ``S^3``; its image on ``SO(3)`` is therefore invariant under the
compositions on either side, hence the Haar measure. The uniform law on ``S^3``
is the direction of a vector of four independent standard normal numbers,
whose law is invariant under every orthogonal map. Likewise
[`random_axis`](@ref) normalizes three standard normal numbers.

Two consequences serve as checks. The angle ``\omega`` of a Haar rotation,
``\operatorname{tr}R = 1 + 2\cos\omega``, has the density
``(1-\cos\omega)/\pi`` on ``[0,\pi]``: small rotations are rare. And the trace
has mean ``0`` and mean square ``1``, the moments of the character of an
irreducible representation.

**Coefficients.** A rotation acts on the Kelvin–Mandel vector of a symmetric
tensor of order 2 by an orthogonal ``6\times 6`` matrix
([Kelvin–Mandel representation](@ref th-kelvin-mandel)), so a law of
Kelvin–Mandel matrices invariant under all orthogonal changes of basis is
invariant under rotations. The Gaussian orthogonal ensemble (symmetric, normal
entries of variance ``1`` on the diagonal and ``1/2`` off it) is one; so is
``Q\,\Lambda\,Q^{\mathsf{T}}`` with ``Q`` of Haar law on ``O(6)`` and
``\Lambda`` any law of eigenvalues, and a matrix of independent standard normal
entries. [`random_tens`](@ref) draws generic tensors from them, and the blocks
of the structured classes from their analogues of size 2 or 3. In Voigt
notation the same draws would not be invariant: the change of basis is not
orthogonal there.

## Recovering angles from a frame

The inverse problem, extracting ``(\theta,\varphi,\psi)`` from an orthonormal
matrix, is `angles(M, Val{3})` in `src/bases.jl`. It is used when a symmetry
frame has been obtained as an eigenvector basis and must be reported as angles.

Two caveats, inherent to any three-angle parametrization of ``SO(3)``:

- **the map is not injective** — ``(\theta,\varphi,\psi)`` and
  ``(-\theta,\varphi+\pi,\psi+\pi)`` give the same ``R``;
- **it degenerates at the poles** (gimbal lock): at ``\theta=0`` the matrix
  reduces to ``R_z(\varphi+\psi)``, so only the *sum* is determined; at
  ``\theta=\pi`` only the *difference* ``\varphi-\psi`` is.

Neither affects the projections: the objective functions of
[Projection onto a symmetry class](@ref th-projection) depend on the frame, not
on the angles chosen to name it, so a degenerate parametrization costs at worst
a redundant starting point.
