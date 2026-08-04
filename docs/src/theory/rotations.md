# [Rotations and Euler angles](@id th-rotations)

Every symmetry class in this library is defined up to an orientation — an axis
for transverse isotropy, a full frame for orthotropy — so the angular convention
is part of the specification, not a detail. This page states the one `TensND`
implements, read off `_rot3_raw` in `src/tens_projection.jl` and `rot2`/`rot3`/
`rot6` in `src/special_tens.jl`.

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
``\varphi`` no longer affects ``\underline{n}``, is why the multi-start grids in
`ext/TensNDNLoptExt.jl` filter those duplicates.

!!! note "Spherical coordinates use the same angles in a different slot order"
    The predefined spherical coordinate system returns its coordinates as
    ``(\theta,\varphi,r)``, not ``(r,\theta,\varphi)``, precisely so that
    ``\theta=\varphi=0`` reproduces the canonical basis in the canonical order.
    See [Curvilinear differential calculus](@ref th-curvilinear).

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
