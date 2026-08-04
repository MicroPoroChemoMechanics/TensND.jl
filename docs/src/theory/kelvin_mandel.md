# [Kelvin–Mandel representation](@id th-kelvin-mandel)

Minor-symmetric order-4 tensors and symmetric order-2 tensors form spaces of
dimension 21 and 6 in 3-D. The Kelvin–Mandel representation makes that explicit
by mapping them onto ``6\times6`` matrices and ``6``-vectors **isometrically**,
which turns tensor algebra into matrix algebra. It goes back to Thomson
[thomson1856](@cite) and Mandel [mandel1965](@cite); its spectral consequences
are developed in [mehrabadi1990](@cite).

## The map

Index pairs are ordered

```math
11\to1,\quad 22\to2,\quad 33\to3,\quad 23\to4,\quad 13\to5,\quad 12\to6,
```

and off-diagonal entries carry a factor ``\sqrt2`` — one per off-diagonal pair,
so ``2`` where both index pairs are off-diagonal:

```math
\mathrm{Mat}(\boldsymbol{a})=
\begin{pmatrix} a_{11}\\ a_{22}\\ a_{33}\\
\sqrt2\,a_{23}\\ \sqrt2\,a_{13}\\ \sqrt2\,a_{12}\end{pmatrix},
\qquad
\mathrm{Mat}(\mathbb{T})_{KL}=
\begin{cases}
T_{ijkl} & K,L\le3\\
\sqrt2\,T_{ijkl} & \text{exactly one of }K,L>3\\
2\,T_{ijkl} & K,L>3 .
\end{cases}
```

These are [`KM`](@ref) and its inverse [`inv_KM`](@ref); the index ordering is
`_KM_COUPLES` in `src/tens_projection.jl`.

## What the ``\sqrt2`` buys

The factors are chosen so that the map is an **isometry** for the Frobenius
scalar product of [Tensor algebra](@ref th-tensor-algebra). Every line below is
an identity, verified in `test/`:

| Tensor statement | Matrix statement |
| :--------------- | :--------------- |
| ``\|\boldsymbol{a}\|^2=a_{ij}a_{ij}`` | ``{}^{t}\mathrm{Mat}(\boldsymbol{a})\,\mathrm{Mat}(\boldsymbol{a})`` |
| ``\mathbb{A}::\mathbb{B}`` | ``\mathrm{tr}\bigl({}^{t}\mathrm{Mat}(\mathbb{A})\,\mathrm{Mat}(\mathbb{B})\bigr)`` |
| ``\mathbb{A}:\mathbb{B}`` | ``\mathrm{Mat}(\mathbb{A})\,\mathrm{Mat}(\mathbb{B})`` |
| ``\mathbb{A}:\boldsymbol{a}`` | ``\mathrm{Mat}(\mathbb{A})\,\mathrm{Mat}(\boldsymbol{a})`` |
| ``\mathbb{A}^{-1}`` | ``\mathrm{Mat}(\mathbb{A})^{-1}`` |
| ``\mathbb{I}`` | ``\boldsymbol{1}_{6\times6}`` |

Double contraction is an ordinary matrix product and inversion an ordinary
matrix inverse. This is why every projection algorithm in this library is
written in the ``6\times6`` picture.

## Contrast with Voigt

Voigt notation uses *engineering* shear strains: a factor ``2`` on the
off-diagonal strain components and ``1`` on the stress ones. The consequences
matter:

| | Kelvin–Mandel | Voigt |
| :--- | :--- | :--- |
| stress and strain treated | identically | differently |
| ``\|\cdot\|_F`` preserved | yes | no |
| ``\mathbb{A}:\mathbb{B}\to`` matrix product | yes | no |
| ``\mathbb{A}^{-1}\to`` matrix inverse | yes | no |
| rotation matrix ``Q`` | **orthogonal** | not orthogonal |

The last line is the decisive one for this library and is developed next.

## Rotation: an orthogonal congruence

Let ``R\in SO(3)`` and let ``\mathbb{R}=R\stackrel{s}{\boxtimes}R`` be the
order-4 tensor rotating minor-symmetric tensors
([Rotations](@ref th-rotations)). Its Kelvin–Mandel matrix

```math
Q=\mathrm{Mat}(R\stackrel{s}{\boxtimes}R)
\qquad\text{satisfies}\qquad
{}^{t}Q\,Q=\boldsymbol{1}_{6\times6},
```

so a change of frame is an **orthogonal congruence**

```math
\mathrm{Mat}(\mathbb{A})\;\longmapsto\;{}^{t}Q\,\mathrm{Mat}(\mathbb{A})\,Q .
```

`TensND` builds it either from Euler angles (`_KM_rotation(θ, ϕ, ψ)`) or
directly from a frame (`_km_congruence`, in `src/tens_walpole.jl`); the latter
is what relates the two views of an orthotropic tensor,

```math
\mathrm{Mat}(\mathbb{C})=Q\;\mathrm{Mat}_{\text{mat}}(\mathbb{C})\;{}^{t}Q ,
```

with ``\mathrm{Mat}_{\text{mat}}`` the block-diagonal matrix in the material
frame ([Orthotropy](@ref th-orthotropy)).

!!! note "Why orthogonality is the whole point"
    Because ``Q`` is orthogonal, the congruence is a *rotation of the
    21-dimensional Euclidean space of elasticity tensors*: it preserves
    ``\|\cdot\|_F``, hence preserves distances to a symmetry class. Every
    projection in [Projection onto a symmetry class](@ref th-projection) may
    therefore be computed in whichever frame is convenient and rotated back,
    and the optimization over orientation is an optimization over a compact
    group acting by isometries. In Voigt notation none of this holds.

    This identity is easy to get wrong and easy to *appear* right: in the
    canonical frame ``Q=\boldsymbol{1}`` and every convention error hides. It is
    therefore pinned by a test on a **rotated** frame
    ([Testing and conventions](@ref dev-testing)).
