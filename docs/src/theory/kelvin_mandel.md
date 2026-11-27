# [Kelvin–Mandel representation](@id th-kelvin-mandel)

!!! tip "See also"
    **Manual**: [Tensors](@ref man-tensors)

    **Tutorials**: [Tensor products, contractions and their identities](@ref tut-tensor-algebra) · [Performance of the structured types](@ref tut-performance)

Minor-symmetric order-4 tensors and symmetric order-2 tensors form spaces of
dimension 21 and 6 in 3-D. The Kelvin–Mandel representation makes that explicit
by mapping them onto ``6\times6`` matrices and ``6``-vectors **isometrically**,
which turns tensor algebra into matrix algebra. It goes back to
[thomson1856](@citet) and [mandel1965](@citet); its spectral consequences are
developed by [mehrabadi1990](@citet).

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

These are [`mandel`](@ref) and its inverse [`from_mandel`](@ref); the index ordering is
`_KM_COUPLES` in `src/tens_projection.jl`.

## Role of the √2 factors

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

## [Contrast with Voigt](@id th-voigt)

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

The Voigt matrices remain the language of the experimental and seismic
literature, where a stiffness is given by its *two-index components*
``C_{IJ}``, in the same order of the pairs (11, 22, 33, 23, 13, 12). With the
weights ``w = (1, 1, 1, \sqrt2, \sqrt2, \sqrt2)`` and ``M`` the Kelvin–Mandel
matrix,

```math
C_{IJ} = \frac{M_{IJ}}{w_I\,w_J}\quad\text{for a stiffness},\qquad
S_{IJ} = M_{IJ}\,w_I\,w_J\quad\text{for a compliance} ,
```

so that ``C_{44} = C_{2323}`` and ``S_{44} = 4S_{2323}``, and
``\sigma_I = C_{IJ}\varepsilon_J``, ``\varepsilon_I = S_{IJ}\sigma_J`` with
``\varepsilon_4 = 2\varepsilon_{23}``. The two weightings differ: the inverse of
the Voigt matrix of a stiffness is the Voigt matrix of the compliance, not of
the inverse read as a stiffness ([`voigt_stiffness`](@ref),
[`voigt_compliance`](@ref)).

In full, the pair ``I \leftrightarrow (ij)`` is ``1 \leftrightarrow 11``,
``2 \leftrightarrow 22``, ``3 \leftrightarrow 33``, ``4 \leftrightarrow 23``,
``5 \leftrightarrow 13``, ``6 \leftrightarrow 12``, and a transversely isotropic
stiffness of axis ``\underline{e}_3`` has the two matrices

```math
\big[C_{IJ}\big] =
\begin{pmatrix}
C_{11} & C_{12} & C_{13} & 0 & 0 & 0\\
C_{12} & C_{11} & C_{13} & 0 & 0 & 0\\
C_{13} & C_{13} & C_{33} & 0 & 0 & 0\\
0 & 0 & 0 & C_{44} & 0 & 0\\
0 & 0 & 0 & 0 & C_{44} & 0\\
0 & 0 & 0 & 0 & 0 & C_{66}
\end{pmatrix},
\qquad
M =
\begin{pmatrix}
C_{11} & C_{12} & C_{13} & 0 & 0 & 0\\
C_{12} & C_{11} & C_{13} & 0 & 0 & 0\\
C_{13} & C_{13} & C_{33} & 0 & 0 & 0\\
0 & 0 & 0 & 2C_{44} & 0 & 0\\
0 & 0 & 0 & 0 & 2C_{44} & 0\\
0 & 0 & 0 & 0 & 0 & 2C_{66}
\end{pmatrix},
```

``C_{66} = (C_{11} - C_{12})/2``, ``C_{44} = C_{2323} = C_{1313}``. A component of a
general stiffness coupling a normal and a shear pair is ``\sqrt2`` times larger
in ``M`` (``M_{14} = \sqrt2\,C_{14}``). **Wherever this documentation writes a
stiffness component with two indices, ``C_{11}``, ``C_{44}``, ``C_{66}``, …, it
is the two-index (Voigt) component ``C_{IJ}`` of the left matrix**, the language
of the literature on waves and anisotropy; the Kelvin–Mandel matrix is the one
the library computes with. The two matrices of a transversely isotropic
stiffness given by its components:

```@example voigt_full
using TensND
C = ti_from_constants(10.0, 4.0, 3.0, 8.0, 2.5, [0, 0, 1]; form = :components)   # C₁₁, C₁₂, C₁₃, C₃₃, C₄₄
voigt_stiffness(C)
```

```@example voigt_full
Matrix(mandel(C))
```

The last line of the table is the decisive one for this library and is
developed next.

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

`TensND` builds it either from Euler angles (`mandel_rotation(θ, ϕ, ψ)`) or
directly from a frame (`_km_congruence`, in `src/tens_anisotropic.jl`); the latter
is what relates the two views of an orthotropic tensor,

```math
\mathrm{Mat}(\mathbb{C})=Q\;\mathrm{Mat}_{\text{mat}}(\mathbb{C})\;{}^{t}Q ,
```

with ``\mathrm{Mat}_{\text{mat}}`` the block-diagonal matrix in the material
frame ([Orthotropy](@ref th-orthotropy)).

!!! note "Consequence of the orthogonality"
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

## In two dimensions

The map is the same in any dimension ``d``: a symmetric tensor of order 2 has
``n = d(d+1)/2`` components, listed with the diagonal pairs first and the
off-diagonal ones from the last column back,

```math
(1,1),\ \ldots,\ (d,d),\ (d-1,d),\ (d-2,d),\ \ldots,\ (1,d),\ (d-2,d-1),\ \ldots,\ (1,2),
```

which is the order above for ``d = 3``, each off-diagonal one weighted by
``\sqrt2``. In the plane, ``n = 3`` and the pairs are ``11\to1``, ``22\to2``,
``12\to3``; the rotation of the basis by the angle ``\theta`` acts by the
orthogonal ``3\times3`` matrix

```math
Q=\begin{pmatrix}
\cos^2\theta & \sin^2\theta & -\sqrt2\,\cos\theta\sin\theta\\
\sin^2\theta & \cos^2\theta & \sqrt2\,\cos\theta\sin\theta\\
\sqrt2\,\cos\theta\sin\theta & -\sqrt2\,\cos\theta\sin\theta & \cos^2\theta-\sin^2\theta
\end{pmatrix} ,
```

and every statement of this page holds with ``6`` replaced by ``n``
([Tensors of the plane, of space and of any dimension](@ref th-dimensions)).
