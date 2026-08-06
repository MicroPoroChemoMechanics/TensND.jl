# [The extended Walpole algebra](@id th-walpole-extended)

The classical Walpole basis of [The Walpole basis](@ref th-walpole) spans the
tensors that are transversely isotropic in the usual sense. That space is *not*
the full space of tensors invariant under rotations about the axis: two
generators are missing. This page derives the eight-dimensional space `TensND`
actually implements (`TensTI{4,T,8}`, `src/tens_anisotropic.jl`), which has no
counterpart in the Echoes manual or, as far as we are aware, in the standard
references.

## Why eight and not six

Fix a unit axis ``\underline{n}`` and let ``SO(2)`` act by rotations about it.
The space of minor-symmetric order-4 tensors is six-dimensional in the
Kelvin–Mandel picture, and the tensors invariant under the action form the
**commutant** of that representation. Decomposing the ``6``-dimensional
Kelvin–Mandel space into isotypic components of ``SO(2)``:

| Harmonic | Subspace | Dimension | Commutant | Parameters |
| :------- | :------- | :-------- | :-------- | :--------- |
| ``m=0`` | axial ``\varepsilon_{nn}`` and in-plane spherical | 2 | ``2\times2`` real matrices | ``\ell_1,\ell_2,\ell_3,\ell_4`` |
| ``m=1`` | axial shears ``(23,13)`` | 2 | ``\cong\mathbb{C}`` | ``z_1=\ell_6+i\,\ell_7`` |
| ``m=2`` | in-plane deviatoric ``(11{-}22,\;12)`` | 2 | ``\cong\mathbb{C}`` | ``z_2=\ell_5+i\,\ell_8`` |

By Schur's lemma the commutant of each isotypic block is a division algebra —
``\mathbb{R}^{2\times2}`` for the trivial (``m=0``) block, ``\mathbb{C}`` for
each two-dimensional real irreducible — giving ``4+2+2=8`` real parameters.
The classical Walpole basis captures the real parts ``\ell_5``, ``\ell_6`` and
drops the imaginary parts ``\ell_7``, ``\ell_8``.

## The two extra generators

Let ``\boldsymbol{w}`` be the in-plane rotation generator,

```math
\boldsymbol{w}\cdot\underline{p}=\underline{n}\times\underline{p},
\qquad
w_{ij}=\varepsilon_{ikj}\,n_k ,
```

which is antisymmetric and odd in ``\underline{n}``. Then

```math
\begin{aligned}
(\mathbb{W}_7)_{ijkl}&=-\tfrac{1}{2}\bigl(
 w_{ik}(\boldsymbol{1}_n)_{jl}+w_{il}(\boldsymbol{1}_n)_{jk}
+w_{jk}(\boldsymbol{1}_n)_{il}+w_{jl}(\boldsymbol{1}_n)_{ik}\bigr)\\[4pt]
(\mathbb{W}_8)_{ijkl}&=+\tfrac{1}{4}\bigl(
 w_{ik}(\boldsymbol{1}_T)_{jl}+w_{il}(\boldsymbol{1}_T)_{jk}
+w_{jk}(\boldsymbol{1}_T)_{il}+w_{jl}(\boldsymbol{1}_T)_{ik}\bigr)
\end{aligned}
```

In the Kelvin–Mandel frame with ``\underline{n}=\underline{e}_3``:

```math
\mathrm{Mat}(\mathbb{W}_7)=
\begin{pmatrix}0&&&&&\\&0&&&&\\&&0&&&\\&&&0&-1&\\&&&1&0&\\&&&&&0\end{pmatrix}
\qquad
\mathrm{Mat}(\mathbb{W}_8)=
\begin{pmatrix}0&0&0&0&0&-\tfrac{1}{\sqrt2}\\
0&0&0&0&0&\tfrac{1}{\sqrt2}\\
0&0&0&0&0&0\\ 0&0&0&0&0&0\\ 0&0&0&0&0&0\\
\tfrac{1}{\sqrt2}&-\tfrac{1}{\sqrt2}&0&0&0&0\end{pmatrix}
```

Both matrices are **antisymmetric**, so

```math
{}^{t}\mathbb{W}_7=-\mathbb{W}_7,
\qquad
{}^{t}\mathbb{W}_8=-\mathbb{W}_8 :
```

they are *major-antisymmetric*. They are therefore invisible to any
major-symmetric description — which is exactly why they are absent from the
classical basis, whose target was the elastic stiffness.

## The algebra

The eight-dimensional space is a commutant, hence an **algebra**: closed under
double contraction and under inversion. In terms of the block and the two
complex numbers,

```math
\mathbb{L}\equiv\bigl(L,\;z_1,\;z_2\bigr),
\qquad
L=\begin{pmatrix}\ell_1&\ell_3\\\ell_4&\ell_2\end{pmatrix},
\quad
z_1=\ell_6+i\ell_7,
\quad
z_2=\ell_5+i\ell_8,
```

```math
\mathbb{L}:\mathbb{M}\equiv\bigl(L\,M,\;z_1 z'_1,\;z_2 z'_2\bigr),
\qquad
\mathbb{L}^{-1}\equiv\Bigl(L^{-1},\;\tfrac{1}{z_1},\;\tfrac{1}{z_2}\Bigr),
```

with an ordinary ``2\times2`` matrix product and inverse, and ordinary
**complex** products and inverses. Setting ``\ell_7=\ell_8=0`` recovers the real
rules of [The Walpole basis](@ref th-walpole) exactly.

## What the extra generators do to order-2 contraction

``\mathbb{W}_7`` and ``\mathbb{W}_8`` annihilate precisely the order-2 tensors
that are themselves axially invariant:

```math
\mathbb{W}_7:\bigl(a\,\boldsymbol{1}_T+b\,\boldsymbol{1}_n\bigr)=
\mathbb{W}_8:\bigl(a\,\boldsymbol{1}_T+b\,\boldsymbol{1}_n\bigr)=0 .
```

!!! warning "They do not annihilate a general symmetric tensor"
    A common shortcut is to assume ``\ell_7`` and ``\ell_8`` never matter for
    order-4 : order-2 contraction. They do. ``\mathbb{W}_7`` maps the axial-shear
    pair ``(\varepsilon_{23},\varepsilon_{13})`` onto itself rotated by a quarter
    turn, and ``\mathbb{W}_8`` does the same for the in-plane deviatoric pair
    ``(\varepsilon_{11}-\varepsilon_{22},\;\varepsilon_{12})``. Both results are
    symmetric and generally nonzero.

    What survives is the weaker statement above, and it is enough for the rules
    that involve only ``\boldsymbol{1}_T`` and ``\boldsymbol{1}_n`` — which is
    why the ``\ell_1..\ell_4`` machinery carries over unchanged to ``N=8``. The
    ``N=8`` contraction methods do account for ``\ell_7,\ell_8``, and are checked
    against the generic dense route.

## Order 2: three generators, not two

The same argument one order down. The space of order-2 tensors invariant under
rotations about ``\underline{n}`` is **three**-dimensional:

```math
\boldsymbol{a}=a\,\boldsymbol{1}_T+b\,\boldsymbol{1}_n+c\,\boldsymbol{w},
```

the third generator being the antisymmetric in-plane rotation ``\boldsymbol{w}``
that a symmetric parametrization cannot represent. With
``\underline{n}=\underline{e}_3``:

```math
\boldsymbol{a}=
\begin{pmatrix}a&-c&0\\ c&a&0\\ 0&0&b\end{pmatrix}.
```

This is `TensTI{2,T,3}`; the familiar symmetric case is `TensTI{2,T,2}` with
``c=0``.

## When this matters

For an elastic stiffness or compliance — major-symmetric by construction — the
extra coefficients vanish and the classical five- or six-parameter description is
complete. They become essential for objects that are **not** major-symmetric,
the archetype being a strain-concentration tensor

```math
\mathbb{A}=\bigl[\mathbb{I}+\mathbb{P}:(\mathbb{C}_1-\mathbb{C}_0)\bigr]^{-1},
```

whose exact average over rotations about an axis lands in the full
eight-dimensional space. Projecting such an object onto the five-parameter
major-symmetric subspace — which is what a *best-fit* projection does, see
[Projection onto a symmetry class](@ref th-projection) — silently discards the
``\ell_3\neq\ell_4`` split together with ``\ell_7`` and ``\ell_8``. Averaging and
projecting are different operations, and this space is where the difference
lives.

Accessors: [`get_ℓ`](@ref) returns the six classical coefficients (dropping
``\ell_7,\ell_8``), [`get_ℓ8`](@ref) always returns the full eight-tuple, padding
with zeros for `N=5` and `N=6` inputs.
