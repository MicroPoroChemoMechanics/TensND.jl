# [The Walpole basis](@id th-walpole)

Transverse isotropy about an axis ``\underline{n}`` is the symmetry class where a
well-chosen basis pays off most: in the Walpole basis [walpole1984](@cite) the
double contraction becomes a ``2\times2`` matrix product plus two scalar
products, and inversion becomes a ``2\times2`` matrix inverse plus two
reciprocals. `TensND` stores such tensors as their coefficients on this basis
([`TensTI`](@ref), `src/tens_walpole.jl`) rather than as 81 components.

## The six tensors

From the axial and transverse projectors

```math
\boldsymbol{1}_n=\underline{n}\otimes\underline{n},
\qquad
\boldsymbol{1}_T=\boldsymbol{1}-\boldsymbol{1}_n ,
```

build [walpole1984](@cite), [walpole1981](@cite):

```math
\begin{aligned}
\mathbb{W}_1&=\boldsymbol{1}_n\otimes\boldsymbol{1}_n
&\mathbb{W}_2&=\tfrac{1}{2}\,\boldsymbol{1}_T\otimes\boldsymbol{1}_T\\[2pt]
\mathbb{W}_3&=\tfrac{1}{\sqrt2}\,\boldsymbol{1}_n\otimes\boldsymbol{1}_T
&\mathbb{W}_4&=\tfrac{1}{\sqrt2}\,\boldsymbol{1}_T\otimes\boldsymbol{1}_n\\[2pt]
\mathbb{W}_5&=\boldsymbol{1}_T\stackrel{s}{\boxtimes}\boldsymbol{1}_T
              -\tfrac{1}{2}\,\boldsymbol{1}_T\otimes\boldsymbol{1}_T
&\mathbb{W}_6&=\boldsymbol{1}_T\stackrel{s}{\boxtimes}\boldsymbol{1}_n
              +\boldsymbol{1}_n\stackrel{s}{\boxtimes}\boldsymbol{1}_T
\end{aligned}
```

and any transversely isotropic order-4 tensor decomposes as

```math
\mathbb{L}=\sum_{i=1}^{6}\ell_i\,\mathbb{W}_i .
```

They are [`tens_W1`](@ref) … [`tens_W6`](@ref), collectively
[`walpole_basis`](@ref).

## Kelvin–Mandel matrices

In the frame where ``\underline{n}=\underline{e}_3``, with the ordering of
[Kelvin–Mandel representation](@ref th-kelvin-mandel):

```math
\mathrm{Mat}(\mathbb{W}_1)=\!\begin{pmatrix}0&&&&&\\&0&&&&\\&&1&&&\\&&&0&&\\&&&&0&\\&&&&&0\end{pmatrix}
\;\;
\mathrm{Mat}(\mathbb{W}_2)=\!\begin{pmatrix}\frac12&\frac12&&&&\\\frac12&\frac12&&&&\\&&0&&&\\&&&0&&\\&&&&0&\\&&&&&0\end{pmatrix}
\;\;
\mathrm{Mat}(\mathbb{W}_3)=\!\begin{pmatrix}0&0&0&&&\\0&0&0&&&\\\frac{1}{\sqrt2}&\frac{1}{\sqrt2}&0&&&\\&&&0&&\\&&&&0&\\&&&&&0\end{pmatrix}
```

```math
\mathrm{Mat}(\mathbb{W}_4)=\!\begin{pmatrix}0&0&\frac{1}{\sqrt2}&&&\\0&0&\frac{1}{\sqrt2}&&&\\0&0&0&&&\\&&&0&&\\&&&&0&\\&&&&&0\end{pmatrix}
\;\;
\mathrm{Mat}(\mathbb{W}_5)=\!\begin{pmatrix}\frac12&-\frac12&&&&\\-\frac12&\frac12&&&&\\&&0&&&\\&&&0&&\\&&&&0&\\&&&&&1\end{pmatrix}
\;\;
\mathrm{Mat}(\mathbb{W}_6)=\!\begin{pmatrix}0&&&&&\\&0&&&&\\&&0&&&\\&&&1&&\\&&&&1&\\&&&&&0\end{pmatrix}
```

Note that ``\mathbb{W}_3`` and ``\mathbb{W}_4`` are transposes of each other and
are **not** major-symmetric individually; every other ``\mathbb{W}_i`` is.

## The synthetic triplet

Gather the first four coefficients into a ``2\times2`` matrix and keep the last
two apart:

```math
\mathbb{L}\equiv\bigl(L,\ \ell_5,\ \ell_6\bigr),
\qquad
L=\begin{pmatrix}\ell_1&\ell_3\\ \ell_4&\ell_2\end{pmatrix}.
```

The Walpole basis is closed under double contraction, and in this notation the
algebra is [walpole1984](@cite)

```math
\mathbb{L}:\mathbb{M}\equiv\bigl(L\,M,\ \ell_5 m_5,\ \ell_6 m_6\bigr),
\qquad
\mathbb{L}^{-1}\equiv\Bigl(L^{-1},\ \tfrac{1}{\ell_5},\ \tfrac{1}{\ell_6}\Bigr).
```

The reason is the multiplication table, which identifies
``(\mathbb{W}_1,\mathbb{W}_2,\mathbb{W}_3,\mathbb{W}_4)`` with the matrix units
``(E_{11},E_{22},E_{12},E_{21})`` of ``2\times2`` matrices, while
``\mathbb{W}_5`` and ``\mathbb{W}_6`` are idempotents spanning two independent
one-dimensional ideals. Every non-vanishing product is:

| ``:`` | ``\mathbb{W}_1`` | ``\mathbb{W}_2`` | ``\mathbb{W}_3`` | ``\mathbb{W}_4`` | ``\mathbb{W}_5`` | ``\mathbb{W}_6`` |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| ``\mathbb{W}_1`` | ``\mathbb{W}_1`` | ``0`` | ``\mathbb{W}_3`` | ``0`` | ``0`` | ``0`` |
| ``\mathbb{W}_2`` | ``0`` | ``\mathbb{W}_2`` | ``0`` | ``\mathbb{W}_4`` | ``0`` | ``0`` |
| ``\mathbb{W}_3`` | ``0`` | ``\mathbb{W}_3`` | ``0`` | ``\mathbb{W}_1`` | ``0`` | ``0`` |
| ``\mathbb{W}_4`` | ``\mathbb{W}_4`` | ``0`` | ``\mathbb{W}_2`` | ``0`` | ``0`` | ``0`` |
| ``\mathbb{W}_5`` | ``0`` | ``0`` | ``0`` | ``0`` | ``\mathbb{W}_5`` | ``0`` |
| ``\mathbb{W}_6`` | ``0`` | ``0`` | ``0`` | ``0`` | ``0`` | ``\mathbb{W}_6`` |

(row ``i``, column ``j`` is ``\mathbb{W}_i:\mathbb{W}_j``; the table is not
symmetric.)

## Orthogonality and norms

For the Frobenius scalar product the basis is **orthogonal but not
orthonormal**:

```math
\langle\mathbb{W}_i,\mathbb{W}_j\rangle=\mathbb{W}_i::\mathbb{W}_j
= g_i\,\delta_{ij},
\qquad
(g_i)=(1,\,1,\,1,\,1,\,2,\,2) .
```

This single fact is what makes every transversely isotropic projection
closed-form: the normal equations of
[Projection onto a symmetry class](@ref th-projection) are diagonal, so each
coefficient is an independent quotient
``\ell_i=\langle\mathbb{A},\mathbb{W}_i\rangle/g_i``.

## The isotropic tensors on this basis

The identity, spherical and deviatoric tensors of
[Isotropic tensors](@ref th-isotropic) are transversely isotropic about *any*
axis, so they have Walpole coefficients — and these are worth recording because
they are easy to get wrong:

```math
\mathbb{I}=\mathbb{W}_1+\mathbb{W}_2+\mathbb{W}_5+\mathbb{W}_6
```

```math
\mathbb{J}=\tfrac{1}{3}\Bigl(\mathbb{W}_1+2\,\mathbb{W}_2
            +\sqrt2\,\mathbb{W}_3+\sqrt2\,\mathbb{W}_4\Bigr)
```

```math
\mathbb{K}=\tfrac{1}{3}\Bigl(2\,\mathbb{W}_1+\mathbb{W}_2
            -\sqrt2\,\mathbb{W}_3-\sqrt2\,\mathbb{W}_4\Bigr)
            +\mathbb{W}_5+\mathbb{W}_6
```

or, in the synthetic notation,

| | ``L`` | ``\ell_5`` | ``\ell_6`` |
| :--- | :--- | :--- | :--- |
| ``\mathbb{I}`` | ``\begin{pmatrix}1&0\\0&1\end{pmatrix}`` | ``1`` | ``1`` |
| ``\mathbb{J}`` | ``\frac13\begin{pmatrix}1&\sqrt2\\ \sqrt2&2\end{pmatrix}`` | ``0`` | ``0`` |
| ``\mathbb{K}`` | ``\frac13\begin{pmatrix}2&-\sqrt2\\ -\sqrt2&1\end{pmatrix}`` | ``1`` | ``1`` |

Both ``\mathbb{J}`` and ``\mathbb{K}`` have singular ``L`` blocks
(``\det L=0``), as they must, being projectors. ``\mathbb{J}+\mathbb{K}=\mathbb{I}``
is immediate on the table. Converting an isotropic tensor to this
representation is [`fromISO`](@ref).

!!! warning "``\mathbb{I}`` is not the sum of the six"
    ``\sum_i\mathbb{W}_i\neq\mathbb{I}``: the sum overshoots by
    ``\mathbb{W}_3+\mathbb{W}_4``, whose norm is ``\sqrt2``. Likewise
    ``\mathbb{J}\neq\mathbb{W}_1+\mathbb{W}_2`` and
    ``\mathbb{K}\neq\mathbb{W}_3+\mathbb{W}_4+\mathbb{W}_5+\mathbb{W}_6``. The
    three correct identities above are pinned by tests
    ([Testing and conventions](@ref dev-testing)).

## Major symmetry and storage

```math
\mathbb{L}={}^{t}\mathbb{L}
\quad\Longleftrightarrow\quad
\ell_3=\ell_4
\quad\Longleftrightarrow\quad
L={}^{t}L .
```

A major-symmetric transversely isotropic tensor — an elastic stiffness, for
instance — therefore needs **five** coefficients, a general one **six**:

| Storage | Coefficients | Case |
| :------ | :----------- | :--- |
| `TensTI{4,T,5}` | ``(\ell_1,\ell_2,\ell_3,\ell_5,\ell_6)`` | major-symmetric, ``\ell_4=\ell_3`` |
| `TensTI{4,T,6}` | ``(\ell_1,\ldots,\ell_6)`` | general |
| `TensTI{4,T,8}` | ``(\ell_1,\ldots,\ell_8)`` | full axially-invariant space, see [The extended Walpole algebra](@ref th-walpole-extended) |

The widening is not cosmetic: the product of two *major-symmetric* Walpole
tensors is generally **not** major-symmetric, because ``L\,M\neq{}^{t}(L\,M)``
unless ``L`` and ``M`` commute. `TensTI{4,T,5} ⊡ TensTI{4,T,5}` therefore
returns an `N=6` container.

The symmetrized basis sometimes met in the literature merges the pair,
``\mathbb{W}^s_3=\mathbb{W}_3+\mathbb{W}_4``, and relabels
``\mathbb{W}^s_4=\mathbb{W}_5``, ``\mathbb{W}^s_5=\mathbb{W}_6``; it spans
exactly the five-dimensional major-symmetric subspace and is
[`walpole_basis_sym`](@ref).

## Order-2 transverse isotropy

The same construction one order down: a transversely isotropic order-2 tensor is

```math
\boldsymbol{a}=a\,\boldsymbol{1}_T+b\,\boldsymbol{1}_n ,
```

with ``a`` the transverse and ``b`` the axial coefficient, stored as
`TensTI{2,T,2}`. Products and inverses are again termwise, since
``\boldsymbol{1}_T`` and ``\boldsymbol{1}_n`` are complementary orthogonal
projectors. The tensor is isotropic exactly when ``a=b``. A third generator
exists and is discussed on
[The extended Walpole algebra](@ref th-walpole-extended).
