# [Toroidal functions](@id th-toroidal-harmonics)

!!! tip "See also"
    **Theory**: [Toroidal coordinates](@ref th-toroidal) · [Curvilinear differential calculus](@ref th-curvilinear)

    **Manual**: [Coordinate systems (symbolic)](@ref man-coorsystems) · [Charts and coordinate conversions](@ref man-charts) · [Toroidal functions](@ref man-toroidal-harmonics)

The harmonic functions that separate in toroidal coordinates
([Toroidal coordinates](@ref th-toroidal)) are built on Legendre functions of
half-integer degree and of argument larger than one, the **toroidal
functions**. This page defines them, gives the relations that the package uses
to evaluate them, and explains how the evaluation stays accurate at every
degree. The numbered equations cited are those of the
[NIST Digital Library of Mathematical Functions](@cite dlmf).

## Definitions

Let ``x=\cosh\alpha\ge1``, ``w=\sinh\alpha=\sqrt{x^2-1}`` and
``\nu=n-\tfrac12`` with ``n`` an integer. ``P_\nu`` is the Legendre function of
the first kind, the solution of Legendre's equation of degree ``\nu`` that
equals 1 at ``x=1``; ``Q_\nu`` is that of the second kind, the solution that
vanishes at infinity. For an integer order ``m\ge0``,

```math
P_\nu^m(x)=w^m\,\frac{\mathrm{d}^mP_\nu}{\mathrm{d}x^m},
\qquad
Q_\nu^m(x)=w^m\,\frac{\mathrm{d}^mQ_\nu}{\mathrm{d}x^m} .
```

``P_\nu^m`` is the function of the same name of the reference. Its function of
the second kind carries the factor ``\mathrm{e}^{\mathrm{i}m\pi}``, absent
here, and its ``\boldsymbol{Q}_\nu^m`` a further division by
``\Gamma(\nu+m+1)`` (14.3.10): at the order 0,
``\boldsymbol{Q}_\nu=Q_\nu/\Gamma(\nu+1)``.

A harmonic function of the toroidal system is a combination of the products

```math
\sqrt{\cosh\alpha-\cos\beta}\;F_{n-1/2}^m(\cosh\alpha)\;
\Bigl\{\cos n\beta,\ \sin n\beta\Bigr\}\;\Bigl\{\cos m\varphi,\ \sin m\varphi\Bigr\},
\qquad F=P\ \text{or}\ Q .
```

## Two integrals

Two integrals over an angle give the functions of the order 0, and serve as the
independent references of the tests. The first is Laplace's (14.12.7 at the
order 0):

```math
P_{n-1/2}(\cosh\alpha)=\frac{1}{\pi}\int_0^\pi\bigl(\cosh\alpha+\sinh\alpha\,\cos t\bigr)^{n-1/2}\,\mathrm{d}t .
```

The second comes from the expansion of the inverse distance (14.19.6 at the
order 0, with ``\boldsymbol{Q}_{n-1/2}=Q_{n-1/2}/\Gamma(n+\tfrac12)``):

```math
\frac{1}{\sqrt{2\cosh\alpha-2\cos\beta}}
=\frac{1}{\pi}\sum_{n=0}^{\infty}\epsilon_n\,Q_{n-1/2}(\cosh\alpha)\,\cos n\beta,
\qquad \epsilon_0=1,\quad \epsilon_{n\ge1}=2 .
```

The functions of the second kind are thus the Fourier coefficients of the left
side,

```math
Q_{n-1/2}(\cosh\alpha)=\int_0^\pi\frac{\cos n\beta}{\sqrt{2\cosh\alpha-2\cos\beta}}\,\mathrm{d}\beta .
```

Read the other way, the expansion says how the constant 1, a harmonic function,
is written in the toroidal system:

```math
1=\frac{\sqrt2}{\pi}\,\sqrt{\cosh\alpha-\cos\beta}\,\sum_{n=0}^{\infty}\epsilon_n\,Q_{n-1/2}(\cosh\alpha)\,\cos n\beta .
```

## Relations in the degree and in the order

Both kinds satisfy the same three relations. In the degree (14.10.3),

```math
(\nu-m+1)\,F_{\nu+1}^m=(2\nu+1)\,x\,F_\nu^m-(\nu+m)\,F_{\nu-1}^m ;
```

for the derivative (14.10.5),

```math
(x^2-1)\,\frac{\mathrm{d}F_\nu^m}{\mathrm{d}x}=\nu\,x\,F_\nu^m-(\nu+m)\,F_{\nu-1}^m ;
```

and in the order. The last one follows from the definition: since
``\mathrm{d}w/\mathrm{d}x=x/w``,

```math
\frac{\mathrm{d}F_\nu^m}{\mathrm{d}x}=\frac{m\,x}{w^2}\,F_\nu^m+\frac{1}{w}\,F_\nu^{m+1}
\qquad\Longrightarrow\qquad
w\,F_\nu^{m+1}=(\nu-m)\,x\,F_\nu^m-(\nu+m)\,F_{\nu-1}^m ,
```

after replacing the derivative by its expression. With
``\mathrm{d}/\mathrm{d}\alpha=w\,\mathrm{d}/\mathrm{d}x``, the derivative with
respect to the coordinate is

```math
\frac{\mathrm{d}}{\mathrm{d}\alpha}F_\nu^m(\cosh\alpha)=\frac{\nu\,\cosh\alpha\,F_\nu^m-(\nu+m)\,F_{\nu-1}^m}{\sinh\alpha} .
```

Every relation calls the degree below. At ``n=0`` it is ``\nu-1=-\tfrac32``,
and the functions are even in ``n``:

```math
F_{-n-1/2}^m=F_{n-1/2}^m .
```

For the first kind this is ``P_{-\nu-1}^m=P_\nu^m`` (14.9.11), true of every
degree. For the second kind it is special to the half-integer degrees:
``Q_{-\nu-1}-Q_\nu`` is proportional to ``\cos(\nu\pi)\,P_\nu`` (14.9.12), and
``\cos\bigl((n-\tfrac12)\pi\bigr)=0``.

## The first two degrees

The recurrence in the degree starts from ``n=0`` and ``n=1``, which are
complete elliptic integrals ``K`` and ``E`` of the modulus written as their
argument (14.5.24 to 14.5.27):

```math
P_{-1/2}(\cosh\alpha)=\frac{2}{\pi\,\cosh\frac{\alpha}{2}}\,K\Bigl(\tanh\frac{\alpha}{2}\Bigr),
\qquad
P_{1/2}(\cosh\alpha)=\frac{2}{\pi}\,\mathrm{e}^{\alpha/2}\,E\Bigl(\sqrt{1-\mathrm{e}^{-2\alpha}}\Bigr),
```

```math
Q_{-1/2}(\cosh\alpha)=2\,\mathrm{e}^{-\alpha/2}\,K\bigl(\mathrm{e}^{-\alpha}\bigr),
\qquad
Q_{1/2}(\cosh\alpha)=2\,\mathrm{e}^{\alpha/2}\,\Bigl[K\bigl(\mathrm{e}^{-\alpha}\bigr)-E\bigl(\mathrm{e}^{-\alpha}\bigr)\Bigr] .
```

The reference gives ``Q_{1/2}`` with the modulus ``\mathrm{sech}(\alpha/2)``,
as the difference of two terms that each grow as ``\mathrm{e}^{\alpha/2}``
while the function decays as ``\mathrm{e}^{-3\alpha/2}``. The form above, with
the modulus of ``Q_{-1/2}``, is the Fourier coefficient ``n=1`` of the previous
section; the tests compare it with that integral.

The integrals are computed by the arithmetic-geometric mean. Starting from
``a_0=1``, ``b_0=k'=\sqrt{1-k^2}`` and ``c_0=k``, the sequences
``a_{j+1}=(a_j+b_j)/2``, ``b_{j+1}=\sqrt{a_jb_j}``, ``c_{j+1}=(a_j-b_j)/2``
converge quadratically, and

```math
K(k)=\frac{\pi}{2\,a_\infty},
\qquad
K(k)-E(k)=K(k)\sum_{j\ge0}2^{\,j-1}\,c_j^2 .
```

The difference ``K-E`` is a sum of positive terms: it is obtained without
subtracting two close numbers, which is what ``Q_{1/2}`` needs when ``\alpha``
is large.

## Which way to run the recurrence

For a given ``\alpha>0`` the recurrence in the degree has two independent
solutions, which behave at large ``n`` as

```math
P_{n-1/2}(\cosh\alpha)\sim\mathrm{e}^{(n-1/2)\alpha},
\qquad
Q_{n-1/2}(\cosh\alpha)\sim\mathrm{e}^{-(n+1/2)\alpha},
```

up to factors that vary as powers of ``n``. Run upward, the recurrence
amplifies any rounding error as the growing solution. This is harmless for
``P``, which is that solution. For ``Q`` an error ``\varepsilon`` made at the
start becomes, relative to ``Q_{n-1/2}``, of the order of
``\varepsilon\,\mathrm{e}^{2n\alpha}``.

The functions of the second kind are therefore computed in two ways:

- where ``\mathrm{e}^{2N\alpha}\le16`` over the degrees asked for, upward from
  the two first degrees, which loses at most four bits;
- elsewhere, through the ratios ``r_n=Q_{n-1/2}/Q_{n-3/2}``, which satisfy

  ```math
  r_n=\frac{\nu}{2n\,x-(\nu+1)\,r_{n+1}},\qquad \nu=n-\tfrac12 .
  ```

  Run downward from ``r=\mathrm{e}^{-\alpha}`` at a degree ``L`` above ``N``,
  this relation damps the error of its start by ``\mathrm{e}^{-2\alpha}`` at
  each step; ``L`` is taken so that it falls below the rounding of the number
  type. The values follow from ``Q_{-1/2}`` by the products of the ratios.

The orders ``m\ge1`` are then reached by the relation in the order, at all the
degrees at once. It divides by ``\sinh\alpha``: for ``P``, whose terms nearly
cancel near the axis, the relative accuracy decreases as
``1/\alpha^2`` when ``\alpha\to0``.

The stopping tests of these iterations bear on values. The functions accept
numbers and dual numbers, whose derivatives converge with the values, and
refuse symbols.

[`toroidal_P`](@ref), [`toroidal_Q`](@ref), [`toroidal_dP`](@ref) and
[`toroidal_dQ`](@ref) return the functions and their derivatives with respect
to ``\alpha`` for all the degrees up to ``N``.
