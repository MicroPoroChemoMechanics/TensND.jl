# [Curvilinear differential calculus](@id th-curvilinear)

!!! tip "See also"
    **Theory**: [Classical coordinate systems](@ref th-predefined-systems) · [Spheroidal coordinates](@ref th-spheroidal) · [Ellipsoidal coordinates](@ref th-ellipsoidal)

    **Manual**: [Coordinate systems (symbolic)](@ref man-coorsystems) · [Coordinate systems (numerical)](@ref man-coorsystems-num) · [Charts and coordinate conversions](@ref man-charts)

    **Tutorials**: [Differential operators, symbolically](@ref tut-symbolic-operators) · [Christoffel symbols](@ref tut-christoffel) · [Differential operators, numerically](@ref tut-numerical-operators) · [Classical identities of the differential operators](@ref tut-operator-identities) · [The operators in the classical systems](@ref tut-operators-catalog) · [The spheroid in spheroidal coordinates](@ref tut-spheroidal-appendix) · [Ellipsoidal coordinates](@ref tut-ellipsoidal)

A coordinate system is built in four steps: a position vector as a function of
the coordinates, its derivatives (the natural basis), their norms and
directions (the Lamé coefficients and the normalized basis), and the rate at
which the basis turns (the Christoffel symbols). The differential operators
follow. This page defines each step and says how `TensND` carries it out,
symbolically (`src/coorsystems.jl`) or by automatic differentiation
(`src/coorsystems_num.jl`); the systems it predefines are on the next three
pages, [Classical coordinate systems](@ref th-predefined-systems),
[Spheroidal coordinates](@ref th-spheroidal) and
[Ellipsoidal coordinates](@ref th-ellipsoidal). A classical treatment of the
underlying tensor analysis is [simmonds1994](@citet).

## Chart, natural basis, metric

A coordinate system is a map ``\underline{OM}(q^1,\ldots,q^d)`` from a chart onto
``\mathbb{R}^d``. Its **natural basis** is the coordinate frame

```math
\underline{a}_i=\frac{\partial\underline{OM}}{\partial q^i},
```

with the dual basis ``\underline{a}^i`` and the metric
``g_{ij}=\underline{a}_i\cdot\underline{a}_j`` of
[Bases and variance](@ref th-bases-variance). Everything on that page applies
here verbatim; a coordinate system is a basis that varies from point to point.

## Lamé coefficients and the normalized basis

The natural vectors are generally neither unit nor orthogonal. Their norms are
the **Lamé coefficients**

```math
\chi_i=\|\underline{a}_i\|,
\qquad
\underline{e}_i=\frac{\underline{a}_i}{\chi_i}
\quad\text{(no summation)},
```

and ``(\underline{e}_i)`` is the **normalized basis** — the frame in which
physical components are read. For an *orthogonal* system the metric is
``g_{ij}=\mathrm{diag}(\chi_i^2)``, the dual vectors are
``\underline{a}^i=\underline{a}_i/\chi_i^2``, and the line element is

```math
\mathrm{d}s^2=\sum_i\chi_i^2\,(\mathrm{d}q^i)^2 .
```

`TensND` stores both bases: [`natural_basis`](@ref) and
[`normalized_basis`](@ref), with [`natvec`](@ref) and [`unitvec`](@ref) for
individual vectors and [`lame_coefficients`](@ref) for the ``\chi_i``.

## Christoffel symbols

The natural basis varies, and the rate at which it does is the connection:

```math
\frac{\partial\underline{a}_j}{\partial q^i}
=\Gamma^k_{ij}\,\underline{a}_k
\qquad\Longleftrightarrow\qquad
\Gamma^k_{ij}=\frac{\partial\underline{a}_j}{\partial q^i}\cdot\underline{a}^k .
```

They are symmetric in ``(i,j)`` because
``\partial_i\partial_j\underline{OM}=\partial_j\partial_i\underline{OM}``. The
**covariant derivative** of a vector field follows,

```math
\nabla_i v^k=\partial_i v^k+\Gamma^k_{ij}v^j ,
```

and it is the ``\Gamma`` terms — not the partial derivatives — that make a
divergence in spherical coordinates differ from a sum of derivatives.

The symbols of the classical systems are listed in
[Classical coordinate systems](@ref th-predefined-systems); all those of the
Cartesian system vanish, which is the definition of a Cartesian chart.

!!! note "Storage convention"
    [`christoffel`](@ref) returns a ``d\times d\times d`` array indexed
    ``\Gamma[i,j,k]=\Gamma^k_{ij}``: the **contravariant index is last**. The
    same convention is used by the `Γ_func` closure of
    the pointwise [`CoordinateSystem`](@ref).

## The five operators

`TensND` defines the differential operators intrinsically, as sums over the
chart directions:

| Operator | Definition | Order |
| :------- | :--------- | :---- |
| [`grad`](@ref) | ``\displaystyle\sum_i\partial_i t\otimes\underline{a}^i`` | ``+1`` |
| [`symgrad`](@ref) | ``\displaystyle\sum_i\partial_i t\stackrel{s}{\otimes}\underline{a}^i`` | ``+1`` |
| [`divergence`](@ref) | ``\displaystyle\sum_i\partial_i t\cdot\underline{a}^i`` | ``-1`` |
| [`laplacian`](@ref) | ``\mathrm{div}\circ\mathrm{grad}`` | ``0`` |
| [`hess`](@ref) | ``\mathrm{grad}\circ\mathrm{grad}`` | ``+2`` |

Here ``\partial_i t`` is the *full* derivative of the tensor field, basis
included — which is where the Christoffel symbols enter.

!!! warning "Index placement is a convention, and this is the one"
    ``\underline{a}^i`` is appended **on the right**, so for a vector field

    ```math
    (\nabla\underline{v})_{ij}=\partial_j v_i ,
    ```

    and `divergence` contracts the **last** index,

    ```math
    (\mathrm{div}\,\boldsymbol{\sigma})_i=\partial_j\sigma_{ij} .
    ```

    A library using the opposite convention differs by a transpose. The choice
    here is what makes ``\nabla^{2}=\mathrm{div}\circ\mathrm{grad}``
    correct as written. For a symmetric field the distinction is immaterial; for
    a general order-2 field it is not, and it is pinned by a test
    ([Testing and conventions](@ref dev-testing)).

Note also that `symgrad` is a primitive, not `(grad + transpose)/2` — applied to
a displacement field it is directly the linearized strain tensor.

## [Orthogonal systems](@id th-orthogonal-operators)

When the natural vectors are mutually orthogonal, as in every predefined
system, ``\underline{a}^i=\underline{e}_i/\chi_i`` and the operators only
involve the Lamé coefficients. For a scalar field ``f`` and a vector field
``\underline{v}=\sum_iv_i\,\underline{e}_i`` of physical components ``v_i``,

```math
\underline{\nabla}f=\sum_i\frac{1}{\chi_i}\,\frac{\partial f}{\partial q^i}\,\underline{e}_i,
\qquad
\mathrm{div}\,\underline{v}=\frac{1}{J}\sum_i\frac{\partial}{\partial q^i}\Bigl(\frac{J}{\chi_i}\,v_i\Bigr),
\qquad
\Delta f=\frac{1}{J}\sum_i\frac{\partial}{\partial q^i}\Bigl(\frac{J}{\chi_i^2}\,\frac{\partial f}{\partial q^i}\Bigr),
```

where ``J=\prod_i\chi_i`` is the Jacobian of the chart, so that the volume
element is ``\mathrm{d}\Omega=J\,\mathrm{d}q^1\cdots\mathrm{d}q^d``. The
divergence formula states that the flux of ``\underline{v}`` through the faces
of the cell ``\mathrm{d}q^1\cdots\mathrm{d}q^d`` is its divergence times the
volume of the cell: the face normal to ``\underline{e}_i`` has the area
``(J/\chi_i)\prod_{j\ne i}\mathrm{d}q^j``. The pages on the predefined systems
write these three formulas out.

## Building a coordinate system

`TensND` builds a system from what is known about it, by one of four routes.

| Given | Constructor | Computed |
| :---- | :---------- | :------- |
| nothing: a predefined system | [`coorsys_spherical`](@ref)`()`, [`coorsys_spheroidal`](@ref)`()`, [`coorsys_ellipsoidal`](@ref)`()`, … | — |
| ``\underline{OM}``, the normalized basis and the ``\chi_i`` in closed form | `CoordinateSystem(OM, coords, ℬ, χ)` | ``\Gamma`` |
| ``\underline{OM}`` alone, as an expression | `CoordinateSystem(OM, coords)` | ``\underline{a}_i``, ``\chi_i``, ``\underline{e}_i``, ``\Gamma`` |
| ``\underline{OM}`` alone, as a function | `CoordinateSystem(OM_func, dim)` | the same, at each point, by automatic differentiation |

The first three return an exact system, whose operators return expressions;
the last one a pointwise system, whose operators return functions of the
point. The predefined systems use the second route, since their bases and Lamé
coefficients are known in a simpler form than a computer algebra system would
find. [Coordinate systems (symbolic)](@ref man-coorsystems) and
[Coordinate systems (numerical)](@ref man-coorsystems-num) show each route.

## Two implementations, one interface

The operators have the same names and meaning whether the derivatives are taken
symbolically or by automatic differentiation.

| | [`CoordinateSystem`](@ref) | the pointwise [`CoordinateSystem`](@ref) |
| :--- | :--- | :--- |
| derivatives | SymPy / Symbolics, exact | `ForwardDiff`, machine precision |
| stores | ``\underline{OM}``, bases, ``\chi_i``, ``\Gamma`` as **expressions** | ``\chi``, ``R``, ``\Gamma`` as **closures**, evaluated pointwise |
| result | a formula, valid everywhere | a number, at one point |
| cost | grows with expression size | constant per point |
| needs | a symbolic ``\underline{OM}`` | any differentiable `OM(x)` |
| best for | deriving and simplifying closed forms | evaluating fields, sweeps, gradients |

A [`PointwiseGeometry`](@ref) stores the geometry as three closures — `χ_func` (Lamé),
`R_func` (the rotation to the normalized frame) and `Γ_func` (Christoffel, with
the ``\Gamma[i,j,k]=\Gamma^k_{ij}`` convention) — so a custom system needs only a
differentiable position map, with no symbolic setup at all.

Because everything downstream is generic in the element type,
`ForwardDiff.Dual` numbers flow through the operators in **both** directions: an
operator result can be differentiated again with respect to the evaluation
point, and with respect to a parameter carried by the field itself. The internal
buffers promote the element type of the geometry (the Lamé coefficients) with
that of the field value, so a `Dual`-valued field evaluated on a `Float64`
coordinate system is handled correctly.

## Simplification machinery

Non-trivial symbolic charts — the spheroidal and the ellipsoidal ones are the
standard examples — produce metric expressions that SymPy will not reduce on
its own. The exact `CoordinateSystem`
therefore accepts:

| Argument | Role |
| :------- | :--- |
| `tmp_coords` | auxiliary variables standing for compound expressions |
| `params` | constants appearing in ``\underline{OM}`` |
| `rules` | rewrite rules applied after each simplification |
| `tmp_var` | substitutions replacing expressions by temporary variables |
| `to_coords` | how to eliminate the temporaries before differentiating |

The spheroidal systems name ``\bar{p}=\sqrt{1-p^2}`` and
``\bar{q}=\sqrt{q^2-1}``, the ellipsoidal one its six square roots: the
position vector is written with these symbols, and `to_coords` restores the
square roots before a derivative is taken. Their use is shown in
[Christoffel symbols](@ref tut-christoffel) among the tutorials.
