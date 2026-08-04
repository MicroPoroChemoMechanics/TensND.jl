# [Curvilinear differential calculus](@id th-curvilinear)

The second pillar of the library, and the second chapter with no counterpart in
the [Echoes manual](https://jfbarthelemy.github.io/echoes/), which works only in
Cartesian coordinates. Everything below is derived from the definitions and
matches `src/coorsystems.jl` (symbolic) and `src/coorsystems_num.jl`
(automatic differentiation). A classical treatment of the underlying tensor
analysis is [simmonds1994](@cite).

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
individual vectors and [`Lame`](@ref) for the ``\chi_i``.

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

!!! note "Storage convention"
    [`Christoffel`](@ref) returns a ``d\times d\times d`` array indexed
    ``\Gamma[i,j,k]=\Gamma^k_{ij}``: the **contravariant index is last**. The
    same convention is used by the `Γ_func` closure of
    [`CoorSystemNum`](@ref).

## The five operators

`TensND` defines the differential operators intrinsically, as sums over the
chart directions:

| Operator | Definition | Order |
| :------- | :--------- | :---- |
| [`GRAD`](@ref) | ``\displaystyle\sum_i\partial_i t\otimes\underline{a}^i`` | ``+1`` |
| [`SYMGRAD`](@ref) | ``\displaystyle\sum_i\partial_i t\stackrel{s}{\otimes}\underline{a}^i`` | ``+1`` |
| [`DIV`](@ref) | ``\displaystyle\sum_i\partial_i t\cdot\underline{a}^i`` | ``-1`` |
| [`LAPLACE`](@ref) | ``\mathrm{DIV}\circ\mathrm{GRAD}`` | ``0`` |
| [`HESS`](@ref) | ``\mathrm{GRAD}\circ\mathrm{GRAD}`` | ``+2`` |

Here ``\partial_i t`` is the *full* derivative of the tensor field, basis
included — which is where the Christoffel symbols enter.

!!! warning "Index placement is a convention, and this is the one"
    ``\underline{a}^i`` is appended **on the right**, so for a vector field

    ```math
    (\nabla\underline{v})_{ij}=\partial_j v_i ,
    ```

    and `DIV` contracts the **last** index,

    ```math
    (\mathrm{DIV}\,\boldsymbol{\sigma})_i=\partial_j\sigma_{ij} .
    ```

    A library using the opposite convention differs by a transpose. The choice
    here is what makes ``\mathrm{LAPLACE}=\mathrm{DIV}\circ\mathrm{GRAD}``
    correct as written. For a symmetric field the distinction is immaterial; for
    a general order-2 field it is not, and it is pinned by a test
    ([Testing and conventions](@ref dev-testing)).

Note also that `SYMGRAD` is a primitive, not `(GRAD + transpose)/2` — applied to
a displacement field it is directly the linearized strain tensor.

## Predefined systems

| System | Coordinates | ``\underline{OM}`` | Lamé ``\chi_i`` |
| :----- | :---------- | :----------------- | :-------------- |
| [`coorsys_cartesian`](@ref) | ``(x,y,z)`` | ``(x,\,y,\,z)`` | ``(1,1,1)`` |
| [`coorsys_polar`](@ref) | ``(r,\theta)`` | ``(r\cos\theta,\;r\sin\theta)`` | ``(1,\,r)`` |
| [`coorsys_cylindrical`](@ref) | ``(r,\theta,z)`` | ``(r\cos\theta,\;r\sin\theta,\;z)`` | ``(1,\,r,\,1)`` |
| [`coorsys_spherical`](@ref) | ``(\theta,\varphi,r)`` | ``(r\sin\theta\cos\varphi,\;r\sin\theta\sin\varphi,\;r\cos\theta)`` | ``(r,\;r\sin\theta,\;1)`` |
| [`coorsys_spheroidal`](@ref) | ``(\varphi,p,q)`` | prolate spheroidal | see the constructor |

Their non-vanishing Christoffel symbols:

```math
\textbf{polar / cylindrical:}\qquad
\Gamma^r_{\theta\theta}=-r,
\qquad
\Gamma^\theta_{r\theta}=\Gamma^\theta_{\theta r}=\frac{1}{r}
```

```math
\textbf{spherical:}\qquad
\Gamma^r_{\theta\theta}=-r,\quad
\Gamma^r_{\varphi\varphi}=-r\sin^2\theta,\quad
\Gamma^\theta_{\varphi\varphi}=-\sin\theta\cos\theta,
```

```math
\Gamma^\theta_{r\theta}=\Gamma^\theta_{\theta r}=\frac{1}{r},\quad
\Gamma^\varphi_{r\varphi}=\Gamma^\varphi_{\varphi r}=\frac{1}{r},\quad
\Gamma^\varphi_{\theta\varphi}=\Gamma^\varphi_{\varphi\theta}=\cot\theta .
```

All Christoffel symbols of the Cartesian system vanish, which is the definition
of a Cartesian chart.

!!! warning "The spherical coordinates are ordered ``(\theta,\varphi,r)``"
    Not ``(r,\theta,\varphi)``. The ordering is chosen so that
    ``\theta=\varphi=0`` reproduces the canonical basis **in the canonical
    order**, which makes the spherical frame a genuine
    [`RotatedBasis`](@ref) — see [Rotations](@ref th-rotations) — and lets
    spherical results be compared with Euler-angle results without a
    permutation. The same ordering applies to
    [`init_spherical`](@ref) and to ``\underline{e}^\theta,\underline{e}^\varphi,
    \underline{e}^r``.

## Two implementations, one interface

The operators have the same names and meaning whether the derivatives are taken
symbolically or by automatic differentiation.

| | [`CoorSystemSym`](@ref) | [`CoorSystemNum`](@ref) |
| :--- | :--- | :--- |
| derivatives | SymPy / Symbolics, exact | `ForwardDiff`, machine precision |
| stores | ``\underline{OM}``, bases, ``\chi_i``, ``\Gamma`` as **expressions** | ``\chi``, ``R``, ``\Gamma`` as **closures**, evaluated pointwise |
| result | a formula, valid everywhere | a number, at one point |
| cost | grows with expression size | constant per point |
| needs | a symbolic ``\underline{OM}`` | any differentiable `OM(x)` |
| best for | deriving and simplifying closed forms | evaluating fields, sweeps, gradients |

`CoorSystemNum` stores its geometry as three closures — `χ_func` (Lamé),
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

Non-trivial symbolic charts — the spheroidal one is the standard example —
produce metric expressions that SymPy will not reduce on its own. `CoorSystemSym`
therefore accepts:

| Argument | Role |
| :------- | :--- |
| `tmp_coords` | auxiliary variables standing for compound expressions |
| `params` | constants appearing in ``\underline{OM}`` |
| `rules` | rewrite rules applied after each simplification |
| `tmp_var` | substitutions replacing expressions by temporary variables |
| `to_coords` | how to eliminate the temporaries before differentiating |

Without them the Christoffel symbols of a spheroidal chart are computable in
principle and unusable in practice. Their use is shown in
[Christoffel symbols](@ref tut-christoffel) among the tutorials.
