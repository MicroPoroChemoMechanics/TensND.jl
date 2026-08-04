# [Submanifolds](@id th-submanifolds)

A hypersurface embedded in ``\mathbb{R}^d`` carries two geometries: the one it
inherits as a metric space, and the one describing how it curves inside the
ambient space. `TensND` represents both with [`SubManifoldSym`](@ref)
(`src/submanifold.jl`), which is a [`CoorSystemSym`](@ref) with one coordinate
fewer than the ambient dimension and a unit normal completing the frame.

## Tangent frame and normal

Let ``\underline{OM}(q^1,\ldots,q^{d-1})`` parametrize the surface. The tangent
vectors, Lamé coefficients and unit tangents are those of
[Curvilinear differential calculus](@ref th-curvilinear),

```math
\underline{a}_i=\partial_i\underline{OM},
\qquad
\chi_i=\|\underline{a}_i\|,
\qquad
\underline{e}_i=\underline{a}_i/\chi_i
\qquad (i=1,\ldots,d-1),
```

and the frame is closed by the **unit normal** ``\underline{n}``, built as the
generalized cross product of the tangent vectors and stored last with
``\chi_d=1``. Reversing two coordinates reverses ``\underline{n}``, and with it
the sign of the second fundamental form below — the orientation is a choice, and
it must be stated whenever a curvature is reported.

## The two fundamental forms

```math
\boxed{\;a_{ij}=\underline{a}_i\cdot\underline{a}_j\;}
\qquad
\boxed{\;b_{ij}=\underline{n}\cdot\partial_i\underline{a}_j
              =-\,\partial_i\underline{n}\cdot\underline{a}_j\;}
```

- ``\boldsymbol{a}`` — the **first fundamental form**, the metric induced by the
  embedding. It measures lengths and angles *within* the surface and knows
  nothing of the ambient space. [`submetric`](@ref).
- ``\boldsymbol{b}`` — the **second fundamental form**, the normal component of
  the variation of the tangent frame. It measures how the surface bends *in* the
  ambient space. [`curvature`](@ref).

The second equality follows from ``\underline{n}\cdot\underline{a}_j=0`` by
differentiation, and is the Weingarten relation.

Both are stored as ``d``-dimensional order-2 tensors with a vanishing last row
and column, so they combine directly with tensors expressed in the full embedded
frame.

## Gauss and Weingarten

Decomposing ``\partial_i\underline{a}_j`` on the embedded frame gives the two
classical equations, tangential and normal:

```math
\partial_i\underline{a}_j=\Gamma^k_{ij}\,\underline{a}_k+b_{ij}\,\underline{n}
\qquad\text{(Gauss)}
```

```math
\partial_i\underline{n}=-\,b_i{}^{j}\,\underline{a}_j
\qquad\text{(Weingarten)}
```

`SubManifoldSym` stores exactly this: the connection array ``\Gamma`` holds the
intrinsic Christoffel symbols in its tangential block, ``\boldsymbol{b}`` in the
slice ``\Gamma[:,:,d]`` (Gauss), and ``-b_i{}^{j}`` in ``\Gamma[:,d,:]``
(Weingarten).

!!! warning "`Riemann` returns Christoffel symbols, not a curvature"
    [`Riemann`](@ref) returns the purely tangential block
    ``\Gamma[1{:}d{-}1,1{:}d{-}1,1{:}d{-}1]`` — the **Christoffel symbols of the
    induced metric**. These are connection coefficients: they are not tensor
    components, they change under a change of chart by an inhomogeneous rule, and
    they can be made to vanish at any single point. A Riemann curvature tensor
    can do none of those things.

    For a sphere of radius ``R`` parametrized by ``(\theta,\varphi)`` the
    function returns

    ```math
    \Gamma^\theta_{\varphi\varphi}=-\sin\theta\cos\theta,
    \qquad
    \Gamma^\varphi_{\theta\varphi}=\Gamma^\varphi_{\varphi\theta}=\cot\theta ,
    ```

    which are indeed the sphere's connection coefficients. The name is
    misleading and is kept only for backward compatibility. The intrinsic
    curvature is recovered from ``\boldsymbol{a}`` and ``\boldsymbol{b}`` by the
    Gauss equation below.

## Curvatures

With the **shape operator** ``b_i{}^{j}=a^{ik}b_{kj}`` — the second form with one
index raised —

| Quantity | Definition | Sphere of radius ``R``, outward normal |
| :------- | :--------- | :------------------------------------- |
| principal curvatures | eigenvalues ``\kappa_i`` of ``b_i{}^{j}`` | ``-1/R,\ -1/R`` |
| mean curvature | ``H=\tfrac{1}{d-1}\,\mathrm{tr}\,b_i{}^{j}`` | ``-1/R`` |
| Gaussian curvature | ``K=\det b_i{}^{j}`` | ``1/R^2`` |

Gauss's *Theorema Egregium* is that ``K`` depends only on ``\boldsymbol{a}``,
even though its definition uses ``\boldsymbol{b}``: it is an intrinsic quantity.
The cylinder is the standard illustration — one principal curvature vanishes, so
``K=0`` and the cylinder is intrinsically flat, which is why it can be unrolled
onto a plane without distortion while a sphere cannot.

## Worked closed forms

For a sphere of radius ``R`` parametrized by ``(\theta,\varphi)``, with the
outward normal ``\underline{n}=\underline{e}^r``:

```math
\boldsymbol{a}=\mathrm{diag}\bigl(R^2,\;R^2\sin^2\theta,\;0\bigr),
\qquad
\boldsymbol{b}=-\frac{\boldsymbol{a}}{R} .
```

The proportionality ``\boldsymbol{b}=-\boldsymbol{a}/R`` is the defining property
of an umbilical surface: every direction curves identically.

| Surface | ``K`` | ``H`` | Comment |
| :------ | :---- | :---- | :------ |
| plane | ``0`` | ``0`` | ``\boldsymbol{b}=0`` |
| sphere, radius ``R`` | ``1/R^2`` | ``-1/R`` | umbilical |
| cylinder, radius ``R`` | ``0`` | ``-1/(2R)`` | intrinsically flat |

## Differential operators on the surface

Because `SubManifoldSym` is an `AbstractCoorSystem`, the operators
[`GRAD`](@ref), [`DIV`](@ref), [`LAPLACE`](@ref) apply to fields defined on the
surface, with the intrinsic connection. The relation that ties the two
geometries together is the Weingarten equation read as a gradient:

```math
\boxed{\;\mathrm{GRAD}(\underline{n})=-\,\boldsymbol{b}\;}
```

It is the cheapest check that an implementation's normal orientation and
curvature sign agree, and it is pinned by a test
([Testing and conventions](@ref dev-testing)).

!!! warning "Compare tensors, not stored arrays"
    This is an identity between **tensors**, and the two sides are not stored
    the same way: ``\mathrm{GRAD}(\underline{n})`` comes out on one basis and
    variance, ``\boldsymbol{b}`` on another. Subtracting their
    [`get_array`](@ref)s therefore gives a nonzero, and rather convincing,
    wrong answer — on a sphere of radius ``R`` the raw arrays differ by a factor
    ``R^2``, so the discrepancy even vanishes for ``R=1``.

    Bring both to a common basis and variance first —
    [`components`](@ref)`(t, ℬ, var)` or [`change_tens`](@ref) — after which the
    difference is exactly zero in *every* variance. This is the practical face
    of [Bases and variance](@ref th-bases-variance): a tensor is not its
    component array, and only components in the same basis may be compared.

!!! note "Declare the parameter ranges, or fight `Abs`"
    Declaring ``\theta`` merely `real` leaves SymPy unable to reduce
    ``|\sin\theta|`` to ``\sin\theta``, and the curvature comes back multiplied
    by ``\sin\theta/|\sin\theta|``. The results are correct but unreadable. Pass
    the `rules` of [`SubManifoldSym`](@ref) — or restrict the symbol's
    assumptions — so that the sign is resolved. This is the same simplification
    machinery described on
    [Curvilinear differential calculus](@ref th-curvilinear).
