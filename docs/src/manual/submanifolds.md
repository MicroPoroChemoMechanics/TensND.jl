# [Submanifolds](@id man-submanifolds)

Hypersurfaces embedded in ``\mathbb{R}^d``: their induced metric, their
curvature, and differential operators on them. The geometry is on
[Submanifolds](@ref th-submanifolds).

## Construction

```julia
SubManifoldSym(OM, coords, tmp_coords = (), params = ();
               rules = Dict(), tmp_var = Dict(), to_coords = Dict())
```

`coords` holds **one coordinate fewer** than the ambient dimension: a surface in
``\mathbb{R}^3`` is parametrized by two. The unit normal is built automatically
and completes the frame, so the stored bases remain ``d``-dimensional.

```@example sm
using TensND, SymPy

R = symbols("R", positive = true)
θ = symbols("θ", positive = true)
ϕ = symbols("ϕ", real = true)

Sphere = SubManifoldSym(
    Tens(R * [sin(θ)cos(ϕ), sin(θ)sin(ϕ), cos(θ)]), (θ, ϕ), (), (R,);
    rules = Dict(abs(sin(θ)) => sin(θ)),
)
tsimplify.(components_canon(normal(Sphere)))
```

!!! note "Use `rules` to resolve sign ambiguities"
    SymPy cannot know that ``\theta\in(0,\pi)``, so ``|\sin\theta|`` survives
    every simplification and multiplies the curvature by
    ``\sin\theta/|\sin\theta|``. The `rules` argument rewrites it away, as
    above. The same machinery is described under
    [Coordinate systems](@ref man-coorsystems).

## The four accessors

| Function | Returns |
| :--- | :--- |
| [`normal`](@ref) | the unit normal ``\underline{n}`` |
| [`submetric`](@ref) | the first fundamental form ``\boldsymbol{a}`` |
| [`curvature`](@ref) | the second fundamental form ``\boldsymbol{b}`` |
| [`connection`](@ref) | the intrinsic connection coefficients ``\Gamma^\gamma_{\alpha\beta}`` |

```@example sm
get_array(submetric(Sphere))
```

```@example sm
tsimplify.(get_array(curvature(Sphere)))
```

Both are returned as ``d``-dimensional order-2 tensors with a vanishing last row
and column, so they combine directly with tensors in the embedded frame.

For a sphere with the outward normal, ``\boldsymbol{b}=-\boldsymbol{a}/R``:

```@example sm
tsimplify.(get_array(curvature(Sphere)) + get_array(submetric(Sphere)) / R)
```

!!! warning "The orientation decides the sign"
    The normal is the generalized cross product of the tangent vectors **in the
    order given**, so swapping two coordinates reverses it — and with it the
    sign of [`curvature`](@ref) and of the mean curvature. State the orientation
    whenever a curvature is reported. The Gaussian curvature, being a
    determinant, is unaffected.

!!! note "`connection` — renamed from `Riemann`"
    [`connection`](@ref) returns the **connection coefficients** (Christoffel
    symbols) of the induced metric, never a Riemann curvature tensor. The old
    name `Riemann` named the wrong object; it still works and forwards, with a
    deprecation warning. For the sphere it
    gives ``\Gamma^\theta_{\varphi\varphi}=-\sin\theta\cos\theta`` and
    ``\Gamma^\varphi_{\theta\varphi}=\cot\theta``. These are connection
    coefficients: not tensor components, and removable at any single point.
    Compute the intrinsic curvature from [`submetric`](@ref) and
    [`curvature`](@ref) instead.

## Curvatures

The shape operator is the second form with one index raised,
``b_\alpha{}^{\beta}=a^{\beta\gamma}b_{\gamma\alpha}``; its determinant is the Gaussian curvature and its
normalized trace the mean curvature.

```@example sm
using LinearAlgebra
a2 = tsimplify.(get_array(submetric(Sphere))[1:2, 1:2])
b2 = tsimplify.(get_array(curvature(Sphere))[1:2, 1:2])
S = tsimplify.(inv(a2) * b2)
(tsimplify(det(S)), tsimplify(tr(S) / 2))
```

| Surface | ``K`` | ``H`` |
| :--- | :--- | :--- |
| plane | ``0`` | ``0`` |
| sphere, radius ``R``, outward normal | ``1/R^2`` | ``-1/R`` |
| cylinder, radius ``R`` | ``0`` | ``-1/(2R)`` |

## Differential operators on the surface

A `SubManifoldSym` is an `AbstractCoorSystem`, so [`GRAD`](@ref),
[`DIV`](@ref) and [`LAPLACE`](@ref) apply to fields defined on it, with the
intrinsic connection. The identity linking the two geometries is the Weingarten
relation read as a gradient,

```math
\mathrm{GRAD}(\underline{n})=-\,\boldsymbol{b}.
```

!!! warning "Compare tensors, not stored arrays"
    The two sides are stored on different bases and variances, so subtracting
    their [`get_array`](@ref)s gives a nonzero — and plausible-looking — wrong
    answer. Bring both to a common basis and variance first:

    ```julia
    components(GRAD(normal(SM), SM), natural_basis(SM), (:cov, :cov)) +
    components(curvature(SM),       natural_basis(SM), (:cov, :cov))
    ```

    which is exactly zero. This is the practical face of
    [Bases and variance](@ref th-bases-variance).

A worked version of all of the above, on the sphere, cylinder, plane and
paraboloid, is in
[Surfaces: fundamental forms and curvature](@ref).
