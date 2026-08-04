```@meta
EditURL = "../../../../scripts/13_submanifolds.jl"
```

# Surfaces: fundamental forms and curvature

A [`SubManifoldSym`](@ref) is a hypersurface embedded in ``\mathbb{R}^d``: one
coordinate fewer than the ambient space, with the unit normal completing the
frame. It carries both geometries — the metric it inherits, and the way it
bends inside the ambient space.

```math
a_{ij}=\underline{a}_i\cdot\underline{a}_j
\qquad\text{(first fundamental form)}
```

```math
b_{ij}=\underline{n}\cdot\partial_i\underline{a}_j
\qquad\text{(second fundamental form)}
```

Theory: [Submanifolds](@ref th-submanifolds).

````@example submanifolds
using TensND
using LinearAlgebra
using SymPy
````

## Resolving ``|\sin\theta|``

SymPy cannot know that ``\theta`` lies in ``(0,\pi)``: declaring the symbol
`real`, or even `positive`, still leaves ``|\sin\theta|`` irreducible, and
every curvature comes back multiplied by ``\sin\theta/|\sin\theta|`` —
correct, but unreadable. The `rules` argument of [`SubManifoldSym`](@ref) is
exactly the mechanism for this: a rewrite applied after each simplification.

````@example submanifolds
R = symbols("R", positive = true)
θ = symbols("θ", positive = true)
ϕ = symbols("ϕ", real = true)
z = symbols("z", real = true)

polar_rules = Dict(abs(sin(θ)) => sin(θ))
````

## The sphere

The classical umbilical surface: every direction curves identically.

````@example submanifolds
Sphere = SubManifoldSym(
    Tens(R * [sin(θ)cos(ϕ), sin(θ)sin(ϕ), cos(θ)]), (θ, ϕ), (), (R,);
    rules = polar_rules,
)
````

The unit normal is the radial vector:

````@example submanifolds
tsimplify.(components_canon(normal(Sphere)))
````

First fundamental form — the familiar ``\mathrm{d}s^2=R^2\mathrm{d}\theta^2
+R^2\sin^2\theta\,\mathrm{d}\varphi^2``:

````@example submanifolds
𝐚 = submetric(Sphere)
get_array(𝐚)
````

Second fundamental form:

````@example submanifolds
𝐛 = curvature(Sphere)
tsimplify.(get_array(𝐛))
````

For a sphere with the **outward** normal, ``\boldsymbol{b}=-\boldsymbol{a}/R``:

````@example submanifolds
tsimplify.(get_array(𝐛) + get_array(𝐚) / R)
````

## Curvatures

The shape operator is the second form with one index raised,
``b_i{}^{j}=a^{ik}b_{kj}``; its determinant is the Gaussian curvature and half
its trace the mean curvature.

````@example submanifolds
function curvatures(SM)
    dim = length(get_array(submetric(SM))[1, :])
    tang = 1:(dim - 1)
    a = tsimplify.(get_array(submetric(SM))[tang, tang])
    b = tsimplify.(get_array(curvature(SM))[tang, tang])
    S = tsimplify.(inv(a) * b)
    return (K = tsimplify(det(S)), H = tsimplify(tr(S) / (dim - 1)), shape = tsimplify.(diag(S)))
end

cs = curvatures(Sphere)
println("sphere : K = ", cs.K, "   H = ", cs.H, "   principal curvatures = ", cs.shape)
````

``K=1/R^2`` and ``H=-1/R``: both principal curvatures equal ``-1/R``, the sign
following the outward orientation of ``\underline{n}``.

## The cylinder

One principal curvature vanishes, so ``K=0``: the cylinder is **intrinsically
flat** and can be unrolled onto a plane without distortion, although it is
visibly curved in space. This is Gauss's *Theorema Egregium* in one line.

````@example submanifolds
Cylinder = SubManifoldSym(Tens([R * cos(ϕ), R * sin(ϕ), z]), (ϕ, z), (), (R,))

tsimplify.(get_array(submetric(Cylinder))), tsimplify.(get_array(curvature(Cylinder)))
````

````@example submanifolds
cc = curvatures(Cylinder)
println("cylinder : K = ", cc.K, "   H = ", cc.H, "   principal curvatures = ", cc.shape)
````

## The plane

The degenerate case: the normal is constant, so ``\boldsymbol{b}=0`` and both
curvatures vanish.

````@example submanifolds
x, y = symbols("x y", real = true)
Plane = SubManifoldSym(Tens([x, y, zero(x)]), (x, y))

tsimplify.(get_array(curvature(Plane)))
````

````@example submanifolds
cp = curvatures(Plane)
println("plane : K = ", cp.K, "   H = ", cp.H)
````

## A paraboloid

A surface whose curvature actually varies with position.

````@example submanifolds
Paraboloid = SubManifoldSym(Tens([x, y, (x^2 + y^2) / (2R)]), (x, y), (), (R,))
pb = curvatures(Paraboloid)
````

The closed forms are large, so they are more usefully read as values. At the
apex ``x=y=0`` the surface osculates a sphere of radius ``R``, giving
``|K|=1/R^2`` and ``|H|=1/R``; both fall off as one moves outwards.

Note the **sign**: the ``(x,y)`` parametrization orients the normal towards
the concave side, opposite to the outward normal used for the sphere above, so
``H`` comes out positive here where the sphere gave ``-1/R``. The Gaussian
curvature ``K``, being a determinant, is insensitive to the orientation.

````@example submanifolds
println("  s = √(x²+y²)      K·R²           H·R")
for s in (0.0, 0.5, 1.0, 2.0, 4.0)
    sub = Dict(x => Sym(s), y => Sym(0), R => Sym(1))
    Kv = Float64(simplify(pb.K.subs(sub)))
    Hv = Float64(simplify(pb.H.subs(sub)))
    println("     ", rpad(s, 12), rpad(round(Kv, digits = 8), 15), round(Hv, digits = 8))
end
````

The apex values, symbolically:

````@example submanifolds
(simplify(pb.K.subs(Dict(x => Sym(0), y => Sym(0)))), simplify(pb.H.subs(Dict(x => Sym(0), y => Sym(0)))))
````

## The Weingarten relation as a gradient

Because a `SubManifoldSym` is a coordinate system, the differential operators
apply to fields defined on the surface. The identity tying the two geometries
together is

```math
\mathrm{GRAD}(\underline{n})=-\,\boldsymbol{b}.
```

````@example submanifolds
G = GRAD(normal(Sphere), Sphere)
````

!!! warning "Compare tensors, not stored arrays"
    The two sides are stored on different bases and variances, so subtracting
    their `get_array`s gives a nonzero — and rather convincing — wrong answer.
    On a sphere of radius ``R`` the raw arrays differ by a factor ``R^2``, so
    the discrepancy even vanishes for ``R=1``.

````@example submanifolds
tsimplify.(get_array(G)), tsimplify.(get_array(𝐛))
````

Brought to a **common** basis and variance the difference is exactly zero, in
every variance:

````@example submanifolds
ℬnat = natural_basis(Sphere)
for var in ((:cov, :cov), (:cont, :cov), (:cov, :cont), (:cont, :cont))
    resid = tsimplify.(components(G, ℬnat, var) + components(𝐛, ℬnat, var))
    println("  variance ", var, " : all zero ? ", all(iszero, resid))
end
````

This is the practical face of [Bases and variance](@ref th-bases-variance): a
tensor is not its component array, and only components in the same basis may
be compared.

## The connection accessor is misnamed

[`Riemann`](@ref) returns the **Christoffel symbols** of the induced metric,
not the Riemann curvature tensor. For the sphere:

````@example submanifolds
Γ = Riemann(Sphere)
println("Γ^θ_ϕϕ = ", tsimplify(Γ[2, 2, 1]))
println("Γ^ϕ_θϕ = ", tsimplify(Γ[1, 2, 2]), "     Γ^ϕ_ϕθ = ", tsimplify(Γ[2, 1, 2]))
````

These are the sphere's connection coefficients. They are not tensor components
— they change inhomogeneously under a change of chart and can be made to
vanish at any single point, which no curvature can. The name is kept for
backward compatibility only; the intrinsic curvature is the ``K`` computed
above from ``\boldsymbol{a}`` and ``\boldsymbol{b}``.

## Summary

| Surface | ``K`` | ``H`` | comment |
|:--|:--|:--|:--|
| plane | ``0`` | ``0`` | ``\boldsymbol{b}=0`` |
| sphere, radius ``R`` | ``1/R^2`` | ``-1/R`` | umbilical |
| cylinder, radius ``R`` | ``0`` | ``-1/(2R)`` | intrinsically flat |
| paraboloid | varies | varies | ``\to1/R^2`` at the apex |

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

