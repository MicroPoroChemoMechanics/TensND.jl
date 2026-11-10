```@meta
EditURL = "../../../../scripts/16_spheroidal_appendix.jl"
```

# [The spheroid in spheroidal coordinates](@id tut-spheroidal-appendix)

!!! tip "See also"
    **Theory**: [Spheroidal coordinates](@ref th-spheroidal) · [Ellipsoidal coordinates](@ref th-ellipsoidal) · [Curvilinear differential calculus](@ref th-curvilinear) · [Spheroidal harmonics and their re-expansion](@ref th-spheroidal-harmonics) · [Submanifolds](@ref th-submanifolds)

    **Manual**: [Coordinate systems (symbolic)](@ref man-coorsystems) · [Charts and coordinate conversions](@ref man-charts) · [Spheroidal harmonics](@ref man-spheroidal-harmonics) · [Submanifolds](@ref man-submanifolds)

    **Tutorials**: [Ellipsoidal coordinates](@ref tut-ellipsoidal) · [Changing the focal distance of spheroidal harmonics](@ref tut-spheroidal-reexpansion) · [Surfaces: fundamental forms and curvature](@ref tut-submanifolds)

The spheroidal coordinates are the system in which a spheroidal inclusion is a
coordinate surface. This tutorial builds the prolate system step by step and
recovers, with `TensND`, the results that a problem set on a spheroid needs:
the Lamé coefficients and the unit vectors, the volume and the surface of a
spheroid, the gradient, the divergence and the Laplacian, then the spheroid
seen as a surface, with its curvatures and its surface operators. They are
the formulas of the appendix on spheroidal coordinates of
[barthelemy2020](@citet), derived on
[Spheroidal coordinates](@ref th-spheroidal). The oblate system follows by a
substitution, checked at the end.

````@example spheroidal_appendix
using TensND
using LinearAlgebra
using SymPyPythonCall
using Plots
````

## The chart

With ``c`` the semi-focal distance, ``-1\le p\le1`` and ``q\ge1``, the point
of coordinates ``(\varphi,p,q)`` is

```math
x_1=c\,\sqrt{1-p^2}\,\sqrt{q^2-1}\,\cos\varphi,\qquad
x_2=c\,\sqrt{1-p^2}\,\sqrt{q^2-1}\,\sin\varphi,\qquad
x_3=c\,p\,q .
```

[`coorsys_spheroidal`](@ref) returns the system. Its position vector is
written with the two auxiliary symbols ``\bar{p}=\sqrt{1-p^2}`` and
``\bar{q}=\sqrt{q^2-1}``, which keep the expressions short:

````@example spheroidal_appendix
Spheroidal = coorsys_spheroidal()
ϕ, p, q = get_coords(Spheroidal)
c = symbols("c", positive = true)
position_vector(Spheroidal)
````

The same map between symbols, without auxiliary symbols, is the chart:

````@example spheroidal_appendix
chart = ProlateSpheroidalChart(c)
x₁, x₂, x₃ = coords_to_cartesian(chart, (ϕ, p, q))
````

The surface ``q=\text{const}`` is a spheroid of semi-axes ``c\,q`` along
``\underline{e}_3`` and ``c\sqrt{q^2-1}`` across:

````@example spheroidal_appendix
simplify((x₁^2 + x₂^2) / (c^2 * (q^2 - 1)) + x₃^2 / (c^2 * q^2))
````

The position vector is also a combination of two spheroidal harmonics of
degree one,

```math
\underline{x}=c\,\bigl(-P_1^1(p)\,P_1^1(q)\,\underline{u}_\varphi+P_1(p)\,P_1(q)\,\underline{e}_3\bigr),
\qquad
\underline{u}_\varphi=\cos\varphi\,\underline{e}_1+\sin\varphi\,\underline{e}_2,
```

with the Legendre functions of [`legendre_P`](@ref) on ``[-1,1]`` and of
[`radial_P`](@ref) beyond. Both return the degrees ``0`` to ``N``, the degree
``n`` at the index ``n+1``:

````@example spheroidal_appendix
P₁¹p = legendre_P(1, 1, p)[2]
P₁¹q = radial_P(chart, 1, 1, q)[2]
P₁¹p, P₁¹q
````

````@example spheroidal_appendix
simplify(x₁ + c * P₁¹p * P₁¹q * cos(ϕ)), simplify(x₃ - c * legendre_P(1, 0, p)[2] * radial_P(chart, 1, 0, q)[2])
````

## Lamé coefficients and unit vectors

The natural vectors are the derivatives of the position vector with respect
to the coordinates; their norms are the Lamé coefficients

```math
\chi_\varphi=c\,\sqrt{1-p^2}\,\sqrt{q^2-1},\qquad
\chi_p=c\,\sqrt{\frac{q^2-p^2}{1-p^2}},\qquad
\chi_q=c\,\sqrt{\frac{q^2-p^2}{q^2-1}} .
```

````@example spheroidal_appendix
χᵩ, χₚ, χ_q = lame_coefficients(Spheroidal)
````

The unit vectors ``\underline{e}^\varphi,\underline{e}^p,\underline{e}^q``
are the columns of the normalized basis, in Cartesian components:

````@example spheroidal_appendix
ℬ = normalized_basis(Spheroidal)
𝐞ᵠ, 𝐞ᵖ, 𝐞q = unitvec(Spheroidal)
Matrix(ℬ)
````

They are orthonormal, and each natural vector is its unit vector times its
Lamé coefficient. The check differentiates the chart directly:

````@example spheroidal_appendix
E = Matrix(ℬ)
J = [diff(xᵢ, ξ) for xᵢ in (x₁, x₂, x₃), ξ in (ϕ, p, q)]
simplify.(transpose(E) * E), simplify.(J - E * Diagonal([χᵩ, χₚ, χ_q]))
````

## Volume and surface of a spheroid

The volume element is the product of the three Lamé coefficients and the
surface element of a spheroid ``q=\text{const}`` the product of the two
tangent ones:

```math
\mathrm{d}\Omega=c^3\,(q^2-p^2)\,\mathrm{d}\varphi\,\mathrm{d}p\,\mathrm{d}q,
\qquad
\mathrm{d}S_q=c^2\,\sqrt{q^2-p^2}\,\sqrt{q^2-1}\,\mathrm{d}\varphi\,\mathrm{d}p .
```

````@example spheroidal_appendix
dΩ = simplify(χᵩ * χₚ * χ_q)
dS = simplify(χᵩ * χₚ)
dΩ, dS
````

Integrating over ``\varphi\in[0,2\pi]`` and ``p\in[-1,1]``, and for the volume
over the radial coordinate from the focal segment ``q=1`` to the spheroid:

````@example spheroidal_appendix
t = symbols("t", positive = true)
V = factor(integrate(integrate(integrate(dΩ(q => t), (ϕ, 0, 2PI)), (p, -1, 1)), (t, 1, q)))
````

````@example spheroidal_appendix
S = simplify(integrate(integrate(dS, (ϕ, 0, 2PI)), (p, -1, 1)))
````

which are ``V=\tfrac43\pi\,c^3\,q\,(q^2-1)`` and
``S=2\pi\,c^2\sqrt{q^2-1}\,\bigl(\sqrt{q^2-1}+q^2\arcsin(1/q)\bigr)``. In terms
of the equatorial radius ``b=c\sqrt{q^2-1}`` and of the aspect ratio
``\omega=q/\sqrt{q^2-1}>1``, the volume is ``\tfrac43\pi\,b^3\,\omega``:

````@example spheroidal_appendix
b, ω = symbols("b ω", positive = true)
simplify(V(c => b / sqrt(q^2 - 1)) / (Sym(4) / 3 * PI * b^3) - q / sqrt(q^2 - 1))
````

The specific surface ``\sigma_v=S/V`` measures how much interface a given
volume of inclusions carries. Compared with that of the sphere of the same
volume, ``3/R`` with ``R=b\,\omega^{1/3}``, it only depends on the aspect
ratio. The oblate branch ``\omega<1`` uses the oblate formulas obtained at the
end of this page.

````@example spheroidal_appendix
σ_prolate(ω) = (q = ω / sqrt(ω^2 - 1); ω^(-2 / 3) / 2 * (1 + q^2 / sqrt(q^2 - 1) * asin(1 / q)))
σ_oblate(ω) = (τ = ω / sqrt(1 - ω^2); ω^(-2 / 3) / 2 * (1 + τ^2 / sqrt(τ^2 + 1) * asinh(1 / τ)))
ωs_oblate = 10 .^ range(-2, -0.01; length = 100)
ωs_prolate = 10 .^ range(0.01, 2; length = 100)
plot(
    ωs_oblate, σ_oblate.(ωs_oblate); label = "oblate", xscale = :log10, lw = 2,
    xlabel = "aspect ratio ω", ylabel = "specific surface / that of the sphere of the same volume",
    legend = :top, size = (700, 420), guidefontsize = 9,
)
plot!(ωs_prolate, σ_prolate.(ωs_prolate); label = "prolate", lw = 2)
scatter!([1.0], [1.0]; label = "sphere", color = :black)
````

A needle (``\omega\gg1``) and a disk (``\omega\ll1``) both expose much more
surface than the sphere of the same volume.

## Gradient, divergence and Laplacian

The operators apply to any field of the coordinates. For a scalar field
``f(\varphi,p,q)``, the physical components of the gradient are the partial
derivatives divided by the Lamé coefficients:

````@example spheroidal_appendix
f = SymFunction("f", real = true)(ϕ, p, q)
simplify.(components(grad(f, Spheroidal), ℬ))
````

The Laplacian is compared with its expression in an orthogonal system,

```math
\Delta f=\frac{1}{c^2\,(q^2-p^2)}\left[
\frac{\partial}{\partial p}\Bigl((1-p^2)\,\frac{\partial f}{\partial p}\Bigr)
+\frac{\partial}{\partial q}\Bigl((q^2-1)\,\frac{\partial f}{\partial q}\Bigr)\right]
+\frac{1}{c^2\,(1-p^2)\,(q^2-1)}\,\frac{\partial^2f}{\partial\varphi^2} :
```

````@example spheroidal_appendix
Δf = (diff((1 - p^2) * diff(f, p), p) + diff((q^2 - 1) * diff(f, q), q)) / (c^2 * (q^2 - p^2)) +
    diff(f, ϕ, 2) / (c^2 * (1 - p^2) * (q^2 - 1))
simplify(laplacian(f, Spheroidal) - Δf)
````

For a vector field of physical components ``(v_\varphi,v_p,v_q)``,

```math
\mathrm{div}\,\underline{v}=\frac{1}{c\,\bar{p}\,\bar{q}}\,\frac{\partial v_\varphi}{\partial\varphi}
+\frac{1}{c\,(q^2-p^2)}\left[
\frac{\partial}{\partial p}\Bigl(\bar{p}\,\sqrt{q^2-p^2}\;v_p\Bigr)
+\frac{\partial}{\partial q}\Bigl(\bar{q}\,\sqrt{q^2-p^2}\;v_q\Bigr)\right] :
```

````@example spheroidal_appendix
vᵩ, vₚ, v_q = (SymFunction("v$i", real = true)(ϕ, p, q) for i in ("ϕ", "p", "q"))
𝐯 = vᵩ * 𝐞ᵠ + vₚ * 𝐞ᵖ + v_q * 𝐞q
D = q^2 - p^2
divv = diff(vᵩ, ϕ) / (c * sqrt(1 - p^2) * sqrt(q^2 - 1)) +
    (diff(sqrt(1 - p^2) * sqrt(D) * vₚ, p) + diff(sqrt(q^2 - 1) * sqrt(D) * v_q, q)) / (c * D)
simplify(divergence(𝐯, Spheroidal) - divv)
````

The Laplace equation separates in these coordinates, and its solutions are
products of Legendre functions of ``p`` and of ``q``. The exterior harmonic of
degree one, ``P_1(p)\,Q_1(q)`` with ``Q_1(q)=q\,\mathrm{arccoth}\,q-1``, is the
potential outside a spheroid placed in a uniform field along its axis:

````@example spheroidal_appendix
simplify(laplacian(p * (q * acoth(q) - 1), Spheroidal))
````

## The spheroid as a coordinate surface

On the spheroid ``q=\text{const}``, the unit normal is ``\underline{e}^q`` and
``(\varphi,p)`` are coordinates of the surface. The gradient of the normal
field, read in the normalized basis, is minus the curvature tensor of the
surface (Weingarten): its diagonal holds the principal curvatures with their
sign changed.

````@example spheroidal_appendix
𝐧 = 𝐞q
∇𝐧 = simplify.(components(grad(𝐧, Spheroidal), ℬ))
````

The principal curvatures along ``\underline{e}^\varphi`` and
``\underline{e}^p`` are therefore

```math
\kappa_\varphi=-\frac{q}{c\,\sqrt{q^2-1}\,\sqrt{q^2-p^2}},
\qquad
\kappa_p=-\frac{q\,\sqrt{q^2-1}}{c\,(q^2-p^2)^{3/2}},
```

negative because the normal points outward:

````@example spheroidal_appendix
κᵩ, κₚ = -∇𝐧[1, 1], -∇𝐧[2, 2]
````

At the equator they are the inverse of the equatorial radius and the curvature
of the meridian ellipse at the end of its minor axis, ``b/a^2`` with
``a=c\,q`` and ``b=c\sqrt{q^2-1}``:

````@example spheroidal_appendix
simplify(κᵩ(p => 0) + 1 / (c * sqrt(q^2 - 1))), simplify(κₚ(p => 0) + c * sqrt(q^2 - 1) / (c * q)^2)
````

The divergence of the normal is minus their sum, twice the mean curvature up
to the sign:

````@example spheroidal_appendix
simplify(divergence(𝐧, Spheroidal)), simplify(divergence(𝐧, Spheroidal) + κᵩ + κₚ)
````

## Surface operators

A field defined on the spheroid depends on ``(\varphi,p)`` only. Its surface
gradient is its gradient, which has no component along the normal:

````@example spheroidal_appendix
g = SymFunction("g", real = true)(ϕ, p)
∇g = grad(g, Spheroidal)
simplify.(components(∇g, ℬ))
````

The surface divergence of a field ``\underline{w}`` removes from the
divergence the stretching along the normal,
``\mathrm{div}_S\,\underline{w}=\mathrm{div}\,\underline{w}-\underline{n}\cdot\nabla\underline{w}\cdot\underline{n}``:

````@example spheroidal_appendix
div_S(𝐰) = divergence(𝐰, Spheroidal) - 𝐧 ⋅ (grad(𝐰, Spheroidal) ⋅ 𝐧)
````

For a tangent field ``\underline{w}=w_\varphi\,\underline{e}^\varphi+w_p\,\underline{e}^p``
it is

```math
\mathrm{div}_S\,\underline{w}=\frac{1}{\chi_\varphi\,\chi_p}\left[
\frac{\partial}{\partial\varphi}\bigl(\chi_p\,w_\varphi\bigr)
+\frac{\partial}{\partial p}\bigl(\chi_\varphi\,w_p\bigr)\right] :
```

````@example spheroidal_appendix
wᵩ, wₚ = (SymFunction("w$i", real = true)(ϕ, p) for i in ("ϕ", "p"))
𝐰 = wᵩ * 𝐞ᵠ + wₚ * 𝐞ᵖ
simplify(div_S(𝐰) - (diff(χₚ * wᵩ, ϕ) + diff(χᵩ * wₚ, p)) / (χᵩ * χₚ))
````

The surface Laplacian is the surface divergence of the surface gradient. It
is the operator that the conditions on an imperfect interface involve:

```math
\Delta_Sg=\frac{1}{c^2}\left[
\frac{1}{(1-p^2)\,(q^2-1)}\,\frac{\partial^2g}{\partial\varphi^2}
+\frac{1}{\sqrt{q^2-p^2}}\,\frac{\partial}{\partial p}
\Bigl(\frac{1-p^2}{\sqrt{q^2-p^2}}\,\frac{\partial g}{\partial p}\Bigr)\right] .
```

````@example spheroidal_appendix
Δ_Sg = (diff(g, ϕ, 2) / ((1 - p^2) * (q^2 - 1)) + diff((1 - p^2) / sqrt(D) * diff(g, p), p) / sqrt(D)) / c^2
simplify(div_S(∇g) - Δ_Sg)
````

## The same spheroid as a submanifold

A [`Submanifold`](@ref) is built from the position vector of a surface alone,
without the system around it. The spheroid ``c=1``, ``q=2`` is parameterized
by the colatitude-like angle ``\vartheta`` such that ``p=\cos\vartheta``, and
by ``\varphi``; in this order the normal points outward.

````@example spheroidal_appendix
ϑ = symbols("ϑ", real = true)
q₀ = Sym(2)
Spheroid = Submanifold(
    Tens([sqrt(q₀^2 - 1) * sin(ϑ) * cos(ϕ), sqrt(q₀^2 - 1) * sin(ϑ) * sin(ϕ), q₀ * cos(ϑ)]),
    (ϑ, ϕ); rules = Dict(abs(sin(ϑ)) => sin(ϑ)),
)
simplify.(canonical_components(normal(Spheroid)))
````

Its first and second fundamental forms:

````@example spheroidal_appendix
𝐚 = simplify.(get_array(submetric(Spheroid))[1:2, 1:2])
𝐛 = simplify.(get_array(curvature(Spheroid))[1:2, 1:2])
𝐚, 𝐛
````

The eigenvalues of ``\boldsymbol{a}^{-1}\boldsymbol{b}`` are the principal
curvatures, here along ``\vartheta`` then ``\varphi``. They are those found
above, at ``c=1``, ``q=2`` and ``p=\cos\vartheta``:

````@example spheroidal_appendix
at(e) = e(c => 1, q => q₀, p => cos(ϑ))
simplify.(diag(inv(𝐚) * 𝐛) - [at(κₚ), at(κᵩ)])
````

The operators of a submanifold are the surface operators. The Laplacian of a
field ``h(\vartheta,\varphi)`` is compared with ``\Delta_S`` above, in which
``\partial/\partial p=-(1/\sin\vartheta)\,\partial/\partial\vartheta``:

````@example spheroidal_appendix
h = SymFunction("h", real = true)(ϑ, ϕ)
D₀ = q₀^2 - cos(ϑ)^2
Δ_Sh = diff(h, ϕ, 2) / (sin(ϑ)^2 * (q₀^2 - 1)) + diff(sin(ϑ) / sqrt(D₀) * diff(h, ϑ), ϑ) / (sqrt(D₀) * sin(ϑ))
simplify(laplacian(h, Spheroid) - Δ_Sh)
````

## The oblate system

The oblate system ``(\varphi,p,\tau)`` revolves the same confocal conics
about the other axis. It is the prolate one in which ``c`` is replaced by
``-\imath\,c`` and ``q`` by ``\imath\,\tau``:

````@example spheroidal_appendix
Oblate = coorsys_spheroidal_oblate()
τ = get_coords(Oblate)[3]
lame_coefficients(Oblate)
````

The substitution is checked on the squares of the Lamé coefficients, which
involve no choice of branch. Since ``\mathrm{d}q=\imath\,\mathrm{d}\tau``, the
line element ``\chi_q^2\,\mathrm{d}q^2`` becomes
``-\chi_q^2\,\mathrm{d}\tau^2``, so that ``\chi_\tau^2`` is the image of
``-\chi_q^2``:

````@example spheroidal_appendix
to_oblate(e) = e(c => -IM * c, q => IM * τ)
simplify.(to_oblate.((χᵩ^2, χₚ^2, -χ_q^2)) .- lame_coefficients(Oblate) .^ 2)
````

and it turns the volume of the prolate spheroid into that of the oblate one,
``\tfrac43\pi\,c^3\,\tau\,(\tau^2+1)``:

````@example spheroidal_appendix
factor(simplify(to_oblate(V)))
````

The operators are obtained in the same way, or directly from the oblate
system. Its Laplacian has ``q^2-1`` replaced by ``\tau^2+1`` and ``q^2-p^2``
by ``\tau^2+p^2``:

````@example spheroidal_appendix
fₒ = SymFunction("f", real = true)(ϕ, p, τ)
Δfₒ = (diff((1 - p^2) * diff(fₒ, p), p) + diff((τ^2 + 1) * diff(fₒ, τ), τ)) / (c^2 * (τ^2 + p^2)) +
    diff(fₒ, ϕ, 2) / (c^2 * (1 - p^2) * (τ^2 + 1))
simplify(laplacian(fₒ, Oblate) - Δfₒ)
````

The surface of the oblate spheroid follows from its own surface element:

````@example spheroidal_appendix
dSₒ = simplify(prod(lame_coefficients(Oblate)[1:2]))
Sₒ = simplify(integrate(integrate(dSₒ, (ϕ, 0, 2PI)), (p, -1, 1)))
````

that is ``2\pi\,c^2\sqrt{\tau^2+1}\,\bigl(\sqrt{\tau^2+1}+\tau^2\,\mathrm{arcsinh}(1/\tau)\bigr)``,
the expression used for the specific surface plotted above.

The harmonics of the two systems, and the change from a system to another of
a different focal distance, are the subject of
[Changing the focal distance of spheroidal harmonics](@ref tut-spheroidal-reexpansion).

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

