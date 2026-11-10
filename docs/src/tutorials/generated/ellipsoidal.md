```@meta
EditURL = "../../../../scripts/17_ellipsoidal.jl"
```

# [Ellipsoidal coordinates](@id tut-ellipsoidal)

!!! tip "See also"
    **Theory**: [Spheroidal coordinates](@ref th-spheroidal) · [Ellipsoidal coordinates](@ref th-ellipsoidal) · [Curvilinear differential calculus](@ref th-curvilinear) · [Changes of coordinates](@ref th-coordinate-changes)

    **Manual**: [Coordinate systems (symbolic)](@ref man-coorsystems) · [Charts and coordinate conversions](@ref man-charts)

    **Tutorials**: [The spheroid in spheroidal coordinates](@ref tut-spheroidal-appendix) · [Changing the focal distance of spheroidal harmonics](@ref tut-spheroidal-reexpansion)

The ellipsoidal coordinates are the system in which a triaxial ellipsoid is a
coordinate surface, as a spheroid is in the spheroidal systems. This tutorial
builds the system of semi-focal distances ``h_2>h_3`` with `TensND`: the
three families of confocal quadrics through a point, the exact system with
its Lamé coefficients and unit vectors, the differential operators, the
ellipsoid seen as a surface, the pointwise system with the inversion of the
chart, and the spheroidal limits. The formulas are those of
[Ellipsoidal coordinates](@ref th-ellipsoidal), in the notation of
[dassios2012](@citet).

````@example ellipsoidal
using TensND
using LinearAlgebra
using SymPyPythonCall
using ForwardDiff
using Printf
using Plots
````

## The confocal quadrics

An ellipsoid of semi-axes ``a_1>a_2>a_3`` along
``\underline{e}_1,\underline{e}_2,\underline{e}_3`` has the semi-focal
distances ``h_2^2=a_1^2-a_3^2`` and ``h_3^2=a_1^2-a_2^2``. Through a point
``\underline{x}`` pass three quadrics confocal with it, whose parameters
``\kappa=\rho^2,\mu^2,\nu^2`` are the roots of

```math
\frac{x_1^2}{\kappa}+\frac{x_2^2}{\kappa-h_3^2}+\frac{x_3^2}{\kappa-h_2^2}=1,
\qquad
0\le\nu^2\le h_3^2\le\mu^2\le h_2^2\le\rho^2 .
```

The coordinates ``(\rho,\mu,\nu)`` are those of
[`EllipsoidalChart`](@ref). For the ellipsoid of semi-axes ``(3,2,1)`` and a
point inside it, [`cartesian_to_coords`](@ref) returns the three roots:

````@example ellipsoidal
chart = EllipsoidalChart(3.0, 2.0, 1.0)
H₂, H₃ = chart.h₂, chart.h₃
x = [1.8, 1.2, 0.5]
ξ₀ = cartesian_to_coords(chart, x)
````

They are ordered as announced, and each one satisfies the equation above:

````@example ellipsoidal
quadric(κ) = x[1]^2 / κ + x[2]^2 / (κ - H₃^2) + x[3]^2 / (κ - H₂^2)
ρ₀, μ₀, ν₀ = ξ₀
(ν₀^2 ≤ H₃^2 ≤ μ₀^2 ≤ H₂^2 ≤ ρ₀^2, all(isapprox.(quadric.((ρ₀^2, μ₀^2, ν₀^2)), 1)))
````

The signs of the denominators distinguish the families. For
``\kappa=\rho^2`` all are positive: the surface is the ellipsoid of semi-axes
``\rho,\sqrt{\rho^2-h_3^2},\sqrt{\rho^2-h_2^2}``. For ``\kappa=\mu^2`` the
denominator of ``x_3^2`` is negative, a hyperboloid of one sheet, and for
``\kappa=\nu^2`` those of ``x_2^2`` and ``x_3^2`` are, a hyperboloid of two
sheets.

In the plane ``x_3=0`` the term in ``x_3`` drops out, and the three families
trace conics with the same foci ``(\pm h_3,0)``: the ellipses
``x_1^2/s^2+x_2^2/(s^2-h_3^2)=1`` with ``s=\rho``, outside the focal ellipse
``s=h_2``, and with ``s=\mu``, inside it, and the hyperbolas
``x_1^2/\nu^2-x_2^2/(h_3^2-\nu^2)=1``. A point of the plane has
``\mu=h_2`` outside the focal ellipse and ``\rho=h_2`` inside:

````@example ellipsoidal
(cartesian_to_coords(chart, [4cos(1.0), sqrt(16 - H₃^2) * sin(1.0), 0.0]),
    cartesian_to_coords(chart, [1.0, 1.0, 0.0]))
````

The figure draws the traces of the three families, the foci being the two
black dots. The ellipse ``\rho=3`` is the section of the ellipsoid of
semi-axes ``(3,2,1)`` by the plane.

````@example ellipsoidal
ellipse(s) = (t = range(0, 2π; length = 200); (s * cos.(t), sqrt(s^2 - H₃^2) * sin.(t)))
function hyperbola(ν)
    T = acosh(5 / abs(ν))
    t = range(-T, T; length = 200)
    return ν * cosh.(t), sqrt(H₃^2 - ν^2) * sinh.(t)
end
ρs = (3.0, 3.6, 4.4, 5.2)
μs = range(H₃, H₂; length = 6)[2:5]
νs = (-1.9, -1.4, -0.9, -0.4, 0.4, 0.9, 1.4, 1.9)
fig = plot(;
    aspect_ratio = :equal, xlims = (-5.2, 5.2), ylims = (-3.8, 3.8), size = (720, 480),
    xlabel = "x₁", ylabel = "x₂", legend = :outertopright, guidefontsize = 9,
)
for (k, s) in enumerate(ρs)
    plot!(fig, ellipse(s)...; color = 1, lw = 1.5, label = k == 1 ? "ρ = const" : "")
end
for (k, s) in enumerate(μs)
    plot!(fig, ellipse(s)...; color = 2, lw = 1.5, label = k == 1 ? "μ = const" : "")
end
for (k, ν) in enumerate(νs)
    plot!(fig, hyperbola(ν)...; color = 3, lw = 1.5, label = k == 1 ? "ν = const" : "")
end
scatter!(fig, [-H₃, H₃], [0, 0]; color = :black, label = "")
````

## The exact system

[`coorsys_ellipsoidal`](@ref) returns the system for symbolic ``h_2`` and
``h_3``. Its position vector is written with six auxiliary symbols, the
square roots that the chart involves,

```math
\bar\rho_i=\sqrt{\rho^2-h_i^2},\qquad
\bar\mu_2=\sqrt{h_2^2-\mu^2},\quad\bar\mu_3=\sqrt{\mu^2-h_3^2},\qquad
\bar\nu_i=\sqrt{h_i^2-\nu^2},
```

so that
``x_1=\rho\mu\nu/(h_2h_3)``, ``x_2=\bar\rho_3\bar\mu_3\bar\nu_3/(h_1h_3)``,
``x_3=\bar\rho_2\bar\mu_2\bar\nu_2/(h_1h_2)`` with
``h_1=\sqrt{h_2^2-h_3^2}``:

````@example ellipsoidal
E = coorsys_ellipsoidal()
ρ, μ, ν = get_coords(E)
h₂, h₃ = symbols("h₂ h₃", positive = true)
ρ̄₂, ρ̄₃, μ̄₂, μ̄₃, ν̄₂, ν̄₃ = symbols("ρ̄₂ ρ̄₃ μ̄₂ μ̄₃ ν̄₂ ν̄₃", positive = true)
position_vector(E)
````

The Lamé coefficients are

```math
\chi_\rho^2=\frac{(\rho^2-\mu^2)(\rho^2-\nu^2)}{(\rho^2-h_3^2)(\rho^2-h_2^2)},
\quad
\chi_\mu^2=\frac{(\rho^2-\mu^2)(\mu^2-\nu^2)}{(\mu^2-h_3^2)(h_2^2-\mu^2)},
\quad
\chi_\nu^2=\frac{(\rho^2-\nu^2)(\mu^2-\nu^2)}{(h_3^2-\nu^2)(h_2^2-\nu^2)} :
```

````@example ellipsoidal
χᵨ, χᵤ, χᵥ = lame_coefficients(E)
χ² = (
    (ρ^2 - μ^2) * (ρ^2 - ν^2) / ((ρ^2 - h₃^2) * (ρ^2 - h₂^2)),
    (ρ^2 - μ^2) * (μ^2 - ν^2) / ((μ^2 - h₃^2) * (h₂^2 - μ^2)),
    (ρ^2 - ν^2) * (μ^2 - ν^2) / ((h₃^2 - ν^2) * (h₂^2 - ν^2)),
)
simplify.((χᵨ, χᵤ, χᵥ) .^ 2 .- χ²)
````

The unit vectors are the columns of the normalized basis. Their expressions
are long, and simplifying a whole matrix of them is slow; the orthonormality
is checked instead at the point ``(h_2,h_3)=(4,2)``,
``(\rho,\mu,\nu)=(5,3,1)``, where every auxiliary symbol is an exact number:

````@example ellipsoidal
pt = [
    ρ => 5, μ => 3, ν => 1, h₂ => 4, h₃ => 2,
    ρ̄₂ => 3, ρ̄₃ => sqrt(Sym(21)), μ̄₂ => sqrt(Sym(7)), μ̄₃ => sqrt(Sym(5)), ν̄₂ => sqrt(Sym(15)), ν̄₃ => sqrt(Sym(3)),
]
at(e) = simplify(e(pt...))
ℬ = normalized_basis(E)
M = at.(Matrix(ℬ))
simplify.(transpose(M) * M)
````

The triad ``(\underline{e}^\rho,\underline{e}^\mu,\underline{e}^\nu)`` is
left-handed, so that its determinant is ``-1``:

````@example ellipsoidal
simplify(det(M))
````

At that point the position vector, evaluated from the auxiliary symbols,
is that of the chart:

````@example ellipsoidal
at.(position_vector(E)) - [Sym(v) for v in coords_to_cartesian(EllipsoidalChart(Sym(4), Sym(2)), (Sym(5), Sym(3), Sym(1)))]
````

## Gradient and Laplacian

For a scalar field ``f(\rho,\mu,\nu)`` the physical components of the gradient
are the partial derivatives divided by the Lamé coefficients:

````@example ellipsoidal
f = SymFunction("f", real = true)(ρ, μ, ν)
∇f = simplify.(components(grad(f, E), ℬ))
simplify.(∇f - [diff(f, ρ) / χᵨ, diff(f, μ) / χᵤ, diff(f, ν) / χᵥ])
````

With the square roots
``R_i=\sqrt{\rho^2-h_i^2}``, ``M_2=\sqrt{h_2^2-\mu^2}``, ``M_3=\sqrt{\mu^2-h_3^2}``,
``N_i=\sqrt{h_i^2-\nu^2}``, the Laplacian is

```math
\begin{aligned}
\Delta f={}&\frac{R_2R_3}{(\rho^2-\mu^2)(\rho^2-\nu^2)}\,
\frac{\partial}{\partial\rho}\Bigl(R_2R_3\,\frac{\partial f}{\partial\rho}\Bigr)
+\frac{M_2M_3}{(\rho^2-\mu^2)(\mu^2-\nu^2)}\,
\frac{\partial}{\partial\mu}\Bigl(M_2M_3\,\frac{\partial f}{\partial\mu}\Bigr)\\
&+\frac{N_2N_3}{(\rho^2-\nu^2)(\mu^2-\nu^2)}\,
\frac{\partial}{\partial\nu}\Bigl(N_2N_3\,\frac{\partial f}{\partial\nu}\Bigr) .
\end{aligned}
```

````@example ellipsoidal
R₂, R₃ = sqrt(ρ^2 - h₂^2), sqrt(ρ^2 - h₃^2)
M₂, M₃ = sqrt(h₂^2 - μ^2), sqrt(μ^2 - h₃^2)
N₂, N₃ = sqrt(h₂^2 - ν^2), sqrt(h₃^2 - ν^2)
Δρ(f) = R₂ * R₃ / ((ρ^2 - μ^2) * (ρ^2 - ν^2)) * diff(R₂ * R₃ * diff(f, ρ), ρ)
Δμ(f) = M₂ * M₃ / ((ρ^2 - μ^2) * (μ^2 - ν^2)) * diff(M₂ * M₃ * diff(f, μ), μ)
Δν(f) = N₂ * N₃ / ((ρ^2 - ν^2) * (μ^2 - ν^2)) * diff(N₂ * N₃ * diff(f, ν), ν)
simplify(laplacian(f, E) - (Δρ(f) + Δμ(f) + Δν(f)))
````

The harmonic polynomials ``x_1^2-x_3^2``, ``x_2x_3`` and ``x_1x_2x_3`` have a
vanishing Laplacian. The position vector provides ``x_1,x_2,x_3``. Simplifying
the Laplacian as a function of the coordinates takes about a minute for each
polynomial, so it is evaluated at two points with exact coordinates, with
``(h_2,h_3)=(4,2)``, one per column of the result:

````@example ellipsoidal
x₁, x₂, x₃ = position_vector(E)
points = (
    [ρ => 5, μ => 3, ν => 1, h₂ => 4, h₃ => 2],
    [ρ => 6, μ => Sym(7) / 2, ν => Sym(3) / 2, h₂ => 4, h₃ => 2],
)
[simplify(laplacian(P, E)(q...)) for P in (x₁^2 - x₃^2, x₂ * x₃, x₁ * x₂ * x₃), q in points]
````

## The potential of a charged ellipsoid

A conducting ellipsoid ``\rho=a_1`` at the potential ``V`` carries an
exterior potential that depends on ``\rho`` only, since ``\rho=\text{const}``
is an equipotential. For a field ``f(\rho)`` the Laplacian reduces to its
first term:

````@example ellipsoidal
fρ = SymFunction("f", real = true)(ρ)
simplify(laplacian(fρ, E) - Δρ(fρ))
````

so that ``\Delta f=0`` is
``\dfrac{\mathrm{d}}{\mathrm{d}\rho}\bigl(R_2R_3\,f'\bigr)=0``, and
``f'=C/(R_2R_3)``. The decreasing solution is the elliptic integral

```math
I(\rho)=\int_\rho^\infty\frac{\mathrm{d}t}{\sqrt{t^2-h_3^2}\,\sqrt{t^2-h_2^2}},
\qquad
I'(\rho)=-\frac{1}{R_2R_3},
```

the derivative of an integral with respect to its lower limit. The Laplacian
of ``I`` follows from that of a function of ``\rho``, with ``I'`` and
``I''`` replaced by their values, without evaluating the integral:

````@example ellipsoidal
I′ = -1 / (R₂ * R₃)
ΔI = laplacian(fρ, E)(diff(fρ, ρ, 2) => diff(I′, ρ), diff(fρ, ρ) => I′)
simplify(ΔI)
````

## The ellipsoid as a surface

On the ellipsoid ``\rho=\text{const}`` the unit normal is
``\underline{e}^\rho``. As for the spheroid, the gradient of the normal field
read in the normalized basis is minus the curvature tensor of the surface,
and its diagonal holds, with the sign changed, the principal curvatures along
``\underline{e}^\mu`` and ``\underline{e}^\nu``:

```math
\kappa_\mu=-\frac{\rho}{\chi_\rho\,(\rho^2-\mu^2)},
\qquad
\kappa_\nu=-\frac{\rho}{\chi_\rho\,(\rho^2-\nu^2)} .
```

````@example ellipsoidal
𝐞ᵨ = unitvec(E)[1]
∇𝐧 = components(grad(𝐞ᵨ, E), ℬ)
κᵤ, κᵥ = -∇𝐧[2, 2], -∇𝐧[3, 3]
simplify(κᵤ + ρ / (χᵨ * (ρ^2 - μ^2))), simplify(κᵥ + ρ / (χᵨ * (ρ^2 - ν^2)))
````

The first row of the gradient vanishes because the normal has a unit norm,
and the off-diagonal entries of the tangent block vanish because the
coordinate lines of the ellipsoid are its lines of curvature:

````@example ellipsoidal
simplify.(∇𝐧[1, :]), simplify(∇𝐧[2, 3]), simplify(∇𝐧[3, 2])
````

The product of the two curvatures is the Gaussian curvature of the ellipsoid,

```math
K=\frac{\rho^2\,(\rho^2-h_3^2)\,(\rho^2-h_2^2)}{(\rho^2-\mu^2)^2\,(\rho^2-\nu^2)^2} :
```

````@example ellipsoidal
simplify(κᵤ * κᵥ - ρ^2 * (ρ^2 - h₃^2) * (ρ^2 - h₂^2) / ((ρ^2 - μ^2)^2 * (ρ^2 - ν^2)^2))
````

At the three vertices of the ellipsoid of semi-axes ``(3,2,1)`` the principal
curvatures are known elementary values: at the vertex on the axis ``i``,
``a_i/a_j^2`` and ``a_i/a_k^2``. The coordinates of the vertices come from
the chart, and the curvatures from the formulas above:

````@example ellipsoidal
a = (3.0, 2.0, 1.0)
curvatures(ρ, μ, ν) = (
    χρ = sqrt((ρ^2 - μ^2) * (ρ^2 - ν^2) / ((ρ^2 - H₃^2) * (ρ^2 - H₂^2)));
    sort([ρ / (χρ * (ρ^2 - μ^2)), ρ / (χρ * (ρ^2 - ν^2))])
)
vertices = ([3.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 1.0])
elementary(i) = sort([a[i] / a[j]^2 for j in 1:3 if j != i])
[curvatures(cartesian_to_coords(chart, v)...) ≈ elementary(i) for (i, v) in enumerate(vertices)]
````

## The pointwise system and the inversion of the chart

`coorsys_ellipsoidal(PointwiseGeometry, h₂, h₃)` evaluates the operators at a
point, with derivatives by automatic differentiation. The field
``F=x_1^2x_2+x_3^4`` is written in Cartesian components and composed with
the chart; its Laplacian in the Cartesian system is
``2x_2+12x_3^2``:

````@example ellipsoidal
Ep = coorsys_ellipsoidal(PointwiseGeometry, H₂, H₃)
F(x) = x[1]^2 * x[2] + x[3]^4
G(ξ) = F(coords_to_cartesian(chart, ξ))
laplacian(G, Ep)(ξ₀), 2x[2] + 12x[3]^2
````

The inverse map is differentiated by `ForwardDiff` too. The rows of its
Jacobian are the gradients of the coordinates, ``\underline{e}^i/\chi_i``,
that is ``\mathrm{diag}(1/\chi)\,R^{\mathsf T}`` with ``R`` the matrix of the
unit vectors:

````@example ellipsoidal
J = ForwardDiff.jacobian(y -> collect(cartesian_to_coords(chart, y)), x)
Jₑ = Diagonal(1 ./ lame_coefficients(Ep, ξ₀)) * transpose(Matrix(normalized_basis(Ep, ξ₀)))
J ≈ Jₑ
````

The chart covers the quarter space ``x_2\ge0``, ``x_3\ge0``. A point of
another quarter has the coordinates of its mirror image in the planes
``x_2=0`` and ``x_3=0``, and ``\nu`` keeps the sign of ``x_1``:

````@example ellipsoidal
xₘ = [1.8, -1.2, -0.5]
cartesian_to_coords(chart, xₘ) ≈ ξ₀, coords_to_cartesian(chart, cartesian_to_coords(chart, xₘ))
````

## The spheroidal limits

When ``h_3\to0`` the ellipsoids ``\rho=\text{const}`` tend to oblate
spheroids of axis ``\underline{e}_3`` and semi-focal distance ``c=h_2``, and
the relations ``\tau=\sqrt{\rho^2-h_2^2}/h_2`` and
``p=\sqrt{h_2^2-\mu^2}/h_2`` connect the coordinates to those of
[`OblateSpheroidalChart`](@ref). The table gives, for the point
``(0.7,0.5,0.6)`` and ``h_2=1``, the differences between these values and
the coordinates of the oblate chart:

````@example ellipsoidal
xₗ = [0.7, 0.5, 0.6]
_, pₒ, τₒ = cartesian_to_coords(OblateSpheroidalChart(1.0), xₗ)
println("      h₃        τ - τₒ        p - pₒ")
for h in (0.4, 0.2, 0.1, 0.05, 0.025)
    ρₗ, μₗ, _ = cartesian_to_coords(EllipsoidalChart(1.0, h), xₗ)
    @printf "  %6.3f   %+.3e   %+.3e\n" h sqrt(ρₗ^2 - 1) - τₒ sqrt(1 - μₗ^2) - pₒ
end
````

When ``h_3\to h_2`` they tend to prolate spheroids, of axis
``\underline{e}_1``. The spheroidal charts have their axis along
``\underline{e}_3``: the components of the point are permuted,
``(x_1,x_2,x_3)\mapsto(x_2,x_3,x_1)``, which keeps the orientation. The
relations are ``q=\rho/h_2`` and ``p=\nu/h_2``:

````@example ellipsoidal
_, pₚ, qₚ = cartesian_to_coords(ProlateSpheroidalChart(1.0), xₗ[[2, 3, 1]])
println("  h₂ - h₃        q - qₚ        p - pₚ")
for ε in (0.4, 0.2, 0.1, 0.05, 0.025)
    ρₗ, _, νₗ = cartesian_to_coords(EllipsoidalChart(1.0, 1 - ε), xₗ)
    @printf "  %6.3f   %+.3e   %+.3e\n" ε ρₗ - qₚ νₗ - pₚ
end
````

The differences vanish as ``h_3^2`` in the oblate case, and as ``h_2-h_3`` in
the prolate one. The chart itself is singular in the limit: a spheroid is
handled by its own system.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

