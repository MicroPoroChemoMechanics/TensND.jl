```@meta
EditURL = "../../../../scripts/23_insulating_torus.jl"
```

# [An insulating torus in a uniform flux](@id tut-insulating-torus)

!!! tip "See also"
    **Theory**: [Toroidal coordinates](@ref th-toroidal) · [Toroidal functions](@ref th-toroidal-harmonics) · [Curvilinear differential calculus](@ref th-curvilinear)

    **Manual**: [Coordinate systems (symbolic)](@ref man-coorsystems) · [Charts and coordinate conversions](@ref man-charts) · [Toroidal functions](@ref man-toroidal-harmonics)

A ring-shaped pore in a conducting medium, under a remote temperature
gradient along its axis: the temperature is harmonic outside the pore and
its flux does not cross the surface. The toroidal chart turns the surface of
the pore into a coordinate surface, and the problem into a tridiagonal
system on the amplitudes of a series of toroidal functions. This tutorial
sets it up and solves it with the chart and the functions of the package,
checks the result by automatic differentiation, and draws the cases of
[radi2016](@cite): the temperature around a torus of ratio ``r/R=0.5`` and
the temperature on its surface for several ratios.

````@example insulating_torus
using TensND
using LinearAlgebra
using ForwardDiff
using Printf
using Plots
gr()
nothing #hide
````

## The torus as a coordinate surface

The torus has the center-line radius ``R`` and the tube radius ``r``. In the
toroidal coordinates ``(\varphi,\alpha,\beta)`` of focal radius
``a=\sqrt{R^2-r^2}``, with ``\delta=\cosh\alpha-\cos\beta``,

```math
\varrho=\frac{a\,\sinh\alpha}{\delta},\qquad x_3=\frac{a\,\sin\beta}{\delta},
```

it is the surface ``\alpha=\alpha_0``, ``\cosh\alpha_0=R/r``; the medium is
``\alpha<\alpha_0``, which contains the axis ``\alpha=0`` and the points at
infinity ``\alpha=\beta=0``.

````@example insulating_torus
R, r = 1.0, 0.5
chart = ToroidalChart(R, r)
a, α₀ = sqrt(R^2 - r^2), acosh(R / r)
coords_to_cartesian(chart, (0.0, α₀, π)), (R - r, 0.0, 0.0)      ## the point of the torus nearest to the axis
````

## The remote field as a series

A function harmonic in the medium and bounded on the axis is a series of the
products ``\sqrt\delta\,P_{n-1/2}(\cosh\alpha)\sin n\beta``, where
``P_{n-1/2}`` is the toroidal function of the first kind
([Separation of the Laplace equation](@ref th-toroidal)). The remote field
``x_3`` of a unit gradient along the axis has the series, of the second kind,

```math
x_3=a\,\sqrt{\delta}\,\sum_{n\ge1}s_n\,Q_{n-1/2}(\cosh\alpha)\,\sin n\beta,
\qquad s_n=\frac{4\sqrt2}{\pi}\,n ,
```

which follows from that of ``1/\sqrt{2\delta}`` by a derivative with respect
to ``\beta``. [`toroidal_Q`](@ref) returns the functions of all degrees up
to ``N`` at once:

````@example insulating_torus
s(n) = 4 * sqrt(2) / π * n
function x₃_series(α, β; N = 60)
    Q = toroidal_Q(N, 0, α)
    return a * sqrt(cosh(α) - cos(β)) * sum(s(n) * Q[n + 1] * sin(n * β) for n in 1:N)
end
x₃_series(1.0, 0.7), a * sin(0.7) / (cosh(1.0) - cos(0.7))
````

## The condition on the torus

The temperature is the remote field plus a perturbation that vanishes at
infinity,

```math
T=a\,\sqrt{\delta}\sum_{n\ge1}f_n(\alpha)\,\sin n\beta,
\qquad f_n=s_n\,Q_{n-1/2}(\cosh\alpha)+b_n\,P_{n-1/2}(\cosh\alpha) ,
```

and the flux does not cross the torus: ``\partial T/\partial\alpha=0`` at
``\alpha_0``. With ``\partial\delta/\partial\alpha=\sinh\alpha``,

```math
\frac{\partial}{\partial\alpha}\Bigl(\sqrt\delta\sum_nf_n\sin n\beta\Bigr)
=\frac{1}{2\sqrt\delta}\Bigl(\sinh\alpha\sum_nf_n\sin n\beta+2\cosh\alpha\sum_nf_n'\sin n\beta-\sum_nf_n'\,2\cos\beta\,\sin n\beta\Bigr) ,
```

and ``2\cos\beta\sin n\beta=\sin(n+1)\beta+\sin(n-1)\beta``, the coefficient
of each ``\sin n\beta`` vanishes when, at ``\alpha_0``,

```math
\sinh\alpha_0\,f_n+2\cosh\alpha_0\,f_n'-f_{n-1}'-f_{n+1}'=0,\qquad n\ge1,\quad f_0'=0 .
```

Three consecutive amplitudes are tied. The unknowns are scaled as
``x_n=b_n\,P_{n-1/2}(\cosh\alpha_0)``, the functions of the first kind
growing as ``\mathrm{e}^{n\alpha_0}``; with ``\ell_n=P_{n-1/2}'/P_{n-1/2}``
the system is tridiagonal,

```math
(\sinh\alpha_0+2\cosh\alpha_0\,\ell_n)\,x_n-\ell_{n-1}\,x_{n-1}-\ell_{n+1}\,x_{n+1}
=-\bigl[s_n(\sinh\alpha_0\,Q_n+2\cosh\alpha_0\,Q_n')-s_{n-1}Q_{n-1}'-s_{n+1}Q_{n+1}'\bigr] ,
```

and is cut at the degree ``N``. [`toroidal_dP`](@ref) and
[`toroidal_dQ`](@ref) give the derivatives with respect to ``\alpha``.

````@example insulating_torus
function amplitudes(α₀; N = 40)
    P, dP = toroidal_P(N + 1, 0, α₀), toroidal_dP(N + 1, 0, α₀)
    Q, dQ = toroidal_Q(N + 1, 0, α₀), toroidal_dQ(N + 1, 0, α₀)
    ℓ = dP ./ P                                                 ## ℓ[n + 1] is ℓₙ
    diagonal = [sinh(α₀) + 2 * cosh(α₀) * ℓ[n + 1] for n in 1:N]
    lower = [-ℓ[n] for n in 2:N]                                ## the factor of xₙ₋₁ in the relation n
    upper = [-ℓ[n + 2] for n in 1:(N - 1)]                      ## the factor of xₙ₊₁
    right = [-(s(n) * (sinh(α₀) * Q[n + 1] + 2 * cosh(α₀) * dQ[n + 1]) - s(n - 1) * dQ[n] - s(n + 1) * dQ[n + 2]) for n in 1:N]
    return (Tridiagonal(lower, diagonal, upper) \ right) ./ P[2:(N + 1)]
end
b = amplitudes(α₀)
b[1:4]
````

The amplitudes decay as ``\mathrm{e}^{-2n\alpha_0}``. The temperature at a
point is the remote field plus the series of the perturbation;
[`cartesian_to_coords`](@ref) gives the coordinates of the point in the
chart.

````@example insulating_torus
function temperature(chart, b, x)
    _, α, β = cartesian_to_coords(chart, x)
    P = toroidal_P(length(b), 0, α)
    return x[3] + chart.a * sqrt(cosh(α) - cos(β)) * sum(b[n] * P[n + 1] * sin(n * β) for n in eachindex(b))
end
T(x) = temperature(chart, b, x)
T([1.8, 0.3, 0.4])
````

## Checks by automatic differentiation

The function is differentiated through the chart and the toroidal functions.
Its Laplacian vanishes in the medium, and its derivative along the normal
vanishes on the torus, here at a few points:

````@example insulating_torus
laplacian_at(x) = tr(ForwardDiff.hessian(T, x))
function normal_flux_at(φ, β)
    x = coords_to_cartesian(chart, (φ, α₀, β))
    ν = ForwardDiff.derivative(α -> coords_to_cartesian(chart, (φ, α, β)), α₀)      ## along 𝐞ᵅ, normal to the torus
    return dot(ForwardDiff.gradient(T, Vector(x)), ν / norm(ν))
end
maximum(abs, laplacian_at.([[1.8, 0.3, 0.4], [0.2, -0.1, 0.9], [0.1, 0.0, -0.2]])),
    maximum(abs, normal_flux_at(φ, β) for φ in (0.0, 1.0), β in (0.3, 1.5, 2.9))
````

## The temperature around the torus

The map of the temperature in a meridian plane, for ``r/R=0.5``, in units of
the remote gradient times ``R``: far from the torus it is ``x_3``, and the
isotherms bend to meet the surface at a right angle.

````@example insulating_torus
ϱs, zs = range(0.0, 2.2; length = 221), range(-1.4, 1.4; length = 281)
inside(ϱ, z) = (ϱ - R)^2 + z^2 < r^2
field = [inside(ϱ, z) ? NaN : T([ϱ, 0.0, z]) for z in zs, ϱ in ϱs]
contourf(ϱs, zs, field; levels = 28, color = :balance, aspect_ratio = 1, xlabel = "ϱ / R", ylabel = "x₃ / R", colorbar_title = "T / (G R)", size = (560, 560), linewidth = 0)
θs = range(0, 2π; length = 181)
plot!(R .+ r .* cos.(θs), r .* sin.(θs); color = :black, lw = 1.5, label = "")
````

## The temperature on the torus

On the surface the temperature is drawn against the angle ``\theta`` around
the tube, from the point farthest from the axis, for three ratios ``r/R``.
Without the pore it would be ``r\sin\theta``; around a thin ring it tends to
``2r\sin\theta``, the temperature on an insulating cylinder across the
gradient, and it loses its symmetry about ``\theta=\pi/2`` as the hole of
the torus closes.

````@example insulating_torus
plot(; xlabel = "θ", ylabel = "T / (G r)", legend = :topright, size = (640, 380))
for ratio in (0.25, 0.5, 0.75)
    torus = ToroidalChart(1.0, ratio)
    amplitude = amplitudes(acosh(1 / ratio); N = 80)
    on_surface(θ) = temperature(torus, amplitude, [1 + 1.000001 * ratio * cos(θ), 0.0, 1.000001 * ratio * sin(θ)]) / ratio
    plot!(θs, on_surface.(θs); lw = 2, label = "r/R = $ratio")
end
plot!(θs, sin.(θs); color = :black, ls = :dash, label = "no pore")
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

