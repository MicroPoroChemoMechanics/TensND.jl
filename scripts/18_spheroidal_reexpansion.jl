# # [Changing the focal distance of spheroidal harmonics](@id tut-spheroidal-reexpansion)
#
# !!! tip "See also"
#     **Theory**: [Changes of coordinates](@ref th-coordinate-changes) · [Spheroidal harmonics and their re-expansion](@ref th-spheroidal-harmonics) · [Spheroidal coordinates](@ref th-spheroidal)
#
#     **Manual**: [Charts and coordinate conversions](@ref man-charts) · [Spheroidal harmonics](@ref man-spheroidal-harmonics)
#
#     **Tutorials**: [Ellipsoidal coordinates](@ref tut-ellipsoidal) · [The spheroid in spheroidal coordinates](@ref tut-spheroidal-appendix)
#
# A field that is harmonic around a spheroid is a series of spheroidal
# harmonics of the system confocal with that spheroid. When a second
# spheroidal surface of the same axis is not confocal with the first one, a
# condition written on it needs the same field in the harmonics of its own
# system. This tutorial computes the harmonics of two coaxial systems, passes
# from one to the other, and measures where and how fast the series converge.
# The formulas are proved on
# [Spheroidal harmonics and their re-expansion](@ref th-spheroidal-harmonics);
# they are those that [riccardi1999](@citet) use for a spheroidal inclusion
# surrounded by a non-confocal spheroidal shell.

import Pkg                                                          #jl
Pkg.activate(joinpath(@__DIR__, "..", "docs"); io = devnull)        #jl
## The `docs` environment declares every dependency the tutorials use     #jl
## (SymPy, Symbolics, BenchmarkTools, ...) and resolves TensND             #jl
## itself through `[sources]`. These four lines are stripped from the      #jl
## generated page and notebook, which already run inside that project.     #jl

using TensND
using LinearAlgebra
using Printf
using Plots

# ## Two systems, one point
#
# The two prolate systems have the axis ``\underline{e}_3`` and the semi-focal
# distances ``c_1=1`` and ``c_2=1.5``. A point has the coordinates
# ``(\varphi,p_1,q_1)`` in the first and ``(\varphi,p_2,q_2)`` in the second:

chart₁ = ProlateSpheroidalChart(1.0)
chart₂ = ProlateSpheroidalChart(1.5)
x = [0.7, -0.4, 2.1]
ξ₁ = cartesian_to_coords(chart₁, x)

#-

ξ₂ = change_coords(chart₂, chart₁, ξ₁)

# The spheroids ``q_1=\text{const}`` and ``q_2=\text{const}`` through the point
# are different surfaces. In a meridian half-plane, with ``\varrho`` the
# distance to the axis:

meridian(chart, q) = (ps = range(-1, 1; length = 201); (chart.c .* sqrt.(1 .- ps .^ 2) .* sqrt(q^2 - 1), chart.c .* ps .* q))
ϱ = hypot(x[1], x[2])
plt = plot(; aspect_ratio = :equal, xlabel = "ϱ", ylabel = "x₃", legend = :outerright, size = (640, 460))
plot!(plt, meridian(chart₁, ξ₁[3])...; lw = 2, label = "q₁ = $(round(ξ₁[3]; digits = 3))")
plot!(plt, meridian(chart₂, ξ₂[3])...; lw = 2, label = "q₂ = $(round(ξ₂[3]; digits = 3))")
plot!(plt, [0, 0], [-1, 1]; lw = 4, color = 1, label = "focal segment, c₁ = 1")
plot!(plt, [0.03, 0.03], [-1.5, 1.5]; lw = 4, color = 2, label = "focal segment, c₂ = 1.5")
scatter!(plt, [ϱ], [x[3]]; color = :black, label = "point")

# ## Harmonics of a system
#
# The harmonics of order ``m`` of a prolate system are
# ``P_n^m(p)\,P_n^m(q)\,\cos m\varphi`` inside and
# ``P_n^m(p)\,Q_n^m(q)\,\cos m\varphi`` outside (and the same with
# ``\sin m\varphi``). [`interior_harmonics`](@ref) and
# [`exterior_harmonics`](@ref) return their meridian parts for the degrees
# ``0`` to ``N``, the degree ``n`` at the index ``n+1``:

interior_harmonics(chart₁, 4, 1, ξ₁)

#-

exterior_harmonics(chart₁, 4, 1, ξ₁)

# Each is harmonic. The pointwise system of the chart checks it: the Laplacian
# of the exterior harmonic ``n=3``, ``m=1`` vanishes to rounding at the point.

CS₁ = coorsys_spheroidal(PointwiseGeometry, 1.0)
G₃¹(ξ) = exterior_harmonics(chart₁, 3, 1, ξ)[4] * cos(ξ[1])
laplacian(G₃¹, CS₁)(Vector(ξ₁)), G₃¹(ξ₁)

# The exterior functions ``Q_n^m`` decay with the degree as fast as the
# interior ones grow, by a factor ``e^{\mp\xi}`` per degree with
# ``\xi=\mathrm{arccosh}\,q``:

q = ξ₁[3]
Qs = radial_Q(chart₁, 40, 0, q)
Ps = radial_P(chart₁, 40, 0, q)
plot(
    0:40, [abs.(Qs) Ps]; yscale = :log10, lw = 2, label = ["Qₙ(q)" "Pₙ(q)"],
    xlabel = "degree n", legend = :right, size = (640, 380),
)
plot!(0:40, exp.(-(0:40) .* acosh(q)); ls = :dash, color = :gray, label = "e^(∓nξ)")
plot!(0:40, exp.((0:40) .* acosh(q)); ls = :dash, color = :gray, label = "")

# ## Interior harmonics: a finite sum
#
# An interior harmonic is a polynomial of the Cartesian coordinates. In the
# other system it is a finite combination of the interior harmonics of degrees
# ``n, n-2,\ldots``,
#
# ```math
# P_n^m(p_1)\,P_n^m(q_1)=\sum_{k=0}^{\lfloor(n-m)/2\rfloor}L_k\,P_{n-2k}^m(p_2)\,P_{n-2k}^m(q_2),
# ```
#
# and [`interior_reexpansion`](@ref) returns the ``L_k``, the target chart
# first:

L = interior_reexpansion(chart₂, chart₁, 5, 1)

# The identity holds at every point, inside the spheroids as well as outside:

function interior_gap(to, from, n, m, x)
    lhs = interior_harmonics(from, n, m, cartesian_to_coords(from, x))[n + 1]
    rhs = sum(interior_reexpansion(to, from, n, m) .* interior_harmonics(to, n, m, cartesian_to_coords(to, x))[(n + 1):-2:(m + 1)])
    return abs(rhs / lhs - 1)
end
[interior_gap(chart₂, chart₁, 5, 1, y) for y in ([0.7, -0.4, 2.1], [0.1, 0.2, 0.3], [3.0, 1.0, -8.0])]

# The coefficients are polynomials of the ratio of the semi-focal distances.
# With symbols, for ``n=2``, ``m=0``:

using SymPyPythonCall
c₁, c₂ = symbols("c₁ c₂", positive = true)
simplify.(interior_reexpansion(ProlateSpheroidalChart(c₂), ProlateSpheroidalChart(c₁), 2, 0))

# ## Exterior harmonics: a series
#
# An exterior harmonic is a series of the exterior harmonics of degrees
# ``n, n+2,\ldots`` of the other system,
#
# ```math
# P_n^m(p_1)\,Q_n^m(q_1)=\sum_{i\ge0}K_i\,P_{n+2i}^m(p_2)\,Q_{n+2i}^m(q_2),
# ```
#
# whose coefficients [`exterior_reexpansion`](@ref) returns up to a given
# index:

K = exterior_reexpansion(chart₂, chart₁, 2, 1, 6)

# The error of the partial sums is measured at the point ``x`` above:

function exterior_errors(to, from, n, m, x, imax)
    exact = exterior_harmonics(from, n, m, cartesian_to_coords(from, x))[n + 1]
    terms = exterior_reexpansion(to, from, n, m, imax) .*
        exterior_harmonics(to, n + 2imax, m, cartesian_to_coords(to, x))[(n + 1):2:end]
    return abs.(cumsum(terms) ./ exact .- 1)
end
errors = exterior_errors(chart₂, chart₁, 2, 1, x, 25)
@printf "terms    error\n"
for N in (1, 3, 5, 10, 15, 20)
    @printf "%5d    %.1e\n" N errors[N]
end

# ## Where the series converges
#
# The series converges outside the smallest spheroid of the target system that
# encloses the focal segment of the source system.
# [`reexpansion_bound`](@ref) returns it as a value of the radial coordinate of
# the target. From ``c_1=1`` to ``c_2=1.5`` the focal segment of the source is
# inside that of the target, and the series converges everywhere off the
# latter; the other way, it converges for ``q_1>c_2/c_1``:

reexpansion_bound(chart₂, chart₁), reexpansion_bound(chart₁, chart₂)

# Beyond the bound the terms decay as ``e^{-2i(\xi-\xi^*)}``, where
# ``\xi=\mathrm{arccosh}\,q`` at the point and ``\xi^*`` at the bound. On the
# axis ``p=1`` of the first system, at three values of ``q_1`` beyond the bound
# and one below it:

q★ = reexpansion_bound(chart₁, chart₂)
plt = plot(;
    yscale = :log10, xlabel = "number of terms", ylabel = "relative error of the partial sum",
    legend = :outerright, size = (720, 420), ylims = (1.0e-16, 1.0e4),
)
for (k, q₁) in enumerate((1.3, 1.6, 2.0, 3.0))
    y = coords_to_cartesian(chart₁, (0.3, 0.6, q₁))
    e = exterior_errors(chart₁, chart₂, 2, 1, y, 40)
    plot!(plt, 1:41, max.(e, 1.0e-17); lw = 2, color = k, label = "q₁ = $q₁")
    rate = exp(-2 * (acosh(q₁) - acosh(q★)))
    plot!(plt, 1:41, e[1] .* rate .^ (0:40); ls = :dash, color = k, label = "")
end
plt

# The dashed lines are the predicted geometric rates. At ``q_1=1.3``, inside the
# spheroid ``q_1=1.5`` through the foci of the second system, the series
# diverges although the harmonic itself is regular there.
#
# The error after 20 terms, in a meridian half-plane, draws the bound: the
# white curve is the spheroid ``q_1=c_2/c_1``.

ϱs = range(0.02, 3.0; length = 90)
zs = range(-3.0, 3.0; length = 180)
err20(ϱ, z) = log10(clamp(exterior_errors(chart₁, chart₂, 2, 1, [ϱ, 0.0, z], 20)[end], 1.0e-16, 1.0e4))
heatmap(
    ϱs, zs, (ϱ, z) -> err20(ϱ, z); aspect_ratio = :equal, c = :viridis, clims = (-16, 4),
    xlabel = "ϱ", ylabel = "x₃", colorbar_title = "log₁₀ error after 20 terms", size = (520, 620),
)
plot!(meridian(chart₁, q★)...; color = :white, lw = 2, label = "")
plot!([0, 0], [-1.5, 1.5]; color = :red, lw = 3, label = "")

# ## Prolate, oblate and spherical
#
# The two systems may be prolate, oblate or spherical in any combination. An
# exterior harmonic of an oblate system (a disk-like inclusion) in the
# harmonics of a prolate one:

oblate = OblateSpheroidalChart(0.8)
y = [1.9, 0.6, 1.4]
reexpansion_bound(chart₁, oblate), cartesian_to_coords(chart₁, y)[3]

#-

exterior_errors(chart₁, oblate, 1, 0, y, 30)[[1, 5, 10, 20, 30]]

# In the spherical chart, the harmonics are the solid harmonics
# ``r^{-n-1}P_n^m(\cos\theta)``: the coefficients are the multipole moments of
# the spheroidal harmonic, and the series converges outside the sphere through
# the foci.

exterior_reexpansion(SphericalChart(), chart₁, 0, 0, 4)

# The first one is ``c_1``: far from the foci, ``Q_0(q)\simeq1/q\simeq c_1/r``.
# The interior harmonics, in the other direction, are finite sums: the solid
# harmonic ``r^4P_4(\cos\theta)`` in the interior harmonics of the oblate
# system,

interior_reexpansion(oblate, SphericalChart(), 4, 0)

# ## Derivatives with respect to the focal distance
#
# The coefficients are differentiable. The sensitivity of the first
# coefficients of an exterior re-expansion to the semi-focal distance of the
# source, by automatic differentiation, against the closed form of
# ``K_0=(c_1/c_2)^{n+1}``:

using ForwardDiff
dK = ForwardDiff.derivative(c -> exterior_reexpansion(chart₂, ProlateSpheroidalChart(c), 2, 1, 3), 1.0)
dK[1], 3 * 1.0^2 / 1.5^3
