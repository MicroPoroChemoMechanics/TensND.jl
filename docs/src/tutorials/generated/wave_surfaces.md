```@meta
EditURL = "../../../../scripts/25_wave_surfaces.jl"
```

# [Slowness and wave surfaces](@id tut-wave-surfaces)

!!! tip "See also"
    **Theory**: [Plane waves](@ref th-plane-waves)

    **Manual**: [Plane waves](@ref man-plane-waves)

    **Tutorials**: [Thomsen's parameters](@ref tut-thomsen) · [Generalized Thomsen parameters](@ref tut-generalized-thomsen) · [Reflection and transmission at an interface](@ref tut-reflection)

In an anisotropic solid the velocity of a wave depends on its direction, and
its energy does not travel along the normal to its fronts. Three surfaces
describe this, one sheet per wave: the phase velocities, the slownesses and
the wave surface, the front that a point source emits. This tutorial
computes them for two rocks of [mensch1997](@cite): a transversely isotropic
sandstone and a limestone without any symmetry.

````@example wave_surfaces
using TensND
using LinearAlgebra
using ForwardDiff
using JSON
using Printf
using Plots
using SymPyPythonCall: symbols
gr()
nothing #hide
````

## Three surfaces

A plane wave ``\underline{u} = \underline{U}\,f(\underline{n}\cdot\underline{x} - Vt)``
of unit direction ``\underline{n}`` travels in the medium of stiffness
``\mathbb{C}`` and density ``\rho`` when

```math
\bigl(\boldsymbol{\Gamma}(\underline{n}) - \rho V^2\,\boldsymbol{1}\bigr)\cdot\underline{U} = \underline{0},
\qquad \boldsymbol{\Gamma}(\underline{n}) = \underline{n}\cdot\mathbb{C}\cdot\underline{n} ,
```

which gives three phase velocities per direction, qP, qS1 and qS2 in
decreasing order ([Plane waves](@ref th-plane-waves)). As ``\underline{n}``
runs over the unit sphere,

- the points ``V\underline{n}`` draw the **phase-velocity surface**;
- the points ``\underline{s} = \underline{n}/V`` draw the **slowness
  surface**: the front of the wave at the time ``t`` is the plane
  ``\underline{s}\cdot\underline{x} = t``;
- the points ``\underline{g}``, the group velocities, draw the **wave
  surface**, the envelope of the fronts ``\underline{n}\cdot\underline{x} = V``
  at unit time: what a pulse emitted at the origin has reached after one
  second. With a unit polarization,

```math
\underline{g} = \frac{(\underline{U}\cdot\mathbb{C}\cdot\underline{U})\cdot\underline{n}}{\rho V},
\qquad \underline{g}\cdot\underline{n} = V,
\qquad \underline{g}\cdot\underline{s} = 1 ,
```

and ``\underline{g}`` is normal to the slowness surface. In an isotropic
medium the three surfaces are spheres and ``\underline{g} = V\underline{n}``.

## The two rocks

The stiffnesses are given by their matrices in the two-index notation, in
GPa, and the densities in kg/m³, stored with their source under
`data/literature`. [`stiffness_from_voigt`](@ref) builds the tensors; the
densities are divided by 1000, so that the velocities come out in km/s and
the slownesses in s/km.

````@example wave_surfaces
data = JSON.parsefile(joinpath(pkgdir(TensND), "data", "literature", "menschGJI1997.json"))
quantity(key) = data["quantities"][key]["value"]
voigt(key) = Float64[row[c] for row in data["tables"][key], c in ("1", "2", "3", "4", "5", "6")]
sandstone = (ℂ = stiffness_from_voigt(voigt("mesaverde_sandstone_stiffness")), ρ = quantity("mesaverde_sandstone_density") / 1000)
limestone = (ℂ = stiffness_from_voigt(voigt("lavoux_stiffness")), ρ = quantity("lavoux_density") / 1000)
voigt_stiffness(sandstone.ℂ)

Main.show_literature("menschGJI1997") # hide
````

The sandstone is transversely isotropic about ``\underline{e}_3``: its matrix
has ``C_{11} = C_{22}``, ``C_{13} = C_{23}``, ``C_{44} = C_{55}``,
``C_{66} = (C_{11} - C_{12})/2`` and no other nonzero entry. The departures
from these equalities:

````@example wave_surfaces
C = voigt_stiffness(sandstone.ℂ)
C[1, 1] - C[2, 2], C[1, 3] - C[2, 3], C[4, 4] - C[5, 5], C[6, 6] - (C[1, 1] - C[1, 2]) / 2
````

## One direction

[`direction_from_angles`](@ref) gives the unit vector of the angle
``\theta`` to ``\underline{e}_3`` and of the azimuth ``\phi``.
[`plane_waves`](@ref) returns the three waves of that direction, with their
polarizations and group velocities, one per column:

````@example wave_surfaces
n = direction_from_angles(deg2rad(50), 0.0)
w = plane_waves(sandstone.ℂ, sandstone.ρ, n)
````

````@example wave_surfaces
w.polarization, w.group
````

The projection of each group velocity on the normal is the phase velocity,
and the ray leaves the normal by an angle that depends on the wave:

````@example wave_surfaces
transpose(w.group) * n - w.velocity, [acosd(dot(g, n) / norm(g)) for g in eachcol(w.group)]
````

## The sandstone in a plane through its axis

In a transversely isotropic medium every plane through the axis is a plane of
symmetry, and the group velocities of the directions of that plane stay in
it. One of the two shear waves, SH, is polarized normal to the plane, along
``\underline{e}_\phi = (-\sin\phi, \cos\phi, 0)``; the other, qSV, in it. The
sorting by velocity mixes them where their velocities cross, so the modes
are named by their polarizations:

````@example wave_surfaces
function ti_modes(w, ϕ)
    eϕ = [-sin(ϕ), cos(ϕ), 0.0]
    sh = abs(dot(w.polarization[:, 2], eϕ)) > abs(dot(w.polarization[:, 3], eϕ)) ? 2 : 3
    return (qP = 1, qSV = 5 - sh, SH = sh)
end
θs = range(0, 2π; length = 721)
meridian = [plane_waves(sandstone.ℂ, sandstone.ρ, direction_from_angles(θ, 0.0)) for θ in θs]
modes = [ti_modes(w, 0.0) for w in meridian]
minimum(abs(w.polarization[2, m.SH]) for (w, m) in zip(meridian, modes))
````

The last number is the smallest component of the SH polarizations along
``\underline{e}_2 = \underline{e}_\phi``. The three sections, in the plane
``(\underline{e}_1, \underline{e}_3)``:

````@example wave_surfaces
function section(point, title, unit)
    panel = plot(; aspect_ratio = :equal, title, xlabel = "x₁ ($unit)", ylabel = "x₃ ($unit)", legend = :topright, titlefontsize = 10)
    for k in (:qP, :qSV, :SH)
        P = [point(w, m[k]) for (w, m) in zip(meridian, modes)]
        plot!(panel, first.(P), last.(P); lw = 1.5, label = string(k))
    end
    return panel
end
plot(
    section((w, i) -> w.velocity[i] * w.direction[[1, 3]], "phase velocity", "km/s"),
    section((w, i) -> w.direction[[1, 3]] / w.velocity[i], "slowness", "s/km"),
    section((w, i) -> w.group[[1, 3], i], "wave surface", "km/s");
    layout = (1, 3), size = (980, 380), left_margin = 3Plots.mm, bottom_margin = 3Plots.mm,
)
````

The slowness curve of the meridian is known in closed form. A slowness
vector ``\underline{s} = (s_1, 0, s_3)`` lies on it when
``\det(\boldsymbol{\Gamma}(\underline{s}) - \rho\boldsymbol{1}) = 0``, and
[`christoffel_invariants`](@ref), which does not normalize its vector, gives
this determinant as ``\rho^3(I_1 - I_2 + I_3 - 1)`` for the invariants of
``\boldsymbol{\Gamma}(\underline{s})/\rho``. With the five components of a
stiffness transversely isotropic about ``\underline{e}_3`` as symbols, it
factors into the ellipse of SH times a quartic, the curves of qP and qSV:

````@example wave_surfaces
C₁₁, C₁₂, C₁₃, C₃₃, C₄₄, ρₛ = symbols("C_11 C_12 C_13 C_33 C_44 ρ", positive = true)
s₁, s₃ = symbols("s_1 s_3", real = true)
I₁, I₂, I₃ = christoffel_invariants(ti_from_constants(C₁₁, C₁₂, C₁₃, C₃₃, C₄₄, [0, 0, 1]), ρₛ, [s₁, 0, s₃])
slowness_curves = tfactor(ρₛ^3 * (I₁ - I₂ + I₃ - 1))
````

The two factors, with ``C_{66} = (C_{11} - C_{12})/2``: the difference
reduces to zero.

````@example wave_surfaces
SH = (C₁₁ - C₁₂) / 2 * s₁^2 + C₄₄ * s₃^2 - ρₛ
qP_qSV = (C₁₁ * s₁^2 + C₄₄ * s₃^2 - ρₛ) * (C₄₄ * s₁^2 + C₃₃ * s₃^2 - ρₛ) - (C₁₃ + C₄₄)^2 * s₁^2 * s₃^2
tsimplify(slowness_curves - SH * qP_qSV)
````

## Cusps of the qSV wave surface

Where the slowness curve of a wave turns from convex to concave, its normal,
that is the ray, turns back: the ray angle ``\psi`` to the axis is no longer
an increasing function of the angle ``\theta`` of the normal, and the wave
surface folds into a triangle with two cusps, crossed three times by the rays
that it covers. The ray angles of the three modes against ``\theta``, on a
quarter of the meridian:

````@example wave_surfaces
θq = range(0.001, π / 2 - 0.001; length = 900)
rays = map(θq) do θ
    w = plane_waves(sandstone.ℂ, sandstone.ρ, direction_from_angles(θ, 0.0))
    m = ti_modes(w, 0.0)
    map(i -> atan(w.group[1, i], w.group[3, i]), (qP = m.qP, qSV = m.qSV, SH = m.SH))
end
plot(
    rad2deg.(θq), [rad2deg.(getindex.(rays, k)) for k in (:qP, :qSV, :SH)]; label = ["qP" "qSV" "SH"], lw = 2,
    xlabel = "angle θ of the normal to the axis (°)", ylabel = "angle ψ of the ray to the axis (°)", legend = :topleft, size = (640, 380)
)
plot!(rad2deg.(θq), rad2deg.(θq); color = :gray, ls = :dash, label = "ψ = θ")
````

The code looks for the angles where ``\psi`` of qSV turns back:

````@example wave_surfaces
ψ = getindex.(rays, :qSV)
turns = [i for i in 2:(length(θq) - 1) if (ψ[i + 1] - ψ[i]) * (ψ[i] - ψ[i - 1]) < 0]
if isempty(turns)
    println("qSV: the ray angle increases with θ over the quadrant, and its wave surface has no cusp")
else
    for i in turns
        @printf "qSV: the ray angle turns back at θ = %.2f°, where ψ = %.2f°\n" rad2deg(θq[i]) rad2deg(ψ[i])
    end
    length(turns) ≥ 2 && @printf "the rays between ψ = %.2f° and %.2f° cross three qSV fronts\n" rad2deg.(extrema(ψ[turns]))...
end
````

## The limestone, without symmetry

The limestone has no plane of symmetry: its group velocities leave any plane
of directions, and its two shear waves are named qS1 and qS2 by their
velocities only. The phase velocities in two planes through
``\underline{e}_3``, of azimuths ``\phi = 0`` and ``\phi = \pi/2``:

````@example wave_surfaces
θl = range(0, π; length = 181)
plot(; xlabel = "angle θ to e₃ (°)", ylabel = "phase velocity (km/s)", legend = :outerright, size = (720, 380))
for (ϕ, style) in ((0.0, :solid), (π / 2, :dash))
    V = [phase_velocities(limestone.ℂ, limestone.ρ, direction_from_angles(θ, ϕ)) for θ in θl]
    for (m, name) in enumerate(("qP", "qS1", "qS2"))
        plot!(rad2deg.(θl), getindex.(V, m); color = m, ls = style, lw = 1.5, label = "$name, ϕ = $(round(rad2deg(ϕ)))°")
    end
end
current()
````

## The three sheets in space

The directions are sampled on a grid of the angles ``(\theta, \phi)``; a
quarter of the azimuths, ``0 < \phi < \pi/2``, is left out so that the inner
sheets show through.

````@example wave_surfaces
θg, ϕg = range(0, π; length = 31), range(π / 2, 2π; length = 37)
grid(medium) = [plane_waves(medium.ℂ, medium.ρ, direction_from_angles(θ, ϕ)) for θ in θg, ϕ in ϕg]
sandstone_grid, limestone_grid = grid(sandstone), grid(limestone)
nothing #hide
````

On every direction of both grids, ``\underline{g}\cdot\underline{n} = V`` and
``\underline{g}\cdot\underline{s} = 1``; the largest angle between the ray and
the normal, for each wave:

````@example wave_surfaces
for (name, waves) in (("sandstone", sandstone_grid), ("limestone", limestone_grid))
    gn = maximum(abs(dot(w.group[:, m], w.direction) - w.velocity[m]) for w in waves, m in 1:3)
    gs = maximum(abs(dot(w.group[:, m], slowness(w)[:, m]) - 1) for w in waves, m in 1:3)
    angle = [maximum(acosd(min(1.0, dot(w.group[:, m], w.direction) / norm(w.group[:, m]))) for w in waves) for m in 1:3]
    @printf "%-9s  max |g·n − V| = %.1e km/s   max |g·s − 1| = %.1e   largest angle: qP %.1f°, qS1 %.1f°, qS2 %.1f°\n" name gn gs angle...
end
````

That ``\underline{g}`` is normal to the slowness surface is checked with the
two tangents ``\partial\underline{s}/\partial\theta`` and
``\partial\underline{s}/\partial\phi``, by ForwardDiff through
`plane_waves`: the largest cosine between the group velocity and a tangent,
at a few directions of the limestone,

````@example wave_surfaces
s_of(m) = a -> (n = direction_from_angles(a[1], a[2]); n / phase_velocities(limestone.ℂ, limestone.ρ, n)[m])
maximum(((θ, ϕ), m) for (θ, ϕ) in ((0.4, 0.3), (1.1, 2.0), (2.0, 4.0)), m in 1:3) do ((θ, ϕ), m)
    t = ForwardDiff.jacobian(s_of(m), [θ, ϕ])
    g = plane_waves(limestone.ℂ, limestone.ρ, direction_from_angles(θ, ϕ)).group[:, m]
    maximum(abs(dot(g, c)) / (norm(g) * norm(c)) for c in eachcol(t))
end
````

The sheets are drawn by Plotly in the browser, in scenes that can be turned
and zoomed. The page receives the points of each sheet and one call to the
drawing function of the documentation theme; colors are named by theme
token, so that the scene reads on the light and the dark scheme alike:

````@example wave_surfaces
js_array(X) = "[" * join(("[" * join(round.(r; sigdigits = 4), ",") * "]" for r in eachrow(X)), ",") * "]"
function sheet(name, color, P)                     ## P: a matrix of points of space
    x, y, z = (js_array(map(q -> q[k], P)) for k in 1:3)
    c = "\"var(--tensnd-fig-$color)\""
    return """{"type":"surface","name":"$name","colorscale":[[0,$c],[1,$c]],"showscale":false,"showlegend":true,"opacity":0.85,"x":$x,"y":$y,"z":$z}"""
end
function scene(uid, title, unit, sheets)
    axis(label) = """{"title":{"text":"$label ($unit)"},"color":"var(--tensnd-fig-ink)"}"""
    layout = """{"title":{"text":"$title"},"height":480,"margin":{"l":0,"r":0,"t":40,"b":0},
    "paper_bgcolor":"rgba(0,0,0,0)","font":{"color":"var(--tensnd-fig-ink)"},"legend":{"x":0,"y":1},"uirevision":"$uid",
    "scene":{"aspectmode":"data","xaxis":$(axis("x₁")),"yaxis":$(axis("x₂")),"zaxis":$(axis("x₃"))}}"""
    return HTML(
        """
        <div id="$uid" style="width:100%;height:480px"></div>
        <script>
        window.tensndPlotly("$uid", [$(join(sheets, ","))], $layout);
        </script>
        """
    )
end
nothing #hide
````

The slowness surface of the limestone: qP inside, qS2 outside. The two shear
sheets touch at the directions where their velocities coincide.

````@example wave_surfaces
colors = ("q1", "q2", "q3")
scene(
    "limestone-slowness", "Limestone: slowness surface", "s/km",
    [sheet(name, colors[m], [slowness(w)[:, m] for w in limestone_grid]) for (m, name) in enumerate(("qP", "qS1", "qS2"))]
)
````

Its wave surface, the front of a pulse after one second:

````@example wave_surfaces
scene(
    "limestone-wave", "Limestone: wave surface", "km/s",
    [sheet(name, colors[m], [w.group[:, m] for w in limestone_grid]) for (m, name) in enumerate(("qP", "qS1", "qS2"))]
)
````

The wave surface of the sandstone, with its modes named by their
polarizations as in the plane above; it has the symmetry of revolution about
``\underline{e}_3``, and a fold of qSV shows as a ring.

````@example wave_surfaces
sandstone_modes = [ti_modes(sandstone_grid[i, j], ϕg[j]) for i in eachindex(θg), j in eachindex(ϕg)]
scene(
    "sandstone-wave", "Sandstone: wave surface", "km/s",
    [sheet(string(k), colors[i], map((w, m) -> w.group[:, m[k]], sandstone_grid, sandstone_modes)) for (i, k) in enumerate((:qP, :qSV, :SH))]
)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

