```@meta
EditURL = "../../../../scripts/26_thomsen_parameters.jl"
```

# [Thomsen's parameters and weak anisotropy](@id tut-thomsen)

!!! tip "See also"
    **Theory**: [Plane waves](@ref th-plane-waves) · [Thomsen's parameters](@ref th-thomsen)

    **Manual**: [Plane waves](@ref man-plane-waves)

    **Tutorials**: [Weak anisotropy of any symmetry](@ref tut-generalized-thomsen) · [Stress-induced anisotropy](@ref tut-stress-anisotropy)

!!! note "Two-index components"
    ``C_{11}``, ``C_{44}``, ``C_{66}``, … are the two-index (Voigt) components
    ``C_{IJ}`` of the stiffness, ``I = 1, \dots, 6`` standing for the pairs 11,
    22, 33, 23, 13, 12 (``C_{44} = C_{2323}``, ``C_{66} = C_{1212}``), not the
    entries of the Kelvin–Mandel matrix, which carry factors ``\sqrt2`` and ``2``.
    Both matrices are written in full in [Contrast with Voigt](@ref th-voigt).

A transversely isotropic medium has five elastic constants. [thomsen1986](@citet)
replaces them by the two velocities along the axis, ``V_{P0}`` and
``V_{S0}``, and three dimensionless numbers ``\varepsilon``, ``\delta`` and
``\gamma`` that vanish for an isotropic medium and are small for most rocks.
This tutorial builds media from the parameters of his Table 1, compares their
exact velocities with the first-order ones, recovers the first-order formulas
from the exact ones with symbols, and shows why ``\delta`` governs the
directions close to the axis.

````@example thomsen_parameters
using TensND
using LinearAlgebra
using ForwardDiff
using JSON
using Printf
using Plots
using SymPyPythonCall: symbols, N
gr()
nothing #hide
````

## The parameters of measured rocks and crystals

Each row of Table 1 of [thomsen1986](@citet) gives ``V_{P0}`` and ``V_{S0}``
(m/s), ``\varepsilon``, ``\delta^*``, ``\delta``, ``\gamma`` and the density
(g/cm³). A selection of rows is stored with its source under
`data/literature`.

````@example thomsen_parameters
data = JSON.parsefile(joinpath(pkgdir(TensND), "data", "literature", "thomsenGEO1986.json"))
rows = data["tables"]["table_1"]
[(r["sample"], r["epsilon"], r["delta"], r["gamma"]) for r in rows]

Main.show_literature("thomsenGEO1986") # hide
````

With the velocities in km/s and the density in g/cm³, the stiffness comes
in GPa. `ti_from_constants(V_P0, V_S0, ε, δ, γ, axis; form = :thomsen,
density = ρ)` returns the transversely isotropic stiffness ``\mathbb{C}`` of
axis ``\underline{e}_3``:

````@example thomsen_parameters
e₃ = [0.0, 0.0, 1.0]
function medium(r)
    ρ = r["rho (g/cm3)"]
    ℂ = ti_from_constants(r["V_P0 (m/s)"] / 1000, r["V_S0 (m/s)"] / 1000, r["epsilon"], r["delta"], r["gamma"], e₃; form = :thomsen, density = ρ)
    return ℂ, ρ
end
table_row(name) = only(filter(r -> r["sample"] == name, rows))
weak_rock, strong_rock = "Mesaverde (6455.1) immature sandstone", "Mesaverde (5501) clayshale"
ℂ, ρ = medium(table_row(weak_rock))
round.(voigt_stiffness(ℂ); digits = 2)
````

The matrix is that of the two-index notation (`voigt_stiffness`), in GPa.
`thomsen_parameters` reads the parameters back from the stiffness:

````@example thomsen_parameters
thomsen_parameters(ℂ, ρ)
````

## Exact and first-order velocities

Along a direction ``\underline{n}`` at the angle ``\theta`` to the axis, the
three waves are a quasi-longitudinal one (qP), a quasi-transverse one
polarized in the plane of ``\underline{n}`` and the axis (qSV), and a
transverse one normal to that plane (SH). `wave_modes` gives their exact
velocities, `weak_wave_modes` the first-order ones

```math
\begin{gathered}
V_{P} \approx V_{P0}\,(1 + \delta\sin^2\theta\cos^2\theta + \varepsilon\sin^4\theta),
\\
V_{SV} \approx V_{S0}\Bigl[1 + \frac{V_{P0}^2}{V_{S0}^2}(\varepsilon - \delta)\sin^2\theta\cos^2\theta\Bigr],
\\
V_{SH} \approx V_{S0}\,(1 + \gamma\sin^2\theta) .
\end{gathered}
```

````@example thomsen_parameters
θs = range(0, π / 2; length = 91)
modes = (:qP, :qSV, :SH)
exact_velocity(ℂ, ρ, θ, k) = wave_modes(ℂ, ρ, direction_from_angles(θ, 0.0))[k][1]
weak_velocity(ℂ, ρ, θ, k) = weak_wave_modes(ℂ, ρ, direction_from_angles(θ, 0.0))[k]

function velocity_panel(name)
    ℂ, ρ = medium(table_row(name))
    V = [exact_velocity(ℂ, ρ, θ, k) for θ in θs, k in modes]
    W = [weak_velocity(ℂ, ρ, θ, k) for θ in θs, k in modes]
    p = plot(; title = name, titlefontsize = 10, xlabel = "angle to the axis θ (°)", ylabel = "phase velocity (km/s)", xticks = 0:15:90, legend = :right)
    plot!(p, rad2deg.(θs), V; lw = 2, color = [1 2 3], label = ["qP" "qSV" "SH"])
    plot!(p, rad2deg.(θs), W; lw = 2, ls = :dash, color = [1 2 3], label = ["first order" "" ""])
    return p
end
plot(velocity_panel(weak_rock), velocity_panel(strong_rock); layout = (1, 2), size = (980, 380), left_margin = 4Plots.mm, bottom_margin = 4Plots.mm)
````

The largest relative error of the first-order velocities over the angles:

````@example thomsen_parameters
for name in (weak_rock, strong_rock)
    ℂ, ρ = medium(table_row(name))
    err = [abs(weak_velocity(ℂ, ρ, θ, k) / exact_velocity(ℂ, ρ, θ, k) - 1) for θ in θs, k in modes]
    @printf "%-38s  qP %.2f %%   qSV %.2f %%   SH %.2f %%\n" name (100 .* maximum(err; dims = 1))...
end
````

## The first-order formulas from the exact ones

The three parameters are scaled by ``\tau`` and the exact velocities expanded
to first order in ``\tau`` with SymPy. ``V_{P0} = V_{S0} + \Delta V`` with
``\Delta V > 0`` tells SymPy that ``V_{P0} > V_{S0}``, so that it takes
``\sqrt{(C_{33} - C_{44})^2} = C_{33} - C_{44}``.

````@example thomsen_parameters
V_S0, ΔV, ρₛ = symbols("V_S0 ΔV ρ", positive = true)
εₛ, δₛ, γₛ, τ, θ = symbols("ε δ γ τ θ", real = true)
ℂτ = ti_from_constants(V_S0 + ΔV, V_S0, τ * εₛ, τ * δₛ, τ * γₛ, [0, 0, 1]; form = :thomsen, density = ρₛ)
n = [sin(θ), 0, cos(θ)]
exact_s = wave_modes(ℂτ, ρₛ, n)
weak_s = weak_wave_modes(ℂτ, ρₛ, n)
first_order(V) = tsubs(V, τ => 0) + τ * tsubs(tdiff(V, τ), τ => 0)
nothing #hide
````

For SH, SymPy reduces the difference between the series and Thomsen's formula
to zero:

````@example thomsen_parameters
tsimplify(first_order(exact_s.SH[1]) - weak_s.SH)
````

For qP and qSV the difference is evaluated at the parameters of each row of
the table and at three angles; it vanishes to rounding:

````@example thomsen_parameters
function series_gap(k)
    Δ = first_order(exact_s[k][1]) - weak_s[k]
    gaps = [
        N(
            tsubs(
                Δ, θ => ϑ, τ => 1, V_S0 => r["V_S0 (m/s)"] / 1000,
                ΔV => (r["V_P0 (m/s)"] - r["V_S0 (m/s)"]) / 1000, εₛ => r["epsilon"], δₛ => r["delta"], γₛ => r["gamma"],
                ρₛ => r["rho (g/cm3)"]
            )
        ) for r in rows, ϑ in (0.3, 0.8, 1.3)
    ]
    return maximum(abs ∘ float, gaps)
end
(qP = series_gap(:qP), qSV = series_gap(:qSV))
````

## Near the axis: δ, not ε

Close to the axis, ``\sin^2\theta\cos^2\theta \approx \theta^2`` while
``\sin^4\theta \approx \theta^4``: the P velocity departs from ``V_{P0}`` at
the rate set by ``\delta``, and ``\varepsilon`` takes over towards the
horizontal, where ``V_P = V_{P0}\sqrt{1 + 2\varepsilon}``. Four rocks with
different signs and sizes of ``\delta`` and ``\varepsilon``:

````@example thomsen_parameters
near = ("Taylor sandstone", "Mesaverde (4903) mudshale", "Dog Creek shale", "Berea sandstone")
p = plot(; xlabel = "angle to the axis θ (°)", ylabel = "V_P(θ) / V_P0", xticks = 0:15:90, legend = :topleft, size = (640, 380), left_margin = 4Plots.mm, bottom_margin = 4Plots.mm)
for name in near
    r = table_row(name)
    ℂ, ρ = medium(r)
    V₀ = r["V_P0 (m/s)"] / 1000
    plot!(p, rad2deg.(θs), [exact_velocity(ℂ, ρ, θ, :qP) / V₀ for θ in θs]; lw = 2, label = @sprintf("%s (ε = %.3f, δ = %.3f)", name, r["epsilon"], r["delta"]))
end
p
````

A reflection survey over a horizontal reflector measures the moveout
velocity, set by the curvature of the phase velocity at the axis,
``V_{\mathrm{nmo}}^2 = V(0)^2\,[1 + V''(0)/V(0)]``. For P it is
``V_{P0}\sqrt{1 + 2\delta}``. `squared_moveout_velocities` gives the squares
for the three modes; the qSV one is negative when ``\varepsilon - \delta``
is large and negative, and the qSV wave has then no moveout velocity:

````@example thomsen_parameters
ℂ, ρ = medium(table_row(weak_rock))
squared_moveout_velocities(ℂ, ρ)
````

The same curvature, taken on the exact P velocity by two nested
ForwardDiff derivatives, for every row; ``\eta = (\varepsilon - \delta)/(1 + 2\delta)``
(`anellipticity`) measures how far the P wave is from an elliptical one, for
which ``V_{\mathrm{nmo}}`` and the horizontal velocity coincide:

````@example thomsen_parameters
function moveout_from_curvature(V)
    V₀ = V(0.0)
    V″ = ForwardDiff.derivative(θ -> ForwardDiff.derivative(V, θ), 0.0)
    return sqrt(V₀^2 * (1 + V″ / V₀))
end
@printf "%-42s %7s %7s %7s %9s %9s\n" "sample" "ε" "δ" "η" "V_nmo" "V(90°)"
for r in rows
    ℂ, ρ = medium(r)
    V = θ -> exact_velocity(ℂ, ρ, θ, :qP)
    @printf "%-42s %7.3f %7.3f %7.3f %9.3f %9.3f\n" r["sample"] r["epsilon"] r["delta"] anellipticity(ℂ) moveout_from_curvature(V) V(π / 2)
end
````

The velocities are in km/s. A survey that reads ``V_{\mathrm{nmo}}`` close
to ``V_{P0}`` learns ``\delta``, and nothing of ``\varepsilon``.

## Two definitions of δ

Thomsen defines a first parameter ``\delta^*`` before replacing it by
``\delta``; both are functions of the same constants, and
``\delta^* = (1 - V_{S0}^2/V_{P0}^2)(2\delta - \varepsilon)`` exactly
([Thomsen's parameters](@ref th-thomsen)). The rows where the printed
``\delta^*`` and the value of this relation differ by more than the rounding
of the table, with the ``\delta`` that the printed ``\delta^*`` would give:

````@example thomsen_parameters
δ_star(r) = (1 - (r["V_S0 (m/s)"] / r["V_P0 (m/s)"])^2) * (2 * r["delta"] - r["epsilon"])
δ_from_star(r) = (r["delta_star"] / (1 - (r["V_S0 (m/s)"] / r["V_P0 (m/s)"])^2) + r["epsilon"]) / 2
@printf "%-42s %9s %9s %9s %9s\n" "sample" "δ* print" "δ* rel." "δ print" "δ rel."
for r in rows
    abs(δ_star(r) - r["delta_star"]) > 0.01 || continue
    @printf "%-42s %9.3f %9.3f %9.3f %9.3f\n" r["sample"] r["delta_star"] δ_star(r) r["delta"] δ_from_star(r)
end
````

These rows are crystals. The computations above use the printed ``\delta``,
the parameter the first-order velocities are written with; the table does not
say which of the two printed values the constants of the source support.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

