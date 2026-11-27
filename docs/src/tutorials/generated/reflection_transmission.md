```@meta
EditURL = "../../../../scripts/28_reflection_transmission.jl"
```

# [Reflection and transmission at an interface](@id tut-reflection)

!!! tip "See also"
    **Theory**: [Waves at an interface](@ref th-interface-waves) · [Plane waves](@ref th-plane-waves) · [Thomsen's parameters](@ref th-thomsen)

    **Manual**: [Waves at an interface](@ref man-interface-waves) · [Plane waves](@ref man-plane-waves)

    **Tutorials**: [Slowness and wave surfaces](@ref tut-wave-surfaces)

A P wave that meets the plane interface between two solids gives rise to
reflected and transmitted P and S waves. This tutorial computes their
directions, their amplitudes and the energy they carry against the angle of
incidence, for isotropic media, for a transversely isotropic one, and across
an interface that slips; then it differentiates the reflection coefficient
with respect to the properties of the second medium and compares the result
with the classical linearized formula.

````@example reflection_transmission
using TensND
using LinearAlgebra
using ForwardDiff
using Printf
using Plots
gr()
nothing #hide
````

## Two half-spaces

An isotropic medium is given by the velocities ``\alpha`` of its P waves and
``\beta`` of its S waves (km/s) and by its density ``\rho`` (g/cm³): its shear
modulus is ``\rho\beta^2`` and its bulk modulus ``\rho(\alpha^2 - 4\beta^2/3)``,
in GPa. The two media are round values of the orders of magnitude of a shale
above a sandstone, an illustration rather than measured rocks. The interface
is the plane ``x_3 = 0`` of normal ``\underline{\nu} = \underline{e}_3``, the
shale on the side ``x_3 < 0``.

````@example reflection_transmission
isotropic(α, β, ρ) = (ℂ = iso_stiffness(ρ * (α^2 - 4β^2 / 3), ρ * β^2), ρ = ρ, α = α, β = β)
shale = isotropic(2.9, 1.4, 2.3)
sandstone = isotropic(3.5, 2.1, 2.4)
ν = [0.0, 0.0, 1.0]
nothing #hide
````

## Snell–Descartes on the slownesses

All the waves excited by an incident one share the tangential part
``\underline{s}_t`` of its slowness: for an incident P wave at the angle
``\theta`` in the plane ``(\underline{e}_1, \underline{e}_3)``,
``\underline{s}_t = p\,\underline{e}_1`` with the horizontal slowness
``p = \sin\theta/\alpha_1``. [`interface_waves`](@ref) returns the six waves of
a medium for that ``p``, three going towards ``+\underline{\nu}`` and three
towards ``-\underline{\nu}``, each by increasing normal slowness ``|\xi|``; the
incident P wave is the first forward wave of the shale.

````@example reflection_transmission
incident_P(m, p) = first(interface_waves(m.ℂ, m.ρ, ν, [p, 0.0, 0.0]))[1]
θ = deg2rad(30)
p = sin(θ) / shale.α
r = reflection_transmission(shale.ℂ, shale.ρ, sandstone.ℂ, sandstone.ρ, ν, incident_P(shale, p))
nothing #hide
````

The reflected waves are the backward waves of the shale and the transmitted
ones the forward waves of the sandstone. Each slowness gives the velocity
``1/|\underline{s}|`` of its wave and its angle to the normal, whose sine over
the velocity is ``p`` for every one of them; the two S waves of an isotropic
medium share their slowness.

````@example reflection_transmission
@printf "incident P: p = %.5f s/km\n" p
for (side, waves) in (("reflected", r.reflected), ("transmitted", r.transmitted)), (w, name) in zip(waves, ("P", "S", "S"))
    s = real.(w.slowness)
    V, angle = 1 / norm(s), acos(abs(s[3]) / norm(s))
    @printf "%-11s %s  V = %.3f km/s  angle to the normal %6.2f°  sin(angle)/V = %.5f s/km\n" side name V rad2deg(angle) sin(angle) / V
end
````

## Amplitudes and energies against the angle of incidence

[`reflection_transmission`](@ref) returns the displacement amplitudes `R` and
`T` of the six waves, on unit polarizations whose sign is arbitrary, and the
fractions `energy` of the incident energy flux that each carries across the
interface. The two S waves of a medium are read together, by the norm of
their total displacement and the sum of their energies. The P-P coefficient
is given a sign by reading the displacement of each P wave along its
direction of propagation ``\underline{s}/|\underline{s}|``:

````@example reflection_transmission
direction(w) = (s = real.(w.slowness); s / norm(s))
along(w) = sum(direction(w) .* w.polarization)
signed_PP(r, incident) = r.R[1] * along(r.reflected[1]) / along(incident)
function coefficients(m₁, m₂, p; kw...)
    incident = incident_P(m₁, p)
    r = reflection_transmission(m₁.ℂ, m₁.ρ, m₂.ℂ, m₂.ρ, ν, incident; kw...)
    S(a, waves) = norm(a[2] * waves[2].polarization + a[3] * waves[3].polarization)
    eR, eT = r.energy
    return (
        PP = signed_PP(r, incident), PS = S(r.R, r.reflected), TP = abs(r.T[1]), TS = S(r.T, r.transmitted),
        energy = (eR[1], eR[2] + eR[3], eT[1], eT[2] + eT[3]), homogeneous = r.transmitted[1].homogeneous,
    )
end
nothing #hide
````

At normal incidence the sign is that of the impedance contrast,
``R = (Z_2 - Z_1)/(Z_2 + Z_1)`` with ``Z = \rho\alpha``:

````@example reflection_transmission
Z₁, Z₂ = shale.ρ * shale.α, sandstone.ρ * sandstone.α
coefficients(shale, sandstone, 0.0).PP, (Z₂ - Z₁) / (Z₂ + Z₁)
````

The transmitted P wave is faster than the incident one: beyond a critical
angle its normal slowness is complex, the wave is evanescent and carries no
energy away. The code locates that angle on the flag `homogeneous`, then by
bisection, and compares it with ``\arcsin(\alpha_1/\alpha_2)``:

````@example reflection_transmission
θs = deg2rad.(0:0.5:89.5)
results = [coefficients(shale, sandstone, sin(θ) / shale.α) for θ in θs]
function critical_angle(a, b)                       ## homogeneous at a, evanescent at b
    for _ in 1:12
        c = (a + b) / 2
        coefficients(shale, sandstone, sin(c) / shale.α).homogeneous ? (a = c) : (b = c)
    end
    return (a + b) / 2
end
k = findfirst(c -> !c.homogeneous, results)
θc = critical_angle(θs[k - 1], θs[k])
rad2deg(θc), asind(shale.α / sandstone.α)
````

The amplitudes, with the P-P coefficient by its modulus past the critical
angle, and the energy fractions with their sum:

````@example reflection_transmission
deg = rad2deg.(θs)
amplitudes = plot(
    deg, [abs.(getfield.(results, :PP)) getfield.(results, :PS) getfield.(results, :TP) getfield.(results, :TS)];
    label = ["reflected P" "reflected S" "transmitted P" "transmitted S"], lw = 2, xlabel = "angle of incidence (°)", ylabel = "amplitude", legend = :topleft
)
vline!(amplitudes, [rad2deg(θc)]; color = :gray, ls = :dash, label = "critical angle")
E = reduce(hcat, [collect(c.energy) for c in results])'
energies = plot(
    deg, [E sum(E; dims = 2)]; label = ["reflected P" "reflected S" "transmitted P" "transmitted S" "sum"],
    lw = 2, ls = [:solid :solid :solid :solid :dash], color = [1 2 3 4 :black], xlabel = "angle of incidence (°)", ylabel = "energy fraction", legend = :left
)
vline!(energies, [rad2deg(θc)]; color = :gray, ls = :dash, label = "")
plot(amplitudes, energies; layout = (1, 2), size = (980, 380), left_margin = 3Plots.mm, bottom_margin = 4Plots.mm)
````

The energy balance, over all the angles:

````@example reflection_transmission
maximum(abs, sum(E; dims = 2) .- 1)
````

## A transversely isotropic shale

A shale is often transversely isotropic about the vertical. Built from
Thomsen's parameters ([Thomsen's parameters](@ref th-thomsen)) with the same
velocities along the axis as the isotropic shale, round values again,

````@example reflection_transmission
vti = (ℂ = ti_from_constants(2.9, 1.4, 0.2, 0.1, 0.15, ν; form = :thomsen, density = 2.3), ρ = 2.3)
thomsen_parameters(vti.ℂ, vti.ρ)
````

the incident qP wave is given by the angle ``\theta`` of its normal: its
horizontal slowness is ``p = \sin\theta/V_{qP}(\theta)``, with the phase
velocity of that direction, and its ray leaves the normal. The signed P-P
coefficient below the critical angle, where it is real, against that of the
isotropic shale:

````@example reflection_transmission
p_vti(θ) = sin(θ) / phase_velocities(vti.ℂ, vti.ρ, direction_from_angles(θ, 0.0))[1]
results_vti = [coefficients(vti, sandstone, p_vti(θ)) for θ in θs]
precritical(c) = c.homogeneous ? real(c.PP) : NaN
plot(
    deg, [precritical.(results) precritical.(results_vti)]; label = ["isotropic shale" "transversely isotropic shale"], lw = 2,
    xlabel = "angle of incidence of the normal (°)", ylabel = "signed P-P coefficient", legend = :topleft, size = (640, 380)
)
````

The energy balance holds for the anisotropic medium as well:

````@example reflection_transmission
maximum(c -> abs(sum(c.energy) - 1), results_vti)
````

## An interface that slips

A fracture or a weak layer is modeled by a linear slip interface
[schoenberg1980](@cite): the traction is continuous and the displacement
jumps by ``\boldsymbol{Z}`` times the traction, ``\boldsymbol{Z}`` the
compliance of the interface, here ``Z_N`` normal to it and ``Z_T`` along it.
The keyword `slip` takes ``\omega\boldsymbol{Z}``, so that the coefficients
depend on the frequency ``f = \omega/2\pi``. With lengths in km and moduli in
GPa, a compliance in m/Pa is multiplied by ``10^6``:

````@example reflection_transmission
Zₙ, Zₜ = 1.0e-10 * 1.0e6, 2.0e-10 * 1.0e6
slip(f) = 2π * f * Diagonal([Zₜ, Zₜ, Zₙ])
fs = range(0, 500; length = 101)
R_slip(m₁, m₂, f) = abs(coefficients(m₁, m₂, 0.0; slip = slip(f)).PP)
nothing #hide
````

At normal incidence, between two identical media, only the slip reflects,
and ``|R| = (\kappa/2)/\sqrt{1 + (\kappa/2)^2}`` with
``\kappa = \omega Z_N\rho\alpha`` [schoenberg1980](@cite). The modulus of the
P-P coefficient against the frequency:

````@example reflection_transmission
κ(f) = 2π * f * Zₙ * sandstone.ρ * sandstone.α
closed(f) = (κ(f) / 2) / sqrt(1 + (κ(f) / 2)^2)
plot(
    fs, [R_slip.(Ref(shale), Ref(sandstone), fs) R_slip.(Ref(sandstone), Ref(sandstone), fs)];
    label = ["shale on sandstone" "sandstone on sandstone"], lw = 2, xlabel = "frequency (Hz)", ylabel = "|R| of the P wave, normal incidence",
    legend = :topleft, size = (640, 380)
)
hline!([abs((Z₂ - Z₁) / (Z₂ + Z₁))]; color = 1, ls = :dot, label = "shale on sandstone, welded")
scatter!(fs[1:10:end], closed.(fs[1:10:end]); color = 2, label = "closed form")
````

The largest difference between the computed coefficient and the closed form:

````@example reflection_transmission
maximum(f -> abs(R_slip(sandstone, sandstone, f) - closed(f)), fs)
````

## The coefficient of a small contrast

For two media that differ little, the P-P coefficient is linear in the
contrasts. The functions of the package take dual numbers, so its gradient
with respect to the density and the velocities ``(\rho_2, \alpha_2, \beta_2)``
of the second medium, at zero contrast (the second medium a copy of the
first), is one call to ForwardDiff:

````@example reflection_transmission
function R_PP(x, θ)                                  ## x = (ρ₂, α₂, β₂)
    m₂ = isotropic(x[2], x[3], x[1])
    incident = incident_P(shale, sin(θ) / shale.α)
    r = reflection_transmission(shale.ℂ, shale.ρ, m₂.ℂ, m₂.ρ, ν, incident)
    return real(signed_PP(r, incident))
end
x₀ = [shale.ρ, shale.α, shale.β]
ForwardDiff.gradient(x -> R_PP(x, deg2rad(20)), x₀)
````

The linearized formula of [aki2002](@citet), for the same convention of
sign, reads with ``p`` the horizontal slowness

```math
R \approx \tfrac12\bigl(1 - 4p^2\beta^2\bigr)\frac{\Delta\rho}{\rho}
  + \frac{\Delta\alpha}{2\alpha\cos^2\theta} - 4p^2\beta^2\,\frac{\Delta\beta}{\beta} ,
```

whose derivatives at zero contrast are the three coefficients of the
contrasts. Both, made dimensionless by ``(\rho, \alpha, \beta)``, against the
angle of incidence:

````@example reflection_transmission
function linearized(θ)
    ρ, α, β = x₀
    q = sin(θ) / α
    return [(1 - 4q^2 * β^2) / (2ρ), 1 / (2α * cos(θ)^2), -4q^2 * β]
end
θd = deg2rad.(0:2:60)
by_duals = [ForwardDiff.gradient(x -> R_PP(x, θ), x₀) .* x₀ for θ in θd]
by_formula = [linearized(θ) .* x₀ for θ in θd]
plot(
    rad2deg.(θd), reduce(hcat, by_duals)'; label = ["ρ ∂R/∂ρ₂" "α ∂R/∂α₂" "β ∂R/∂β₂"], lw = 2,
    xlabel = "angle of incidence (°)", ylabel = "sensitivity of R", legend = :topleft, size = (640, 380)
)
scatter!(rad2deg.(θd), reduce(hcat, by_formula)'; color = [1 2 3], label = "")
````

The markers are the formula; the largest difference between the two:

````@example reflection_transmission
maximum(maximum(abs, a - b) for (a, b) in zip(by_duals, by_formula))
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

