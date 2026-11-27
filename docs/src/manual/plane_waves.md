# [Plane waves](@id man-plane-waves)

!!! tip "See also"
    **Theory**: [Plane waves](@ref th-plane-waves) · [Thomsen's parameters](@ref th-thomsen) · [Weak anisotropy of any symmetry](@ref th-weak-anisotropy) · [Shear-wave splitting](@ref th-splitting)

    **Manual**: [Waves at an interface](@ref man-interface-waves) · [Third-order elasticity](@ref man-third-order) · [Elastic moduli](@ref man-elastic-moduli) · [Symbolic tensors of a material symmetry](@ref man-symbolic-materials)

    **Tutorials**: [Slowness and wave surfaces](@ref tut-wave-surfaces) · [Thomsen's parameters](@ref tut-thomsen) · [Generalized Thomsen parameters](@ref tut-generalized-thomsen)

A plane wave ``\underline{u} = \underline{U}\,f(\underline{n}\cdot\underline{x} - Vt)``
of unit direction ``\underline{n}`` travels in a homogeneous medium of stiffness
``\mathbb{C}`` and mass density ``\rho`` when its phase velocity ``V`` and its
polarization ``\underline{U}`` solve the Christoffel problem

```math
\bigl(\boldsymbol{\Gamma}(\underline{n}) - \rho V^2\,\boldsymbol{1}\bigr)\cdot\underline{U} = \underline{0},
\qquad \boldsymbol{\Gamma}(\underline{n}) = \underline{n}\cdot\mathbb{C}\cdot\underline{n} .
```

Three waves go along each direction: the fastest, quasi-longitudinal (qP),
and two quasi-transverse ones (qS1 faster than qS2). This page gives the
syntax; the derivations are on [Plane waves](@ref th-plane-waves).

| to … | write |
| :--- | :--- |
| solve the Christoffel problem | `plane_waves(ℂ, ρ, n)` |
| read the phase velocities only | `phase_velocities(ℂ, ρ, n)` |
| read the slowness vectors ``\underline{n}/V`` | `slowness(plane_waves(ℂ, ρ, n))` |
| get closed forms, isotropic or transversely isotropic | `wave_modes(ℂ, ρ, n)` |
| get the characteristic polynomial of ``\boldsymbol{\Gamma}/\rho`` | `christoffel_invariants(ℂ, ρ, n)` |
| read Thomsen's parameters | `thomsen_parameters(ℂ, ρ)`, `ti_constants(ℂ; form = :thomsen, density = ρ)` |
| build a stiffness from them | `ti_from_constants(V_P0, V_S0, ε, δ, γ, axis; form = :thomsen, density = ρ)` |
| velocities to first order in Thomsen's parameters | `weak_wave_modes(ℂ, ρ, n)`, `anellipticity(ℂ)` |
| squared moveout velocities of a horizontal reflector | `squared_moveout_velocities(ℂ, ρ)` |
| velocities to first order, any symmetry | `weak_plane_waves(ℂ, ρ, n)`, `generalized_thomsen_parameters(ℂ)` |
| measure the splitting of the shear waves | `shear_splitting(ℂ, ρ, n)` |

## The three waves of a direction

[`plane_waves`](@ref) takes the stiffness, the density and a direction, any
nonzero vector of three components, normalized on the way. With the stiffness
in GPa and the density in g/cm³, the velocities come out in km/s. The example is
a transversely isotropic stiffness of axis ``\underline{e}_3``, given by its
components ``(C_{1111}, C_{1122}, C_{1133}, C_{3333}, C_{2323})``
([Parametrizations and conversions](@ref man-parametrizations)):

```@example waves
using TensND, LinearAlgebra

ℂ = ti_from_constants(52.9, 15.7, 21.8, 47.8, 16.4, [0.0, 0.0, 1.0])
ρ = 2.45
n = [1.0, 0.0, 1.0]
w = plane_waves(ℂ, ρ, n)
```

The result is a [`PlaneWaves`](@ref): its fields `direction` (unit),
`velocity` (qP, qS1, qS2), `polarization` and `group` hold one wave per column.

```@example waves
w.polarization
```

[`phase_velocities`](@ref) returns the field `velocity` alone, and
[`slowness`](@ref) the vectors ``\underline{s} = \underline{n}/V``, one per
column:

```@example waves
phase_velocities(ℂ, ρ, n) == w.velocity, slowness(w)
```

## The group velocity

The energy of a wave travels at its group velocity
``\underline{g} = (\underline{U}\cdot\mathbb{C}\cdot\underline{U})\cdot\underline{n}/(\rho V)``,
which is not along ``\underline{n}`` in an anisotropic medium. Its projection on
``\underline{n}`` is the phase velocity:

```@example waves
transpose(w.group) * w.direction - w.velocity
```

The angle between the ray and the wave normal, in degrees, for each of the
three waves:

```@example waves
[acosd(dot(g, w.direction) / norm(g)) for g in eachcol(w.group)]
```

## Closed forms

[`wave_modes`](@ref) gives the waves of an isotropic or of a transversely
isotropic medium in closed form, under the name of their mode: `P` and `S`, or
`qP`, `qSV` (polarized in the plane of the axis and of ``\underline{n}``) and
`SH` (normal to it). Each mode is a pair `(velocity, polarization)`; no
velocity is compared with another, so that a mode keeps its name where two
velocities cross.

```@example waves
modes = wave_modes(ℂ, ρ, n)
```

```@example waves
map(first, modes)
```

The polarizations of qP and qSV are given by their components along the
transverse direction and along the axis, in that order, for a direction on the
side of the axis, ``\underline{n}\cdot\underline{a} \geq 0``; the direction
``-\underline{n}`` carries the same waves. An isotropic medium has one shear
velocity for every polarization normal to ``\underline{n}``, and the
polarization of `S` is `nothing`:

```@example waves
wave_modes(iso_stiffness(20.0, 12.0), 2.5, [0.0, 0.0, 1.0])
```

No closed form needs a comparison, so the constants may be symbols:

```@example waves
using SymPyPythonCall

C₁₁, C₁₂, C₁₃, C₃₃, C₄₄, ϱ = symbols("C_11 C_12 C_13 C_33 C_44 rho", positive = true)
ϑ = symbols("theta", real = true)
ℂₛ = ti_from_constants(C₁₁, C₁₂, C₁₃, C₃₃, C₄₄, [0, 0, 1])
simplify(wave_modes(ℂₛ, ϱ, [sin(ϑ), 0, cos(ϑ)]).SH[1]^2)
```

[`christoffel_invariants`](@ref) returns the trace ``I_1``, the trace of the
adjugate ``I_2`` and the determinant ``I_3`` of
``\boldsymbol{\Gamma}(\underline{n})/\rho``, polynomial in the stiffness and in
the components of ``\underline{n}``. For a unit ``\underline{n}`` the squared
velocities are the roots of ``v^3 - I_1v^2 + I_2v - I_3``:

```@example waves
I₁, I₂, I₃ = christoffel_invariants(ℂ, ρ, w.direction)
[v^3 - I₁ * v^2 + I₂ * v - I₃ for v in w.velocity .^ 2]
```

```@example waves
k, μ = symbols("k mu", positive = true)
simplify.(christoffel_invariants(iso_stiffness(k, μ), ϱ, [0, 0, 1]))
```

## Thomsen's parameters

A transversely isotropic medium is often described by the velocities of its
P and S waves along the axis, ``V_{P0}`` and ``V_{S0}``, and three
dimensionless numbers ``\varepsilon``, ``\delta``, ``\gamma`` that vanish for an
isotropic medium [thomsen1986](@cite) ([Thomsen's parameters](@ref th-thomsen)).
[`thomsen_parameters`](@ref) reads them, by their exact definitions:

```@example waves
p = thomsen_parameters(ℂ, ρ)
```

On `ℂₛ`, the symbolic stiffness of the previous section, the definitions
themselves come out:

```@example waves
map(simplify, thomsen_parameters(ℂₛ, ϱ))
```

They are also the form `:thomsen` of [`ti_constants`](@ref) and
[`ti_from_constants`](@ref), with the density as a keyword; the converse gives
the stiffness back:

```@example waves
ℂₜ = ti_from_constants(values(p)..., [0.0, 0.0, 1.0]; form = :thomsen, density = ρ)
ti_constants(ℂ; form = :thomsen, density = ρ) == values(p), ti_constants(ℂₜ) .- ti_constants(ℂ)
```

``\delta`` fixes ``(C_{1133} + C_{2323})^2`` only: the converse keeps
``C_{1133} + C_{2323} > 0``. Without `density`, the form `:thomsen` throws an
`ArgumentError`.

The anellipticity ``\eta = (\varepsilon - \delta)/(1 + 2\delta)`` vanishes when
the slowness curve of the P wave is an ellipse; for `ℂ` and for `ℂₛ`:

```@example waves
anellipticity(ℂ), simplify(anellipticity(ℂₛ))
```

[`weak_wave_modes`](@ref) gives the velocities to first order in
``(\varepsilon, \delta, \gamma)``, to compare with the exact ones:

```@example waves
weak_wave_modes(ℂ, ρ, n), map(first, wave_modes(ℂ, ρ, n))
```

With Thomsen's parameters as symbols, the stiffness of the form `:thomsen`
gives back the first-order formulas of
[Thomsen's parameters](@ref th-thomsen), at the angle ``\vartheta`` to the axis:

```@example waves
V_P0, V_S0 = symbols("V_P0 V_S0", positive = true)
εₜ, δₜ, γₜ = symbols("epsilon delta gamma", real = true)
ℂᵥ = ti_from_constants(V_P0, V_S0, εₜ, δₜ, γₜ, [0, 0, 1]; form = :thomsen, density = ϱ)
wₜ = weak_wave_modes(ℂᵥ, ϱ, [sin(ϑ), 0, cos(ϑ)])
```

The parameters of ``\mathbb{C}`` and the angle ``\pi/4`` of ``\underline{n}``
substituted in these formulas give the numeric result above:

```@example waves
at = (V_P0 => p.V_P0, V_S0 => p.V_S0, εₜ => p.ε, δₜ => p.δ, γₜ => p.γ, ϱ => ρ, ϑ => PI / 4)
[N(subs(v, at...)) for v in wₜ] .- collect(weak_wave_modes(ℂ, ρ, n))
```

[`squared_moveout_velocities`](@ref) gives the squares of the short-spread
moveout velocities of a horizontal reflector in a medium of vertical axis,
exact for any strength of the anisotropy. The squares are returned because the
qSV one is negative in many rocks, where no moveout velocity exists:

```@example waves
squared_moveout_velocities(ℂ, ρ)
```

In Thomsen's parameters:

```@example waves
map(simplify, squared_moveout_velocities(ℂᵥ, ϱ))
```

## Any symmetry

The stiffness of a rock without any symmetry, the Lavoux limestone
[mensch1997](@cite), is read from its matrix in the two-index notation
([`stiffness_from_voigt`](@ref)), in GPa, and its density in kg/m³:

```@example waves
using JSON

data = JSON.parsefile(joinpath(pkgdir(TensND), "data", "literature", "menschGJI1997.json"))
rows = data["tables"]["lavoux_stiffness"]
ℂₗ = stiffness_from_voigt(Float64[row[c] for row in rows, c in ("1", "2", "3", "4", "5", "6")])
ρₗ = data["quantities"]["lavoux_density"]["value"] / 1000       # g/cm³
nothing # hide
```

```@example waves
Main.show_literature("menschGJI1997") # hide
```

[`weak_plane_waves`](@ref) gives the three waves to first order in the
deviation of the stiffness from isotropy, each as a pair
`(velocity, polarization)`, without any reference medium
([Weak anisotropy of any symmetry](@ref th-weak-anisotropy)):

```@example waves
nₗ = [0.3, -0.5, 0.8]
map(first, weak_plane_waves(ℂₗ, ρₗ, nₗ)), phase_velocities(ℂₗ, ρₗ, nₗ)
```

[`generalized_thomsen_parameters`](@ref) returns the seventeen dimensionless
parameters that set these velocities about the vertical
``\underline{e}_3`` [mensch1997](@cite):

```@example waves
map(x -> round(x; digits = 4), generalized_thomsen_parameters(ℂₗ))
```

For a transversely isotropic medium of axis ``\underline{e}_3`` they reduce to
Thomsen's: ``\varepsilon_x = \varepsilon_y = \varepsilon`` and
``\gamma_y = \gamma`` exactly, ``\delta_x = \delta_y`` equal to ``\delta`` to
first order.

```@example waves
g = generalized_thomsen_parameters(ℂ)
(g.εx, g.εy, p.ε), (g.γy, p.γ), (g.δx, g.δy, p.δ)
```

For an orthotropic stiffness in its axes, given by its nine components
``C_{IJ}``, the parameters are ratios of these components, and those that
couple the axes vanish:

```@example waves
c₁₁, c₂₂, c₃₃, c₁₂, c₁₃, c₂₃, c₄₄, c₅₅, c₆₆ = symbols("C_11 C_22 C_33 C_12 C_13 C_23 C_44 C_55 C_66", positive = true)
ℂₒ = TensOrtho(c₁₁, c₂₂, c₃₃, c₁₂, c₁₃, c₂₃, c₄₄, c₅₅, c₆₆, CanonicalBasis{3, Sym}())
map(simplify, generalized_thomsen_parameters(ℂₒ))
```

In the plane ``(\underline{e}_1, \underline{e}_3)``, a plane of symmetry, the
first-order modulus of the qP wave is
``\rho V_P^2 \simeq \underline{n}\cdot\boldsymbol{\Gamma}(\underline{n})\cdot\underline{n} = C_{11}\sin^4\vartheta + C_{33}\cos^4\vartheta + (2C_{13} + 4C_{55})\sin^2\vartheta\cos^2\vartheta``.
Both differences reduce to zero:

```@example waves
nₛ = [sin(ϑ), 0, cos(ϑ)]
ρV² = ϱ * weak_plane_waves(ℂₒ, ϱ, nₛ).qP[1]^2
nΓn = transpose(nₛ) * Array(canonical_components(acoustic_tensor(ℂₒ, nₛ))) * nₛ
simplify(ρV² - nΓn), simplify(ρV² - (c₁₁ * sin(ϑ)^4 + c₃₃ * cos(ϑ)^4 + (2c₁₃ + 4c₅₅) * sin(ϑ)^2 * cos(ϑ)^2))
```

[`shear_splitting`](@ref) measures how far apart the two shear waves of a
direction travel: the relative difference of their velocities, the
polarization of the fast one, and the delay of the slow one per unit length of
path (s/km here) ([Shear-wave splitting](@ref th-splitting)):

```@example waves
shear_splitting(ℂₗ, ρₗ, [0.0, 0.0, 1.0])
```

## Numbers, dual numbers, symbols

Every function above takes numbers and the dual numbers of ForwardDiff: the
derivative of a velocity with respect to the density, or to a component of the
stiffness, is one call.

```@example waves
using ForwardDiff

ForwardDiff.derivative(r -> phase_velocities(ℂ, r, n)[1], ρ), -w.velocity[1] / (2ρ)
```

```@example waves
ForwardDiff.derivative(c -> phase_velocities(ti_from_constants(c, 15.7, 21.8, 47.8, 16.4, [0.0, 0.0, 1.0]), ρ, n)[1], 52.9)
```

`wave_modes`, `christoffel_invariants`, `thomsen_parameters`, `anellipticity`,
`weak_wave_modes`, `squared_moveout_velocities`, `weak_plane_waves` and
`generalized_thomsen_parameters` also take symbols, as the examples above show. `plane_waves` and
`shear_splitting` sort the waves by their velocities, a comparison that symbols
do not allow: they throw an `ArgumentError` on a symbolic stiffness.

A direction in which an eigenvalue of ``\boldsymbol{\Gamma}/\rho`` is not
positive carries no wave: the stiffness is not positive definite there, or the
density is not positive, and `plane_waves` throws a `DomainError`. With
``C_{1212} < 0``, the transverse wave along ``\underline{e}_1`` polarized along
``\underline{e}_2`` does not exist:

```@example waves
ℂ₋ = from_mandel(Matrix(Diagonal([1.0, 1, 1, 1, 1, -1])))
try
    plane_waves(ℂ₋, 1.0, [1.0, 0.0, 0.0])
catch e
    e
end
```
