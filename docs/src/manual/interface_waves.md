# [Waves at an interface](@id man-interface-waves)

!!! tip "See also"
    **Theory**: [Waves at an interface](@ref th-interface-waves) · [Characteristic roots](@ref th-characteristic-roots) · [Plane waves](@ref th-plane-waves)

    **Manual**: [Plane waves](@ref man-plane-waves) · [Elastic moduli](@ref man-elastic-moduli)

    **Tutorials**: [Reflection and transmission at an interface](@ref tut-reflection)

A wave ``\underline{u} = a\,\underline{U}\exp[\mathrm{i}\omega(\underline{s}\cdot\underline{x} - t)]``
meets a plane interface of unit normal ``\underline{\nu}``. All the waves it
excites share the tangential part ``\underline{s}_t`` of its slowness
(Snell–Descartes), and in each medium the normal slowness ``\xi`` of
``\underline{s} = \underline{s}_t + \xi\,\underline{\nu}`` solves a Christoffel
problem with six roots, real or in conjugate pairs
([Waves at an interface](@ref th-interface-waves)). This page gives the syntax.

| to … | write |
| :--- | :--- |
| the six waves of a medium for a tangential slowness | `forward, backward = interface_waves(ℂ, ρ, ν, sₜ)` |
| the reflected and transmitted waves of an incident one | `reflection_transmission(ℂ₁, ρ₁, ℂ₂, ρ₂, ν, incident)` |
| the same across a linear slip interface | `reflection_transmission(…; slip = ω * Z)` |
| the characteristic roots of a stiffness of the plane | `characteristic_pair(ℂ)`, `characteristic_roots(ℂ)` |
| its compliances in the engineering form | `engineering_compliances(ℂ)` |

## The six waves of a medium

The two media of the examples are isotropic, built from their P and S
velocities (km/s) and their densities (g/cm³), so that slownesses are in s/km:

```@example interface
using TensND, LinearAlgebra

iso(α, β, ρ) = iso_stiffness(ρ * (α^2 - 4β^2 / 3), ρ * β^2)       # bulk and shear moduli, GPa
ℂ₁, ρ₁ = iso(2.9, 1.4, 2.3), 2.3
ℂ₂, ρ₂ = iso(3.5, 2.1, 2.4), 2.4
ν = [0.0, 0.0, 1.0]
p = sind(30) / 2.9                     # tangential slowness of a P wave at 30^\circ in the medium 1
forward, backward = interface_waves(ℂ₁, ρ₁, ν, [p, 0.0, 0.0])
nothing # hide
```

[`interface_waves`](@ref) returns two tuples of three [`InterfaceWave`](@ref):
the `forward` waves carry energy, or decay, towards ``+\underline{\nu}``, the
`backward` ones towards ``-\underline{\nu}``; each tuple is sorted by increasing
``|\xi|``, the quasi-longitudinal wave first. The part of `sₜ` along
``\underline{\nu}`` is removed. Each wave has the fields `normal_slowness`
(``\xi``, complex), `slowness`, `polarization` (unit, complex), `traction`
(``\underline{\tau} = (\underline{\nu}\cdot\mathbb{C}\cdot\underline{s})\cdot\underline{U}``),
`flux` (the energy flux through the interface per unit
``\tfrac12\omega^2|a|^2``) and `homogeneous`:

```@example interface
[(w.normal_slowness, w.flux, w.homogeneous) for w in forward]
```

The two shear waves of an isotropic medium share their normal slowness; their
polarizations are then any orthonormal pair of their plane. A real ``\xi`` is a
homogeneous wave; past its critical angle a wave turns evanescent, with a
complex ``\xi`` of positive imaginary part for a wave that decays towards
``+\underline{\nu}``. In the medium 2 the P wave is evanescent for a tangential
slowness above ``1/3.5`` s/km:

```@example interface
[(w.normal_slowness, w.homogeneous) for w in first(interface_waves(ℂ₂, ρ₂, ν, [0.32, 0.0, 0.0]))]
```

## Reflection and transmission

The incident wave is a forward wave of the medium 1, which lies on the side
``\underline{\nu}\cdot\underline{x} < 0``. [`reflection_transmission`](@ref)
returns the three `reflected` waves of the medium 1, the three `transmitted`
waves of the medium 2, their complex amplitudes `R` and `T` per unit incident
amplitude, and `energy`, the fractions of the incident energy flux carried by
each wave, as a pair (reflected, transmitted):

```@example interface
incident = forward[1]
r = reflection_transmission(ℂ₁, ρ₁, ℂ₂, ρ₂, ν, incident)
abs.(r.R), abs.(r.T)
```

```@example interface
r.energy, sum(r.energy[1]) + sum(r.energy[2])
```

The amplitudes are those of the displacement on unit polarizations, whose sign
or phase is arbitrary: their moduli, the displacement ``R\,\underline{U}`` and
the energy fractions are meaningful. Two waves of the same normal slowness (the
shear waves of an isotropic medium) are read together, by the sum of their
displacements or of their energies. An evanescent wave carries no energy
across the interface.

The interface above is welded: displacement and traction are continuous. With
`slip`, a ``3\times3`` matrix ``\omega\boldsymbol{Z}`` (angular frequency times
the compliance of the interface, in the unit of a length per stress), the
displacement jumps by ``\boldsymbol{Z}`` times the traction
[schoenberg1980](@cite):

```@example interface
Z = Diagonal([2.0e-4, 2.0e-4, 1.0e-4])                  # tangential and normal compliances, km/GPa
rₛ = reflection_transmission(ℂ₁, ρ₁, ℂ₂, ρ₂, ν, incident; slip = 2π * 50.0 * Z)
abs.(rₛ.R), sum(rₛ.energy[1]) + sum(rₛ.energy[2])
```

Both functions take dual numbers, through the derivative of each root: the
derivative of a coefficient with respect to the density of the medium 2,

```@example interface
using ForwardDiff

ForwardDiff.derivative(ϱ -> abs(reflection_transmission(ℂ₁, ρ₁, ℂ₂, ϱ, ν, incident).R[1]), ρ₂)
```

The roots are found numerically and sorted by comparisons: symbols are refused
with an `ArgumentError`.

## Characteristic roots of a stiffness of the plane

Without inertia, the same problem in the plane gives the displacements
``\underline{U}\,f(x_1 + \mu x_2)`` in equilibrium in a medium of stiffness
``\mathbb{C}``. With the compliances of the plane in the engineering form,
``\varepsilon_{11} = a_{11}\sigma_{11} + a_{12}\sigma_{22} + a_{16}\sigma_{12}``,
``\varepsilon_{22} = a_{12}\sigma_{11} + a_{22}\sigma_{22} + a_{26}\sigma_{12}``,
``2\varepsilon_{12} = a_{16}\sigma_{11} + a_{26}\sigma_{22} + a_{66}\sigma_{12}``,
the characteristic roots ``\mu`` solve [lekhnitskii1981](@cite)

```math
a_{11}\mu^4 - 2a_{16}\mu^3 + (2a_{12} + a_{66})\mu^2 - 2a_{26}\mu + a_{22} = 0 ,
```

which has no real root for a positive definite stiffness: two roots
``\mu_1, \mu_2`` of positive imaginary part and their conjugates
([Characteristic roots](@ref th-characteristic-roots)). The example is an
orthotropic stiffness of the plane built from its engineering constants
``(E_1, E_2, \nu_{12}, G_{12})`` in its material frame, here the canonical one
([`ortho_from_constants`](@ref)):

```@example interface
ℂₚ = ortho_from_constants(40.0, 10.0, 0.25, 5.0, Basis(0.0); role = :stiffness)
a = engineering_compliances(ℂₚ)                 # (a₁₁, a₂₂, a₁₂, a₆₆, a₁₆, a₂₆)
```

[`characteristic_pair`](@ref) returns the sum ``s = \mu_1 + \mu_2`` and the
product ``p = \mu_1\mu_2``, [`characteristic_roots`](@ref) the two roots:

```@example interface
characteristic_pair(ℂₚ), characteristic_roots(ℂₚ)
```

```@example interface
a₁₁, a₂₂, a₁₂, a₆₆, a₁₆, a₂₆ = a
[a₁₁ * μ^4 - 2a₁₆ * μ^3 + (2a₁₂ + a₆₆) * μ^2 - 2a₂₆ * μ + a₂₂ for μ in characteristic_roots(ℂₚ)]
```

The roots of an isotropic tensor merge at ``\mathrm{i}``, where they stop being
differentiable functions of the stiffness while their sum and product stay
smooth. An isotropic, a square or an orthotropic tensor has its pair in closed
form, in its own axes turned to the canonical ones; any other tensor has its
quartic solved numerically. The same material turned by ``\pi/6``, as an
orthotropic tensor and as a general one:

```@example interface
ℂᵣ = ortho_from_constants(40.0, 10.0, 0.25, 5.0, Basis(π / 6); role = :stiffness)
characteristic_pair(ℂᵣ), characteristic_pair(Tens(Array(canonical_components(ℂᵣ))))
```

```@example interface
characteristic_roots(iso_stiffness(Val(2), 5.0, 2.0))
```
