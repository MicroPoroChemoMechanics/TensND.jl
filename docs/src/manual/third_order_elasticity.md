# [Third-order elasticity](@id man-third-order)

!!! tip "See also"
    **Theory**: [Third-order elasticity](@ref th-third-order) · [Plane waves](@ref th-plane-waves)

    **Manual**: [Plane waves](@ref man-plane-waves) · [Elastic moduli](@ref man-elastic-moduli)

    **Tutorials**: [Stress-induced anisotropy](@ref tut-stress-anisotropy)

A solid isotropic at orders two and three stores, per unit reference volume,
the strain energy

```math
W = \tfrac12\lambda\,(\operatorname{tr}\boldsymbol{E})^2 + \mu\operatorname{tr}\boldsymbol{E}^2
  + \tfrac13(l - m)\,(\operatorname{tr}\boldsymbol{E})^3 + m\operatorname{tr}\boldsymbol{E}\operatorname{tr}\boldsymbol{E}^2 + n\det\boldsymbol{E}
```

in the Green–Lagrange strain
``\boldsymbol{E} = \tfrac12(\boldsymbol{F}^{\mathsf T}\boldsymbol{F} - \boldsymbol{1})``,
with Lamé's constants ``(\lambda, \mu)`` and Murnaghan's ``(l, m, n)``. Under a
static deformation of gradient ``\boldsymbol{F}`` its small waves depend on the
stress: the solid becomes anisotropic. This page gives the syntax; the
derivations are on [Third-order elasticity](@ref th-third-order).

| to … | write |
| :--- | :--- |
| define the solid | `s = ThirdOrderIsotropic(λ, μ, l, m, n)` |
| read its constants in other notations | `third_order_constants(s)` |
| the energy, the stress, the tangent | `strain_energy(s, E)`, `cauchy_stress(s, F)`, `tangent_stiffness(s, E)` |
| the tensor of the small motions under a deformation | `acoustoelastic_tensor(s, F)` |
| the waves under a deformation | `prestressed_waves(s, ρ₀, F, N)` |
| the first-order effect of a stress on the velocities | `acoustoelastic_coefficients(s)` |

## The solid and its constants

[`ThirdOrderIsotropic`](@ref) takes the five constants, here in GPa, of the
order of those of a sandstone:

```@example toe
using TensND, LinearAlgebra

s = ThirdOrderIsotropic(8.3, 11.9, -3.4e3, -6.8e3, -6.7e3)
```

[`third_order_constants`](@ref) returns the same third-order constants in the
two other usual notations: the components ``C_{111}, C_{112}, C_{123}`` of
``\partial^3W/\partial\boldsymbol{E}^3`` in the two-index notation, and Landau's
``A, B, C``:

```@example toe
third_order_constants(s)
```

## Energy, stress and tangent

The strain is a symmetric ``3\times3`` matrix of canonical components; the
deformation gradient a ``3\times3`` matrix. [`strain_energy`](@ref) returns
``W``, [`cauchy_stress`](@ref) the Cauchy stress
``\boldsymbol{\sigma} = \boldsymbol{F}\cdot\boldsymbol{S}\cdot\boldsymbol{F}^{\mathsf T}/\det\boldsymbol{F}``,
with ``\boldsymbol{S} = \partial W/\partial\boldsymbol{E}`` the second
Piola–Kirchhoff stress:

```@example toe
F = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 0.999]           # a shortening of 10⁻³ along e₃
E = (transpose(F) * F - I) / 2
strain_energy(s, E), cauchy_stress(s, F)
```

[`tangent_stiffness`](@ref) returns
``\mathbb{C}(\boldsymbol{E}) = \partial^2W/\partial\boldsymbol{E}^2``, a tensor with
the minor and major symmetries; at ``\boldsymbol{E} = \boldsymbol{0}`` it is the
isotropic stiffness of ``(\lambda, \mu)``. Its matrix in the two-index notation
([`voigt_stiffness`](@ref)):

```@example toe
round.(voigt_stiffness(tangent_stiffness(s, E)); digits = 3)
```

## The tensor of the small motions

A small motion superposed on the static deformation obeys, in the reference
coordinates, ``\rho_0\ddot u_i = \partial_{X_j}(A_{ijkl}\,\partial_{X_l}u_k)``,
with

```math
A_{ijkl} = \frac{\partial^2W}{\partial F_{ij}\,\partial F_{kl}} = \delta_{ik}\,S_{jl} + F_{ip}\,F_{kq}\,\mathbb{C}_{pjql}(\boldsymbol{E}) .
```

[`acoustoelastic_tensor`](@ref) returns it as a tensor of order four
(`AbstractTens{4,3}`), whose components in the canonical basis
[`canonical_components`](@ref) reads. It has the major symmetry, and under
stress not the minor ones:

```@example toe
A = acoustoelastic_tensor(s, F)
Aᵢⱼₖₗ = Array(canonical_components(A))
Aᵢⱼₖₗ ≈ permutedims(Aᵢⱼₖₗ, (3, 4, 1, 2)), Aᵢⱼₖₗ ≈ permutedims(Aᵢⱼₖₗ, (2, 1, 3, 4))
```

## Waves under a deformation

[`prestressed_waves`](@ref) solves the Christoffel problem of ``A`` along a
direction ``\underline{N}`` of the reference configuration, for the reference
density ``\rho_0`` (g/cm³ here, velocities in km/s), and returns a
[`PlaneWaves`](@ref) sorted by decreasing velocity:

```@example toe
ρ₀ = 2.3
N = [0.0, 0.0, 1.0]
w = prestressed_waves(s, ρ₀, F, N)
```

Its velocities are natural velocities: a length of the reference configuration
over a travel time, what a measurement between two faces of a sample gives
when the reference length is kept. The true velocity, a length of the deformed
solid over the same time, is the natural one times
``|\boldsymbol{F}\cdot\underline{N}|``. Against the velocities of the unstressed
solid:

```@example toe
λ, μ = s.λ, s.μ
w.velocity .* norm(F * N), phase_velocities(iso_stiffness(λ + 2μ / 3, μ), ρ₀, N)
```

A deformation that makes an eigenvalue of the acoustic tensor nonpositive
throws a `DomainError`, as `plane_waves` does.

## The first-order effect of a stress

[`acoustoelastic_coefficients`](@ref) returns, in closed form
([johnson1996](@citet), Eqs. (18)–(20), after [hughes1953](@citet)), the derivatives of ``\rho_0v^2`` (``v`` the true velocity)
with respect to a compressive stress ``\sigma`` at ``\sigma = 0``: uniaxial
along the propagation (`P_parallel`, `S_parallel`), uniaxial normal to it
(`P_perpendicular`, and `S_perpendicular_parallel` or
`S_perpendicular_perpendicular` for a shear wave polarized along the stress or
normal to it), hydrostatic (`P_hydrostatic`, `S_hydrostatic`); and the two
coefficients ``a_P`` and ``b_S`` of the anisotropy that a uniaxial stress
induces [johnson1996](@cite):

```@example toe
c = acoustoelastic_coefficients(s)
```

The same derivative by ForwardDiff on `prestressed_waves`, for a P wave along
a compression ``\sigma`` applied along ``\underline{e}_3``, the deformation being
taken to first order, ``\boldsymbol{F} = \boldsymbol{1} - \mathbb{S}:(\sigma\,\underline{e}_3\otimes\underline{e}_3)``:

```@example toe
using ForwardDiff

Eₘ, νₘ = young_poisson(iso_stiffness(λ + 2μ / 3, μ))
Fσ(σ) = [1 + νₘ * σ / Eₘ 0 0; 0 1 + νₘ * σ / Eₘ 0; 0 0 1 - σ / Eₘ]
true_P(σ) = norm(Fσ(σ) * N) * prestressed_waves(s, ρ₀, Fσ(σ), N).velocity[1]
ForwardDiff.derivative(σ -> ρ₀ * true_P(σ)^2, 0.0), c.P_parallel
```

## A background anisotropic at rest

[`strain_stiffening`](@ref) returns the change of stiffness that a strain brings
about through the third-order constants, ``\partial^3W/\partial\boldsymbol{E}^3 : \boldsymbol{E}``.
Added to the stiffness at rest of an anisotropic rock, it is the model of
[prioul2004](@citet) ([Third-order elasticity](@ref th-third-order), section 5).
[`third_order_from_constants`](@ref) builds the solid from the constants of the
two-index notation, as they are usually published; ``\lambda`` and ``\mu`` play
no part in the stiffening:

```@example toe
sₚ = third_order_from_constants(0.0, 0.0, -11300.0, -4800.0, 5800.0)       # GPa
ℂ⁰ = ti_from_constants(36.5, 14.9, 15.7, 24.6, 5.9, [0.0, 0.0, 1.0]; form = :components, role = :stiffness)
Eₕ = Array(canonical_components(inv(ℂ⁰) ⊡ Tens(-0.01 * Matrix(1.0I, 3, 3))))   # 10 MPa more of pressure
round.(voigt_stiffness(ℂ⁰ + strain_stiffening(sₚ, Eₕ)); digits = 2)
```

The tangent of the isotropic solid is the isotropic stiffness of ``\lambda``,
``\mu`` plus this stiffening:

```@example toe
tangent_stiffness(s, Matrix(E)) ≈ TensISO{3}(3λ + 2μ, 2μ) + strain_stiffening(s, Matrix(E))
```

## With symbols

The constants in the other notations, the energy, the stress, the tangent and
the coefficients are polynomial or rational in the five constants, so they
take symbols. The other notations first:

```@example toe
using SymPyPythonCall

λₛ, μₛ, lₛ, mₛ, nₛ = symbols("lambda mu l m n", real = true)
sₛ = ThirdOrderIsotropic(λₛ, μₛ, lₛ, mₛ, nₛ)
third_order_constants(sₛ)
```

Under a hydrostatic strain ``\boldsymbol{E} = e\,\boldsymbol{1}`` the tangent
stays isotropic, of Lamé constants ``\lambda + (6l - 2m + n)\,e`` and
``\mu + (3m - n/2)\,e``. The entries ``C_{12}`` and ``C_{44}`` of its two-index
matrix, and their differences with these values:

```@example toe
e = symbols("e", real = true)
Cₑ = simplify.(voigt_stiffness(tangent_stiffness(sₛ, [i == j ? e : zero(e) for i in 1:3, j in 1:3])))
Cₑ[1, 2], Cₑ[4, 4]
```

```@example toe
simplify(Cₑ[1, 2] - (λₛ + (6lₛ - 2mₛ + nₛ) * e)), simplify(Cₑ[4, 4] - (μₛ + (3mₛ - nₛ / 2) * e))
```

The coefficients of the first-order effect of a stress, among them ``a_P`` and
``b_S``:

```@example toe
cₛ = acoustoelastic_coefficients(sₛ)
simplify(cₛ.P_hydrostatic), cₛ.a_P, cₛ.b_S
```

The constants of `s` substituted in ``a_P`` and ``b_S`` give the numeric values
above:

```@example toe
at = (λₛ => s.λ, μₛ => s.μ, lₛ => s.l, mₛ => s.m, nₛ => s.n)
SymPyPythonCall.N(subs(cₛ.a_P, at...)) - c.a_P, SymPyPythonCall.N(subs(cₛ.b_S, at...)) - c.b_S                # N names the direction here
```
