# [Elastic moduli](@id man-elastic-moduli)

<!-- crosslinks: written by .github/scripts/crosslinks.py from docs/crosslinks.toml; do not edit -->
!!! tip "See also"
    **Theory**: [Isotropic tensors](@ref th-isotropic)

    **Manual**: [Structured tensors](@ref man-structured) · [Symbolic and numeric](@ref man-symbolic-numeric) · [Parametrizations and conversions](@ref man-parametrizations)

    **Tutorials**: [Elastic Green's functions](@ref tut-green-function) · [The elastic sphere, solved symbolically](@ref tut-sphere-problems) · [Differentiation and interoperability](@ref tut-ad-interop) · [Performance of the structured types](@ref tut-performance) · [Differential operators, symbolically](@ref tut-symbolic-operators)
<!-- /crosslinks -->

An isotropic linear elastic material has two independent constants, and six of
them are in common use: the bulk modulus ``k``, Young's modulus ``E``, Lamé's
first parameter ``\lambda``, the shear modulus ``\mu``, Poisson's ratio ``\nu`` and
the P-wave (oedometric) modulus ``M``. Any two of them fix the others. TensND
makes the tensor the meeting point: the stiffness, or the compliance, is built
from any pair, and any modulus is read from either, so that a conversion is one
function inside another.

| to … | write |
| :--- | :--- |
| build a stiffness from two moduli | `iso_stiffness_<pair>(a, b)`, e.g. `iso_stiffness_young_poisson(E, ν)` or `iso_stiffness_E_ν(E, ν)` |
| build a compliance from two moduli | `iso_compliance_<pair>(a, b)`, e.g. `iso_compliance_young_poisson(E, ν)` |
| read one modulus | `young_modulus(ℂ)`, `poisson_ratio(ℂ)`, …, or `young_modulus(𝕊; role = :compliance)` |
| read two moduli | `young_poisson(ℂ)` (`E_ν(ℂ)`), `lame_shear(ℂ)` (`λ_μ(ℂ)`), …, with the same `role` |
| work in 2-D (plane strain) | the same builders with `Val(2)` first |
| switch between plane strain and plane stress | `plane_strain_stiffness`, `plane_stress_stiffness`, `plane_stress_to_plane_strain`, `plane_strain_to_plane_stress` |

## The stiffness and its moduli

In dimension 3 the stiffness is
``\mathbb{C} = 3k\,\mathbb{J} + 2\mu\,\mathbb{K} = (3\lambda + 2\mu)\,\mathbb{J} + 2\mu\,\mathbb{K}``,
a `TensISO{4,3}`. The other moduli follow from ``(\lambda, \mu)``:

```math
k = \lambda + \tfrac{2}{3}\mu, \qquad
E = \frac{\mu(3\lambda + 2\mu)}{\lambda + \mu}, \qquad
\nu = \frac{\lambda}{2(\lambda + \mu)}, \qquad
M = \lambda + 2\mu .
```

Each has its reader, [`bulk_modulus`](@ref), [`young_modulus`](@ref),
[`lame_lambda`](@ref), [`shear_modulus`](@ref), [`poisson_ratio`](@ref) and
[`pwave_modulus`](@ref), and each pair too, named after its two moduli in the
order ``k, E, \lambda, \mu, \nu, M``: [`bulk_shear`](@ref), [`young_poisson`](@ref),
[`lame_shear`](@ref), [`young_pwave`](@ref), …

```@example moduli
using TensND

ℂ = iso_stiffness_young_poisson(70.0, 0.33)   # an aluminum, GPa
(k = bulk_modulus(ℂ), E = young_modulus(ℂ), λ = lame_lambda(ℂ),
 μ = shear_modulus(ℂ), ν = poisson_ratio(ℂ), M = pwave_modulus(ℂ))
```

## Stiffness or compliance

The compliance is ``\mathbb{S} = \mathbb{C}^{-1} = a\,\mathbb{J} + b\,\mathbb{K}`` with
``a = 1/(3k)`` and ``b = 1/(2\mu)``. Every reader takes the keyword `role`, which
says which of the two it is given: `:stiffness`, the default, or `:compliance`.
Each modulus then has a formula of its own on each, written directly on the two
coefficients, never through an inverse:

| modulus | on ``\mathbb{C} = \alpha\,\mathbb{J} + \beta\,\mathbb{K}`` | on ``\mathbb{S} = a\,\mathbb{J} + b\,\mathbb{K}`` |
| :--- | :--- | :--- |
| ``k`` | ``\alpha/3`` | ``1/(3a)`` |
| ``\mu`` | ``\beta/2`` | ``1/(2b)`` |
| ``\lambda`` | ``(\alpha - \beta)/3`` | ``(b - a)/(3ab)`` |
| ``E`` | ``9k\mu/(3k + \mu)`` | ``3/(a + 2b)`` |
| ``\nu`` | ``(3k - 2\mu)/(2(3k + \mu))`` | ``(b - a)/(a + 2b)`` |
| ``M`` | ``\lambda + 2\mu`` | ``(b + 2a)/(3ab)`` |

```@example moduli
𝕊 = iso_compliance_young_poisson(70.0, 0.33)   # the same aluminum
young_poisson(𝕊; role = :compliance)
```

The compliance builders write ``\mathbb{S}`` directly as well: from ``(E, \nu)``,
``\mathbb{S} = (1 - 2\nu)/E\;\mathbb{J} + (1 + \nu)/E\;\mathbb{K}``. Nothing in a
`TensISO` tells a stiffness from a compliance: the role is the caller's to give,
and a compliance read as a stiffness gives the moduli of another material.

## Building from any pair

The stiffness of a pair is `iso_stiffness_<pair>`, and its compliance
`iso_compliance_<pair>` ([`iso_compliance`](@ref) for ``(k, \mu)``), with the same
names. The
fifteen pairs, and how each gives ``(\lambda, \mu)`` (the formulas are those of
the classical [conversion table](https://en.wikipedia.org/wiki/Elastic_modulus),
in this notation):

| pair | builder | ``\lambda`` | ``\mu`` |
| :--- | :--- | :--- | :--- |
| ``(k, E)`` | [`iso_stiffness_bulk_young`](@ref) | ``3k(3k - E)/(9k - E)`` | ``3kE/(9k - E)`` |
| ``(k, \lambda)`` | [`iso_stiffness_bulk_lame`](@ref) | ``\lambda`` | ``3(k - \lambda)/2`` |
| ``(k, \mu)`` | [`iso_stiffness`](@ref) | ``k - 2\mu/3`` | ``\mu`` |
| ``(k, \nu)`` | [`iso_stiffness_bulk_poisson`](@ref) | ``3k\nu/(1 + \nu)`` | ``3k(1 - 2\nu)/(2(1 + \nu))`` |
| ``(k, M)`` | [`iso_stiffness_bulk_pwave`](@ref) | ``(3k - M)/2`` | ``3(M - k)/4`` |
| ``(E, \lambda)`` | [`iso_stiffness_young_lame`](@ref) | ``\lambda`` | ``(E - 3\lambda + R)/4`` |
| ``(E, \mu)`` | [`iso_stiffness_young_shear`](@ref) | ``\mu(E - 2\mu)/(3\mu - E)`` | ``\mu`` |
| ``(E, \nu)`` | [`iso_stiffness_young_poisson`](@ref) | ``E\nu/((1 + \nu)(1 - 2\nu))`` | ``E/(2(1 + \nu))`` |
| ``(E, M)`` | [`iso_stiffness_young_pwave`](@ref) | ``(M - E + S)/4`` | ``(3M + E - S)/8`` |
| ``(\lambda, \mu)`` | [`iso_stiffness_lame_shear`](@ref) | ``\lambda`` | ``\mu`` |
| ``(\lambda, \nu)`` | [`iso_stiffness_lame_poisson`](@ref) | ``\lambda`` | ``\lambda(1 - 2\nu)/(2\nu)`` |
| ``(\lambda, M)`` | [`iso_stiffness_lame_pwave`](@ref) | ``\lambda`` | ``(M - \lambda)/2`` |
| ``(\mu, \nu)`` | [`iso_stiffness_shear_poisson`](@ref) | ``2\mu\nu/(1 - 2\nu)`` | ``\mu`` |
| ``(\mu, M)`` | [`iso_stiffness_shear_pwave`](@ref) | ``M - 2\mu`` | ``\mu`` |
| ``(\nu, M)`` | [`iso_stiffness_poisson_pwave`](@ref) | ``M\nu/(1 - \nu)`` | ``M(1 - 2\nu)/(2(1 - \nu))`` |

with ``R = \sqrt{E^2 + 9\lambda^2 + 2E\lambda}`` and
``S = \sqrt{E^2 + 9M^2 - 10EM}``. These two pairs solve a quadratic equation.
From ``(E, \lambda)`` the root kept is the only one with ``-1 < \nu < 1/2``. From
``(E, M)`` both roots are materials, one with ``\nu \geq 0`` and one with
``\nu \leq 0``, and ``E`` and ``M`` cannot tell them apart; the one kept is
``\nu \geq 0``.

The pair ``(\lambda, \nu)`` says nothing of ``\mu`` at ``\nu = 0``, where ``\lambda = 0``
as well. A conversion reads the result on the stiffness:

```@example moduli
E, ν = young_poisson(iso_stiffness_lame_shear(40.0, 26.0))
```

```@example moduli
young_shear(iso_stiffness_bulk_poisson(70.0, 0.25))
```

## Shorter names

The pairs without the P-wave modulus also have the names of their symbols, in
ASCII and with Greek letters, as exact aliases: `k_mu` and `k_μ` are
[`bulk_shear`](@ref), `E_nu` and `E_ν` are [`young_poisson`](@ref), `λ_μ` and
`lambda_mu` are [`lame_shear`](@ref), and the builders likewise,
`iso_stiffness_E_nu`, `iso_stiffness_E_ν`, `iso_stiffness_λ_μ`, … and
`iso_compliance_E_nu`, `iso_compliance_E_ν`, `iso_compliance_λ_μ`, …. The pairs
with ``M`` have none: `k_M` or `μ_M` name the moduli of a matrix in
MeanFieldHomogenization, and the two meanings would be read one for the other.

```@example moduli
k, μ = k_μ(iso_stiffness_E_ν(70.0, 0.33))
```

## Symbolic moduli

No formula depends on the sign of a value, so the moduli may be symbols. With
SymPy, the closed forms come out exactly:

```@example moduli
using SymPyPythonCall
λ, μ = symbols("λ μ", positive = true)
ℂ = iso_stiffness_lame_shear(λ, μ)
simplify(young_modulus(ℂ)), simplify(poisson_ratio(ℂ))
```

They go through ForwardDiff as well: a derivative of a modulus with respect to
another pair is that of the composition.

```@example moduli
using ForwardDiff
ForwardDiff.derivative(ν -> shear_modulus(iso_stiffness_young_poisson(70.0, ν)), 0.33)
```

## Two dimensions

A `TensISO{4,2}` is the in-plane stiffness of a 3-D material under plane
strain, ``\mathbb{C} = 2(\lambda + \mu)\,\mathbb{J} + 2\mu\,\mathbb{K}``. Its ``E``,
``\lambda``, ``\mu``, ``\nu`` and ``M`` are those of the 3-D material, while its bulk
modulus is the in-plane one, ``k = \lambda + \mu``, the factor of ``\mathbb{J}``
over the dimension, as for [`bulk_shear`](@ref). Every builder takes the
dimension as a first argument `Val(2)` (`Val(3)` is the default):

```@example moduli
ℂ₂ = iso_stiffness_young_poisson(Val(2), 70.0, 0.33)
bulk_shear(ℂ₂), young_poisson(ℂ₂)
```

On a 2-D compliance ``\mathbb{S} = a\,\mathbb{J} + b\,\mathbb{K}``, with
``a = 1/(2(\lambda + \mu))`` and ``b = 1/(2\mu)``, the readers use

```math
k = \frac{1}{2a}, \quad \mu = \frac{1}{2b}, \quad \lambda = \frac{b - a}{2ab}, \quad
E = \frac{3b - a}{2b^2}, \quad \nu = \frac{b - a}{2b}, \quad M = \frac{a + b}{2ab},
```

and the compliance builders take `Val(2)` as the stiffness builders do:

```@example moduli
𝕊₂ = iso_compliance_young_poisson(Val(2), 70.0, 0.33)
young_poisson(𝕊₂; role = :compliance)
```

The pairs with ``k`` then use the in-plane bulk modulus: ``\mu = k - \lambda``,
``\mu = k(1 - 2\nu)``, ``\mu = M - k``, and, from ``(k, E)``,
``\mu = (3k - \sqrt{9k^2 - 4Ek})/2``. That last quadratic has two roots of
physical materials, of Poisson's ratios ``\nu`` and ``-1/2 - \nu``; the one kept
is ``\nu \geq -1/4``.

## Plane strain and plane stress

A 2-D computation on a 3-D material assumes either plane strain,
``\varepsilon_{i3} = 0``, or plane stress, ``\sigma_{i3} = 0``.
[`plane_strain_stiffness`](@ref) and [`plane_stress_stiffness`](@ref) give the
in-plane stiffness of each. Plane strain keeps the in-plane block of
``\mathbb{C}``; plane stress eliminates ``\varepsilon_{33}`` from
``\sigma_{33} = 0``, which keeps ``\mu`` and replaces ``\lambda`` by

```math
\lambda^* = \frac{2\lambda\mu}{\lambda + 2\mu} = \frac{E\nu}{1 - \nu^2} .
```

```@example moduli
ℂ = iso_stiffness_young_poisson(70.0, 0.33)
lame_shear(plane_strain_stiffness(ℂ)), lame_shear(plane_stress_stiffness(ℂ))
```

Both are `TensISO{4,2}`, so a plane-stress state is computed with the
plane-strain formulas of another material. [`plane_stress_to_plane_strain`](@ref)
gives its moduli and [`plane_strain_to_plane_stress`](@ref) goes back:

```math
(E, \nu) \mapsto \left(\frac{E(1 + 2\nu)}{(1 + \nu)^2},\ \frac{\nu}{1 + \nu}\right),
\qquad
(E, \nu) \mapsto \left(\frac{E}{1 - \nu^2},\ \frac{\nu}{1 - \nu}\right).
```

```@example moduli
E′, ν′ = plane_stress_to_plane_strain(70.0, 0.33)
young_poisson(plane_stress_stiffness(ℂ)) .≈ (E′, ν′)
```
