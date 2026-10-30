# [Elastic moduli](@id man-elastic-moduli)

An isotropic linear elastic material has two independent constants, and six of
them are in common use: the bulk modulus ``k``, Young's modulus ``E``, Lamé's
first parameter ``\lambda``, the shear modulus ``\mu``, Poisson's ratio ``\nu`` and
the P-wave (oedometric) modulus ``M``. Any two of them fix the others. TensND
makes the stiffness tensor the meeting point: it is built from any pair and any
modulus is read from it, so that a conversion is one function inside another.

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

A compliance is the inverse of a stiffness, exactly in the class, so its moduli
are those of `inv(𝕊)`.

## Building from any pair

The stiffness of a pair is `iso_stiffness_<pair>`, with the same names. The
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
`iso_stiffness_E_nu`, `iso_stiffness_E_ν`, `iso_stiffness_λ_μ`, …. The pairs
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
