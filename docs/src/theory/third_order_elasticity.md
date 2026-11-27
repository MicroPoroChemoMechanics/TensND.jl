# [Third-order elasticity](@id th-third-order)

!!! tip "See also"
    **Manual**: [Third-order elasticity](@ref man-third-order)

    **Tutorials**: [Stress-induced anisotropy](@ref tut-stress-anisotropy)

A solid under stress transmits small waves at velocities that depend on the
stress: a rock compressed along one direction becomes anisotropic. Linear
elasticity cannot describe it, since its stiffness does not depend on the
state. This page writes the energy of a solid isotropic at orders two and
three, derives the tensor that governs a small wave superposed on a static
deformation, and obtains the first-order effect of a stress on the
velocities.

## 1. The energy and its derivatives

Let ``\boldsymbol{F}`` be the gradient of a deformation from a reference
configuration of density ``\rho_0``, and
``\boldsymbol{E} = \tfrac12(\boldsymbol{F}^{\mathsf T}\cdot\boldsymbol{F} - \boldsymbol{1})``
the Green–Lagrange strain. An isotropic energy per unit reference volume is a
function of the invariants of ``\boldsymbol{E}``; cut at the third order, it
has five constants, Lamé's ``\lambda``, ``\mu`` and Murnaghan's ``l``, ``m``,
``n``:

```math
W = \tfrac12\lambda\,(\operatorname{tr}\boldsymbol{E})^2 + \mu\operatorname{tr}\boldsymbol{E}^2
  + \tfrac13(l - m)\,(\operatorname{tr}\boldsymbol{E})^3 + m\operatorname{tr}\boldsymbol{E}\operatorname{tr}\boldsymbol{E}^2 + n\det\boldsymbol{E} .
```

The derivative of ``\det\boldsymbol{E}`` is its cofactor,
``\boldsymbol{E}^2 - \operatorname{tr}\boldsymbol{E}\,\boldsymbol{E} + \tfrac12[(\operatorname{tr}\boldsymbol{E})^2 - \operatorname{tr}\boldsymbol{E}^2]\,\boldsymbol{1}``
(Cayley–Hamilton), so that the second Piola–Kirchhoff stress is

```math
\boldsymbol{S} = \frac{\partial W}{\partial\boldsymbol{E}}
= \bigl[\lambda\operatorname{tr}\boldsymbol{E} + (l - m)(\operatorname{tr}\boldsymbol{E})^2 + m\operatorname{tr}\boldsymbol{E}^2\bigr]\boldsymbol{1}
+ 2(\mu + m\operatorname{tr}\boldsymbol{E})\,\boldsymbol{E} + n\operatorname{cof}\boldsymbol{E} ,
```

and its tangent ``\mathbb{C}(\boldsymbol{E}) = \partial\boldsymbol{S}/\partial\boldsymbol{E}``
is the isotropic stiffness of ``\lambda``, ``\mu`` plus a term linear in
``\boldsymbol{E}``. The third derivative of ``W`` is a constant isotropic
tensor of order six; its two-index components are
``C_{111} = 2l + 4m``, ``C_{112} = 2l``, ``C_{123} = 2l - 2m + n``, and Landau's
constants ``A = n``, ``B = m - n/2``, ``C = l - m + n/2`` are a third
notation of the same energy ([`third_order_constants`](@ref)). The Cauchy
stress is ``\boldsymbol{\sigma} = \boldsymbol{F}\cdot\boldsymbol{S}\cdot\boldsymbol{F}^{\mathsf T}/\det\boldsymbol{F}``.

## 2. A small wave on a deformed solid

Let the solid be held in the homogeneous deformation ``\boldsymbol{F}``, and
add a small displacement ``\underline{u}`` written as a function of the
reference position ``\underline{X}``. The first Piola–Kirchhoff stress
``\boldsymbol{P} = \partial W/\partial\boldsymbol{F}`` varies by
``\delta P_{ij} = A_{ijkl}\,\partial u_k/\partial X_l`` with

```math
A_{ijkl} = \frac{\partial^2 W}{\partial F_{ij}\,\partial F_{kl}}
= \delta_{ik}\,S_{jl} + F_{ip}\,F_{kq}\,\mathbb{C}_{pjql}(\boldsymbol{E}) ,
```

obtained by differentiating ``\boldsymbol{P} = \boldsymbol{F}\cdot\boldsymbol{S}``
with ``\partial E_{pq}/\partial F_{kl} = \tfrac12(\delta_{ql}F_{kp} + \delta_{pl}F_{kq})``.
The equation of motion in reference coordinates,
``\rho_0\,\ddot u_i = \partial\,\delta P_{ij}/\partial X_j``, has the form of
the [Christoffel problem](@ref th-plane-waves) with ``\mathbb{C}`` replaced by
``\mathbb{A}``: a plane wave of *reference* direction ``\underline{N}`` has

```math
\bigl(\boldsymbol{\Gamma} - \rho_0 W^2\,\boldsymbol{1}\bigr)\cdot\underline{U} = \underline{0},
\qquad
\Gamma_{ik} = N_j\,A_{ijkl}\,N_l .
```

``\mathbb{A}`` has the major symmetry, so ``\boldsymbol{\Gamma}`` is symmetric,
but not the minor ones: the stress enters through
``\delta_{ik}S_{jl}``. ``W`` is the *natural velocity*, a reference length
over a travel time: what is measured between two faces of a sample when the
distance between them is taken in the unstressed state. The *true velocity*
``v``, a length of the deformed solid over the same time, is
``v = |\boldsymbol{F}\cdot\underline{N}|\,W``. The group velocity follows as in
[Plane waves](@ref th-plane-waves), with ``\mathbb{A}`` and ``\rho_0``
([`acoustoelastic_tensor`](@ref), [`prestressed_waves`](@ref)).

## 3. Hydrostatic pressure

Take ``\boldsymbol{F} = (1 + e)\,\boldsymbol{1}``, to first order in ``e``. Then
``\boldsymbol{E} \simeq e\,\boldsymbol{1}``,
``\boldsymbol{S} \simeq (3\lambda + 2\mu)\,e\,\boldsymbol{1} = 3K e\,\boldsymbol{1}``
with ``K = \lambda + 2\mu/3``, and the pressure is ``p = -3Ke``. The
cofactor of ``e\,\boldsymbol{1}`` being ``e^2\,\boldsymbol{1}``, its derivative
at that point is
``e\,(\boldsymbol{1}\otimes\boldsymbol{1} - \mathbb{I})`` to first order, with
``\mathbb{I}`` the identity of the symmetric tensors, and the tangent is
isotropic,

```math
\mathbb{C} \simeq \lambda'\,\boldsymbol{1}\otimes\boldsymbol{1} + 2\mu'\,\mathbb{I},
\qquad
\lambda' = \lambda + (6l - 2m + n)\,e,\qquad
\mu' = \mu + (3m - n/2)\,e .
```

Along ``\underline{N} = \underline{e}_1``, the P wave has
``\rho_0 W_P^2 = A_{1111} = 3Ke + (1 + e)^2(\lambda' + 2\mu')``, and an S wave
polarized along ``\underline{e}_2`` has
``\rho_0 W_S^2 = A_{2121} = 3Ke + (1 + e)^2\mu'``. With ``v = (1 + e)W`` and
``e = -p/(3K)``,

```math
\begin{gathered}
\rho_0 v_P^2 = \lambda + 2\mu - \frac{p}{3K}\,(6l + 4m + 7\lambda + 10\mu),
 \\
\rho_0 v_S^2 = \mu - \frac{p}{3K}\,\Bigl(3m - \frac{n}{2} + 3\lambda + 6\mu\Bigr),
\end{gathered}
```

the expressions (20a) and (20b) of [johnson1996](@citet), who refer them to
[hughes1953](@citet).

## 4. Uniaxial stress and the anisotropy it induces

A compression ``\sigma`` along the unit vector ``\underline{m}`` produces, to
first order, the strain
``\boldsymbol{E} = -\sigma\,[(1 + \nu)\,\underline{m}\otimes\underline{m} - \nu\,\boldsymbol{1}]/E_Y``
of linear elasticity (``E_Y`` and ``\nu`` Young's modulus and Poisson's
ratio). The same computation as above, with this strain, gives the slopes
``\mathrm{d}(\rho_0 v^2)/\mathrm{d}\sigma`` of the P and S waves that travel
along ``\underline{m}`` or across it, the expressions (18)–(19) of
[johnson1996](@citet); [`acoustoelastic_coefficients`](@ref) returns them,
and the tests of the package compare each one with the derivative of the
exact tensor ``\mathbb{A}``.

The dependence on the angle ``\theta`` between the propagation and
``\underline{m}`` follows from the structure of ``\mathbb{A}``. To first order
in ``\sigma`` it is the isotropic stiffness plus three contributions, each
linear in ``\boldsymbol{E}``, hence in the tensors ``\boldsymbol{1}`` and
``\underline{m}\otimes\underline{m}``:

- the third-order term of ``\mathbb{C}(\boldsymbol{E})``, built from
  ``\boldsymbol{1}\otimes\boldsymbol{1}``, ``\mathbb{I}``,
  ``\boldsymbol{E}\otimes\boldsymbol{1} + \boldsymbol{1}\otimes\boldsymbol{E}``
  and the symmetrized product of ``\boldsymbol{E}`` and ``\boldsymbol{1}``;
- ``\delta_{ik}S_{jl}``, which adds ``\underline{N}\cdot\boldsymbol{S}\cdot\underline{N}\,\boldsymbol{1}``
  to ``\boldsymbol{\Gamma}``, the same for the three waves;
- ``F_{ip}F_{kq}\,\mathbb{C}_{pjql}``, which adds
  ``\boldsymbol{H}\cdot\boldsymbol{\Gamma}^0 + \boldsymbol{\Gamma}^0\cdot\boldsymbol{H}``
  with ``\boldsymbol{F} = \boldsymbol{1} + \boldsymbol{H}``.

None of them contains ``\underline{m}\otimes\underline{m}\otimes\underline{m}\otimes\underline{m}``,
which would need ``\boldsymbol{E}`` twice. The stressed solid is
transversely isotropic of axis ``\underline{m}``, and in
[Thomsen's parameters](@ref th-thomsen) each of these tensors has
``\varepsilon = \delta``: for
``\underline{m}\otimes\underline{m}\otimes\boldsymbol{1} + \boldsymbol{1}\otimes\underline{m}\otimes\underline{m}``,
``\delta C_{33} = 2``, ``\delta C_{13} = 1`` and
``\delta C_{11} = \delta C_{44} = 0``, so that ``\varepsilon = \delta = -1/C_{33}``;
for the symmetrized product of ``\underline{m}\otimes\underline{m}`` and
``\boldsymbol{1}``, ``\delta C_{33} = 1``, ``\delta C_{44} = 1/4`` and
``\delta C_{11} = \delta C_{13} = 0``, so that ``\varepsilon = \delta = -1/(2C_{33})``.
The first-order P velocity of an elliptical medium is
``\rho V_P^2 \simeq C_{33}(1 + 2\varepsilon\sin^2\theta)``, and the
difference of the two S moduli is ``(C_{66} - C_{44})\sin^2\theta``; the
last contribution adds to the SV modulus a term in
``(\underline{U}\cdot\underline{m})^2 = \sin^2\theta``. Both anisotropies go
as ``\sin^2\theta``:

```math
\frac{v_{P\parallel} - v_P(\theta)}{v_{P\parallel}} \simeq -\sigma\,a_P\sin^2\theta,
\qquad
\frac{v_{S_1}(\theta) - v_{S_2}(\theta)}{v_{S_1}(\theta)} \simeq -\sigma\,b_S\sin^2\theta ,
```

``S_1`` being polarized in the plane of ``\underline{m}`` and of the
propagation, ``S_2`` normal to it, and

```math
a_P = \frac{2\lambda + 5\mu + 2m}{2\mu\,(\lambda + 2\mu)},\qquad
b_S = \frac{4\mu + n}{8\mu^2} ,
```

read on the slopes at ``\theta = 90^\circ``. In rocks, where the third-order
constants are larger than the second-order ones by orders of magnitude
(cracks closing under stress), ``a_P`` and ``b_S`` are large and negative: the
velocities are highest along the compression.
