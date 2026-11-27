# [Waves at an interface](@id th-interface-waves)

!!! tip "See also"
    **Manual**: [Waves at an interface](@ref man-interface-waves)

    **Tutorials**: [Reflection and transmission at an interface](@ref tut-reflection)

A plane wave that meets the plane separating two media gives rise to reflected
and transmitted waves. This page derives which ones (the law of
Snell-Descartes), how their amplitudes follow from the conditions on the
plane, and why the energy is shared among them. The same problem without
inertia gives the characteristic roots of plane anisotropic elasticity.

## 1. The law of Snell-Descartes

Let the plane ``\underline{\nu}\cdot\underline{x} = 0``, of unit normal
``\underline{\nu}``, separate the medium 1 (``\underline{\nu}\cdot\underline{x} < 0``)
from the medium 2. A harmonic wave of slowness ``\underline{s}``
([Plane waves](@ref th-plane-waves)) is

```math
\underline{u} = a\,\underline{U}\,\exp\bigl[\mathrm{i}\omega(\underline{s}\cdot\underline{x} - t)\bigr],
\qquad
\underline{s} = \underline{s}_t + \xi\,\underline{\nu},\quad \underline{s}_t\cdot\underline{\nu} = 0 .
```

On the plane it varies as
``\exp[\mathrm{i}\omega(\underline{s}_t\cdot\underline{x} - t)]``. The conditions
that link the two media hold at every point of the plane and at every time:
all the waves involved have the frequency and the *tangential slowness*
``\underline{s}_t`` of the incident one. This is the law of Snell-Descartes.
For a wave of velocity ``V`` whose normal makes the angle ``\theta`` with
``\underline{\nu}``, ``|\underline{s}_t| = \sin\theta/V``: in isotropic media
it is the familiar ``\sin\theta_1/V_1 = \sin\theta_2/V_2``. In an anisotropic
medium ``V`` depends on the direction and the angle is not known in advance;
what is known is ``\underline{s}_t``, and the unknown is the *normal slowness*
``\xi``.

## 2. The six waves of a tangential slowness

The Christoffel problem
``(\boldsymbol{\Gamma}(\underline{s}) - \rho\,\boldsymbol{1})\cdot\underline{U} = \underline{0}``
with ``\underline{s} = \underline{s}_t + \xi\,\underline{\nu}`` is quadratic in
``\xi``:

```math
\bigl[\boldsymbol{Q} + \xi\,(\boldsymbol{R} + \boldsymbol{R}^{\mathsf T}) + \xi^2\,\boldsymbol{T} - \rho\,\boldsymbol{1}\bigr]\cdot\underline{U} = \underline{0},
```

```math
\boldsymbol{Q} = \underline{s}_t\cdot\mathbb{C}\cdot\underline{s}_t,\qquad
\boldsymbol{R} = \underline{s}_t\cdot\mathbb{C}\cdot\underline{\nu},\qquad
\boldsymbol{T} = \underline{\nu}\cdot\mathbb{C}\cdot\underline{\nu},
\qquad (\underline{a}\cdot\mathbb{C}\cdot\underline{b})_{ik} = a_j\,C_{ijkl}\,b_l .
```

Its determinant is a polynomial of degree six in ``\xi`` with real
coefficients: six roots, real or in conjugate pairs. Geometrically, they are
the intersections of the slowness surface with the line parallel to
``\underline{\nu}`` through ``\underline{s}_t``.

- A **real root** is a *homogeneous* wave, an ordinary plane wave. It carries
  energy along its group velocity, towards the medium 2 or away from it.
- A **complex root** is an *evanescent* wave: its amplitude varies as
  ``\exp(-\omega\operatorname{Im}\xi\;\underline{\nu}\cdot\underline{x})``. It
  decays away from the plane on the side where
  ``\operatorname{Im}\xi\;\underline{\nu}\cdot\underline{x} > 0`` and carries
  no energy across it (section 4).

Three of the six waves are *forward*: they carry energy, or decay, towards
``+\underline{\nu}``; the three others are *backward*. The direction of a
homogeneous wave is that of its energy flux, not of its slowness: in an
anisotropic medium a wave with ``\xi > 0`` may carry energy towards
``-\underline{\nu}``.

When ``|\underline{s}_t|`` grows, the line leaves a sheet of the slowness
surface: two real roots merge and become a conjugate pair. The incidence at
which this happens for a wave is its *critical angle*.

The roots are computed as the eigenvalues of a real ``6\times6`` matrix. The
*traction* that a wave of unit amplitude exerts on the plane is
``\mathrm{i}\omega\,\underline{\tau}`` with

```math
\underline{\tau} = (\underline{\nu}\cdot\mathbb{C}\cdot\underline{s})\cdot\underline{U}
= (\boldsymbol{R}^{\mathsf T} + \xi\,\boldsymbol{T})\cdot\underline{U} ,
```

since ``\boldsymbol{\sigma}\cdot\underline{\nu} = \mathbb{C}:\operatorname{grad}\underline{u}\cdot\underline{\nu}``
and ``\operatorname{grad}\underline{u} = \mathrm{i}\omega\,\underline{u}\otimes\underline{s}``.
Writing ``\tilde{\boldsymbol{Q}} = \boldsymbol{Q} - \rho\,\boldsymbol{1}``, the
quadratic problem is
``\tilde{\boldsymbol{Q}}\cdot\underline{U} + \xi\,\boldsymbol{R}\cdot\underline{U} + \xi\,\underline{\tau} = \underline{0}``,
and with ``\xi\,\underline{U} = \boldsymbol{T}^{-1}\cdot(\underline{\tau} - \boldsymbol{R}^{\mathsf T}\cdot\underline{U})``
from the definition of the traction, the pair
``\underline{\eta} = (\underline{U}, \underline{\tau})`` solves a linear
eigenvalue problem,

```math
\boldsymbol{N}\cdot\underline{\eta} = \xi\,\underline{\eta},
\qquad
\boldsymbol{N} =
\begin{pmatrix}
-\boldsymbol{T}^{-1}\boldsymbol{R}^{\mathsf T} & \boldsymbol{T}^{-1}\\
\boldsymbol{R}\,\boldsymbol{T}^{-1}\boldsymbol{R}^{\mathsf T} - \tilde{\boldsymbol{Q}} & -\boldsymbol{R}\,\boldsymbol{T}^{-1}
\end{pmatrix} .
```

[`interface_waves`](@ref) returns the six waves, sorted in forward and
backward ones.

## 3. Reflection and transmission

An incident wave, a forward wave of the medium 1, excites the three backward
waves of the medium 1, *reflected*, with amplitudes ``R_m``, and the three
forward waves of the medium 2, *transmitted*, with amplitudes ``T_m``: the
waves that go or decay away from the plane. A **welded** interface keeps the
displacement and the traction,

```math
\underline{U}^{\mathrm{inc}} + \sum_{m=1}^{3} R_m\,\underline{U}_m^{\mathrm{r}} = \sum_{m=1}^{3} T_m\,\underline{U}_m^{\mathrm{t}},
\qquad
\underline{\tau}^{\mathrm{inc}} + \sum_{m=1}^{3} R_m\,\underline{\tau}_m^{\mathrm{r}} = \sum_{m=1}^{3} T_m\,\underline{\tau}_m^{\mathrm{t}} ,
```

six equations for the six amplitudes, with unit polarizations. At normal
incidence on two isotropic media, a P wave gives
``R = (Z_1 - Z_2)/(Z_1 + Z_2)`` and ``T = 2Z_1/(Z_1 + Z_2)`` in absolute
value, with the impedances ``Z_i = \rho_i V_{P,i}``.

A thin compliant layer, a fracture or an imperfect bond is modeled by a
**linear-slip** interface [schoenberg1983, schoenberg1995](@cite): the traction stays
continuous and the displacement jumps in proportion to it,

```math
\underline{u}_2 - \underline{u}_1 = \boldsymbol{Z}\cdot\boldsymbol{\sigma}\cdot\underline{\nu} ,
```

with ``\boldsymbol{Z}`` the compliance of the interface (a length divided by
a stiffness). For the harmonic waves above, the first condition becomes

```math
\sum_m T_m\,\underline{U}_m^{\mathrm{t}} - \underline{U}^{\mathrm{inc}} - \sum_m R_m\,\underline{U}_m^{\mathrm{r}}
= \mathrm{i}\,\omega\boldsymbol{Z}\cdot\sum_m T_m\,\underline{\tau}_m^{\mathrm{t}} ,
```

in which the frequency enters through the product ``\omega\boldsymbol{Z}``
only: a slip interface reflects the high frequencies more. For a wave at
normal incidence between two identical media of impedance ``Z``, with the
compliance ``Z_N`` or ``Z_T`` of its polarization and
``\kappa = \omega Z_{N,T}\,Z``,

```math
|R|^2 = \frac{\kappa^2}{\kappa^2 + 4},\qquad |T|^2 = \frac{4}{\kappa^2 + 4} :
```

the two conditions are ``1 + R = T`` for the traction and
``T - (1 - R) = -\mathrm{i}\kappa\,T`` for the displacement, up to the signs
that the orientation of the polarizations fixes.

[`reflection_transmission`](@ref) solves the six equations, welded or with
`slip` ``= \omega\boldsymbol{Z}``.

## 4. Energy

The mean energy flux of a field through the plane, towards
``+\underline{\nu}``, is
``-\tfrac12\operatorname{Re}(\boldsymbol{\sigma}\cdot\underline{\nu}\cdot\overline{\dot{\underline{u}}})``.
For one wave of amplitude ``a`` it is
``\tfrac12\omega^2|a|^2\operatorname{Re}(\underline{\tau}\cdot\overline{\underline{U}})``;
for a sum of waves, products of two different waves appear as well. They
cancel. The matrix ``\boldsymbol{N}`` of section 2 is such that
``\boldsymbol{J}\boldsymbol{N}`` is symmetric, with ``\boldsymbol{J}`` the
exchange of the two blocks,

```math
\boldsymbol{J}\boldsymbol{N} =
\begin{pmatrix}
\boldsymbol{R}\,\boldsymbol{T}^{-1}\boldsymbol{R}^{\mathsf T} - \tilde{\boldsymbol{Q}} & -\boldsymbol{R}\,\boldsymbol{T}^{-1}\\
-\boldsymbol{T}^{-1}\boldsymbol{R}^{\mathsf T} & \boldsymbol{T}^{-1}
\end{pmatrix},
\qquad
\boldsymbol{J} = \begin{pmatrix}\boldsymbol{0} & \boldsymbol{1}\\ \boldsymbol{1} & \boldsymbol{0}\end{pmatrix} .
```

For two eigenvectors, ``\underline{\eta}_a\cdot\boldsymbol{J}\boldsymbol{N}\cdot\underline{\eta}_b``
equals both ``\xi_b\,\underline{\eta}_a\cdot\boldsymbol{J}\cdot\underline{\eta}_b``
and ``\xi_a\,\underline{\eta}_a\cdot\boldsymbol{J}\cdot\underline{\eta}_b``,
hence

```math
\underline{\tau}_a\cdot\underline{U}_b + \underline{U}_a\cdot\underline{\tau}_b = 0
\qquad\text{whenever } \xi_a \ne \xi_b .
```

The matrix being real, the conjugate of an eigenvector is the eigenvector of
the conjugate root. Applied to a wave ``a`` and to the conjugate of a wave
``b`` with ``\xi_a \ne \overline{\xi_b}``, the relation gives
``\underline{\tau}_a\cdot\overline{\underline{U}_b} = -\overline{\underline{\tau}_b\cdot\overline{\underline{U}_a}}``.
In the flux of the sum of the two waves, of amplitudes ``a_a`` and ``a_b``,
the mixed term is
``\operatorname{Re}(c\,\underline{\tau}_a\cdot\overline{\underline{U}_b} + \bar c\,\underline{\tau}_b\cdot\overline{\underline{U}_a})``
with ``c = a_a\overline{a_b}``, the real part of a number minus its conjugate:
it vanishes. In particular:

- an evanescent wave carries no mean flux through the plane
  (``a = b``, ``\xi_a \ne \overline{\xi_a}``);
- the fluxes of the waves of one side add, the homogeneous ones only.

When two real roots coincide (the two S waves of an isotropic medium), the
relation is not automatic: the two polarizations are then chosen in their
common plane so that it holds.

The displacement and the traction being continuous across a welded
interface, so is the flux of the total field, and the same holds with a real
symmetric ``\boldsymbol{Z}``, which stores energy without dissipating it. The
flux of the incident wave is therefore shared: with
``F = \operatorname{Re}(\underline{\tau}\cdot\overline{\underline{U}})`` for
each wave (negative for a reflected one),

```math
\sum_m |R_m|^2\,\frac{-F_m^{\mathrm{r}}}{F^{\mathrm{inc}}} + \sum_m |T_m|^2\,\frac{F_m^{\mathrm{t}}}{F^{\mathrm{inc}}} = 1 ,
```

each term being the fraction of the incident energy carried by one wave,
zero for an evanescent one.

## [5. Characteristic roots of a plane stiffness](@id th-characteristic-roots)

Without inertia, the same construction solves plane anisotropic elasticity.
A displacement of the plane
``\underline{u} = \underline{U}\,f(x_1 + \mu\,x_2)`` is in equilibrium in a
medium of stiffness ``\mathbb{C}`` when

```math
\begin{gathered}
\bigl[\boldsymbol{Q} + \mu\,(\boldsymbol{R} + \boldsymbol{R}^{\mathsf T}) + \mu^2\,\boldsymbol{T}\bigr]\cdot\underline{U} = \underline{0},
 \\
\boldsymbol{Q} = \underline{e}_1\cdot\mathbb{C}\cdot\underline{e}_1, \\
\boldsymbol{R} = \underline{e}_1\cdot\mathbb{C}\cdot\underline{e}_2, \\
\boldsymbol{T} = \underline{e}_2\cdot\mathbb{C}\cdot\underline{e}_2 ,
\end{gathered}
```

the problem of section 2 with ``\rho = 0``, ``\underline{s}_t = \underline{e}_1``
and ``\underline{\nu} = \underline{e}_2``, in dimension two. The matrix is
positive definite for a real ``\mu`` (it is the acoustic tensor of the
direction ``\underline{e}_1 + \mu\,\underline{e}_2``), so that the determinant,
a polynomial of degree four with real coefficients, has no real root: two
roots ``\mu_1``, ``\mu_2`` of positive imaginary part and their conjugates,
the *characteristic roots*. Every plane elastic field is a sum of two
functions of the complex variables ``z_k = x_1 + \mu_k x_2`` and of their
conjugates.

With the compliances of the plane in the engineering form,

```math
\begin{gathered}
\varepsilon_{11} = a_{11}\sigma_{11} + a_{12}\sigma_{22} + a_{16}\sigma_{12}, \\
\varepsilon_{22} = a_{12}\sigma_{11} + a_{22}\sigma_{22} + a_{26}\sigma_{12}, \\
2\varepsilon_{12} = a_{16}\sigma_{11} + a_{26}\sigma_{22} + a_{66}\sigma_{12} ,
\end{gathered}
```

the polynomial is, up to a factor, that of [lekhnitskii1981](@citet),

```math
a_{11}\,\mu^4 - 2a_{16}\,\mu^3 + (2a_{12} + a_{66})\,\mu^2 - 2a_{26}\,\mu + a_{22} = 0 .
```

It follows from the stress function ``F`` such that
``\sigma_{11} = \partial^2F/\partial x_2^2``,
``\sigma_{22} = \partial^2F/\partial x_1^2``,
``\sigma_{12} = -\partial^2F/\partial x_1\partial x_2``, which satisfies the
equilibrium: for ``F = \phi(x_1 + \mu x_2)`` the stresses are
``(\mu^2, 1, -\mu)\,\phi''``, and the compatibility
``\partial^2\varepsilon_{11}/\partial x_2^2 + \partial^2\varepsilon_{22}/\partial x_1^2 = 2\,\partial^2\varepsilon_{12}/\partial x_1\partial x_2``
of the strains gives the quartic.

In the axes of an orthotropic medium ``a_{16} = a_{26} = 0`` and the equation
is biquadratic. With ``s = \mu_1 + \mu_2`` and ``p = \mu_1\mu_2``, the quartic
divided by ``a_{11}`` is ``(\mu^2 - s\mu + p)(\mu^2 - \bar s\mu + \bar p)``,
whence

```math
p = -\sqrt{\frac{a_{22}}{a_{11}}},\qquad
s = \mathrm{i}\sqrt{\frac{2a_{12} + a_{66}}{a_{11}} + 2\sqrt{\frac{a_{22}}{a_{11}}}} .
```

An isotropic medium has ``a_{11} = a_{22}`` and ``2a_{12} + a_{66} = 2a_{11}``:
``p = -1``, ``s = 2\mathrm{i}``, the double root ``\mu_1 = \mu_2 = \mathrm{i}``,
and ``z_k = x_1 + \mathrm{i}x_2`` is the complex variable of the isotropic
plane. Near it the two roots are not differentiable functions of the
stiffness, while their sum and their product are: a computation that must
remain regular through isotropy is written with the pair ``(s, p)``
([`characteristic_pair`](@ref), [`characteristic_roots`](@ref)).

The stiffness seen by a small wave depends on the stress the medium is under:
it is the subject of [Third-order elasticity](@ref th-third-order).
