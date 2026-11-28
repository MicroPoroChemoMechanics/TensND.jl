# [Plane waves](@id th-plane-waves)

!!! tip "See also"
    **Manual**: [Plane waves](@ref man-plane-waves)

    **Tutorials**: [Slowness and wave surfaces](@ref tut-wave-surfaces)

!!! note "Two-index components"
    ``C_{11}``, ``C_{44}``, ``C_{66}``, … are the two-index (Voigt) components
    ``C_{IJ}`` of the stiffness, ``I = 1, \dots, 6`` standing for the pairs 11,
    22, 33, 23, 13, 12 (``C_{44} = C_{2323}``, ``C_{66} = C_{1212}``), not the
    entries of the Kelvin–Mandel matrix, which carry factors ``\sqrt2`` and ``2``.
    Both matrices are written in full in [Contrast with Voigt](@ref th-voigt).

A homogeneous elastic medium carries, in each direction, three plane waves.
(The Christoffel *problem* of this page, named after E. B. Christoffel as the
[Christoffel symbols](@ref th-curvilinear) of a chart, has nothing else in
common with them.)
This page derives the eigenvalue problem that gives their velocities and
polarizations, the velocity at which each one carries energy, and the closed
forms of the isotropic and transversely isotropic classes.

## 1. The Christoffel problem

A medium of stiffness ``\mathbb{C}`` and mass density ``\rho`` moves, without
body force, according to

```math
\rho\,\ddot{\underline{u}} = \operatorname{div}\bigl(\mathbb{C}:\operatorname{grad}\underline{u}\bigr),
\qquad
\rho\,\ddot u_i = C_{ijkl}\,\frac{\partial^2 u_k}{\partial x_j\,\partial x_l} .
```

A *plane wave* is a displacement that depends on the position through its
projection on a unit vector ``\underline{n}`` only, and travels along it at the
velocity ``V``:

```math
\underline{u}(\underline{x},t) = \underline{U}\,f(\underline{n}\cdot\underline{x} - V t) .
```

``\underline{U}`` is the *polarization*, the direction in which the particles
move. Its second derivatives are
``\ddot{\underline{u}} = V^2\,\underline{U}\,f''`` and
``\partial^2 u_k/\partial x_j\partial x_l = n_j n_l\,U_k\,f''``, so that the
equation of motion becomes

```math
\bigl(\boldsymbol{\Gamma}(\underline{n}) - \rho V^2\,\boldsymbol{1}\bigr)\cdot\underline{U} = \underline{0},
\qquad
\Gamma_{ik}(\underline{n}) = n_j\,C_{ijkl}\,n_l ,
```

the *Christoffel problem*: ``\rho V^2`` is an eigenvalue of the *acoustic
tensor* ``\boldsymbol{\Gamma}(\underline{n}) = \underline{n}\cdot\mathbb{C}\cdot\underline{n}``
and ``\underline{U}`` the eigenvector that goes with it.

The major symmetry of ``\mathbb{C}`` makes ``\boldsymbol{\Gamma}`` symmetric. Its
minor symmetries give
``\underline{U}\cdot\boldsymbol{\Gamma}\cdot\underline{U} = (\underline{U}\stackrel{s}{\otimes}\underline{n}):\mathbb{C}:(\underline{U}\stackrel{s}{\otimes}\underline{n})``,
which is positive for a positive definite stiffness. The three eigenvalues are
therefore real and positive, and the three polarizations orthogonal: three
waves travel along ``\underline{n}``. In an isotropic medium one is polarized
along ``\underline{n}`` and two across it (section 4); in an anisotropic one
the polarizations are neither along nor across ``\underline{n}``, and the
waves are named by continuity, *quasi-longitudinal* (qP, the fastest) and
*quasi-transverse* (qS1 faster than qS2).

[`plane_waves`](@ref) solves the problem for any stiffness, the three waves
sorted by decreasing velocity.

## 2. Slowness

The *slowness vector* of a wave is ``\underline{s} = \underline{n}/V``: the
wave is ``\underline{U}\,g(\underline{s}\cdot\underline{x} - t)``. The acoustic
tensor being quadratic in its argument,
``\boldsymbol{\Gamma}(\underline{s}) = \boldsymbol{\Gamma}(\underline{n})/V^2``,
and the Christoffel problem reads

```math
\bigl(\boldsymbol{\Gamma}(\underline{s}) - \rho\,\boldsymbol{1}\bigr)\cdot\underline{U} = \underline{0} .
```

The slowness vectors of all the waves of the medium lie on the surface
``\det(\boldsymbol{\Gamma}(\underline{s}) - \rho\,\boldsymbol{1}) = 0``, of degree six
in ``\underline{s}``: the *slowness surface*, with three sheets, one per wave.
It is the natural object at an interface, where a component of the slowness is
imposed ([Waves at an interface](@ref th-interface-waves)).

The squared velocities ``v = V^2`` along ``\underline{n}`` are the roots of the
characteristic polynomial of ``\boldsymbol{\Gamma}/\rho``,

```math
v^3 - I_1\,v^2 + I_2\,v - I_3 = 0 ,
```

with ``I_1`` its trace, ``I_2`` the trace of its adjugate and ``I_3`` its
determinant, polynomials in the stiffness and in the direction
([`christoffel_invariants`](@ref)).

## 3. Group velocity

A harmonic wave of wave vector ``\underline{k} = k\,\underline{n}`` and angular
frequency ``\omega = kV`` is
``\underline{U}\exp[\mathrm{i}(\underline{k}\cdot\underline{x} - \omega t)]``.
With a unit polarization, the Christoffel problem multiplied by
``\underline{U}`` gives the *dispersion relation*

```math
\rho\,\omega^2 = \underline{U}\cdot\boldsymbol{\Gamma}(\underline{k})\cdot\underline{U} = C_{ijkl}\,U_i\,k_j\,U_k\,k_l .
```

A packet of waves of neighboring wave vectors travels at the *group velocity*
``\underline{g} = \partial\omega/\partial\underline{k}``, the gradient of
``\omega(\underline{k})``. In the dispersion relation both the wave vector and the
polarization depend on ``\underline{k}``: ``\underline{U}(\underline{k})`` is the
unit eigenvector of ``\boldsymbol{\Gamma}(\underline{k})`` of the mode followed,

```math
\boldsymbol{\Gamma}(\underline{k})\cdot\underline{U}(\underline{k}) = \rho\,\omega(\underline{k})^2\,\underline{U}(\underline{k}),
\qquad
\underline{U}\cdot\underline{U} = 1,
\qquad
\rho\,\omega^2 = C_{ijkl}\,U_i\,k_j\,U_k\,k_l .
```

The derivative of the last expression with respect to ``k_m`` has two kinds of
terms, those where ``\underline{U}`` is differentiated and those where the
``\underline{k}`` written in ``\boldsymbol{\Gamma}`` is:

```math
2\rho\,\omega\,\frac{\partial\omega}{\partial k_m}
  = \underbrace{C_{ijkl}\Big(\frac{\partial U_i}{\partial k_m}\,k_j\,U_k\,k_l + U_i\,k_j\,\frac{\partial U_k}{\partial k_m}\,k_l\Big)}_{\text{polarization}}
  + \underbrace{C_{ijkl}\,U_i\,U_k\,\big(\delta_{jm}\,k_l + k_j\,\delta_{lm}\big)}_{\text{wave vector in }\boldsymbol{\Gamma}} .
```

The two polarization terms are equal by the major symmetry, and their sum is
``2\,(\partial\underline{U}/\partial k_m)\cdot\boldsymbol{\Gamma}(\underline{k})\cdot\underline{U}``.
The eigenvalue equation replaces ``\boldsymbol{\Gamma}(\underline{k})\cdot\underline{U}``
by ``\rho\omega^2\underline{U}``, and the sum becomes
``2\rho\omega^2\,(\partial\underline{U}/\partial k_m)\cdot\underline{U} = \rho\omega^2\,\partial(\underline{U}\cdot\underline{U})/\partial k_m = 0``,
the polarization staying a unit vector. The change of the wave vector inside
``\boldsymbol{\Gamma}`` is thus all that remains: its two terms,
``C_{imkl}U_iU_kk_l`` and ``C_{ijkm}U_ik_jU_k``, are equal by the major
symmetry, and with ``\omega = kV`` and ``\underline{k} = k\,\underline{n}``

```math
\underline{g} = \frac{(\underline{U}\cdot\mathbb{C}\cdot\underline{U})\cdot\underline{n}}{\rho V},
\qquad
g_m = \frac{C_{imkl}\,U_i\,U_k\,n_l}{\rho V} .
```

It is the first-order perturbation of an eigenvalue: the eigenvalue of a
symmetric matrix that depends on a parameter changes as the matrix does along
its unit eigenvector, ``\delta(\rho\omega^2) = \underline{U}\cdot\delta\boldsymbol{\Gamma}\cdot\underline{U}``,
the change of the eigenvector contributing nothing at first order.

Three consequences follow.

- ``\underline{g}\cdot\underline{n} = \underline{U}\cdot\boldsymbol{\Gamma}(\underline{n})\cdot\underline{U}/(\rho V) = V``:
  the group velocity projects on the wave normal as the phase velocity. It is
  at least as large, and it is not along ``\underline{n}`` unless the medium
  is isotropic or the direction special.
- ``\underline{g}\cdot\underline{s} = 1``, and a variation
  ``\delta\underline{s}`` on the slowness surface keeps the frequency, so that
  ``\underline{g}\cdot\delta\underline{s} = 0``: the group velocity is normal
  to the slowness surface.
- The mean energy flux of the harmonic wave of amplitude ``a`` is
  ``-\tfrac12\operatorname{Re}(\boldsymbol{\sigma}\cdot\overline{\dot{\underline{u}}})``,
  of components ``\tfrac12\omega^2|a|^2\,C_{ijkl}\,U_i\,U_k\,n_l/V``, and its
  mean energy density ``\tfrac12\rho\,\omega^2|a|^2``, the kinetic and the
  elastic parts being equal. Their ratio is ``\underline{g}``: the energy
  travels at the group velocity.

The surface described by ``\underline{g}`` when ``\underline{n}`` covers the
sphere is the *wave surface*: the shape, after a unit time, of the disturbance
sent by a point source. A sheet of the slowness surface that is not convex
gives a folded wave surface, with several arrivals in some directions.

## 4. Isotropic medium

With ``\mathbb{C} = 3k\,\mathbb{J} + 2\mu\,\mathbb{K}`` and
``\lambda = k - 2\mu/3``, ``C_{ijkl} = \lambda\,\delta_{ij}\delta_{kl} + \mu(\delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk})``
and

```math
\boldsymbol{\Gamma}(\underline{n}) = (\lambda + \mu)\,\underline{n}\otimes\underline{n} + \mu\,\boldsymbol{1} .
```

One wave is polarized along ``\underline{n}``, the P wave, with
``\rho V_P^2 = \lambda + 2\mu``; every vector normal to ``\underline{n}`` is a
polarization of the S wave, ``\rho V_S^2 = \mu``. The velocities do not depend
on the direction and the group velocity is ``V\underline{n}``.

## 5. Transversely isotropic medium

Let ``\underline{a}`` be the axis of the medium and ``\theta`` the angle
between ``\underline{n}`` and ``\underline{a}``. In a frame whose third vector
is ``\underline{a}``, the stiffness has the five components
``C_{11}``, ``C_{33}``, ``C_{13}``, ``C_{44}`` and
``C_{66} = (C_{11} - C_{12})/2`` of the two-index notation
([The two-index matrix](@ref th-voigt)). Take the first vector in the plane
of ``\underline{a}`` and ``\underline{n}``, so that
``\underline{n} = (\sin\theta, 0, \cos\theta)``. The acoustic tensor is then

```math
\boldsymbol{\Gamma} =
\begin{pmatrix}
C_{11}\sin^2\theta + C_{44}\cos^2\theta & 0 & (C_{13} + C_{44})\sin\theta\cos\theta\\
0 & C_{66}\sin^2\theta + C_{44}\cos^2\theta & 0\\
(C_{13} + C_{44})\sin\theta\cos\theta & 0 & C_{44}\sin^2\theta + C_{33}\cos^2\theta
\end{pmatrix} .
```

The second direction, normal to the plane of ``\underline{a}`` and
``\underline{n}``, is a polarization: the *SH* wave,

```math
\rho V_{SH}^2 = C_{66}\sin^2\theta + C_{44}\cos^2\theta .
```

The two others are polarized in that plane and solve a quadratic equation:

```math
2\rho V_{qP,\,qSV}^2 = C_{11}\sin^2\theta + C_{33}\cos^2\theta + C_{44} \pm \sqrt{D} ,
```

```math
D = \bigl[(C_{11} - C_{44})\sin^2\theta - (C_{33} - C_{44})\cos^2\theta\bigr]^2 + 4\,(C_{13} + C_{44})^2\sin^2\theta\cos^2\theta .
```

Along the axis, ``\rho V_P^2 = C_{33}`` and the two S waves share
``\rho V_S^2 = C_{44}``; across it, ``\rho V_P^2 = C_{11}``,
``\rho V_{SV}^2 = C_{44}`` and ``\rho V_{SH}^2 = C_{66}``.

[`wave_modes`](@ref) returns these closed forms mode by mode, without
comparing a velocity with another, so that they hold for symbols; they are the
starting point of [Thomsen's parameters](@ref th-thomsen).
