# [Thomsen's parameters](@id th-thomsen)

!!! tip "See also"
    **Manual**: [Plane waves](@ref man-plane-waves)

    **Tutorials**: [Thomsen's parameters and weak anisotropy](@ref tut-thomsen)

The five constants of a transversely isotropic stiffness do not tell at a
glance how its waves differ from those of an isotropic medium.
[thomsen1986](@citet) replaces them, with the density, by two velocities and
three dimensionless numbers that do. This page defines them, shows what each
one controls, and derives the velocities of a weakly anisotropic medium from
the exact ones of [Plane waves](@ref th-plane-waves).

## 1. Definitions

In a frame whose third vector is the axis, with the components of the
two-index notation and ``C_{66} = (C_{11} - C_{12})/2``:

```math
V_{P0} = \sqrt{\frac{C_{33}}{\rho}},\qquad
V_{S0} = \sqrt{\frac{C_{44}}{\rho}},\qquad
\varepsilon = \frac{C_{11} - C_{33}}{2C_{33}},\qquad
\gamma = \frac{C_{66} - C_{44}}{2C_{44}},
```

```math
\delta = \frac{(C_{13} + C_{44})^2 - (C_{33} - C_{44})^2}{2C_{33}\,(C_{33} - C_{44})} .
```

``V_{P0}`` and ``V_{S0}`` are the velocities of the P and S waves along the
axis. The definitions are exact: they hold for any strength of the
anisotropy, and they are inverted by

```math
\begin{gathered}
C_{33} = \rho V_{P0}^2, \\
C_{44} = \rho V_{S0}^2, \\
C_{11} = C_{33}(1 + 2\varepsilon), \\
C_{66} = C_{44}(1 + 2\gamma), \\
(C_{13} + C_{44})^2 = (C_{33} - C_{44})\,\bigl[C_{33}(1 + 2\delta) - C_{44}\bigr] .
\end{gathered}
```

The last relation fixes ``C_{13} + C_{44}`` up to its sign, which the
velocities do not see (it enters ``D`` by its square); the positive root is
kept. An isotropic medium has ``\varepsilon = \delta = \gamma = 0``.
[`thomsen_parameters`](@ref) reads them, and the form `:thomsen` of
[`ti_from_constants`](@ref) builds the stiffness.

## 2. What each parameter controls

Across the axis, the exact velocities of
[the transversely isotropic medium](@ref th-plane-waves) are

```math
V_P(90^\circ) = V_{P0}\sqrt{1 + 2\varepsilon},\qquad
V_{SH}(90^\circ) = V_{S0}\sqrt{1 + 2\gamma},\qquad
V_{SV}(90^\circ) = V_{S0} :
```

``\varepsilon`` and ``\gamma`` measure the difference between the horizontal
and the vertical velocities of the P and SH waves, for a vertical axis.

``\delta`` governs the P wave near the axis. Let
``f = 1 - C_{44}/C_{33} = 1 - V_{S0}^2/V_{P0}^2``. The definitions give
``C_{11} - C_{44} = C_{33}(f + 2\varepsilon)`` and
``(C_{13} + C_{44})^2 = C_{33}^2\,f\,(f + 2\delta)``, so that the quantity
``D`` of the exact velocities, with ``s = \sin\theta`` and ``c = \cos\theta``,
is

```math
\begin{aligned}
\frac{D}{C_{33}^2} &= \bigl[f\,(s^2 - c^2) + 2\varepsilon s^2\bigr]^2 + 4f\,(f + 2\delta)\,s^2c^2\\
&= f^2 + 8f\,\delta\,s^2c^2 + 4f\,\varepsilon\,s^2(s^2 - c^2) + 4\varepsilon^2 s^4 ,
\end{aligned}
```

an exact expression (``f^2(s^2 - c^2)^2 + 4f^2s^2c^2 = f^2``). Near the axis,
``s^2 \simeq \theta^2`` and ``c^2 \simeq 1``:
``\sqrt{D}/C_{33} = f + (4\delta - 2\varepsilon)\,\theta^2 + O(\theta^4)``, and

```math
\rho V_P^2 = C_{33}\,(1 + 2\delta\,\theta^2) + O(\theta^4),
\qquad
\rho V_{SV}^2 = C_{44} + 2C_{33}\,(\varepsilon - \delta)\,\theta^2 + O(\theta^4) ,
```

whatever the size of the parameters: the curvature of the P velocity at the
axis is ``\delta``, not ``\varepsilon``. A P wave reflected at a small offset
feels ``\delta`` only (section 4).

## 3. Weak anisotropy

To first order in ``\varepsilon`` and ``\delta``, the exact expression of
``D`` gives
``\sqrt{D}/C_{33} \simeq f + 4\delta\,s^2c^2 + 2\varepsilon\,s^2(s^2 - c^2)``.
The sum ``C_{11}s^2 + C_{33}c^2 + C_{44}`` being
``C_{33}(1 + 2\varepsilon s^2) + C_{44}``, the two roots are

```math
\rho V_P^2 \simeq C_{33}\,\bigl(1 + 2\delta\,s^2c^2 + 2\varepsilon\,s^4\bigr),
\qquad
\rho V_{SV}^2 \simeq C_{44} + 2C_{33}\,(\varepsilon - \delta)\,s^2c^2 ,
```

and the SH wave has
``\rho V_{SH}^2 = C_{44}(1 + 2\gamma\,s^2)`` exactly. Taking the square roots
to the same order,

```math
V_P(\theta) \simeq V_{P0}\,\bigl(1 + \delta\sin^2\theta\cos^2\theta + \varepsilon\sin^4\theta\bigr),
```

```math
V_{SV}(\theta) \simeq V_{S0}\Bigl[1 + \frac{V_{P0}^2}{V_{S0}^2}\,(\varepsilon - \delta)\sin^2\theta\cos^2\theta\Bigr],
\qquad
V_{SH}(\theta) \simeq V_{S0}\,\bigl(1 + \gamma\sin^2\theta\bigr) ,
```

the formulas of [thomsen1986](@citet), computed by [`weak_wave_modes`](@ref).
The SV velocity is the same along and across the axis, and its variation in
between is driven by ``\varepsilon - \delta``, magnified by
``V_{P0}^2/V_{S0}^2``.

Two related quantities are met in the literature. The first is

```math
\delta^* = \frac{2(C_{13} + C_{44})^2 - (C_{33} - C_{44})(C_{11} + C_{33} - 2C_{44})}{2C_{33}^2}
= f\,(2\delta - \varepsilon) ,
```

the equality being exact, by the two relations of section 2. The second is the
*anellipticity*

```math
\eta = \frac{\varepsilon - \delta}{1 + 2\delta} ,
```

which vanishes when the slowness curve of the P wave is an ellipse: for
``\delta = \varepsilon`` the exact ``D`` is the square
``C_{33}^2(f + 2\varepsilon s^2)^2``, hence
``\rho V_P^2 = C_{11}\sin^2\theta + C_{33}\cos^2\theta`` and
``\rho V_{SV}^2 = C_{44}`` ([`anellipticity`](@ref)).

## 4. Moveout velocity

A source and a receiver at the surface, a distance ``x`` apart, record the
wave reflected by a horizontal plane at the depth ``z`` in a medium of
vertical axis. For small offsets the travel time ``t`` follows
``t^2 \simeq t_0^2 + x^2/V_{\mathrm{nmo}}^2``, which defines the *moveout
velocity* ``V_{\mathrm{nmo}}``.

The energy travels along the ray, at the group velocity of magnitude ``g``
and angle ``\psi`` to the vertical: ``x = 2z\tan\psi`` and
``t = 2z/(g\cos\psi)``. Write ``g(\psi) \simeq g_0(1 + \tfrac12 g_2\psi^2)``
near the axis. Then ``x^2 \simeq 4z^2\psi^2`` and
``t^2 \simeq (4z^2/g_0^2)\,[1 + (1 - g_2)\,\psi^2]``, so that

```math
V_{\mathrm{nmo}}^2 = \frac{g_0^2}{1 - g_2} .
```

It remains to relate the group velocity to the phase velocity
``V(\theta) \simeq V_0(1 + \tfrac12 v_2\theta^2)``. In the plane of the axis
and of the wave normal, the group velocity has the component ``V`` along the
normal ([``\underline{g}\cdot\underline{n} = V``](@ref th-plane-waves)) and
``\mathrm{d}V/\mathrm{d}\theta`` across it: it is normal to the slowness curve
``\underline{n}(\theta)/V(\theta)``, whose tangent is
``\underline{n}'/V - \underline{n}\,V'/V^2`` with ``\underline{n}'`` the unit
vector across the normal. Hence
``g^2 = V^2 + V'^2`` and ``\tan(\psi - \theta) = V'/V``, that is
``\psi \simeq (1 + v_2)\,\theta`` and
``g \simeq V_0[1 + \tfrac12 v_2(1 + v_2)\,\theta^2]``: ``g_0 = V_0`` and
``g_2 = v_2/(1 + v_2)``. Finally

```math
V_{\mathrm{nmo}}^2 = V_0^2\,(1 + v_2),\qquad v_2 = \frac{V''(0)}{V(0)} .
```

The expansions of section 2 give ``v_2`` exactly, ``2\delta`` for P,
``2(V_{P0}/V_{S0})^2(\varepsilon - \delta)`` for SV and ``2\gamma`` for SH:

```math
\begin{gathered}
\bigl(V_{\mathrm{nmo}}^{P}\bigr)^2 = V_{P0}^2\,(1 + 2\delta), \\
\bigl(V_{\mathrm{nmo}}^{SV}\bigr)^2 = V_{S0}^2 + 2V_{P0}^2\,(\varepsilon - \delta), \\
\bigl(V_{\mathrm{nmo}}^{SH}\bigr)^2 = V_{S0}^2\,(1 + 2\gamma) ,
\end{gathered}
```

with no assumption on the strength of the anisotropy
([`squared_moveout_velocities`](@ref)). When ``\varepsilon - \delta`` is
large and negative, the square of the SV one is negative: the travel time
decreases with the offset near the axis, and no moveout velocity exists.

The next page, [Weak anisotropy of any symmetry](@ref th-weak-anisotropy),
extends the first-order velocities to a stiffness without any symmetry.
