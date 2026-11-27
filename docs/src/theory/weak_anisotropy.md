# [Weak anisotropy of any symmetry](@id th-weak-anisotropy)

!!! tip "See also"
    **Manual**: [Plane waves](@ref man-plane-waves)

    **Tutorials**: [Weak anisotropy of any symmetry](@ref tut-generalized-thomsen)

The velocities of a transversely isotropic medium have closed forms, from
which [Thomsen's first-order formulas](@ref th-thomsen) follow. A stiffness
without symmetry has none. This page obtains its velocities to first order in
the deviation from isotropy by a perturbation of the
[Christoffel problem](@ref th-plane-waves), and reads on the result the
parameters that generalize Thomsen's.

## 1. Perturbation of the Christoffel problem

Let ``\mathbb{C} = \mathbb{C}^0 + \delta\mathbb{C}`` with ``\mathbb{C}^0``
isotropic, of Lamé constants ``\lambda``, ``\mu``, and ``\delta\mathbb{C}``
small. The acoustic tensor along a unit vector ``\underline{n}`` is
``\boldsymbol{\Gamma} = \boldsymbol{\Gamma}^0 + \delta\boldsymbol{\Gamma}`` with

```math
\boldsymbol{\Gamma}^0 = (\lambda + \mu)\,\underline{n}\otimes\underline{n} + \mu\,\boldsymbol{1},
\qquad
\delta\boldsymbol{\Gamma} = \underline{n}\cdot\delta\mathbb{C}\cdot\underline{n} .
```

``\boldsymbol{\Gamma}^0`` has the simple eigenvalue ``\lambda + 2\mu`` on
``\underline{n}`` and the double eigenvalue ``\mu`` on the plane normal to
``\underline{n}``. The first-order perturbation of a symmetric matrix
[jech1989](@cite) treats them differently.

**The simple eigenvalue** moves by the mean value of the perturbation on its
eigenvector:

```math
\rho V_P^2 \simeq \lambda + 2\mu + \underline{n}\cdot\delta\boldsymbol{\Gamma}\cdot\underline{n}
= \underline{n}\cdot\boldsymbol{\Gamma}(\underline{n})\cdot\underline{n}
= (\underline{n}\otimes\underline{n}):\mathbb{C}:(\underline{n}\otimes\underline{n}) .
```

**The double eigenvalue** splits. Any vector of the plane is an eigenvector of
``\boldsymbol{\Gamma}^0``; the perturbation chooses among them. On an
orthonormal pair ``(\underline{e}_1, \underline{e}_2)`` normal to
``\underline{n}``, the two eigenvalues become those of the ``2\times2`` matrix
``\mu\,\delta_{ab} + \underline{e}_a\cdot\delta\boldsymbol{\Gamma}\cdot\underline{e}_b = \Gamma_{ab}``:

```math
\rho V_{S1,S2}^2 \simeq \tfrac12\Bigl[\Gamma_{11} + \Gamma_{22} \pm \sqrt{(\Gamma_{11} - \Gamma_{22})^2 + 4\,\Gamma_{12}^2}\Bigr],
\qquad
\Gamma_{ab} = \underline{e}_a\cdot\boldsymbol{\Gamma}(\underline{n})\cdot\underline{e}_b ,
```

and the two polarizations, to the lowest order, are its eigenvectors.

Both results are written with the acoustic tensor of ``\mathbb{C}`` itself:
the isotropic tensor ``\mathbb{C}^0`` has disappeared. The first-order
velocities do not depend on the reference medium, and they amount to this:
*keep the longitudinal direction and the transverse plane of an isotropic
medium, and restrict the Christoffel problem to each*. The neglected terms
are the couplings ``\underline{e}_a\cdot\boldsymbol{\Gamma}\cdot\underline{n}``
between the two, which enter the eigenvalues by their squares divided by
``\lambda + \mu``: the error is of second order. [`weak_plane_waves`](@ref)
computes these velocities; they take no comparison and hold for symbols.

## 2. The quartic form of the P wave

The first-order modulus of the P wave is a homogeneous polynomial of degree
four in the components of ``\underline{n}``. With the two-index components
``C_{IJ}`` ([The two-index matrix](@ref th-voigt)) and the vector
``\underline{m} = (n_1^2, n_2^2, n_3^2, 2n_2n_3, 2n_1n_3, 2n_1n_2)`` of the
components of ``\underline{n}\otimes\underline{n}`` with engineering shears,
``(\underline{n}\otimes\underline{n}):\mathbb{C}:(\underline{n}\otimes\underline{n}) = m_I\,C_{IJ}\,m_J``.
Subtracting ``M\,(n_1^2 + n_2^2 + n_3^2)^2 = M`` for a modulus ``M`` and
dividing by it,

```math
\begin{aligned}
\frac{\rho V_P^2}{M} \simeq 1 &+ 2\varepsilon_z^*\,n_3^4 + 2\varepsilon_y\,n_1^4 + 2\varepsilon_x\,n_2^4
 + 2\delta_z\,n_1^2n_2^2 + 2\delta_y\,n_1^2n_3^2 + 2\delta_x\,n_2^2n_3^2\\
&+ 4\bigl(\varepsilon_{15}\,n_1^3n_3 + \varepsilon_{16}\,n_1^3n_2 + \varepsilon_{24}\,n_2^3n_3 + \varepsilon_{26}\,n_2^3n_1 + \varepsilon_{34}\,n_3^3n_2 + \varepsilon_{35}\,n_3^3n_1\bigr)\\
&+ 4\bigl(\chi_x\,n_1^2n_2n_3 + \chi_y\,n_1n_2^2n_3 + \chi_z\,n_1n_2n_3^2\bigr) ,
\end{aligned}
```

an identity, in which the fifteen coefficients are the *generalized
parameters* of [mensch1997](@citet):

```math
\varepsilon_x = \frac{C_{22} - M}{2M},\quad
\varepsilon_y = \frac{C_{11} - M}{2M},\quad
\varepsilon_z^* = \frac{C_{33} - M}{2M},
```

```math
\delta_x = \frac{C_{23} + 2C_{44} - M}{M},\quad
\delta_y = \frac{C_{13} + 2C_{55} - M}{M},\quad
\delta_z = \frac{C_{12} + 2C_{66} - M}{M},
```

```math
\chi_x = \frac{C_{14} + 2C_{56}}{M},\quad
\chi_y = \frac{C_{25} + 2C_{46}}{M},\quad
\chi_z = \frac{C_{36} + 2C_{45}}{M},\qquad
\varepsilon_{IJ} = \frac{C_{IJ}}{M}\ \ (IJ = 15, 16, 24, 26, 34, 35) .
```

The index ``x``, ``y`` or ``z`` names the plane normal to that axis. Of the
twenty-one constants of the stiffness, the P wave sees to first order these
fifteen combinations only: the count of the fully symmetric tensors of order
four. With the modulus of the vertical P wave as reference, ``M = C_{33}``,
``\varepsilon_z^*`` vanishes and fourteen remain; it is the choice of
[`generalized_thomsen_parameters`](@ref).

In a symmetry plane the form reduces to Thomsen's. For
``\underline{n} = (\sin\theta, 0, \cos\theta)`` in a medium that has the
horizontal plane or the plane ``yz`` as a mirror (``C_{15} = C_{35} = 0``),
and ``M = C_{33}``,

```math
V_P(\theta) \simeq V_{P0}\,\bigl(1 + \delta_y\sin^2\theta\cos^2\theta + \varepsilon_y\sin^4\theta\bigr) .
```

For a transversely isotropic medium of axis ``\underline{e}_3``,
``\varepsilon_x = \varepsilon_y = \varepsilon`` and
``\delta_x = \delta_y = (C_{13} + 2C_{44} - C_{33})/C_{33}``, which is not
Thomsen's ``\delta`` but agrees with it to first order [sayers1994](@cite):
with ``a = C_{13} + C_{44}`` and ``b = C_{33} - C_{44}``, ``a - b = C_{33}\,\delta_x`` and

```math
\delta = \frac{a^2 - b^2}{2C_{33}\,b} = \delta_x + \frac{\delta_x^2}{2\,(1 - C_{44}/C_{33})} .
```

## [3. Shear-wave splitting](@id th-splitting)

The two S waves of an isotropic medium have the same velocity; a weak
anisotropy separates them, by

```math
\rho\,(V_{S1}^2 - V_{S2}^2) \simeq \sqrt{(\Gamma_{11} - \Gamma_{22})^2 + 4\,\Gamma_{12}^2} ,
```

a quantity of first order which vanishes in the directions where the
``2\times2`` problem is a multiple of the identity. A shear wave that enters
such a medium splits into a fast and a slow wave, polarized along the two
eigenvectors, and the two arrive after a path of length ``L`` with the delay
``L\,(1/V_{S2} - 1/V_{S1})``. [`shear_splitting`](@ref) returns the relative
difference ``(V_{S1} - V_{S2})/V_{S1}`` of the exact velocities, the fast
polarization and the delay per unit length.

Along the vertical, ``\Gamma_{11} = C_{55}``, ``\Gamma_{22} = C_{44}`` and
``\Gamma_{12} = C_{45}``. In the frame where ``C_{45} = 0`` the two vertical S
waves are polarized along ``\underline{e}_1`` and ``\underline{e}_2``, with
the moduli ``C_{55}`` and ``C_{44}``, the references of the three parameters
of the S waves [mensch1997](@cite),

```math
\gamma_x = \frac{C_{44} - C_{55}}{2C_{55}},\qquad
\gamma_y = \frac{C_{66} - C_{44}}{2C_{44}},\qquad
\gamma_z = \frac{C_{55} - C_{44}}{2C_{44}} .
```

Writing ``V_{kl}`` for the velocity of the S wave that travels along
``\underline{e}_k`` polarized along ``\underline{e}_l``
(``\rho V_{xz}^2 = C_{55}``, ``\rho V_{yz}^2 = C_{44}``,
``\rho V_{xy}^2 = C_{66}``), they are to first order
``\gamma_x \simeq (V_{yz} - V_{xz})/V_{xz}``,
``\gamma_y \simeq (V_{xy} - V_{zy})/V_{zy}`` and
``\gamma_z \simeq (V_{xz} - V_{yz})/V_{yz}``: ``\gamma_z`` is the splitting of
the vertical waves, and ``\gamma_y`` reduces to Thomsen's ``\gamma`` for a
transversely isotropic medium of vertical axis, where
``\gamma_x = \gamma_z = 0``.

The waves considered so far travel in an unbounded medium. The next page,
[Waves at an interface](@ref th-interface-waves), follows them across a plane
that separates two media.
