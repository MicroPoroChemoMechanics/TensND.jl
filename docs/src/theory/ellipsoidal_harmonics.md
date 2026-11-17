# [Ellipsoidal harmonics](@id th-ellipsoidal-harmonics)

This page defines the harmonics that separate in the coordinates of an
[`EllipsoidalChart`](@ref): the Lamé functions of the first and second kinds
and the products they form [dassios2012](@cite). They expand a harmonic field
inside and outside a triaxial ellipsoid, as the
[spheroidal harmonics](@ref th-spheroidal-harmonics) do for a spheroid.

## 1. Lamé functions

Let ``h_2 > h_3 > 0`` be the semi-focal distances of the chart,
``h_1^2 = h_2^2 - h_3^2``, and ``(\rho, \mu, \nu)`` the coordinates of a point,
``\nu^2 \le h_3^2 \le \mu^2 \le h_2^2 \le \rho^2``. A harmonic function separates
as ``E(\rho)\,E(\mu)\,E(\nu)``, ``E`` a Lamé function of degree ``n``: a solution of

```math
(x^2 - h_3^2)(x^2 - h_2^2)\,E'' + x\,(2x^2 - h_3^2 - h_2^2)\,E'
  + \big[(h_3^2 + h_2^2)\,p - n(n+1)\,x^2\big]\,E = 0
```

of the form ``\psi(x)\,P(x)``, with ``P`` a polynomial of the powers of one
parity and ``\psi`` one of ``1``, ``\sqrt{|x^2 - h_3^2|}``,
``\sqrt{|x^2 - h_2^2|}`` and their product: the classes K, L, M and N. Writing
``P = \sum_k c_k\,x^{d-2k}``, ``d`` the degree of ``P``, turns the equation into a
three-term recurrence on the ``c_k``; the recurrence stops at ``k = \lfloor d/2\rfloor``
only for the values of ``p`` that are the eigenvalues of a tridiagonal matrix.
Each class gives as many functions as its matrix has eigenvalues, ``2n + 1``
in all, normalized by ``c_0 = 1``; [`lame_functions`](@ref) returns them,
computed in floating point from the primal values of ``h_2`` and ``h_3``.

## 2. Interior and exterior harmonics

The products ``E(\rho)\,E(\mu)\,E(\nu)`` are polynomials of degree ``n`` in the
Cartesian coordinates, harmonic in the whole space once the square roots of
each class take the signs of the octant of the point: ``\nu`` carries the sign
of ``x_1``, ``\sqrt{h_3^2-\nu^2}`` that of ``x_2`` and ``\sqrt{h_2^2-\mu^2}`` that of
``x_3`` [dassios2012](@cite), which [`ellipsoidal_signed_coordinates`](@ref)
returns with the coordinates and [`lame_value`](@ref) uses. The harmonics that
vanish at infinity replace ``E(\rho)`` by the function of the second kind

```math
F(\rho) = (2n+1)\,E(\rho)\int_\rho^\infty
  \frac{\mathrm{d}t}{E(t)^2\sqrt{t^2-h_3^2}\,\sqrt{t^2-h_2^2}} ,
```

which decays as ``\rho^{-n-1}``. [`lame_second_kind`](@ref) computes it by
Gauss–Legendre quadrature in ``s = \rho/t``, the interval being split at
``s = 1/2`` and its upper half mapped by ``1 - s^2 = (\kappa/h_2)^2\sinh^2\tau``,
``\kappa = \sqrt{\rho^2 - h_2^2}``, which keeps the integrand smooth when
``\kappa`` is small, that is for an elongated ellipsoid.
