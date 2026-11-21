# [Toroidal functions](@id man-toroidal-harmonics)

!!! tip "See also"
    **Theory**: [Toroidal coordinates](@ref th-toroidal) · [Toroidal functions](@ref th-toroidal-harmonics) · [Curvilinear differential calculus](@ref th-curvilinear)

    **Manual**: [Coordinate systems (symbolic)](@ref man-coorsystems) · [Charts and coordinate conversions](@ref man-charts)

How to evaluate the functions on which the harmonic functions of a toroidal
chart are built. The definitions and the method are on
[Toroidal functions](@ref th-toroidal-harmonics).

## The functions and their degrees

Every function returns the values for all the degrees up to `N` at once: the
entry `n + 1` is the function of degree ``n-\tfrac12``, at the order `m` and at
the coordinate ``\alpha``.

```@example tor
using TensND

toroidal_P(3, 0, 1.0)          # P₋₁⸝₂, P₁⸝₂, P₃⸝₂, P₅⸝₂ at cosh α, α = 1
```

```@example tor
toroidal_Q(3, 0, 1.0)
```

The first kind is bounded on the axis of the chart, ``\alpha=0``, where every
degree equals 1, and grows toward the focal circle; the second kind is singular
on the axis and decays toward the focal circle.

```@example tor
toroidal_P(3, 0, 0.0)
```

The order is the second argument:

```@example tor
toroidal_Q(3, 2, 1.0)
```

## Derivatives

[`toroidal_dP`](@ref) and [`toroidal_dQ`](@ref) return the derivatives with
respect to ``\alpha``, which a condition on a torus ``\alpha=\alpha_0`` needs.

```@example tor
toroidal_dQ(3, 0, 1.0)
```

They agree with the derivative of the values by `ForwardDiff`, which the
functions accept:

```@example tor
using ForwardDiff
ForwardDiff.derivative(α -> toroidal_Q(3, 0, α), 1.0)
```

## A harmonic function of the chart

The product of ``\sqrt{\cosh\alpha-\cos\beta}``, of a toroidal function and of
a sine or cosine of ``n\beta`` and of ``m\varphi`` is harmonic. Its Laplacian,
computed by the pointwise toroidal system, vanishes to rounding:

```@example tor
CS = coorsys_toroidal(PointwiseGeometry, 1.0)
n, m = 2, 1
h(ξ) = sqrt(cosh(ξ[2]) - cos(ξ[3])) * toroidal_Q(n, m, ξ[2])[n + 1] * cos(n * ξ[3]) * cos(m * ξ[1])
ξ = [0.3, 0.8, 1.1]          # (φ, α, β)
h(ξ), laplacian(h, CS)(ξ)
```

## The constant as a series

Outside the axis, the constant 1 is the series of the functions of the second
kind

```math
1=\frac{\sqrt2}{\pi}\,\sqrt{\cosh\alpha-\cos\beta}\,\sum_{n\ge0}\epsilon_n\,Q_{n-1/2}(\cosh\alpha)\,\cos n\beta,
\qquad \epsilon_0=1,\ \epsilon_{n\ge1}=2 ,
```

whose terms decay as ``\mathrm{e}^{-n\alpha}``:

```@example tor
α, β, N = 0.8, 1.1, 40
Q = toroidal_Q(N, 0, α)
sqrt(2) / π * sqrt(cosh(α) - cos(β)) * sum((n == 0 ? 1 : 2) * Q[n + 1] * cos(n * β) for n in 0:N)
```

## Number types

The functions take numbers and dual numbers. They are computed by iterations
that stop on a value, and refuse a symbol:

```@example tor
using Symbolics
@variables a
try
    toroidal_Q(2, 0, a)
catch err
    println(sprint(showerror, err))
end
```
