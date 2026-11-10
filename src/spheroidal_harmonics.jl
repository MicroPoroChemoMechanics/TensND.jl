# ── Spheroidal harmonics and their re-expansion ──────────────────────────────
#
# The harmonic functions that separate in a spheroidal chart are products of
# associated Legendre functions of the two meridian coordinates. This file
# evaluates them and re-expands a harmonic of one chart in the harmonics of a
# coaxial chart of another semi-focal distance, prolate, oblate or spherical.
#
# Conventions. With `P_n` the Legendre polynomial and `Q_0(x) = arccoth x`,
#
#     P_n^m(x) = (-√(1-x²))^m dᵐP_n/dxᵐ           (|x| ≤ 1)
#     P_n^m(x) = ( √(x²-1))^m dᵐP_n/dxᵐ           (x > 1)
#     Q_n^m(x) = ( √(x²-1))^m dᵐQ_n/dxᵐ           (x > 1)
#
# and, on the oblate side, the real functions of `τ > 0`
#
#     𝒫_n^m(τ) = i⁻ⁿ P_n^m(iτ),    𝒬_n^m(τ) = iⁿ⁺¹ Q_n^m(iτ),    √((iτ)²-1) = i √(τ²+1).
#
# The interior harmonic of degree `n` and order `m ≥ 0` is
# `P_n^m(p) 𝖯_n^m(ζ) e^{imϕ}`, regular on the focal set; the exterior one is
# `P_n^m(p) 𝖰_n^m(ζ) e^{imϕ}`, which vanishes at infinity; `(𝖯, 𝖰, ζ)` stands
# for `(P, Q, q)` or `(𝒫, 𝒬, τ)`. In the spherical chart they are the solid
# harmonics `rⁿ P_n^m(cos θ) e^{imϕ}` and `r⁻ⁿ⁻¹ P_n^m(cos θ) e^{imϕ}`.

function _check_degree_order(N::Integer, m::Integer)
    m ≥ 0 || throw(DomainError(m, "the order m must be non-negative"))
    N ≥ 0 || throw(DomainError(N, "the degree must be non-negative"))
    return nothing
end

# `(n-m+1) F_{n+1} = (2n+1) x F_n + σ (n+m) F_{n-1}` upward from `F_m^m`, into
# `v[n+1]`, `m ≤ n ≤ N`. Each step reads the two previous ones: a recurrence.
function _legendre_forward(N::Integer, m::Integer, x, Fmm, σ::Integer)
    v = zeros(typeof(Fmm), N + 1)
    m > N && return v
    v[m + 1] = Fmm
    m < N && (v[m + 2] = (2m + 1) * x * Fmm)
    for n in (m + 1):(N - 1)
        v[n + 2] = ((2n + 1) * x * v[n + 1] + σ * (n + m) * v[n]) / (n - m + 1)
    end
    return v
end

# `(2m-1)!! sᵐ`, the value of `F_m^m` up to its sign.
_double_factorial_power(s, m::Integer) = prod((2k - 1) * s for k in 1:m; init = one(s))

"""
    legendre_P(N, m, x) -> Vector

Associated Legendre functions of the first kind on the cut, `|x| ≤ 1`:
`v[n+1] = P_n^m(x) = (-√(1-x²))^m dᵐP_n/dxᵐ` for `m ≤ n ≤ N`, and zero for
`n < m`. The upward recurrence in the degree is stable there. The argument may
be a number, a dual or a symbol.

# Examples
```jldoctest
julia> legendre_P(2, 1, 0.5) ≈ [0, -√0.75, -1.5 * √0.75]
true
```

See also [`radial_P`](@ref), [`radial_Q`](@ref).
"""
function legendre_P(N::Integer, m::Integer, x)
    _check_degree_order(N, m)
    s = sqrt((1 - x) * (1 + x))
    return _legendre_forward(N, m, x, _double_factorial_power(-s, m), -1)
end

"""
    radial_P(chart, N, m, ζ) -> Vector

Radial factor of the interior harmonics of `chart`, `v[n+1]` for `m ≤ n ≤ N`
(zero for `n < m`):

- [`ProlateSpheroidalChart`](@ref): `P_n^m(q) = (√(q²-1))^m dᵐP_n/dqᵐ`, `q ≥ 1`;
- [`OblateSpheroidalChart`](@ref): `𝒫_n^m(τ) = i⁻ⁿ P_n^m(iτ)`, `τ ≥ 0`, real;
- [`SphericalChart`](@ref): `rⁿ`.

These are the dominant solutions of their recurrences, which are run upward.
The argument may be a number, a dual or a symbol.

See also [`radial_Q`](@ref), [`interior_harmonics`](@ref).
"""
function radial_P(chart::_SpheroidalChart, N::Integer, m::Integer, ζ)
    _check_degree_order(N, m)
    s = _spheroidal_sign(chart)
    # `𝒫` obeys the recurrence of `P` with the sign of its last term reversed.
    return _legendre_forward(N, m, ζ, _double_factorial_power(sqrt(ζ^2 - s), m), -s)
end
function radial_P(::SphericalChart, N::Integer, m::Integer, r)
    _check_degree_order(N, m)
    return [n < m ? zero(r) : r^n for n in 0:N]
end

# `e^{-ξ}` with `ξ = arccosh q` or `arcsinh τ`: the ratio `𝖰_n/𝖰_{n-1}` at large
# degree, and the rate at which the exterior functions decay with the degree.
_decay_ratio(chart::_SpheroidalChart, ζ) = inv(ζ + sqrt(ζ^2 - _spheroidal_sign(chart)))

"""
    radial_Q(chart, N, m, ζ) -> Vector

Radial factor of the exterior harmonics of `chart`, `v[n+1]` for `m ≤ n ≤ N`
(zero for `n < m`):

- [`ProlateSpheroidalChart`](@ref): `Q_n^m(q) = (√(q²-1))^m dᵐQ_n/dqᵐ`, `q > 1`;
- [`OblateSpheroidalChart`](@ref): `𝒬_n^m(τ) = iⁿ⁺¹ Q_n^m(iτ)`, `τ > 0`, real;
- [`SphericalChart`](@ref): `r⁻ⁿ⁻¹`.

`𝖰_n^m` is the minimal solution of its recurrence: run upward, the recurrence
loses it to rounding. The ratios `𝖰_n^m/𝖰_{n-1}^m` are therefore computed
downward, from their limit `e^{-ξ}` far above `N`, and the first value follows
from the Casoratian with the dominant solution at the lowest degree,

    𝖯_{m+1}^m 𝖰_m^m - s 𝖯_m^m 𝖰_{m+1}^m = (-1)^m (2m)!,

with `s = 1` (prolate) or `s = -1` (oblate). The result is accurate to rounding
at every degree and differentiable (duals in `ζ`). It is singular on the focal
set (`q = 1`, `τ = 0`).

See also [`radial_P`](@ref), [`exterior_harmonics`](@ref).
"""
function radial_Q(chart::_SpheroidalChart, N::Integer, m::Integer, ζ)
    _check_degree_order(N, m)
    s = _spheroidal_sign(chart)
    ratio = _decay_ratio(chart, ζ)
    T = typeof(ratio)
    # The downward run is an iteration sized on a value: numbers and duals
    # only, a symbol has no such size.
    ξ = -log(primal(ratio))
    ξ > zero(ξ) || throw(
        DomainError(ζ, "the exterior harmonics are singular on the focal set of the chart")
    )
    # A start `L` degrees above `N` leaves an error `e^{-2Lξ}` on the ratio at
    # `N`: `L` is taken so that it falls below the rounding of the number type.
    L = ceil(Int, precision(typeof(ξ)) * log(2) / (2 * ξ)) + 20
    v = zeros(T, N + 1)                      # holds 𝖰_n/𝖰_{n-1} first, for m < n ≤ N
    for n in (N + 1 + L):-1:(m + 1)
        ratio = (n + m) / ((2n + 1) * ζ - s * (n - m + 1) * ratio)
        n ≤ N && (v[n + 1] = ratio)
        n == m + 1 && (v[m + 1] = ratio)
    end
    # The Casoratian at the lowest degree fixes 𝖰_m, the ratios give the rest:
    # the dominant solution is needed at the degrees m and m+1 only, where it
    # cannot overflow.
    P = radial_P(chart, m + 1, m, ζ)
    v[m + 1] = (-1)^m * _factorial(T, 2m) / (P[m + 2] - s * P[m + 1] * v[m + 1])
    for n in (m + 1):N                       # a recurrence
        v[n + 1] *= v[n]
    end
    return v
end
function radial_Q(::SphericalChart, N::Integer, m::Integer, r)
    _check_degree_order(N, m)
    return [n < m ? zero(inv(r)) : inv(r)^(n + 1) for n in 0:N]
end

# The two meridian coordinates of a point, angular then radial, as the
# harmonics read them: `(p, ζ)`, or `(cos θ, r)` in the spherical chart.
_meridian(::_SpheroidalChart, ξ) = (ξ[2], ξ[3])
_meridian(::SphericalChart, ξ) = (cos(ξ[1]), ξ[3])

"""
    interior_harmonics(chart, N, m, ξ) -> Vector

Meridian part of the interior harmonics of order `m ≥ 0` of `chart` at the
point of coordinates `ξ`: `v[n+1] = P_n^m(p) 𝖯_n^m(ζ)` for `m ≤ n ≤ N`, or
`rⁿ P_n^m(cos θ)` in the [`SphericalChart`](@ref). Multiplied by `cos mϕ` or
`sin mϕ`, each is harmonic and regular on the focal set of the chart.

See also [`exterior_harmonics`](@ref), [`interior_reexpansion`](@ref).
"""
function interior_harmonics(chart::AbstractChart{3}, N::Integer, m::Integer, ξ)
    p, ζ = _meridian(chart, ξ)
    return legendre_P(N, m, p) .* radial_P(chart, N, m, ζ)
end

"""
    exterior_harmonics(chart, N, m, ξ) -> Vector

Meridian part of the exterior harmonics of order `m ≥ 0` of `chart` at the
point of coordinates `ξ`: `v[n+1] = P_n^m(p) 𝖰_n^m(ζ)` for `m ≤ n ≤ N`, or
`r⁻ⁿ⁻¹ P_n^m(cos θ)` in the [`SphericalChart`](@ref). Multiplied by `cos mϕ` or
`sin mϕ`, each is harmonic off the focal set of the chart and vanishes at
infinity.

See also [`interior_harmonics`](@ref), [`exterior_reexpansion`](@ref).
"""
function exterior_harmonics(chart::AbstractChart{3}, N::Integer, m::Integer, ξ)
    p, ζ = _meridian(chart, ξ)
    return legendre_P(N, m, p) .* radial_Q(chart, N, m, ζ)
end

# ── Re-expansion coefficients ────────────────────────────────────────────────

# The signed square of the semi-focal distance, `c²` (prolate) or `-c²`
# (oblate): the only parameter of a chart that the re-expansion sees. The
# formulas of one case are those of the other continued in this parameter.
_signed_focal_square(chart::_SpheroidalChart) = _spheroidal_sign(chart) * chart.c^2

"""
    reexpansion_polynomials(n, imax, x) -> Vector

The polynomials of the re-expansion of spheroidal harmonics,
`Φ[i+1] = Φᵢ⁽ⁿ⁾(x)` for `0 ≤ i ≤ imax`:

    Φᵢ⁽ⁿ⁾(x) = (2n+4i+1)/(2n+1) · (n+½)ᵢ/i! · ₂F₁(-i, n+i+½; n+3/2; x)
             = (2n+4i+1)/(2n+2i+1) · Pᵢ^(n+½, -1)(1-2x),

with `Pᵢ^(α, β)` the Jacobi polynomial: `Φ₀ = 1`, `Φ₁ = (2n+5)(1-x)/2`, and
`Φᵢ(1) = 0` for `i ≥ 1`. They are evaluated by the three-term recurrence of
the Jacobi polynomials; the hypergeometric sum cancels for `0 < x < 1`.
"""
function reexpansion_polynomials(n::Integer, imax::Integer, x)
    T = typeof(x / 2)
    α = T(2n + 1) / 2
    σ = α - 1                                # α + β, with β = -1
    y = 1 - 2 * x
    J = zeros(T, imax + 1)
    J[1] = one(T)
    imax ≥ 1 && (J[2] = (α + 1) * (1 - x))
    for i in 1:(imax - 1)                    # a recurrence
        A = (2i + σ + 1) * (2i + σ + 2) / (2 * (i + 1) * (i + σ + 1))
        B = (α^2 - 1) * (2i + σ + 1) / (2 * (i + 1) * (i + σ + 1) * (2i + σ))
        C = (i + α) * (i - 1) * (2i + σ + 2) / ((i + 1) * (i + σ + 1) * (2i + σ))
        J[i + 2] = (A * y + B) * J[i + 1] - C * J[i]
    end
    return [T(2n + 4i + 1) / (2n + 2i + 1) * J[i + 1] for i in 0:imax]
end

# `g_a/g_b` with `g_n = (n+m)!/(n-m)!`, for `a ≤ b`.
_order_ratio(::Type{T}, a::Integer, b::Integer, m::Integer) where {T} =
    inv(prod(T(k + m) / (k - m) for k in (a + 1):b; init = one(T)))

_factorial(::Type{T}, k::Integer) where {T} = prod(T, 1:k; init = one(T))
# `(a)ₖ/k!`, factor by factor: neither the Pochhammer symbol nor the factorial
# is formed, which would overflow long before their ratio does.
_pochhammer_ratio(a, k::Integer) = prod((a + j) / (j + 1) for j in 0:(k - 1); init = one(a))
# `(-1)^m (n+m)!/(2n+1)!!`: the factor of `r^{-n-1} P_n^m(cos θ)` in
# `c^{-n-1} P_n^m(p) Q_n^m(q)` far from the foci.
_far_field_factor(::Type{T}, n::Integer, m::Integer) where {T} =
    (-1)^m * prod(T(k) / (2k + 1) for k in 1:n; init = one(T)) * prod(T, (n + 1):(n + m); init = one(T))
# `(2n)!/(2ⁿ n! (n-m)!) = (2n-1)!!/(n-m)!`, the leading coefficient of `dᵐP_n/dxᵐ`.
_leading_coefficient(::Type{T}, n::Integer, m::Integer) where {T} =
    prod(T(2k - 1) / k for k in 1:(n - m); init = one(T)) * prod(T(2k - 1) for k in (n - m + 1):n; init = one(T))

# A sequence from its first term and the ratio of two consecutive ones.
function _from_ratios(first, ratio::F, kmax::Integer) where {F}
    v = fill(first, kmax + 1)
    for k in 0:(kmax - 1)                    # a recurrence
        v[k + 2] = v[k + 1] * ratio(k)
    end
    return v
end

"""
    exterior_reexpansion(to, from, n, m, imax) -> Vector

Coefficients `K[i+1]`, `0 ≤ i ≤ imax`, of the exterior harmonic of degree `n`
and order `m` of the chart `from` in the exterior harmonics of the coaxial
chart `to`:

    P_n^m(p₁) 𝖰_n^m(ζ₁) = Σᵢ Kᵢ P_{n+2i}^m(p₂) 𝖰_{n+2i}^m(ζ₂).

Between two spheroidal charts of semi-focal distances `c₁`, `c₂` and signs
`s₁`, `s₂` (`1` prolate, `-1` oblate),

    Kᵢ = (-s₂)ⁱ (c₁/c₂)ⁿ⁺¹ · (n+m)! (n+2i-m)! / ((n-m)! (n+2i+m)!) · Φᵢ⁽ⁿ⁾(s₁c₁²/(s₂c₂²)),

with `Φᵢ⁽ⁿ⁾` of [`reexpansion_polynomials`](@ref). A [`SphericalChart`](@ref)
on either side gives the expansion in, or of, the solid harmonic
`r⁻ⁿ⁻¹ P_n^m(cos θ)`. The series converges outside the smallest coordinate
surface of `to` that encloses the focal set of `from`, see
[`reexpansion_bound`](@ref).

# Examples
```jldoctest
julia> from, to = ProlateSpheroidalChart(1.0), ProlateSpheroidalChart(1.5); x = [0.7, -0.4, 2.1];

julia> K = exterior_reexpansion(to, from, 2, 1, 40);

julia> lhs = exterior_harmonics(from, 2, 1, cartesian_to_coords(from, x))[3];

julia> rhs = sum(K .* exterior_harmonics(to, 82, 1, cartesian_to_coords(to, x))[3:2:end]);

julia> isapprox(lhs, rhs; rtol = 1e-12)
true
```

See also [`interior_reexpansion`](@ref).
"""
function exterior_reexpansion(
        to::_SpheroidalChart, from::_SpheroidalChart, n::Integer, m::Integer, imax::Integer
    )
    _check_reexpansion(n, m)
    Φ = reexpansion_polynomials(n, imax, _signed_focal_square(from) / _signed_focal_square(to))
    T = eltype(Φ)
    λ, s = from.c / to.c, _spheroidal_sign(to)
    return [(-s)^i * λ^(n + 1) * _order_ratio(T, n, n + 2i, m) * Φ[i + 1] for i in 0:imax]
end
# Spheroidal exterior harmonic in the solid harmonics `r^{-l-1} P_l^m`, `r > c`:
#
#     aₖ = sᵏ c^{n+2k+1} (-1)^m (n+m)!/(n-m)! · 2ⁿ (n+2k-m)! (n+k)! / (k! (2n+2k+1)!).
function exterior_reexpansion(
        ::SphericalChart, from::_SpheroidalChart, n::Integer, m::Integer, imax::Integer
    )
    _check_reexpansion(n, m)
    c, s = from.c, _spheroidal_sign(from)
    T = typeof(c / 2)
    ratio(k) = s * c^2 * T((n + 2k - m + 1) * (n + 2k - m + 2)) / (2 * (k + 1) * (2n + 2k + 3))
    return _from_ratios(c^(n + 1) * _far_field_factor(T, n, m), ratio, imax)
end
# Solid harmonic `r^{-n-1} P_n^m` in the spheroidal exterior harmonics.
function exterior_reexpansion(
        to::_SpheroidalChart, ::SphericalChart, n::Integer, m::Integer, imax::Integer
    )
    _check_reexpansion(n, m)
    c, s = to.c, _spheroidal_sign(to)
    T = typeof(c / 2)
    e = c^(n + 1) * _far_field_factor(T, n, m)
    # `Φᵢ⁽ⁿ⁾(0)`, the value of the polynomial when the source has no focal distance.
    Φ₀(i) = T(2n + 4i + 1) / (2n + 1) * _pochhammer_ratio(T(2n + 1) / 2, i)
    return [(-s)^i * Φ₀(i) * _order_ratio(T, n, n + 2i, m) / e for i in 0:imax]
end

"""
    interior_reexpansion(to, from, n, m) -> Vector

Coefficients `L[k+1]`, `0 ≤ k ≤ ⌊(n-m)/2⌋`, of the interior harmonic of degree
`n` and order `m` of the chart `from` in the interior harmonics of the coaxial
chart `to`, of degrees `j = n - 2k`:

    P_n^m(p₁) 𝖯_n^m(ζ₁) = Σₖ Lₖ P_j^m(p₂) 𝖯_j^m(ζ₂).

The sum is finite and the identity holds at every point: both sides are the
same polynomial of the Cartesian coordinates. Between two spheroidal charts,

    Lₖ = (-s₁)ᵏ (c₂/c₁)ʲ · (n+m)! (j-m)! / ((n-m)! (j+m)!) · (2j+1)/(2n+1) · Φₖ⁽ʲ⁾(s₂c₂²/(s₁c₁²)),

with `Φₖ⁽ʲ⁾` of [`reexpansion_polynomials`](@ref). A [`SphericalChart`](@ref)
on either side gives the expansion in, or of, the solid harmonic
`rⁿ P_n^m(cos θ)`.

# Examples
```jldoctest
julia> from, to = OblateSpheroidalChart(1.0), ProlateSpheroidalChart(0.6); x = [0.7, -0.4, 0.3];

julia> L = interior_reexpansion(to, from, 4, 2);

julia> lhs = interior_harmonics(from, 4, 2, cartesian_to_coords(from, x))[5];

julia> rhs = sum(L .* interior_harmonics(to, 4, 2, cartesian_to_coords(to, x))[5:-2:3]);

julia> lhs ≈ rhs
true
```

See also [`exterior_reexpansion`](@ref).
"""
function interior_reexpansion(to::_SpheroidalChart, from::_SpheroidalChart, n::Integer, m::Integer)
    _check_reexpansion(n, m)
    x = _signed_focal_square(to) / _signed_focal_square(from)
    T = typeof(x / 2)
    μ, s = to.c / from.c, _spheroidal_sign(from)
    return [
        (-s)^k * μ^(n - 2k) * inv(_order_ratio(T, n - 2k, n, m)) * T(2 * (n - 2k) + 1) / (2n + 1) *
            reexpansion_polynomials(n - 2k, k, x)[k + 1]
            for k in 0:((n - m) ÷ 2)
    ]
end
# Spheroidal interior harmonic in the solid harmonics `r^l P_l^m`:
#
#     aₖ = (-s)ᵏ (n+m)!/(n-m)! · (2n-2k)! / (2ⁿ k! (n-k)! (n-2k+m)!) · c^{-(n-2k)}.
function interior_reexpansion(::SphericalChart, from::_SpheroidalChart, n::Integer, m::Integer)
    _check_reexpansion(n, m)
    c, s = from.c, _spheroidal_sign(from)
    T = typeof(c / 2)
    ratio(k) = -s * c^2 * T((n - 2k + m) * (n - 2k + m - 1)) / (2 * (2n - 2k - 1) * (k + 1))
    return _from_ratios(_leading_coefficient(T, n, m) / c^n, ratio, (n - m) ÷ 2)
end
# Solid harmonic `r^n P_n^m` in the spheroidal interior harmonics.
function interior_reexpansion(to::_SpheroidalChart, ::SphericalChart, n::Integer, m::Integer)
    _check_reexpansion(n, m)
    c, s = to.c, _spheroidal_sign(to)
    T = typeof(c / 2)
    d = _leading_coefficient(T, n, m)
    return [
        s^k * c^n / d * inv(_order_ratio(T, n - 2k, n, m)) * T(2 * (n - 2k) + 1) / (2n - 2k + 1) *
            _pochhammer_ratio(T(2n - 2k + 1) / 2, k)
            for k in 0:((n - m) ÷ 2)
    ]
end

function _check_reexpansion(n::Integer, m::Integer)
    0 ≤ m ≤ n || throw(DomainError((n, m), "a harmonic of degree n and order m needs 0 ≤ m ≤ n"))
    return nothing
end

"""
    reexpansion_bound(to, from)

The value of the radial coordinate of `to` beyond which the exterior
re-expansion of [`exterior_reexpansion`](@ref) converges: the coordinate
surface of `to` through the edge of the focal set of `from`. Below it the
series diverges. The terms of the series decay as `e^{-2i(ξ-ξ*)}`, with
`ξ = arccosh q` or `arcsinh τ` and `ξ*` its value at the bound.

| `from` → `to` | bound |
| :------------ | :---- |
| prolate `c₁` → prolate `c₂` | `q₂ > max(1, c₁/c₂)` |
| oblate `c₁` → oblate `c₂` | `τ₂ > √(max(0, c₁²/c₂² - 1))` |
| prolate `c₁` → oblate `c₂` | `τ₂ > c₁/c₂` |
| oblate `c₁` → prolate `c₂` | `q₂ > √(1 + c₁²/c₂²)` |
| spheroidal `c₁` → spherical | `r > c₁` |
| spherical → spheroidal | the focal set of `to` |
"""
reexpansion_bound(to::ProlateSpheroidalChart, from::ProlateSpheroidalChart) =
    max(one(from.c / to.c), from.c / to.c)
reexpansion_bound(to::OblateSpheroidalChart, from::OblateSpheroidalChart) =
    sqrt(max(zero(from.c / to.c), (from.c / to.c)^2 - 1))
reexpansion_bound(to::OblateSpheroidalChart, from::ProlateSpheroidalChart) = from.c / to.c
reexpansion_bound(to::ProlateSpheroidalChart, from::OblateSpheroidalChart) = sqrt(1 + (from.c / to.c)^2)
reexpansion_bound(::SphericalChart, from::_SpheroidalChart) = from.c
reexpansion_bound(to::ProlateSpheroidalChart, ::SphericalChart) = one(to.c)
reexpansion_bound(to::OblateSpheroidalChart, ::SphericalChart) = zero(to.c)

export legendre_P, radial_P, radial_Q, interior_harmonics, exterior_harmonics
export reexpansion_polynomials, exterior_reexpansion, interior_reexpansion, reexpansion_bound
