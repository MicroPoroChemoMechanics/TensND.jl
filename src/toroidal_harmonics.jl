# ── Toroidal functions ───────────────────────────────────────────────────────
#
# The harmonic functions that separate in a toroidal chart are, up to the
# factor `√(cosh α - cos β)`, products of `e^{inβ}`, `e^{imϕ}` and a Legendre
# function of half-integer degree `n - 1/2`, order `m` and argument
# `cosh α ≥ 1`. This file evaluates those functions of `α` and their
# derivatives.
#
# Conventions. With `x = cosh α`, `w = sinh α = √(x²-1)` and `ν = n - 1/2`,
#
#     P_ν^m(x) = wᵐ dᵐP_ν/dxᵐ,      Q_ν^m(x) = wᵐ dᵐQ_ν/dxᵐ        (x > 1),
#
# `P_ν` being the Legendre function of the first kind, equal to 1 at `x = 1`,
# and `Q_ν` that of the second kind, which vanishes at infinity. Both kinds
# obey, in the degree and in the order,
#
#     (ν-m+1) F_{ν+1}^m = (2ν+1) x F_ν^m - (ν+m) F_{ν-1}^m,
#     w F_ν^{m+1} = (ν-m) x F_ν^m - (ν+m) F_{ν-1}^m,
#     w dF_ν^m/dx = [ν x F_ν^m - (ν+m) F_{ν-1}^m] / w,
#
# and are even in `n`: `F_{-n-1/2}^m = F_{n-1/2}^m`, which gives the degree
# below `n = 0`.

# The complete elliptic integral `K(k)` and the difference `K(k) - E(k)`, from
# the modulus `k` and its complement `k′ = √(1-k²)` given separately, by the
# arithmetic-geometric mean: `K = π/(2 a_∞)` and `K - E = K Σ 2^{n-1} c_n²`, a
# sum of positive terms. The iteration stops on a value (numbers and duals,
# not symbols); one more step after the test brings the derivatives carried by
# duals to the same accuracy, the convergence being quadratic.
function _complete_elliptic(k, k′)
    a, b, c = one(k), k′, k
    weight = one(k) / 2
    S = weight * c^2
    tol = eps(typeof(primal(a)))
    converged = false
    for _ in 1:64
        a, b, c = (a + b) / 2, sqrt(a * b), (a - b) / 2
        weight *= 2
        S += weight * c^2
        converged && return (π / (2a), π * S / (2a))
        converged = abs(primal(c)) ≤ tol * abs(primal(a))
    end
    # Reached only for a complementary modulus that is zero or not a number.
    throw(DomainError(k′, "the arithmetic-geometric mean did not converge for this complementary modulus"))
end

# The toroidal functions are computed by iterations sized on the value of `α`:
# numbers and duals, not symbols.
function _check_toroidal_argument(α, lowest)
    is_hard_numeric(typeof(α)) || throw(
        ArgumentError(
            "the toroidal functions are evaluated by iterations that stop on the value of α: " *
                "a number or a dual is needed, got a $(typeof(α))",
        )
    )
    primal(α) ≥ lowest || throw(DomainError(α, "the toroidal coordinate α is non-negative"))
    return nothing
end

# `P_{-1/2}(cosh α) = (2/π) sech(α/2) K(tanh(α/2))` and
# `P_{1/2}(cosh α) = (2/π) e^{α/2} E(√(1-e^{-2α}))`.
function _toroidal_P_seeds(α)
    k′ = sech(α / 2)
    K, _ = _complete_elliptic(tanh(α / 2), k′)
    q = exp(-α)
    K₁, KE₁ = _complete_elliptic(sqrt(-expm1(-2α)), q)
    return (2 * k′ * K / π, 2 * (K₁ - KE₁) / (π * sqrt(q)))
end

# `Q_{-1/2}(cosh α) = 2 e^{-α/2} K(e^{-α})` and
# `Q_{1/2}(cosh α) = 2 e^{α/2} [K(e^{-α}) - E(e^{-α})]`: the difference comes
# from the mean as a sum, with no cancellation at large `α`.
function _toroidal_Q_seeds(α)
    q = exp(-α)
    K, KE = _complete_elliptic(q, sqrt(-expm1(-2α)))
    return (2 * sqrt(q) * K, 2 * KE / sqrt(q))
end

# The degree below `n - 1/2`, in a vector of the degrees `0 ≤ n`: that of
# `n = 1` under `n = 0`, the functions being even in `n`.
_toroidal_below(F::AbstractVector, n::Integer) = F[iszero(n) ? 2 : n]

# From the order `k` to the order `k+1`, at every degree:
# `w F_ν^{k+1} = (ν-k) x F_ν^k - (ν+k) F_{ν-1}^k`.
_toroidal_raise_order(F::AbstractVector, k::Integer, x, w) =
    [((n - 1 // 2 - k) * x * F[n + 1] - (n - 1 // 2 + k) * _toroidal_below(F, n)) / w for n in 0:(length(F) - 1)]

# `dF_ν^m/dα = [ν x F_ν^m - (ν+m) F_{ν-1}^m] / w`, both kinds.
_toroidal_derivative(F::AbstractVector, m::Integer, x, w) =
    [((n - 1 // 2) * x * F[n + 1] - (n - 1 // 2 + m) * _toroidal_below(F, n)) / w for n in 0:(length(F) - 1)]

# The functions of the first kind at the degrees `n - 1/2`, `0 ≤ n ≤ N` with
# `N ≥ 1`: the dominant solution of the recurrence in the degree, run upward.
function _toroidal_P(N::Integer, m::Integer, α)
    x, w = cosh(α), sinh(α)
    F = _legendre_forward(N, 0, x, _toroidal_P_seeds(α)..., -1, -1 // 2)
    return foldl((G, k) -> _toroidal_raise_order(G, k, x, w), 0:(m - 1); init = F)
end

# The functions of the second kind: the minimal solution. Run upward, the
# recurrence loses it to rounding by the factor `e^{2αn}` at the degree
# `n - 1/2`. Where that factor stays below 16 over the degrees asked for, the
# recurrence is run upward from the two seeds; elsewhere the ratios
# `Q_{n-1/2}/Q_{n-3/2}` are computed downward from their limit `e^{-α}` far
# above `N`, and the first seed gives the values.
function _toroidal_Q(N::Integer, m::Integer, α)
    x, w = cosh(α), sinh(α)
    Q₀, Q₁ = _toroidal_Q_seeds(α)
    F = 2 * primal(α) * N ≤ 4 * log(oftype(primal(α), 2)) ?
        _legendre_forward(N, 0, x, Q₀, Q₁, -1, -1 // 2) : _toroidal_Q_downward(N, x, α, Q₀)
    return foldl((G, k) -> _toroidal_raise_order(G, k, x, w), 0:(m - 1); init = F)
end

function _toroidal_Q_downward(N::Integer, x, α, Q₀)
    ratio = exp(-α)
    ξ = primal(α)
    # A start `L` degrees above `N` leaves an error `e^{-2Lα}` on the ratio at
    # `N`: `L` is taken so that it falls below the rounding of the number type.
    L = ceil(Int, precision(typeof(ξ)) * log(2) / (2 * ξ)) + 20
    v = zeros(typeof(ratio * Q₀), N + 1)     # holds Q_{n-1/2}/Q_{n-3/2} first, for 1 ≤ n ≤ N
    for n in (N + L):-1:1
        ν = n - 1 // 2
        ratio = ν / (2n * x - (ν + 1) * ratio)
        n ≤ N && (v[n + 1] = ratio)
    end
    v[1] = Q₀
    return cumprod(v)
end

"""
    toroidal_P(N, m, α) -> Vector

Toroidal functions of the first kind: `v[n+1] = P_{n-1/2}^m(cosh α)` for
`0 ≤ n ≤ N`, with `P_ν^m(x) = (√(x²-1))^m dᵐP_ν/dxᵐ` and `m ≥ 0`. They are
bounded on the axis of the toroidal chart, `α = 0`, where `P_{n-1/2} = 1`, and
grow toward its focal circle as `e^{(n-1/2)α}`.

`α` is a number or a dual, `α ≥ 0`; at an order `m ≥ 1` it must be positive,
the functions being then computed as a quotient by `sinh α`.

# Examples
```jldoctest
julia> toroidal_P(2, 0, 0.0) ≈ [1, 1, 1]
true
```

See also [`toroidal_Q`](@ref), [`toroidal_dP`](@ref), [`ToroidalChart`](@ref).
"""
function toroidal_P(N::Integer, m::Integer, α)
    _check_degree_order(N, m)
    _check_toroidal_argument(α, zero(primal(α)))
    _check_off_axis(m, α)
    return resize!(_toroidal_P(max(N, 1), m, α), N + 1)
end

_check_off_axis(m::Integer, α) = iszero(m) || primal(α) > 0 || throw(
    DomainError(α, "at an order m ≥ 1 the toroidal functions are evaluated off the axis, α > 0")
)

"""
    toroidal_Q(N, m, α) -> Vector

Toroidal functions of the second kind: `v[n+1] = Q_{n-1/2}^m(cosh α)` for
`0 ≤ n ≤ N`, with `Q_ν^m(x) = (√(x²-1))^m dᵐQ_ν/dxᵐ` and `m ≥ 0`. They decay
toward the focal circle of the toroidal chart as `e^{-(n+1/2)α}` and are
singular on its axis, `α = 0`.

`Q_{n-1/2}` is the minimal solution of the recurrence in the degree: it is
computed upward where the loss `e^{2αN}` of that run stays below 16, and
elsewhere from the ratios of consecutive degrees, computed downward. The result
is accurate to rounding at every degree and differentiable (duals in `α`).

# Examples
```jldoctest
julia> toroidal_Q(0, 0, 20.0)[1] ≈ π * exp(-10.0)
true
```

See also [`toroidal_P`](@ref), [`toroidal_dQ`](@ref), [`ToroidalChart`](@ref).
"""
function toroidal_Q(N::Integer, m::Integer, α)
    _check_degree_order(N, m)
    _check_toroidal_argument(α, zero(primal(α)))
    _check_off_singular_axis(α)
    return resize!(_toroidal_Q(max(N, 1), m, α), N + 1)
end

_check_off_singular_axis(α) = primal(α) > 0 || throw(
    DomainError(α, "the toroidal functions of the second kind are singular on the axis of the chart, α = 0")
)

"""
    toroidal_dP(N, m, α) -> Vector

Derivatives with respect to `α` of the functions of [`toroidal_P`](@ref),
`v[n+1] = d P_{n-1/2}^m(cosh α)/dα` for `0 ≤ n ≤ N`, from
`sinh α dF_ν^m/dα = ν cosh α F_ν^m - (ν+m) F_{ν-1}^m`. `α` is positive.
"""
function toroidal_dP(N::Integer, m::Integer, α)
    _check_degree_order(N, m)
    _check_toroidal_argument(α, zero(primal(α)))
    _check_off_axis(1, α)
    return resize!(_toroidal_derivative(_toroidal_P(max(N, 1), m, α), m, cosh(α), sinh(α)), N + 1)
end

"""
    toroidal_dQ(N, m, α) -> Vector

Derivatives with respect to `α` of the functions of [`toroidal_Q`](@ref),
`v[n+1] = d Q_{n-1/2}^m(cosh α)/dα` for `0 ≤ n ≤ N`, from the same relation as
[`toroidal_dP`](@ref).
"""
function toroidal_dQ(N::Integer, m::Integer, α)
    _check_degree_order(N, m)
    _check_toroidal_argument(α, zero(primal(α)))
    _check_off_singular_axis(α)
    return resize!(_toroidal_derivative(_toroidal_Q(max(N, 1), m, α), m, cosh(α), sinh(α)), N + 1)
end

export toroidal_P, toroidal_Q, toroidal_dP, toroidal_dQ
