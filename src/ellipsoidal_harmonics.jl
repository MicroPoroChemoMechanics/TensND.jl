# ── Lamé functions and ellipsoidal harmonics ─────────────────────────────────
#
# In the ellipsoidal coordinates (ρ, μ, ν) of an `EllipsoidalChart` of
# semi-focal distances h₂ > h₃ (h₁² = h₂² - h₃², ν² ≤ h₃² ≤ μ² ≤ h₂² ≤ ρ²), a
# harmonic function separates as E(ρ) E(μ) E(ν), E being a Lamé function of
# degree n: one of the 2n + 1 solutions of
#
#     (x² - h₃²)(x² - h₂²) E″ + x (2x² - h₃² - h₂²) E′ + [(h₃² + h₂²) p - n(n+1) x²] E = 0
#
# of the form ψ(x) P(x), P a polynomial of every other power and ψ one of 1,
# √|x² - h₃²|, √|x² - h₂²| and their product: the classes K, L, M, N
# (dassios2012, chapter 3). Writing P = Σ_k c_k x^{d-2k} turns the equation into a
# three-term recurrence on the c_k whose truncation is a tridiagonal
# eigenproblem for p (dassios2012, Eqs. (3.22), (3.48), (3.81)); its eigenvalues
# give the functions of the class, with c₀ = 1.
#
# The interior harmonics E(ρ) E(μ) E(ν) are polynomials of degree n in x; the
# exterior ones are F(ρ) E(μ) E(ν), with the Lamé function of the second kind
#
#     F(ρ) = (2n + 1) E(ρ) I(ρ) ,   I(ρ) = ∫_ρ^∞ dt / [E(t)² √(t² - h₃²) √(t² - h₂²)] ,
#
# which decays as ρ^{-n-1}. The geometry enters the Lamé functions only through
# numbers compared with one another (the eigenvalues of the recurrence, sorted):
# they are computed in `Float64`, while the harmonics built on them may carry
# duals in the point and in ρ.

"""
    LameFunction

Lamé function ``E(x) = \\psi(x) \\sum_k c_k\\,x^{d - 2k}`` of degree `n` of an
[`EllipsoidalChart`](@ref), with ``\\psi`` the product of
``\\sqrt{|x^2 - h_3^2|}`` (if `a == 1`) and ``\\sqrt{|x^2 - h_2^2|}`` (if `b == 1`),
``d = n - a - b``, `p` its separation constant and `coefficients` the ``c_k``,
``c_0 = 1``.

See also [`lame_functions`](@ref), [`lame_value`](@ref), [`lame_second_kind`](@ref).
"""
struct LameFunction
    n::Int
    a::Int
    b::Int
    p::Float64
    coefficients::Vector{Float64}
end

"""
    lame_functions(chart::EllipsoidalChart, n) -> Vector{LameFunction}

The ``2n + 1`` Lamé functions of degree `n` of `chart`, by class (K, L, M, N)
and, within a class, by increasing separation constant [dassios2012](@cite). They
are computed in `Float64` from the primal values of the semi-focal distances:
they come from eigenvalues sorted by value.

# Examples
```jldoctest
julia> length(lame_functions(EllipsoidalChart(1.2, 0.8), 3))
7
```
"""
function lame_functions(chart::EllipsoidalChart, n::Integer)
    n ≥ 0 || throw(DomainError(n, "the degree of a Lamé function must be non-negative"))
    h₂², h₃² = Float64(primal(chart.h₂))^2, Float64(primal(chart.h₃))^2
    α, β = h₃² + h₂², h₃² * h₂²
    out = LameFunction[]
    for (a, b) in ((0, 0), (1, 0), (0, 1), (1, 1))
        d = n - a - b
        d ≥ 0 || continue
        r = d ÷ 2 + 1                   # coefficients of x^d, x^{d-2}, …, x^{d mod 2}
        # 2(k+1)(2n-2k-1) c_{k+1} = [α (p - s_k²) - e_k] c_k + β f_k c_{k-1}, that
        # is α p c = A c with A tridiagonal, c_r = 0 closing it
        A = zeros(r, r)
        for k in 0:(r - 1)
            s, e, f = _lame_recurrence(n, a, b, k, h₂², h₃²)
            A[k + 1, k + 1] = α * s^2 + e
            k + 1 < r && (A[k + 1, k + 2] = 2 * (k + 1) * (2n - 2k - 1))
            k > 0 && (A[k + 1, k] = -β * f)
        end
        F = eigen(A ./ α)
        for j in sortperm(real.(F.values))
            c = real.(F.vectors[:, j])
            push!(out, LameFunction(n, a, b, real(F.values[j]), c ./ c[1]))
        end
    end
    return out
end

# The terms of the recurrence of the class (a, b) at the step k: the square s²
# of the shifted degree, the diagonal term e of the classes L and M, and the
# factor f of c_{k-1} (dassios2012, Eqs. (3.22), (3.48), (3.81)).
function _lame_recurrence(n, a, b, k, h₂², h₃²)
    a == b == 0 && return (n - 2k, 0.0, (n - 2k + 2) * (n - 2k + 1))
    if a + b == 1
        h² = a == 1 ? h₂² : h₃²          # class L carries h₂² in its diagonal, M h₃²
        return (n - 2k - 1, (2n - 4k - 1) * h², (n - 2k + 1) * (n - 2k))
    end
    return (n - 2k - 1, 0.0, (n - 2k) * (n - 2k - 1))
end

# The polynomial part of E and its derivative at x.
function _lame_polynomial(E::LameFunction, x)
    d = E.n - E.a - E.b
    P, dP = zero(x), zero(x)
    for (k, c) in enumerate(E.coefficients)
        j = d - 2 * (k - 1)
        P += c * x^j
        j > 0 && (dP += c * j * x^(j - 1))
    end
    return P, dP
end

"""
    lame_value(E, x, s₃, s₂, σ₃, σ₂) -> (value, derivative)
    lame_value(E, c, i) -> (value, derivative)

The Lamé function `E` and its derivative at `x`, the square roots of its class
being given as ``s_3 = \\sqrt{|x^2 - h_3^2|}`` and ``s_2 = \\sqrt{|x^2 - h_2^2|}``
with the sign of the coordinate's branch, their derivatives following from
``s\\,s' = \\sigma x``, ``\\sigma`` (`σ₃`, `σ₂`) the sign of ``x^2 - h^2``. The second
form takes the coordinates `c` of [`ellipsoidal_signed_coordinates`](@ref) and
the coordinate `i` = 1, 2, 3 (``\\rho``, ``\\mu``, ``\\nu``).
"""
function lame_value(E::LameFunction, x, s₃, s₂, σ₃::Int, σ₂::Int)
    P, dP = _lame_polynomial(E, x)
    ψ, dψ = one(x), zero(x)
    if E.a == 1
        ψ, dψ = ψ * s₃, dψ * s₃ + ψ * σ₃ * x / s₃
    end
    if E.b == 1
        ψ, dψ = ψ * s₂, dψ * s₂ + ψ * σ₂ * x / s₂
    end
    return ψ * P, dψ * P + ψ * dP
end
function lame_value(E::LameFunction, c::NamedTuple, i::Integer)
    i == 1 && return lame_value(E, c.ρ, c.ρ3, c.ρ2, 1, 1)
    i == 2 && return lame_value(E, c.μ, c.μ3, c.μ2, 1, -1)
    i == 3 || throw(DomainError(i, "the ellipsoidal coordinate is 1 (ρ), 2 (μ) or 3 (ν)"))
    return lame_value(E, c.ν, c.ν3, c.ν2, -1, -1)
end

"""
    ellipsoidal_signed_coordinates(chart::EllipsoidalChart, x) -> NamedTuple

The ellipsoidal coordinates ``(\\rho, \\mu, \\nu)`` of the point `x`, from
[`cartesian_to_coords`](@ref), and the six square roots
``\\sqrt{|\\xi^2 - h_k^2|}`` that enter the Lamé functions, signed by the octant
of `x` [dassios2012; § 1.2](@cite): ``\\nu`` carries the sign of ``x_1``,
``\\sqrt{h_3^2 - \\nu^2}`` that of ``x_2`` and ``\\sqrt{h_2^2 - \\mu^2}`` that of
``x_3``. Fields `ρ, μ, ν, ρ3, ρ2, μ3, μ2, ν3, ν2`. With them the products
``E(\\rho)E(\\mu)E(\\nu)`` are harmonic in the whole space, whereas the chart alone
covers the quarter space ``x_2, x_3 \\ge 0``.
"""
function ellipsoidal_signed_coordinates(chart::EllipsoidalChart, x)
    ρ, μ, ν = cartesian_to_coords(chart, x)
    h₂², h₃² = chart.h₂^2, chart.h₃^2
    return (;
        ρ, μ, ν,
        ρ3 = sqrt(ρ^2 - h₃²), ρ2 = sqrt(ρ^2 - h₂²),
        μ3 = sqrt(μ^2 - h₃²), μ2 = _octant_sign(x[3]) * sqrt(h₂² - μ^2),
        ν3 = _octant_sign(x[2]) * sqrt(h₃² - ν^2), ν2 = sqrt(h₂² - ν^2),
    )
end
# The sign of a Cartesian component for the octants, +1 on its zero plane.
_octant_sign(v) = primal(v) < 0 ? -1 : 1

"""
    lame_second_kind(chart::EllipsoidalChart, E, ρ[, nodes]) -> (F, dF)

The Lamé function of the second kind ``F(\\rho) = (2n+1)\\,E(\\rho)\\,I(\\rho)`` and
its derivative, with
``I(\\rho) = \\int_\\rho^\\infty \\mathrm{d}t / [E(t)^2 \\sqrt{t^2-h_3^2}\\sqrt{t^2-h_2^2}]``,
computed as
``\\rho \\int_0^1 \\mathrm{d}s / [E(\\rho/s)^2 \\sqrt{\\rho^2-h_3^2 s^2}\\sqrt{\\rho^2-h_2^2 s^2}]``.
The factor ``\\sqrt{\\rho^2-h_2^2 s^2}`` falls to ``\\kappa = \\sqrt{\\rho^2-h_2^2}`` at
``s = 1``: for an elongated ellipsoid the integrand varies there over a length
``(\\kappa/h_2)^2``. The interval is split at ``s = 1/2``; on ``[1/2, 1]`` the
variable ``\\tau``, with ``1 - s^2 = (\\kappa/h_2)^2 \\sinh^2\\tau``, turns that
factor into ``\\kappa\\cosh\\tau`` and leaves a smooth integrand. Both parts are
integrated by Gauss–Legendre on `nodes` (nodes and weights on ``[0, 1]``,
``n + 24`` of them by default; a caller evaluating many functions passes them
once). `ρ` may carry a derivative.
"""
function lame_second_kind(chart::EllipsoidalChart, E::LameFunction, ρ, nodes = _gauss_legendre_01(E.n + 24))
    h₂, h₃ = chart.h₂, chart.h₃
    xs, ws = nodes
    κ = sqrt(ρ^2 - h₂^2)
    τ₁ = asinh(sqrt(oftype(κ, 3)) * h₂ / (2κ))     # s = 1/2
    Eρs(s) = first(lame_value(E, ρ / s, sqrt((ρ / s)^2 - h₃^2), sqrt((ρ / s)^2 - h₂^2), 1, 1))
    I = zero(κ)
    for (x, w) in zip(xs, ws)
        s = x / 2                                             # on [0, 1/2]
        I += (w / 2) / (Eρs(s)^2 * sqrt(ρ^2 - h₃^2 * s^2) * sqrt(ρ^2 - h₂^2 * s^2))
        τ = τ₁ * x                                            # on [1/2, 1]
        sh = sinh(τ)
        s = sqrt(1 - (κ * sh / h₂)^2)
        I += τ₁ * w * κ * sh / (h₂^2 * s * Eρs(s)^2 * sqrt(ρ^2 - h₃^2 * s^2))
    end
    I *= ρ
    s₃ = sqrt(ρ^2 - h₃^2)
    Eρ, dEρ = lame_value(E, ρ, s₃, κ, 1, 1)
    dI = -1 / (Eρ^2 * s₃ * κ)
    return (2E.n + 1) * Eρ * I, (2E.n + 1) * (dEρ * I + Eρ * dI)
end

# Gauss–Legendre rule of n nodes on [0, 1], from the Jacobi matrix
# (Golub–Welsch): its eigenvalues are the nodes on [-1, 1], the squares of the
# first components of its eigenvectors the weights over 2.
function _gauss_legendre_01(n::Integer)
    β = [k / sqrt(4k^2 - 1) for k in 1:(n - 1)]
    F = eigen(SymTridiagonal(zeros(n), β))
    return (1 .+ F.values) ./ 2, F.vectors[1, :] .^ 2
end

export LameFunction, lame_functions, lame_value, lame_second_kind, ellipsoidal_signed_coordinates
