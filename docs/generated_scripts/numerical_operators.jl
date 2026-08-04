import Pkg
Pkg.activate(joinpath(@__DIR__, "..", "docs"); io = devnull)
# The `docs` environment declares every dependency the tutorials use
# (SymPy, Symbolics, NLopt, BenchmarkTools, ...) and resolves TensND
# itself through `[sources]`. These four lines are stripped from the
# generated page and notebook, which already run inside that project.

using TensND
using LinearAlgebra
using Printf
using ForwardDiff

CS_polar = coorsys_polar_num()
r₀, θ₀ = 2.5, 0.7

println("  n  m       numerical        exact")
for (n, m) in [(2, 1), (3, 1), (3, 3), (2, 2)]
    f = x -> x[1]^n * sin(m * x[2])
    num = LAPLACE(f, CS_polar)([r₀, θ₀])
    exact = (n^2 - m^2) * r₀^(n - 2) * sin(m * θ₀)
    @printf "  %d  %d   %+.10f   %+.10f\n" n m num exact
end

CS_sph = coorsys_spherical_num()
x₀ = [π / 4, π / 3, 3.0]

println("   n     ∇[e_r]         exact          ∇²             exact")
for n in [1, 2, 3, -1]
    f = x -> x[3]^n
    g = GRAD(f, CS_sph)(x₀)
    lap = LAPLACE(f, CS_sph)(x₀)
    @printf "  %2d  %+.10f  %+.10f  %+.10f  %+.10f\n" n g[3] n * x₀[3]^(n - 1) lap n * (n + 1) * x₀[3]^(n - 2)
end

f_ang = x -> x[3]^2 * sin(x[1]) * cos(x[2])
g_num = GRAD(f_ang, CS_sph)(x₀)
g_exact = [
    x₀[3] * cos(x₀[1]) * cos(x₀[2]),      # e_θ
    -x₀[3] * sin(x₀[2]),                    # e_φ
    2x₀[3] * sin(x₀[1]) * cos(x₀[2]),     # e_r
]
(round.(g_num, digits = 10), round.(g_exact, digits = 10))

a, b = 1.0, 3.0
p_i = 1.0
κ, μ = 2.0, 1.0
λ = κ - 2μ / 3

A = p_i * a^3 / (3κ * (b^3 - a^3))
B = p_i * a^3 * b^3 / (4μ * (b^3 - a^3))

u_func = x -> [0.0, 0.0, A * x[3] + B / x[3]^2]          # (θ, φ, r) components
ε_func = x -> SYMGRAD(u_func, CS_sph)(x)
function σ_func(x)
    ε = ε_func(x)
    trε = sum(ε[i, i] for i in 1:3)
    return [λ * trε * (i == j) + 2μ * ε[i, j] for i in 1:3, j in 1:3]
end

θc, ϕc = π / 3, π / 4
println("   r      ε_rr        exact       ε_θθ        exact       σ_rr        exact      |div σ|")
for rv in [1.1, 1.5, 2.0, 2.5, 2.9]
    x = [θc, ϕc, rv]
    ε, σ = ε_func(x), σ_func(x)
    divσ = norm(DIV(σ_func, CS_sph)(x))
    @printf "  %.1f  %+.7f  %+.7f  %+.7f  %+.7f  %+.7f  %+.7f  %.1e\n" rv ε[3, 3] (A - 2B / rv^3) ε[1, 1] (A + B / rv^3) σ[3, 3] (3κ * A - 4μ * B / rv^3) divσ
end

(σ_func([θc, ϕc, a])[3, 3], σ_func([θc, ϕc, b])[3, 3])

∇lap = ForwardDiff.gradient(x -> LAPLACE(y -> y[3]^3, CS_sph)(x), x₀)
round.(∇lap, digits = 10)

∇ε = ForwardDiff.gradient(x -> SYMGRAD(y -> [0.0, 0.0, y[3]^2], CS_sph)(x)[3, 3], x₀)
round.(∇ε, digits = 10)

ForwardDiff.derivative(α -> LAPLACE(x -> α * x[3]^2, CS_sph)(x₀), 3.0)

ForwardDiff.derivative(α -> SYMGRAD(x -> [zero(α), zero(α), α * x[3]], CS_sph)(x₀)[3, 3], 2.0)

dαdiv = ForwardDiff.derivative(α -> DIV(x -> α * [0.0 0 0; 0 0 0; 0 0 1], CS_sph)(x₀)[3], 2.0)
(dαdiv, 2 / x₀[3])

function εθθ_of_p(pv)
    Av = pv * a^3 / (3κ * (b^3 - a^3))
    Bv = pv * a^3 * b^3 / (4μ * (b^3 - a^3))
    u = x -> [zero(pv), zero(pv), Av * x[3] + Bv / x[3]^2]
    return SYMGRAD(u, CS_sph)([θc, ϕc, 2.0])[1, 1]
end

dε_ad = ForwardDiff.derivative(εθθ_of_p, p_i)
@printf "dε_θθ/dp : AD = %+.12f   exact = %+.12f\n" dε_ad (εθθ_of_p(p_i) / p_i)

a_ell = 2.0
CS_ell = CoorSystemNum(x -> [a_ell * cosh(x[1]) * cos(x[2]), a_ell * sinh(x[1]) * sin(x[2])], 2)
x_ell = [1.0, 0.8]

(Lame(CS_ell, x_ell), a_ell * sqrt(sinh(x_ell[1])^2 + sin(x_ell[2])^2))

for (name, f) in ("μ" => (x -> x[1]), "cosh μ cos ν" => (x -> cosh(x[1]) * cos(x[2])))
    @printf "  ∇²(%-12s) = %+.3e\n" name LAPLACE(f, CS_ell)(x_ell)
end

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
