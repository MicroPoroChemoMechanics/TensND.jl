# ============================================================================
#  Numerical differential operators with CoorSystemNum
#  (automatic differentiation via ForwardDiff)
#
#  This script illustrates:
#   1. Basic scalar/vector operators in polar and spherical coordinates
#   2. Lamé problem: hollow sphere under internal pressure
#      - numerical strain via SYMGRAD
#      - numerical stress divergence via DIV (rank-2)
#      - parametric study over radius
#   3. Custom coordinate system from an OM(x) function (elliptic coords)
# ============================================================================

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)

using TensND, LinearAlgebra, Printf

# ─────────────────────────────────────────────────────────────────────────────
# 1.  Polar 2D — Laplacian of r^n sin(mθ)
# ─────────────────────────────────────────────────────────────────────────────
#
# Exact result: ∇²(r^n f(θ)) = [n(n-1) + n - m²] r^(n-2) f(θ)  = (n²-m²) r^(n-2) f(θ)
# For n=2, m=1:  ∇²(r² sinθ) = (4-1) r⁰ sinθ = 3 sinθ
# For n=3, m=1:  ∇²(r³ sinθ) = (9-1) r sinθ = 8r sinθ

CS_polar = coorsys_polar_num()

println("── Polar coordinate system ──────────────────────────────")
for (n, m) in [(2,1), (3,1), (3,3)]
    f = x -> x[1]^n * sin(m * x[2])   # f(r,θ) = rⁿ sin(mθ)
    r₀, θ₀ = 2.5, 0.7
    lap = LAPLACE(f, CS_polar)([r₀, θ₀])
    exact = (n^2 - m^2) * r₀^(n-2) * sin(m * θ₀)
    @printf("  ∇²(r^%d sin(%dθ)) at (r=%.1f,θ=%.2f):  num = %+.6f   exact = %+.6f\n",
            n, m, r₀, θ₀, lap, exact)
end

# ─────────────────────────────────────────────────────────────────────────────
# 2.  Spherical 3D — Gradient and Laplacian of classical functions
# ─────────────────────────────────────────────────────────────────────────────
#
# Coordinates: (θ, ϕ, r)  — matches CoorSystemSym convention
# Lamé: χ = (r, r sinθ, 1)
#
# Known results (physical components in the normalized basis):
#   ∇(r^n) = n r^(n-1) 𝐞ʳ
#   ∇²(r^n) = n(n+1) r^(n-2)   (for n ≠ 0)

CS_sph = coorsys_spherical_num()

println("\n── Spherical coordinate system ──────────────────────────")
θ₀, ϕ₀, r₀ = π/4, π/3, 3.0
x₀ = [θ₀, ϕ₀, r₀]

for n in [1, 2, 3, -1]
    f = x -> x[3]^n   # f(θ,ϕ,r) = rⁿ
    g = GRAD(f, CS_sph)(x₀)
    lap = LAPLACE(f, CS_sph)(x₀)
    # physical gradient: only e_r component ≠ 0
    grad_r_exact = n * r₀^(n-1)
    lap_exact    = n * (n + 1) * r₀^(n-2)
    @printf("  r^%-2d :  ∇r[eᵣ] = %+.6f (exact %+.6f),  ∇² = %+.6f (exact %+.6f)\n",
            n, g[3], grad_r_exact, lap, lap_exact)
end

# Gradient of a scalar with angular dependence: f = r² sinθ cosϕ
f_ang = x -> x[3]^2 * sin(x[1]) * cos(x[2])
g_ang = GRAD(f_ang, CS_sph)(x₀)
# Exact: ∇f = 2r sinθ cosϕ 𝐞ʳ + r cosθ cosϕ 𝐞ᶿ - r sinϕ 𝐞ᵠ  (physical)
g_exact = [r₀*cos(θ₀)*cos(ϕ₀),   # eθ component
           -r₀*sin(ϕ₀),           # eϕ component
           2r₀*sin(θ₀)*cos(ϕ₀)]   # er component
println("\n  ∇(r² sinθ cosϕ) at (θ=π/4, ϕ=π/3, r=3):")
println("    num   = $([round(g_ang[i],sigdigits=7) for i in 1:3])")
println("    exact = $([round(g_exact[i],sigdigits=7) for i in 1:3])")

# ─────────────────────────────────────────────────────────────────────────────
# 3.  Lamé problem — hollow elastic sphere under internal pressure
# ─────────────────────────────────────────────────────────────────────────────
#
# Geometry:  inner radius a, outer radius b
# Loading:   internal pressure p (outer surface free)
# Material:  isotropic, bulk modulus κ, shear modulus μ
#
# Exact solution (Lamé):
#   u_r(r) = A·r + B/r²
#   A = p a³ / [3κ(b³-a³)]
#   B = p a³ b³ / [4μ(b³-a³)]
#
# Non-zero physical strain components:
#   εᵣᵣ = A - 2B/r³
#   εθθ = εϕϕ = A + B/r³
#
# Stress (isotropic Hooke):   σ = λ tr(ε) 1 + 2μ ε,  λ = κ - 2μ/3
#   σᵣᵣ(a) = -p,  σᵣᵣ(b) = 0,  div σ = 0

println("\n── Lamé problem: hollow sphere under internal pressure ──")

a, b   = 1.0, 3.0     # inner / outer radius
p_i    = 1.0           # internal pressure
κ_val  = 2.0           # bulk modulus
μ_val  = 1.0           # shear modulus
λ_val  = κ_val - 2μ_val/3

A = p_i * a^3 / (3κ_val * (b^3 - a^3))
B = p_i * a^3 * b^3 / (4μ_val * (b^3 - a^3))

# Radial displacement field: u = (0, 0, u_r(r)) in the (θ,ϕ,r) normalized basis
u_func = x -> [0.0, 0.0, A*x[3] + B/x[3]^2]

# Strain and stress functions
ε_func = x -> SYMGRAD(u_func, CS_sph)(x)
σ_func = x -> begin
    ε = ε_func(x)
    tr_ε = sum(ε[i,i] for i in 1:3)
    [λ_val*tr_ε*(i==j ? 1.0 : 0.0) + 2μ_val*ε[i,j] for i in 1:3, j in 1:3]
end

θ₀, ϕ₀ = π/3, π/4

println("\n  r     ε_rr(num)   ε_rr(exact)  ε_tt(num)  ε_tt(exact)  σ_rr(num)  σ_rr(exact)  |div σ|")
println("  ", "-"^90)
for r in [1.1, 1.5, 2.0, 2.5, 2.9]
    x = [θ₀, ϕ₀, r]
    ε = ε_func(x)
    σ = σ_func(x)
    divσ_norm = norm(DIV(σ_func, CS_sph)(x))

    ε_rr_ex = A - 2B/r^3
    ε_tt_ex = A + B/r^3
    σ_rr_ex = (3κ_val*A - 4μ_val*B/r^3)

    @printf("  %.1f  %+.6f  %+.6f   %+.6f  %+.6f   %+.6f  %+.6f   %.2e\n",
            r, ε[3,3], ε_rr_ex, ε[1,1], ε_tt_ex, σ[3,3], σ_rr_ex, divσ_norm)
end

# Verify boundary conditions
println("\n  Boundary conditions:")
σ_inner = σ_func([θ₀, ϕ₀, a])
σ_outer = σ_func([θ₀, ϕ₀, b])
@printf("    σ_rr(r=a) = %.6f  (should be -%.6f)\n", σ_inner[3,3], p_i)
@printf("    σ_rr(r=b) = %.6f  (should be  0.000000)\n", σ_outer[3,3])

# ─────────────────────────────────────────────────────────────────────────────
# 4.  Generic constructor from OM function — toroidal-like coordinate
#     (prolate spheroidal in 2D as demonstration)
# ─────────────────────────────────────────────────────────────────────────────
#
# Use elliptic coords: OM(μ,ν) = (a cosh(μ)cos(ν), a sinh(μ)sin(ν))
# χ₁ = χ₂ = a√(sinh²μ + sin²ν)
# ∇²f = [∂²f/∂μ² + ∂²f/∂ν²] / (a²(sinh²μ + sin²ν))

println("\n── Elliptic coordinates (generic OM constructor) ────────")
a_ell = 2.0
OM_elliptic = x -> [a_ell * cosh(x[1]) * cos(x[2]),
                    a_ell * sinh(x[1]) * sin(x[2])]
CS_ell = CoorSystemNum(OM_elliptic, 2)

μ₀, ν₀ = 1.0, 0.8
x₀ = [μ₀, ν₀]
h² = a_ell^2 * (sinh(μ₀)^2 + sin(ν₀)^2)

# Test Laplacian of f = μ (linear in μ → ∂²/∂μ²=0, ∂²/∂ν²=0) → ∇²μ = 0
f_mu = x -> x[1]
lap_mu = LAPLACE(f_mu, CS_ell)(x₀)
@printf("  ∇²μ = %.2e  (should be 0)\n", lap_mu)

# Test: ∇²(cosh μ cos ν) = 0  (harmonic function)
f_harm = x -> cosh(x[1]) * cos(x[2])
lap_harm = LAPLACE(f_harm, CS_ell)(x₀)
@printf("  ∇²(cosh μ cos ν) = %.2e  (should be 0)\n", lap_harm)
