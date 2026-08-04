import Pkg
Pkg.activate(joinpath(@__DIR__, "..", "docs"); io = devnull)
# The `docs` environment declares every dependency the tutorials use
# (SymPy, Symbolics, NLopt, BenchmarkTools, ...) and resolves TensND
# itself through `[sources]`. These four lines are stripped from the
# generated page and notebook, which already run inside that project.

using TensND
using LinearAlgebra
using SymPy
using Printf

Spherical = coorsys_spherical()
coords = getcoords(Spherical)
Γ = Christoffel(Spherical)
d = 3

residuals = [
    tsimplify(
        ∂(natvec(Spherical, j, :cov), coords[i], Spherical) ⋅ natvec(Spherical, k, :cont)
            - Γ[i, j, k],
    )
        for i in 1:d, j in 1:d, k in 1:d
]
println("all Γ recovered from the definition ? ", all(iszero, residuals))

tsimplify(maximum(abs.(Γ - permutedims(Γ, (2, 1, 3)))))

function christoffel_table(name, CS)
    c = getcoords(CS)
    G = Christoffel(CS)
    dim = length(c)
    println("── ", name, "   coords = ", c, "   Lamé = ", Lame(CS))
    any_nz = false
    for i in 1:dim, j in 1:dim, k in 1:dim
        v = tsimplify(G[i, j, k])
        if !iszero(v)
            @printf "     Γ^%s_%s%s = %s\n" string(c[k]) string(c[i]) string(c[j]) string(v)
            any_nz = true
        end
    end
    any_nz || println("     (all vanish — this is a Cartesian chart)")
    return println()
end

christoffel_table("cartesian", coorsys_cartesian())
christoffel_table("polar", coorsys_polar())
christoffel_table("cylindrical", coorsys_cylindrical())
christoffel_table("spherical", coorsys_spherical())

θ, ϕ, r = getcoords(Spherical)
𝐞ᶿ, 𝐞ᵠ, 𝐞ʳ = unitvec(Spherical)
@set_coorsys Spherical

get_array(𝐞ʳ), tsimplify(DIV(𝐞ʳ))

CSnum = coorsys_spherical_num()
x₀ = [0.7, 1.1, 2.3]                       # (θ, ϕ, r)

Γnum = Christoffel(CSnum, x₀)
Γsym = [Float64(tsimplify(Γ[i, j, k]).subs(Dict(θ => x₀[1], ϕ => x₀[2], r => x₀[3]))) for i in 1:3, j in 1:3, k in 1:3]

println("‖Γ_num − Γ_sym‖ at (θ,ϕ,r) = ", x₀, " : ", norm(Γnum - Γsym))

(Lame(CSnum, x₀), (x₀[3], x₀[3] * sin(x₀[1]), 1.0))

Spheroidal = coorsys_spheroidal()
getcoords(Spheroidal)

components_canon(getOM(Spheroidal))

Lame(Spheroidal)

simplify(LAPLACE(getOM(Spheroidal)[1]^2, Spheroidal))

ϕs, ps, qs = getcoords(Spheroidal)
P = sympy.assoc_legendre

for (n, m) in ((2, 0), (2, 1), (3, 2))
    T = P(n, m, ps) * P(n, m, qs) * cos(m * ϕs)
    lap = simplify(LAPLACE(T, Spheroidal))
    @printf "  ∇²[P_%d^%d(p) P_%d^%d(q) cos(%dφ)] = %s\n" n m n m m string(lap)
end

a_ell = 2.0
OM_elliptic = x -> [a_ell * cosh(x[1]) * cos(x[2]), a_ell * sinh(x[1]) * sin(x[2])]
CS_ell = CoorSystemNum(OM_elliptic, 2)

x_ell = [1.0, 0.8]
χ_ell = Lame(CS_ell, x_ell)
χ_exact = a_ell * sqrt(sinh(x_ell[1])^2 + sin(x_ell[2])^2)
println("Lamé  : ", round.(χ_ell, digits = 10), "   exact (both equal) : ", round(χ_exact, digits = 10))

lap_harm = LAPLACE(x -> cosh(x[1]) * cos(x[2]), CS_ell)(x_ell)
@printf "∇²(cosh μ cos ν) = %.3e   (should vanish)\n" lap_harm

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
