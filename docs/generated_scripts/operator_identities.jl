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

u, v, w = symbols("u v w", real = true)

charts = [
    "cartesian" => coorsys_cartesian(),
    "polar" => coorsys_polar(),
    "spherical" => coorsys_spherical(),
    "non-orthogonal" => CoorSystemSym(Tens([u + v^2, v]), (u, v)),
]

for (name, CS) in charts
    ortho = isorthonormal(normalized_basis(CS))
    @printf "  %-16s dim %d   coords %-14s orthogonal: %s\n" name get_dim(CS) string(getcoords(CS)) ortho
end

zeroresidual(x::Number) = tsimplify(x)
zeroresidual(x::AbstractArray) = tsimplify(maximum(abs, tsimplify.(x)))
zeroresidual(t::AbstractTens) = zeroresidual(get_array(t))

function check(label, chart, residual)
    r = zeroresidual(residual)
    ok = iszero(r)
    @printf "  %-14s %-42s %s\n" chart label (ok ? "0  ✓" : "→ $r")
    return ok
end

println("\n── scalar identities ──")
for (name, CS) in charts
    q = getcoords(CS)
    f = SymFunction("f", real = true)(q...)
    g = SymFunction("g", real = true)(q...)

    check("∇²f = DIV(GRAD f)", name, LAPLACE(f, CS) - DIV(GRAD(f, CS), CS))
    check("∇²f = tr(HESS f)", name, LAPLACE(f, CS) - tr(HESS(f, CS)))
    check("HESS f symmetric", name, get_array(HESS(f, CS)) - transpose(get_array(HESS(f, CS))))
    check("GRAD(fg) = f GRAD g + g GRAD f", name, GRAD(f * g, CS) - (f * GRAD(g, CS) + g * GRAD(f, CS)))
    check("∇²(fg) = f∇²g + 2∇f·∇g + g∇²f", name,
        LAPLACE(f * g, CS) - (f * LAPLACE(g, CS) + 2 * (GRAD(f, CS) ⋅ GRAD(g, CS)) + g * LAPLACE(f, CS)))
end

println("\n── vector identities ──")
for (name, CS) in charts
    d = get_dim(CS)
    q = getcoords(CS)
    f = SymFunction("f", real = true)(q...)
    𝐯 = sum(SymFunction("v$i", real = true)(q...) * unitvec(CS, i) for i in 1:d)
    𝟏 = tens_Id2(Val(d), Val(Sym))

    check("DIV(f𝐯) = GRAD f ⋅ 𝐯 + f DIV 𝐯", name,
        DIV(f * 𝐯, CS) - (GRAD(f, CS) ⋅ 𝐯 + f * DIV(𝐯, CS)))
    check("DIV(f𝟏) = GRAD f", name, DIV(f * 𝟏, CS) - GRAD(f, CS))
    check("DIV(GRAD 𝐯) = ∇²𝐯", name, DIV(GRAD(𝐯, CS), CS) - LAPLACE(𝐯, CS))
    check("tr(GRAD 𝐯) = DIV 𝐯", name, tr(GRAD(𝐯, CS)) - DIV(𝐯, CS))
end

println("\n── Leibniz rule for ∂ ──")
for (name, CS) in charts
    d = get_dim(CS)
    q = getcoords(CS)
    𝐚 = sum(SymFunction("a$i", real = true)(q...) * unitvec(CS, i) for i in 1:d)
    𝐛 = sum(SymFunction("b$i", real = true)(q...) * unitvec(CS, i) for i in 1:d)
    i = 1

    check("∂(𝐚⊗𝐛) = ∂𝐚⊗𝐛 + 𝐚⊗∂𝐛", name,
        ∂(𝐚 ⊗ 𝐛, i, CS) - (∂(𝐚, i, CS) ⊗ 𝐛 + 𝐚 ⊗ ∂(𝐛, i, CS)))
    check("∂(𝐚⊗ˢ𝐛) = ∂𝐚⊗ˢ𝐛 + 𝐚⊗ˢ∂𝐛", name,
        ∂(𝐚 ⊗ˢ 𝐛, i, CS) - (∂(𝐚, i, CS) ⊗ˢ 𝐛 + 𝐚 ⊗ˢ ∂(𝐛, i, CS)))
    check("∂(𝐚⋅𝐛) = ∂𝐚⋅𝐛 + 𝐚⋅∂𝐛", name,
        ∂(𝐚 ⋅ 𝐛, i, CS) - (∂(𝐚, i, CS) ⋅ 𝐛 + 𝐚 ⋅ ∂(𝐛, i, CS)))
end

println("\n── order-2 fields ──")
for (name, CS) in charts
    d = get_dim(CS)
    q = getcoords(CS)
    𝟏 = tens_Id2(Val(d), Val(Sym))
    𝐯 = sum(SymFunction("v$i", real = true)(q...) * unitvec(CS, i) for i in 1:d)
    𝛔 = sum(
        SymFunction("s$(min(i, j))$(max(i, j))", real = true)(q...) *
            unitvec(CS, i) ⊗ˢ unitvec(CS, j) for i in 1:d, j in 1:d
    )

    check("∂𝟏 = 0  (metric compatibility)", name, ∂(𝟏, 1, CS))
    check("DIV(f𝟏) with f=1 ⟹ 0", name, DIV(𝟏, CS))
    check("DIV(𝛔⋅𝐯) = DIV𝛔⋅𝐯 + 𝛔:GRAD𝐯", name,
        DIV(𝛔 ⋅ 𝐯, CS) - (DIV(𝛔, CS) ⋅ 𝐯 + 𝛔 ⊡ GRAD(𝐯, CS)))
end

println("\n── harmonic functions survive any chart ──")
harmonics = [
    "cartesian" => (coorsys_cartesian(), (q -> q[1]^2 - q[2]^2)),
    "polar" => (coorsys_polar(), (q -> q[1]^2 * cos(2q[2]))),
    "spherical" => (coorsys_spherical(), (q -> 1 / q[3])),
    "non-orthogonal" => (charts[4][2], (q -> (q[1] + q[2]^2)^2 - q[2]^2)),
]
for (name, (CS, f)) in harmonics
    check("∇²(harmonic) = 0", name, LAPLACE(f(getcoords(CS)), CS))
end

Spheroidal = coorsys_spheroidal()
ϕs, ps, qs = getcoords(Spheroidal)
P = sympy.assoc_legendre

println("\n── spheroidal harmonics ──")
for (n, m) in ((2, 0), (2, 1))
    T = P(n, m, ps) * P(n, m, qs) * cos(m * ϕs)
    check("∇²[P_$(n)^$(m)(p)P_$(n)^$(m)(q)cos($(m)φ)] = 0", "spheroidal",
        simplify(LAPLACE(T, Spheroidal)))
end

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
