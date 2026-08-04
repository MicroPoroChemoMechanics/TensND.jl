import Pkg
Pkg.activate(joinpath(@__DIR__, "..", "docs"); io = devnull)
# The `docs` environment declares every dependency the tutorials use
# (SymPy, Symbolics, NLopt, BenchmarkTools, ...) and resolves TensND
# itself through `[sources]`. These four lines are stripped from the
# generated page and notebook, which already run inside that project.

using TensND
using LinearAlgebra
using SymPy

R = symbols("R", positive = true)
θ = symbols("θ", positive = true)
ϕ = symbols("ϕ", real = true)
z = symbols("z", real = true)

polar_rules = Dict(abs(sin(θ)) => sin(θ))

Sphere = SubManifoldSym(
    Tens(R * [sin(θ)cos(ϕ), sin(θ)sin(ϕ), cos(θ)]), (θ, ϕ), (), (R,);
    rules = polar_rules,
)

tsimplify.(components_canon(normal(Sphere)))

𝐚 = submetric(Sphere)
get_array(𝐚)

𝐛 = curvature(Sphere)
tsimplify.(get_array(𝐛))

tsimplify.(get_array(𝐛) + get_array(𝐚) / R)

function curvatures(SM)
    dim = length(get_array(submetric(SM))[1, :])
    tang = 1:(dim - 1)
    a = tsimplify.(get_array(submetric(SM))[tang, tang])
    b = tsimplify.(get_array(curvature(SM))[tang, tang])
    S = tsimplify.(inv(a) * b)
    return (K = tsimplify(det(S)), H = tsimplify(tr(S) / (dim - 1)), shape = tsimplify.(diag(S)))
end

cs = curvatures(Sphere)
println("sphere : K = ", cs.K, "   H = ", cs.H, "   principal curvatures = ", cs.shape)

Cylinder = SubManifoldSym(Tens([R * cos(ϕ), R * sin(ϕ), z]), (ϕ, z), (), (R,))

tsimplify.(get_array(submetric(Cylinder))), tsimplify.(get_array(curvature(Cylinder)))

cc = curvatures(Cylinder)
println("cylinder : K = ", cc.K, "   H = ", cc.H, "   principal curvatures = ", cc.shape)

x, y = symbols("x y", real = true)
Plane = SubManifoldSym(Tens([x, y, zero(x)]), (x, y))

tsimplify.(get_array(curvature(Plane)))

cp = curvatures(Plane)
println("plane : K = ", cp.K, "   H = ", cp.H)

Paraboloid = SubManifoldSym(Tens([x, y, (x^2 + y^2) / (2R)]), (x, y), (), (R,))
pb = curvatures(Paraboloid)

println("  s = √(x²+y²)      K·R²           H·R")
for s in (0.0, 0.5, 1.0, 2.0, 4.0)
    sub = Dict(x => Sym(s), y => Sym(0), R => Sym(1))
    Kv = Float64(simplify(pb.K.subs(sub)))
    Hv = Float64(simplify(pb.H.subs(sub)))
    println("     ", rpad(s, 12), rpad(round(Kv, digits = 8), 15), round(Hv, digits = 8))
end

(simplify(pb.K.subs(Dict(x => Sym(0), y => Sym(0)))), simplify(pb.H.subs(Dict(x => Sym(0), y => Sym(0)))))

G = GRAD(normal(Sphere), Sphere)

tsimplify.(get_array(G)), tsimplify.(get_array(𝐛))

ℬnat = natural_basis(Sphere)
for var in ((:cov, :cov), (:cont, :cov), (:cov, :cont), (:cont, :cont))
    resid = tsimplify.(components(G, ℬnat, var) + components(𝐛, ℬnat, var))
    println("  variance ", var, " : all zero ? ", all(iszero, resid))
end

Γ = connection(Sphere)
println("Γ^θ_ϕϕ = ", tsimplify(Γ[2, 2, 1]))
println("Γ^ϕ_θϕ = ", tsimplify(Γ[1, 2, 2]), "     Γ^ϕ_ϕθ = ", tsimplify(Γ[2, 1, 2]))

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
