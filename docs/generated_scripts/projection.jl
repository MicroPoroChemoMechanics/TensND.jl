import Pkg
Pkg.activate(joinpath(@__DIR__, "..", "docs"); io = devnull)
# The `docs` environment declares every dependency the tutorials use
# (SymPy, Symbolics, NLopt, BenchmarkTools, ...) and resolves TensND
# itself through `[sources]`. These four lines are stripped from the
# generated page and notebook, which already run inside that project.

using TensND
using LinearAlgebra
using Printf
using Random

Random.seed!(20260804)
arr(x) = get_array(x)
𝕀, 𝕁, 𝕂 = ISO(Val(3), Val(Float64))

nv = [0.0, 0.0, 1.0]
C_ti = tens_TI(10.0, 3.0, 2.5, 12.0, 2.0, nv)

B, d, drel = proj_tens(:TI, arr(C_ti), nv)
@printf "projecting a TI tensor onto TI (right axis): d = %.3e, drel = %.3e\n" d drel
println("returned type: ", typeof(B))

for α in (0.0, 15.0, 30.0, 45.0, 90.0)
    m = [sin(deg2rad(α)), 0.0, cos(deg2rad(α))]
    _, _, dr = proj_tens(:TI, arr(C_ti), m)
    @printf "  axis tilted by %5.1f°  →  drel = %.4e\n" α dr
end

θ, ϕ = 0.6, 1.1
n_tilt = [sin(θ)cos(ϕ), sin(θ)sin(ϕ), cos(θ)]
C_tilt = tens_TI(10.0, 3.0, 2.5, 12.0, 2.0, n_tilt)

_, _, drel_canon = proj_tens(:TI, arr(C_tilt), [0.0, 0.0, 1.0])
_, _, drel_true = proj_tens(:TI, arr(C_tilt), n_tilt)
@printf "assuming the canonical axis : drel = %.4e\n" drel_canon
@printf "using the true axis         : drel = %.4e\n" drel_true

using NLopt

B_opt, d_opt, drel_opt = proj_tens(:TI, arr(C_tilt))
@printf "optimized axis              : drel = %.4e\n" drel_opt
println("recovered axis : ", round.(collect(axis(B_opt)), digits = 8))
println("true axis      : ", round.(n_tilt, digits = 8))

min(
    norm(collect(axis(B_opt)) - n_tilt),
    norm(collect(axis(B_opt)) + n_tilt),
)

results = [proj_tens(:TI, arr(C_tilt))[3] for _ in 1:5]
all(r === results[1] for r in results), results[1]

ℬmat = Basis(0.4, 0.9, 0.3)
t_ortho = TensOrtho(10.0, 8.0, 12.0, 3.0, 2.5, 1.5, 2.0, 3.0, 3.5, ℬmat)

_, _, dr_fixed = proj_tens(:ORTHO, arr(t_ortho), CanonicalBasis{3, Float64}())
Bo, _, dr_opt = proj_tens(:ORTHO, arr(t_ortho))
@printf "ORTHO, canonical frame assumed : drel = %.4e\n" dr_fixed
@printf "ORTHO, frame optimized         : drel = %.4e\n" dr_opt

generic(t) = Tens(arr(t))

materials = [
    "isotropic" => generic(3 * 20.0 * 𝕁 + 2 * 8.0 * 𝕂),
    "transversely isotropic (e₃)" => generic(C_ti),
    "transversely isotropic (tilted)" => generic(C_tilt),
    "orthotropic (rotated frame)" => generic(t_ortho),
    "fully anisotropic" => generic(inv_KM((m -> (m + m') / 2)(rand(6, 6)))),
]

println("with angle optimization:")
for (name, C) in materials
    _, _, dr, sym = best_sym_tens(C; optimize_angles = true)
    @printf "  %-34s → %-6s (drel = %.3e)\n" name string(sym) dr
end

println("\nwithout angle optimization (cheap path):")
for (name, C) in materials
    _, _, dr, sym = best_sym_tens(C)
    @printf "  %-34s → %-6s (drel = %.3e)\n" name string(sym) dr
end

using BenchmarkTools

C_generic = generic(C_tilt)
t_cheap = @belapsed best_sym_tens($C_generic)
t_opt = @belapsed best_sym_tens($C_generic; optimize_angles = true)
@printf "cheap path     : %8.1f µs\n" 1.0e6 * t_cheap
@printf "optimized path : %8.1f µs   (×%.0f)\n" 1.0e6 * t_opt t_opt / t_cheap

for (name, C) in materials
    @printf "  %-34s ISO=%-5s TI=%-5s ORTHO=%-5s\n" name string(is_ISO(C)) string(is_TI(C)) string(is_ORTHO(C))
end

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
