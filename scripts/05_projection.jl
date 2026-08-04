# # Projection and symmetry detection
#
# Given a tensor, what symmetry does it have — and if it has one, in what
# orientation? This tutorial runs [`proj_tens`](@ref) and
# [`best_sym_tens`](@ref) on tensors of known symmetry in known frames,
# with the axis fixed and with the axis optimized, and measures what each mode
# costs.
#
# Theory: [Projection onto a symmetry class](@ref th-projection).

import Pkg                                                          #jl
Pkg.activate(joinpath(@__DIR__, "..", "docs"); io = devnull)        #jl
## The `docs` environment declares every dependency the tutorials use     #jl
## (SymPy, Symbolics, NLopt, BenchmarkTools, ...) and resolves TensND      #jl
## itself through `[sources]`. These four lines are stripped from the      #jl
## generated page and notebook, which already run inside that project.     #jl

using TensND
using LinearAlgebra
using Printf
using Random

Random.seed!(20260804)
arr(x) = get_array(x)
𝕀, 𝕁, 𝕂 = ISO(Val(3), Val(Float64))

# ## The return contract
#
# [`proj_tens`](@ref) returns `(B, d, drel)`: the projected tensor, the absolute
# Frobenius distance and the relative one. The relative distance is what makes a
# tolerance dimensionless and independent of the units of the moduli.

nv = [0.0, 0.0, 1.0]
C_ti = tens_TI(10.0, 3.0, 2.5, 12.0, 2.0, nv)

B, d, drel = proj_tens(:TI, arr(C_ti), nv)
@printf "projecting a TI tensor onto TI (right axis): d = %.3e, drel = %.3e\n" d drel
println("returned type: ", typeof(B))

# Projecting a tensor onto a class it already belongs to is exact — the
# projection is the identity on the subspace.

# ## The wrong axis costs something
#
# The same tensor projected onto TI about axes progressively further from the
# true one:

for α in (0.0, 15.0, 30.0, 45.0, 90.0)
    m = [sin(deg2rad(α)), 0.0, cos(deg2rad(α))]
    _, _, dr = proj_tens(:TI, arr(C_ti), m)
    @printf "  axis tilted by %5.1f°  →  drel = %.4e\n" α dr
end

# At 90° the axis lies in the isotropy plane and the fit is worst; at 0° it is
# exact. This function of the orientation is what the optimizer minimizes.

# ## Fixed frame versus optimized frame
#
# Build a genuinely transversely isotropic tensor about a **tilted** axis, then
# hide that fact from the projector.

θ, ϕ = 0.6, 1.1
n_tilt = [sin(θ)cos(ϕ), sin(θ)sin(ϕ), cos(θ)]
C_tilt = tens_TI(10.0, 3.0, 2.5, 12.0, 2.0, n_tilt)

_, _, drel_canon = proj_tens(:TI, arr(C_tilt), [0.0, 0.0, 1.0])
_, _, drel_true = proj_tens(:TI, arr(C_tilt), n_tilt)
@printf "assuming the canonical axis : drel = %.4e\n" drel_canon
@printf "using the true axis         : drel = %.4e\n" drel_true

# Without `NLopt` the axis must be supplied. Loading it activates the extension
# and enables the no-axis methods, which search over the orientation.

using NLopt

B_opt, d_opt, drel_opt = proj_tens(:TI, arr(C_tilt))
@printf "optimized axis              : drel = %.4e\n" drel_opt
println("recovered axis : ", round.(collect(axis(B_opt)), digits = 8))
println("true axis      : ", round.(n_tilt, digits = 8))

# The axis is recovered up to a sign, which is immaterial: ``\underline{n}`` and
# ``-\underline{n}`` define the same symmetry.

min(
    norm(collect(axis(B_opt)) - n_tilt),
    norm(collect(axis(B_opt)) + n_tilt),
)

# ## Determinism
#
# The multi-start is deterministic — no clock-seeded global stage — so repeated
# calls return bit-identical results. This matters: the previous stochastic
# strategy produced intermittent CI failures.

results = [proj_tens(:TI, arr(C_tilt))[3] for _ in 1:5]
all(r === results[1] for r in results), results[1]

# ## Orthotropy in a rotated frame

ℬmat = Basis(0.4, 0.9, 0.3)
t_ortho = TensOrtho(10.0, 8.0, 12.0, 3.0, 2.5, 1.5, 2.0, 3.0, 3.5, ℬmat)

_, _, dr_fixed = proj_tens(:ORTHO, arr(t_ortho), CanonicalBasis{3, Float64}())
Bo, _, dr_opt = proj_tens(:ORTHO, arr(t_ortho))
@printf "ORTHO, canonical frame assumed : drel = %.4e\n" dr_fixed
@printf "ORTHO, frame optimized         : drel = %.4e\n" dr_opt

# The optimized frame recovers the material frame up to the ordering and signs
# of its axes — an orthotropic frame is defined only up to the 24 rotations
# mapping the three orthogonal planes onto themselves.

# ## The detection cascade
#
# [`best_sym_tens`](@ref) tries ISO → TI → ORTHO and returns the first class
# whose relative error falls below ``\varepsilon``.

# Each is wrapped as a **generic** tensor first, so that the detection has to
# discover the symmetry from the components rather than read it off a
# structured container.

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

# Without optimization the axis and frame are guessed from the Kelvin–Mandel
# eigenstructure, which is exact when the tensor really does have the symmetry —
# so the cheap path finds the tilted cases too, and needs no `NLopt`:

println("\nwithout angle optimization (cheap path):")
for (name, C) in materials
    _, _, dr, sym = best_sym_tens(C)
    @printf "  %-34s → %-6s (drel = %.3e)\n" name string(sym) dr
end

# ## Cost
#
# The cheap path is a handful of closed-form averages; the optimized one runs a
# multi-start local search over two or three angles.

using BenchmarkTools

C_generic = generic(C_tilt)
t_cheap = @belapsed best_sym_tens($C_generic)
t_opt = @belapsed best_sym_tens($C_generic; optimize_angles = true)
@printf "cheap path     : %8.1f µs\n" 1.0e6 * t_cheap
@printf "optimized path : %8.1f µs   (×%.0f)\n" 1.0e6 * t_opt t_opt / t_cheap

# Use the optimized path when the orientation is genuinely unknown and the
# tensor is not exactly of the class sought — that is, when the eigenstructure
# candidate is only a guess. Otherwise the cheap path gives the same answer for
# a fraction of the cost.

# ## Predicates
#
# The same computation with a boolean answer:

for (name, C) in materials
    @printf "  %-34s ISO=%-5s TI=%-5s ORTHO=%-5s\n" name string(is_ISO(C)) string(is_TI(C)) string(is_ORTHO(C))
end

# Note the hierarchy: an isotropic tensor satisfies all three predicates, since
# ISO ⊂ TI ⊂ ORTHO.
