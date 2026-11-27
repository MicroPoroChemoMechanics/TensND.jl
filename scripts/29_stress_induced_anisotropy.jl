# # [Stress-induced anisotropy](@id tut-stress-anisotropy)
#
# !!! tip "See also"
#     **Theory**: [Third-order elasticity](@ref th-third-order) · [Plane waves](@ref th-plane-waves)
#
#     **Manual**: [Plane waves](@ref man-plane-waves)
#
#     **Tutorials**: [Thomsen's parameters and weak anisotropy](@ref tut-thomsen) · [Weak anisotropy of any symmetry](@ref tut-generalized-thomsen)
#
# An isotropic rock under a uniaxial stress carries waves whose velocities
# depend on their direction: the stress makes it anisotropic. The effect is
# governed by the third-order elastic constants, which in rocks are hundreds
# or thousands of times larger than the second-order ones
# [johnson1996](@cite). This tutorial computes, for the materials of their
# Table 1, the first-order coefficients of the anisotropy, then the velocities
# of Barre granite under a uniaxial compression, and checks one of the slopes
# of [hughes1953](@citet) with symbols.

import Pkg                                                          #jl
Pkg.activate(joinpath(@__DIR__, "..", "docs"); io = devnull)        #jl

using TensND
using LinearAlgebra
using JSON
using Printf
using Plots
using SymPyPythonCall: symbols
gr()
nothing #hide

# ## A solid isotropic at orders two and three
#
# The strain energy per unit reference volume of an isotropic solid, cut at
# the third order in the Green–Lagrange strain
# ``\boldsymbol{E} = \tfrac12(\boldsymbol{F}^{\mathsf T}\cdot\boldsymbol{F} - \boldsymbol{1})``,
# has Lamé's constants ``\lambda``, ``\mu`` and Murnaghan's ``l``, ``m``, ``n``:
#
# ```math
# W = \tfrac12\lambda\,(\operatorname{tr}\boldsymbol{E})^2 + \mu\operatorname{tr}\boldsymbol{E}^2
#   + \tfrac13(l - m)\,(\operatorname{tr}\boldsymbol{E})^3 + m\operatorname{tr}\boldsymbol{E}\operatorname{tr}\boldsymbol{E}^2 + n\det\boldsymbol{E}
# ```
#
# ([Third-order elasticity](@ref th-third-order)). Table 1 of
# [johnson1996](@citet) gives, for metals, a glass and rocks, the bulk modulus
# ``K = \lambda + 2\mu/3``, the shear modulus ``\mu`` and ``l``, ``m``, ``n``,
# in GPa. It is stored with its source under `data/literature`.

data = JSON.parsefile(joinpath(pkgdir(TensND), "data", "literature", "johnsonJGR1996.json"))
table = data["tables"]["table_1"]
function solid(r)
    K, μ = r["K (GPa)"], r["mu (GPa)"]
    return ThirdOrderIsotropic(K - 2μ / 3, μ, r["l (GPa)"], r["m (GPa)"], r["n (GPa)"])
end
row(name) = only(filter(r -> r["material"] == name, table))
barre = solid(row("Barre granite"))

#md Main.show_literature("johnsonJGR1996") # hide

# The same constants in the two other usual notations, the components
# ``C_{111}``, ``C_{112}``, ``C_{123}`` and Landau's ``A``, ``B``, ``C``:

third_order_constants(barre)

# ## The coefficients of the induced anisotropy
#
# Under a uniaxial compression ``\sigma`` (positive in compression), at the
# angle ``\theta`` between propagation and stress, the P velocity and the
# splitting of the two S waves follow, to first order in ``\sigma``,
#
# ```math
# \frac{v_{P}(0) - v_P(\theta)}{v_{P}(0)} \simeq -\sigma\,a_P\sin^2\theta,
# \qquad
# \frac{v_{S_1}(\theta) - v_{S_2}(\theta)}{v_{S_1}(\theta)} \simeq -\sigma\,b_S\sin^2\theta ,
# ```
#
# ``S_1`` being polarized in the plane of the stress and of the propagation.
# `acoustoelastic_coefficients` gives ``a_P``, ``b_S`` and the slopes of
# [hughes1953](@citet). They are computed here from the constants of each row,
# next to the values predicted in the table, in GPa⁻¹, with their relative
# gap:

@printf "%-40s %8s %8s %7s %8s %8s %7s\n" "material" "b_S" "table" "gap" "a_P" "table" "gap"
for r in table
    c = acoustoelastic_coefficients(solid(r))
    b, a = r["b_S predicted (1/GPa)"], r["a_P predicted (1/GPa)"]
    @printf "%-40s %8.3g %8.3g %6.1f%% %8.3g %8.3g %6.1f%%\n" r["material"] c.b_S b 100 * (c.b_S / b - 1) c.a_P a 100 * (c.a_P / a - 1)
end

# The rows agree within the rounding of the printed moduli, except that of the
# water-saturated, thermally cracked sandstone: there the closed forms,
# applied to the shear modulus ``\mu`` printed in the same row, give the
# values of the columns `b_S` and `a_P`, which differ from the predictions
# printed next to them. The measured ``b_S`` of that row, for reference:

row("Water-saturated F32, thermally cracked")["b_S measured (1/GPa)"]

# ## Natural and true velocities
#
# A wave of reference direction ``\underline{N}`` in a solid deformed by the
# gradient ``\boldsymbol{F}`` has a natural velocity ``W``, the reference
# length of a path over its travel time, which `prestressed_waves` computes
# from the tensor ``A_{ijkl} = \partial^2 W/\partial F_{ij}\partial F_{kl}``.
# A segment of reference direction ``\underline{N}`` is stretched by
# ``|\boldsymbol{F}\cdot\underline{N}|``, so that the true velocity, the deformed
# length over the same time, is ``v = |\boldsymbol{F}\cdot\underline{N}|\,W``.
#
# ## Barre granite under a uniaxial compression
#
# To first order the compression ``\sigma`` along ``\underline{e}_3`` strains
# the solid by Hooke's law,
# ``\boldsymbol{F} = \boldsymbol{1} + \mathbb{S} : (-\sigma\,\underline{e}_3\otimes\underline{e}_3)``,
# that is ``F_{33} = 1 - \sigma/E`` and ``F_{11} = F_{22} = 1 + \nu\sigma/E``
# with Young's modulus ``E`` and Poisson's ratio ``\nu``:

young_poisson(λ, μ) = (μ * (3λ + 2μ) / (λ + μ), λ / (2 * (λ + μ)))
function uniaxial_gradient(σ, E, ν)
    stretch = (1 + ν * σ / E, 1 + ν * σ / E, 1 - σ / E)
    return [i == j ? stretch[i] : zero(stretch[i]) for i in 1:3, j in 1:3]
end
E, ν = young_poisson(barre.λ, barre.μ)
F(σ) = uniaxial_gradient(σ, E, ν)
cauchy_stress(barre, F(0.01))

# At ``\sigma`` = 10 MPa = 0.01 GPa the Cauchy stress is the uniaxial
# compression, up to terms of second order. The true velocities against the
# stress, for four angles between propagation and stress (density in g/cm³,
# velocities in km/s):

ρ₀ = row("Barre granite")["rho (kg/m3)"] / 1000
function true_velocities(σ, θ)
    N = direction_from_angles(θ, 0.0)
    return norm(F(σ) * N) .* prestressed_waves(barre, ρ₀, F(σ), N).velocity
end
σs = range(0, 0.02; length = 41)
angles = (0, 30, 60, 90)
p = plot(; xlabel = "uniaxial compression σ (MPa)", ylabel = "true P velocity (km/s)", legend = :topleft, size = (640, 380), left_margin = 4Plots.mm, bottom_margin = 4Plots.mm)
for θ in angles
    plot!(p, 1000 .* σs, [true_velocities(σ, deg2rad(θ))[1] for σ in σs]; lw = 2, label = "θ = $(θ)°")
end
p

# The P anisotropy and the S splitting against the angle, at two stresses,
# next to the ``\sin^2\theta`` laws (dashed). The splitting compares two
# velocities of the same direction, where the factor
# ``|\boldsymbol{F}\cdot\underline{N}|`` cancels; the P anisotropy compares two
# directions, where it does not.

c = acoustoelastic_coefficients(barre)
θs = range(0, π / 2; length = 46)
pP = plot(; xlabel = "θ (°)", ylabel = "(v_P(0) − v_P(θ)) / v_P(0)", xticks = 0:15:90, legend = :topleft)
pS = plot(; xlabel = "θ (°)", ylabel = "(v_S1 − v_S2) / v_S1", xticks = 0:15:90, legend = :topleft)
for (i, σ) in enumerate((0.005, 0.02))
    v₀ = true_velocities(σ, 0.0)[1]
    v = [true_velocities(σ, θ) for θ in θs]
    label = "σ = $(round(Int, 1000σ)) MPa"
    plot!(pP, rad2deg.(θs), [(v₀ - w[1]) / v₀ for w in v]; lw = 2, color = i, label)
    plot!(pP, rad2deg.(θs), -σ * c.a_P .* sin.(θs) .^ 2; lw = 2, ls = :dash, color = i, label = "")
    plot!(pS, rad2deg.(θs), [(w[2] - w[3]) / w[2] for w in v]; lw = 2, color = i, label)
    plot!(pS, rad2deg.(θs), -σ * c.b_S .* sin.(θs) .^ 2; lw = 2, ls = :dash, color = i, label = "")
end
plot(pP, pS; layout = (1, 2), size = (980, 380), left_margin = 4Plots.mm, bottom_margin = 4Plots.mm)

# The first-order laws hold while the stress is small against the scales
# ``\mu(\lambda + \mu)/|m|`` and ``\mu^2/|n|`` set by the third-order
# constants, a few tens of MPa for this granite (values in GPa):

(scale_P = barre.μ * (barre.λ + barre.μ) / abs(barre.m), scale_S = barre.μ^2 / abs(barre.n))

# At 20 MPa the exact anisotropies already depart from them.

# ## A slope of Hughes and Kelly with symbols
#
# The same functions take SymPy symbols. The true P modulus
# ``\rho_0 v^2`` along the stress is the component ``A_{3333}`` times
# ``|\boldsymbol{F}\cdot\underline{e}_3|^2``; its derivative at
# ``\sigma = 0``, from `acoustoelastic_tensor`, is the slope `P_parallel` of
# `acoustoelastic_coefficients`:

λₛ, μₛ, lₛ, mₛ, nₛ, σₛ = symbols("λ μ l m n σ", real = true)
sₛ = ThirdOrderIsotropic(λₛ, μₛ, lₛ, mₛ, nₛ)
Fₛ = uniaxial_gradient(σₛ, young_poisson(λₛ, μₛ)...)
A = acoustoelastic_tensor(sₛ, Fₛ)
slope = tsimplify(tsubs(tdiff(A[3, 3, 3, 3] * Fₛ[3, 3]^2, σₛ), σₛ => 0))

# and the difference with the closed form reduces to zero:

tsimplify(slope - acoustoelastic_coefficients(sₛ).P_parallel)
