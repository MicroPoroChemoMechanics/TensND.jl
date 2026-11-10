# # [The operators in the classical systems](@id tut-operators-catalog)
#
# !!! tip "See also"
#     **Theory**: [Curvilinear differential calculus](@ref th-curvilinear) · [Classical coordinate systems](@ref th-predefined-systems)
#
#     **Manual**: [Coordinate systems (symbolic)](@ref man-coorsystems) · [Coordinate systems (numerical)](@ref man-coorsystems-num)
#
#     **Tutorials**: [Differential operators, symbolically](@ref tut-symbolic-operators) · [Christoffel symbols](@ref tut-christoffel) · [Differential operators, numerically](@ref tut-numerical-operators) · [Classical identities of the differential operators](@ref tut-operator-identities)
#
# A reader who needs the expression of an operator in a given system finds it
# here, computed by `TensND` on arbitrary fields and compared with the textbook
# formula written just above it. The formulas are those of
# [Classical coordinate systems](@ref th-predefined-systems), derived from the
# general ones of an orthogonal system
# ([Curvilinear differential calculus](@ref th-orthogonal-operators)); each
# comparison returns exactly zero. The systems are the polar, cylindrical and
# spherical ones; the last section evaluates the same operators at a point, by
# automatic differentiation, and compares them with the closed forms.

import Pkg                                                          #jl
Pkg.activate(joinpath(@__DIR__, "..", "docs"); io = devnull)        #jl
## The `docs` environment declares every dependency the tutorials use     #jl
## (SymPy, Symbolics, BenchmarkTools, ...) and resolves TensND             #jl
## itself through `[sources]`. These four lines are stripped from the      #jl
## generated page and notebook, which already run inside that project.     #jl

using TensND
using LinearAlgebra
using SymPyPythonCall
using Plots

# ## Polar coordinates
#
# The coordinates are ``(r,\theta)``, the unit vectors ``\underline{e}^r`` and
# ``\underline{e}^\theta``, the Lamé coefficients ``(1,r)``. A scalar field
# ``f(r,\theta)`` and a vector field
# ``\underline{v}=v_r\,\underline{e}^r+v_\theta\,\underline{e}^\theta`` of
# physical components ``v_r(r,\theta)`` and ``v_\theta(r,\theta)`` are arbitrary
# functions:

Polar = coorsys_polar()
r, θ = get_coords(Polar)
𝐞ʳ, 𝐞ᶿ = unitvec(Polar)
ℬ = normalized_basis(Polar)
f = SymFunction("f", real = true)(r, θ)
vᵣ, vθ = (SymFunction("v$i", real = true)(r, θ) for i in ("r", "θ"))
𝐯 = vᵣ * 𝐞ʳ + vθ * 𝐞ᶿ

# The gradient of ``f`` has the physical components
#
# ```math
# \underline{\nabla}f=\frac{\partial f}{\partial r}\,\underline{e}^r
# +\frac{1}{r}\,\frac{\partial f}{\partial\theta}\,\underline{e}^\theta :
# ```

∇f = components(grad(f, Polar), ℬ)

#-

simplify.(∇f - [diff(f, r), diff(f, θ) / r])

# The divergence of ``\underline{v}``,
#
# ```math
# \mathrm{div}\,\underline{v}=\frac{1}{r}\,\frac{\partial(r\,v_r)}{\partial r}
# +\frac{1}{r}\,\frac{\partial v_\theta}{\partial\theta} :
# ```

simplify(divergence(𝐯, Polar) - (diff(r * vᵣ, r) + diff(vθ, θ)) / r)

# The Laplacian,
#
# ```math
# \Delta f=\frac{1}{r}\,\frac{\partial}{\partial r}\Bigl(r\,\frac{\partial f}{\partial r}\Bigr)
# +\frac{1}{r^2}\,\frac{\partial^2f}{\partial\theta^2} :
# ```

Δf = simplify(laplacian(f, Polar))

#-

simplify(Δf - (diff(r * diff(f, r), r) / r + diff(f, θ, 2) / r^2))

# The strain tensor ``\boldsymbol{\varepsilon}=\mathrm{sym}\,\underline{\nabla}\underline{v}``
# has the physical components
#
# ```math
# \varepsilon_{rr}=\frac{\partial v_r}{\partial r},
# \qquad
# \varepsilon_{\theta\theta}=\frac{1}{r}\,\frac{\partial v_\theta}{\partial\theta}+\frac{v_r}{r},
# \qquad
# \varepsilon_{r\theta}=\frac{1}{2}\Bigl(\frac{1}{r}\,\frac{\partial v_r}{\partial\theta}
# +\frac{\partial v_\theta}{\partial r}-\frac{v_\theta}{r}\Bigr) :
# ```

ε = components(symgrad(𝐯, Polar), ℬ)
ε_formula = [
    diff(vᵣ, r) (diff(vᵣ, θ) / r + diff(vθ, r) - vθ / r) / 2
    (diff(vᵣ, θ) / r + diff(vθ, r) - vθ / r) / 2 diff(vθ, θ) / r + vᵣ / r
]
simplify.(ε - ε_formula)

# ## Cylindrical coordinates
#
# The coordinates are ``(r,\theta,z)`` and the Lamé coefficients ``(1,r,1)``.
# The fields depend on the three coordinates:

Cyl = coorsys_cylindrical()
rc, θc, zc = get_coords(Cyl)
𝐞ʳᶜ, 𝐞ᶿᶜ, 𝐞ᶻᶜ = unitvec(Cyl)
ℬᶜ = normalized_basis(Cyl)
fᶜ = SymFunction("f", real = true)(rc, θc, zc)
uᵣ, uθ, u_z = (SymFunction("u$i", real = true)(rc, θc, zc) for i in ("r", "θ", "z"))
𝐮 = uᵣ * 𝐞ʳᶜ + uθ * 𝐞ᶿᶜ + u_z * 𝐞ᶻᶜ

# The gradient, the divergence and the Laplacian are the polar ones completed
# by the axial terms,
#
# ```math
# \underline{\nabla}f=\frac{\partial f}{\partial r}\,\underline{e}^r
# +\frac{1}{r}\,\frac{\partial f}{\partial\theta}\,\underline{e}^\theta
# +\frac{\partial f}{\partial z}\,\underline{e}^z,
# \qquad
# \mathrm{div}\,\underline{u}=\frac{1}{r}\,\frac{\partial(r\,u_r)}{\partial r}
# +\frac{1}{r}\,\frac{\partial u_\theta}{\partial\theta}+\frac{\partial u_z}{\partial z},
# ```
#
# ```math
# \Delta f=\frac{1}{r}\,\frac{\partial}{\partial r}\Bigl(r\,\frac{\partial f}{\partial r}\Bigr)
# +\frac{1}{r^2}\,\frac{\partial^2f}{\partial\theta^2}+\frac{\partial^2f}{\partial z^2} :
# ```

∇fᶜ = components(grad(fᶜ, Cyl), ℬᶜ)
simplify.(∇fᶜ - [diff(fᶜ, rc), diff(fᶜ, θc) / rc, diff(fᶜ, zc)])

#-

div_formula = (diff(rc * uᵣ, rc) + diff(uθ, θc)) / rc + diff(u_z, zc)
simplify(divergence(𝐮, Cyl) - div_formula)

#-

Δ_formula = diff(rc * diff(fᶜ, rc), rc) / rc + diff(fᶜ, θc, 2) / rc^2 + diff(fᶜ, zc, 2)
simplify(laplacian(fᶜ, Cyl) - Δ_formula)

# The six components of the strain tensor of the displacement
# ``\underline{u}`` are
#
# ```math
# \begin{aligned}
# \varepsilon_{rr}&=\frac{\partial u_r}{\partial r},&
# \varepsilon_{\theta\theta}&=\frac{1}{r}\,\frac{\partial u_\theta}{\partial\theta}+\frac{u_r}{r},&
# \varepsilon_{zz}&=\frac{\partial u_z}{\partial z},\\
# \varepsilon_{r\theta}&=\frac{1}{2}\Bigl(\frac{1}{r}\,\frac{\partial u_r}{\partial\theta}
# +\frac{\partial u_\theta}{\partial r}-\frac{u_\theta}{r}\Bigr),&
# \varepsilon_{\theta z}&=\frac{1}{2}\Bigl(\frac{1}{r}\,\frac{\partial u_z}{\partial\theta}
# +\frac{\partial u_\theta}{\partial z}\Bigr),&
# \varepsilon_{rz}&=\frac{1}{2}\Bigl(\frac{\partial u_r}{\partial z}
# +\frac{\partial u_z}{\partial r}\Bigr) :
# \end{aligned}
# ```

εᶜ = components(symgrad(𝐮, Cyl), ℬᶜ)
εᶜ_rθ = (diff(uᵣ, θc) / rc + diff(uθ, rc) - uθ / rc) / 2
εᶜ_θz = (diff(u_z, θc) / rc + diff(uθ, zc)) / 2
εᶜ_rz = (diff(uᵣ, zc) + diff(u_z, rc)) / 2
εᶜ_formula = [
    diff(uᵣ, rc) εᶜ_rθ εᶜ_rz
    εᶜ_rθ diff(uθ, θc) / rc + uᵣ / rc εᶜ_θz
    εᶜ_rz εᶜ_θz diff(u_z, zc)
]
simplify.(εᶜ - εᶜ_formula)

# ## Spherical coordinates
#
# The coordinates are ordered ``(\theta,\varphi,r)``, the colatitude first and
# the radius last, so that the components of a vector are
# ``(v_\theta,v_\varphi,v_r)``. The Lamé coefficients are
# ``(r,\,r\sin\theta,\,1)``.

Sph = coorsys_spherical()
θs, ϕs, rs = get_coords(Sph)
𝐞ᶿˢ, 𝐞ᵠˢ, 𝐞ʳˢ = unitvec(Sph)
ℬˢ = normalized_basis(Sph)
fˢ = SymFunction("f", real = true)(θs, ϕs, rs)
wθ, wϕ, wᵣ = (SymFunction("w$i", real = true)(θs, ϕs, rs) for i in ("θ", "ϕ", "r"))
𝐰 = wθ * 𝐞ᶿˢ + wϕ * 𝐞ᵠˢ + wᵣ * 𝐞ʳˢ

# The gradient, in the order of the components,
#
# ```math
# \underline{\nabla}f=\frac{1}{r}\,\frac{\partial f}{\partial\theta}\,\underline{e}^\theta
# +\frac{1}{r\sin\theta}\,\frac{\partial f}{\partial\varphi}\,\underline{e}^\varphi
# +\frac{\partial f}{\partial r}\,\underline{e}^r :
# ```

∇fˢ = components(grad(fˢ, Sph), ℬˢ)
simplify.(∇fˢ - [diff(fˢ, θs) / rs, diff(fˢ, ϕs) / (rs * sin(θs)), diff(fˢ, rs)])

# The divergence,
#
# ```math
# \mathrm{div}\,\underline{w}=\frac{1}{r\sin\theta}\,\frac{\partial(\sin\theta\,w_\theta)}{\partial\theta}
# +\frac{1}{r\sin\theta}\,\frac{\partial w_\varphi}{\partial\varphi}
# +\frac{1}{r^2}\,\frac{\partial(r^2\,w_r)}{\partial r} :
# ```

divw = diff(sin(θs) * wθ, θs) / (rs * sin(θs)) + diff(wϕ, ϕs) / (rs * sin(θs)) + diff(rs^2 * wᵣ, rs) / rs^2
simplify(divergence(𝐰, Sph) - divw)

# The Laplacian,
#
# ```math
# \Delta f=\frac{1}{r^2\sin\theta}\,\frac{\partial}{\partial\theta}\Bigl(\sin\theta\,\frac{\partial f}{\partial\theta}\Bigr)
# +\frac{1}{r^2\sin^2\theta}\,\frac{\partial^2f}{\partial\varphi^2}
# +\frac{1}{r^2}\,\frac{\partial}{\partial r}\Bigl(r^2\,\frac{\partial f}{\partial r}\Bigr) :
# ```

Δˢ = diff(sin(θs) * diff(fˢ, θs), θs) / (rs^2 * sin(θs)) + diff(fˢ, ϕs, 2) / (rs^2 * sin(θs)^2) +
    diff(rs^2 * diff(fˢ, rs), rs) / rs^2
simplify(laplacian(fˢ, Sph) - Δˢ)

# The strain tensor of ``\underline{w}`` has the six physical components
#
# ```math
# \begin{aligned}
# \varepsilon_{\theta\theta}&=\frac{1}{r}\,\frac{\partial w_\theta}{\partial\theta}+\frac{w_r}{r},\qquad
# \varepsilon_{rr}=\frac{\partial w_r}{\partial r},\\
# \varepsilon_{\varphi\varphi}&=\frac{1}{r\sin\theta}\,\frac{\partial w_\varphi}{\partial\varphi}
# +\frac{w_\theta\cot\theta}{r}+\frac{w_r}{r},\\
# \varepsilon_{\theta\varphi}&=\frac{1}{2}\Bigl(\frac{1}{r\sin\theta}\,\frac{\partial w_\theta}{\partial\varphi}
# +\frac{1}{r}\,\frac{\partial w_\varphi}{\partial\theta}-\frac{w_\varphi\cot\theta}{r}\Bigr),\\
# \varepsilon_{\theta r}&=\frac{1}{2}\Bigl(\frac{1}{r}\,\frac{\partial w_r}{\partial\theta}
# +\frac{\partial w_\theta}{\partial r}-\frac{w_\theta}{r}\Bigr),\qquad
# \varepsilon_{\varphi r}=\frac{1}{2}\Bigl(\frac{1}{r\sin\theta}\,\frac{\partial w_r}{\partial\varphi}
# +\frac{\partial w_\varphi}{\partial r}-\frac{w_\varphi}{r}\Bigr) .
# \end{aligned}
# ```
#
# The components of a displacement that depends on the three coordinates are
# simplified one at a time, which keeps each block short:

εˢ = components(symgrad(𝐰, Sph), ℬˢ)
εˢ_θθ = diff(wθ, θs) / rs + wᵣ / rs
εˢ_ϕϕ = diff(wϕ, ϕs) / (rs * sin(θs)) + wθ * cot(θs) / rs + wᵣ / rs
εˢ_rr = diff(wᵣ, rs)
simplify(εˢ[1, 1] - εˢ_θθ), simplify(εˢ[2, 2] - εˢ_ϕϕ), simplify(εˢ[3, 3] - εˢ_rr)

#-

εˢ_θϕ = (diff(wθ, ϕs) / (rs * sin(θs)) + diff(wϕ, θs) / rs - wϕ * cot(θs) / rs) / 2
εˢ_θr = (diff(wᵣ, θs) / rs + diff(wθ, rs) - wθ / rs) / 2
εˢ_ϕr = (diff(wᵣ, ϕs) / (rs * sin(θs)) + diff(wϕ, rs) - wϕ / rs) / 2
simplify(εˢ[1, 2] - εˢ_θϕ), simplify(εˢ[1, 3] - εˢ_θr), simplify(εˢ[2, 3] - εˢ_ϕr)

# The equilibrium equations of a body in spherical coordinates, without body
# force, state that the divergence of the symmetric stress tensor
# ``\boldsymbol{\sigma}`` vanishes. For an axisymmetric stress state, whose six
# physical components depend on ``(\theta,r)`` only, they read
#
# ```math
# \begin{aligned}
# (\mathrm{div}\,\boldsymbol{\sigma})_r&=\frac{\partial\sigma_{rr}}{\partial r}
# +\frac{1}{r}\,\frac{\partial\sigma_{r\theta}}{\partial\theta}
# +\frac{1}{r}\bigl(2\sigma_{rr}-\sigma_{\theta\theta}-\sigma_{\varphi\varphi}+\sigma_{r\theta}\cot\theta\bigr),\\
# (\mathrm{div}\,\boldsymbol{\sigma})_\theta&=\frac{\partial\sigma_{r\theta}}{\partial r}
# +\frac{1}{r}\,\frac{\partial\sigma_{\theta\theta}}{\partial\theta}
# +\frac{1}{r}\bigl((\sigma_{\theta\theta}-\sigma_{\varphi\varphi})\cot\theta+3\sigma_{r\theta}\bigr),\\
# (\mathrm{div}\,\boldsymbol{\sigma})_\varphi&=\frac{\partial\sigma_{r\varphi}}{\partial r}
# +\frac{1}{r}\,\frac{\partial\sigma_{\theta\varphi}}{\partial\theta}
# +\frac{1}{r}\bigl(3\sigma_{r\varphi}+2\sigma_{\theta\varphi}\cot\theta\bigr) .
# \end{aligned}
# ```
#
# The components of the divergence involve ``\sin2\theta`` and ``\tan\theta``,
# which `simplify` reduces only after the trigonometric functions of multiple
# angles are expanded:

σθθ, σϕϕ, σrr, σθϕ, σθr, σϕr = (
    SymFunction("σ$i", real = true)(θs, rs) for i in ("θθ", "ϕϕ", "rr", "θϕ", "θr", "ϕr")
)
𝛔 = σθθ * 𝐞ᶿˢ ⊗ 𝐞ᶿˢ + σϕϕ * 𝐞ᵠˢ ⊗ 𝐞ᵠˢ + σrr * 𝐞ʳˢ ⊗ 𝐞ʳˢ +
    σθϕ * (𝐞ᶿˢ ⊗ 𝐞ᵠˢ + 𝐞ᵠˢ ⊗ 𝐞ᶿˢ) + σθr * (𝐞ᶿˢ ⊗ 𝐞ʳˢ + 𝐞ʳˢ ⊗ 𝐞ᶿˢ) + σϕr * (𝐞ᵠˢ ⊗ 𝐞ʳˢ + 𝐞ʳˢ ⊗ 𝐞ᵠˢ)
div𝛔 = components(divergence(𝛔, Sph), ℬˢ)
eq_r = diff(σrr, rs) + diff(σθr, θs) / rs + (2σrr - σθθ - σϕϕ + σθr * cot(θs)) / rs
eq_θ = diff(σθr, rs) + diff(σθθ, θs) / rs + ((σθθ - σϕϕ) * cot(θs) + 3σθr) / rs
eq_ϕ = diff(σϕr, rs) + diff(σθϕ, θs) / rs + (3σϕr + 2σθϕ * cot(θs)) / rs
simplify.(expand.(div𝛔 - [eq_θ, eq_ϕ, eq_r], trig = true))

# ## The Hessian in polar coordinates
#
# The Hessian of a scalar field, ``\mathrm{grad}\,\mathrm{grad}\,f``, has the
# physical components
#
# ```math
# H_{rr}=\frac{\partial^2f}{\partial r^2},
# \qquad
# H_{r\theta}=\frac{1}{r}\,\frac{\partial^2f}{\partial r\,\partial\theta}-\frac{1}{r^2}\,\frac{\partial f}{\partial\theta},
# \qquad
# H_{\theta\theta}=\frac{1}{r}\,\frac{\partial f}{\partial r}+\frac{1}{r^2}\,\frac{\partial^2f}{\partial\theta^2} :
# ```

H = components(hess(f, Polar), ℬ)
H_formula = [
    diff(f, r, 2) diff(f, r, θ) / r - diff(f, θ) / r^2
    diff(f, r, θ) / r - diff(f, θ) / r^2 diff(f, r) / r + diff(f, θ, 2) / r^2
]
simplify.(H - H_formula)

# Its trace is the Laplacian:

simplify(tr(H) - laplacian(f, Polar))

# ## The same operators at a point
#
# A pointwise system, `coorsys_spherical(PointwiseGeometry)`, evaluates the
# operators of a field given as a function of the coordinates at one point,
# with derivatives by automatic differentiation. The field
# ``f=e^{-r^2}\,(1+\cos\theta+\sin^2\theta\cos2\varphi)`` is not harmonic. Its
# Laplacian is computed in closed form by the exact system, then evaluated at
# ``(\theta,\varphi,r)=(\pi/3,\pi/4,0.8)``, and by the pointwise system at
# the same point:

g = exp(-rs^2) * (1 + cos(θs) + sin(θs)^2 * cos(2ϕs))
Δg = simplify(laplacian(g, Sph))
x₀ = (π / 3, π / 4, 8 // 10)
Δg_exact = N(Δg(θs => x₀[1], ϕs => x₀[2], rs => x₀[3]))

#-

SphP = coorsys_spherical(PointwiseGeometry)
g_num(x) = exp(-x[3]^2) * (1 + cos(x[1]) + sin(x[1])^2 * cos(2x[2]))
Δg_point = laplacian(g_num, SphP)([π / 3, π / 4, 0.8])
Δg_point - Δg_exact

# The gradient evaluated in the same way gives the three physical components:

∇g_point = grad(g_num, SphP)([π / 3, π / 4, 0.8])
∇g_exact = [N(e(θs => x₀[1], ϕs => x₀[2], rs => x₀[3])) for e in simplify.(components(grad(g, Sph), ℬˢ))]
maximum(abs, ∇g_point - ∇g_exact)

# The Laplacian along the radius ``\theta=\pi/3``, ``\varphi=\pi/4``: the
# markers are the pointwise values, the line the closed form.

Δg_fun = lambdify(Δg, (θs, ϕs, rs))
radii = range(0.2, 2.0; length = 80)
rm = range(0.2, 2.0; length = 10)
plot(
    radii, [Δg_fun(π / 3, π / 4, ρ) for ρ in radii]; lw = 2, label = "closed form",
    xlabel = "r", ylabel = "Δf", size = (700, 420), guidefontsize = 9, legend = :bottomright,
)
scatter!(rm, [laplacian(g_num, SphP)([π / 3, π / 4, ρ]) for ρ in rm]; label = "pointwise", color = :black)
