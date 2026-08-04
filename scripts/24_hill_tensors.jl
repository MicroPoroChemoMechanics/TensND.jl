# # Acoustic and Hill tensors
#
# From the stiffness ``\mathbb{C}`` and a direction ``\underline{\xi}``, the
# **acoustic tensor** ``\boldsymbol{K}=\underline{\xi}\cdot\mathbb{C}\cdot
# \underline{\xi}`` controls wave propagation and, through its inverse, the
# Green operator
#
# ```math
# \mathbb{\Gamma}=\underline{\xi}\stackrel{s}{\otimes}\boldsymbol{K}^{-1}
#                 \stackrel{s}{\otimes}\underline{\xi},
# \qquad
# \mathbb{\Lambda}=\mathbb{C}:\mathbb{\Gamma}:\mathbb{C},
# ```
#
# the integrand of every Hill polarization tensor and the kernel of
# stress-intensity-factor integrals on crack fronts.
#
# The point of this tutorial is that the construction is written **once** and
# runs unchanged on an isotropic and on a transversely isotropic stiffness — the
# structured types of [The Walpole basis](@ref th-walpole) doing the work.
# Background: [mura1987](@cite), [hoenig1978](@cite).

import Pkg                                                          #jl
Pkg.activate(joinpath(@__DIR__, "..", "docs"); io = devnull)        #jl
## The `docs` environment declares every dependency the tutorials use     #jl
## (SymPy, Symbolics, NLopt, BenchmarkTools, ...) and resolves TensND      #jl
## itself through `[sources]`. These four lines are stripped from the      #jl
## generated page and notebook, which already run inside that project.     #jl

using TensND
using LinearAlgebra
using SymPy

# ## Setting up
#
# The direction ``\underline{\xi}`` is taken as the position vector of a
# spherical chart, so that ``\underline{\xi}=\underline{e}^r`` up to its norm
# and the angular dependence is explicit.

Cartesian = coorsys_cartesian(symbols("x y z", real = true))
𝐞₁, 𝐞₂, 𝐞₃ = unitvec(Cartesian)

Spherical = coorsys_spherical((symbols("θ ϕ", real = true)..., symbols("ξ", positive = true)))
θ, ϕ, ξ = getcoords(Spherical)
𝐞ᶿ, 𝐞ᵠ, 𝐞ʳ = unitvec(Spherical)

𝕀, 𝕁, 𝕂 = iso_projectors(Val(3), Val(Sym))
𝟏 = tens_Id2(Val(3), Val(Sym))
𝛏 = getOM(Spherical)

# ## The isotropic case
#
# With ``\mathbb{C}=3\lambda\mathbb{J}+2\mu\mathbb{I}``:

λ = symbols("λ", real = true)
μ = symbols("μ", positive = true)
ℂ = 3λ * 𝕁 + 2μ * 𝕀

𝐊 = 𝛏 ⋅ ℂ ⋅ 𝛏

# The acoustic tensor is transversely isotropic about ``\underline{\xi}``: one
# longitudinal eigenvalue and a doubly degenerate transverse one.

tsimplify(𝐊)

# The Green operator and the Hill-type kernel:

ℾ = 𝛏 ⊗ˢ 𝐊^(-1) ⊗ˢ 𝛏
𝚲 = tsimplify(ℂ ⊡ ℾ ⊡ ℂ)

# It must reproduce the closed form
#
# ```math
# \mathbb{\Lambda}=\frac{\lambda^{2}}{\lambda+2\mu}\,\boldsymbol{1}\otimes\boldsymbol{1}
# +\frac{2\lambda\mu}{\lambda+2\mu}\bigl(\boldsymbol{1}\otimes\underline{e}^r\otimes\underline{e}^r
#  +\underline{e}^r\otimes\underline{e}^r\otimes\boldsymbol{1}\bigr)
# +4\mu\Bigl(\underline{e}^r\stackrel{s}{\otimes}\boldsymbol{1}\stackrel{s}{\otimes}\underline{e}^r
#  -\frac{\lambda+\mu}{\lambda+2\mu}\,\underline{e}^r{}^{\otimes4}\Bigr)
# ```

𝚲₂ = tsimplify(
    λ^2 / (λ + 2μ) * 𝟏 ⊗ 𝟏
        + 2λ * μ / (λ + 2μ) * (𝟏 ⊗ 𝐞ʳ ⊗ 𝐞ʳ + 𝐞ʳ ⊗ 𝐞ʳ ⊗ 𝟏)
        + 4μ * (𝐞ʳ ⊗ˢ 𝟏 ⊗ˢ 𝐞ʳ - (λ + μ) / (λ + 2μ) * 𝐞ʳ ⊗ 𝐞ʳ ⊗ 𝐞ʳ ⊗ 𝐞ʳ)
)

intrinsic(tsimplify(𝚲 - 𝚲₂), Spherical)

# Identically zero.

# ## The transversely isotropic case
#
# The same three lines, with a stiffness built by [`tens_TI`](@ref) instead.
# Nothing in the construction changes: `⋅`, `inv` and `⊡` dispatch on the
# structured type.
#
# !!! note "This replaces a legacy hand-rolled Walpole basis"
#     An earlier research script re-defined the Walpole tensors and the TI
#     constructors locally. Everything it needed is now in the library —
#     [`walpole_basis`](@ref), [`tens_TI`](@ref), [`tens_TI_eng`](@ref),
#     [`tens_TI_Hoenig`](@ref) — so the local definitions are gone.

C₁₁₁₁, C₁₁₂₂, C₁₁₃₃, C₃₃₃₃, C₂₃₂₃ = symbols("C₁₁₁₁ C₁₁₂₂ C₁₁₃₃ C₃₃₃₃ C₂₃₂₃", positive = true)
n = 𝐞₃
ℂᵗⁱ = tens_TI(C₁₁₁₁, C₁₁₂₂, C₁₁₃₃, C₃₃₃₃, C₂₃₂₃, [Sym(0), Sym(0), Sym(1)])

typeof(ℂᵗⁱ), get_ℓ(ℂᵗⁱ)

# The acoustic tensor along the symmetry axis. Taking
# ``\underline{\xi}=\underline{e}_3`` makes it diagonal, with the longitudinal
# modulus ``C_{3333}`` and the doubly degenerate shear modulus ``C_{2323}``:

𝐊ᵃˣ = tsimplify(𝐞₃ ⋅ ℂᵗⁱ ⋅ 𝐞₃)
get_array(𝐊ᵃˣ)

# In the isotropy plane, ``\underline{\xi}=\underline{e}_1``, the three
# eigenvalues are all distinct — the anisotropy is fully visible:

𝐊ᵗ = tsimplify(𝐞₁ ⋅ ℂᵗⁱ ⋅ 𝐞₁)
get_array(𝐊ᵗ)

# The two coincide only when the material is isotropic. Substituting the
# isotropic relations ``C_{1111}=C_{3333}=\lambda+2\mu``,
# ``C_{1122}=C_{1133}=\lambda``, ``C_{2323}=\mu``:

iso_subs = Dict(
    C₁₁₁₁ => λ + 2μ, C₃₃₃₃ => λ + 2μ,
    C₁₁₂₂ => λ, C₁₁₃₃ => λ, C₂₃₂₃ => μ,
)
(tsimplify(subs.(get_array(𝐊ᵃˣ), iso_subs...)), tsimplify(subs.(get_array(𝐊ᵗ), iso_subs...)))

# ## The Green operator for the TI medium
#
# Along the symmetry axis the acoustic tensor is diagonal, so its inverse is
# immediate and ``\mathbb{\Lambda}`` follows:

ℾᵗⁱ = 𝐞₃ ⊗ˢ inv(𝐊ᵃˣ) ⊗ˢ 𝐞₃
𝚲ᵗⁱ = tsimplify(ℂᵗⁱ ⊡ ℾᵗⁱ ⊡ ℂᵗⁱ)

get_array(𝚲ᵗⁱ)[3, 3, 3, 3]

# The ``3333`` component is ``C_{3333}`` itself: along the symmetry axis the
# Green operator exactly undoes the stiffness, as it must for a longitudinal
# wave.

# ## Checking against the isotropic limit
#
# Substituting the isotropic moduli into the TI result must reproduce the
# isotropic ``\mathbb{\Lambda}`` evaluated at ``\underline{\xi}=\underline{e}_3``:

Λᵗⁱ_iso = tsimplify(subs.(get_array(𝚲ᵗⁱ), iso_subs...))
Λ_iso_axis = tsimplify(subs.(get_array(𝚲₂), θ => Sym(0), ϕ => Sym(0)))

tsimplify(Λᵗⁱ_iso - Λ_iso_axis)

# Zero: the transversely isotropic construction degenerates correctly.
