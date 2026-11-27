# ─────────────────────────────────────────────────────────────────────────────
#  third_order_elasticity.jl — a hyperelastic solid isotropic at orders two
#  and three, and the small waves it carries under a static deformation.
#
#  The strain energy per unit reference volume is a function of the
#  Green-Lagrange strain E = ½(FᵀF − 𝟏), cut at the third order (Murnaghan):
#
#       W = ½λ (tr E)² + μ tr E² + ⅓(l − m)(tr E)³ + m tr E tr E² + n det E ,
#
#  with Lamé's (λ, μ) and Murnaghan's (l, m, n). The second Piola-Kirchhoff
#  stress and its tangent follow by differentiation,
#
#       S = ∂W/∂E ,   ℂ(E) = ∂S/∂E ,
#
#  written here in closed form, polynomial in E (cof E, the derivative of
#  det E, is E² − (tr E) E + ½[(tr E)² − tr E²] 𝟏).
#
#  A small motion u superposed on a static deformation of gradient F obeys, in
#  the reference coordinates X,
#
#       ρ₀ ü_i = ∂/∂X_j (A_ijkl ∂u_k/∂X_l) ,
#       A_ijkl = ∂²W/∂F_ij ∂F_kl = δ_ik S_jl + F_ip F_kq ℂ_pjql(E) ,
#
#  so that a plane wave of reference direction N has ρ₀ W² U = Γ·U with
#  Γ_ik = N_j A_ijkl N_l. W is the natural velocity: a length of the
#  reference configuration over a travel time, what a measurement between two
#  faces of a sample gives when the reference length is kept. A has the major
#  symmetry, hence a symmetric Γ, but not the minor ones: the stress enters
#  through δ_ik S_jl.
# ─────────────────────────────────────────────────────────────────────────────

"""
    ThirdOrderIsotropic(λ, μ, l, m, n)

A hyperelastic solid isotropic at orders two and three: Lamé's constants
`λ`, `μ` and Murnaghan's third-order constants `l`, `m`, `n`, of the strain
energy

```math
W = \\tfrac12\\lambda\\,(\\operatorname{tr}\\boldsymbol{E})^2 + \\mu\\operatorname{tr}\\boldsymbol{E}^2
  + \\tfrac13(l - m)\\,(\\operatorname{tr}\\boldsymbol{E})^3 + m\\operatorname{tr}\\boldsymbol{E}\\operatorname{tr}\\boldsymbol{E}^2 + n\\det\\boldsymbol{E}
```

in the Green–Lagrange strain ``\\boldsymbol{E}``
([Third-order elasticity](@ref th-third-order)). Numbers, dual numbers and
symbols. See [`third_order_constants`](@ref) for the other notations.
"""
struct ThirdOrderIsotropic{T}
    λ::T
    μ::T
    l::T
    m::T
    n::T
end
ThirdOrderIsotropic(λ, μ, l, m, n) = ThirdOrderIsotropic(promote(λ, μ, l, m, n)...)

"""
    third_order_constants(s::ThirdOrderIsotropic) → (C111, C112, C123, A, B, C)

The third-order constants of `s` in the two other usual notations: the
components ``C_{111} = 2l + 4m``, ``C_{112} = 2l``, ``C_{123} = 2l - 2m + n`` of
the sixth-order tensor ``\\partial^3W/\\partial\\boldsymbol{E}^3`` in the two-index
notation, and Landau's ``A = n``, ``B = m - n/2``, ``C = l - m + n/2``.
"""
third_order_constants(s::ThirdOrderIsotropic) =
    (C111 = 2 * s.l + 4 * s.m, C112 = 2 * s.l, C123 = 2 * s.l - 2 * s.m + s.n, A = s.n, B = s.m - s.n / 2, C = s.l - s.m + s.n / 2)

"""
    strain_energy(s::ThirdOrderIsotropic, E) → W

The strain energy of `s` per unit reference volume for the Green–Lagrange
strain `E`, a symmetric 3×3 matrix of canonical components.
"""
strain_energy(s::ThirdOrderIsotropic, E::AbstractMatrix) = _energy(s, SymmetricTensor{2, 3}(E))

# The kernels work on the tensors of Tensors.jl, in the canonical basis, whose
# types are known statically; the public functions take and return matrices or
# `Tens`. With E symmetric, E·E is `tdot(E)` and tr(E²) is E ⊡ E.
function _energy(s::ThirdOrderIsotropic, E::SymmetricTensor{2, 3})
    I₁, I₂ = tr(E), E ⊡ E
    return s.λ * I₁^2 / 2 + s.μ * I₂ + (s.l - s.m) * I₁^3 / 3 + s.m * I₁ * I₂ + s.n * det(E)
end

# The derivative of det E, by Cayley-Hamilton.
_cofactor(E::SymmetricTensor{2, 3}) = tdot(E) - tr(E) * E + (tr(E)^2 - E ⊡ E) / 2 * one(E)

# S = ∂W/∂E.
function _second_piola(s::ThirdOrderIsotropic, E::SymmetricTensor{2, 3})
    I₁, I₂ = tr(E), E ⊡ E
    return (s.λ * I₁ + (s.l - s.m) * I₁^2 + s.m * I₂) * one(E) + 2 * (s.μ + s.m * I₁) * E + s.n * _cofactor(E)
end

# The Green-Lagrange strain of a deformation gradient: ½(Fᵀ·F − 𝟏).
_strain(F::Tensor{2, 3}) = (tdot(F) - one(tdot(F))) / 2

# ℂ = ∂S/∂E: the isotropic stiffness of λ, μ and the term linear in E, the
# sixth-order tensor ∂³W/∂E³ contracted with E (`strain_stiffening`).
_tangent(s::ThirdOrderIsotropic, E::SymmetricTensor{2, 3}) =
    s.λ * (one(E) ⊗ one(E)) + 2 * s.μ * otimesul(one(E), one(E)) + _stiffening(s, E)

function _stiffening(s::ThirdOrderIsotropic, E::SymmetricTensor{2, 3})
    𝟏, I₁ = one(E), tr(E)
    a = 2 * s.m - s.n
    return (2 * (s.l - s.m) + s.n) * I₁ * (𝟏 ⊗ 𝟏) + a * (𝟏 ⊗ E + E ⊗ 𝟏) + a * I₁ * otimesul(𝟏, 𝟏) + s.n * (otimesul(𝟏, E) + otimesul(E, 𝟏))
end

# 𝔸 = 𝟏 ⊠ S + (F ⊠ 𝟏) : ℂ : (Fᵀ ⊠ 𝟏).
function _acoustoelastic(s::ThirdOrderIsotropic, F::Tensor{2, 3})
    E = _strain(F)
    𝟏 = one(F)
    return otimesu(𝟏, _second_piola(s, E)) + otimesu(F, 𝟏) ⊡ _tangent(s, E) ⊡ otimesu(transpose(F), 𝟏)
end

"""
    tangent_stiffness(s::ThirdOrderIsotropic, E) → AbstractTens{4,3}

The tangent ``\\mathbb{C}(\\boldsymbol{E}) = \\partial^2 W/\\partial\\boldsymbol{E}^2`` of
the second Piola–Kirchhoff stress of `s` at the Green–Lagrange strain `E` (a
symmetric 3×3 matrix of canonical components):

```math
\\mathbb{C} = c_1\\,\\boldsymbol{1}\\otimes\\boldsymbol{1}
 + c_2\\,(\\boldsymbol{1}\\otimes\\boldsymbol{E} + \\boldsymbol{E}\\otimes\\boldsymbol{1})
 + c_3\\,\\boldsymbol{1}\\stackrel{s}{\\boxtimes}\\boldsymbol{1}
 + n\\,(\\boldsymbol{1}\\stackrel{s}{\\boxtimes}\\boldsymbol{E} + \\boldsymbol{E}\\stackrel{s}{\\boxtimes}\\boldsymbol{1}),
```

with ``c_1 = \\lambda + (2l - 2m + n)\\operatorname{tr}\\boldsymbol{E}``,
``c_2 = 2m - n`` and ``c_3 = 2\\mu + (2m - n)\\operatorname{tr}\\boldsymbol{E}``: the
isotropic stiffness of `λ`, `μ` plus a term linear in `E`. It has the minor
and major symmetries ([Third-order elasticity](@ref th-third-order)).
"""
tangent_stiffness(s::ThirdOrderIsotropic, E::AbstractMatrix) = Tens(_tangent(s, SymmetricTensor{2, 3}(E)))

"""
    strain_stiffening(s::ThirdOrderIsotropic, E) → AbstractTens{4,3}

The change of stiffness that the strain `E` (a symmetric 3×3 matrix of
canonical components) brings about through the third-order constants of `s`,
``\\partial^3 W/\\partial\\boldsymbol{E}^3 : \\boldsymbol{E}``:

```math
\\begin{aligned}
\\Delta\\mathbb{C} &= (2l - 2m + n)\\operatorname{tr}\\boldsymbol{E}\\;\\boldsymbol{1}\\otimes\\boldsymbol{1}
 + (2m - n)\\,(\\boldsymbol{1}\\otimes\\boldsymbol{E} + \\boldsymbol{E}\\otimes\\boldsymbol{1})\\\\
&\\quad + (2m - n)\\operatorname{tr}\\boldsymbol{E}\\;\\boldsymbol{1}\\stackrel{s}{\\boxtimes}\\boldsymbol{1}
 + n\\,(\\boldsymbol{1}\\stackrel{s}{\\boxtimes}\\boldsymbol{E} + \\boldsymbol{E}\\stackrel{s}{\\boxtimes}\\boldsymbol{1}) ,
\\end{aligned}
```

that is ``\\Delta C_{11} = C_{111}E_{11} + C_{112}(E_{22} + E_{33})`` and its
companions in the two-index notation for principal strains. It does not
depend on `λ` and `μ`: added to the stiffness of an anisotropic background, it
is the model of [prioul2004](@citet) for a rock whose anisotropy at rest is
not isotropic while its nonlinearity is ([Third-order elasticity](@ref th-third-order)).
"""
strain_stiffening(s::ThirdOrderIsotropic, E::AbstractMatrix) = Tens(_stiffening(s, SymmetricTensor{2, 3}(E)))

"""
    third_order_from_constants(λ, μ, C₁₁₁, C₁₁₂, C₁₂₃) → ThirdOrderIsotropic

The solid of Lamé's constants `λ`, `μ` and of the third-order constants given
in the two-index notation: Murnaghan's ``l = C_{112}/2``,
``m = (C_{111} - C_{112})/4`` and ``n = C_{123} - C_{112} + (C_{111} - C_{112})/2``,
the converse of [`third_order_constants`](@ref).
"""
third_order_from_constants(λ, μ, C₁₁₁, C₁₁₂, C₁₂₃) =
    ThirdOrderIsotropic(λ, μ, C₁₁₂ / 2, (C₁₁₁ - C₁₁₂) / 4, C₁₂₃ - C₁₁₂ + (C₁₁₁ - C₁₁₂) / 2)

"""
    acoustoelastic_tensor(s::ThirdOrderIsotropic, F) → AbstractTens{4,3}

The tensor ``\\mathbb{A} = \\partial^2W/\\partial\\boldsymbol{F}^2`` that governs a small
motion superposed on the static deformation of gradient `F` (a 3×3 matrix of
canonical components), in the reference coordinates:

```math
\\mathbb{A} = \\boldsymbol{1}\\boxtimes\\boldsymbol{S}
 + (\\boldsymbol{F}\\boxtimes\\boldsymbol{1}) : \\mathbb{C}(\\boldsymbol{E}) : (\\boldsymbol{F}^{\\mathsf T}\\boxtimes\\boldsymbol{1}),
\\qquad
A_{ijkl} = \\delta_{ik}\\,S_{jl} + F_{ip}\\,F_{kq}\\,\\mathbb{C}_{pjql} ,
```

with ``\\boldsymbol{E} = \\tfrac12(\\boldsymbol{F}^{\\mathsf T}\\cdot\\boldsymbol{F} - \\boldsymbol{1})``,
``\\boldsymbol{S}`` the second Piola–Kirchhoff stress and ``\\mathbb{C}`` its
tangent ([`tangent_stiffness`](@ref)). It has the major symmetry and, under
stress, not the minor ones. For ``\\boldsymbol{F} = \\boldsymbol{1}`` it is the
isotropic stiffness of `λ`, `μ`.
"""
acoustoelastic_tensor(s::ThirdOrderIsotropic, F::AbstractMatrix) = Tens(_acoustoelastic(s, Tensor{2, 3}(F)))

"""
    cauchy_stress(s::ThirdOrderIsotropic, F) → SMatrix{3,3}

The Cauchy stress ``\\boldsymbol{\\sigma} = \\boldsymbol{F}\\cdot\\boldsymbol{S}\\cdot\\boldsymbol{F}^{\\mathsf T}/\\det\\boldsymbol{F}``
of `s` under the deformation gradient `F`.
"""
function cauchy_stress(s::ThirdOrderIsotropic, F::AbstractMatrix)
    𝐅 = Tensor{2, 3}(F)
    return SMatrix{3, 3}(𝐅 ⋅ _second_piola(s, _strain(𝐅)) ⋅ transpose(𝐅) / det(𝐅))
end

# The two contractions of a tensor of order 4 without the minor symmetries that
# the waves need: (N·𝔸·N)_ik = N_j A_ijkl N_l, the acoustic tensor, and
# (U·𝔸·U)_jl = U_i A_ijkl U_k, whose product with N is ρ₀W times the group
# velocity. `dotdot` of Tensors.jl contracts the first and fourth indices,
# which is the same only for a minor-symmetric tensor.
_contract24(A::Tensor{4, 3}, n) = SMatrix{3, 3}(ntuple(q -> ((i, k) = (mod1(q, 3), cld(q, 3)); sum(n[j] * A[i, j, k, l] * n[l] for j in 1:3, l in 1:3)), Val(9)))
_contract13(A::Tensor{4, 3}, u) = SMatrix{3, 3}(ntuple(q -> ((j, l) = (mod1(q, 3), cld(q, 3)); sum(u[i] * A[i, j, k, l] * u[k] for i in 1:3, k in 1:3)), Val(9)))

"""
    prestressed_waves(s::ThirdOrderIsotropic, ρ₀, F, N) → PlaneWaves

The three small plane waves of reference direction `N` in the solid `s`, of
reference density `ρ₀`, under the static deformation of gradient `F`: natural
velocities (a reference length over a travel time), sorted by decreasing
value, and polarizations, the eigenpairs of
``\\boldsymbol{\\Gamma} = \\underline{N}\\cdot\\mathbb{A}\\cdot\\underline{N}``
([`acoustoelastic_tensor`](@ref)). The field `group` holds the natural group
velocities ``g_j = U_i\\,A_{ijkl}\\,U_k\\,N_l/(\\rho_0 W)``. Numbers and dual
numbers.
"""
function prestressed_waves(s::ThirdOrderIsotropic, ρ₀, F::AbstractMatrix, N)
    A = _acoustoelastic(s, Tensor{2, 3}(F))
    T = float_like(promote_type(eltype(A), typeof(ρ₀), map(typeof, _direction_components(N))...))
    n = _unit_direction(N, T)
    Γ = _contract24(A, n)
    λ, U = _eigen3(Γ / convert(T, ρ₀))
    primal(λ[1]) > 0 || _throw_unstable(λ[1])
    V = SVector{3, T}(sqrt(λ[3]), sqrt(λ[2]), sqrt(λ[1]))
    P = SMatrix{3, 3, T}(U[:, 3]..., U[:, 2]..., U[:, 1]...)
    g = map(m -> _contract13(A, P[:, m]) * n / (ρ₀ * V[m]), SVector(1, 2, 3))
    return PlaneWaves{T}(n, V, P, SMatrix{3, 3, T}(g[1]..., g[2]..., g[3]...))
end

"""
    acoustoelastic_coefficients(s::ThirdOrderIsotropic)

The first-order effect of a stress on the velocities of `s`, in closed form
(Johnson & Rasolofosaon 1996, Eqs. (18)–(20), after Hughes & Kelly 1953):
the derivatives of ``\\rho_0 v^2`` (reference density, true
velocity ``v``, a length of the deformed solid over a travel time) with respect
to a **compressive** stress ``\\sigma``, taken at ``\\sigma = 0``,

- `P_parallel`, `S_parallel`: uniaxial stress along the direction of propagation;
- `P_perpendicular`, `S_perpendicular_parallel`, `S_perpendicular_perpendicular`:
  uniaxial stress normal to it, the shear wave being polarized along the
  stress or normal to it;
- `P_hydrostatic`, `S_hydrostatic`: pressure ``\\sigma``;

and the two coefficients of the anisotropy that a uniaxial stress induces,
at the angle ``\\theta`` between propagation and stress,

```math
\\frac{v_{P\\parallel} - v_P(\\theta)}{v_{P\\parallel}} \\simeq -\\sigma\\,a_P\\sin^2\\theta,
\\qquad
\\frac{v_{S_1}(\\theta) - v_{S_2}(\\theta)}{v_{S_1}(\\theta)} \\simeq -\\sigma\\,b_S\\sin^2\\theta ,
```

`a_P` and `b_S` (``S_1`` polarized in the plane of the stress and of the
propagation, ``S_2`` normal to it). Derivation in
[Third-order elasticity](@ref th-third-order). Numbers, dual numbers and symbols.
"""
function acoustoelastic_coefficients(s::ThirdOrderIsotropic)
    λ, μ, l, m, n = s.λ, s.μ, s.l, s.m, s.n
    K3 = 3 * λ + 2 * μ               # 3K; p: stress along the propagation, q: normal to it
    Pp = -(2 * l + λ + (λ + μ) / μ * (4 * m + 4 * λ + 10 * μ)) / K3
    Sp = -(m + λ * n / (4 * μ) + 4 * λ + 4 * μ) / K3
    Pq = -(2 * l - 2 * λ / μ * (m + λ + 2 * μ)) / K3
    Sqp = -(m + λ * n / (4 * μ) + λ + 2 * μ) / K3
    Sqq = -(m - (λ + μ) / (2 * μ) * n - 2 * λ) / K3
    return (
        P_parallel = Pp, S_parallel = Sp, P_perpendicular = Pq,
        S_perpendicular_parallel = Sqp, S_perpendicular_perpendicular = Sqq,
        P_hydrostatic = Pp + 2 * Pq, S_hydrostatic = Sp + Sqp + Sqq,
        a_P = (2 * λ + 5 * μ + 2 * m) / (2 * μ * (λ + 2 * μ)), b_S = (4 * μ + n) / (8 * μ^2),
    )
end

export ThirdOrderIsotropic, third_order_constants, third_order_from_constants, strain_stiffening, strain_energy, tangent_stiffness, acoustoelastic_tensor, cauchy_stress, prestressed_waves, acoustoelastic_coefficients
