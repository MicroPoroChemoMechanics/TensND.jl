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

_tr(E) = E[1, 1] + E[2, 2] + E[3, 3]
_tr2(E) = sum(E[i, j] * E[j, i] for i in 1:3, j in 1:3)
# x δ_ij, by selection: a Bool factor is not a number for every symbolic type.
_δ(i, j, x) = i == j ? x : zero(x)

"""
    strain_energy(s::ThirdOrderIsotropic, E) → W

The strain energy of `s` per unit reference volume for the Green–Lagrange
strain `E`, a symmetric 3×3 matrix of canonical components.
"""
function strain_energy(s::ThirdOrderIsotropic, E::AbstractMatrix)
    I₁, I₂ = _tr(E), _tr2(E)
    return s.λ * I₁^2 / 2 + s.μ * I₂ + (s.l - s.m) * I₁^3 / 3 + s.m * I₁ * I₂ + s.n * det(SMatrix{3, 3}(E))
end

# S_jl = ∂W/∂E_jl and ℂ_pjql = ∂S_pj/∂E_ql for a symmetric E, entry by entry.
function _second_piola(s::ThirdOrderIsotropic, E, j, l)
    I₁, I₂ = _tr(E), _tr2(E)
    E² = sum(E[j, a] * E[a, l] for a in 1:3)
    cof = E² - I₁ * E[j, l] + _δ(j, l, (I₁^2 - I₂) / 2)
    return _δ(j, l, s.λ * I₁ + (s.l - s.m) * I₁^2 + s.m * I₂) + 2 * (s.μ + s.m * I₁) * E[j, l] + s.n * cof
end

function _tangent(s::ThirdOrderIsotropic, E, p, j, q, l)
    I₁ = _tr(E)
    G = 2 * (s.μ + s.m * I₁) - s.n * I₁                               # factor of the symmetric identity
    𝕀 = (_δ(p, q, _δ(j, l, G)) + _δ(p, l, _δ(j, q, G))) / 2
    # derivative of E², symmetrized on (q, l); that of cof E follows
    dE² = (_δ(p, q, E[l, j]) + _δ(p, l, E[q, j]) + _δ(j, l, E[p, q]) + _δ(j, q, E[p, l])) / 2
    dcof = dE² - _δ(q, l, E[p, j]) + _δ(p, j, _δ(q, l, I₁) - E[q, l])
    return _δ(p, j, _δ(q, l, s.λ + 2 * (s.l - s.m) * I₁)) + 2 * s.m * (_δ(p, j, E[q, l]) + _δ(q, l, E[p, j])) + 𝕀 + s.n * dcof
end

const _VOIGT_PAIRS = ((1, 1), (2, 2), (3, 3), (2, 3), (1, 3), (1, 2))

"""
    tangent_stiffness(s::ThirdOrderIsotropic, E) → AbstractTens{4,3}

The tangent ``\\mathbb{C}(\\boldsymbol{E}) = \\partial^2 W/\\partial\\boldsymbol{E}^2`` of
the second Piola–Kirchhoff stress of `s` at the Green–Lagrange strain `E` (a
symmetric 3×3 matrix of canonical components): the isotropic stiffness of
`λ`, `μ` plus a term linear in `E`, the third-order tensor
``\\partial^3 W/\\partial\\boldsymbol{E}^3`` contracted with `E`. It has the minor
and major symmetries.
"""
function tangent_stiffness(s::ThirdOrderIsotropic, E::AbstractMatrix)
    T = promote_type(typeof(s.λ), eltype(E))
    M = SMatrix{6, 6, T}(
        ntuple(Val(36)) do q
            (p, j), (r, l) = _VOIGT_PAIRS[mod1(q, 6)], _VOIGT_PAIRS[cld(q, 6)]
            _tangent(s, E, p, j, r, l)
        end
    )
    return stiffness_from_voigt(M)
end

"""
    acoustoelastic_tensor(s::ThirdOrderIsotropic, F) → Array{T,4}

The tensor ``A_{ijkl} = \\partial^2W/\\partial F_{ij}\\partial F_{kl}`` that governs
a small motion superposed on the static deformation of gradient `F` (a 3×3
matrix of canonical components), in the reference coordinates:

```math
A_{ijkl} = \\delta_{ik}\\,S_{jl} + F_{ip}\\,F_{kq}\\,\\mathbb{C}_{pjql}(\\boldsymbol{E}),
\\qquad \\boldsymbol{E} = \\tfrac12(\\boldsymbol{F}^{\\mathsf T}\\boldsymbol{F} - \\boldsymbol{1}) ,
```

with ``\\boldsymbol{S}`` the second Piola–Kirchhoff stress and ``\\mathbb{C}`` its
tangent. It has the major symmetry and, under stress, not the minor ones. For
``\\boldsymbol{F} = \\boldsymbol{1}`` it is the isotropic stiffness of `λ`, `μ`.
"""
function acoustoelastic_tensor(s::ThirdOrderIsotropic, F::AbstractMatrix)
    E = (transpose(F) * F - I) / 2
    return [
        _δ(i, k, _second_piola(s, E, j, l)) + sum(F[i, p] * F[k, q] * _tangent(s, E, p, j, q, l) for p in 1:3, q in 1:3)
            for i in 1:3, j in 1:3, k in 1:3, l in 1:3
    ]
end

"""
    cauchy_stress(s::ThirdOrderIsotropic, F) → SMatrix{3,3}

The Cauchy stress ``\\boldsymbol{\\sigma} = \\boldsymbol{F}\\cdot\\boldsymbol{S}\\cdot\\boldsymbol{F}^{\\mathsf T}/\\det\\boldsymbol{F}``
of `s` under the deformation gradient `F`.
"""
function cauchy_stress(s::ThirdOrderIsotropic, F::AbstractMatrix)
    E = (transpose(F) * F - I) / 2
    S = SMatrix{3, 3}(ntuple(q -> _second_piola(s, E, mod1(q, 3), cld(q, 3)), Val(9)))
    Fs = SMatrix{3, 3}(F)
    return Fs * S * transpose(Fs) / det(Fs)
end

"""
    prestressed_waves(s::ThirdOrderIsotropic, ρ₀, F, N) → PlaneWaves

The three small plane waves of reference direction `N` in the solid `s`, of
reference density `ρ₀`, under the static deformation of gradient `F`: natural
velocities (a reference length over a travel time), sorted by decreasing
value, and polarizations, the eigenpairs of
``\\Gamma_{ik} = N_j A_{ijkl} N_l`` ([`acoustoelastic_tensor`](@ref)). The field
`group` holds the natural group velocities. Numbers and dual numbers.
"""
function prestressed_waves(s::ThirdOrderIsotropic, ρ₀, F::AbstractMatrix, N)
    A = acoustoelastic_tensor(s, F)
    T = float_like(promote_type(eltype(A), typeof(ρ₀), map(typeof, _direction_components(N))...))
    n = _unit_direction(N, T)
    Γ = SMatrix{3, 3, T}(ntuple(q -> ((i, k) = (mod1(q, 3), cld(q, 3)); sum(n[j] * A[i, j, k, l] * n[l] for j in 1:3, l in 1:3)), Val(9)))
    λ, U = _eigen3(Γ / convert(T, ρ₀))
    primal(λ[1]) > 0 || _throw_unstable(λ[1])
    V = SVector{3, T}(sqrt(λ[3]), sqrt(λ[2]), sqrt(λ[1]))
    P = SMatrix{3, 3, T}(U[:, 3]..., U[:, 2]..., U[:, 1]...)
    # g_j = A_ijkl U_i U_k N_l / (ρ₀ V)
    g = map(SVector(1, 2, 3)) do w
        SVector{3, T}(ntuple(j -> sum(A[i, j, k, l] * P[i, w] * P[k, w] * n[l] for i in 1:3, k in 1:3, l in 1:3) / (ρ₀ * V[w]), Val(3)))
    end
    return PlaneWaves{T}(n, V, P, SMatrix{3, 3, T}(g[1]..., g[2]..., g[3]...))
end

"""
    acoustoelastic_coefficients(s::ThirdOrderIsotropic)

The first-order effect of a stress on the velocities of `s`, in closed form
(Hughes & Kelly): the derivatives of ``\\rho_0 v^2`` (reference density, true
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

export ThirdOrderIsotropic, third_order_constants, strain_energy, tangent_stiffness, acoustoelastic_tensor, cauchy_stress, prestressed_waves, acoustoelastic_coefficients
