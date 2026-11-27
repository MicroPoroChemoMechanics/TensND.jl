# ─────────────────────────────────────────────────────────────────────────────
#  ortho_constants.jl — an orthotropic elasticity tensor from its engineering
#  constants, and back.
#
#  In the material frame (e₁, e₂, e₃) the compliance of an orthotropic solid
#  reads, with the engineering shears γ_ij = 2ε_ij,
#
#       ε₁₁ = σ₁₁/E₁ − ν₂₁σ₂₂/E₂ − ν₃₁σ₃₃/E₃   (and the two others) ,
#       γ₂₃ = σ₂₃/G₂₃ ,  γ₁₃ = σ₁₃/G₁₃ ,  γ₁₂ = σ₁₂/G₁₂ ,
#
#  with ν_ij/E_i = ν_ji/E_j (major symmetry): ν_ij is the contraction along
#  e_j under a tension along e_i. Nine constants in space,
#  (E₁, E₂, E₃, ν₁₂, ν₁₃, ν₂₃, G₂₃, G₁₃, G₁₂); four in the plane,
#  (E₁, E₂, ν₁₂, G₁₂). As tensor components, S₁₁₂₂ = −ν₁₂/E₁ and
#  S₂₃₂₃ = 1/(4G₂₃).
# ─────────────────────────────────────────────────────────────────────────────

"""
    ortho_from_constants(E₁, E₂, E₃, ν₁₂, ν₁₃, ν₂₃, G₂₃, G₁₃, G₁₂, frame; role = :compliance) -> TensOrtho
    ortho_from_constants(E₁, E₂, ν₁₂, G₁₂, frame; role = :compliance) -> TensOrtho

The orthotropic tensor of the engineering constants in the material frame
`frame`, of space (nine constants) or of the plane (four): Young's moduli
``E_i``, Poisson's ratios ``\\nu_{ij}`` (contraction along ``\\underline{e}_j``
under a tension along ``\\underline{e}_i``, so that
``\\nu_{ij}/E_i = \\nu_{ji}/E_j``) and shear moduli ``G_{ij}``. It is the
**compliance** of these constants,

```math
S_{iiii} = \\frac{1}{E_i},\\qquad S_{iijj} = -\\frac{\\nu_{ij}}{E_i},\\qquad S_{ijij} = \\frac{1}{4G_{ij}}\\quad(i \\ne j) ,
```

or, with `role = :stiffness`, the stiffness of the same material, its inverse.
Numbers, dual numbers and symbols. See [`ortho_constants`](@ref).
"""
ortho_from_constants(E₁, E₂, E₃, ν₁₂, ν₁₃, ν₂₃, G₂₃, G₁₃, G₁₂, frame::OrthonormalBasis{3}; role::Symbol = :compliance) =
    _of_role(
    _role(role),
    TensOrtho(inv(E₁), inv(E₂), inv(E₃), -ν₁₂ / E₁, -ν₁₃ / E₁, -ν₂₃ / E₂, inv(4 * G₂₃), inv(4 * G₁₃), inv(4 * G₁₂), frame)
)
ortho_from_constants(E₁, E₂, ν₁₂, G₁₂, frame::OrthonormalBasis{2}; role::Symbol = :compliance) =
    _of_role(_role(role), TensOrtho(inv(E₁), inv(E₂), -ν₁₂ / E₁, inv(4 * G₁₂), frame))

# The tensor of the wanted role from the compliance, and the compliance from it.
_of_role(::Val{:compliance}, S) = S
_of_role(::Val{:stiffness}, S) = inv(S)

"""
    ortho_constants(t::TensOrtho{4}; role = :compliance)

The engineering constants of the major-symmetric orthotropic tensor `t`, read
as a compliance or, with `role = :stiffness`, as a stiffness:
`(E₁, E₂, E₃, ν₁₂, ν₁₃, ν₂₃, G₂₃, G₁₃, G₁₂)` in space, `(E₁, E₂, ν₁₂, G₁₂)` in
the plane. The inverse of [`ortho_from_constants`](@ref).
"""
ortho_constants(t::TensOrtho{4}; role::Symbol = :compliance) = _ortho_engineering(_of_role(_role(role), t))

function _ortho_engineering(S::TensOrtho{4, 3, <:Any, 9})
    S₁₁, S₂₂, S₃₃, S₁₂, S₁₃, S₂₃, S₄₄, S₅₅, S₆₆ = get_data(S)
    return (inv(S₁₁), inv(S₂₂), inv(S₃₃), -S₁₂ / S₁₁, -S₁₃ / S₁₁, -S₂₃ / S₂₂, inv(4 * S₄₄), inv(4 * S₅₅), inv(4 * S₆₆))
end
function _ortho_engineering(S::TensOrtho{4, 2, <:Any, 4})
    S₁₁, S₂₂, S₁₂, S₆₆ = get_data(S)
    return (inv(S₁₁), inv(S₂₂), -S₁₂ / S₁₁, inv(4 * S₆₆))
end
_ortho_engineering(::TensOrtho{4}) = throw(
    ArgumentError("ortho_constants: the engineering constants are those of a major-symmetric tensor (9 constants in space, 4 in the plane)")
)

export ortho_from_constants, ortho_constants
