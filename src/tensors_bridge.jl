# =============================================================================
#  tensors_bridge.jl — conversions to and from Tensors.jl, and the Voigt forms.
#
#  A Tensors.jl tensor has no basis: it is a bare component array, and in a
#  finite-element code those components are in the global frame. A TensND
#  tensor reads its components in its own basis (`get_array`, `t[i, j]`), so
#  every conversion below goes through the canonical basis explicitly; handing
#  the own-basis components of a tensor written in a rotated basis to an
#  assembler would rotate the material without a word.
# =============================================================================

export to_tensors, from_tensors, voigt_stress, voigt_strain, stress_from_voigt, strain_from_voigt

"""
    to_tensors(t::AbstractTens{2,dim}) -> Tensors.SymmetricTensor{2,dim}
    to_tensors(t::AbstractTens{4,dim}) -> Tensors.SymmetricTensor{4,dim}

The `Tensors.jl` counterpart of `t`, **in the canonical basis**: the components a
finite-element code (Ferrite among others) contracts with its shape-function
gradients.

A second-order tensor is symmetrized. A fourth-order one is assumed to have the
minor symmetries, as stiffnesses, compliances and concentration tensors do; the
major symmetry is not assumed, `SymmetricTensor{4,3}` holding 36 independent
components.

!!! warning "Not `get_array`"
    `get_array(t)` and `t[i, j, k, l]` are the components of `t` in its own
    basis. For a tensor written in a rotated basis they are not the canonical
    ones; `to_tensors` performs the change of basis.
"""
function to_tensors(t::AbstractTens{2, dim}) where {dim}
    A = components_canon(t)
    return SymmetricTensor{2, dim}((i, j) -> (A[i, j] + A[j, i]) / 2)
end

function to_tensors(t::AbstractTens{4, dim}) where {dim}
    A = components_canon(t)
    return SymmetricTensor{4, dim}((i, j, k, l) -> A[i, j, k, l])
end

"""
    from_tensors(t::Tensors.SymmetricTensor) -> Tens

The TensND tensor whose components in the canonical basis are those of `t`:
the inverse of [`to_tensors`](@ref), for a strain or a stress coming out of a
finite-element code.
"""
from_tensors(t::SymmetricTensor) = Tens(t)

"""
    voigt_stress(σ::AbstractTens{2,3}) -> NTuple{6}
    voigt_strain(ε::AbstractTens{2,3}) -> NTuple{6}

Canonical components of a symmetric second-order tensor in the **Voigt** order
`(11, 22, 33, 12, 13, 23)` of Abaqus-style interfaces.

The two differ by the engineering shears: `voigt_strain` doubles the
off-diagonal terms (``\\gamma_{12} = 2\\varepsilon_{12}``), `voigt_stress` does
not. Two functions rather than one with a flag, because a single `voigt` is the
classic way to lose a factor of two between the strain that goes in and the
stress that comes out. Voigt is not an isometry, unlike the Kelvin–Mandel form
of [`mandel`](@ref), and belongs only at the boundary with codes that demand it.
"""
voigt_stress(σ::AbstractTens{2, 3}) = _voigt(components_canon(σ), 1)

voigt_strain(ε::AbstractTens{2, 3}) = _voigt(components_canon(ε), 2)

@doc (@doc voigt_stress) voigt_strain

_voigt(A, s) = (A[1, 1], A[2, 2], A[3, 3], s * A[1, 2], s * A[1, 3], s * A[2, 3])

"""
    stress_from_voigt(v) -> Tens{2,3}
    strain_from_voigt(v) -> Tens{2,3}

Inverses of [`voigt_stress`](@ref) and [`voigt_strain`](@ref): the tensor, in the
canonical basis, whose Voigt components `(11, 22, 33, 12, 13, 23)` are `v`, the
engineering shears of a strain being halved.
"""
stress_from_voigt(v) = _from_voigt(v, 1)

strain_from_voigt(v) = _from_voigt(v, 2)

@doc (@doc stress_from_voigt) strain_from_voigt

function _from_voigt(v, s)
    length(v) == 6 || throw(DimensionMismatch("a Voigt vector has 6 components; got $(length(v))"))
    # Voigt position of the off-diagonal pair {i, j}: 12 → 4, 13 → 5, 23 → 6.
    return Tens(SymmetricTensor{2, 3}((i, j) -> i == j ? v[i] : v[i + j + 1] / s))
end
