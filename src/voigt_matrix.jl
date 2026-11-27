# ─────────────────────────────────────────────────────────────────────────────
#  voigt_matrix.jl — the 6×6 matrices of a stiffness and of a compliance in
#  the two-index (Voigt) notation.
#
#  The pairs are in the classical order 11, 22, 33, 23, 13, 12, that of the
#  Kelvin-Mandel notation of this package (`kelvin_mandel.jl`); only the
#  weights differ. With w = (1, 1, 1, √2, √2, √2) and M the Kelvin-Mandel
#  matrix,
#
#       C_IJ = M_IJ / (w_I w_J)      for a stiffness  (C₄₄ = C₂₃₂₃) ,
#       S_IJ = M_IJ · (w_I w_J)      for a compliance (S₄₄ = 4 S₂₃₂₃) ,
#
#  so that σ_I = C_IJ ε_J and ε_I = S_IJ σ_J hold with engineering shears
#  (ε₄ = 2ε₂₃) and plain stresses. The two matrices are inverses of each other
#  when the tensors are. Two functions rather than one with a flag: the
#  factors are those of the role, not of the tensor. The Voigt vectors of
#  order 2 of `tensors_bridge.jl` use another order of the shears, (12, 13,
#  23), that of the finite-element interfaces they serve.
# ─────────────────────────────────────────────────────────────────────────────

_voigt_weights(::Type{T}) where {T} = mandel_weights(Val(3), T)

# The Kelvin-Mandel matrix, in the canonical basis, of the minor-symmetric
# part of a tensor of order 4.
_canonical_mandel(t::AbstractTens{4, 3}) =
    SMatrix{6, 6}(_mandel_minor(canonical_components(t)))

"""
    voigt_stiffness(ℂ::AbstractTens{4,3}) → SMatrix{6,6}
    voigt_compliance(𝕊::AbstractTens{4,3}) → SMatrix{6,6}

The matrix of a stiffness, or of a compliance, in the two-index notation, with
the pairs in the order `(11, 22, 33, 23, 13, 12)` and canonical components:

```math
\\sigma_I = C_{IJ}\\,\\varepsilon_J, \\qquad \\varepsilon_I = S_{IJ}\\,\\sigma_J,
\\qquad \\varepsilon_4 = 2\\varepsilon_{23},\\ \\varepsilon_5 = 2\\varepsilon_{13},\\ \\varepsilon_6 = 2\\varepsilon_{12} ,
```

that is ``C_{IJ} = C_{ijkl}`` and ``S_{IJ} = S_{ijkl}`` times 1, 2 or 4 according
to the number of shear pairs among ``I`` and ``J``. They differ from the
Kelvin–Mandel matrix of [`mandel`](@ref) by these weights only, and unlike it
they are not the matrix of the tensor in an orthonormal basis: use them at the
boundary with a source written in this notation. The components read are those
of the minor-symmetric part of the tensor.

See also [`stiffness_from_voigt`](@ref), [`compliance_from_voigt`](@ref).
"""
function voigt_stiffness(ℂ::AbstractTens{4, 3})
    M = _canonical_mandel(ℂ)
    w = SVector(_voigt_weights(eltype(M)))
    return M ./ (w * transpose(w))
end

function voigt_compliance(𝕊::AbstractTens{4, 3})
    M = _canonical_mandel(𝕊)
    w = SVector(_voigt_weights(eltype(M)))
    return M .* (w * transpose(w))
end

@doc (@doc voigt_stiffness) voigt_compliance

"""
    stiffness_from_voigt(C::AbstractMatrix) → Tens{4,3}
    compliance_from_voigt(S::AbstractMatrix) → Tens{4,3}

Inverses of [`voigt_stiffness`](@ref) and [`voigt_compliance`](@ref): the
tensor, in the canonical basis, whose 6×6 matrix in the two-index notation
`(11, 22, 33, 23, 13, 12)` is given.
"""
stiffness_from_voigt(C::AbstractMatrix) = _from_voigt_matrix(C, (m, ww) -> m * ww)
compliance_from_voigt(S::AbstractMatrix) = _from_voigt_matrix(S, (m, ww) -> m / ww)

@doc (@doc stiffness_from_voigt) compliance_from_voigt

function _from_voigt_matrix(V::AbstractMatrix, weigh::F) where {F}
    size(V) == (6, 6) || throw(DimensionMismatch("a matrix in the two-index notation is 6×6; got $(size(V))"))
    T = sqrt2_eltype(eltype(V))
    w = SVector(_voigt_weights(T))
    return from_mandel(SMatrix{6, 6, T}(weigh.(V, w * transpose(w))))
end

export voigt_stiffness, voigt_compliance, stiffness_from_voigt, compliance_from_voigt
