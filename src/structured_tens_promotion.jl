##############################################################################
# Cross-type promotion & dispatch between TensISO, TensTI{4}, TensTI,     #
# TensOrtho.                                                                 #
#                                                                            #
# Rule: a binary operation on two structured tensors returns the             #
# highest-symmetry type compatible with both operands' references            #
# (axis / material frame).  When references are incompatible the operation   #
# falls back to the generic `Tens` route via `get_array`.                     #
#                                                                            #
# Symmetry lattice at 4th order (major-symmetric):                           #
#                                                                            #
#   TensISO{4}  ⊂  TensTI{4,N=5}  ⊂  TensOrtho    (with aligned axis)     #
#                                                                            #
# At 2nd order:                                                              #
#                                                                            #
#   TensISO{2}  ⊂  TensTI{2}                                                 #
#                                                                            #
# Products (double contraction) of two major-symmetric tensors are not       #
# always major-symmetric, so not every `⊡` lifts to the richer container.    #
# See implementation comments for each method for the mathematical           #
# justification.                                                             #
##############################################################################

# ──────────────────────────────────────────────────────────────────────────────
# Promotion helpers
# ──────────────────────────────────────────────────────────────────────────────

"""
    iso_to_ortho(A::TensISO{4,3,T}, frame::OrthonormalBasis{3}) → TensOrtho{T}

Convert an isotropic 4th-order tensor `α·𝕁 + β·𝕂` into a `TensOrtho` stored in
the given material frame. An isotropic tensor is orthotropic in *any* frame;
the orthotropic coefficients are

    C₁₁ = C₂₂ = C₃₃ = (α + 2β) / 3
    C₁₂ = C₁₃ = C₂₃ = (α − β)  / 3
    C₄₄ = C₅₅ = C₆₆ = β / 2

This is the obvious promotion for operations that combine `TensISO{4}` with a
`TensOrtho`, e.g. `TensISO + TensOrtho → TensOrtho`.

# Examples
```julia
julia> I4 = TensISO{3}(2.0, 3.0);   # α=2, β=3

julia> O = iso_to_ortho(I4, CanonicalBasis{3,Float64}());

julia> typeof(O) === TensOrtho{Float64}
true

julia> get_data(O)[1], get_data(O)[4], get_data(O)[7]
(2.6666666666666665, -0.3333333333333333, 1.5)
```

See also [`walpole_to_ortho`](@ref), [`fromISO`](@ref).
"""
function iso_to_ortho(A::TensISO{4, 3, T}, frame::OrthonormalBasis{3}) where {T}
    α, β = get_data(A)
    C_diag = (α + 2β) / 3          # C₁₁ = C₂₂ = C₃₃
    C_off = (α - β) / 3            # C₁₂ = C₁₃ = C₂₃
    C_sh = β / 2                   # C₄₄ = C₅₅ = C₆₆
    return TensOrtho(
        C_diag, C_diag, C_diag,
        C_off, C_off, C_off,
        C_sh, C_sh, C_sh,
        frame,
    )
end

"""
    _axis_on_frame_index(n, frame::OrthonormalBasis{3}; atol=1e-10) → Int

Return the index `k ∈ {1,2,3}` such that `n` is parallel to the `k`-th axis
of `frame` (i.e. `|n ⋅ eₖ| ≈ 1`), or `0` if `n` is not aligned with any axis.

Used to decide whether a `TensTI{4}` (TI with axis `n`) and a `TensOrtho`
(with a material frame) are compatible for symmetry-preserving arithmetic.

# Examples
```julia
julia> TensND._axis_on_frame_index((0.0, 0.0, 1.0), CanonicalBasis{3,Float64}())
3

julia> TensND._axis_on_frame_index((1.0, 0.0, 0.0), CanonicalBasis{3,Float64}())
1

julia> TensND._axis_on_frame_index((1.0, 1.0, 0.0) ./ √2, CanonicalBasis{3,Float64}())
0
```
"""
function _axis_on_frame_index(n, frame::OrthonormalBasis{3}; atol = 1.0e-10)
    nv = _extract_vec(n)
    E = vecbasis(frame, :cov)
    for k in 1:3
        dot_val = nv[1] * E[1, k] + nv[2] * E[2, k] + nv[3] * E[3, k]
        if abs(abs(dot_val) - one(dot_val)) < atol
            return k
        end
    end
    return 0
end

"""
    walpole_to_ortho(A::TensTI{4, T, 5}, frame::OrthonormalBasis{3}, axis_idx::Int) → TensOrtho{T}

Convert a major-symmetric `TensTI{4, T, 5}` into a `TensOrtho` stored in the
given material frame, assuming the Walpole axis `A.n` is aligned with axis
`axis_idx ∈ {1,2,3}` of the frame. A TI tensor is a special case of an
orthotropic tensor (with the 1–2 equivalence about the TI axis).

Mapping in the rotated material frame (TI axis = axis k):

    Ctrtr     = (ℓ₂ + ℓ₅) / 2    (transverse–transverse C in KM)
    Cmix      = (ℓ₂ − ℓ₅) / 2    (between the two transverse axes)
    Cax_tr    = ℓ₃ / √2          (axial ↔ transverse)
    Cax       = ℓ₁               (axial–axial)
    Cshear_ax = ℓ₆ / 2           (shear involving the axial axis)
    Cshear_tr = ℓ₅ / 2           (shear in the transverse plane)

The 9 orthotropic coefficients are then permuted according to `axis_idx`.

Restricted to `N=5` (major-symmetric Walpole): a non-major-symmetric TI
(`N=6`) has no orthotropic counterpart in `TensOrtho` (which carries only
9 major-symmetric constants).

See also [`iso_to_ortho`](@ref), [`_axis_on_frame_index`](@ref).
"""
function walpole_to_ortho(A::TensTI{4, T, 5}, frame::OrthonormalBasis{3}, axis_idx::Int) where {T}
    ℓ₁, ℓ₂, ℓ₃, _, ℓ₅, ℓ₆ = get_ℓ(A)
    sq2 = sqrt(T(2))
    Ctrtr = (ℓ₂ + ℓ₅) / 2
    Cmix = (ℓ₂ - ℓ₅) / 2
    Cax_tr = ℓ₃ / sq2
    Cax = ℓ₁
    Cshear_ax = ℓ₆ / 2
    Cshear_tr = ℓ₅ / 2

    if axis_idx == 3
        return TensOrtho(
            Ctrtr, Ctrtr, Cax,
            Cmix, Cax_tr, Cax_tr,
            Cshear_ax, Cshear_ax, Cshear_tr,
            frame,
        )
    elseif axis_idx == 1
        # TI axis = e₁; transverse plane = (e₂, e₃)
        return TensOrtho(
            Cax, Ctrtr, Ctrtr,
            Cax_tr, Cax_tr, Cmix,
            Cshear_tr, Cshear_ax, Cshear_ax,
            frame,
        )
    elseif axis_idx == 2
        # TI axis = e₂; transverse plane = (e₁, e₃)
        return TensOrtho(
            Ctrtr, Cax, Ctrtr,
            Cax_tr, Cmix, Cax_tr,
            Cshear_ax, Cshear_tr, Cshear_ax,
            frame,
        )
    else
        error("axis_idx must be 1, 2 or 3; got $axis_idx")
    end
end

"""
    _lift_walpole_N6(A::TensTI{4, T, 5}) → TensTI{4, T, 6}

Lift a major-symmetric `N=5` Walpole tensor to the general `N=6` form by
duplicating the shared coefficient `ℓ₃ = ℓ₄`.
"""
function _lift_walpole_N6(A::TensTI{4, T, 5}) where {T}
    ℓ₁, ℓ₂, ℓ₃, _, ℓ₅, ℓ₆ = get_ℓ(A)
    return TensTI{4, T, 6}((ℓ₁, ℓ₂, ℓ₃, ℓ₃, ℓ₅, ℓ₆), axis(A))
end

# ──────────────────────────────────────────────────────────────────────────────
# Arithmetic (+, −) — cross-type
# ──────────────────────────────────────────────────────────────────────────────
#
# Rule: the sum of two structured tensors of compatible references is
# structured (unlike products, there is no asymmetry issue).  Promote the
# lower-symmetry operand to the richer container, then delegate to the
# existing same-type +/- method.

# ── TensISO{4,3} ± TensOrtho ────────────────────────────────────────────────

for OP in (:+, :-)
    @eval @inline function Base.$OP(A::TensISO{4, 3}, B::TensOrtho)
        return $OP(iso_to_ortho(A, frame(B)), B)
    end
    @eval @inline function Base.$OP(A::TensOrtho, B::TensISO{4, 3})
        return $OP(A, iso_to_ortho(B, frame(A)))
    end
end

# ── TensTI{4,N=5} ± TensOrtho (aligned axis) ─────────────────────────────

for OP in (:+, :-)
    @eval function Base.$OP(A::TensTI{4, <:Any, 5}, B::TensOrtho)
        k = _axis_on_frame_index(axis(A), frame(B))
        k == 0 && throw(
            AssertionError(
                "TensTI{4,N=5} $($(string(OP))) TensOrtho requires the Walpole axis " *
                    "to be aligned with one of the frame axes",
            ),
        )
        return $OP(walpole_to_ortho(A, frame(B), k), B)
    end
    @eval function Base.$OP(A::TensOrtho, B::TensTI{4, <:Any, 5})
        k = _axis_on_frame_index(axis(B), frame(A))
        k == 0 && throw(
            AssertionError(
                "TensOrtho $($(string(OP))) TensTI{4,N=5} requires the Walpole axis " *
                    "to be aligned with one of the frame axes",
            ),
        )
        return $OP(A, walpole_to_ortho(B, frame(A), k))
    end
end

# ── TensTI{4,N=5} ± TensTI{4,N=6} (same axis) ─────────────────────────

for OP in (:+, :-)
    @eval function Base.$OP(A::TensTI{4, <:Any, 5}, B::TensTI{4, <:Any, 6})
        return $OP(_lift_walpole_N6(A), B)
    end
    @eval function Base.$OP(A::TensTI{4, <:Any, 6}, B::TensTI{4, <:Any, 5})
        return $OP(A, _lift_walpole_N6(B))
    end
end

# ── TensISO{4,3} ± TensTI{4,N} (any axis) ──────────────────────────────────
# An iso 4-tensor is a special case of TI(axis) for any axis. Promote the
# iso operand to TI using `fromISO` (lift to N=6 if the TI partner is N=6)
# so the sum stays in the TI Walpole basis instead of falling through to
# `TensCanonical` via the unstructured `+ ::AbstractTens` route.

for OP in (:+, :-)
    @eval function Base.$OP(A::TensISO{4, 3}, B::TensTI{4, <:Any, 5})
        return $OP(fromISO(A, axis(B)), B)
    end
    @eval function Base.$OP(A::TensTI{4, <:Any, 5}, B::TensISO{4, 3})
        return $OP(A, fromISO(B, axis(A)))
    end
    @eval function Base.$OP(A::TensISO{4, 3}, B::TensTI{4, <:Any, 6})
        return $OP(_lift_walpole_N6(fromISO(A, axis(B))), B)
    end
    @eval function Base.$OP(A::TensTI{4, <:Any, 6}, B::TensISO{4, 3})
        return $OP(A, _lift_walpole_N6(fromISO(B, axis(A))))
    end
end

# ── TensISO{2,3} ± TensTI{2} ────────────────────────────────────────────────
# TensISO{2,3}(λ) ≡ λ·𝟏 ≡ TensTI{2}(λ, λ, n) for any axis n.

for OP in (:+, :-)
    @eval function Base.$OP(A::TensISO{2, 3}, B::TensTI{2, <:Any, 2})
        λ = get_data(A)[1]
        a, b = get_data(B)
        return TensTI{2}($OP(λ, a), $OP(λ, b), axis(B))
    end
    @eval function Base.$OP(A::TensTI{2, <:Any, 2}, B::TensISO{2, 3})
        a, b = get_data(A)
        λ = get_data(B)[1]
        return TensTI{2}($OP(a, λ), $OP(b, λ), axis(A))
    end
end

# ──────────────────────────────────────────────────────────────────────────────
# Cross-order double contraction (4th order ⊡ 2nd order → 2nd order)
# ──────────────────────────────────────────────────────────────────────────────
#
# Mathematical foundation (all derivations in the Walpole basis with
# `t = a·nT + b·nₙ`):
#
#   W₁ ⊡ t = b·nₙ         W₂ ⊡ t = a·nT
#   W₃ ⊡ t = √2 a · nₙ    W₄ ⊡ t = (b/√2)·nT
#   W₅ ⊡ t = 0            W₆ ⊡ t = 0
#
# And (left contraction) t ⊡ Wᵢ is the same as Wᵢ ⊡ t for W₁,W₂,W₅,W₆ but
# swaps W₃ and W₄ (non-major-symmetric case).

# ── TensISO{4,3} ⊡ TensTI{2} / TensTI{2} ⊡ TensISO{4,3} ─────────────────────
#
# (α𝕁 + β𝕂) ⊡ (a·nT + b·nₙ) = β·t + ((α−β)/3)·(tr t)·𝟏
# where tr t = 2a + b.  Expanding in the (nT, nₙ) basis:
#
#   new_a = (α − β)·(2a + b)/3 + β·a
#   new_b = (α − β)·(2a + b)/3 + β·b

function Tensors.dcontract(A::TensISO{4, 3}, B::TensTI{2, <:Any, 2})
    α, β = get_data(A)
    a, b = get_data(B)
    sph = (α - β) * (2a + b) / 3
    return TensTI{2}(sph + β * a, sph + β * b, axis(B))
end

function Tensors.dcontract(A::TensTI{2, <:Any, 2}, B::TensISO{4, 3})
    # (α𝕁 + β𝕂) is major-symmetric, so A⊡B = B⊡A as 2nd-order output.
    return Tensors.dcontract(B, A)
end

# ── TensTI{4} ⊡ TensTI{2} / TensTI{2} ⊡ TensTI{4} (same axis) ──────────

function Tensors.dcontract(A::TensTI{4}, B::TensTI{2, <:Any, 2})
    @assert axis(A) == axis(B) "dcontract(TensTI{4}, TensTI{2}) requires the same axis"
    T = promote_type(eltype(A), eltype(B))
    ℓ₁, ℓ₂, ℓ₃, ℓ₄, _, _ = get_ℓ(A)
    a, b = get_data(B)
    sq2 = sqrt(T(2))
    new_a = ℓ₂ * a + ℓ₄ * b / sq2
    new_b = ℓ₁ * b + sq2 * ℓ₃ * a
    return TensTI{2}(new_a, new_b, axis(A))
end

function Tensors.dcontract(A::TensTI{2, <:Any, 2}, B::TensTI{4})
    @assert axis(A) == axis(B) "dcontract(TensTI{2}, TensTI{4}) requires the same axis"
    T = promote_type(eltype(A), eltype(B))
    ℓ₁, ℓ₂, ℓ₃, ℓ₄, _, _ = get_ℓ(B)
    a, b = get_data(A)
    sq2 = sqrt(T(2))
    # For `t ⊡ W`, the contraction is on the FIRST pair of W, swapping the
    # roles of ℓ₃ (W₃) and ℓ₄ (W₄) in the axial/transverse mix.
    new_a = ℓ₂ * a + ℓ₃ * b / sq2
    new_b = ℓ₁ * b + sq2 * ℓ₄ * a
    return TensTI{2}(new_a, new_b, axis(A))
end

# ──────────────────────────────────────────────────────────────────────────────
# Single contraction (2nd order · 2nd order)
# ──────────────────────────────────────────────────────────────────────────────
#
# For TI 2nd-order tensors `a·nT + b·nₙ`, single contraction is diagonal in
# the (nT, nₙ) decomposition because nT and nₙ are orthogonal projectors:
#
#   (a₁·nT + b₁·nₙ) · (a₂·nT + b₂·nₙ)
#   = a₁a₂·nT·nT + b₁b₂·nₙ·nₙ
#   = a₁a₂·nT + b₁b₂·nₙ

function LinearAlgebra.dot(A::TensTI{2, <:Any, 2}, B::TensTI{2, <:Any, 2})
    @assert axis(A) == axis(B) "dot(TensTI{2}, TensTI{2}) requires the same axis"
    a₁, b₁ = get_data(A)
    a₂, b₂ = get_data(B)
    return TensTI{2}(a₁ * a₂, b₁ * b₂, axis(A))
end

# ── TensTI{2} · TensISO{2,3} ────────────────────────────────────────────────
# TensISO{2,3}(λ) ≡ λ·𝟏, so A·B = λ·A (TI output preserved).

function LinearAlgebra.dot(A::TensTI{2, <:Any, 2}, B::TensISO{2, 3})
    a, b = get_data(A)
    λ = get_data(B)[1]
    return TensTI{2}(a * λ, b * λ, axis(A))
end

function LinearAlgebra.dot(A::TensISO{2, 3}, B::TensTI{2, <:Any, 2})
    λ = get_data(A)[1]
    a, b = get_data(B)
    return TensTI{2}(λ * a, λ * b, axis(B))
end

# ──────────────────────────────────────────────────────────────────────────────
# Exports
# ──────────────────────────────────────────────────────────────────────────────

export iso_to_ortho, walpole_to_ortho
