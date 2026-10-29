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
#   TensISO{4}  ⊂  TensCubic     ⊂  TensOrtho    (on the same frame)      #
#                                                                            #
# TensCubic and TensTI are INCOMPARABLE: their intersection is the isotropic #
# class and neither contains the other; a TI axis along a cube axis leaves   #
# them a tetragonal, hence orthotropic, sum; otherwise the generic route.    #
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
    iso_to_ortho(A::TensISO{4,3,T}, frame::OrthonormalBasis{3}) → TensOrtho{4,T,9}

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

julia> O isa TensOrtho{4, Float64, 9}
true

julia> get_data(O)[1], get_data(O)[4], get_data(O)[7]
(2.6666666666666665, -0.3333333333333333, 1.5)
```

See also [`walpole_to_ortho`](@ref), [`iso_to_ti`](@ref).
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
    walpole_to_ortho(A::TensTI{4, T, 5 or 6}, frame::OrthonormalBasis{3}, axis_idx::Int) → TensOrtho{4, T, 9 or 12}

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

A general `TensTI{4, T, 6}` (`ℓ₃ ≠ ℓ₄`) maps the same way onto the
twelve-constant form: `ℓ₃/√2` sits in the axial row of the normal block (axial
stress, transverse strain), `ℓ₄/√2` in the axial column.

See also [`iso_to_ortho`](@ref), `_axis_on_frame_index`.
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

function walpole_to_ortho(A::TensTI{4, T, 6}, frame::OrthonormalBasis{3}, axis_idx::Int) where {T}
    ℓ₁, ℓ₂, ℓ₃, ℓ₄, ℓ₅, ℓ₆ = get_ℓ(A)
    sq2 = sqrt(T(2))
    Ctrtr = (ℓ₂ + ℓ₅) / 2
    Cmix = (ℓ₂ - ℓ₅) / 2
    Cax_tr = ℓ₃ / sq2   # axial row
    Ctr_ax = ℓ₄ / sq2   # axial column
    Cshear_ax = ℓ₆ / 2
    Cshear_tr = ℓ₅ / 2
    # (C₁₁,C₂₂,C₃₃,C₁₂,C₁₃,C₂₃,C₄₄,C₅₅,C₆₆,C₂₁,C₃₁,C₃₂), Cₘₙ in row m, column n
    if axis_idx == 3
        return TensOrtho(
            Ctrtr, Ctrtr, ℓ₁, Cmix, Ctr_ax, Ctr_ax,
            Cshear_ax, Cshear_ax, Cshear_tr, Cmix, Cax_tr, Cax_tr, frame,
        )
    elseif axis_idx == 1
        return TensOrtho(
            ℓ₁, Ctrtr, Ctrtr, Cax_tr, Cax_tr, Cmix,
            Cshear_tr, Cshear_ax, Cshear_ax, Ctr_ax, Ctr_ax, Cmix, frame,
        )
    elseif axis_idx == 2
        return TensOrtho(
            Ctrtr, ℓ₁, Ctrtr, Ctr_ax, Cmix, Cax_tr,
            Cshear_ax, Cshear_tr, Cshear_ax, Cax_tr, Cmix, Ctr_ax, frame,
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

"""
    iso_to_cubic(A::TensISO{4,3,T}, frame::OrthonormalBasis{3}) → TensCubic{T}

Convert an isotropic 4th-order tensor `α𝕁 + β𝕂` into a `TensCubic` on `frame`.

Exact and frame-independent: `𝕂 = 𝔼 + 𝕋`, so the cubic coefficients are simply
`(α, β, β)` — and they are the same whatever frame is chosen, which is the
formal statement that an isotropic tensor is cubic about *every* cube.

Used by the mixed arithmetic below so that `TensISO + TensCubic` stays a
`TensCubic` instead of collapsing to the unstructured route.

# Examples
```julia
julia> C = iso_to_cubic(TensISO{3}(18.0, 6.0), CanonicalBasis{3,Float64}());

julia> cubic_anisotropy(C)          # isotropic, so exactly zero
0.0
```
"""
function iso_to_cubic(A::TensISO{4, 3, T}, frame::OrthonormalBasis{3}) where {T}
    α, β = get_data(A)
    return TensCubic{T}((α, β, β), frame)
end

"""
    cubic_to_ortho(A::TensCubic{T}) → TensOrtho{4,T,9}

Widen a cubic tensor to the orthotropic class on its own cube frame:
`C₁₁ = C₂₂ = C₃₃`, `C₁₂ = C₁₃ = C₂₃`, `C₄₄ = C₅₅ = C₆₆`.

Exact — no projection, no loss. It is the embedding that makes
`TensCubic ± TensOrtho` well defined on a shared frame.
"""
cubic_to_ortho(A::TensCubic) = _cubic_as_ortho(A, frame(A))

# A cubic tensor written on a frame of the same cube (its own, or one of its
# axes permuted and signed): the same coefficients, the octahedral group
# leaving it unchanged.
function _cubic_as_ortho(A::TensCubic, ℱ)
    C₁₁, C₁₂, C₄₄ = arg_cubic(A)
    return TensOrtho(C₁₁, C₁₁, C₁₁, C₁₂, C₁₂, C₁₂, C₄₄, C₄₄, C₄₄, ℱ)
end

# ── TensISO{4,3} ± TensCubic ────────────────────────────────────────────────
# An isotropic tensor is cubic about every cube, so no frame check is needed:
# the iso operand is lifted onto whatever frame the cubic one carries.

for OP in (:+, :-)
    @eval @inline function Base.$OP(A::TensISO{4, 3}, B::TensCubic)
        return _same_frame_binary($OP, iso_to_cubic(A, frame(B)), B)
    end
    @eval @inline function Base.$OP(A::TensCubic, B::TensISO{4, 3})
        return _same_frame_binary($OP, A, iso_to_cubic(B, frame(A)))
    end
end

# ── TensCubic ± TensOrtho (same frame) ──────────────────────────────────────
# Cubic is orthotropic with the three families equal, so on a frame of the same
# cube (the orthotropic one, up to a signed permutation of the axes) the sum is
# orthotropic. Otherwise it is generally fully anisotropic: the generic route,
# as the rule at the top of this file says.

for OP in (:+, :-)
    @eval @inline function Base.$OP(A::TensCubic, B::TensOrtho{4})
        frame(A) == frame(B) || _permuted_frames(frame(A), frame(B)) || return _structured_sum($OP, A, B)
        return _same_frame_binary($OP, _cubic_as_ortho(A, frame(B)), B)
    end
    @eval @inline function Base.$OP(A::TensOrtho{4}, B::TensCubic)
        frame(A) == frame(B) || _permuted_frames(frame(A), frame(B)) || return _structured_sum($OP, A, B)
        return _same_frame_binary($OP, A, _cubic_as_ortho(B, frame(A)))
    end
end

# ── Sums on a reference known to be shared ──────────────────────────────────
# Coefficient by coefficient, a concrete type: `+` checks the reference at run
# time and returns a union, which the conversions below would carry for nothing.

for (order, N) in ((2, 3), (4, 9), (4, 12))
    @eval @inline _same_frame_binary(op::F, A::TensOrtho{$order, <:Any, $N}, B::TensOrtho{$order, <:Any, $N}) where {F} =
        _rebuild(A, map(op, get_data(A), get_data(B)))
end
@inline _same_frame_binary(op::F, A::TensOrtho{4}, B::TensOrtho{4}) where {F} =
    _same_frame_binary(op, _ortho_general(A), _ortho_general(B))
@inline _same_frame_binary(op::F, A::TensCubic, B::TensCubic) where {F} = _rebuild(A, map(op, get_data(A), get_data(B)))

# ── The generic route between structured fourth-order tensors ──────────────
# When two structured tensors share no class, their sum is the general tensor of
# canonical Kelvin–Mandel matrix the sum of theirs, each formed in closed form:
# the congruence of its sparse material matrix for a tensor with a frame, its
# entries for a TI one. Measured 60 times faster than summing their 81
# canonical components (163 ns against 10 µs for two cubic tensors), and the
# result is in canonical components, which is what its users read next
# (working in the frame of one operand saved 14 ns on the sum and cost
# 3.4 times more once its canonical components were asked for).

const _KMStructured = Union{TensOrtho{4}, TensCubic, TensTI{4}}
@noinline _structured_sum(op::F, A::_KMStructured, B::_KMStructured) where {F} = from_mandel(op(_km_canonical(A), _km_canonical(B)))
@inline _km_canonical(t::Union{TensOrtho{4}, TensCubic}) = _km_in_frame(t)
@inline _km_canonical(t::TensTI{4}) = _km_by_entries(t)

# At order 2, the tensor of the canonical components: symmetric for the D∞h and
# orthotropic shapes, stored so from the start, a concrete type the caller
# splits without boxing; the shapes of `N = 3` have a skew part, and the storage
# of their sum is decided by its value.
const _SymmetricOrder2 = Union{TensTI{2, <:Any, 2}, TensOrtho{2}}
@noinline function _structured_sum(op::F, A::_SymmetricOrder2, B::_SymmetricOrder2) where {F}
    a, b = components_canon(A), components_canon(B)
    return Tens(Tensors.SymmetricTensor{2, 3}((i, j) -> op(a[i, j], b[i, j])))
end
@noinline _structured_sum(op::F, A::Union{TensTI{2}, TensOrtho{2}}, B::Union{TensTI{2}, TensOrtho{2}}) where {F} =
    Tens(Tensors.Tensor{2, 3}(broadcast(op, components_canon(A), components_canon(B))))

# ── TensCubic ± TensTI{4,N=5|6} (axis along a cube axis) ────────────────────
# Cubic and TI are incomparable classes, but a TI axis along an axis of the cube
# leaves them the tetragonal symmetry of that axis, orthotropic on the frame of
# the cube. Otherwise the generic route.

for OP in (:+, :-), NW in (5, 6)
    @eval function Base.$OP(A::TensCubic, B::TensTI{4, <:Any, $NW})
        k = _axis_on_frame_index(axis(B), frame(A))
        k == 0 && return _structured_sum($OP, A, B)
        return _same_frame_binary($OP, cubic_to_ortho(A), walpole_to_ortho(B, frame(A), k))
    end
    @eval function Base.$OP(A::TensTI{4, <:Any, $NW}, B::TensCubic)
        k = _axis_on_frame_index(axis(A), frame(B))
        k == 0 && return _structured_sum($OP, A, B)
        return _same_frame_binary($OP, walpole_to_ortho(A, frame(B), k), cubic_to_ortho(B))
    end
end

# ── TensISO{4,3} ± TensOrtho ────────────────────────────────────────────────

for OP in (:+, :-)
    @eval @inline function Base.$OP(A::TensISO{4, 3}, B::TensOrtho{4})
        return _same_frame_binary($OP, iso_to_ortho(A, frame(B)), B)
    end
    @eval @inline function Base.$OP(A::TensOrtho{4}, B::TensISO{4, 3})
        return _same_frame_binary($OP, A, iso_to_ortho(B, frame(A)))
    end
end

# ── TensTI{4,N=5|6} ± TensOrtho (aligned axis) ───────────────────────────
# A six-parameter TI tensor lands on the twelve-constant form.

for OP in (:+, :-), NW in (5, 6)
    @eval function Base.$OP(A::TensTI{4, <:Any, $NW}, B::TensOrtho{4})
        k = _axis_on_frame_index(axis(A), frame(B))
        k == 0 && return _structured_sum($OP, A, B)
        return _same_frame_binary($OP, walpole_to_ortho(A, frame(B), k), B)
    end
    @eval function Base.$OP(A::TensOrtho{4}, B::TensTI{4, <:Any, $NW})
        k = _axis_on_frame_index(axis(B), frame(A))
        k == 0 && return _structured_sum($OP, A, B)
        return _same_frame_binary($OP, A, walpole_to_ortho(B, frame(A), k))
    end
end

# ── TensTI{4,N=5} ± TensTI{4,N=6} ─────────────────────────────────────────
# The five-constant operand lifted to six, then the TI sum, which checks the axes.

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
# iso operand to TI using `iso_to_ti` (lift to N=6 if the TI partner is N=6)
# so the sum stays in the TI Walpole basis instead of falling through to
# `TensCanonical` via the unstructured `+ ::AbstractTens` route.

for OP in (:+, :-)
    @eval function Base.$OP(A::TensISO{4, 3}, B::TensTI{4, <:Any, 5})
        return _coaxial_binary($OP, iso_to_ti(A, axis(B)), B)
    end
    @eval function Base.$OP(A::TensTI{4, <:Any, 5}, B::TensISO{4, 3})
        return _coaxial_binary($OP, A, iso_to_ti(B, axis(A)))
    end
    @eval function Base.$OP(A::TensISO{4, 3}, B::TensTI{4, <:Any, 6})
        return _coaxial_binary($OP, _lift_walpole_N6(iso_to_ti(A, axis(B))), B)
    end
    @eval function Base.$OP(A::TensTI{4, <:Any, 6}, B::TensISO{4, 3})
        return _coaxial_binary($OP, A, _lift_walpole_N6(iso_to_ti(B, axis(A))))
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

# ── TensISO{2,3} ⊡ TensTI{4} / TensTI{4} ⊡ TensISO{2,3} ─────────────────────
#
# `TensISO{2,3}(λ) ≡ λ·𝟏 ≡ TensTI{2}(λ, λ, n)` for any axis, so the whole rule
# is a conversion onto the axis of the 4th-order operand followed by the
# TensTI{4} ⊡ TensTI{2} rule below. Without these two methods the contraction
# of the second-order identity with a transversely isotropic 4th-order tensor
# — `𝟏 : 𝔸`, the trace of a concentration tensor, i.e. exactly what a Levin or
# a Biot post-processing forms — fell through to the generic 81-component
# array path and came back as an unstructured `Tens`.

function Tensors.dcontract(A::TensISO{2, 3}, B::TensTI{4})
    λ = get_data(A)[1]
    return Tensors.dcontract(TensTI{2}(λ, λ, axis(B)), B)
end

function Tensors.dcontract(A::TensTI{4}, B::TensISO{2, 3})
    λ = get_data(B)[1]
    return Tensors.dcontract(A, TensTI{2}(λ, λ, axis(A)))
end

# ── TensTI{4} ⊡ TensTI{2} / TensTI{2} ⊡ TensTI{4} (same axis) ──────────

function Tensors.dcontract(A::TensTI{4}, B::TensTI{2, <:Any, 2})
    axis(A) == axis(B) ||
        return Tensors.dcontract(generic_tens(A), generic_tens(B))
    T = promote_type(eltype(A), eltype(B))
    ℓ₁, ℓ₂, ℓ₃, ℓ₄, _, _ = get_ℓ(A)
    a, b = get_data(B)
    sq2 = sqrt(T(2))
    new_a = ℓ₂ * a + ℓ₄ * b / sq2
    new_b = ℓ₁ * b + sq2 * ℓ₃ * a
    return TensTI{2}(new_a, new_b, axis(A))
end

function Tensors.dcontract(A::TensTI{2, <:Any, 2}, B::TensTI{4})
    axis(A) == axis(B) ||
        return Tensors.dcontract(generic_tens(A), generic_tens(B))
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
    axis(A) == axis(B) ||
        return LinearAlgebra.dot(generic_tens(A), generic_tens(B))
    a₁, b₁ = get_data(A)
    a₂, b₂ = get_data(B)
    return TensTI{2}(a₁ * a₂, b₁ * b₂, axis(A))
end

# ── TensTI{2} · vector ──────────────────────────────────────────────────────
# A vector given by its canonical components (a `Vector`, a `Vec`), in closed
# form: `a (v − (n⋅v) n) + b (n⋅v) n + c n × v`. Without it the call fell
# through to `LinearAlgebra.dot` of two arrays, an inner product of a 3×3
# with a 3-vector.

function LinearAlgebra.dot(A::TensTI{2}, v::AbstractVector)
    length(v) == 3 || throw(DimensionMismatch("a vector of 3 components is contracted with a TensTI{2}; got $(length(v))"))
    a, b, c = _ti2_abc(A)
    n = axis(A)
    nv = n[1] * v[1] + n[2] * v[2] + n[3] * v[3]
    nxv = (n[2] * v[3] - n[3] * v[2], n[3] * v[1] - n[1] * v[3], n[1] * v[2] - n[2] * v[1])
    return Vec{3}(ntuple(i -> a * (v[i] - nv * n[i]) + b * nv * n[i] + c * nxv[i], 3))
end
# A TensND vector is an `AbstractVector` too: its canonical components, and a
# tensor back, as the generic contraction gives.
LinearAlgebra.dot(A::TensTI{2}, v::AbstractTens{1, 3}) = Tens(A ⋅ Vec{3}(Tuple(components_canon(v))))
_ti2_abc(A::TensTI{2, T, 2}) where {T} = (get_data(A)..., zero(T))
_ti2_abc(A::TensTI{2, T, 3}) where {T} = get_data(A)

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

# ── TensISO{4,3} ⊡ TensCubic / TensCubic ⊡ TensISO{4,3} ─────────────────────
# Both operands are then cubic on the same cube, so the componentwise rule of
# `dcontract(::TensCubic, ::TensCubic)` applies and the product stays in the
# class -- which the dense route would not have preserved.

@inline Tensors.dcontract(A::TensISO{4, 3}, B::TensCubic) =
    Tensors.dcontract(iso_to_cubic(A, frame(B)), B)
@inline Tensors.dcontract(A::TensCubic, B::TensISO{4, 3}) =
    Tensors.dcontract(A, iso_to_cubic(B, frame(A)))

# ── TensCubic ⊡ TensOrtho / TensOrtho ⊡ TensCubic ───────────────────────────
# Widen to orthotropy and let its block route take over: on a shared frame the
# result is the twelve-constant orthotropic form of
# `dcontract(::TensOrtho, ::TensOrtho)`, on different frames a generic tensor.

@inline Tensors.dcontract(A::TensCubic, B::TensOrtho{4}) =
    frame(A) == frame(B) ? Tensors.dcontract(cubic_to_ortho(A), B) :
    Tensors.dcontract(generic_tens(A), generic_tens(B))
@inline Tensors.dcontract(A::TensOrtho{4}, B::TensCubic) =
    frame(A) == frame(B) ? Tensors.dcontract(A, cubic_to_ortho(B)) :
    Tensors.dcontract(generic_tens(A), generic_tens(B))

# ── TensISO{4,3} ⊡ TensOrtho / TensOrtho ⊡ TensISO{4,3} ─────────────────────
# An isotropic tensor is orthotropic in *every* frame, so promoting it onto the
# other operand's material frame is exact and lets the closed form of
# `dcontract(::TensOrtho, ::TensOrtho)` take over instead of the dense route:
# the product is orthotropic, in the twelve-constant form.
# Mirrors the ± promotion above.

@inline Tensors.dcontract(A::TensISO{4, 3}, B::TensOrtho{4}) =
    Tensors.dcontract(iso_to_ortho(A, frame(B)), B)
@inline Tensors.dcontract(A::TensOrtho{4}, B::TensISO{4, 3}) =
    Tensors.dcontract(A, iso_to_ortho(B, frame(A)))

# ── TensTI{4,N=5|6} ⊡ TensOrtho / TensOrtho ⊡ TensTI{4,N=5|6} ───────────────
# A transversely isotropic tensor whose axis is one of the frame axes is
# orthotropic on that frame: the block product applies. Any other axis gives a
# generally anisotropic product, by the generic route.

for NW in (5, 6)
    @eval function Tensors.dcontract(A::TensTI{4, <:Any, $NW}, B::TensOrtho{4})
        k = _axis_on_frame_index(axis(A), frame(B))
        k == 0 && return Tensors.dcontract(generic_tens(A), generic_tens(B))
        return Tensors.dcontract(walpole_to_ortho(A, frame(B), k), B)
    end
    @eval function Tensors.dcontract(A::TensOrtho{4}, B::TensTI{4, <:Any, $NW})
        k = _axis_on_frame_index(axis(B), frame(A))
        k == 0 && return Tensors.dcontract(generic_tens(A), generic_tens(B))
        return Tensors.dcontract(A, walpole_to_ortho(B, frame(A), k))
    end
end

# ──────────────────────────────────────────────────────────────────────────────
# Order 2: TensOrtho{2} with TensISO{2,3}, TensTI{2,T,2} and the order-4 classes
# ──────────────────────────────────────────────────────────────────────────────

"""
    iso_to_ortho(A::TensISO{2,3,T}, frame::OrthonormalBasis{3}) → TensOrtho{2,T,3}

`λ𝟏 = λP₁ + λP₂ + λP₃` on any frame: the order-2 counterpart of the order-4
promotion.
"""
function iso_to_ortho(A::TensISO{2, 3}, frame::OrthonormalBasis{3})
    λ = get_data(A)[1]
    return TensOrtho{2}(λ, λ, λ, frame)
end

"""
    walpole_to_ortho(A::TensTI{2, T, 2}, frame::OrthonormalBasis{3}, axis_idx::Int) → TensOrtho{2}

`a nT + b nₙ` on a frame whose axis `axis_idx` is the TI axis: `b` on that axis,
`a` on the two others.
"""
function walpole_to_ortho(A::TensTI{2, T, 2}, frame::OrthonormalBasis{3}, axis_idx::Int) where {T}
    a, b = get_data(A)
    1 ≤ axis_idx ≤ 3 || error("axis_idx must be 1, 2 or 3; got $axis_idx")
    return TensOrtho{2}(ntuple(m -> m == axis_idx ? b : a, Val(3))..., frame)
end

# ± : an isotropic operand on the frame of the other; a TI operand when its
# axis is one of the frame axes, the generic route otherwise (as at order 4).
for OP in (:+, :-)
    @eval @inline Base.$OP(A::TensISO{2, 3}, B::TensOrtho{2}) = _same_frame_binary($OP, iso_to_ortho(A, frame(B)), B)
    @eval @inline Base.$OP(A::TensOrtho{2}, B::TensISO{2, 3}) = _same_frame_binary($OP, A, iso_to_ortho(B, frame(A)))
    @eval function Base.$OP(A::TensTI{2, <:Any, 2}, B::TensOrtho{2})
        k = _axis_on_frame_index(axis(A), frame(B))
        k == 0 && return _structured_sum($OP, A, B)
        return _same_frame_binary($OP, walpole_to_ortho(A, frame(B), k), B)
    end
    @eval function Base.$OP(A::TensOrtho{2}, B::TensTI{2, <:Any, 2})
        k = _axis_on_frame_index(axis(B), frame(A))
        k == 0 && return _structured_sum($OP, A, B)
        return _same_frame_binary($OP, A, walpole_to_ortho(B, frame(A), k))
    end
end

# ── Single contraction ───────────────────────────────────────────────────────

"""
    dot(A::TensOrtho{2}, B::TensOrtho{2}) → TensOrtho{2}

`(Σ aₘPₘ)·(Σ bₘPₘ) = Σ aₘbₘPₘ` on a shared frame; the generic product otherwise.
"""
function LinearAlgebra.dot(A::TensOrtho{2}, B::TensOrtho{2})
    frame(A) == frame(B) || return dot(generic_tens(A), generic_tens(B))
    return _rebuild(A, map(*, get_data(A), get_data(B)))
end
function LinearAlgebra.dot(A::TensTI{2, <:Any, 2}, B::TensOrtho{2})
    k = _axis_on_frame_index(axis(A), frame(B))
    k == 0 && return dot(generic_tens(A), generic_tens(B))
    return dot(walpole_to_ortho(A, frame(B), k), B)
end
function LinearAlgebra.dot(A::TensOrtho{2}, B::TensTI{2, <:Any, 2})
    k = _axis_on_frame_index(axis(B), frame(A))
    k == 0 && return dot(generic_tens(A), generic_tens(B))
    return dot(A, walpole_to_ortho(B, frame(A), k))
end
LinearAlgebra.dot(A::TensOrtho{2}, B::TensISO{2, 3}) = _rebuild(A, get_data(A) .* get_data(B)[1])
LinearAlgebra.dot(A::TensISO{2, 3}, B::TensOrtho{2}) = _rebuild(B, get_data(A)[1] .* get_data(B))

"""
    dot(A::TensOrtho{2}, v) → Vec{3}

`𝐀·v = Σ aₘ (eₘ·v) eₘ`, from the canonical components of `v`.
"""
function LinearAlgebra.dot(A::TensOrtho{2}, v::AbstractVector)
    length(v) == 3 || throw(DimensionMismatch("a vector of 3 components is contracted with a TensOrtho{2}; got $(length(v))"))
    E = SMatrix{3, 3}(vecbasis(frame(A), :cov))
    return Vec{3}(Tuple(E * (SVector(get_data(A)) .* (transpose(E) * SVector{3}(v[1], v[2], v[3])))))
end
LinearAlgebra.dot(A::TensOrtho{2}, v::AbstractTens{1, 3}) = Tens(A ⋅ Vec{3}(Tuple(components_canon(v))))

# ── Double contraction with an order-4 orthotropic tensor ────────────────────
# On a shared frame `ℂ : 𝐀 = Σₘ (Σₙ Cₘₙ aₙ) Pₘ` and `𝐀 : ℂ = Σₙ (Σₘ aₘ Cₘₙ) Pₙ`:
# the normal block acts on the coefficients, the shear terms vanish on a
# diagonal tensor. The order-4 operand may be any class that is orthotropic on
# that frame.

function Tensors.dcontract(A::TensOrtho{4}, B::TensOrtho{2})
    frame(A) == frame(B) || return Tensors.dcontract(generic_tens(A), generic_tens(B))
    return _rebuild(B, Tuple(_ortho_normal_block(A) * SVector(get_data(B))))
end
function Tensors.dcontract(A::TensOrtho{2}, B::TensOrtho{4})
    frame(A) == frame(B) || return Tensors.dcontract(generic_tens(A), generic_tens(B))
    return _rebuild(A, Tuple(transpose(_ortho_normal_block(B)) * SVector(get_data(A))))
end
@inline Tensors.dcontract(A::TensISO{4, 3}, B::TensOrtho{2}) = Tensors.dcontract(iso_to_ortho(A, frame(B)), B)
@inline Tensors.dcontract(A::TensOrtho{2}, B::TensISO{4, 3}) = Tensors.dcontract(A, iso_to_ortho(B, frame(A)))
@inline Tensors.dcontract(A::TensISO{2, 3}, B::TensOrtho{4}) = Tensors.dcontract(iso_to_ortho(A, frame(B)), B)
@inline Tensors.dcontract(A::TensOrtho{4}, B::TensISO{2, 3}) = Tensors.dcontract(A, iso_to_ortho(B, frame(A)))
@inline Tensors.dcontract(A::TensCubic, B::TensOrtho{2}) =
    frame(A) == frame(B) ? Tensors.dcontract(cubic_to_ortho(A), B) : Tensors.dcontract(generic_tens(A), generic_tens(B))
@inline Tensors.dcontract(A::TensOrtho{2}, B::TensCubic) =
    frame(A) == frame(B) ? Tensors.dcontract(A, cubic_to_ortho(B)) : Tensors.dcontract(generic_tens(A), generic_tens(B))
for NW in (5, 6)
    @eval function Tensors.dcontract(A::TensTI{4, <:Any, $NW}, B::TensOrtho{2})
        k = _axis_on_frame_index(axis(A), frame(B))
        k == 0 && return Tensors.dcontract(generic_tens(A), generic_tens(B))
        return Tensors.dcontract(walpole_to_ortho(A, frame(B), k), B)
    end
    @eval function Tensors.dcontract(A::TensOrtho{2}, B::TensTI{4, <:Any, $NW})
        k = _axis_on_frame_index(axis(B), frame(A))
        k == 0 && return Tensors.dcontract(generic_tens(A), generic_tens(B))
        return Tensors.dcontract(A, walpole_to_ortho(B, frame(A), k))
    end
end

# ── Tensor product ───────────────────────────────────────────────────────────
# `𝐀 ⊗ 𝐁 = Σₘₙ aₘbₙ Pₘ⊗Pₙ`: the normal block `a bᵀ`, no shear. A self product
# is major-symmetric (nine constants), a product of two tensors generally not
# (twelve), as `TensTI{2} ⊗ TensTI{2}` gives N=5 and N=6.

function Tensors.otimes(A::TensOrtho{2})
    a = SVector(get_data(A))
    z = zero(eltype(a))
    O = _ortho_from_blocks(a * transpose(a), (z, z, z), frame(A))
    return Tensors.majorsymmetric(O)
end
function Tensors.otimes(A::TensOrtho{2}, B::TensOrtho{2})
    frame(A) == frame(B) || return Tensors.otimes(generic_tens(A), generic_tens(B))
    a, b = SVector(get_data(A)), SVector(get_data(B))
    z = zero(promote_type(eltype(a), eltype(b)))
    return _ortho_from_blocks(a * transpose(b), (z, z, z), frame(A))
end
Tensors.otimes(A::TensISO{2, 3}, B::TensOrtho{2}) = Tensors.otimes(iso_to_ortho(A, frame(B)), B)
Tensors.otimes(A::TensOrtho{2}, B::TensISO{2, 3}) = Tensors.otimes(A, iso_to_ortho(B, frame(A)))

# ──────────────────────────────────────────────────────────────────────────────
# Exports
# ──────────────────────────────────────────────────────────────────────────────

export iso_to_ortho, walpole_to_ortho, iso_to_cubic, cubic_to_ortho
