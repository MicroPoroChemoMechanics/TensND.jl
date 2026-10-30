##############################################################################
# tens_anisotropic.jl — the structured ANISOTROPIC tensor types, the           #
# counterpart of `tens_isotropic.jl`:                                          #
#                                                                              #
#   TensTI     transversely isotropic tensors, stored in the Walpole basis     #
#   TensOrtho  orthotropic tensors, nine or twelve constants plus a frame      #
##############################################################################

# ─────────────────────────────────────────────────────────────────────────────
# TensTI{4}
# ─────────────────────────────────────────────────────────────────────────────
#
# A transversely isotropic (TI) 4th-order tensor with symmetry axis n can be
# written in the Walpole basis {W₁,…,W₆} as
#
#   L = ℓ₁W₁ + ℓ₂W₂ + ℓ₃W₃ + ℓ₄W₄ + ℓ₅W₅ + ℓ₆W₆
#
# where (nₙ = n⊗n, nT = 1 − nₙ):
#   W₁ = nₙ⊗nₙ
#   W₂ = (nT⊗nT)/2
#   W₃ = (nₙ⊗nT)/√2
#   W₄ = (nT⊗nₙ)/√2
#   W₅ = nT⊠ˢnT − (nT⊗nT)/2
#   W₆ = nT⊠ˢnₙ + nₙ⊠ˢnT
#
# For major-symmetric tensors ℓ₃ = ℓ₄ → stored with N=5 data scalars.
# General (non-major-sym) tensors use N=6.
#
# Synthetic notation: L ≡ ([[ℓ₁,ℓ₃],[ℓ₄,ℓ₂]], ℓ₅, ℓ₆)
#   Product:  (L⊡M)_mat = L_mat × M_mat  ,  (L⊡M)₅ = ℓ₅m₅ , (L⊡M)₆ = ℓ₆m₆
#   Inverse:  (L⁻¹)_mat = (L_mat)⁻¹      ,  1/ℓ₅           , 1/ℓ₆
# ─────────────────────────────────────────────────────────────────────────────

# ──────────────────────────────────────────────────────────────────────────────
# TensTI — parametric transversely isotropic tensor (order 2 or 4)
# ──────────────────────────────────────────────────────────────────────────────

"""
    TensTI{order, T, N} <: AbstractTens{order, 3, T}

Transversely isotropic tensor of order `order` (always `dim=3`) with symmetry
axis `n`, parametrized like `TensISO{order, dim, T, N}`: `N` is the number of
stored coefficients, and the invariance group fixes which tensors a shape holds.

Five concrete shapes are supported:

| Parametrization   | Role                                   | Invariance group | Stored coefficients              |
|-------------------|----------------------------------------|------------------|----------------------------------|
| `TensTI{2, T, 2}` | 2nd-order TI                           | `D∞h`            | `data = (a, b)`, `n::NTuple{3}`  |
| `TensTI{2, T, 3}` | 2nd-order, axially invariant           | `SO(2)`          | `data = (a, b, c)`, `n`          |
| `TensTI{4, T, 5}` | 4th-order TI, major-symmetric          | `D∞h`            | `data = (ℓ₁,ℓ₂,ℓ₃,ℓ₅,ℓ₆)`, `n`   |
| `TensTI{4, T, 6}` | 4th-order TI, general                  | `D∞h`            | `data = (ℓ₁,…,ℓ₆)`, `n`          |
| `TensTI{4, T, 8}` | 4th-order, axially invariant           | `SO(2)`          | `data = (ℓ₁,…,ℓ₈)`, `n`          |

`D∞h` is the full transverse-isotropy group (rotations about `n`, mirrors
through and across it); `SO(2)` keeps only the rotations about `n`, and with
them the terms odd in `n` (`c·w`, `ℓ₇W₇`, `ℓ₈W₈`, built on the in-plane
rotation generator `w`). The two orders correspond: `N=2` ↔ `N=5, 6`, `N=3` ↔
`N=8` (the contraction of an order-4 shape with an order-2 one of the same
group stays in it).

- **Order 2** (`N=2`): `𝐀 = a·nT + b·nₙ` where `nₙ = n⊗n`, `nT = 𝟏 − nₙ`.
  `a` is the transverse coefficient, `b` the axial one.  When `a = b`, the
  tensor is isotropic and equivalent to `TensISO{2,3,T}(a)`.

- **Order 2, `N=3`**: `𝐀 = a·nT + b·nₙ + c·w`, the skew part `c·w` invariant
  under the rotations about `n` but not under the mirrors.

- **Order 4**: stored in the Walpole basis `{W₁,…,W₆}`, with
  `W₁ = nₙ⊗nₙ`, `W₂ = (nT⊗nT)/2`, `W₃ = (nₙ⊗nT)/√2`,
  `W₄ = (nT⊗nₙ)/√2`, `W₅ = nT⊠ˢnT − (nT⊗nT)/2`, `W₆ = nT⊠ˢnₙ + nₙ⊠ˢnT`.
  Major-symmetric tensors have `ℓ₃ = ℓ₄` and are stored under `N=5`.
  Synthetic notation: `L ≡ ([[ℓ₁,ℓ₃],[ℓ₄,ℓ₂]], ℓ₅, ℓ₆)`. `N=8` adds the two
  major-antisymmetric generators `W₇`, `W₈` of the `SO(2)` commutant.

The product of two shapes of the same order returns the smallest shape holding
it (`N=5 ⊡ N=5 → N=6`), the inverse keeps the shape, the major transpose swaps
`ℓ₃` and `ℓ₄` (and changes the sign of `ℓ₇`, `ℓ₈`), and
`Tensors.majorsymmetric` returns `N=5`.

See also [`TensOrtho`](@ref), [`ti_from_constants`](@ref).
"""
struct TensTI{order, T, N} <: AbstractTens{order, 3, T}
    data::NTuple{N, T}
    n::NTuple{3, T}       # symmetry axis (assumed unit vector)
    function TensTI{order, T, N}(data::NTuple{N, T}, n::NTuple{3, T}) where {order, T, N}
        order == 4 && T !== sqrt2_eltype(T) && throw(
            ArgumentError(
                "a fourth-order TensTI stores its Walpole coefficients in the type of its " *
                    "components, which carry 1/√2: $T cannot hold them; use " *
                    "$(sqrt2_eltype(T)), which the constructors `TensTI{4}(ℓ…, n)` choose"
            )
        )
        return new{order, T, N}(data, n)
    end
end

# Mixed-eltype constructor.  The symmetry axis is a piece of GEOMETRY and is
# routinely a plain `Float64` (`(0.0, 0.0, 1.0)`) while the Walpole
# coefficients get promoted — most importantly to `ForwardDiff.Dual` when a
# scheme is differentiated with respect to a modulus.  Requiring a single `T`
# for both made every such call fail with
# `MethodError: no method matching TensTI{4, Dual, N}(::NTuple{N, Dual}, ::NTuple{3, Float64})`.
# Convert both tuples to `T` instead.
TensTI{order, T, N}(data::NTuple{N, Any}, n::NTuple{3, Any}) where {order, T, N} =
    TensTI{order, T, N}(map(T, data), map(T, n))

# ── Traits ────────────────────────────────────────────────────────────────────

Base.eltype(::Type{TensTI{order, T, N}}) where {order, T, N} = T
Base.length(::TensTI{order}) where {order} = 3^order
Base.size(::TensTI{order}) where {order} = ntuple(_ -> 3, Val(order))

get_basis(::TensTI{order, T}) where {order, T} = CanonicalBasis{3, T}()
get_var(::TensTI{order}) where {order} = ntuple(_ -> :cont, Val(order))
get_var(::TensTI, ::Integer) = :cont
get_data(t::TensTI) = t.data
"""
    axis(t::TensTI)

Return the symmetry axis of a transversely isotropic tensor.
"""
axis(t::TensTI) = t.n

# ── Rebuild helper (used by symbolic ops) ─────────────────────────────────────
_rebuild(t::TensTI{order}, new_data) where {order} =
    TensTI{order, _ti_eltype(Val(order), eltype(new_data)), length(new_data)}(new_data, axis(t))

# The element type of a TensTI of order `order` holding coefficients of type `T`:
# those of order 4 carry 1/√2 in their components (see `sqrt2_eltype`).
_ti_eltype(::Val{4}, ::Type{T}) where {T} = sqrt2_eltype(T)
_ti_eltype(::Val, ::Type{T}) where {T} = T

# ── 4th-order TI accessors ───────────────────────────────────────────────────

"""
    walpole_coeffs(t::TensTI{4,T,N}) → NTuple{6,T}

Always returns a 6-tuple `(ℓ₁, ℓ₂, ℓ₃, ℓ₄, ℓ₅, ℓ₆)` of Walpole coefficients.
For `N=5` (major-symmetric), `ℓ₃ = ℓ₄` is stored once and duplicated on read.
"""
walpole_coeffs(t::TensTI{4, T, 5}) where {T} =
    (t.data[1], t.data[2], t.data[3], t.data[3], t.data[4], t.data[5])
walpole_coeffs(t::TensTI{4, T, 6}) where {T} = t.data


# ── Constructors ──────────────────────────────────────────────────────────────

# TensTI{4}(ℓ₁,ℓ₂,ℓ₃,ℓ₄,ℓ₅,ℓ₆, n) → TensTI{4, T, 6}
# General (not necessarily major-symmetric) TI 4th-order tensor with axis `n`.
function TensTI{4}(ℓ₁, ℓ₂, ℓ₃, ℓ₄, ℓ₅, ℓ₆, n)
    T = sqrt2_eltype(
        promote_type(
            typeof(ℓ₁), typeof(ℓ₂), typeof(ℓ₃), typeof(ℓ₄),
            typeof(ℓ₅), typeof(ℓ₆), eltype(n)
        )
    )
    nv = _extract_vec(n)
    return TensTI{4, T, 6}(
        (T(ℓ₁), T(ℓ₂), T(ℓ₃), T(ℓ₄), T(ℓ₅), T(ℓ₆)),
        (T(nv[1]), T(nv[2]), T(nv[3]))
    )
end

# TensTI{4}(ℓ₁,ℓ₂,ℓ₃,ℓ₅,ℓ₆, n) → TensTI{4, T, 5}
# Major-symmetric TI 4th-order tensor (ℓ₃ = ℓ₄), 5 independent scalars.
function TensTI{4}(ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆, n)
    T = sqrt2_eltype(
        promote_type(
            typeof(ℓ₁), typeof(ℓ₂), typeof(ℓ₃),
            typeof(ℓ₅), typeof(ℓ₆), eltype(n)
        )
    )
    nv = _extract_vec(n)
    return TensTI{4, T, 5}(
        (T(ℓ₁), T(ℓ₂), T(ℓ₃), T(ℓ₅), T(ℓ₆)),
        (T(nv[1]), T(nv[2]), T(nv[3]))
    )
end

# Extract a plain 3-vector from various input types
_extract_vec(n::NTuple{3}) = n
# Mixed element types (`(sin(θ), 0, cos(θ))`): promoted to one.
_extract_vec(n::Tuple{Number, Number, Number}) = promote(n...)
_extract_vec(n::AbstractVector) = (n[1], n[2], n[3])
# Canonical components: an axis is canonical, whatever basis the vector that
# gives it is written in.
_extract_vec(n::AbstractTens) = _extract_vec(canonical_components(n))
_extract_vec(n::Vec{3}) = (n[1], n[2], n[3])
_extract_vec(n::AbstractArray) = (n[1], n[2], n[3])

# ── Basis tensors Wᵢ ─────────────────────────────────────────────────────────

"""
    tens_W1(n) → TensTI{4, T, 6}   (W₁ = nₙ⊗nₙ, coeffs (1,0,0,0,0,0))
"""
tens_W1(n) = TensTI{4}(
    one(eltype_of(n)), zero(eltype_of(n)), zero(eltype_of(n)),
    zero(eltype_of(n)), zero(eltype_of(n)), zero(eltype_of(n)), n
)

"""
    tens_W2(n) → TensTI{4, T, 6}   (W₂ = (nT⊗nT)/2, coeffs (0,1,0,0,0,0))
"""
tens_W2(n) = TensTI{4}(
    zero(eltype_of(n)), one(eltype_of(n)), zero(eltype_of(n)),
    zero(eltype_of(n)), zero(eltype_of(n)), zero(eltype_of(n)), n
)

"""
    tens_W3(n) → TensTI{4, T, 6}   (W₃ = (nₙ⊗nT)/√2, coeffs (0,0,1,0,0,0))
"""
tens_W3(n) = TensTI{4}(
    zero(eltype_of(n)), zero(eltype_of(n)), one(eltype_of(n)),
    zero(eltype_of(n)), zero(eltype_of(n)), zero(eltype_of(n)), n
)

"""
    tens_W4(n) → TensTI{4, T, 6}   (W₄ = (nT⊗nₙ)/√2, coeffs (0,0,0,1,0,0))
"""
tens_W4(n) = TensTI{4}(
    zero(eltype_of(n)), zero(eltype_of(n)), zero(eltype_of(n)),
    one(eltype_of(n)), zero(eltype_of(n)), zero(eltype_of(n)), n
)

"""
    tens_W5(n) → TensTI{4, T, 6}   (W₅ = nT⊠ˢnT − (nT⊗nT)/2, coeffs (0,0,0,0,1,0))
"""
tens_W5(n) = TensTI{4}(
    zero(eltype_of(n)), zero(eltype_of(n)), zero(eltype_of(n)),
    zero(eltype_of(n)), one(eltype_of(n)), zero(eltype_of(n)), n
)

"""
    tens_W6(n) → TensTI{4, T, 6}   (W₆ = nT⊠ˢnₙ + nₙ⊠ˢnT, coeffs (0,0,0,0,0,1))
"""
tens_W6(n) = TensTI{4}(
    zero(eltype_of(n)), zero(eltype_of(n)), zero(eltype_of(n)),
    zero(eltype_of(n)), zero(eltype_of(n)), one(eltype_of(n)), n
)

# Helper: get element type from various axis representations
eltype_of(::AbstractArray{T}) where {T} = T
# `NTuple{N, T}` also matches the empty tuple, where `T` is bound by nothing.
eltype_of(::Tuple{T, Vararg{T}}) where {T} = T
eltype_of(::AbstractTens{1, 3, T}) where {T} = T

"""
    walpole_basis(n) → (W₁, W₂, W₃, W₄, W₅, W₆)

Return the six general (`N=6`) Walpole basis tensors for the symmetry axis
`n`. These span the full TI 4th-order tensor space, including the
non-major-symmetric components `W₃ ≠ W₄`.
"""
walpole_basis(n) = (tens_W1(n), tens_W2(n), tens_W3(n), tens_W4(n), tens_W5(n), tens_W6(n))

"""
    walpole_basis_sym(n) → (W₁ˢ, W₂ˢ, W₃ˢ, W₄ˢ, W₅ˢ)

Return the five major-symmetric (`N=5`) Walpole basis tensors for the symmetry
axis `n`, where `W₃ˢ = W₃ + W₄`. Use this for building stiffness / compliance
TI tensors, which are always major-symmetric.
"""
function walpole_basis_sym(n)
    T = eltype_of(n)
    o, z = one(T), zero(T)
    W1s = TensTI{4}(o, z, z, z, z, n)         # ℓ₁=1
    W2s = TensTI{4}(z, o, z, z, z, n)         # ℓ₂=1
    W3s = TensTI{4}(z, z, o, z, z, n)         # ℓ₃=1  (W₃+W₄)
    W4s = TensTI{4}(z, z, z, o, z, n)         # ℓ₅=1
    W5s = TensTI{4}(z, z, z, z, o, n)         # ℓ₆=1
    return W1s, W2s, W3s, W4s, W5s
end

# ── get_array ─────────────────────────────────────────────────────────────────

"""
    get_array(t::TensTI{4, T}) → Array{T,4}

Compute the 3×3×3×3 component array from the Walpole coefficients and axis.
"""
function get_array(t::TensTI{4, T}) where {T}
    ℓ₁, ℓ₂, ℓ₃, ℓ₄, ℓ₅, ℓ₆ = walpole_coeffs(t)
    n = t.n
    sq2 = sqrt(T(2))
    δ(i, j) = i == j ? one(T) : zero(T)
    nn(i, j) = n[i] * n[j]
    nT(i, j) = δ(i, j) - nn(i, j)
    result = Array{T, 4}(undef, 3, 3, 3, 3)
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        W1 = nn(i, j) * nn(k, l)
        W2 = nT(i, j) * nT(k, l) / 2
        W3 = nn(i, j) * nT(k, l) / sq2
        W4 = nT(i, j) * nn(k, l) / sq2
        W5 = (nT(i, k) * nT(j, l) + nT(i, l) * nT(j, k)) / 2 - nT(i, j) * nT(k, l) / 2
        W6 = (nT(i, k) * nn(j, l) + nT(i, l) * nn(j, k) + nn(i, k) * nT(j, l) + nn(i, l) * nT(j, k)) / 2
        result[i, j, k, l] = ℓ₁ * W1 + ℓ₂ * W2 + ℓ₃ * W3 + ℓ₄ * W4 + ℓ₅ * W5 + ℓ₆ * W6
    end
    return result
end

function Base.getindex(t::TensTI{4}, i::Integer, j::Integer, k::Integer, l::Integer)
    ℓ₁, ℓ₂, ℓ₃, ℓ₄, ℓ₅, ℓ₆ = walpole_coeffs(t)
    n = t.n
    T = eltype(t)
    sq2 = sqrt(T(2))
    δ(a, b) = a == b ? one(T) : zero(T)
    nn(a, b) = n[a] * n[b]
    nT(a, b) = δ(a, b) - nn(a, b)
    W1 = nn(i, j) * nn(k, l)
    W2 = nT(i, j) * nT(k, l) / 2
    W3 = nn(i, j) * nT(k, l) / sq2
    W4 = nT(i, j) * nn(k, l) / sq2
    W5 = (nT(i, k) * nT(j, l) + nT(i, l) * nT(j, k)) / 2 - nT(i, j) * nT(k, l) / 2
    W6 = (nT(i, k) * nn(j, l) + nT(i, l) * nn(j, k) + nn(i, k) * nT(j, l) + nn(i, l) * nT(j, k)) / 2
    return ℓ₁ * W1 + ℓ₂ * W2 + ℓ₃ * W3 + ℓ₄ * W4 + ℓ₅ * W5 + ℓ₆ * W6
end

# ── Kelvin-Mandel matrix ──────────────────────────────────────────────────────

"""
    _mandel(t::TensTI{4})

Kelvin-Mandel (6×6) matrix of the Walpole tensor.
"""
_mandel(t::TensTI{4}) = Matrix(_km_by_entries(t))

# ── Arithmetic ────────────────────────────────────────────────────────────────
# Scalar ops (-, α*A, A*α, A/α) defined in structured_tens_ops.jl.
# Binary ± (same-N, mixed-N, mixed-axis fallback) implemented in the unified
# block further below (see "Arithmetic — axis-aware ±").

# ── Double contraction (Walpole product rule) ─────────────────────────────────

"""
    dcontract(A::TensTI{4}, B::TensTI{4}) → TensTI{4, T, 6}

Product rule via 2×2 matrix product + scalar products for ℓ₅, ℓ₆.
Always returns N=6 since the product of two symmetric tensors need not be symmetric.

If the axes differ, the operation falls back to the generic (unstructured)
`Tens` route via `get_array` — the product of two TI tensors with different
axes is generally fully anisotropic.
"""
function Tensors.dcontract(A::TensTI{4}, B::TensTI{4})
    A.n == B.n || return Tensors.dcontract(generic_tens(A), generic_tens(B))
    return _dcontract_coaxial(A, B)
end

# The Walpole product of two TI tensors sharing their axis (not checked): the
# callers that build one operand on the other's axis call it directly, so that
# the result is inferred as a `TensTI` rather than a union with the generic route.
function _dcontract_coaxial(A::TensTI{4}, B::TensTI{4})
    ℓA₁, ℓA₂, ℓA₃, ℓA₄, ℓA₅, ℓA₆ = walpole_coeffs(A)
    ℓB₁, ℓB₂, ℓB₃, ℓB₄, ℓB₅, ℓB₆ = walpole_coeffs(B)
    # 2×2 matrix rule: M_A × M_B where M = [[ℓ₁,ℓ₃],[ℓ₄,ℓ₂]]
    n₁ = ℓA₁ * ℓB₁ + ℓA₃ * ℓB₄
    n₃ = ℓA₁ * ℓB₃ + ℓA₃ * ℓB₂
    n₄ = ℓA₄ * ℓB₁ + ℓA₂ * ℓB₄
    n₂ = ℓA₄ * ℓB₃ + ℓA₂ * ℓB₂
    n₅ = ℓA₅ * ℓB₅
    n₆ = ℓA₆ * ℓB₆
    T = promote_type(eltype(A), eltype(B))
    return TensTI{4, T, 6}((T(n₁), T(n₂), T(n₃), T(n₄), T(n₅), T(n₆)), A.n)
end

# ── Inverse ───────────────────────────────────────────────────────────────────

"""
    inv(t::TensTI{4, T, 5}) → TensTI{4, T, 5}
    inv(t::TensTI{4, T, 6}) → TensTI{4, T, 6}

Inverse via the 2×2 Walpole matrix and scalar inverses for ℓ₅, ℓ₆.
"""
function Base.inv(t::TensTI{4, T, 5}) where {T}
    ℓ₁, ℓ₂, ℓ₃, _, ℓ₅, ℓ₆ = walpole_coeffs(t)   # ℓ₄=ℓ₃ for N=5
    det = ℓ₁ * ℓ₂ - ℓ₃ * ℓ₃
    return TensTI{4, T, 5}((ℓ₂ / det, ℓ₁ / det, -ℓ₃ / det, one(T) / ℓ₅, one(T) / ℓ₆), t.n)
end

function Base.inv(t::TensTI{4, T, 6}) where {T}
    ℓ₁, ℓ₂, ℓ₃, ℓ₄, ℓ₅, ℓ₆ = walpole_coeffs(t)
    det = ℓ₁ * ℓ₂ - ℓ₃ * ℓ₄
    return TensTI{4, T, 6}((ℓ₂ / det, ℓ₁ / det, -ℓ₃ / det, -ℓ₄ / det, one(T) / ℓ₅, one(T) / ℓ₆), t.n)
end

@inline Base.literal_pow(::typeof(^), A::TensTI{4}, ::Val{-1}) = inv(A)

# ── Symmetry tests ────────────────────────────────────────────────────────────

LinearAlgebra.issymmetric(::TensTI{4, T, 5}) where {T} = true
LinearAlgebra.issymmetric(t::TensTI{4, T, 6}) where {T} = isequal(t.data[3], t.data[4])
Tensors.isminorsymmetric(::TensTI{4}) = true
Tensors.ismajorsymmetric(::TensTI{4, T, 5}) where {T} = true
Tensors.ismajorsymmetric(t::TensTI{4, T, 6}) where {T} = isequal(t.data[3], t.data[4])

# ── iso_to_ti ───────────────────────────────────────────────────────────────────

"""
    iso_to_ti(A::TensISO{4,3}, n) → TensTI{4, T, 5}

Convert an isotropic 4th-order tensor `αJ + βK` into its Walpole representation.

Formulas: ℓ₁=(α+2β)/3, ℓ₂=(2α+β)/3 (note: dim=3 → these are (3k,2μ) related),
          ℓ₃=ℓ₄=√2(α−β)/3, ℓ₅=ℓ₆=β.
Here `α` = data[1] and `β` = data[2] in TensISO (coefficients of J and K).
"""
function iso_to_ti(A::TensISO{4, 3, T}, n) where {T}
    α, β = get_data(A)    # A = α*J + β*K
    sq2 = sqrt(T(2))
    ℓ₁ = (α + 2β) / 3
    ℓ₂ = (2α + β) / 3   # Note: for 3D, 1-1/dim = 2/3 and 1/dim = 1/3
    ℓ₃ = sq2 * (α - β) / 3
    ℓ₅ = β
    ℓ₆ = β
    return TensTI{4}(ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆, n)
end

"""
    dcontract(A::TensTI{4}, B::TensISO{4,3}) → TensTI{4, T, 6}
    dcontract(A::TensISO{4,3}, B::TensTI{4}) → TensTI{4, T, 6}
"""
function Tensors.dcontract(A::TensTI{4}, B::TensISO{4, 3})
    return _dcontract_coaxial(A, iso_to_ti(B, A.n))
end
function Tensors.dcontract(A::TensISO{4, 3}, B::TensTI{4})
    return _dcontract_coaxial(iso_to_ti(A, B.n), B)
end

# ── TI constants ─────────────────────────────────────────────────────────────

"""
    ti_from_constants(c₁, c₂, c₃, c₄, c₅, n; form = :components) → TensTI{4, T, 5}

The major-symmetric transversely isotropic tensor of axis `n` given by five
constants, read according to `form`:

- `:components` (the default): its components
  `(C₁₁₁₁, C₁₁₂₂, C₁₁₃₃, C₃₃₃₃, C₂₃₂₃)` in a frame whose third axis is `n`, for
  a stiffness or a compliance alike. Its Walpole coefficients are
  `ℓ₁ = C₃₃₃₃`, `ℓ₂ = C₁₁₁₁ + C₁₁₂₂`, `ℓ₃ = √2 C₁₁₃₃`, `ℓ₅ = C₁₁₁₁ − C₁₁₂₂`,
  `ℓ₆ = 2 C₂₃₂₃`.
- `:engineering`: the **compliance** of the engineering constants
  `(E₁, E₃, ν₁₂, ν₃₁, G₃₁)`: the transverse and axial Young's moduli, the
  in-plane and axial-transverse Poisson's ratios (`ν₃₁/E₃ = ν₁₃/E₁`) and the
  axial shear modulus, so that `S₁₁₁₁ = 1/E₁`, `S₃₃₃₃ = 1/E₃`,
  `S₁₁₂₂ = −ν₁₂/E₁`, `S₁₁₃₃ = −ν₃₁/E₃`, `S₂₃₂₃ = 1/(4G₃₁)`.
- `:hoenig`: the **compliance** of Hoenig's parameters `(E, ν₁, ν₂, H, Γ)`
  (Hoenig, 1978): `S₁₁₁₁ = 1/E`, `S₁₁₂₂ = −ν₁/E`, `S₁₁₃₃ = −ν₂/E`,
  `S₃₃₃₃ = 1/(E H)`, `S₂₃₂₃ = (1 + ν₁)/(2 E Γ)`.

The stiffness of the last two is `inv` of the result. [`ti_constants`](@ref)
reads the constants back.
"""
function ti_from_constants(c₁, c₂, c₃, c₄, c₅, n; form::Symbol = :components)
    form === :components && return _ti_from_components(c₁, c₂, c₃, c₄, c₅, n)
    form === :engineering && return _ti_from_components(_compliance_engineering(c₁, c₂, c₃, c₄, c₅)..., n)
    form === :hoenig && return _ti_from_components(_compliance_hoenig(c₁, c₂, c₃, c₄, c₅)..., n)
    throw(_unknown_ti_form(form))
end

"""
    ti_constants(t::TensTI{4}; form = :components) → NTuple{5}

The five constants of a transversely isotropic tensor, the reciprocal of
[`ti_from_constants`](@ref) for the same `form`: its components
`(C₁₁₁₁, C₁₁₂₂, C₁₁₃₃, C₃₃₃₃, C₂₃₂₃)` in its axis frame, read from the stored
Walpole coefficients (`C₃₃₃₃ = ℓ₁`, `C₁₁₁₁ = (ℓ₂ + ℓ₅)/2`,
`C₁₁₂₂ = (ℓ₂ − ℓ₅)/2`, `C₁₁₃₃ = ℓ₃/√2`, `C₂₃₂₃ = ℓ₆/2`); with
`form = :engineering`, `(E₁, E₃, ν₁₂, ν₃₁, G₃₁)` of the **compliance** `t`; with
`form = :hoenig`, `(E, ν₁, ν₂, H, Γ)` of the **compliance** `t`.

A six- or eight-parameter tensor (a product of TI tensors lands there) has
these five constants only when it is major-symmetric, ``\\ell_3 = \\ell_4`` and
``\\ell_7 = \\ell_8 = 0``: exactly for symbolic coefficients, to rounding for
numbers. Otherwise it is refused, rather than read as if `ℓ₄`, `ℓ₇` and `ℓ₈` were
not there; take `Tensors.majorsymmetric(t)` first if its symmetric part is meant.
"""
function ti_constants(t::TensTI{4}; form::Symbol = :components)
    form === :components && return _ti_components(t)
    form === :engineering && return _engineering_of_compliance(_ti_components(t)...)
    form === :hoenig && return _hoenig_of_compliance(_ti_components(t)...)
    throw(_unknown_ti_form(form))
end

_unknown_ti_form(form) =
    ArgumentError("unknown form `:$form` of the TI constants; use `:components`, `:engineering` or `:hoenig`")

function _ti_from_components(C₁₁₁₁, C₁₁₂₂, C₁₁₃₃, C₃₃₃₃, C₂₃₂₃, n)
    T = promote_type(
        typeof(C₁₁₁₁), typeof(C₁₁₂₂), typeof(C₁₁₃₃),
        typeof(C₃₃₃₃), typeof(C₂₃₂₃)
    )
    sq2 = sqrt(T(2))
    ℓ₁ = C₃₃₃₃
    ℓ₂ = C₁₁₁₁ + C₁₁₂₂
    ℓ₃ = sq2 * C₁₁₃₃
    ℓ₅ = C₁₁₁₁ - C₁₁₂₂
    ℓ₆ = 2 * C₂₃₂₃
    return TensTI{4}(ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆, n)
end

function _ti_components(t::TensTI{4})
    ℓ₁, ℓ₂, ℓ₃, ℓ₄, ℓ₅, ℓ₆, ℓ₇, ℓ₈ = walpole_coeffs(t, Val(8))
    _major_symmetric_ti(ℓ₃ - ℓ₄, ℓ₇, ℓ₈, max(abs(ℓ₁), abs(ℓ₂), abs(ℓ₃), abs(ℓ₅), abs(ℓ₆))) || throw(
        ArgumentError(
            "ti_constants: the tensor is not major-symmetric (ℓ₃ ≠ ℓ₄, or ℓ₇, ℓ₈ ≠ 0), so it has " *
                "no five TI constants; take `Tensors.majorsymmetric(t)` first if its " *
                "symmetric part is meant"
        )
    )
    return _ti_components(TensTI{4, eltype(t), 5}((ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆), axis(t)))
end

function _ti_components(t::TensTI{4, T, 5}) where {T}
    ℓ₁, ℓ₂, ℓ₃, _, ℓ₅, ℓ₆ = walpole_coeffs(t)
    sq2 = sqrt(T(2))
    C₃₃₃₃ = ℓ₁
    C₁₁₁₁ = (ℓ₂ + ℓ₅) / 2
    C₁₁₂₂ = (ℓ₂ - ℓ₅) / 2
    C₁₁₃₃ = ℓ₃ / sq2
    C₂₃₂₃ = ℓ₆ / 2
    return (C₁₁₁₁, C₁₁₂₂, C₁₁₃₃, C₃₃₃₃, C₂₃₂₃)
end

# Whether the major-antisymmetric generators vanish: to rounding against the
# scale `s` for comparable numbers, structurally for symbolic ones.
_major_symmetric_ti(d, ℓ₇, ℓ₈, s) = _major_symmetric_ti(Val(is_hard_numeric(typeof(d))), d, ℓ₇, ℓ₈, s)
_major_symmetric_ti(::Val{true}, d, ℓ₇, ℓ₈, s) =
    max(abs(d), abs(ℓ₇), abs(ℓ₈)) <= 16 * eps(float(primal(s)))      # 16 ulps of the scale
_major_symmetric_ti(::Val{false}, d, ℓ₇, ℓ₈, s) =
    all(x -> isequal(tsimplify(x), zero(x)), (d, ℓ₇, ℓ₈))

# The compliance components of the engineering constants, and back.
_compliance_engineering(E₁, E₃, ν₁₂, ν₃₁, G₃₁) = (inv(E₁), -ν₁₂ / E₁, -ν₃₁ / E₃, inv(E₃), inv(4 * G₃₁))

function _engineering_of_compliance(S₁₁₁₁, S₁₁₂₂, S₁₁₃₃, S₃₃₃₃, S₂₃₂₃)
    E₁ = inv(S₁₁₁₁)
    E₃ = inv(S₃₃₃₃)
    return (E₁, E₃, -E₁ * S₁₁₂₂, -E₃ * S₁₁₃₃, inv(4 * S₂₃₂₃))
end

# The compliance components of Hoenig's parameters, and back.
_compliance_hoenig(E, ν₁, ν₂, H, Γ) = (inv(E), -ν₁ / E, -ν₂ / E, inv(E * H), (1 + ν₁) / (2 * E * Γ))

function _hoenig_of_compliance(S₁₁₁₁, S₁₁₂₂, S₁₁₃₃, S₃₃₃₃, S₂₃₂₃)
    E = inv(S₁₁₁₁)
    ν₁ = -E * S₁₁₂₂
    return (E, ν₁, -E * S₁₁₃₃, inv(S₃₃₃₃ * E), (1 + ν₁) / (2 * E * S₂₃₂₃))
end

# ── isisotropic / istransverselyisotropic ─────────────────────────────────────────────────────────────

"""
    istransverselyisotropic(A)

Return `true` if `A` is a `TensTI{4}`, indicating transverse isotropy.
"""
istransverselyisotropic(::TensTI{4}) = true
istransverselyisotropic(::Any) = false
isisotropic(::TensTI{4}) = false
isorthotropic(::TensTI{4}) = false

# Symbolic helpers (tsimplify, tsubs, …) defined in structured_tens_ops.jl

# ── Display ───────────────────────────────────────────────────────────────────

function Base.show(io::IO, A::TensTI{4, <:Any, 5})
    ℓ₁, ℓ₂, ℓ₃, _, ℓ₅, ℓ₆ = walpole_coeffs(A)
    print(
        io, "(", ℓ₁, ") W₁ˢ + (", ℓ₂, ") W₂ˢ + (", ℓ₃,
        ") W₃ˢ + (", ℓ₅, ") W₄ˢ + (", ℓ₆, ") W₅ˢ"
    )
    return print(io, "\n  axis n = ", A.n)
end
function Base.show(io::IO, A::TensTI{4, <:Any, 6})
    ℓ₁, ℓ₂, ℓ₃, ℓ₄, ℓ₅, ℓ₆ = walpole_coeffs(A)
    print(
        io, "(", ℓ₁, ") W₁ + (", ℓ₂, ") W₂ + (", ℓ₃,
        ") W₃ + (", ℓ₄, ") W₄ + (", ℓ₅, ") W₅ + (", ℓ₆, ") W₆"
    )
    return print(io, "\n  axis n = ", A.n)
end

function pprint(A::TensTI{4, <:Any, 5})
    ℓ₁, ℓ₂, ℓ₃, _, ℓ₅, ℓ₆ = walpole_coeffs(A)
    println(
        "(", ℓ₁, ") W₁ˢ + (", ℓ₂, ") W₂ˢ + (", ℓ₃,
        ") W₃ˢ + (", ℓ₅, ") W₄ˢ + (", ℓ₆, ") W₅ˢ"
    )
    return println("  axis n = ", A.n)
end
function pprint(A::TensTI{4, <:Any, 6})
    ℓ₁, ℓ₂, ℓ₃, ℓ₄, ℓ₅, ℓ₆ = walpole_coeffs(A)
    println(
        "(", ℓ₁, ") W₁ + (", ℓ₂, ") W₂ + (", ℓ₃,
        ") W₃ + (", ℓ₄, ") W₄ + (", ℓ₅, ") W₅ + (", ℓ₆, ") W₆"
    )
    return println("  axis n = ", A.n)
end


##############################################################################
# TensTI — transversely isotropic tensor (parametric on order)              #
##############################################################################
#
# Follows the same parametric-order pattern as TensISO{order,dim,T,N}:
#
# TensTI{order,T,N} <: AbstractTens{order,3,T}
#
# Order 2 (N=2):  𝐀 = a·nT + b·nₙ  where nₙ=n⊗n, nT=𝟏−nₙ
#   data = (a, b), a = transverse coeff, b = axial coeff
#   When a = b, isotropic: 𝐀 = a·𝟏 (equiv. TensISO{2,3,T}(a))
#
# Order 4 is handled by TensTI{4, T, N} (Walpole basis, N=5 or 6).
# Future unification TensTI{4} → TensTI{4,T,N} is possible.
# ─────────────────────────────────────────────────────────────────────────────

# (The TensTI struct + traits are defined at the top of this file.)

# ── Convenience constructors ─────────────────────────────────────────────────

# TensTI{2}(a, b, n) → TensTI{2,T,2}
#
# Construct a TI 2nd-order tensor `a·nT + b·nₙ` with symmetry axis `n`.
#
# Examples:
#   n = [0., 0., 1.]
#   A = TensTI{2}(5.0, 8.0, n)
#   get_array(A) → [5 0 0; 0 5 0; 0 0 8]
function TensTI{2}(a, b, n)
    T = promote_type(typeof(a), typeof(b), eltype(n))
    nv = _extract_vec(n)
    return TensTI{2, T, 2}((T(a), T(b)), (T(nv[1]), T(nv[2]), T(nv[3])))
end

# ── get_array (order 2) ──────────────────────────────────────────────────────

"""
    get_array(t::TensTI{2,T,2}) → Array{T,2}

Compute the 3×3 component array: `a*(δᵢⱼ − nᵢnⱼ) + b*nᵢnⱼ`.

# Examples
```julia
julia> A = TensTI{2}(5.0, 8.0, [0., 0., 1.]);

julia> get_array(A)
3×3 Matrix{Float64}:
 5.0  0.0  0.0
 0.0  5.0  0.0
 0.0  0.0  8.0
```
"""
function get_array(t::TensTI{2, T, 2}) where {T}
    a, b = t.data
    n = SVector(t.n)
    # Static algebra copied once into the result: faster than filling it entry
    # by entry, and `(b n) nᵀ` keeps the association of the component formula.
    return copyto!(Matrix{T}(undef, 3, 3), a * (I - n * transpose(n)) + (b * n) * transpose(n))
end

function Base.getindex(t::TensTI{2, T, 2}, i::Integer, j::Integer) where {T}
    a, b = t.data
    n = t.n
    δij = i == j ? one(T) : zero(T)
    nnij = n[i] * n[j]
    return a * (δij - nnij) + b * nnij
end

# ── KM ───────────────────────────────────────────────────────────────────────

_mandel(t::TensTI{2}) = tomandel(_static_tensor(t))

# ── Arithmetic — axis-aware ± ────────────────────────────────────────────────
# Scalar ops (-, α*A, A*α, A/α) defined in structured_tens_ops.jl.
#
# Binary ± between two TensTI of the same order:
#   • same axis, same N     → structured result (data-wise ±)
#   • same axis, mixed N    → lift to the richer parametrization (see the
#                              mixed-N methods in the N=8 / N=3 section below
#                              and in structured_tens_promotion.jl)
#   • different axes        → fall back to the generic `Tens` route: the sum
#                              of two TI tensors with different axes has no TI
#                              structure. This replaces the former hard
#                              assertion, enabling multi-axis accumulation in
#                              scheme kernels (e.g. self-consistent estimates
#                              with several inclusion-family axes).
# Note: axes are compared for strict equality; `n` and `−n` are treated as
# different axes (the antisymmetric components ℓ₇, ℓ₈ and the order-2 `c`
# coefficient are odd in `n`), falling back to the generic route.

"""
    generic_tens(t::AbstractTens)

The tensor `t` as a generic `Tens` of its canonical components: a structured
tensor (`TensISO`, `TensTI`, `TensOrtho`, `TensCubic`) forgets its class, so that
it accepts any update — `zero(generic_tens(t))` is a zero of the same order,
dimension and element type that a non-isotropic term can be added to, where
`zero(t)` of a `TensTI` is a `TensISO`.
"""
generic_tens(t::AbstractTens) = Tens(_static_tensor(t))

@inline _generic_binary(op, A::AbstractTens{order, dim}, B::AbstractTens{order, dim}) where {order, dim} =
    Tens(_static_tensor(broadcast(op, get_array(A), get_array(B)), Val(dim), Val(order)))

for OP in (:+, :-)
    @eval @inline function Base.$OP(A::TensTI{order, <:Any, N}, B::TensTI{order, <:Any, N}) where {order, N}
        axis(A) == axis(B) && return _coaxial_binary($OP, A, B)
        c = _axis_cosine(axis(A), axis(B))
        c === nothing && return _structured_sum($OP, A, B)
        abs(abs(c) - one(c)) < 1.0e-10 && return _coaxial_binary($OP, A, _about_axis(B, axis(A), c))
        abs(c) < 1.0e-10 && _has_orthogonal_ortho(A) && return _orthogonal_ti_binary($OP, A, B)
        return _structured_sum($OP, A, B)
    end
end

# Two TI tensors on different axes keep the symmetries the two share: about one
# axis when the axes are parallel (opposite, or equal within 1e-10), the second
# operand being written about the first's axis; orthotropic on the frame
# (n_A, n_B, n_A × n_B) when they are orthogonal and both of the D∞h shapes,
# which have the three mirrors of that frame; nothing else otherwise, the
# generic route. On a symbolic axis only equal axes are recognized: the cosine
# is `nothing`.
@inline _axis_cosine(nA, nB) = is_hard_numeric(real(promote_type(eltype(nA), eltype(nB)))) ?
    nA[1] * nB[1] + nA[2] * nB[2] + nA[3] * nB[3] : nothing

# B written about n, parallel to its own axis: the terms odd in the axis (`c·w`
# of `N = 3`, `ℓ₇`, `ℓ₈` of `N = 8`) change sign when n is opposite.
_about_axis(B::TensTI{order, T, N}, n, c) where {order, T, N} =
    TensTI{order, promote_type(T, eltype(n)), N}(c > 0 ? get_data(B) : _odd_in_axis_flipped(get_data(B)), n)
_odd_in_axis_flipped(d::NTuple{3}) = (d[1], d[2], -d[3])
_odd_in_axis_flipped(d::NTuple{8}) = (d[1], d[2], d[3], d[4], d[5], d[6], -d[7], -d[8])
_odd_in_axis_flipped(d::NTuple) = d

# Orthogonal axes: the D∞h shapes are orthotropic on (n_A, n_B, n_A × n_B), n_B
# made exactly orthogonal to n_A; the SO(2) shapes share no symmetry there.
_has_orthogonal_ortho(::Union{TensTI{2, <:Any, 2}, TensTI{4, <:Any, 5}, TensTI{4, <:Any, 6}}) = true
_has_orthogonal_ortho(::TensTI) = false
@noinline function _orthogonal_ti_binary(op::F, A::TensTI, B::TensTI) where {F}
    nA, nB = SVector(axis(A)), SVector(axis(B))
    e₂ = normalize(nB - (transpose(nA) * nB) * nA)
    ℱ = RotatedBasis(hcat(nA, e₂, cross(nA, e₂)), Val(3), Val(:typed))
    return _same_frame_binary(op, walpole_to_ortho(A, ℱ, 1), walpole_to_ortho(B, ℱ, 2))
end

# The data-wise sum or difference of two TI tensors of the same parametrization
# sharing their axis (not checked): an isotropic operand lifted onto the other's
# axis goes straight here, so that the result is inferred as a `TensTI`.
@inline _coaxial_binary(op::F, A::TensTI{order, <:Any, N}, B::TensTI{order, <:Any, N}) where {F, order, N} =
    _rebuild(A, broadcast(op, get_data(A), get_data(B)))

# ── Inverse (order 2) ───────────────────────────────────────────────────────

"""
    inv(t::TensTI{2,T,2}) → TensTI{2,T,2}

Inverse: `(a·nT + b·nₙ)⁻¹ = (1/a)·nT + (1/b)·nₙ`.

# Examples
```julia
julia> A = TensTI{2}(5.0, 8.0, [0., 0., 1.]);

julia> inv(A).data
(0.2, 0.125)
```
"""
# In the type of 1/a, as `inv` of an integer matrix is a float one.
@inline Base.inv(t::TensTI{2, T, 2}) where {T} = TensTI{2}(inv(t.data[1]), inv(t.data[2]), t.n)
@inline Base.literal_pow(::typeof(^), A::TensTI{2}, ::Val{-1}) = inv(A)

# ── Trace (order 2) ─────────────────────────────────────────────────────────

"""
    tr(t::TensTI{2}) → scalar

Trace: `tr(a·nT + b·nₙ) = 2a + b`.
"""
LinearAlgebra.tr(t::TensTI{2}) = 2 * t.data[1] + t.data[2]

# ── Symmetry ─────────────────────────────────────────────────────────────────

LinearAlgebra.issymmetric(::TensTI{2}) = true
isisotropic(t::TensTI{2}) = t.data[1] == t.data[2]
istransverselyisotropic(::TensTI) = true
isorthotropic(::TensTI) = false

# Symbolic helpers (tsimplify, tsubs, …) defined in structured_tens_ops.jl

# ── Display ──────────────────────────────────────────────────────────────────

function Base.show(io::IO, A::TensTI{2})
    a, b = get_data(A)
    print(io, "(", a, ") nT + (", b, ") nₙ")
    return print(io, "\n  axis n = ", A.n)
end

function pprint(A::TensTI{2})
    a, b = get_data(A)
    println("(", a, ") nT + (", b, ") nₙ")
    return println("  axis n = ", A.n)
end


# ── express / components for TensTI ──────────────────────────────────────
#
# A structured tensor keeps its Walpole/parametric components in the CANONICAL
# frame — `get_basis(::TensTI) === CanonicalBasis` — with the symmetry axis
# stored separately. Re-expressing it in another orthonormal basis therefore
# needs the very same multi-index rotation as any other tensor. Wrapping the
# stored components in the new basis *without* rotating them relabels them, and
# silently returns a DIFFERENT physical tensor: that is only correct when the
# target basis is canonical, or for `TensISO`, which is genuinely invariant.
#
# The rotation itself is not re-derived here: the components are wrapped as a
# plain `Tens` in the canonical basis and handed to the generic
# `TensOrthonormal` path, which is the one already exercised by the rest of the
# package. `express` short-circuits on an equal basis, so the canonical
# target stays allocation-free apart from the wrapper.

_as_canonical_tens(t::AbstractTens{order, 3, T}) where {order, T} =
    Tens(_static_tensor(t), CanonicalBasis{3, T}())

express(t::TensTI{2, T}, ℬ::OrthonormalBasis{3, T}) where {T} =
    express(_as_canonical_tens(t), ℬ)
components(t::TensTI{2, T}, ℬ::OrthonormalBasis{3, T}, v::NTuple{2, Symbol}) where {T} =
    components(_as_canonical_tens(t), ℬ, v)
components(t::TensTI{2}) = get_array(t)
components(t::TensTI{2}, ::NTuple{2, Symbol}) = get_array(t)

# ── otimes specializations (TensTI{2} → TensTI{4}) ───────────────────────

"""
    otimes(A::TensTI{2}) → TensTI{4, T, 5}

Self tensor product of a TI 2nd-order tensor.  The result is always
major-symmetric (ℓ₃ = ℓ₄) and lives in the Walpole basis with N=5.

    (a·nT + b·nₙ) ⊗ (a·nT + b·nₙ)
    = b²W₁ + 2a²W₂ + √2·ab·(W₃+W₄)
"""
function Tensors.otimes(A::TensTI{2, T, 2}) where {T}
    a, b = A.data
    sq2 = sqrt(T(2))
    return TensTI{4, T, 5}((b * b, T(2) * a * a, sq2 * a * b, zero(T), zero(T)), A.n)
end

"""
    otimes(A::TensTI{2}, B::TensTI{2}) → TensTI{4, T, 6}

Tensor product of two TI 2nd-order tensors with the same axis.
Falls back to generic `otimes` if axes differ.

    (a₁·nT + b₁·nₙ) ⊗ (a₂·nT + b₂·nₙ)
    = b₁b₂·W₁ + 2a₁a₂·W₂ + √2·b₁a₂·W₃ + √2·a₁b₂·W₄
"""
function Tensors.otimes(A::TensTI{2, T1, 2}, B::TensTI{2, T2, 2}) where {T1, T2}
    if A.n != B.n
        return invoke(Tensors.otimes, Tuple{AbstractTens{2, 3}, AbstractTens{2, 3}}, A, B)
    end
    T = promote_type(T1, T2)
    a₁, b₁ = A.data
    a₂, b₂ = B.data
    sq2 = sqrt(T(2))
    return TensTI{4, T, 6}(
        (
            T(b₁ * b₂), T(2) * a₁ * a₂, sq2 * T(b₁ * a₂), sq2 * T(a₁ * b₂),
            zero(T), zero(T),
        ), A.n
    )
end

"""
    otimes(A::TensISO{2,3}, B::TensTI{2}) → TensTI{4, T, 6}

Tensor product of a 3D isotropic 2nd-order tensor with a TI 2nd-order tensor.
The isotropic tensor `λ·𝟏` is treated as `TensTI{2}(λ,λ,n)` with the axis of B.
"""
function Tensors.otimes(A::TensISO{2, 3}, B::TensTI{2, T2, 2}) where {T2}
    T = promote_type(eltype(A), T2)
    λ = A.data[1]
    a₂, b₂ = B.data
    sq2 = sqrt(T(2))
    return TensTI{4, T, 6}(
        (
            T(λ * b₂), T(2) * λ * a₂, sq2 * T(λ * a₂), sq2 * T(λ * b₂),
            zero(T), zero(T),
        ), B.n
    )
end

"""
    otimes(A::TensTI{2}, B::TensISO{2,3}) → TensTI{4, T, 6}

Tensor product of a TI 2nd-order tensor with a 3D isotropic 2nd-order tensor.
The isotropic tensor `λ·𝟏` is treated as `TensTI{2}(λ,λ,n)` with the axis of A.
"""
function Tensors.otimes(A::TensTI{2, T1, 2}, B::TensISO{2, 3}) where {T1}
    T = promote_type(T1, eltype(B))
    a₁, b₁ = A.data
    λ = B.data[1]
    sq2 = sqrt(T(2))
    return TensTI{4, T, 6}(
        (
            T(b₁ * λ), T(2) * a₁ * λ, sq2 * T(b₁ * λ), sq2 * T(a₁ * λ),
            zero(T), zero(T),
        ), A.n
    )
end


##############################################################################
# TensOrtho — orthotropic tensors, order 4 (9 or 12 constants) and 2 (3)
##############################################################################
#
# In the material frame (e₁,e₂,e₃) with Pₘ = eₘ⊗eₘ:
#
#   ℂ = C₁₁P₁⊗P₁ + C₂₂P₂⊗P₂ + C₃₃P₃⊗P₃
#     + C₁₂(P₁⊗P₂+P₂⊗P₁) + C₁₃(P₁⊗P₃+P₃⊗P₁) + C₂₃(P₂⊗P₃+P₃⊗P₂)
#     + 2C₄₄(P₂⊠ˢP₃) + 2C₅₅(P₁⊠ˢP₃) + 2C₆₆(P₁⊠ˢP₂)
#
# where C₄₄=C₂₃₂₃, C₅₅=C₁₃₁₃, C₆₆=C₁₂₁₂.
#
# KM in the material frame (Kelvin-Mandel, ordering 11,22,33,23,13,12):
#
#   [[C₁₁,C₁₂,C₁₃, 0,  0,  0 ],
#    [C₁₂,C₂₂,C₂₃, 0,  0,  0 ],
#    [C₁₃,C₂₃,C₃₃, 0,  0,  0 ],
#    [ 0,  0,  0, 2C₄₄, 0,  0 ],
#    [ 0,  0,  0,  0, 2C₅₅, 0 ],
#    [ 0,  0,  0,  0,  0, 2C₆₆]]
# ─────────────────────────────────────────────────────────────────────────────

"""
    TensOrtho{order, T, N, B<:OrthonormalBasis{3}} <: AbstractTens{order, 3, T}

Orthotropic tensor of order `order` (always `dim=3`) with material frame
`(e₁,e₂,e₃)` of type `B`, parametrized like [`TensTI`](@ref): `N` is the number
of stored coefficients, with `Pₘ = eₘ⊗eₘ`.

Three concrete shapes are supported:

| Parametrization           | Role                                   | Invariance group | Stored coefficients                                    |
|---------------------------|----------------------------------------|------------------|--------------------------------------------------------|
| `TensOrtho{2, T, 3, B}`   | 2nd-order orthotropic                  | `D₂h`            | `data = (a₁, a₂, a₃)`, `frame::B`                      |
| `TensOrtho{4, T, 9, B}`   | 4th-order orthotropic, major-symmetric | `D₂h`            | `data = (C₁₁,C₂₂,C₃₃,C₁₂,C₁₃,C₂₃,C₄₄,C₅₅,C₆₆)`, `frame` |
| `TensOrtho{4, T, 12, B}`  | 4th-order orthotropic, general         | `D₂h`            | `data = (C₁₁,…,C₆₆, C₂₁,C₃₁,C₃₂)`, `frame`             |

`D₂h` is generated by the three mirror planes of the frame. At even order its
rotation subgroup `D₂` (the half-turns about the frame axes) has the same
invariant tensors, the inversion acting trivially: unlike transverse isotropy,
orthotropy has no chiral shape, and these three are the whole family. The two
orders correspond: `N=3` ↔ `N=9, 12`.

- **Order 2** (`N=3`): `𝐀 = a₁P₁ + a₂P₂ + a₃P₃`, diagonal in the frame. A
  second-order tensor invariant under the three mirror planes is necessarily of
  this form, hence symmetric: there is no general order-2 counterpart.

- **Order 4**: in the frame,

      ℂ = Σₘₙ Cₘₙ Pₘ⊗Pₙ + 2C₄₄(P₂⊠ˢP₃) + 2C₅₅(P₁⊠ˢP₃) + 2C₆₆(P₁⊠ˢP₂)

  (`C₄₄=C₂₃₂₃`, `C₅₅=C₁₃₁₃`, `C₆₆=C₁₂₁₂`), so that the Kelvin-Mandel matrix in
  the material frame is block-diagonal:

      [[C₁₁,C₁₂,C₁₃, 0,   0,   0  ],
       [C₂₁,C₂₂,C₂₃, 0,   0,   0  ],
       [C₃₁,C₃₂,C₃₃, 0,   0,   0  ],
       [ 0,  0,  0,  2C₄₄, 0,   0  ],
       [ 0,  0,  0,   0,  2C₅₅, 0  ],
       [ 0,  0,  0,   0,   0,  2C₆₆]]

  Major-symmetric tensors have `Cₙₘ = Cₘₙ` and are stored under `N=9`.
  Synthetic notation: `ℂ ≡ (C, C₄₄, C₅₅, C₆₆)` with `C` the 3×3 normal block.
  Product:  `(ℂ⊡𝔻)_C = C D`, `(ℂ⊡𝔻)ₘₘ = 2CₘₘDₘₘ` for the shears;
  inverse: `(ℂ⁻¹)_C = C⁻¹`, `1/(4Cₘₘ)`.

The product of two shapes of the same order returns the smallest shape holding
it (`N=9 ⊡ N=9 → N=12`), the inverse keeps the shape, the major transpose
transposes the normal block, and `Tensors.majorsymmetric` returns `N=9`.

The frame's element type is a free type parameter, independent of the data
eltype `T`: differentiating w.r.t. the coefficients (`T = ForwardDiff.Dual`)
does not require a `Dual`-typed geometric frame. Keeping `B` as a concrete type
parameter (rather than erasing it to the abstract `OrthonormalBasis{3}`) is what
lets `t.frame` be inferred and stored unboxed.

See also [`iso_to_ortho`](@ref), [`walpole_to_ortho`](@ref), [`cubic_to_ortho`](@ref).
"""
struct TensOrtho{order, T, N, B <: OrthonormalBasis{3}} <: AbstractTens{order, 3, T}
    data::NTuple{N, T}    # (a₁,a₂,a₃) or (C₁₁,C₂₂,C₃₃,C₁₂,C₁₃,C₂₃,C₄₄,C₅₅,C₆₆[,C₂₁,C₃₁,C₃₂])
    frame::B               # material frame (e₁,e₂,e₃)
    function TensOrtho{order, T, N, B}(data::NTuple{N, T}, frame::B) where {order, T, N, B <: OrthonormalBasis{3}}
        (order, N) in ((2, 3), (4, 9), (4, 12)) || throw(
            ArgumentError(
                "a TensOrtho of order 2 stores 3 coefficients and one of order 4 stores " *
                    "9 or 12; got order $order with $N"
            )
        )
        return new{order, T, N, B}(data, frame)
    end
end

# `TensOrtho{order, T}(data, frame)` and `TensOrtho{order, T, N}(data, frame)`:
# the number of coefficients and the frame type are read from the arguments, the
# coefficients converted to `T` (the frame keeps its own element type, see above).
TensOrtho{order, T}(data::NTuple{N, Any}, frame::B) where {order, T, N, B <: OrthonormalBasis{3}} =
    TensOrtho{order, T, N, B}(map(x -> convert(T, x), data), frame)
TensOrtho{order, T, N}(data::NTuple{N, Any}, frame::B) where {order, T, N, B <: OrthonormalBasis{3}} =
    TensOrtho{order, T, N, B}(map(x -> convert(T, x), data), frame)

# ── Traits ────────────────────────────────────────────────────────────────────

Base.eltype(::Type{<:TensOrtho{order, T}}) where {order, T} = T
Base.length(::TensOrtho{order}) where {order} = 3^order
Base.size(::TensOrtho{order}) where {order} = ntuple(_ -> 3, Val(order))

get_basis(::TensOrtho{order, T}) where {order, T} = CanonicalBasis{3, T}()
get_var(::TensOrtho{order}) where {order} = ntuple(_ -> :cont, Val(order))
get_var(::TensOrtho, ::Integer) = :cont
get_data(t::TensOrtho) = t.data
"""
    frame(t::TensOrtho)

Return the material frame of an orthotropic tensor.
"""
frame(t::TensOrtho) = t.frame

# ── Rebuild helper (used by symbolic ops) ─────────────────────────────────────
_rebuild(t::TensOrtho{order}, new_data) where {order} = TensOrtho{order, eltype(new_data)}(new_data, frame(t))

# ── Constructors ──────────────────────────────────────────────────────────────

# TensOrtho{4}(C₁₁,C₂₂,C₃₃,C₁₂,C₁₃,C₂₃,C₄₄,C₅₅,C₆₆, frame) → TensOrtho{4, T, 9}
# Major-symmetric orthotropic 4th-order tensor, 9 independent constants.
function TensOrtho{4}(C11, C22, C33, C12, C13, C23, C44, C55, C66, frame::OrthonormalBasis{3})
    T = promote_type(
        typeof(C11), typeof(C22), typeof(C33), typeof(C12), typeof(C13), typeof(C23),
        typeof(C44), typeof(C55), typeof(C66), eltype(frame)
    )
    return TensOrtho{4, T, 9}((C11, C22, C33, C12, C13, C23, C44, C55, C66), frame)
end

# TensOrtho{4}(C₁₁,…,C₆₆, C₂₁,C₃₁,C₃₂, frame) → TensOrtho{4, T, 12}
# General (not necessarily major-symmetric) orthotropic 4th-order tensor.
function TensOrtho{4}(
        C11, C22, C33, C12, C13, C23, C44, C55, C66, C21, C31, C32, frame::OrthonormalBasis{3}
    )
    T = promote_type(
        typeof(C11), typeof(C22), typeof(C33), typeof(C12), typeof(C13), typeof(C23),
        typeof(C44), typeof(C55), typeof(C66), typeof(C21), typeof(C31), typeof(C32), eltype(frame)
    )
    return TensOrtho{4, T, 12}((C11, C22, C33, C12, C13, C23, C44, C55, C66, C21, C31, C32), frame)
end

# TensOrtho{2}(a₁, a₂, a₃, frame) → TensOrtho{2, T, 3}
# `a₁P₁ + a₂P₂ + a₃P₃`, diagonal in the frame.
function TensOrtho{2}(a1, a2, a3, frame::OrthonormalBasis{3})
    T = promote_type(typeof(a1), typeof(a2), typeof(a3), eltype(frame))
    return TensOrtho{2, T, 3}((a1, a2, a3), frame)
end

"""
    TensOrtho(C11,C22,C33,C12,C13,C23,C44,C55,C66, frame)
    TensOrtho(C11,C22,C33,C12,C13,C23,C44,C55,C66,C21,C31,C32, frame)

Shorthand for `TensOrtho{4}(…)`: the orthotropic 4th-order tensor of the 9
constants of a major-symmetric one, or the 12 of a general one (`C₂₁, C₃₁, C₃₂`
last), in the material frame `frame`. The order-2 tensor is `TensOrtho{2}(a₁,
a₂, a₃, frame)`.
"""
TensOrtho(C11, C22, C33, C12, C13, C23, C44, C55, C66, frame::OrthonormalBasis{3}) =
    TensOrtho{4}(C11, C22, C33, C12, C13, C23, C44, C55, C66, frame)
TensOrtho(C11, C22, C33, C12, C13, C23, C44, C55, C66, C21, C31, C32, frame::OrthonormalBasis{3}) =
    TensOrtho{4}(C11, C22, C33, C12, C13, C23, C44, C55, C66, C21, C31, C32, frame)

"""
    TensOrtho(KMmat::AbstractMatrix, frame)
    TensOrtho(KMmat::AbstractMatrix, frame, Val(12))

Build a `TensOrtho` from a 6×6 Kelvin-Mandel matrix expressed in the material frame.
The matrix must have the block-diagonal orthotropic structure:
upper-left 3×3 for normal stresses and lower-right 3×3 diagonal for shear.
The nine-constant form reads the upper triangle of the normal block; with
`Val(12)` the whole block is read, for a tensor without major symmetry.
"""
TensOrtho(KMmat::AbstractMatrix, frame::OrthonormalBasis{3}, ::Val{9}) = TensOrtho(KMmat, frame)
function TensOrtho(KMmat::AbstractMatrix, frame::OrthonormalBasis{3}, ::Val{12})
    O = TensOrtho(KMmat, frame)
    return TensOrtho{4, eltype(O)}((get_data(O)..., KMmat[2, 1], KMmat[3, 1], KMmat[3, 2]), frame)
end
function TensOrtho(KMmat::AbstractMatrix, frame::OrthonormalBasis{3})
    T = eltype(KMmat)
    C11 = KMmat[1, 1]; C22 = KMmat[2, 2]; C33 = KMmat[3, 3]
    C12 = KMmat[1, 2]; C13 = KMmat[1, 3]; C23 = KMmat[2, 3]
    C44 = KMmat[4, 4] / 2
    C55 = KMmat[5, 5] / 2
    C66 = KMmat[6, 6] / 2
    return TensOrtho{4, T}(
        (
            T(C11), T(C22), T(C33), T(C12), T(C13), T(C23),
            T(C44), T(C55), T(C66),
        ), frame
    )
end

# ── get_array ─────────────────────────────────────────────────────────────────

# One entry of the canonical-frame component array.  Factored out of
# `get_array` so that `getindex` can evaluate a SINGLE component without
# materializing all 81 — see the comment on `getindex(::TensOrtho, …)`.
#
# `E` is `vecbasis(frame, :cov)`, whose column `m` is the frame vector `eₘ`;
# `Pₘ[i,j] = E[i,m] * E[j,m]`.  The factor 2 in `2Cₘₘ(Pₘ⊠ˢPₙ + Pₙ⊠ˢPₘ)` is the
# standard Voigt-to-tensor conversion for the shear moduli.
@inline function _ortho_entry(
        E, C11, C22, C33, C12, C13, C23, C44, C55, C66,
        i::Integer, j::Integer, k::Integer, l::Integer
    )
    @inbounds begin
        P(m, a, b) = E[a, m] * E[b, m]
        return (
            C11 * P(1, i, j) * P(1, k, l)
                + C22 * P(2, i, j) * P(2, k, l)
                + C33 * P(3, i, j) * P(3, k, l)
                + C12 * (P(1, i, j) * P(2, k, l) + P(2, i, j) * P(1, k, l))
                + C13 * (P(1, i, j) * P(3, k, l) + P(3, i, j) * P(1, k, l))
                + C23 * (P(2, i, j) * P(3, k, l) + P(3, i, j) * P(2, k, l))
                + C44 * (
                E[i, 2] * E[k, 3] * E[j, 3] * E[l, 2] + E[i, 2] * E[l, 3] * E[j, 3] * E[k, 2] +
                    E[j, 2] * E[k, 3] * E[i, 3] * E[l, 2] + E[j, 2] * E[l, 3] * E[i, 3] * E[k, 2] +
                    E[i, 3] * E[k, 2] * E[j, 2] * E[l, 3] + E[i, 3] * E[l, 2] * E[j, 2] * E[k, 3] +
                    E[j, 3] * E[k, 2] * E[i, 2] * E[l, 3] + E[j, 3] * E[l, 2] * E[i, 2] * E[k, 3]
            ) / 2
                + C55 * (
                E[i, 1] * E[k, 3] * E[j, 3] * E[l, 1] + E[i, 1] * E[l, 3] * E[j, 3] * E[k, 1] +
                    E[j, 1] * E[k, 3] * E[i, 3] * E[l, 1] + E[j, 1] * E[l, 3] * E[i, 3] * E[k, 1] +
                    E[i, 3] * E[k, 1] * E[j, 1] * E[l, 3] + E[i, 3] * E[l, 1] * E[j, 1] * E[k, 3] +
                    E[j, 3] * E[k, 1] * E[i, 1] * E[l, 3] + E[j, 3] * E[l, 1] * E[i, 1] * E[k, 3]
            ) / 2
                + C66 * (
                E[i, 1] * E[k, 2] * E[j, 2] * E[l, 1] + E[i, 1] * E[l, 2] * E[j, 2] * E[k, 1] +
                    E[j, 1] * E[k, 2] * E[i, 2] * E[l, 1] + E[j, 1] * E[l, 2] * E[i, 2] * E[k, 1] +
                    E[i, 2] * E[k, 1] * E[j, 1] * E[l, 2] + E[i, 2] * E[l, 1] * E[j, 1] * E[k, 2] +
                    E[j, 2] * E[k, 1] * E[i, 1] * E[l, 2] + E[j, 2] * E[l, 1] * E[i, 1] * E[k, 2]
            ) / 2
        )
    end
end

# The general form: the major-symmetric entry of the symmetric part of the
# normal block, plus the skew part `aₘₙ (Pₘ⊗Pₙ − Pₙ⊗Pₘ)`, `aₘₙ = (Cₘₙ − Cₙₘ)/2`.
# Written this way the nine-constant entry above is reused as it is, and a
# tensor whose skew part vanishes has exactly the components of its symmetric
# form.
@inline function _ortho_entry(
        E, C11, C22, C33, C12, C13, C23, C44, C55, C66, C21, C31, C32,
        i::Integer, j::Integer, k::Integer, l::Integer
    )
    @inbounds begin
        P(m, a, b) = E[a, m] * E[b, m]
        skew(m, n) = P(m, i, j) * P(n, k, l) - P(n, i, j) * P(m, k, l)
        return _ortho_entry(
            E, C11, C22, C33, (C12 + C21) / 2, (C13 + C31) / 2, (C23 + C32) / 2,
            C44, C55, C66, i, j, k, l
        ) + (C12 - C21) / 2 * skew(1, 2) + (C13 - C31) / 2 * skew(1, 3) + (C23 - C32) / 2 * skew(2, 3)
    end
end

"""
    get_array(t::TensOrtho{4,T}) → Array{T,4}

Compute the 3×3×3×3 component array in the canonical frame.
"""
# Its 36 independent components, through `_static_tensor`: a third faster than
# the 81 evaluated one by one, and exactly minor-symmetric. (The TI forms keep
# their loops, measured faster: `getindex` repeats a setup the loop shares.)
get_array(t::TensOrtho{4}) = Array(_static_tensor(t))

"""
    getindex(t::TensOrtho, i, j, k, l)

Single component in closed form, without building the 81 components:
`TensOrtho <: AbstractArray`, so every generic algorithm (`Array(t)`, `norm`,
broadcasting, `show`, `tomandel`, …) walks it element by element, and a dense
access per element would make an O(81) traversal cost O(81²).
"""
Base.@propagate_inbounds Base.getindex(
    t::TensOrtho{4}, i::Integer, j::Integer, k::Integer, l::Integer
) = _ortho_entry(vecbasis(t.frame, :cov), get_data(t)..., i, j, k, l)

Base.IndexStyle(::Type{<:TensOrtho}) = IndexCartesian()

# ── KM in the material frame ──────────────────────────────────────────────────

"""
    mandel(t::TensOrtho)

Returns the 6×6 Kelvin-Mandel matrix in the **canonical** frame.
Use `mandel_material(t)` for the block-diagonal form in the material frame.
"""
_mandel(t::TensOrtho{4}) = Matrix(_km_in_frame(t))

"""
    mandel_material(t::TensOrtho)

Returns the 6×6 Kelvin-Mandel matrix in the material frame (block-diagonal).
"""
function mandel_material(t::TensOrtho{4, T}) where {T}
    N = _ortho_normal_block(t)
    C44, C55, C66 = _ortho_shears(t)
    z = zero(T)
    return @SMatrix [
        N[1, 1]  N[1, 2]  N[1, 3]   z    z    z  ;
        N[2, 1]  N[2, 2]  N[2, 3]   z    z    z  ;
        N[3, 1]  N[3, 2]  N[3, 3]   z    z    z  ;
        z        z        z      2C44   z    z  ;
        z        z        z        z  2C55   z  ;
        z        z        z        z    z  2C66
    ]
end

# ── The two blocks of the material-frame Kelvin-Mandel matrix ─────────────────
#
# The normal 3×3 block `Cₘₙ` (row: stress, column: strain) and the three shear
# moduli `(C₄₄, C₅₅, C₆₆)`; every closed form below works on them, whatever the
# number of stored constants.
function _ortho_normal_block(t::TensOrtho{4, T, 9}) where {T}
    C11, C22, C33, C12, C13, C23 = get_data(t)
    return SMatrix{3, 3, T}(C11, C12, C13, C12, C22, C23, C13, C23, C33)
end
function _ortho_normal_block(t::TensOrtho{4, T, 12}) where {T}
    C11, C22, C33, C12, C13, C23, _, _, _, C21, C31, C32 = get_data(t)
    return SMatrix{3, 3, T}(C11, C21, C31, C12, C22, C32, C13, C23, C33)
end
_ortho_shears(t::TensOrtho{4}) = (get_data(t)[7], get_data(t)[8], get_data(t)[9])

# The general tensor of normal block `M` and shear moduli `s` on `frame`.
# Its element type spans the frame's as well, as the positional constructor's
# does: the components are products of both.
_ortho_from_blocks(M::SMatrix{3, 3}, s, frame) =
    TensOrtho{4, promote_type(eltype(M), eltype(frame))}((M[1, 1], M[2, 2], M[3, 3], M[1, 2], M[1, 3], M[2, 3], s..., M[2, 1], M[3, 1], M[3, 2]), frame)

# Each branch of a sum on different references returns a concrete type, the
# heavy ones out of line: the union of the result is then split by the caller
# without boxing, where a branch returning a union itself boxed the result of
# the sum on a shared frame too (112 bytes and +60 % measured).

# An orthotropic tensor written on the frame ℱ whose axis π[j] is its own axis j
# up to the sign: an orthotropic tensor is even in each axis, so the signs do
# not enter, and with ρ = π⁻¹ the normal block becomes N[ρ(i), ρ(k)] and the
# shear modulus of the axis i (the one its plane leaves out) that of ρ(i).
@noinline function _ortho_on_frame(B::TensOrtho{4, T, 9}, π, ℱ) where {T}
    ρ, d = invperm(π), get_data(B)
    n(a, b) = _ortho_normal(d, a, b)
    return TensOrtho{4, T}(
        (n(ρ[1], ρ[1]), n(ρ[2], ρ[2]), n(ρ[3], ρ[3]), n(ρ[1], ρ[2]), n(ρ[1], ρ[3]), n(ρ[2], ρ[3]), d[6 + ρ[1]], d[6 + ρ[2]], d[6 + ρ[3]]), ℱ
    )
end
@noinline function _ortho_on_frame(B::TensOrtho{4, T, 12}, π, ℱ) where {T}
    ρ, d = invperm(π), get_data(B)
    n(a, b) = _ortho_normal(d, a, b)
    return TensOrtho{4, T}(
        (
            n(ρ[1], ρ[1]), n(ρ[2], ρ[2]), n(ρ[3], ρ[3]), n(ρ[1], ρ[2]), n(ρ[1], ρ[3]), n(ρ[2], ρ[3]),
            d[6 + ρ[1]], d[6 + ρ[2]], d[6 + ρ[3]], n(ρ[2], ρ[1]), n(ρ[3], ρ[1]), n(ρ[3], ρ[2]),
        ), ℱ
    )
end
_ortho_on_frame(B::TensOrtho{2, T, 3}, π, ℱ) where {T} = (ρ = invperm(π); d = get_data(B); TensOrtho{2, T}((d[ρ[1]], d[ρ[2]], d[ρ[3]]), ℱ))

# The entry (a, b) of the normal block: C₁₂ at 4, C₁₃ at 5, C₂₃ at 6 (a + b + 1),
# the twelve-constant form holding the lower triangle three places further on.
_ortho_normal(d::NTuple{9}, a, b) = a == b ? d[a] : d[a + b + 1]
_ortho_normal(d::NTuple{12}, a, b) = a == b ? d[a] : d[a + b + 1 + (a > b ? 6 : 0)]

# The twelve-constant form of a nine-constant tensor, `Cₙₘ = Cₘₙ` written out.
_ortho_general(t::TensOrtho{4, <:Any, 12}) = t
_ortho_general(t::TensOrtho{4, T, 9}) where {T} =
    (d = get_data(t); TensOrtho{4, T}((d..., d[4], d[5], d[6]), frame(t)))

# ── Kelvin-Mandel congruence of the material frame ────────────────────────────
#
# `mandel(t)` is the Kelvin-Mandel matrix in the CANONICAL frame, `mandel_material(t)`
# the block-diagonal one in the material frame.  The two are related by the
# congruence
#
#     mandel(t) == Q * mandel_material(t) * transpose(Q)
#
# where `Q` is the 6×6 Kelvin-Mandel representation of the 4th-order rotation
# operator `R ⊠ˢ R` (i.e. `mandel(rotation_tensor4(θ,ϕ,ψ))`, see `rotation_tensor4`), built here directly
# from the frame vectors so that no angles need to be recovered.  `Q` is
# orthogonal, which is what makes Kelvin-Mandel the right normalization for
# this: in Voigt notation the analogous matrix is not.
#
# `tomandel` orders the pairs (11, 22, 33, 32, 31, 21), and the off-diagonal
# Kelvin-Mandel basis tensors carry a √2, which is where the three cases below
# come from.  The identity above is asserted in
# `test/test_tens_anisotropic.jl`, on a rotated frame — do not "simplify" the √2
# bookkeeping without re-running it.
const _KM_PAIR = ((1, 1), (2, 2), (3, 3), (3, 2), (3, 1), (2, 1))

@inline function _km_entry(E, s2, I::Int, J::Int)
    (i, j) = _KM_PAIR[I]
    (a, b) = _KM_PAIR[J]
    return @inbounds if i == j && a == b
        E[i, a] * E[j, a]
    elseif i == j
        s2 * E[i, a] * E[j, b]
    elseif a == b
        s2 * E[i, a] * E[j, a]
    else
        E[i, a] * E[j, b] + E[i, b] * E[j, a]
    end
end

@inline function _km_congruence(E::AbstractMatrix{T}) where {T}
    s2 = sqrt(T(2))
    return SMatrix{6, 6, T}(_km_entry(E, s2, I, J) for I in 1:6, J in 1:6)
end

# ── Kelvin-Mandel matrix of a structured tensor, without its 81 components ────
#
# `tomandel(_static_tensor(t))` evaluated all 81 components into an array,
# checked its symmetry, built a `SymmetricTensor` and read 36 entries back out.
#
# `_km_by_entries` evaluates the 36 entries directly, with the scaling and the
# operation order of Tensors.jl's `tomandel` (`s * x` with one off-diagonal
# pair, `s * s * x` with two, `s = √(2one(T))`): the same expressions as before,
# symbolic ones included, 4 to 30 times faster.
#
# `_km_in_frame` goes further for a tensor that carries its material frame
# (`TensOrtho`, `TensCubic`): the congruence of its sparse material-frame matrix,
# about 20 times faster than the 81 components. On a symbolic element type it
# takes the entry route instead, which keeps the expressions a symbolic user
# had (the congruence regroups them into equal but different forms).
@inline function _km_scaled_entry(t, s, I::Int, J::Int)
    (i, j) = _KM_PAIR[I]
    (k, l) = _KM_PAIR[J]
    x = @inbounds t[i, j, k, l]
    return I > 3 && J > 3 ? s * s * x : I > 3 || J > 3 ? s * x : x
end

function _km_by_entries(t::AbstractTens{4, 3, T}) where {T}
    s = √(2one(T))
    return SMatrix{6, 6, T}(_km_scaled_entry(t, s, I, J) for I in 1:6, J in 1:6)
end

_km_in_frame(t) = _km_in_frame(Val(is_hard_numeric(eltype(t))), t)
_km_in_frame(::Val{false}, t) = _km_by_entries(t)
function _km_in_frame(::Val{true}, t)
    E = vecbasis(frame(t), :cov)
    T = promote_type(eltype(t), eltype(E))
    Q = _km_congruence(SMatrix{3, 3, T}(E))
    return Q * mandel_material(t) * transpose(Q)
end

# ── Arithmetic ────────────────────────────────────────────────────────────────
# Scalar ops (-, α*A, A*α, A/α) defined in
# structured_tens_ops.jl

# Two forms of the same length on the same frame add coefficient by coefficient,
# on frames equal up to a signed permutation of the axes once the second is
# written on the first's frame, and take the generic route otherwise; a nine-constant
# operand is written in the twelve-constant form against a general one.
# (One method per shape: a free `N` shared by both arguments would make the
# method less specific than the mixed one, which would then call itself.)
for OP in (:+, :-)
    for (order, N) in ((2, 3), (4, 9), (4, 12))
        @eval @inline function Base.$OP(A::TensOrtho{$order, <:Any, $N}, B::TensOrtho{$order, <:Any, $N})
            frame(A) == frame(B) && return _same_frame_binary($OP, A, B)
            π = _frame_permutation(frame(A), frame(B))
            π === nothing && return _structured_sum($OP, A, B)
            return _same_frame_binary($OP, A, _ortho_on_frame(B, π, frame(A)))
        end
    end
    @eval @inline Base.$OP(A::TensOrtho{4}, B::TensOrtho{4}) = $OP(_ortho_general(A), _ortho_general(B))
end

# ── Inverse ───────────────────────────────────────────────────────────────────

"""
    inv(t::TensOrtho) → TensOrtho

Inverse in closed form: the upper 3×3 block `[[C₁₁,C₁₂,C₁₃],[C₂₁,C₂₂,C₂₃],
[C₃₁,C₃₂,C₃₃]]` is inverted via its scalar adjugate/determinant, and each shear
term inverts as `Cₘₘ' = 1/(4Cₘₘ)` (since the KM diagonal entry is `2Cₘₘ`). No
heap allocation, no LU factorization — matches the block structure documented
on `TensOrtho`. The inverse keeps the form: nine constants (the inverse of a
symmetric block is symmetric) or twelve.
"""
function Base.inv(t::TensOrtho{4, T, 9}) where {T}
    a, b, c, d, e, f = get_data(t)[1:6]   # a=C₁₁, b=C₂₂, c=C₃₃, d=C₁₂, e=C₁₃, f=C₂₃
    C44, C55, C66 = get_data(t)[7:9]
    det = a * b * c + 2 * d * e * f - a * f * f - b * e * e - c * d * d
    inv_det = one(T) / det
    C11_inv = (b * c - f * f) * inv_det
    C22_inv = (a * c - e * e) * inv_det
    C33_inv = (a * b - d * d) * inv_det
    C12_inv = (e * f - c * d) * inv_det
    C13_inv = (d * f - b * e) * inv_det
    C23_inv = (d * e - a * f) * inv_det
    C44_inv = inv(4 * C44)
    C55_inv = inv(4 * C55)
    C66_inv = inv(4 * C66)
    return TensOrtho(
        C11_inv, C22_inv, C33_inv, C12_inv, C13_inv, C23_inv,
        C44_inv, C55_inv, C66_inv, t.frame
    )
end

function Base.inv(t::TensOrtho{4, T, 12}) where {T}
    M = _ortho_normal_block(t)
    # Adjugate over determinant: the cofactors in closed form, no pivoting, so
    # that a symbolic block inverts as a numeric one.
    a, d, g, b, e, h, c, f, k = M   # M = [a b c; d e f; g h k], column-major
    A = e * k - f * h
    B = f * g - d * k
    C = d * h - e * g
    inv_det = one(T) / (a * A + b * B + c * C)
    Minv = SMatrix{3, 3, T}(
        A * inv_det, B * inv_det, C * inv_det,
        (c * h - b * k) * inv_det, (a * k - c * g) * inv_det, (b * g - a * h) * inv_det,
        (b * f - c * e) * inv_det, (c * d - a * f) * inv_det, (a * e - b * d) * inv_det,
    )
    return _ortho_from_blocks(Minv, map(x -> inv(4 * x), _ortho_shears(t)), frame(t))
end

@inline Base.literal_pow(::typeof(^), A::TensOrtho, ::Val{-1}) = inv(A)

# ── Major transpose ───────────────────────────────────────────────────────────
#
# `tᵀ_ijkl = t_klij`: the normal block is transposed, the shear terms are their
# own transposes. A nine-constant tensor is its own transpose.

@inline Base.transpose(A::TensOrtho{4, <:Any, 9}) = A
@inline function Base.transpose(A::TensOrtho{4, <:Any, 12})
    C11, C22, C33, C12, C13, C23, C44, C55, C66, C21, C31, C32 = get_data(A)
    return _rebuild(A, (C11, C22, C33, C21, C31, C32, C44, C55, C66, C12, C13, C23))
end
@inline Base.adjoint(A::TensOrtho) = transpose(conj(A))

# ── Double contraction ────────────────────────────────────────────────────────

"""
    dcontract(A::TensOrtho{4}, B::TensOrtho{4}) → TensOrtho{4, T, 12}

Closed form via the material-frame Kelvin-Mandel blocks.

Double contraction *is* the matrix product in Kelvin-Mandel, and in the shared
material frame both operands are block-diagonal, `[3×3] ⊕ diag(2C₄₄,2C₅₅,2C₆₆)`,
so the product is one 3×3 product plus three scalar products:

    mandel_material(A ⊡ B) == mandel_material(A) * mandel_material(B)

The result is orthotropic on the same frame, in the **twelve-constant** form:
the product of two symmetric 3×3 blocks is not symmetric unless they commute,
so `A ⊡ B` has no major symmetry in general, the same widening that makes
`dcontract(::TensTI{4}, ::TensTI{4})` return `N=6` rather than `N=5`. Take
`Tensors.majorsymmetric` of the result for its nine-constant symmetric part.

If the two material frames differ the product is generally fully anisotropic
and the generic route is taken instead.
"""
# Inlined, so that a caller splits the union on the spot instead of receiving
# a boxed result; the generic route is kept out of line.
@inline function Tensors.dcontract(A::TensOrtho{4}, B::TensOrtho{4})
    frame(A) == frame(B) || return _generic_dcontract(A, B)
    # The diagonal Kelvin-Mandel entries are `2Cₘₘ`: `2Cₐ · 2C_b = 2 (2CₐC_b)`.
    s = map((a, b) -> 2 * a * b, _ortho_shears(A), _ortho_shears(B))
    return _ortho_from_blocks(_ortho_normal_block(A) * _ortho_normal_block(B), s, frame(A))
end

@noinline _generic_dcontract(A, B) = Tensors.dcontract(generic_tens(A), generic_tens(B))

# ── Symmetry ──────────────────────────────────────────────────────────────────

LinearAlgebra.issymmetric(::TensOrtho{4, <:Any, 9}) = true
# The twelve-constant form is major-symmetric when its normal block is symmetric.
function LinearAlgebra.issymmetric(t::TensOrtho{4, <:Any, 12})
    _, _, _, C12, C13, C23, _, _, _, C21, C31, C32 = get_data(t)
    return isequal(C12, C21) && isequal(C13, C31) && isequal(C23, C32)
end
Tensors.isminorsymmetric(::TensOrtho{4}) = true
Tensors.ismajorsymmetric(t::TensOrtho{4}) = issymmetric(t)

# ── isisotropic / istransverselyisotropic / isorthotropic ───────────────────────────────────────────────────

isisotropic(::TensOrtho{4}) = false
istransverselyisotropic(::TensOrtho) = false
isorthotropic(::TensOrtho) = true
isorthotropic(::Any) = false   # universal fallback

# ──────────────────────────────────────────────────────────────────────────────
# Unified symmetry accessors
# ──────────────────────────────────────────────────────────────────────────────

"""
    symmetry(t) -> Symbol

Return the material symmetry class imposed by the container type of `t`:
`:ISO`, `:TI`, `:ORTHO`, or `:ANISO` (default for any unstructured tensor).

This is a *type-level* query — it tells you what symmetry the storage
guarantees, not whether the numerical components happen to satisfy a tighter
symmetry.  For value-level detection use `best_sym_tens(t)`.

# Examples
```julia
julia> symmetry(TensISO{3}(2.0, 3.0))
:ISO

julia> symmetry(ti_from_constants(10., 3., 2.5, 12., 2., [0., 0., 1.]))
:TI
```

See also [`reference`](@ref), [`isisotropic`](@ref), [`istransverselyisotropic`](@ref),
[`isorthotropic`](@ref), [`best_sym_tens`](@ref).
"""
symmetry(::TensISO) = :ISO
symmetry(::TensTI) = :TI
symmetry(::TensOrtho) = :ORTHO
symmetry(::Any) = :ANISO

"""
    reference(t)

Return the geometric reference that parametrises the material symmetry of `t`:
the symmetry axis `NTuple{3}` for a `TensTI`, the material frame
`OrthonormalBasis{3}` for a `TensOrtho`, and `nothing` for `TensISO` or any
tensor without a structured reference.

# Examples
```julia
julia> reference(ti_from_constants(10., 3., 2.5, 12., 2., [0., 0., 1.]))
(0.0, 0.0, 1.0)

julia> reference(TensISO{3}(2.0, 3.0)) === nothing
true
```

See also [`axis`](@ref), [`frame`](@ref), [`symmetry`](@ref).
"""
reference(t::TensTI) = axis(t)
reference(t::TensOrtho) = frame(t)
reference(::TensISO) = nothing
reference(::Any) = nothing

# Symbolic helpers (tsimplify, tsubs, …) defined in structured_tens_ops.jl

# ── Display ───────────────────────────────────────────────────────────────────

function Base.show(io::IO, A::TensOrtho)
    _show_ortho_terms(io, A)
    return print(io, "\n  frame: ", vecbasis(A.frame, :cov))
end

function pprint(A::TensOrtho)
    _show_ortho_terms(stdout, A)
    return println("\n  frame: ", vecbasis(A.frame, :cov))
end

function _show_ortho_terms(io::IO, A::TensOrtho{4, <:Any, 9})
    C11, C22, C33, C12, C13, C23, C44, C55, C66 = get_data(A)
    print(io, "(", C11, ") P₁⊗P₁ + (", C22, ") P₂⊗P₂ + (", C33, ") P₃⊗P₃")
    print(io, "\n  + (", C12, ")(P₁⊗P₂+P₂⊗P₁) + (", C13, ")(P₁⊗P₃+P₃⊗P₁) + (", C23, ")(P₂⊗P₃+P₃⊗P₂)")
    return print(io, "\n  + 2(", C44, ")(P₂⊠ˢP₃) + 2(", C55, ")(P₁⊠ˢP₃) + 2(", C66, ")(P₁⊠ˢP₂)")
end
function _show_ortho_terms(io::IO, A::TensOrtho{4, <:Any, 12})
    C11, C22, C33, C12, C13, C23, C44, C55, C66, C21, C31, C32 = get_data(A)
    print(io, "(", C11, ") P₁⊗P₁ + (", C22, ") P₂⊗P₂ + (", C33, ") P₃⊗P₃")
    print(io, "\n  + (", C12, ") P₁⊗P₂ + (", C21, ") P₂⊗P₁ + (", C13, ") P₁⊗P₃ + (", C31, ") P₃⊗P₁")
    print(io, "\n  + (", C23, ") P₂⊗P₃ + (", C32, ") P₃⊗P₂")
    return print(io, "\n  + 2(", C44, ")(P₂⊠ˢP₃) + 2(", C55, ")(P₁⊠ˢP₃) + 2(", C66, ")(P₁⊠ˢP₂)")
end


# ─────────────────────────────────────────────────────────────────────────────
# TensOrtho{2} — orthotropic 2nd-order tensor, the counterpart of TensTI{2,T,2}
# ─────────────────────────────────────────────────────────────────────────────
#
#   𝐀 = a₁P₁ + a₂P₂ + a₃P₃,   Pₘ = eₘ⊗eₘ,
#
# diagonal in the material frame. Product, inverse and trace act on the three
# coefficients; with an order-4 orthotropic tensor on the same frame the double
# contraction is the normal block acting on (a₁, a₂, a₃), the shear terms
# vanishing on a diagonal tensor.

# Component (i, j) in the frame whose vectors are the columns of `E`.
@inline _ortho2_entry(E, a1, a2, a3, i::Integer, j::Integer) =
    @inbounds a1 * E[i, 1] * E[j, 1] + a2 * E[i, 2] * E[j, 2] + a3 * E[i, 3] * E[j, 3]

Base.@propagate_inbounds Base.getindex(t::TensOrtho{2}, i::Integer, j::Integer) =
    _ortho2_entry(vecbasis(t.frame, :cov), get_data(t)..., i, j)

@inline _static_tensor(t::TensOrtho{2}) = SymmetricTensor{2, 3, eltype(t)}((i, j) -> @inbounds t[i, j])

"""
    get_array(t::TensOrtho{2,T}) → Array{T,2}

The 3×3 component array in the canonical frame, `Σₘ aₘ eₘ⊗eₘ`.
"""
get_array(t::TensOrtho{2}) = Array(_static_tensor(t))

_mandel(t::TensOrtho{2}) = tomandel(_static_tensor(t))

# The same tensor on another frame (geometry only): what the projections rotate.
_rebuild_frame(t::TensOrtho{order}, frame::OrthonormalBasis{3}) where {order} =
    TensOrtho{order, eltype(t)}(get_data(t), frame)

"""
    inv(t::TensOrtho{2,T,3}) → TensOrtho{2,T,3}

Inverse: `(a₁P₁ + a₂P₂ + a₃P₃)⁻¹ = (1/a₁)P₁ + (1/a₂)P₂ + (1/a₃)P₃`.
"""
@inline Base.inv(t::TensOrtho{2}) = _rebuild(t, map(inv, get_data(t)))

"""
    tr(t::TensOrtho{2}) → scalar

Trace: `a₁ + a₂ + a₃`.
"""
LinearAlgebra.tr(t::TensOrtho{2}) = sum(get_data(t))

LinearAlgebra.issymmetric(::TensOrtho{2}) = true
@inline Base.transpose(A::TensOrtho{2}) = A
isisotropic(t::TensOrtho{2}) = ((a1, a2, a3) = get_data(t); a1 == a2 && a2 == a3)

# The components in an orthonormal basis ℬ in closed form: the same
# coefficients with the frame written in ℬ, as for the order-4 tensor.
function express(t::TensOrtho{2, T}, ℬ::OrthonormalBasis{3, T}) where {T}
    F = vecbasis(frame(t), :cov)
    Tf = promote_type(T, eltype(F))
    E = transpose(SMatrix{3, 3, Tf}(vecbasis(ℬ, :cov))) * SMatrix{3, 3, Tf}(F)
    S = SymmetricTensor{2, 3}((i, j) -> _ortho2_entry(E, get_data(t)..., i, j))
    return Tens(S, ℬ, get_var(t))
end
components(t::TensOrtho{2, T}, ℬ::OrthonormalBasis{3, T}, v::NTuple{2, Symbol}) where {T} =
    components(_as_canonical_tens(t), ℬ, v)

function _show_ortho_terms(io::IO, A::TensOrtho{2})
    a1, a2, a3 = get_data(A)
    return print(io, "(", a1, ") P₁ + (", a2, ") P₂ + (", a3, ") P₃")
end

##############################################################################
# Shared express / components for TensTI{4} and TensOrtho
# (both are 3D order-4 tensors stored in the canonical frame)
##############################################################################

# See the note above `_as_canonical_tens`: the components must be rotated, not
# relabeled. They are, in closed form: in an orthonormal basis ℬ the tensor has
# the same coefficients with its axis written in ℬ, so its components there come
# from the Walpole formula itself, exact for every element type, instead of the
# 81 canonical ones rotated index by index. The eight-coefficient form keeps the
# generic route: its ℓ₇ and ℓ₈ terms change sign with the orientation of ℬ.
express(t::TensTI{4, T}, ℬ::OrthonormalBasis{3, T}) where {T} =
    Tens(_static_tensor(_with_axis_in(t, ℬ)), ℬ, get_var(t))
express(t::TensTI{4, T, 8}, ℬ::OrthonormalBasis{3, T}) where {T} =
    express(_as_canonical_tens(t), ℬ)

# The same tensor written about the components of its axis in ℬ.
_with_axis_in(t::TensTI{4, T, N}, ℬ::OrthonormalBasis{3, T}) where {T, N} =
    TensTI{4, T, N}(get_data(t), Tuple(transpose(SMatrix{3, 3, T}(vecbasis(ℬ, :cov))) * SVector(axis(t))))
components(t::TensTI{4, T}, ℬ::OrthonormalBasis{3, T}, v::NTuple{4, Symbol}) where {T} =
    components(_as_canonical_tens(t), ℬ, v)
components(t::TensTI{4}) = get_array(t)
components(t::TensTI{4}, ::NTuple{4, Symbol}) = get_array(t)

# TensOrtho — same reasoning: stored in the canonical frame, so a rotated basis
# needs a genuine change of components (`express` in `tens_cubic.jl`, shared
# with `TensCubic`).
components(t::TensOrtho{4, T}, ℬ::OrthonormalBasis{3, T}, v::NTuple{4, Symbol}) where {T} =
    components(_as_canonical_tens(t), ℬ, v)
components(t::TensOrtho{order}) where {order} = get_array(t)
components(t::TensOrtho{order}, ::NTuple{order, Symbol}) where {order} = get_array(t)

##############################################################################
# TensTI{4,T,8} — full axially-invariant (azimuthal-average) 4th-order space #
# TensTI{2,T,3} — full axially-invariant 2nd-order space                     #
##############################################################################
#
# The space of minor-symmetric 4th-order tensors invariant under all rotations
# about an axis n is EIGHT-dimensional — it is the commutant of the SO(2)
# action on the 6-dim Kelvin-Mandel space, which decomposes into
#
#   m=0 (invariants)      : {ε_nn-axial, ε-in-plane-spherical} → full 2×2 block
#                            (4 parameters: ℓ₁, ℓ₂, ℓ₃, ℓ₄ — Walpole W₁..W₄)
#   m=1 (axial shears)    : commutant ≅ ℂ → z₁ = ℓ₆ + i ℓ₇
#   m=2 (in-plane devia.) : commutant ≅ ℂ → z₂ = ℓ₅ + i ℓ₈
#
# The two extra generators beyond the classical Walpole basis are the
# antisymmetric (major-antisymmetric) couplings that appear e.g. in the exact
# azimuthal average of strain-concentration tensors:
#
#   W₇[i,j,k,l] = −(1/2)( w[i,k]nₙ[j,l] + w[i,l]nₙ[j,k]
#                        + w[j,k]nₙ[i,l] + w[j,l]nₙ[i,k] )
#   W₈[i,j,k,l] = +(1/4)( w[i,k]nT[j,l] + w[i,l]nT[j,k]
#                        + w[j,k]nT[i,l] + w[j,l]nT[i,k] )
#
# where w is the in-plane rotation generator w·p = n × p (w[i,j] = ε[i,k,j]n[k],
# odd in n).  In the Kelvin-Mandel frame with n = e₃ (ordering 11,22,33,23,13,12):
#
#   W₇ : M₄₅ = −1, M₅₄ = +1                       (m=1 antisymmetric coupling)
#   W₈ : M₆₁ = +1/√2, M₆₂ = −1/√2, M₁₆ = −1/√2, M₂₆ = +1/√2   (m=2)
#
# Because the 8-dim space is a commutant ALGEBRA, it is closed under double
# contraction and inversion, with the cheap product rule
#
#   block 2×2 :  [[ℓ₁,ℓ₃],[ℓ₄,ℓ₂]]  → matrix product / matrix inverse
#   z₁ = ℓ₆ + i ℓ₇                   → complex product / complex inverse
#   z₂ = ℓ₅ + i ℓ₈                   → complex product / complex inverse
#
# W₇ and W₈ annihilate, under double contraction, exactly those 2nd-order
# tensors that are themselves axially invariant about n — i.e. of the form
# a·nT + b·nₙ. (They do NOT annihilate a general symmetric 2nd-order tensor:
# W₇ maps the (23,13) axial-shear pair onto itself rotated by a quarter turn,
# and W₈ does the same for the (11-22, 12) in-plane deviatoric pair. Both
# results are symmetric and generally nonzero.)
#
# That weaker statement is all the ℓ₁..ℓ₄-based rules need, since those concern
# contraction with nT and nₙ; they therefore remain valid for N=8. Contraction
# with a *general* symmetric 2nd-order tensor does pick up ℓ₇ and ℓ₈, and the
# N=8 `dcontract` methods account for them — verified against the generic dense
# route in `test/test_tens_anisotropic.jl`.
#
# Similarly, the space of 2nd-order tensors invariant under rotations about n
# is THREE-dimensional: a·nT + b·nₙ + c·w (the antisymmetric in-plane part c·w
# is what a plain symmetric TI parametrization cannot represent).
# ─────────────────────────────────────────────────────────────────────────────

# ── Constructors ─────────────────────────────────────────────────────────────

# TensTI{4}(ℓ₁,…,ℓ₈, n) → TensTI{4, T, 8}
# Full axially-invariant 4th-order tensor (see block comment above).
function TensTI{4}(ℓ₁, ℓ₂, ℓ₃, ℓ₄, ℓ₅, ℓ₆, ℓ₇, ℓ₈, n)
    T = sqrt2_eltype(
        promote_type(
            typeof(ℓ₁), typeof(ℓ₂), typeof(ℓ₃), typeof(ℓ₄),
            typeof(ℓ₅), typeof(ℓ₆), typeof(ℓ₇), typeof(ℓ₈), eltype_of(n)
        )
    )
    nv = _extract_vec(n)
    return TensTI{4, T, 8}(
        (T(ℓ₁), T(ℓ₂), T(ℓ₃), T(ℓ₄), T(ℓ₅), T(ℓ₆), T(ℓ₇), T(ℓ₈)),
        (T(nv[1]), T(nv[2]), T(nv[3]))
    )
end

# TensTI{2}(a, b, c, n) → TensTI{2, T, 3}
# Full axially-invariant 2nd-order tensor a·nT + b·nₙ + c·w.
function TensTI{2}(a, b, c, n)
    T = promote_type(typeof(a), typeof(b), typeof(c), eltype_of(n))
    nv = _extract_vec(n)
    return TensTI{2, T, 3}((T(a), T(b), T(c)), (T(nv[1]), T(nv[2]), T(nv[3])))
end

"""
    tens_W7(n) → TensTI{4, T, 8}   (m=1 antisymmetric generator)
"""
tens_W7(n) = TensTI{4}(
    zero(eltype_of(n)), zero(eltype_of(n)), zero(eltype_of(n)), zero(eltype_of(n)),
    zero(eltype_of(n)), zero(eltype_of(n)), one(eltype_of(n)), zero(eltype_of(n)), n
)

"""
    tens_W8(n) → TensTI{4, T, 8}   (m=2 antisymmetric generator)
"""
tens_W8(n) = TensTI{4}(
    zero(eltype_of(n)), zero(eltype_of(n)), zero(eltype_of(n)), zero(eltype_of(n)),
    zero(eltype_of(n)), zero(eltype_of(n)), zero(eltype_of(n)), one(eltype_of(n)), n
)

# ── Accessors ────────────────────────────────────────────────────────────────

# For N=8, `walpole_coeffs` returns the six Walpole coefficients (dropping ℓ₇, ℓ₈).
# This keeps every ℓ₁..ℓ₄-based 4th⊡2nd contraction rule valid (W₅..W₈
# annihilate symmetric 2nd-order tensors); use `walpole_coeffs(t, Val(8))` when the
# antisymmetric couplings matter.
walpole_coeffs(t::TensTI{4, T, 8}) where {T} =
    (t.data[1], t.data[2], t.data[3], t.data[4], t.data[5], t.data[6])

"""
    walpole_coeffs(t::TensTI{4,T,N}, ::Val{8}) → NTuple{8,T}

Always returns the 8-tuple `(ℓ₁, …, ℓ₆, ℓ₇, ℓ₈)` of coefficients in the full
axially-invariant basis `{W₁,…,W₈}`.  For `N=5`/`N=6` the antisymmetric
couplings `ℓ₇ = ℓ₈ = 0`.
"""
walpole_coeffs(t::TensTI{4, T, 5}, ::Val{8}) where {T} =
    (t.data[1], t.data[2], t.data[3], t.data[3], t.data[4], t.data[5], zero(T), zero(T))
walpole_coeffs(t::TensTI{4, T, 6}, ::Val{8}) where {T} =
    (t.data[1], t.data[2], t.data[3], t.data[4], t.data[5], t.data[6], zero(T), zero(T))
walpole_coeffs(t::TensTI{4, T, 8}, ::Val{8}) where {T} = t.data

"""
    _lift_walpole_N8(A::TensTI{4, T, N}) → TensTI{4, T, 8}

Lift a Walpole tensor (`N=5` or `N=6`) to the full axially-invariant `N=8`
form with vanishing antisymmetric couplings.
"""
_lift_walpole_N8(A::TensTI{4, T, 5}) where {T} = TensTI{4, T, 8}(walpole_coeffs(A, Val(8)), axis(A))
_lift_walpole_N8(A::TensTI{4, T, 6}) where {T} = TensTI{4, T, 8}(walpole_coeffs(A, Val(8)), axis(A))

# ── Major transpose ───────────────────────────────────────────────────────────
#
# `tᵀ_ijkl = t_klij` (order 4), `tᵀ_ij = t_ji` (order 2), shape kept: in the
# Walpole basis the transpose swaps ℓ₃ and ℓ₄ and changes the sign of the
# major-antisymmetric ℓ₇, ℓ₈; at order 2 it changes the sign of the skew c·w.
# The symmetric shapes are their own transposes.

@inline Base.transpose(A::Union{TensTI{2, <:Any, 2}, TensTI{4, <:Any, 5}}) = A
@inline function Base.transpose(A::TensTI{4, <:Any, 6})
    ℓ₁, ℓ₂, ℓ₃, ℓ₄, ℓ₅, ℓ₆ = get_data(A)
    return _rebuild(A, (ℓ₁, ℓ₂, ℓ₄, ℓ₃, ℓ₅, ℓ₆))
end
@inline function Base.transpose(A::TensTI{4, <:Any, 8})
    ℓ₁, ℓ₂, ℓ₃, ℓ₄, ℓ₅, ℓ₆, ℓ₇, ℓ₈ = get_data(A)
    return _rebuild(A, (ℓ₁, ℓ₂, ℓ₄, ℓ₃, ℓ₅, ℓ₆, -ℓ₇, -ℓ₈))
end
@inline function Base.transpose(A::TensTI{2, <:Any, 3})
    a, b, c = get_data(A)
    return _rebuild(A, (a, b, -c))
end
@inline Base.adjoint(A::TensTI) = transpose(conj(A))
_lift_walpole_N8(A::TensTI{4, T, 8}) where {T} = A

# Lift a 2nd-order symmetric TI (N=2) to the full axially-invariant N=3 form.
_lift_ti2_N3(A::TensTI{2, T, 2}) where {T} =
    TensTI{2, T, 3}((A.data[1], A.data[2], zero(T)), axis(A))
_lift_ti2_N3(A::TensTI{2, T, 3}) where {T} = A

# ── get_array ────────────────────────────────────────────────────────────────

# `[n]ₓ`, the matrix of `p ↦ n × p`: the in-plane rotation generator of a TI
# tensor about `n` (`transpose`-free, so a symbolic axis is not conjugated).
@inline function _cross_matrix(n::NTuple{3, T}) where {T}
    n₁, n₂, n₃ = n
    return @SMatrix [zero(T) -n₃ n₂; n₃ zero(T) -n₁; -n₂ n₁ zero(T)]
end

function get_array(t::TensTI{4, T, 8}) where {T}
    ℓ₁, ℓ₂, ℓ₃, ℓ₄, ℓ₅, ℓ₆, ℓ₇, ℓ₈ = t.data
    n = t.n
    sq2 = sqrt(T(2))
    δ(i, j) = i == j ? one(T) : zero(T)
    nn(i, j) = n[i] * n[j]
    nT(i, j) = δ(i, j) - nn(i, j)
    # In-plane rotation generator w·p = n × p  (odd in n)
    W = _cross_matrix(n)
    w(i, j) = W[i, j]
    result = Array{T, 4}(undef, 3, 3, 3, 3)
    # Fill the canonical (i ≤ j, k ≤ l) entries only, then mirror — this keeps
    # the minor symmetry EXACT in floating point (summation order would
    # otherwise differ between mirrored entries, breaking the
    # `tensor_or_array` SymmetricTensor detection and the 6×6 Mandel routes).
    for i in 1:3, j in i:3, k in 1:3, l in k:3
        W1 = nn(i, j) * nn(k, l)
        W2 = nT(i, j) * nT(k, l) / 2
        W3 = nn(i, j) * nT(k, l) / sq2
        W4 = nT(i, j) * nn(k, l) / sq2
        W5 = (nT(i, k) * nT(j, l) + nT(i, l) * nT(j, k)) / 2 - nT(i, j) * nT(k, l) / 2
        W6 = (nT(i, k) * nn(j, l) + nT(i, l) * nn(j, k) + nn(i, k) * nT(j, l) + nn(i, l) * nT(j, k)) / 2
        W7 = -(w(i, k) * nn(j, l) + w(i, l) * nn(j, k) + w(j, k) * nn(i, l) + w(j, l) * nn(i, k)) / 2
        W8 = (w(i, k) * nT(j, l) + w(i, l) * nT(j, k) + w(j, k) * nT(i, l) + w(j, l) * nT(i, k)) / 4
        val =
            ℓ₁ * W1 + ℓ₂ * W2 + ℓ₃ * W3 + ℓ₄ * W4 +
            ℓ₅ * W5 + ℓ₆ * W6 + ℓ₇ * W7 + ℓ₈ * W8
        result[i, j, k, l] = val
        result[j, i, k, l] = val
        result[i, j, l, k] = val
        result[j, i, l, k] = val
    end
    return result
end

# Canonicalizes (i,j) and (k,l) before evaluating — mirrors `get_array`'s own
# canonical-fill-then-mirror strategy, so `getindex` matches it bit-for-bit
# (the W7/W8 sums are not literally invariant to summation order under an
# (i,j) or (k,l) swap, unlike the simple products in W1-W4).
function Base.getindex(t::TensTI{4, T, 8}, i::Integer, j::Integer, k::Integer, l::Integer) where {T}
    ii, jj = i <= j ? (i, j) : (j, i)
    kk, ll = k <= l ? (k, l) : (l, k)
    ℓ₁, ℓ₂, ℓ₃, ℓ₄, ℓ₅, ℓ₆, ℓ₇, ℓ₈ = t.data
    n = t.n
    sq2 = sqrt(T(2))
    δ(a, b) = a == b ? one(T) : zero(T)
    nn(a, b) = n[a] * n[b]
    nT(a, b) = δ(a, b) - nn(a, b)
    W = _cross_matrix(n)
    w(a, b) = W[a, b]
    W1 = nn(ii, jj) * nn(kk, ll)
    W2 = nT(ii, jj) * nT(kk, ll) / 2
    W3 = nn(ii, jj) * nT(kk, ll) / sq2
    W4 = nT(ii, jj) * nn(kk, ll) / sq2
    W5 = (nT(ii, kk) * nT(jj, ll) + nT(ii, ll) * nT(jj, kk)) / 2 - nT(ii, jj) * nT(kk, ll) / 2
    W6 = (nT(ii, kk) * nn(jj, ll) + nT(ii, ll) * nn(jj, kk) + nn(ii, kk) * nT(jj, ll) + nn(ii, ll) * nT(jj, kk)) / 2
    W7 = -(w(ii, kk) * nn(jj, ll) + w(ii, ll) * nn(jj, kk) + w(jj, kk) * nn(ii, ll) + w(jj, ll) * nn(ii, kk)) / 2
    W8 = (w(ii, kk) * nT(jj, ll) + w(ii, ll) * nT(jj, kk) + w(jj, kk) * nT(ii, ll) + w(jj, ll) * nT(ii, kk)) / 4
    return ℓ₁ * W1 + ℓ₂ * W2 + ℓ₃ * W3 + ℓ₄ * W4 + ℓ₅ * W5 + ℓ₆ * W6 + ℓ₇ * W7 + ℓ₈ * W8
end

function get_array(t::TensTI{2, T, 3}) where {T}
    a, b, c = t.data
    n = SVector(t.n)
    nn = n * transpose(n)
    return copyto!(Matrix{T}(undef, 3, 3), a * (I - nn) + b * nn + c * _cross_matrix(t.n))
end

function Base.getindex(t::TensTI{2, T, 3}, i::Integer, j::Integer) where {T}
    a, b, c = t.data
    n = t.n
    δij = i == j ? one(T) : zero(T)
    wij = _cross_matrix(n)[i, j]
    nnij = n[i] * n[j]
    return a * (δij - nnij) + b * nnij + c * wij
end

# ── Symmetry queries ─────────────────────────────────────────────────────────

LinearAlgebra.issymmetric(t::TensTI{4, T, 8}) where {T} =
    isequal(t.data[3], t.data[4]) && iszero(t.data[7]) && iszero(t.data[8])
Tensors.ismajorsymmetric(t::TensTI{4, T, 8}) where {T} = issymmetric(t)
LinearAlgebra.issymmetric(t::TensTI{2, T, 3}) where {T} = iszero(t.data[3])

LinearAlgebra.tr(t::TensTI{2, T, 3}) where {T} = 2 * t.data[1] + t.data[2]

# ── Double contraction (commutant-algebra product rule) ──────────────────────

"""
    dcontract(A::TensTI{4,T,8}, B::TensTI{4,T,8}) → TensTI{4, T, 8}

Closed product rule in the 8-dim commutant algebra:
2×2 block product for (ℓ₁..ℓ₄), complex products for z₁ = ℓ₆ + iℓ₇ (m=1)
and z₂ = ℓ₅ + iℓ₈ (m=2).
"""
function Tensors.dcontract(A::TensTI{4, <:Any, 8}, B::TensTI{4, <:Any, 8})
    A.n == B.n || return Tensors.dcontract(generic_tens(A), generic_tens(B))
    return _dcontract_coaxial(A, B)
end

function _dcontract_coaxial(A::TensTI{4, <:Any, 8}, B::TensTI{4, <:Any, 8})
    ℓA₁, ℓA₂, ℓA₃, ℓA₄, ℓA₅, ℓA₆, ℓA₇, ℓA₈ = A.data
    ℓB₁, ℓB₂, ℓB₃, ℓB₄, ℓB₅, ℓB₆, ℓB₇, ℓB₈ = B.data
    n₁ = ℓA₁ * ℓB₁ + ℓA₃ * ℓB₄
    n₃ = ℓA₁ * ℓB₃ + ℓA₃ * ℓB₂
    n₄ = ℓA₄ * ℓB₁ + ℓA₂ * ℓB₄
    n₂ = ℓA₄ * ℓB₃ + ℓA₂ * ℓB₂
    n₆ = ℓA₆ * ℓB₆ - ℓA₇ * ℓB₇
    n₇ = ℓA₆ * ℓB₇ + ℓA₇ * ℓB₆
    n₅ = ℓA₅ * ℓB₅ - ℓA₈ * ℓB₈
    n₈ = ℓA₅ * ℓB₈ + ℓA₈ * ℓB₅
    T = promote_type(eltype(A), eltype(B))
    return TensTI{4, T, 8}(
        (T(n₁), T(n₂), T(n₃), T(n₄), T(n₅), T(n₆), T(n₇), T(n₈)), A.n
    )
end

Tensors.dcontract(A::TensTI{4, <:Any, 8}, B::TensTI{4}) =
    Tensors.dcontract(A, _lift_walpole_N8(B))
Tensors.dcontract(A::TensTI{4}, B::TensTI{4, <:Any, 8}) =
    Tensors.dcontract(_lift_walpole_N8(A), B)
_dcontract_coaxial(A::TensTI{4, <:Any, 8}, B::TensTI{4}) = _dcontract_coaxial(A, _lift_walpole_N8(B))
_dcontract_coaxial(A::TensTI{4}, B::TensTI{4, <:Any, 8}) = _dcontract_coaxial(_lift_walpole_N8(A), B)

# ── Inverse ──────────────────────────────────────────────────────────────────

"""
    inv(t::TensTI{4, T, 8}) → TensTI{4, T, 8}

Inverse in the commutant algebra: 2×2 block inverse + two complex inverses.
"""
function Base.inv(t::TensTI{4, T, 8}) where {T}
    ℓ₁, ℓ₂, ℓ₃, ℓ₄, ℓ₅, ℓ₆, ℓ₇, ℓ₈ = t.data
    det = ℓ₁ * ℓ₂ - ℓ₃ * ℓ₄
    d₁ = ℓ₆ * ℓ₆ + ℓ₇ * ℓ₇     # |z₁|²
    d₂ = ℓ₅ * ℓ₅ + ℓ₈ * ℓ₈     # |z₂|²
    return TensTI{4, T, 8}(
        (
            ℓ₂ / det, ℓ₁ / det, -ℓ₃ / det, -ℓ₄ / det,
            ℓ₅ / d₂, ℓ₆ / d₁, -ℓ₇ / d₁, -ℓ₈ / d₂,
        ), t.n
    )
end

@inline Base.literal_pow(::typeof(^), A::TensTI{4, <:Any, 8}, ::Val{-1}) = inv(A)

"""
    inv(t::TensTI{2, T, 3}) → TensTI{2, T, 3}

Inverse of `a·nT + b·nₙ + c·w`: the in-plane part is the complex number
`a + ic` (since `w² = −nT`), the axial part is `b`:
`(a + ic)⁻¹ ⊕ b⁻¹` → `(a/(a²+c²), 1/b, −c/(a²+c²))`.
"""
function Base.inv(t::TensTI{2, T, 3}) where {T}
    a, b, c = t.data
    d = a * a + c * c
    return TensTI{2}(a / d, inv(b), -c / d, t.n)
end

@inline Base.literal_pow(::typeof(^), A::TensTI{2, <:Any, 3}, ::Val{-1}) = inv(A)

# ── Mixed-N ± (same order, lift to the richer parametrization) ───────────────

for OP in (:+, :-)
    @eval Base.$OP(A::TensTI{4, <:Any, 8}, B::TensTI{4, <:Any, 5}) =
        $OP(A, _lift_walpole_N8(B))
    @eval Base.$OP(A::TensTI{4, <:Any, 5}, B::TensTI{4, <:Any, 8}) =
        $OP(_lift_walpole_N8(A), B)
    @eval Base.$OP(A::TensTI{4, <:Any, 8}, B::TensTI{4, <:Any, 6}) =
        $OP(A, _lift_walpole_N8(B))
    @eval Base.$OP(A::TensTI{4, <:Any, 6}, B::TensTI{4, <:Any, 8}) =
        $OP(_lift_walpole_N8(A), B)
    @eval Base.$OP(A::TensTI{2, <:Any, 3}, B::TensTI{2, <:Any, 2}) =
        $OP(A, _lift_ti2_N3(B))
    @eval Base.$OP(A::TensTI{2, <:Any, 2}, B::TensTI{2, <:Any, 3}) =
        $OP(_lift_ti2_N3(A), B)
    @eval function Base.$OP(A::TensISO{4, 3}, B::TensTI{4, <:Any, 8})
        return _coaxial_binary($OP, _lift_walpole_N8(iso_to_ti(A, axis(B))), B)
    end
    @eval function Base.$OP(A::TensTI{4, <:Any, 8}, B::TensISO{4, 3})
        return _coaxial_binary($OP, A, _lift_walpole_N8(iso_to_ti(B, axis(A))))
    end
    @eval function Base.$OP(A::TensISO{2, 3}, B::TensTI{2, <:Any, 3})
        λ = get_data(A)[1]
        a, b, c = get_data(B)
        return TensTI{2}($OP(λ, a), $OP(λ, b), $OP(zero(λ), c), axis(B))
    end
    @eval function Base.$OP(A::TensTI{2, <:Any, 3}, B::TensISO{2, 3})
        a, b, c = get_data(A)
        λ = get_data(B)[1]
        return TensTI{2}($OP(a, λ), $OP(b, λ), c, axis(A))
    end
end

# ── 2nd-order products (dot) in the N=3 space ────────────────────────────────
#
# With nT, nₙ orthogonal projectors and w the in-plane rotation generator:
#   nT·nT = nT, nₙ·nₙ = nₙ, nT·nₙ = 0, w·nT = nT·w = w, w·nₙ = nₙ·w = 0,
#   w·w = −nT
# so the in-plane part (a, c) multiplies like the complex number a + ic and
# the axial part like the scalar b.

function LinearAlgebra.dot(A::TensTI{2, <:Any, 3}, B::TensTI{2, <:Any, 3})
    axis(A) == axis(B) || return LinearAlgebra.dot(generic_tens(A), generic_tens(B))
    a₁, b₁, c₁ = get_data(A)
    a₂, b₂, c₂ = get_data(B)
    return TensTI{2}(a₁ * a₂ - c₁ * c₂, b₁ * b₂, a₁ * c₂ + c₁ * a₂, axis(A))
end

LinearAlgebra.dot(A::TensTI{2, <:Any, 3}, B::TensTI{2, <:Any, 2}) =
    LinearAlgebra.dot(A, _lift_ti2_N3(B))
LinearAlgebra.dot(A::TensTI{2, <:Any, 2}, B::TensTI{2, <:Any, 3}) =
    LinearAlgebra.dot(_lift_ti2_N3(A), B)

function LinearAlgebra.dot(A::TensTI{2, <:Any, 3}, B::TensISO{2, 3})
    a, b, c = get_data(A)
    λ = get_data(B)[1]
    return TensTI{2}(a * λ, b * λ, c * λ, axis(A))
end
function LinearAlgebra.dot(A::TensISO{2, 3}, B::TensTI{2, <:Any, 3})
    λ = get_data(A)[1]
    a, b, c = get_data(B)
    return TensTI{2}(λ * a, λ * b, λ * c, axis(B))
end

# ── 4th ⊡ 2nd with the N=3 space ─────────────────────────────────────────────
# A minor-symmetric 4th-order tensor annihilates the antisymmetric part c·w,
# so the contraction reduces to the symmetric (a, b) part.

Tensors.dcontract(A::TensTI{4}, B::TensTI{2, <:Any, 3}) =
    Tensors.dcontract(A, TensTI{2}(B.data[1], B.data[2], axis(B)))
Tensors.dcontract(A::TensTI{2, <:Any, 3}, B::TensTI{4}) =
    Tensors.dcontract(TensTI{2}(A.data[1], A.data[2], axis(A)), B)

# ── isisotropic for the new shapes ────────────────────────────────────────────────

isisotropic(t::TensTI{2, <:Any, 3}) = t.data[1] == t.data[2] && iszero(t.data[3])

# ── Display ──────────────────────────────────────────────────────────────────

function Base.show(io::IO, A::TensTI{4, <:Any, 8})
    ℓ₁, ℓ₂, ℓ₃, ℓ₄, ℓ₅, ℓ₆, ℓ₇, ℓ₈ = A.data
    print(
        io, "(", ℓ₁, ") W₁ + (", ℓ₂, ") W₂ + (", ℓ₃, ") W₃ + (", ℓ₄,
        ") W₄ + (", ℓ₅, ") W₅ + (", ℓ₆, ") W₆ + (", ℓ₇, ") W₇ + (", ℓ₈, ") W₈"
    )
    return print(io, "\n  axis n = ", A.n)
end

function pprint(A::TensTI{4, <:Any, 8})
    ℓ₁, ℓ₂, ℓ₃, ℓ₄, ℓ₅, ℓ₆, ℓ₇, ℓ₈ = A.data
    println(
        "(", ℓ₁, ") W₁ + (", ℓ₂, ") W₂ + (", ℓ₃, ") W₃ + (", ℓ₄,
        ") W₄ + (", ℓ₅, ") W₅ + (", ℓ₆, ") W₆ + (", ℓ₇, ") W₇ + (", ℓ₈, ") W₈"
    )
    return println("  axis n = ", A.n)
end

function Base.show(io::IO, A::TensTI{2, <:Any, 3})
    a, b, c = get_data(A)
    print(io, "(", a, ") nT + (", b, ") nₙ + (", c, ") w")
    return print(io, "\n  axis n = ", A.n)
end

function pprint(A::TensTI{2, <:Any, 3})
    a, b, c = get_data(A)
    println("(", a, ") nT + (", b, ") nₙ + (", c, ") w")
    return println("  axis n = ", A.n)
end

##############################################################################
# Exports
##############################################################################

export TensTI, TensOrtho
export tens_W1, tens_W2, tens_W3, tens_W4, tens_W5, tens_W6, tens_W7, tens_W8
export walpole_basis, walpole_basis_sym
export walpole_coeffs, axis, frame, reference, symmetry
export iso_to_ti, istransverselyisotropic, isorthotropic
export ti_from_constants, ti_constants
export mandel_material
