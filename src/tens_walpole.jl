##############################################################################
# TensTI{4} — transversely isotropic 4th-order tensors (Walpole basis)    #
# TensOrtho  — orthotropic 4th-order tensors                                #
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
axis `n`, parametrised like `TensISO{order, dim, T, N}`.

Three concrete shapes are supported:

| Parametrization          | Role                               | Stored coefficients              |
|--------------------------|------------------------------------|----------------------------------|
| `TensTI{2, T, 2}`        | 2nd-order TI                       | `data = (a, b)`, `n::NTuple{3}`  |
| `TensTI{4, T, 5}`        | 4th-order TI, major-symmetric      | `data = (ℓ₁,ℓ₂,ℓ₃,ℓ₅,ℓ₆)`, `n`  |
| `TensTI{4, T, 6}`        | 4th-order TI, general              | `data = (ℓ₁,…,ℓ₆)`, `n`          |

- **Order 2** (`N=2`): `𝐀 = a·nT + b·nₙ` where `nₙ = n⊗n`, `nT = 𝟏 − nₙ`.
  `a` is the transverse coefficient, `b` the axial one.  When `a = b`, the
  tensor is isotropic and equivalent to `TensISO{2,3,T}(a)`.

- **Order 4**: stored in the Walpole basis `{W₁,…,W₆}`, with
  `W₁ = nₙ⊗nₙ`, `W₂ = (nT⊗nT)/2`, `W₃ = (nₙ⊗nT)/√2`,
  `W₄ = (nT⊗nₙ)/√2`, `W₅ = nT⊠ˢnT − (nT⊗nT)/2`, `W₆ = nT⊠ˢnₙ + nₙ⊠ˢnT`.
  Major-symmetric tensors have `ℓ₃ = ℓ₄` and are stored under `N=5`.
  Synthetic notation: `L ≡ ([[ℓ₁,ℓ₃],[ℓ₄,ℓ₂]], ℓ₅, ℓ₆)`.

See also [`tens_TI`](@ref), [`tens_TI_eng`](@ref), [`tens_TI_Hoenig`](@ref).
"""
struct TensTI{order, T, N} <: AbstractTens{order, 3, T}
    data::NTuple{N, T}
    n::NTuple{3, T}       # symmetry axis (assumed unit vector)
    TensTI{order, T, N}(data::NTuple{N, T}, n::NTuple{3, T}) where {order, T, N} =
        new{order, T, N}(data, n)
end

# ── Traits ────────────────────────────────────────────────────────────────────

@pure Base.eltype(::Type{TensTI{order, T, N}}) where {order, T, N} = T
@pure Base.length(::TensTI{order}) where {order} = 3^order
@pure Base.size(::TensTI{order}) where {order} = ntuple(_ -> 3, Val(order))

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
    TensTI{order, eltype(new_data), length(new_data)}(new_data, axis(t))

# ── 4th-order TI accessors ───────────────────────────────────────────────────

"""
    get_ℓ(t::TensTI{4,T,N}) → NTuple{6,T}

Always returns a 6-tuple `(ℓ₁, ℓ₂, ℓ₃, ℓ₄, ℓ₅, ℓ₆)` of Walpole coefficients.
For `N=5` (major-symmetric), `ℓ₃ = ℓ₄` is stored once and duplicated on read.
"""
get_ℓ(t::TensTI{4, T, 5}) where {T} =
    (t.data[1], t.data[2], t.data[3], t.data[3], t.data[4], t.data[5])
get_ℓ(t::TensTI{4, T, 6}) where {T} = t.data

# Helper: 2×2 Walpole matrix [[ℓ₁,ℓ₃],[ℓ₄,ℓ₂]]
function _walpole_mat(t::TensTI{4})
    ℓ₁, ℓ₂, ℓ₃, ℓ₄ = get_ℓ(t)[1:4]
    return SMatrix{2, 2}(ℓ₁, ℓ₄, ℓ₃, ℓ₂)   # column-major: [col1, col2] = [[ℓ₁,ℓ₄],[ℓ₃,ℓ₂]]
end

# ── Constructors ──────────────────────────────────────────────────────────────

# TensTI{4}(ℓ₁,ℓ₂,ℓ₃,ℓ₄,ℓ₅,ℓ₆, n) → TensTI{4, T, 6}
# General (not necessarily major-symmetric) TI 4th-order tensor with axis `n`.
function TensTI{4}(ℓ₁, ℓ₂, ℓ₃, ℓ₄, ℓ₅, ℓ₆, n)
    T = promote_type(
        typeof(ℓ₁), typeof(ℓ₂), typeof(ℓ₃), typeof(ℓ₄),
        typeof(ℓ₅), typeof(ℓ₆), eltype(n)
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
    T = promote_type(
        typeof(ℓ₁), typeof(ℓ₂), typeof(ℓ₃),
        typeof(ℓ₅), typeof(ℓ₆), eltype(n)
    )
    nv = _extract_vec(n)
    return TensTI{4, T, 5}(
        (T(ℓ₁), T(ℓ₂), T(ℓ₃), T(ℓ₅), T(ℓ₆)),
        (T(nv[1]), T(nv[2]), T(nv[3]))
    )
end

# Extract a plain 3-vector from various input types
_extract_vec(n::NTuple{3}) = n
_extract_vec(n::AbstractVector) = (n[1], n[2], n[3])
_extract_vec(n::AbstractTens) = _extract_vec(get_array(n))
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
eltype_of(::NTuple{N, T}) where {N, T} = T
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

"""
    Walpole(n; sym::Bool = false)

Legacy entry point: dispatches to [`walpole_basis`](@ref) (6-tuple, default)
or [`walpole_basis_sym`](@ref) (5-tuple, `sym=true`). Kept for backward
compatibility with older scripts; new code should prefer the dedicated
functions whose return arity is deterministic from the name.
"""
Walpole(n; sym::Bool = false) = sym ? walpole_basis_sym(n) : walpole_basis(n)

# ── get_array ─────────────────────────────────────────────────────────────────

"""
    get_array(t::TensTI{4, T}) → Array{T,4}

Compute the 3×3×3×3 component array from the Walpole coefficients and axis.
"""
function get_array(t::TensTI{4, T}) where {T}
    ℓ₁, ℓ₂, ℓ₃, ℓ₄, ℓ₅, ℓ₆ = get_ℓ(t)
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

Base.getindex(t::TensTI{4}, i::Integer, j::Integer, k::Integer, l::Integer) =
    get_array(t)[i, j, k, l]

# ── Kelvin-Mandel matrix ──────────────────────────────────────────────────────

"""
    KM(t::TensTI{4})

Kelvin-Mandel (6×6) matrix of the Walpole tensor.
"""
KM(t::TensTI{4}) = tomandel(tensor_or_array(get_array(t)))

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
    A.n == B.n || return Tensors.dcontract(_generic_tens(A), _generic_tens(B))
    ℓA₁, ℓA₂, ℓA₃, ℓA₄, ℓA₅, ℓA₆ = get_ℓ(A)
    ℓB₁, ℓB₂, ℓB₃, ℓB₄, ℓB₅, ℓB₆ = get_ℓ(B)
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
    ℓ₁, ℓ₂, ℓ₃, _, ℓ₅, ℓ₆ = get_ℓ(t)   # ℓ₄=ℓ₃ for N=5
    det = ℓ₁ * ℓ₂ - ℓ₃ * ℓ₃
    return TensTI{4, T, 5}((ℓ₂ / det, ℓ₁ / det, -ℓ₃ / det, one(T) / ℓ₅, one(T) / ℓ₆), t.n)
end

function Base.inv(t::TensTI{4, T, 6}) where {T}
    ℓ₁, ℓ₂, ℓ₃, ℓ₄, ℓ₅, ℓ₆ = get_ℓ(t)
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

# ── fromISO ───────────────────────────────────────────────────────────────────

"""
    fromISO(A::TensISO{4,3}, n) → TensTI{4, T, 5}

Convert an isotropic 4th-order tensor `αJ + βK` into its Walpole representation.

Formulas: ℓ₁=(α+2β)/3, ℓ₂=(2α+β)/3 (note: dim=3 → these are (3k,2μ) related),
          ℓ₃=ℓ₄=√2(α−β)/3, ℓ₅=ℓ₆=β.
Here `α` = data[1] and `β` = data[2] in TensISO (coefficients of J and K).
"""
function fromISO(A::TensISO{4, 3, T}, n) where {T}
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
    return Tensors.dcontract(A, fromISO(B, A.n))
end
function Tensors.dcontract(A::TensISO{4, 3}, B::TensTI{4})
    return Tensors.dcontract(fromISO(A, B.n), B)
end

# ── TI convenience constructors ──────────────────────────────────────────────

"""
    tens_TI(C₁₁₁₁, C₁₁₂₂, C₁₁₃₃, C₃₃₃₃, C₂₃₂₃, n) → TensTI{4, T, 5}

Construct a major-symmetric TI 4th-order tensor from its 5 independent
components and symmetry axis `n`.  Works for both stiffness and compliance
tensors (the formula is the same).

Walpole coefficients:
- `ℓ₁ = C₃₃₃₃`
- `ℓ₂ = C₁₁₁₁ + C₁₁₂₂`
- `ℓ₃ = √2 C₁₁₃₃`
- `ℓ₅ = C₁₁₁₁ − C₁₁₂₂`
- `ℓ₆ = 2 C₂₃₂₃`

See also [`arg_TI`](@ref), [`tens_TI_eng`](@ref).
"""
function tens_TI(C₁₁₁₁, C₁₁₂₂, C₁₁₃₃, C₃₃₃₃, C₂₃₂₃, n)
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

"""
    arg_TI(t::TensTI{4}) → (C₁₁₁₁, C₁₁₂₂, C₁₁₃₃, C₃₃₃₃, C₂₃₂₃)

Extract the 5 independent TI components from a Walpole tensor,
directly from the stored coefficients (no array materialisation).

Inverse of [`tens_TI`](@ref):
- `C₃₃₃₃ = ℓ₁`
- `C₁₁₁₁ = (ℓ₂ + ℓ₅)/2`
- `C₁₁₂₂ = (ℓ₂ − ℓ₅)/2`
- `C₁₁₃₃ = ℓ₃/√2`
- `C₂₃₂₃ = ℓ₆/2`

See also [`arg_TI_eng`](@ref).
"""
function arg_TI(t::TensTI{4})
    ℓ₁, ℓ₂, ℓ₃, _, ℓ₅, ℓ₆ = get_ℓ(t)
    T = eltype(t)
    sq2 = sqrt(T(2))
    C₃₃₃₃ = ℓ₁
    C₁₁₁₁ = (ℓ₂ + ℓ₅) / 2
    C₁₁₂₂ = (ℓ₂ - ℓ₅) / 2
    C₁₁₃₃ = ℓ₃ / sq2
    C₂₃₂₃ = ℓ₆ / 2
    return (C₁₁₁₁, C₁₁₂₂, C₁₁₃₃, C₃₃₃₃, C₂₃₂₃)
end

"""
    tens_TI_eng(E₁, E₃, ν₁₂, ν₃₁, G₃₁, n) → TensTI{4, T, 5}

Construct the TI **compliance** tensor from 5 engineering constants
and symmetry axis `n`.

- `E₁` : transverse Young's modulus (isotropic plane)
- `E₃` : axial Young's modulus (symmetry axis)
- `ν₁₂`: in-plane Poisson's ratio
- `ν₃₁`: axial-transverse Poisson's ratio  (`ν₃₁/E₃ = ν₁₃/E₁`)
- `G₃₁`: axial shear modulus

To obtain the stiffness tensor, invert the result: `inv(tens_TI_eng(…))`.

See also [`arg_TI_eng`](@ref), [`tens_TI`](@ref).
"""
function tens_TI_eng(E₁, E₃, ν₁₂, ν₃₁, G₃₁, n)
    S₁₁₁₁ = inv(E₁)
    S₃₃₃₃ = inv(E₃)
    S₁₁₂₂ = -ν₁₂ / E₁
    S₁₁₃₃ = -ν₃₁ / E₃
    S₂₃₂₃ = inv(4 * G₃₁)
    return tens_TI(S₁₁₁₁, S₁₁₂₂, S₁₁₃₃, S₃₃₃₃, S₂₃₂₃, n)
end

"""
    arg_TI_eng(𝕊::TensTI{4}) → (E₁, E₃, ν₁₂, ν₃₁, G₃₁)

Extract engineering constants from a TI **compliance** tensor.

See also [`tens_TI_eng`](@ref), [`arg_TI`](@ref).
"""
function arg_TI_eng(𝕊::TensTI{4})
    S₁₁₁₁, S₁₁₂₂, S₁₁₃₃, S₃₃₃₃, S₂₃₂₃ = arg_TI(𝕊)
    E₁ = inv(S₁₁₁₁)
    E₃ = inv(S₃₃₃₃)
    ν₁₂ = -E₁ * S₁₁₂₂
    ν₃₁ = -E₃ * S₁₁₃₃
    G₃₁ = inv(4 * S₂₃₂₃)
    return (E₁, E₃, ν₁₂, ν₃₁, G₃₁)
end

"""
    tens_TI_Hoenig(E, ν₁, ν₂, H, Γ, n) → TensTI{4, T, 5}

Construct the TI **compliance** tensor from 5 Hoenig parameters
(Hoenig, 1978) and symmetry axis `n`.

- `E`  : transverse Young's modulus (`= 1/S₁₁₁₁`)
- `ν₁` : in-plane Poisson's ratio (`= −E S₁₁₂₂`)
- `ν₂` : axial-transverse Poisson's ratio (`= −E S₁₁₃₃`)
- `H`  : axial-to-transverse modulus ratio (`= 1/(E S₃₃₃₃)`)
- `Γ`  : shear anisotropy parameter (`= (1+ν₁)/(2 E S₂₃₂₃)`)

Compliance components:
- `S₁₁₁₁ = 1/E`
- `S₁₁₂₂ = −ν₁/E`
- `S₁₁₃₃ = −ν₂/E`
- `S₃₃₃₃ = 1/(E H)`
- `S₂₃₂₃ = (1+ν₁)/(2 E Γ)`

To obtain the stiffness tensor, invert the result: `inv(tens_TI_Hoenig(…))`.

See also [`arg_TI_Hoenig`](@ref), [`tens_TI_eng`](@ref), [`tens_TI`](@ref).
"""
function tens_TI_Hoenig(E, ν₁, ν₂, H, Γ, n)
    S₁₁₁₁ = inv(E)
    S₃₃₃₃ = inv(E * H)
    S₁₁₂₂ = -ν₁ / E
    S₁₁₃₃ = -ν₂ / E
    S₂₃₂₃ = (1 + ν₁) / (2 * E * Γ)
    return tens_TI(S₁₁₁₁, S₁₁₂₂, S₁₁₃₃, S₃₃₃₃, S₂₃₂₃, n)
end

"""
    arg_TI_Hoenig(𝕊::TensTI{4}) → (E, ν₁, ν₂, H, Γ)

Extract the 5 Hoenig parameters from a TI **compliance** tensor.

See also [`tens_TI_Hoenig`](@ref), [`arg_TI_eng`](@ref).
"""
function arg_TI_Hoenig(𝕊::TensTI{4})
    S₁₁₁₁, S₁₁₂₂, S₁₁₃₃, S₃₃₃₃, S₂₃₂₃ = arg_TI(𝕊)
    E = inv(S₁₁₁₁)
    ν₁ = -E * S₁₁₂₂
    ν₂ = -E * S₁₁₃₃
    H = inv(S₃₃₃₃ * E)
    Γ = (1 + ν₁) / (2 * E * S₂₃₂₃)
    return (E, ν₁, ν₂, H, Γ)
end

# ── is_ISO / is_TI ─────────────────────────────────────────────────────────────

"""
    is_TI(A)

Return `true` if `A` is a `TensTI{4}`, indicating transverse isotropy.
"""
is_TI(::TensTI{4}) = true
is_TI(::Any) = false
is_ISO(::TensTI{4}) = false
is_ORTHO(::TensTI{4}) = false

# Symbolic helpers (tsimplify, tsubs, …) defined in structured_tens_ops.jl

# ── Display ───────────────────────────────────────────────────────────────────

function Base.show(io::IO, A::TensTI{4, <:Any, 5})
    ℓ₁, ℓ₂, ℓ₃, _, ℓ₅, ℓ₆ = get_ℓ(A)
    print(
        io, "(", ℓ₁, ") W₁ˢ + (", ℓ₂, ") W₂ˢ + (", ℓ₃,
        ") W₃ˢ + (", ℓ₅, ") W₄ˢ + (", ℓ₆, ") W₅ˢ"
    )
    return print(io, "\n  axis n = ", A.n)
end
function Base.show(io::IO, A::TensTI{4, <:Any, 6})
    ℓ₁, ℓ₂, ℓ₃, ℓ₄, ℓ₅, ℓ₆ = get_ℓ(A)
    print(
        io, "(", ℓ₁, ") W₁ + (", ℓ₂, ") W₂ + (", ℓ₃,
        ") W₃ + (", ℓ₄, ") W₄ + (", ℓ₅, ") W₅ + (", ℓ₆, ") W₆"
    )
    return print(io, "\n  axis n = ", A.n)
end

function intrinsic(A::TensTI{4, <:Any, 5})
    ℓ₁, ℓ₂, ℓ₃, _, ℓ₅, ℓ₆ = get_ℓ(A)
    println(
        "(", ℓ₁, ") W₁ˢ + (", ℓ₂, ") W₂ˢ + (", ℓ₃,
        ") W₃ˢ + (", ℓ₅, ") W₄ˢ + (", ℓ₆, ") W₅ˢ"
    )
    return println("  axis n = ", A.n)
end
function intrinsic(A::TensTI{4, <:Any, 6})
    ℓ₁, ℓ₂, ℓ₃, ℓ₄, ℓ₅, ℓ₆ = get_ℓ(A)
    println(
        "(", ℓ₁, ") W₁ + (", ℓ₂, ") W₂ + (", ℓ₃,
        ") W₃ + (", ℓ₄, ") W₄ + (", ℓ₅, ") W₅ + (", ℓ₆, ") W₆"
    )
    return println("  axis n = ", A.n)
end

for OP in (:show, :print, :display)
    @eval function Base.$OP(A::TensTI{4})
        $OP(typeof(A))
        print("→ decomposition: ")
        intrinsic(A)
        print("→ KM: ")
        return $OP(KM(A))
    end
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
    n = t.n
    δ(i, j) = i == j ? one(T) : zero(T)
    result = Array{T, 2}(undef, 3, 3)
    for i in 1:3, j in 1:3
        result[i, j] = a * (δ(i, j) - n[i] * n[j]) + b * n[i] * n[j]
    end
    return result
end

Base.getindex(t::TensTI{2}, i::Integer, j::Integer) =
    get_array(t)[i, j]

# ── KM ───────────────────────────────────────────────────────────────────────

KM(t::TensTI{2}) = tomandel(tensor_or_array(get_array(t)))

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

# Generic (unstructured) view of a structured tensor, used by all fallbacks.
_generic_tens(t::AbstractTens) = Tens(tensor_or_array(get_array(t)))

@inline _generic_binary(op, A::AbstractTens, B::AbstractTens) =
    Tens(tensor_or_array(broadcast(op, get_array(A), get_array(B))))

@inline function Base.:+(A::TensTI{order, <:Any, N}, B::TensTI{order, <:Any, N}) where {order, N}
    axis(A) == axis(B) || return _generic_binary(+, A, B)
    return _rebuild(A, get_data(A) .+ get_data(B))
end
@inline function Base.:-(A::TensTI{order, <:Any, N}, B::TensTI{order, <:Any, N}) where {order, N}
    axis(A) == axis(B) || return _generic_binary(-, A, B)
    return _rebuild(A, get_data(A) .- get_data(B))
end

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
@inline Base.inv(t::TensTI{2, T, 2}) where {T} =
    TensTI{2, T, 2}((one(T) / t.data[1], one(T) / t.data[2]), t.n)
@inline Base.literal_pow(::typeof(^), A::TensTI{2}, ::Val{-1}) = inv(A)

# ── Trace (order 2) ─────────────────────────────────────────────────────────

"""
    tr(t::TensTI{2}) → scalar

Trace: `tr(a·nT + b·nₙ) = 2a + b`.
"""
LinearAlgebra.tr(t::TensTI{2}) = 2 * t.data[1] + t.data[2]

# ── Symmetry ─────────────────────────────────────────────────────────────────

LinearAlgebra.issymmetric(::TensTI{2}) = true
is_ISO(t::TensTI{2}) = t.data[1] == t.data[2]
is_TI(::TensTI) = true
is_ORTHO(::TensTI) = false

# Symbolic helpers (tsimplify, tsubs, …) defined in structured_tens_ops.jl

# ── Display ──────────────────────────────────────────────────────────────────

function Base.show(io::IO, A::TensTI{2})
    a, b = get_data(A)
    print(io, "(", a, ") nT + (", b, ") nₙ")
    return print(io, "\n  axis n = ", A.n)
end

function intrinsic(A::TensTI{2})
    a, b = get_data(A)
    println("(", a, ") nT + (", b, ") nₙ")
    return println("  axis n = ", A.n)
end

for OP in (:show, :print, :display)
    @eval function Base.$OP(A::TensTI{2})
        $OP(typeof(A))
        print("→ decomposition: ")
        return intrinsic(A)
    end
end

# ── change_tens / components for TensTI ──────────────────────────────────────

change_tens(t::TensTI{2, T}, ℬ::OrthonormalBasis{3, T}) where {T} =
    Tens(tensor_or_array(get_array(t)), ℬ)
components(t::TensTI{2, T}, ::OrthonormalBasis{3, T}, ::NTuple{2, Symbol}) where {T} =
    get_array(t)
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
# TensOrtho — orthotropic 4th-order tensor
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
    TensOrtho{T} <: AbstractTens{4,3,T}

Orthotropic 4th-order tensor with material frame `(e₁,e₂,e₃)` and 9 independent
elastic constants `(C₁₁,C₂₂,C₃₃,C₁₂,C₁₃,C₂₃,C₄₄,C₅₅,C₆₆)` where
`C₄₄=C₂₃₂₃`, `C₅₅=C₁₃₁₃`, `C₆₆=C₁₂₁₂`:

    ℂ = C₁₁P₁⊗P₁ + C₂₂P₂⊗P₂ + C₃₃P₃⊗P₃
      + C₁₂(P₁⊗P₂+P₂⊗P₁) + C₁₃(P₁⊗P₃+P₃⊗P₁) + C₂₃(P₂⊗P₃+P₃⊗P₂)
      + 2C₄₄(P₂⊠ˢP₃) + 2C₅₅(P₁⊠ˢP₃) + 2C₆₆(P₁⊠ˢP₂)

with `Pₘ = eₘ⊗eₘ`. The Kelvin-Mandel matrix in the material frame is block-diagonal:

    [[C₁₁,C₁₂,C₁₃, 0,   0,   0  ],
     [C₁₂,C₂₂,C₂₃, 0,   0,   0  ],
     [C₁₃,C₂₃,C₃₃, 0,   0,   0  ],
     [ 0,  0,  0,  2C₄₄, 0,   0  ],
     [ 0,  0,  0,   0,  2C₅₅, 0  ],
     [ 0,  0,  0,   0,   0,  2C₆₆]]
"""
struct TensOrtho{T} <: AbstractTens{4, 3, T}
    data::NTuple{9, T}            # (C₁₁,C₂₂,C₃₃,C₁₂,C₁₃,C₂₃,C₄₄,C₅₅,C₆₆)
    frame::OrthonormalBasis{3}    # material frame (e₁,e₂,e₃)
    # The frame eltype is intentionally decoupled from the data eltype `T`:
    # differentiating w.r.t. the elastic constants (T = ForwardDiff.Dual)
    # must not require a Dual-typed geometric frame.  `get_array` promotes
    # the (Float64) frame against the Dual data as needed.
end

# ── Traits ────────────────────────────────────────────────────────────────────

@pure Base.eltype(::Type{TensOrtho{T}}) where {T} = T
@pure Base.length(::TensOrtho) = 81
@pure Base.size(::TensOrtho) = (3, 3, 3, 3)

get_basis(::TensOrtho{T}) where {T} = CanonicalBasis{3, T}()
get_var(::TensOrtho) = (:cont, :cont, :cont, :cont)
get_var(::TensOrtho, ::Integer) = :cont
get_data(t::TensOrtho) = t.data
"""
    frame(t::TensOrtho)

Return the material frame of an orthotropic tensor.
"""
frame(t::TensOrtho) = t.frame

# ── Rebuild helper (used by symbolic ops) ─────────────────────────────────────
_rebuild(t::TensOrtho, new_data) = TensOrtho{eltype(new_data)}(new_data, frame(t))

# ── Constructors ──────────────────────────────────────────────────────────────

"""
    TensOrtho(C11,C22,C33,C12,C13,C23,C44,C55,C66, frame)

Orthotropic tensor from the 9 elastic constants in the material frame `frame`.
"""
function TensOrtho(
        C11, C22, C33, C12, C13, C23, C44, C55, C66,
        frame::OrthonormalBasis{3}
    )
    T = promote_type(
        typeof(C11), typeof(C22), typeof(C33),
        typeof(C12), typeof(C13), typeof(C23),
        typeof(C44), typeof(C55), typeof(C66), eltype(frame)
    )
    return TensOrtho{T}(
        (
            T(C11), T(C22), T(C33), T(C12), T(C13), T(C23),
            T(C44), T(C55), T(C66),
        ), frame
    )
end

"""
    TensOrtho(KMmat::AbstractMatrix, frame)

Build a `TensOrtho` from a 6×6 Kelvin-Mandel matrix expressed in the material frame.
The matrix must have the block-diagonal orthotropic structure:
upper-left 3×3 for normal stresses and lower-right 3×3 diagonal for shear.
"""
function TensOrtho(KMmat::AbstractMatrix, frame::OrthonormalBasis{3})
    T = eltype(KMmat)
    C11 = KMmat[1, 1]; C22 = KMmat[2, 2]; C33 = KMmat[3, 3]
    C12 = KMmat[1, 2]; C13 = KMmat[1, 3]; C23 = KMmat[2, 3]
    C44 = KMmat[4, 4] / 2
    C55 = KMmat[5, 5] / 2
    C66 = KMmat[6, 6] / 2
    return TensOrtho{T}(
        (
            T(C11), T(C22), T(C33), T(C12), T(C13), T(C23),
            T(C44), T(C55), T(C66),
        ), frame
    )
end

# ── get_array ─────────────────────────────────────────────────────────────────

"""
    get_array(t::TensOrtho{T}) → Array{T,4}

Compute the 3×3×3×3 component array in the canonical frame.
"""
function get_array(t::TensOrtho{T}) where {T}
    C11, C22, C33, C12, C13, C23, C44, C55, C66 = get_data(t)
    # Frame vectors as columns of vecbasis(frame, :cov) → e[m] = frame vector m
    E = vecbasis(t.frame, :cov)   # 3×3 matrix, column m = eₘ
    result = Array{T, 4}(undef, 3, 3, 3, 3)
    # Pₘ[i,j] = E[i,m]*E[j,m]
    P(m, i, j) = E[i, m] * E[j, m]
    # (A ⊠ˢ B)[i,j,k,l] = (A[i,k]*B[j,l] + A[i,l]*B[j,k] + A[j,k]*B[i,l] + A[j,l]*B[i,k])/4
    # Note: the factor 2C in the formula accounts for the 2 in "2Cₘₘ(Pₘ⊠ˢPₙ + Pₙ⊠ˢPₘ)"
    # which is the standard Voigt-to-tensor conversion for shear moduli.
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        val = (
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
        result[i, j, k, l] = val
    end
    return result
end

Base.getindex(t::TensOrtho, i::Integer, j::Integer, k::Integer, l::Integer) =
    get_array(t)[i, j, k, l]

# ── KM in the material frame ──────────────────────────────────────────────────

"""
    KM(t::TensOrtho)

Returns the 6×6 Kelvin-Mandel matrix in the **canonical** frame.
Use `KM_material(t)` for the block-diagonal form in the material frame.
"""
KM(t::TensOrtho) = tomandel(tensor_or_array(get_array(t)))

"""
    KM_material(t::TensOrtho)

Returns the 6×6 Kelvin-Mandel matrix in the material frame (block-diagonal).
"""
function KM_material(t::TensOrtho{T}) where {T}
    C11, C22, C33, C12, C13, C23, C44, C55, C66 = get_data(t)
    z = zero(T)
    return [
        C11  C12  C13   z    z    z  ;
        C12  C22  C23   z    z    z  ;
        C13  C23  C33   z    z    z  ;
        z    z    z  2C44   z    z  ;
        z    z    z    z  2C55   z  ;
        z    z    z    z    z  2C66
    ]
end

# ── Arithmetic ────────────────────────────────────────────────────────────────
# Scalar ops (-, α*A, A*α, A/α) and _check_same_reference defined in
# structured_tens_ops.jl

@inline function Base.:+(A::TensOrtho, B::TensOrtho)
    _check_same_reference(A, B)
    return _rebuild(A, get_data(A) .+ get_data(B))
end
@inline function Base.:-(A::TensOrtho, B::TensOrtho)
    _check_same_reference(A, B)
    return _rebuild(A, get_data(A) .- get_data(B))
end

# ── Inverse ───────────────────────────────────────────────────────────────────

"""
    inv(t::TensOrtho) → TensOrtho

Inverse via the KM matrix in the material frame (block-diagonal, efficiently invertible).
"""
function Base.inv(t::TensOrtho{T}) where {T}
    Km = KM_material(t)
    Km_inv = inv(Km)
    return TensOrtho(Km_inv, t.frame)
end

@inline Base.literal_pow(::typeof(^), A::TensOrtho, ::Val{-1}) = inv(A)

# ── Symmetry ──────────────────────────────────────────────────────────────────

LinearAlgebra.issymmetric(::TensOrtho) = true
Tensors.isminorsymmetric(::TensOrtho) = true
Tensors.ismajorsymmetric(::TensOrtho) = true

# ── is_ISO / is_TI / is_ORTHO ───────────────────────────────────────────────────

is_ISO(::TensOrtho) = false
is_TI(::TensOrtho) = false
is_ORTHO(::TensOrtho) = true
is_ORTHO(::Any) = false   # universal fallback

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

julia> symmetry(tens_TI(10., 3., 2.5, 12., 2., [0., 0., 1.]))
:TI
```

See also [`reference`](@ref), [`is_ISO`](@ref), [`is_TI`](@ref),
[`is_ORTHO`](@ref), [`best_sym_tens`](@ref).
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
julia> reference(tens_TI(10., 3., 2.5, 12., 2., [0., 0., 1.]))
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
    C11, C22, C33, C12, C13, C23, C44, C55, C66 = get_data(A)
    print(io, "(", C11, ") P₁⊗P₁ + (", C22, ") P₂⊗P₂ + (", C33, ") P₃⊗P₃")
    print(io, "\n  + (", C12, ")(P₁⊗P₂+P₂⊗P₁) + (", C13, ")(P₁⊗P₃+P₃⊗P₁) + (", C23, ")(P₂⊗P₃+P₃⊗P₂)")
    print(io, "\n  + 2(", C44, ")(P₂⊠ˢP₃) + 2(", C55, ")(P₁⊠ˢP₃) + 2(", C66, ")(P₁⊠ˢP₂)")
    return print(io, "\n  frame: ", vecbasis(A.frame, :cov))
end

function intrinsic(A::TensOrtho)
    C11, C22, C33, C12, C13, C23, C44, C55, C66 = get_data(A)
    println("(", C11, ") P₁⊗P₁ + (", C22, ") P₂⊗P₂ + (", C33, ") P₃⊗P₃")
    println("  + (", C12, ")(P₁⊗P₂+P₂⊗P₁) + (", C13, ")(P₁⊗P₃+P₃⊗P₁) + (", C23, ")(P₂⊗P₃+P₃⊗P₂)")
    println("  + 2(", C44, ")(P₂⊠ˢP₃) + 2(", C55, ")(P₁⊠ˢP₃) + 2(", C66, ")(P₁⊠ˢP₂)")
    return println("  frame: ", vecbasis(A.frame, :cov))
end

for OP in (:show, :print, :display)
    @eval function Base.$OP(A::TensOrtho)
        $OP(typeof(A))
        print("→ decomposition: ")
        intrinsic(A)
        print("→ KM (material frame): ")
        $OP(KM_material(A))
        print("→ KM (canonical frame): ")
        return $OP(KM(A))
    end
end

##############################################################################
# Shared change_tens / components for TensTI{4} and TensOrtho
# (both are 3D order-4 tensors stored in the canonical frame)
##############################################################################

# TensTI{4}: T used to link tensor eltype with basis eltype
change_tens(t::TensTI{4, T}, ℬ::OrthonormalBasis{3, T}) where {T} =
    Tens(tensor_or_array(get_array(t)), ℬ)
components(t::TensTI{4, T}, ::OrthonormalBasis{3, T}, ::NTuple{4, Symbol}) where {T} =
    get_array(t)
components(t::TensTI{4}) = get_array(t)
components(t::TensTI{4}, ::NTuple{4, Symbol}) = get_array(t)

# TensOrtho
change_tens(t::TensOrtho{T}, ℬ::OrthonormalBasis{3, T}) where {T} =
    Tens(tensor_or_array(get_array(t)), ℬ)
components(t::TensOrtho{T}, ::OrthonormalBasis{3, T}, ::NTuple{4, Symbol}) where {T} =
    get_array(t)
components(t::TensOrtho) = get_array(t)
components(t::TensOrtho, ::NTuple{4, Symbol}) = get_array(t)

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
# Both W₇ and W₈ annihilate every symmetric 2nd-order tensor under double
# contraction (their minor-symmetrized structure cancels), so all existing
# 4th⊡2nd rules based on ℓ₁..ℓ₄ remain valid for N=8.
#
# Similarly, the space of 2nd-order tensors invariant under rotations about n
# is THREE-dimensional: a·nT + b·nₙ + c·w (the antisymmetric in-plane part c·w
# is what a plain symmetric TI parametrization cannot represent).
# ─────────────────────────────────────────────────────────────────────────────

# ── Constructors ─────────────────────────────────────────────────────────────

# TensTI{4}(ℓ₁,…,ℓ₈, n) → TensTI{4, T, 8}
# Full axially-invariant 4th-order tensor (see block comment above).
function TensTI{4}(ℓ₁, ℓ₂, ℓ₃, ℓ₄, ℓ₅, ℓ₆, ℓ₇, ℓ₈, n)
    T = promote_type(
        typeof(ℓ₁), typeof(ℓ₂), typeof(ℓ₃), typeof(ℓ₄),
        typeof(ℓ₅), typeof(ℓ₆), typeof(ℓ₇), typeof(ℓ₈), eltype_of(n)
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

# For N=8, `get_ℓ` returns the six Walpole coefficients (dropping ℓ₇, ℓ₈).
# This keeps every ℓ₁..ℓ₄-based 4th⊡2nd contraction rule valid (W₅..W₈
# annihilate symmetric 2nd-order tensors); use `get_ℓ8` when the
# antisymmetric couplings matter.
get_ℓ(t::TensTI{4, T, 8}) where {T} =
    (t.data[1], t.data[2], t.data[3], t.data[4], t.data[5], t.data[6])

"""
    get_ℓ8(t::TensTI{4,T,N}) → NTuple{8,T}

Always returns the 8-tuple `(ℓ₁, …, ℓ₆, ℓ₇, ℓ₈)` of coefficients in the full
axially-invariant basis `{W₁,…,W₈}`.  For `N=5`/`N=6` the antisymmetric
couplings `ℓ₇ = ℓ₈ = 0`.
"""
get_ℓ8(t::TensTI{4, T, 5}) where {T} =
    (t.data[1], t.data[2], t.data[3], t.data[3], t.data[4], t.data[5], zero(T), zero(T))
get_ℓ8(t::TensTI{4, T, 6}) where {T} =
    (t.data[1], t.data[2], t.data[3], t.data[4], t.data[5], t.data[6], zero(T), zero(T))
get_ℓ8(t::TensTI{4, T, 8}) where {T} = t.data

"""
    _lift_walpole_N8(A::TensTI{4, T, N}) → TensTI{4, T, 8}

Lift a Walpole tensor (`N=5` or `N=6`) to the full axially-invariant `N=8`
form with vanishing antisymmetric couplings.
"""
_lift_walpole_N8(A::TensTI{4, T, 5}) where {T} = TensTI{4, T, 8}(get_ℓ8(A), axis(A))
_lift_walpole_N8(A::TensTI{4, T, 6}) where {T} = TensTI{4, T, 8}(get_ℓ8(A), axis(A))
_lift_walpole_N8(A::TensTI{4, T, 8}) where {T} = A

# Lift a 2nd-order symmetric TI (N=2) to the full axially-invariant N=3 form.
_lift_ti2_N3(A::TensTI{2, T, 2}) where {T} =
    TensTI{2, T, 3}((A.data[1], A.data[2], zero(T)), axis(A))
_lift_ti2_N3(A::TensTI{2, T, 3}) where {T} = A

# ── get_array ────────────────────────────────────────────────────────────────

function get_array(t::TensTI{4, T, 8}) where {T}
    ℓ₁, ℓ₂, ℓ₃, ℓ₄, ℓ₅, ℓ₆, ℓ₇, ℓ₈ = t.data
    n = t.n
    sq2 = sqrt(T(2))
    δ(i, j) = i == j ? one(T) : zero(T)
    nn(i, j) = n[i] * n[j]
    nT(i, j) = δ(i, j) - nn(i, j)
    # In-plane rotation generator w·p = n × p  (odd in n)
    ε(i, j, k) =
        (i, j, k) in ((1, 2, 3), (2, 3, 1), (3, 1, 2)) ? one(T) :
        (i, j, k) in ((3, 2, 1), (1, 3, 2), (2, 1, 3)) ? -one(T) : zero(T)
    w(i, j) = ε(i, 1, j) * n[1] + ε(i, 2, j) * n[2] + ε(i, 3, j) * n[3]
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

function get_array(t::TensTI{2, T, 3}) where {T}
    a, b, c = t.data
    n = t.n
    δ(i, j) = i == j ? one(T) : zero(T)
    ε(i, j, k) =
        (i, j, k) in ((1, 2, 3), (2, 3, 1), (3, 1, 2)) ? one(T) :
        (i, j, k) in ((3, 2, 1), (1, 3, 2), (2, 1, 3)) ? -one(T) : zero(T)
    w(i, j) = ε(i, 1, j) * n[1] + ε(i, 2, j) * n[2] + ε(i, 3, j) * n[3]
    result = Array{T, 2}(undef, 3, 3)
    for i in 1:3, j in 1:3
        nnij = n[i] * n[j]
        result[i, j] = a * (δ(i, j) - nnij) + b * nnij + c * w(i, j)
    end
    return result
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
    A.n == B.n || return Tensors.dcontract(_generic_tens(A), _generic_tens(B))
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
    return TensTI{2, T, 3}((a / d, one(T) / b, -c / d), t.n)
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
        return $OP(_lift_walpole_N8(fromISO(A, axis(B))), B)
    end
    @eval function Base.$OP(A::TensTI{4, <:Any, 8}, B::TensISO{4, 3})
        return $OP(A, _lift_walpole_N8(fromISO(B, axis(A))))
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
    axis(A) == axis(B) || return LinearAlgebra.dot(_generic_tens(A), _generic_tens(B))
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

# ── is_ISO for the new shapes ────────────────────────────────────────────────

is_ISO(t::TensTI{2, <:Any, 3}) = t.data[1] == t.data[2] && iszero(t.data[3])

# ── Display ──────────────────────────────────────────────────────────────────

function Base.show(io::IO, A::TensTI{4, <:Any, 8})
    ℓ₁, ℓ₂, ℓ₃, ℓ₄, ℓ₅, ℓ₆, ℓ₇, ℓ₈ = A.data
    print(
        io, "(", ℓ₁, ") W₁ + (", ℓ₂, ") W₂ + (", ℓ₃, ") W₃ + (", ℓ₄,
        ") W₄ + (", ℓ₅, ") W₅ + (", ℓ₆, ") W₆ + (", ℓ₇, ") W₇ + (", ℓ₈, ") W₈"
    )
    return print(io, "\n  axis n = ", A.n)
end

function intrinsic(A::TensTI{4, <:Any, 8})
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

function intrinsic(A::TensTI{2, <:Any, 3})
    a, b, c = get_data(A)
    println("(", a, ") nT + (", b, ") nₙ + (", c, ") w")
    return println("  axis n = ", A.n)
end

##############################################################################
# Exports
##############################################################################

export TensTI, TensOrtho
export tens_W1, tens_W2, tens_W3, tens_W4, tens_W5, tens_W6, tens_W7, tens_W8
export Walpole, walpole_basis, walpole_basis_sym
export get_ℓ, get_ℓ8, axis, frame, reference, symmetry
export fromISO, is_TI, is_ORTHO
export tens_TI, arg_TI, tens_TI_eng, arg_TI_eng, tens_TI_Hoenig, arg_TI_Hoenig
export KM_material
