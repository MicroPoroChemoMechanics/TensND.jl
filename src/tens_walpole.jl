##############################################################################
# TensWalpole — transversely isotropic 4th-order tensors (Walpole basis)    #
# TensOrtho  — orthotropic 4th-order tensors                                #
##############################################################################

# ─────────────────────────────────────────────────────────────────────────────
# TensWalpole
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

"""
    TensWalpole{T,N} <: AbstractTens{4,3,T}

Transversely isotropic 4th-order tensor stored in the Walpole basis {W₁,…,W₆}
with symmetry axis `n` (assumed unit vector):

    L = ℓ₁W₁ + ℓ₂W₂ + ℓ₃W₃ + ℓ₄W₄ + ℓ₅W₅ + ℓ₆W₆

where (`nₙ = n⊗n`, `nT = 𝟏 − nₙ`):

| Tensor | Expression |
|--------|-----------|
| W₁ | `nₙ⊗nₙ` |
| W₂ | `(nT⊗nT)/2` |
| W₃ | `(nₙ⊗nT)/√2` |
| W₄ | `(nT⊗nₙ)/√2` |
| W₅ | `nT⊠ˢnT − (nT⊗nT)/2` |
| W₆ | `nT⊠ˢnₙ + nₙ⊠ˢnT` |

`N=5` (major-symmetric, `ℓ₃=ℓ₄`): `data=(ℓ₁,ℓ₂,ℓ₃,ℓ₅,ℓ₆)`.
`N=6` (general): `data=(ℓ₁,ℓ₂,ℓ₃,ℓ₄,ℓ₅,ℓ₆)`.

Synthetic notation: `L ≡ ([[ℓ₁,ℓ₃],[ℓ₄,ℓ₂]], ℓ₅, ℓ₆)`.
- Double contraction: `(L⊡M)_mat = L_mat × M_mat`, `(L⊡M)₅ = ℓ₅m₅`, `(L⊡M)₆ = ℓ₆m₆`
- Inverse: `(L⁻¹)_mat = (L_mat)⁻¹`, `1/ℓ₅`, `1/ℓ₆`
"""
struct TensWalpole{T, N} <: AbstractTens{4, 3, T}
    data::NTuple{N, T}   # N=5: (ℓ₁,ℓ₂,ℓ₃,ℓ₅,ℓ₆)  N=6: (ℓ₁,ℓ₂,ℓ₃,ℓ₄,ℓ₅,ℓ₆)
    n::NTuple{3, T}      # symmetry axis (assumed to be a unit vector)
end

# ── Traits ────────────────────────────────────────────────────────────────────

@pure Base.eltype(::Type{TensWalpole{T, N}}) where {T, N} = T
@pure Base.length(::TensWalpole) = 81   # 3^4
@pure Base.size(::TensWalpole) = (3, 3, 3, 3)

getbasis(::TensWalpole{T}) where {T} = CanonicalBasis{3, T}()
getvar(::TensWalpole) = (:cont, :cont, :cont, :cont)
getvar(::TensWalpole, ::Integer) = :cont
getdata(t::TensWalpole) = t.data

# ── Rebuild helper (used by symbolic ops) ─────────────────────────────────────
_rebuild(t::TensWalpole, new_data) =
    TensWalpole{eltype(new_data), length(new_data)}(new_data, getaxis(t))

# ── Accessors ─────────────────────────────────────────────────────────────────

"""
    get_ℓ(t::TensWalpole) → NTuple{6}

Always returns a 6-tuple `(ℓ₁,ℓ₂,ℓ₃,ℓ₄,ℓ₅,ℓ₆)`.
For N=5 (symmetric), ℓ₃ = ℓ₄ is stored once so data[3] is duplicated.
"""
get_ℓ(t::TensWalpole{T, 5}) where {T} =
    (t.data[1], t.data[2], t.data[3], t.data[3], t.data[4], t.data[5])
get_ℓ(t::TensWalpole{T, 6}) where {T} = t.data

"""
    getaxis(t::TensWalpole) → NTuple{3}

Returns the symmetry axis as a 3-tuple.
"""
getaxis(t::TensWalpole) = t.n

# Helper: 2×2 Walpole matrix [[ℓ₁,ℓ₃],[ℓ₄,ℓ₂]]
function _walpole_mat(t::TensWalpole)
    ℓ₁, ℓ₂, ℓ₃, ℓ₄ = get_ℓ(t)[1:4]
    return SMatrix{2, 2}(ℓ₁, ℓ₄, ℓ₃, ℓ₂)   # column-major: [col1, col2] = [[ℓ₁,ℓ₄],[ℓ₃,ℓ₂]]
end

# ── Constructors ──────────────────────────────────────────────────────────────

"""
    TensWalpole(ℓ₁,ℓ₂,ℓ₃,ℓ₄,ℓ₅,ℓ₆, n) → TensWalpole{T,6}

General (not necessarily major-symmetric) Walpole tensor with axis `n`.
"""
function TensWalpole(ℓ₁, ℓ₂, ℓ₃, ℓ₄, ℓ₅, ℓ₆, n)
    T = promote_type(
        typeof(ℓ₁), typeof(ℓ₂), typeof(ℓ₃), typeof(ℓ₄),
        typeof(ℓ₅), typeof(ℓ₆), eltype(n)
    )
    nv = _extract_vec(n)
    return TensWalpole{T, 6}(
        (T(ℓ₁), T(ℓ₂), T(ℓ₃), T(ℓ₄), T(ℓ₅), T(ℓ₆)),
        (T(nv[1]), T(nv[2]), T(nv[3]))
    )
end

"""
    TensWalpole(ℓ₁,ℓ₂,ℓ₃,ℓ₅,ℓ₆, n) → TensWalpole{T,5}

Major-symmetric Walpole tensor (ℓ₃ = ℓ₄), 5 independent scalars, with axis `n`.
"""
function TensWalpole(ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆, n)
    T = promote_type(
        typeof(ℓ₁), typeof(ℓ₂), typeof(ℓ₃),
        typeof(ℓ₅), typeof(ℓ₆), eltype(n)
    )
    nv = _extract_vec(n)
    return TensWalpole{T, 5}(
        (T(ℓ₁), T(ℓ₂), T(ℓ₃), T(ℓ₅), T(ℓ₆)),
        (T(nv[1]), T(nv[2]), T(nv[3]))
    )
end

# Extract a plain 3-vector from various input types
_extract_vec(n::NTuple{3}) = n
_extract_vec(n::AbstractVector) = (n[1], n[2], n[3])
_extract_vec(n::AbstractTens) = _extract_vec(getarray(n))
_extract_vec(n::Vec{3}) = (n[1], n[2], n[3])
_extract_vec(n::AbstractArray) = (n[1], n[2], n[3])

# ── Basis tensors Wᵢ ─────────────────────────────────────────────────────────

"""
    tensW1(n) → TensWalpole{T,6}   (W₁ = nₙ⊗nₙ, coeffs (1,0,0,0,0,0))
"""
tensW1(n) = TensWalpole(
    one(eltype_of(n)), zero(eltype_of(n)), zero(eltype_of(n)),
    zero(eltype_of(n)), zero(eltype_of(n)), zero(eltype_of(n)), n
)

"""
    tensW2(n) → TensWalpole{T,6}   (W₂ = (nT⊗nT)/2, coeffs (0,1,0,0,0,0))
"""
tensW2(n) = TensWalpole(
    zero(eltype_of(n)), one(eltype_of(n)), zero(eltype_of(n)),
    zero(eltype_of(n)), zero(eltype_of(n)), zero(eltype_of(n)), n
)

"""
    tensW3(n) → TensWalpole{T,6}   (W₃ = (nₙ⊗nT)/√2, coeffs (0,0,1,0,0,0))
"""
tensW3(n) = TensWalpole(
    zero(eltype_of(n)), zero(eltype_of(n)), one(eltype_of(n)),
    zero(eltype_of(n)), zero(eltype_of(n)), zero(eltype_of(n)), n
)

"""
    tensW4(n) → TensWalpole{T,6}   (W₄ = (nT⊗nₙ)/√2, coeffs (0,0,0,1,0,0))
"""
tensW4(n) = TensWalpole(
    zero(eltype_of(n)), zero(eltype_of(n)), zero(eltype_of(n)),
    one(eltype_of(n)), zero(eltype_of(n)), zero(eltype_of(n)), n
)

"""
    tensW5(n) → TensWalpole{T,6}   (W₅ = nT⊠ˢnT − (nT⊗nT)/2, coeffs (0,0,0,0,1,0))
"""
tensW5(n) = TensWalpole(
    zero(eltype_of(n)), zero(eltype_of(n)), zero(eltype_of(n)),
    zero(eltype_of(n)), one(eltype_of(n)), zero(eltype_of(n)), n
)

"""
    tensW6(n) → TensWalpole{T,6}   (W₆ = nT⊠ˢnₙ + nₙ⊠ˢnT, coeffs (0,0,0,0,0,1))
"""
tensW6(n) = TensWalpole(
    zero(eltype_of(n)), zero(eltype_of(n)), zero(eltype_of(n)),
    zero(eltype_of(n)), zero(eltype_of(n)), one(eltype_of(n)), n
)

# Helper: get element type from various axis representations
eltype_of(::AbstractArray{T}) where {T} = T
eltype_of(::NTuple{N, T}) where {N, T} = T
eltype_of(::AbstractTens{1, 3, T}) where {T} = T

"""
    Walpole(n)           → (W₁,W₂,W₃,W₄,W₅,W₆)
    Walpole(n; sym=true) → (W₁ˢ,W₂ˢ,W₃ˢ,W₄ˢ,W₅ˢ) where W₃ˢ = W₃+W₄
"""
function Walpole(n; sym::Bool = false)
    if sym
        T = eltype_of(n)
        o, z = one(T), zero(T)
        W1s = TensWalpole(o, z, z, z, z, n)         # N=5: ℓ₁=1
        W2s = TensWalpole(z, o, z, z, z, n)         # N=5: ℓ₂=1
        W3s = TensWalpole(z, z, o, z, z, n)         # N=5: ℓ₃=1  (W₃+W₄)
        W4s = TensWalpole(z, z, z, o, z, n)         # N=5: ℓ₅=1
        W5s = TensWalpole(z, z, z, z, o, n)         # N=5: ℓ₆=1
        return W1s, W2s, W3s, W4s, W5s
    else
        return tensW1(n), tensW2(n), tensW3(n), tensW4(n), tensW5(n), tensW6(n)
    end
end

# ── getarray ─────────────────────────────────────────────────────────────────

"""
    getarray(t::TensWalpole{T}) → Array{T,4}

Compute the 3×3×3×3 component array from the Walpole coefficients and axis.
"""
function getarray(t::TensWalpole{T}) where {T}
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

Base.getindex(t::TensWalpole, i::Integer, j::Integer, k::Integer, l::Integer) =
    getarray(t)[i, j, k, l]

# ── Kelvin-Mandel matrix ──────────────────────────────────────────────────────

"""
    KM(t::TensWalpole)

Kelvin-Mandel (6×6) matrix of the Walpole tensor.
"""
KM(t::TensWalpole) = tomandel(tensor_or_array(getarray(t)))

# ── Arithmetic ────────────────────────────────────────────────────────────────
# Scalar ops (-, α*A, A*α, A/α) and _check_same_reference defined in
# structured_tens_ops.jl

@inline function Base.:+(A::TensWalpole{<:Any, N}, B::TensWalpole{<:Any, N}) where {N}
    _check_same_reference(A, B)
    return _rebuild(A, getdata(A) .+ getdata(B))
end
@inline function Base.:-(A::TensWalpole{<:Any, N}, B::TensWalpole{<:Any, N}) where {N}
    _check_same_reference(A, B)
    return _rebuild(A, getdata(A) .- getdata(B))
end

# ── Double contraction (Walpole product rule) ─────────────────────────────────

"""
    dcontract(A::TensWalpole, B::TensWalpole) → TensWalpole{T,6}

Product rule via 2×2 matrix product + scalar products for ℓ₅, ℓ₆.
Always returns N=6 since the product of two symmetric tensors need not be symmetric.
"""
function Tensors.dcontract(A::TensWalpole, B::TensWalpole)
    @assert A.n == B.n "dcontract(TensWalpole,TensWalpole) requires the same axis"
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
    return TensWalpole{T, 6}((T(n₁), T(n₂), T(n₃), T(n₄), T(n₅), T(n₆)), A.n)
end

# ── Inverse ───────────────────────────────────────────────────────────────────

"""
    inv(t::TensWalpole{T,5}) → TensWalpole{T,5}
    inv(t::TensWalpole{T,6}) → TensWalpole{T,6}

Inverse via the 2×2 Walpole matrix and scalar inverses for ℓ₅, ℓ₆.
"""
function Base.inv(t::TensWalpole{T, 5}) where {T}
    ℓ₁, ℓ₂, ℓ₃, _, ℓ₅, ℓ₆ = get_ℓ(t)   # ℓ₄=ℓ₃ for N=5
    det = ℓ₁ * ℓ₂ - ℓ₃ * ℓ₃
    return TensWalpole{T, 5}((ℓ₂ / det, ℓ₁ / det, -ℓ₃ / det, one(T) / ℓ₅, one(T) / ℓ₆), t.n)
end

function Base.inv(t::TensWalpole{T, 6}) where {T}
    ℓ₁, ℓ₂, ℓ₃, ℓ₄, ℓ₅, ℓ₆ = get_ℓ(t)
    det = ℓ₁ * ℓ₂ - ℓ₃ * ℓ₄
    return TensWalpole{T, 6}((ℓ₂ / det, ℓ₁ / det, -ℓ₃ / det, -ℓ₄ / det, one(T) / ℓ₅, one(T) / ℓ₆), t.n)
end

@inline Base.literal_pow(::typeof(^), A::TensWalpole, ::Val{-1}) = inv(A)

# ── Symmetry tests ────────────────────────────────────────────────────────────

LinearAlgebra.issymmetric(::TensWalpole{T, 5}) where {T} = true
LinearAlgebra.issymmetric(t::TensWalpole{T, 6}) where {T} = isequal(t.data[3], t.data[4])
Tensors.isminorsymmetric(::TensWalpole) = true
Tensors.ismajorsymmetric(::TensWalpole{T, 5}) where {T} = true
Tensors.ismajorsymmetric(t::TensWalpole{T, 6}) where {T} = isequal(t.data[3], t.data[4])

# ── fromISO ───────────────────────────────────────────────────────────────────

"""
    fromISO(A::TensISO{4,3}, n) → TensWalpole{T,5}

Convert an isotropic 4th-order tensor `αJ + βK` into its Walpole representation.

Formulas: ℓ₁=(α+2β)/3, ℓ₂=(2α+β)/3 (note: dim=3 → these are (3k,2μ) related),
          ℓ₃=ℓ₄=√2(α−β)/3, ℓ₅=ℓ₆=β.
Here `α` = data[1] and `β` = data[2] in TensISO (coefficients of J and K).
"""
function fromISO(A::TensISO{4, 3, T}, n) where {T}
    α, β = getdata(A)    # A = α*J + β*K
    sq2 = sqrt(T(2))
    ℓ₁ = (α + 2β) / 3
    ℓ₂ = (2α + β) / 3   # Note: for 3D, 1-1/dim = 2/3 and 1/dim = 1/3
    ℓ₃ = sq2 * (α - β) / 3
    ℓ₅ = β
    ℓ₆ = β
    return TensWalpole(ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆, n)
end

"""
    dcontract(A::TensWalpole, B::TensISO{4,3}) → TensWalpole{T,6}
    dcontract(A::TensISO{4,3}, B::TensWalpole) → TensWalpole{T,6}
"""
function Tensors.dcontract(A::TensWalpole, B::TensISO{4, 3})
    return Tensors.dcontract(A, fromISO(B, A.n))
end
function Tensors.dcontract(A::TensISO{4, 3}, B::TensWalpole)
    return Tensors.dcontract(fromISO(A, B.n), B)
end

# ── TI convenience constructors ──────────────────────────────────────────────

"""
    tensTI(C₁₁₁₁, C₁₁₂₂, C₁₁₃₃, C₃₃₃₃, C₂₃₂₃, n) → TensWalpole{T,5}

Construct a major-symmetric TI 4th-order tensor from its 5 independent
components and symmetry axis `n`.  Works for both stiffness and compliance
tensors (the formula is the same).

Walpole coefficients:
- `ℓ₁ = C₃₃₃₃`
- `ℓ₂ = C₁₁₁₁ + C₁₁₂₂`
- `ℓ₃ = √2 C₁₁₃₃`
- `ℓ₅ = C₁₁₁₁ − C₁₁₂₂`
- `ℓ₆ = 2 C₂₃₂₃`

See also [`argTI`](@ref), [`tensTI_eng`](@ref).
"""
function tensTI(C₁₁₁₁, C₁₁₂₂, C₁₁₃₃, C₃₃₃₃, C₂₃₂₃, n)
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
    return TensWalpole(ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆, n)
end

"""
    argTI(t::TensWalpole) → (C₁₁₁₁, C₁₁₂₂, C₁₁₃₃, C₃₃₃₃, C₂₃₂₃)

Extract the 5 independent TI components from a Walpole tensor,
directly from the stored coefficients (no array materialisation).

Inverse of [`tensTI`](@ref):
- `C₃₃₃₃ = ℓ₁`
- `C₁₁₁₁ = (ℓ₂ + ℓ₅)/2`
- `C₁₁₂₂ = (ℓ₂ − ℓ₅)/2`
- `C₁₁₃₃ = ℓ₃/√2`
- `C₂₃₂₃ = ℓ₆/2`

See also [`argTI_eng`](@ref).
"""
function argTI(t::TensWalpole)
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
    tensTI_eng(E₁, E₃, ν₁₂, ν₃₁, G₃₁, n) → TensWalpole{T,5}

Construct the TI **compliance** tensor from 5 engineering constants
and symmetry axis `n`.

- `E₁` : transverse Young's modulus (isotropic plane)
- `E₃` : axial Young's modulus (symmetry axis)
- `ν₁₂`: in-plane Poisson's ratio
- `ν₃₁`: axial-transverse Poisson's ratio  (`ν₃₁/E₃ = ν₁₃/E₁`)
- `G₃₁`: axial shear modulus

To obtain the stiffness tensor, invert the result: `inv(tensTI_eng(…))`.

See also [`argTI_eng`](@ref), [`tensTI`](@ref).
"""
function tensTI_eng(E₁, E₃, ν₁₂, ν₃₁, G₃₁, n)
    S₁₁₁₁ = inv(E₁)
    S₃₃₃₃ = inv(E₃)
    S₁₁₂₂ = -ν₁₂ / E₁
    S₁₁₃₃ = -ν₃₁ / E₃
    S₂₃₂₃ = inv(4 * G₃₁)
    return tensTI(S₁₁₁₁, S₁₁₂₂, S₁₁₃₃, S₃₃₃₃, S₂₃₂₃, n)
end

"""
    argTI_eng(𝕊::TensWalpole) → (E₁, E₃, ν₁₂, ν₃₁, G₃₁)

Extract engineering constants from a TI **compliance** tensor.

See also [`tensTI_eng`](@ref), [`argTI`](@ref).
"""
function argTI_eng(𝕊::TensWalpole)
    S₁₁₁₁, S₁₁₂₂, S₁₁₃₃, S₃₃₃₃, S₂₃₂₃ = argTI(𝕊)
    E₁ = inv(S₁₁₁₁)
    E₃ = inv(S₃₃₃₃)
    ν₁₂ = -E₁ * S₁₁₂₂
    ν₃₁ = -E₃ * S₁₁₃₃
    G₃₁ = inv(4 * S₂₃₂₃)
    return (E₁, E₃, ν₁₂, ν₃₁, G₃₁)
end

"""
    tensTI_Hoenig(E, ν₁, ν₂, H, Γ, n) → TensWalpole{T,5}

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

To obtain the stiffness tensor, invert the result: `inv(tensTI_Hoenig(…))`.

See also [`argTI_Hoenig`](@ref), [`tensTI_eng`](@ref), [`tensTI`](@ref).
"""
function tensTI_Hoenig(E, ν₁, ν₂, H, Γ, n)
    S₁₁₁₁ = inv(E)
    S₃₃₃₃ = inv(E * H)
    S₁₁₂₂ = -ν₁ / E
    S₁₁₃₃ = -ν₂ / E
    S₂₃₂₃ = (1 + ν₁) / (2 * E * Γ)
    return tensTI(S₁₁₁₁, S₁₁₂₂, S₁₁₃₃, S₃₃₃₃, S₂₃₂₃, n)
end

"""
    argTI_Hoenig(𝕊::TensWalpole) → (E, ν₁, ν₂, H, Γ)

Extract the 5 Hoenig parameters from a TI **compliance** tensor.

See also [`tensTI_Hoenig`](@ref), [`argTI_eng`](@ref).
"""
function argTI_Hoenig(𝕊::TensWalpole)
    S₁₁₁₁, S₁₁₂₂, S₁₁₃₃, S₃₃₃₃, S₂₃₂₃ = argTI(𝕊)
    E = inv(S₁₁₁₁)
    ν₁ = -E * S₁₁₂₂
    ν₂ = -E * S₁₁₃₃
    H = inv(S₃₃₃₃ * E)
    Γ = (1 + ν₁) / (2 * E * S₂₃₂₃)
    return (E, ν₁, ν₂, H, Γ)
end

# ── isISO / isTI ─────────────────────────────────────────────────────────────

"""
    isTI(A)

Return `true` if `A` is a `TensWalpole`, indicating transverse isotropy.
"""
isTI(::TensWalpole) = true
isTI(::Any) = false
isISO(::TensWalpole) = false
isOrtho(::TensWalpole) = false

# Symbolic helpers (tsimplify, tsubs, …) defined in structured_tens_ops.jl

# ── Display ───────────────────────────────────────────────────────────────────

function Base.show(io::IO, A::TensWalpole{<:Any, 5})
    ℓ₁, ℓ₂, ℓ₃, _, ℓ₅, ℓ₆ = get_ℓ(A)
    print(
        io, "(", ℓ₁, ") W₁ˢ + (", ℓ₂, ") W₂ˢ + (", ℓ₃,
        ") W₃ˢ + (", ℓ₅, ") W₄ˢ + (", ℓ₆, ") W₅ˢ"
    )
    return print(io, "\n  axis n = ", A.n)
end
function Base.show(io::IO, A::TensWalpole{<:Any, 6})
    ℓ₁, ℓ₂, ℓ₃, ℓ₄, ℓ₅, ℓ₆ = get_ℓ(A)
    print(
        io, "(", ℓ₁, ") W₁ + (", ℓ₂, ") W₂ + (", ℓ₃,
        ") W₃ + (", ℓ₄, ") W₄ + (", ℓ₅, ") W₅ + (", ℓ₆, ") W₆"
    )
    return print(io, "\n  axis n = ", A.n)
end

function intrinsic(A::TensWalpole{<:Any, 5})
    ℓ₁, ℓ₂, ℓ₃, _, ℓ₅, ℓ₆ = get_ℓ(A)
    println(
        "(", ℓ₁, ") W₁ˢ + (", ℓ₂, ") W₂ˢ + (", ℓ₃,
        ") W₃ˢ + (", ℓ₅, ") W₄ˢ + (", ℓ₆, ") W₅ˢ"
    )
    return println("  axis n = ", A.n)
end
function intrinsic(A::TensWalpole{<:Any, 6})
    ℓ₁, ℓ₂, ℓ₃, ℓ₄, ℓ₅, ℓ₆ = get_ℓ(A)
    println(
        "(", ℓ₁, ") W₁ + (", ℓ₂, ") W₂ + (", ℓ₃,
        ") W₃ + (", ℓ₄, ") W₄ + (", ℓ₅, ") W₅ + (", ℓ₆, ") W₆"
    )
    return println("  axis n = ", A.n)
end

for OP in (:show, :print, :display)
    @eval function Base.$OP(A::TensWalpole)
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
# Order 4 is handled by TensWalpole{T,N} (Walpole basis, N=5 or 6).
# Future unification TensWalpole → TensTI{4,T,N} is possible.
# ─────────────────────────────────────────────────────────────────────────────

"""
    TensTI{order,T,N} <: AbstractTens{order,3,T}

Transversely isotropic tensor of order `order` (always dim=3) with symmetry
axis `n`, parametrised like `TensISO{order,dim,T,N}`.

**Order 2** (`N=2`): `data = (a, b)` → `𝐀 = a·nT + b·nₙ`
where `nₙ = n⊗n`, `nT = 𝟏 − nₙ`.
- `a`: transverse coefficient (plane ⊥ n)
- `b`: axial coefficient (along n)
- When `a = b`: isotropic, equivalent to `TensISO{2,3,T}(a)`

**Order 4** is handled separately by [`TensWalpole{T,N}`](@ref) (Walpole basis
with 5 or 6 coefficients).

# Constructor

    TensTI{2}(a, b, n) → TensTI{2,T,2}

Construct a TI 2nd-order tensor `a·nT + b·nₙ` with symmetry axis `n`.
The axis `n` can be a `Vector`, `NTuple{3}`, or any `AbstractTens` of order 1.

# Examples
```julia
julia> n = [0., 0., 1.];

julia> A = TensTI{2}(5.0, 8.0, n);

julia> getarray(A)
3×3 Matrix{Float64}:
 5.0  0.0  0.0
 0.0  5.0  0.0
 0.0  0.0  8.0

julia> tr(A)
18.0

julia> isISO(A)
false

julia> inv(A).data
(0.2, 0.125)

julia> B = TensTI{2}(5.0, 5.0, n); isISO(B)
true
```
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

getbasis(::TensTI{order, T}) where {order, T} = CanonicalBasis{3, T}()
getvar(::TensTI{order}) where {order} = ntuple(_ -> :cont, Val(order))
getvar(::TensTI, ::Integer) = :cont
getdata(t::TensTI) = t.data
getaxis(t::TensTI) = t.n

# ── Rebuild helper (used by symbolic ops) ─────────────────────────────────────
_rebuild(t::TensTI{order}, new_data) where {order} =
    TensTI{order, eltype(new_data), length(new_data)}(new_data, getaxis(t))

# ── Convenience constructors ─────────────────────────────────────────────────

# TensTI{2}(a, b, n) → TensTI{2,T,2}
#
# Construct a TI 2nd-order tensor `a·nT + b·nₙ` with symmetry axis `n`.
#
# Examples:
#   n = [0., 0., 1.]
#   A = TensTI{2}(5.0, 8.0, n)
#   getarray(A) → [5 0 0; 0 5 0; 0 0 8]
function TensTI{2}(a, b, n)
    T = promote_type(typeof(a), typeof(b), eltype(n))
    nv = _extract_vec(n)
    return TensTI{2, T, 2}((T(a), T(b)), (T(nv[1]), T(nv[2]), T(nv[3])))
end

# ── getarray (order 2) ──────────────────────────────────────────────────────

"""
    getarray(t::TensTI{2,T,2}) → Array{T,2}

Compute the 3×3 component array: `a*(δᵢⱼ − nᵢnⱼ) + b*nᵢnⱼ`.

# Examples
```julia
julia> A = TensTI{2}(5.0, 8.0, [0., 0., 1.]);

julia> getarray(A)
3×3 Matrix{Float64}:
 5.0  0.0  0.0
 0.0  5.0  0.0
 0.0  0.0  8.0
```
"""
function getarray(t::TensTI{2, T, 2}) where {T}
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
    getarray(t)[i, j]

# ── KM ───────────────────────────────────────────────────────────────────────

KM(t::TensTI{2}) = tomandel(tensor_or_array(getarray(t)))

# ── Arithmetic ───────────────────────────────────────────────────────────────
# Scalar ops (-, α*A, A*α, A/α) and _check_same_reference defined in
# structured_tens_ops.jl

@inline function Base.:+(A::TensTI{order, <:Any, N}, B::TensTI{order, <:Any, N}) where {order, N}
    _check_same_reference(A, B)
    return _rebuild(A, getdata(A) .+ getdata(B))
end
@inline function Base.:-(A::TensTI{order, <:Any, N}, B::TensTI{order, <:Any, N}) where {order, N}
    _check_same_reference(A, B)
    return _rebuild(A, getdata(A) .- getdata(B))
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
isISO(t::TensTI{2}) = t.data[1] == t.data[2]
isTI(::TensTI) = true
isOrtho(::TensTI) = false

# Symbolic helpers (tsimplify, tsubs, …) defined in structured_tens_ops.jl

# ── Display ──────────────────────────────────────────────────────────────────

function Base.show(io::IO, A::TensTI{2})
    a, b = getdata(A)
    print(io, "(", a, ") nT + (", b, ") nₙ")
    return print(io, "\n  axis n = ", A.n)
end

function intrinsic(A::TensTI{2})
    a, b = getdata(A)
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
    Tens(tensor_or_array(getarray(t)), ℬ)
components(t::TensTI{2, T}, ::OrthonormalBasis{3, T}, ::NTuple{2, Symbol}) where {T} =
    getarray(t)
components(t::TensTI{2}) = getarray(t)
components(t::TensTI{2}, ::NTuple{2, Symbol}) = getarray(t)

# ── otimes specializations (TensTI{2} → TensWalpole) ───────────────────────

"""
    otimes(A::TensTI{2}) → TensWalpole{T,5}

Self tensor product of a TI 2nd-order tensor.  The result is always
major-symmetric (ℓ₃ = ℓ₄) and lives in the Walpole basis with N=5.

    (a·nT + b·nₙ) ⊗ (a·nT + b·nₙ)
    = b²W₁ + 2a²W₂ + √2·ab·(W₃+W₄)
"""
function Tensors.otimes(A::TensTI{2, T, 2}) where {T}
    a, b = A.data
    sq2 = sqrt(T(2))
    return TensWalpole{T, 5}((b * b, T(2) * a * a, sq2 * a * b, zero(T), zero(T)), A.n)
end

"""
    otimes(A::TensTI{2}, B::TensTI{2}) → TensWalpole{T,6}

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
    return TensWalpole{T, 6}(
        (
            T(b₁ * b₂), T(2) * a₁ * a₂, sq2 * T(b₁ * a₂), sq2 * T(a₁ * b₂),
            zero(T), zero(T),
        ), A.n
    )
end

"""
    otimes(A::TensISO{2,3}, B::TensTI{2}) → TensWalpole{T,6}

Tensor product of a 3D isotropic 2nd-order tensor with a TI 2nd-order tensor.
The isotropic tensor `λ·𝟏` is treated as `TensTI{2}(λ,λ,n)` with the axis of B.
"""
function Tensors.otimes(A::TensISO{2, 3}, B::TensTI{2, T2, 2}) where {T2}
    T = promote_type(eltype(A), T2)
    λ = A.data[1]
    a₂, b₂ = B.data
    sq2 = sqrt(T(2))
    return TensWalpole{T, 6}(
        (
            T(λ * b₂), T(2) * λ * a₂, sq2 * T(λ * a₂), sq2 * T(λ * b₂),
            zero(T), zero(T),
        ), B.n
    )
end

"""
    otimes(A::TensTI{2}, B::TensISO{2,3}) → TensWalpole{T,6}

Tensor product of a TI 2nd-order tensor with a 3D isotropic 2nd-order tensor.
The isotropic tensor `λ·𝟏` is treated as `TensTI{2}(λ,λ,n)` with the axis of A.
"""
function Tensors.otimes(A::TensTI{2, T1, 2}, B::TensISO{2, 3}) where {T1}
    T = promote_type(T1, eltype(B))
    a₁, b₁ = A.data
    λ = B.data[1]
    sq2 = sqrt(T(2))
    return TensWalpole{T, 6}(
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
    frame::OrthonormalBasis{3, T} # material frame (e₁,e₂,e₃)
end

# ── Traits ────────────────────────────────────────────────────────────────────

@pure Base.eltype(::Type{TensOrtho{T}}) where {T} = T
@pure Base.length(::TensOrtho) = 81
@pure Base.size(::TensOrtho) = (3, 3, 3, 3)

getbasis(::TensOrtho{T}) where {T} = CanonicalBasis{3, T}()
getvar(::TensOrtho) = (:cont, :cont, :cont, :cont)
getvar(::TensOrtho, ::Integer) = :cont
getdata(t::TensOrtho) = t.data
getframe(t::TensOrtho) = t.frame

# ── Rebuild helper (used by symbolic ops) ─────────────────────────────────────
_rebuild(t::TensOrtho, new_data) = TensOrtho{eltype(new_data)}(new_data, getframe(t))

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

# ── getarray ─────────────────────────────────────────────────────────────────

"""
    getarray(t::TensOrtho{T}) → Array{T,4}

Compute the 3×3×3×3 component array in the canonical frame.
"""
function getarray(t::TensOrtho{T}) where {T}
    C11, C22, C33, C12, C13, C23, C44, C55, C66 = getdata(t)
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
    getarray(t)[i, j, k, l]

# ── KM in the material frame ──────────────────────────────────────────────────

"""
    KM(t::TensOrtho)

Returns the 6×6 Kelvin-Mandel matrix in the **canonical** frame.
Use `KM_material(t)` for the block-diagonal form in the material frame.
"""
KM(t::TensOrtho) = tomandel(tensor_or_array(getarray(t)))

"""
    KM_material(t::TensOrtho)

Returns the 6×6 Kelvin-Mandel matrix in the material frame (block-diagonal).
"""
function KM_material(t::TensOrtho{T}) where {T}
    C11, C22, C33, C12, C13, C23, C44, C55, C66 = getdata(t)
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
    return _rebuild(A, getdata(A) .+ getdata(B))
end
@inline function Base.:-(A::TensOrtho, B::TensOrtho)
    _check_same_reference(A, B)
    return _rebuild(A, getdata(A) .- getdata(B))
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

# ── isISO / isTI / isOrtho ───────────────────────────────────────────────────

isISO(::TensOrtho) = false
isTI(::TensOrtho) = false
isOrtho(::TensOrtho) = true
isOrtho(::Any) = false   # universal fallback

# Symbolic helpers (tsimplify, tsubs, …) defined in structured_tens_ops.jl

# ── Display ───────────────────────────────────────────────────────────────────

function Base.show(io::IO, A::TensOrtho)
    C11, C22, C33, C12, C13, C23, C44, C55, C66 = getdata(A)
    print(io, "(", C11, ") P₁⊗P₁ + (", C22, ") P₂⊗P₂ + (", C33, ") P₃⊗P₃")
    print(io, "\n  + (", C12, ")(P₁⊗P₂+P₂⊗P₁) + (", C13, ")(P₁⊗P₃+P₃⊗P₁) + (", C23, ")(P₂⊗P₃+P₃⊗P₂)")
    print(io, "\n  + 2(", C44, ")(P₂⊠ˢP₃) + 2(", C55, ")(P₁⊠ˢP₃) + 2(", C66, ")(P₁⊠ˢP₂)")
    return print(io, "\n  frame: ", vecbasis(A.frame, :cov))
end

function intrinsic(A::TensOrtho)
    C11, C22, C33, C12, C13, C23, C44, C55, C66 = getdata(A)
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
# Shared change_tens / components for TensWalpole and TensOrtho
# (both are 3D order-4 tensors stored in the canonical frame)
##############################################################################

for TT in (:TensWalpole, :TensOrtho)
    # T used to link tensor eltype with basis eltype:
    @eval change_tens(t::$TT{T}, ℬ::OrthonormalBasis{3, T}) where {T} =
        Tens(tensor_or_array(getarray(t)), ℬ)
    @eval components(t::$TT{T}, ::OrthonormalBasis{3, T}, ::NTuple{4, Symbol}) where {T} =
        getarray(t)
    # T not needed for these:
    @eval components(t::$TT) = getarray(t)
    @eval components(t::$TT, ::NTuple{4, Symbol}) = getarray(t)
end

##############################################################################
# Exports
##############################################################################

export TensWalpole, TensTI, TensOrtho
export tensW1, tensW2, tensW3, tensW4, tensW5, tensW6, Walpole
export get_ℓ, getaxis, getframe
export fromISO, isTI, isOrtho
export tensTI, argTI, tensTI_eng, argTI_eng, tensTI_Hoenig, argTI_Hoenig
export KM_material
