# ── Contraction kernels ──────────────────────────────────────────────────────
#
# Every contraction TensND performs has one of three shapes, at any order and
# for any element type (floating-point, dual, rational, symbolic):
#
# - one index of an array changed by a matrix — a change of basis or of
#   variance, a Christoffel correction: `_contract_index`;
# - the last `k` indices of one array against the first `k` of another — `⋅`,
#   `⊡`, `⊙`, `dotdot`: `_contract_trailing`;
# - two indices of one array traced against each other: `contract`.
#
# In column-major order the last `k` indices of an array vary slowest and the
# first `k` fastest, so the second shape is one matrix product of two reshapes
# (BLAS for floating-point arrays, the generic product otherwise), and the first
# is a loop over the two blocks of indices the changed one separates. `einsum`
# did the same arithmetic through an `EinCode` built anew, and dispatched at run
# time, on every call.

"""
    _contract_index(m::AbstractArray, c::AbstractMatrix, k::Integer) -> Array

The index `k` of `m` changed by the matrix `c`:
`out[…, a, …] = Σᵢ m[…, i, …] c[i, a]`, with `a` in position `k`.

The element type is that of the products, `promote_type` of the two. The sum
is accumulated with `muladd`: a fused multiply-add on floating-point numbers,
as in a BLAS kernel, and `x*y + s` on other types.
"""
function _contract_index(
        m::AbstractArray{Tm, N}, c::AbstractMatrix{Tc}, k::Integer,
    ) where {Tm, Tc, N}
    T = promote_type(Tm, Tc)
    sz = size(m)
    n = sz[k]
    size(c, 1) == n ||
        throw(DimensionMismatch("index $k of size $n changed by a matrix of size $(size(c))"))
    na = size(c, 2)
    L = prod(ntuple(j -> j < k ? sz[j] : 1, Val(N)))
    R = prod(ntuple(j -> j > k ? sz[j] : 1, Val(N)))
    M = reshape(m, L, n, R)
    out = Array{T, 3}(undef, L, na, R)
    @inbounds for r in 1:R, a in 1:na, l in 1:L
        s = zero(T)
        for i in 1:n
            s = muladd(M[l, i, r], c[i, a], s)
        end
        out[l, a, r] = s
    end
    return reshape(out, ntuple(j -> j == k ? na : sz[j], Val(N)))
end

"""
    _contract_trailing(t1::AbstractArray, t2::AbstractArray, ::Val{k}) -> Array

The last `k` indices of `t1` contracted with the first `k` of `t2`:
`out[i…, j…] = Σ_p t1[i…, p…] t2[p…, j…]`. With nothing left on either side the
result is a zero-dimensional array.
"""
function _contract_trailing(
        t1::AbstractArray{T1, N1}, t2::AbstractArray{T2, N2}, ::Val{k},
    ) where {T1, T2, N1, N2, k}
    s1, s2 = size(t1), size(t2)
    inner = ntuple(j -> s1[N1 - k + j], Val(k))
    inner == ntuple(j -> s2[j], Val(k)) || throw(
        DimensionMismatch("contracting the last $k indices of size $s1 with the first $k of size $s2")
    )
    p = prod(inner)
    C = reshape(_dense(t1), :, p) * reshape(_dense(t2), p, :)
    return reshape(C, (ntuple(j -> s1[j], Val(N1 - k))..., ntuple(j -> s2[k + j], Val(N2 - k))...))
end

# A contiguous array to reshape into a matrix: the BLAS product needs one, and
# a Tensors.jl array is not.
_dense(a::Array) = a
_dense(a::AbstractArray) = Array(a)

# The array itself when it already has the element type asked, a converted copy
# otherwise.
_narrow_eltype(::Type{T}, a::AbstractArray{T}) where {T} = a
_narrow_eltype(::Type{T}, a::AbstractArray) where {T} = T.(a)
