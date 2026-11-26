# =============================================================================
#  kelvin_mandel.jl — the Kelvin-Mandel notation, in any dimension.
#
#  A symmetric tensor of order 2 of dimension d has n = d(d+1)/2 independent
#  components, listed here in the order of Tensors.jl: the diagonal ones, then
#  the pairs (i, j), i < j, from the last column to the second one, each column
#  from bottom to top: (2,3), (1,3), (1,2) for d = 3. The off-diagonal ones
#  carry √2, so that the vector has the norm of the tensor, and a tensor of
#  order 4 with the minor symmetries is the n×n matrix of the map it is.
#
#  In dimensions 2 and 3 the components are stored in Tensors.jl types, whose
#  `tomandel` / `frommandel` are used: nothing changes there. In dimension 4 and
#  above, where Tensors.jl has no type, the same notation is written here on
#  plain arrays (`_mandel_array`, `_from_mandel_array`).
# =============================================================================

# ── Index pairs and weights of a dimension ───────────────────────────────────

"""
    mandel_pairs(Val(d)) -> NTuple{d(d+1)/2, NTuple{2, Int}}
    mandel_pairs(d::Integer) -> Vector{NTuple{2, Int}}

The pair of indices of each Kelvin-Mandel slot in dimension `d`, in the order
of Tensors.jl: the diagonal, then the off-diagonal pairs from the last column
back, `((1,1), (2,2), (3,3), (2,3), (1,3), (1,2))` for `d = 3` and
`((1,1), (2,2), (1,2))` for `d = 2`. With `Val(d)` the result is a tuple, a
constant of the compiled code; with an integer, a vector.

See also [`mandel_length`](@ref), [`mandel_weights`](@ref), [`mandel_dimension`](@ref).
"""
@generated function mandel_pairs(::Val{d}) where {d}
    pairs = Tuple{Int, Int}[(i, i) for i in 1:d]
    append!(pairs, (i, j) for j in d:-1:2 for i in (j - 1):-1:1)
    return Tuple(pairs)
end

# The same pairs for a dimension known at run time (dimension 4 and above,
# where the components are a plain array): a vector, so that what is built
# from it is inferred.
function mandel_pairs(d::Integer)
    pairs = [(i, i) for i in 1:d]
    return append!(pairs, (i, j) for j in d:-1:2 for i in (j - 1):-1:1)
end

"""
    mandel_length(d::Integer) -> Int

The number ``d(d+1)/2`` of Kelvin-Mandel components of a symmetric tensor of
order 2 in dimension `d`, the size of the Kelvin-Mandel matrix of a tensor of
order 4: 3 in the plane, 6 in space. [`mandel_dimension`](@ref) is its inverse.
"""
mandel_length(d::Integer) = d * (d + 1) ÷ 2

"""
    mandel_dimension(n::Integer) -> Int

The dimension `d` whose Kelvin-Mandel size is `n = d(d+1)/2`: the inverse of
[`mandel_length`](@ref). A size that is that of no dimension is refused.
"""
function mandel_dimension(n::Integer)
    d = (isqrt(8 * n + 1) - 1) ÷ 2
    mandel_length(d) == n || throw(DimensionMismatch("$n is the Kelvin-Mandel size d(d+1)/2 of no dimension d"))
    return d
end

"""
    mandel_weights(Val(d), T = Float64) -> NTuple{d(d+1)/2, T}

The weight of each Kelvin-Mandel slot in dimension `d`, in the order of
[`mandel_pairs`](@ref): `1` on a diagonal pair, `√2` on an off-diagonal one, so
that the Euclidean norm of the components is that of the tensor.
"""
mandel_weights(::Val{d}, ::Type{T} = Float64) where {d, T} = map(p -> _mandel_weight(T, p), mandel_pairs(Val(d)))


# ── The minor-symmetric part, read on any component array ────────────────────
#
# The Kelvin-Mandel components of the minor-symmetric part of a component array,
# read as it is, with no assumption of major symmetry (a concentration tensor
# has none), for `mandel(A; symmetrize = true)`. The four entries of a pair are
# summed two by two, so that on a minor-symmetric array the average is exact:
# `((a + a) + (a + a)) / 4 == a` to the last bit.
_minor_part(A, i, j, k, l) = ((A[i, j, k, l] + A[j, i, k, l]) + (A[i, j, l, k] + A[j, i, l, k])) / 4
_minor_part(A, i, j) = (A[i, j] + A[j, i]) / 2

# The dimension is read from the size. Dimensions 3 and 2 name their `Val`,
# whose pairs are a tuple; another dimension takes them as a vector. Every
# branch returns the same array type, so that the result is inferred.
function _mandel_minor(A::AbstractArray)
    d = size(A, 1)
    d == 3 && return _mandel_minor(A, mandel_pairs(Val(3)))
    d == 2 && return _mandel_minor(A, mandel_pairs(Val(2)))
    return _mandel_minor(A, mandel_pairs(d))
end

# The weight of a slot: one on a diagonal pair, √2 on another.
_mandel_weight(::Type{T}, (i, j)) where {T} = i == j ? one(T) : sqrt(T(2))

_mandel_minor(A::AbstractArray{T, 4}, P) where {T} =
    [_mandel_weight(T, P[I]) * _mandel_weight(T, P[J]) * _minor_part(A, P[I]..., P[J]...) for I in eachindex(P), J in eachindex(P)]
_mandel_minor(A::AbstractArray{T, 2}, P) where {T} = [_mandel_weight(T, p) * _minor_part(A, p...) for p in P]

# ── Plain arrays: dimension 4 and above ──────────────────────────────────────
#
# Tensors.jl stores a tensor with its minor symmetries (`SymmetricTensor`) or
# without (`Tensor`), and writes the first as an n×n matrix with the √2, the
# second as a d²×d² matrix of its bare components, the pairs (i, j), i > j,
# following the others. A plain array has no such type: it is written in the
# first form when its components have the minor symmetries, to the tolerance
# `_store_symmetric` applies when a tensor of dimension 2 or 3 is built, and in
# the second form otherwise.

# Whether a component array has the minor symmetries (order 4) or is symmetric
# (order 2): to 16 ulps of its largest component for floating-point element
# types, duals included, and structurally for the exact ones.
_minor_symmetric(A::AbstractArray{T}) where {T} = _minor_symmetric(Val(T <: ApproxType), A)
_minor_symmetric(::Val{true}, A::AbstractArray{T}) where {T} = _minor_gap(A) <= 16 * _approx_eps(T) * maximum(abs, A)
_minor_symmetric(::Val{false}, A::AbstractArray{T, 4}) where {T} = all(
    _sym_eq(A[i, j, k, l], A[j, i, k, l]) && _sym_eq(A[i, j, k, l], A[i, j, l, k])
        for i in axes(A, 1), j in axes(A, 2), k in axes(A, 3), l in axes(A, 4)
)
_minor_symmetric(::Val{false}, A::AbstractArray{T, 2}) where {T} =
    all(_sym_eq(A[i, j], A[j, i]) for i in axes(A, 1), j in axes(A, 2))

# The largest minor-antisymmetric component.
_minor_gap(A::AbstractArray{T, 4}) where {T} = maximum(
    max(abs(A[i, j, k, l] - A[j, i, k, l]), abs(A[i, j, k, l] - A[i, j, l, k]))
        for i in axes(A, 1), j in axes(A, 2), k in axes(A, 3), l in axes(A, 4)
)
_minor_gap(A::AbstractArray{T, 2}) where {T} = maximum(abs(A[i, j] - A[j, i]) for i in axes(A, 1), j in axes(A, 2))

# The pairs of the form without symmetry: those of the symmetric one, then
# their transposes off the diagonal.
function _full_pairs(d::Integer)
    P = mandel_pairs(d)
    return append!(P, map(reverse, P[(d + 1):end]))
end

_mandel_full(A::AbstractArray{T, 4}, P) where {T} = [A[P[I]..., P[J]...] for I in eachindex(P), J in eachindex(P)]
_mandel_full(A::AbstractArray{T, 2}, P) where {T} = [A[p...] for p in P]

"""
    _mandel_array(A) -> Matrix or Vector

The Kelvin-Mandel form of a plain component array of order 2 or 4, of any
dimension: the form of a symmetric tensor when `A` has the minor symmetries,
the form without symmetry otherwise.
"""
function _mandel_array(A::AbstractArray)
    d = size(A, 1)
    return _minor_symmetric(A) ? _mandel_minor(A, mandel_pairs(d)) : _mandel_full(A, _full_pairs(d))
end

# The dimension of a Kelvin-Mandel size `n`, and whether it is that of the
# symmetric form: `n = d(d+1)/2` is tried first, `n = d²` next. The two meet
# (36 is both 8·9/2 and 6²): `dim` settles it.
function _mandel_dimension(n::Integer, dim::Union{Nothing, Integer})
    if dim === nothing
        d = (isqrt(8 * n + 1) - 1) ÷ 2
        mandel_length(d) == n && return d, true
        d = isqrt(n)
        d^2 == n && return d, false
    else
        mandel_length(dim) == n && return Int(dim), true
        dim^2 == n && return Int(dim), false
    end
    throw(DimensionMismatch("no Kelvin-Mandel tensor of $(dim === nothing ? "any dimension" : "dimension $dim") has $n components"))
end

"""
    _from_mandel_array(v; dim = nothing) -> Array

The component array whose Kelvin-Mandel vector or matrix is `v`, the reciprocal
of [`_mandel_array`](@ref), for the sizes Tensors.jl has no type for.
"""
function _from_mandel_array(v::AbstractVector{T}; dim = nothing) where {T}
    d, symmetric = _mandel_dimension(length(v), dim)
    S = symmetric ? sqrt2_eltype(T) : T
    A = Array{S, 2}(undef, d, d)
    if symmetric
        for (I, p) in enumerate(mandel_pairs(d))
            A[p[1], p[2]] = A[p[2], p[1]] = v[I] / _mandel_weight(S, p)
        end
    else
        for (I, (i, j)) in enumerate(_full_pairs(d))
            A[i, j] = v[I]
        end
    end
    return A
end

function _from_mandel_array(M::AbstractMatrix{T}; dim = nothing) where {T}
    size(M, 1) == size(M, 2) || throw(DimensionMismatch("no Kelvin-Mandel tensor type for a $(size(M, 1))×$(size(M, 2)) matrix"))
    d, symmetric = _mandel_dimension(size(M, 1), dim)
    S = symmetric ? sqrt2_eltype(T) : T
    A = Array{S, 4}(undef, d, d, d, d)
    if symmetric
        P = mandel_pairs(d)
        for (J, q) in enumerate(P), (I, p) in enumerate(P)
            (i, j), (k, l) = p, q
            x = M[I, J] / (_mandel_weight(S, p) * _mandel_weight(S, q))
            A[i, j, k, l] = A[j, i, k, l] = A[i, j, l, k] = A[j, i, l, k] = x
        end
    else
        P = _full_pairs(d)
        for (J, (k, l)) in enumerate(P), (I, (i, j)) in enumerate(P)
            A[i, j, k, l] = M[I, J]
        end
    end
    return A
end

# ── Kelvin-Mandel notation ────────────────────────────────────────────────────

"""
    mandel(t; symmetrize = false, kwargs...)
    mandel(t, b::AbstractBasis; symmetrize = false, kwargs...)
    mandel(t, b::AbstractBasis, var::NTuple{order,Symbol}; symmetrize = false, kwargs...)

Write the components of a second or fourth order tensor in Kelvin-Mandel notation,
in its own basis or in the basis `b` (with the variances `var`, in the order of
[`components`](@ref)). `t` is a TensND tensor, a Tensors.jl one, or an array of
components. A `TensISO`, `TensTI`, `TensOrtho` or `TensCubic` is written from its
coefficients, in the canonical basis.

The size follows the symmetries of the storage: in dimension `d`, with
`n = d(d+1)/2` (6 in space, 3 in the plane), a matrix `n×n` for a fourth-order
tensor stored with its minor symmetries, `d²×d²` otherwise, and a vector of `n`
or `d²` components at order 2. In dimension 4 and above, where the components
are a plain array, the first form is taken when they have the minor symmetries. With `symmetrize = true` the
result is the Kelvin-Mandel matrix or vector of the **minor-symmetric part** of
the components, whatever the storage, with no assumption of major symmetry: the
form a concentration tensor is handled in. On components that have the minor
symmetries it reads them exactly, to the last bit.

The other keyword arguments go to `Tensors.tomandel`, in dimensions 2 and 3.

# Examples
```jldoctest
julia> σ = Tens(SymmetricTensor{2,3}((i, j) -> symbols(string("σ", i, j), real = true))) ;

julia> mandel(σ)
6-element Vector{Sym{PythonCall.Py}}:
         σ₁₁
         σ₂₂
         σ₃₃
 √2⋅σ₃₂
 √2⋅σ₃₁
 √2⋅σ₂₁

julia> C = Tens(SymmetricTensor{4,3}((i, j, k, l) -> symbols(string("C", i, j, k, l), real = true))) ;

julia> mandel(C)
6×6 Matrix{Sym{PythonCall.Py}}:
         C₁₁₁₁          C₁₁₂₂          C₁₁₃₃  …  √2⋅C₁₁₃₁  √2⋅C₁₁₂₁
         C₂₂₁₁          C₂₂₂₂          C₂₂₃₃     √2⋅C₂₂₃₁  √2⋅C₂₂₂₁
         C₃₃₁₁          C₃₃₂₂          C₃₃₃₃     √2⋅C₃₃₃₁  √2⋅C₃₃₂₁
 √2⋅C₃₂₁₁  √2⋅C₃₂₂₂  √2⋅C₃₂₃₃           2⋅C₃₂₃₁        2⋅C₃₂₂₁
 √2⋅C₃₁₁₁  √2⋅C₃₁₂₂  √2⋅C₃₁₃₃           2⋅C₃₁₃₁        2⋅C₃₁₂₁
 √2⋅C₂₁₁₁  √2⋅C₂₁₂₂  √2⋅C₂₁₃₃  …        2⋅C₂₁₃₁        2⋅C₂₁₂₁
```
"""
mandel(t::Tensors.AllTensors; symmetrize::Bool = false, kwargs...) =
    symmetrize ? _mandel_minor(t) : tomandel(t; kwargs...)
mandel(t::AbstractTens; symmetrize::Bool = false, kwargs...) =
    symmetrize ? _mandel_minor(get_array(t)) : _mandel(t; kwargs...)

# Without symmetrization, a structured tensor has a method of its own, written
# from its coefficients; any other is read from its stored components.
_mandel(t::AbstractTens; kwargs...) = _tomandel(get_array(t); kwargs...)

# Tensors.jl's notation for its own types, the same one written here for the
# plain arrays of dimension 4 and above.
_tomandel(a::Tensors.AllTensors; kwargs...) = tomandel(a; kwargs...)
_tomandel(a::AbstractArray) = _mandel_array(a)

function mandel(
        t::AbstractTens{order, dim},
        b::AbstractBasis{dim},
        var::NTuple{order, Symbol};
        symmetrize::Bool = false,
        kwargs...,
    ) where {order, dim}
    c = components(t, b, var)
    return symmetrize ? _mandel_minor(c) : _tomandel(_static_tensor(c, Val(dim), Val(order)); kwargs...)
end

function mandel(t::AbstractTens{order, dim}, b::AbstractBasis{dim}; symmetrize::Bool = false, kwargs...) where {order, dim}
    c = components(t, b)
    return symmetrize ? _mandel_minor(c) : _tomandel(_static_tensor(c, Val(dim), Val(order)); kwargs...)
end

# Constant-propagated, so that `symmetrize = true` written at the call infers
# the matrix it returns rather than the union with the generic branch.
Base.@constprop :aggressive mandel(t::AbstractArray; symmetrize::Bool = false, kwargs...) =
    symmetrize ? _mandel_minor(t) : mandel(Tens(t); kwargs...)
mandel(t::AbstractArray, b::AbstractBasis; kwargs...) = mandel(Tens(t), b; kwargs...)


"""
    from_mandel(v::AbstractVecOrMat[, b::AbstractBasis]; kwargs...)
    from_mandel(TT::Type{<:Tensors.AllTensors}, v[, b]; kwargs...)

The tensor whose Kelvin-Mandel vector or matrix is `v`, in the basis `b` (the
canonical one by default). The size of `v` gives the Tensors.jl type `TT` —
inferred when `v` is a static array — unless it is given.

The size also gives the dimension `d`: `d(d+1)/2` for the form of a symmetric
tensor, `d²` for the form without symmetry, in any dimension. A size that is
both (36 is `8·9/2` and `6²`) is read as the first; the keyword `dim` gives the
dimension otherwise.
"""
from_mandel(TT::Type{<:Tensors.AllTensors}, v::AbstractVecOrMat; kwargs...) =
    Tens(frommandel(TT, v; kwargs...))
from_mandel(TT::Type{<:Tensors.AllTensors}, v::AbstractVecOrMat, b::AbstractBasis; kwargs...) =
    Tens(frommandel(TT, v; kwargs...), b)
# The shape of `v` gives the tensor type: 6×6 → SymmetricTensor{4,3}, 9×9 →
# Tensor{4,3}, 3×3 → SymmetricTensor{4,2}, 4×4 → Tensor{4,2}, and the vectors of
# the same lengths → the order-2 types.
#
# A static array carries its size in its type, and the tensor type is chosen by
# dispatch: inferred, and the fastest route. A plain array knows its size only
# at run time, so the type is chosen by tests on it: a dispatch on
# `Val(size(v))` would be a dynamic one, measured twelve times slower. Each
# branch names its type literally, so the result is a union of four concrete
# types, which Julia splits.
from_mandel(v::StaticMatrix{6, 6}, b::AbstractBasis...; kwargs...) = from_mandel(SymmetricTensor{4, 3}, v, b...; kwargs...)
from_mandel(v::StaticMatrix{9, 9}, b::AbstractBasis...; kwargs...) = from_mandel(Tensor{4, 3}, v, b...; kwargs...)
from_mandel(v::StaticMatrix{3, 3}, b::AbstractBasis...; kwargs...) = from_mandel(SymmetricTensor{4, 2}, v, b...; kwargs...)
from_mandel(v::StaticMatrix{4, 4}, b::AbstractBasis...; kwargs...) = from_mandel(Tensor{4, 2}, v, b...; kwargs...)
from_mandel(v::StaticVector{6}, b::AbstractBasis...; kwargs...) = from_mandel(SymmetricTensor{2, 3}, v, b...; kwargs...)
from_mandel(v::StaticVector{9}, b::AbstractBasis...; kwargs...) = from_mandel(Tensor{2, 3}, v, b...; kwargs...)
from_mandel(v::StaticVector{3}, b::AbstractBasis...; kwargs...) = from_mandel(SymmetricTensor{2, 2}, v, b...; kwargs...)
from_mandel(v::StaticVector{4}, b::AbstractBasis...; kwargs...) = from_mandel(Tensor{2, 2}, v, b...; kwargs...)

function from_mandel(v::AbstractMatrix, b::AbstractBasis...; kwargs...)
    m, n = size(v)
    (m, n) == (6, 6) && return from_mandel(SymmetricTensor{4, 3}, v, b...; kwargs...)
    (m, n) == (9, 9) && return from_mandel(Tensor{4, 3}, v, b...; kwargs...)
    (m, n) == (3, 3) && return from_mandel(SymmetricTensor{4, 2}, v, b...; kwargs...)
    (m, n) == (4, 4) && return from_mandel(Tensor{4, 2}, v, b...; kwargs...)
    return Tens(_from_mandel_array(v; kwargs...), b...)
end

function from_mandel(v::AbstractVector, b::AbstractBasis...; kwargs...)
    n = length(v)
    n == 6 && return from_mandel(SymmetricTensor{2, 3}, v, b...; kwargs...)
    n == 9 && return from_mandel(Tensor{2, 3}, v, b...; kwargs...)
    n == 3 && return from_mandel(SymmetricTensor{2, 2}, v, b...; kwargs...)
    n == 4 && return from_mandel(Tensor{2, 2}, v, b...; kwargs...)
    return Tens(_from_mandel_array(v; kwargs...), b...)
end


export mandel, from_mandel, mandel_pairs, mandel_length, mandel_dimension, mandel_weights
