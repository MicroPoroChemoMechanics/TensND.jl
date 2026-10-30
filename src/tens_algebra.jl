# ── Tensor algebra: arithmetic, products and contractions ────────────────────

for OP in (:+, :-)
    @eval @inline function Base.$OP(
            t1::AbstractTens{order, dim},
            t2::AbstractTens{order, dim},
        ) where {order, dim}
        nt1, nt2 = same_basis_same_var(t1, t2)
        return Tens($OP(get_array(nt1), get_array(nt2)), get_basis(nt1), get_var(nt1))
    end
    @eval @inline function Base.$OP(
            t1::AllTensOrthogonal{order, dim, T},
            t2::UniformScaling{T},
        ) where {order, dim, T <: Sym}
        nt1 = TensOrthonormal(t1)
        return Tens($OP(get_array(nt1), t2), get_basis(nt1), get_var(nt1))
    end
    @eval @inline function Base.$OP(
            t1::UniformScaling{T},
            t2::AllTensOrthogonal{order, dim, T},
        ) where {order, dim, T <: Sym}
        nt2 = TensOrthonormal(t2)
        return Tens($OP(t1, get_array(nt2)), get_basis(nt2), get_var(nt2))
    end
end

# `-` rather than `.-`: on a Tensors.jl array the broadcast returns a plain
# `Array`, whose symmetries the constructor would check back over every
# component; the negation keeps the storage, as `*` and `/` below do.
@inline Base.:-(t::AbstractTens) = Tens(-get_array(t), get_basis(t), get_var(t))
@inline Base.:*(α::Number, t::AbstractTens) = Tens(α * get_array(t), get_basis(t), get_var(t))
@inline Base.:*(t::AbstractTens, α::Number) = Tens(α * get_array(t), get_basis(t), get_var(t))
@inline Base.:/(t::AbstractTens, α::Number) = Tens(get_array(t) / α, get_basis(t), get_var(t))
# The same tie with Symbolics as for the structured tensors
# (`structured_tens_ops.jl`); the components are divided elementwise, since a
# real array divided by a `Num` would otherwise go to Symbolics' own method.
@inline Base.:/(t::AbstractTens{order, dim, <:Real}, α::Symbolics.Num) where {order, dim} =
    Tens(get_array(t) ./ α, get_basis(t), get_var(t))

@inline Base.inv(t::AbstractTens{2}) =
    Tens(inv(get_array(t)), get_basis(t), (invvar(get_var(t, 2)), invvar(get_var(t, 1))))
@inline Base.inv(t::AbstractTens{4}) = Tens(
    _inv_array(get_array(t)),
    get_basis(t),
    (
        invvar(get_var(t, 3)),
        invvar(get_var(t, 4)),
        invvar(get_var(t, 1)),
        invvar(get_var(t, 2)),
    ),
)

# The components of the inverse of a 4th-order tensor: the matrix inverse of its
# Mandel (minor symmetries) or Voigt (none) form, as Tensors.jl takes it.
#
# A tensor left in full storage by an operation although its components have
# the minor symmetries (to the tolerance of `_store_symmetric`, which decides
# the storage when a tensor is built) is singular as a 9×9 matrix: its three
# minor-antisymmetric directions are null. `inv(𝕀 + ℙ : δℂ)` then failed, or
# returned the inverse of round-off, depending on the basis the operands were
# written in. Such a tensor is inverted on the symmetric tensors, as it would
# be had it been stored symmetric.
_inv_array(A) = inv(A)
_inv_array(A::SymmetricTensor{4, dim, <:ForwardDiff.Dual}) where {dim} =
    frommandel(SymmetricTensor{4, dim}, _inv_matrix(tomandel(SMatrix, A)))
_inv_array(A::Tensor{4, dim}) where {dim} =
    _store_symmetric(A) ? Tensor{4, dim}(_inv_array(Tensors.symmetric(A))) : _inv_full(A)

_inv_full(A::Tensor{4}) = inv(A)
_inv_full(A::Tensor{4, dim, <:ForwardDiff.Dual}) where {dim} =
    fromvoigt(Tensor{4, dim}, _inv_matrix(tovoigt(SMatrix, A)))

# The inverse of a matrix of duals by the implicit-function theorem,
# d(M⁻¹) = -M⁻¹ dM M⁻¹: one inverse of the values, then two products per
# partial, recursively for nested duals. Factorizing the duals themselves would
# unroll StaticArrays' LU over the dual type, a much larger function, which
# LLVM 23 (Julia 1.14) miscompiles on AVX-512 targets (a crash in its register
# coalescer).
_inv_matrix(M::StaticMatrix) = inv(M)
function _inv_matrix(M::StaticMatrix{N, N, ForwardDiff.Dual{Tg, V, P}}) where {N, Tg, V, P}
    Mi = _inv_matrix(ForwardDiff.value.(M))
    dMi = ntuple(k -> -(Mi * ForwardDiff.partials.(M, k) * Mi), Val(P))
    return map((v, d...) -> ForwardDiff.Dual{Tg}(v, d...), Mi, dMi...)
end

"""
    otimes(t1::AbstractTens{order1,dim}, t2::AbstractTens{order2,dim})

Define a tensor product between two tensors

`(aⁱeᵢ) ⊗ (bʲeⱼ) = aⁱbʲ eᵢ⊗eⱼ`
"""
function Tensors.otimes(
        t1::AbstractTens{order1, dim},
        t2::AbstractTens{order2, dim},
    ) where {order1, order2, dim}
    if isequal(t1, t2)
        return otimes(t1)
    else
        nt1, nt2 = same_basis(t1, t2)
        data = otimes(get_array(nt1), get_array(nt2))
        var = (get_var(nt1)..., get_var(nt2)...)
        return Tens(data, get_basis(nt1), var)
    end
end

function Tensors.otimes(
        t1::TensOrthonormal{order1, dim},
        t2::TensOrthonormal{order2, dim},
    ) where {order1, order2, dim}
    if isequal(t1, t2)
        return otimes(t1)
    else
        nt1, nt2 = same_basis(t1, t2)
        data = otimes(get_array(nt1), get_array(nt2))
        return Tens(data, get_basis(nt1))
    end
end

# The product of a tensor with itself. Tensors.jl stores `v ⊗ v` and `S ⊗ S`
# (`S` symmetric) with their symmetries; any other array, a non-symmetric
# order-2 tensor or a vector beyond dimension 3, goes the general way.
_self_otimes(a::Vec) = otimes(a)
_self_otimes(a::SymmetricTensor{2}) = otimes(a)
_self_otimes(a::AbstractArray) = otimes(a, a)

Tensors.otimes(v::AbstractTens{1, dim}) where {dim} =
    Tens(_self_otimes(get_array(v)), get_basis(v), (get_var(v)..., get_var(v)...))

Tensors.otimes(v::TensOrthonormal{1, dim}) where {dim} =
    Tens(_self_otimes(get_array(v)), get_basis(v))

@inline function Tensors.otimes(S::SymmetricTensor{2, dim}) where {dim}
    return SymmetricTensor{4, dim}(
        @inline function (i, j, k, l)
            return @inbounds S[i, j] * S[k, l]
        end
    )
end

Tensors.otimes(t::AbstractTens{2, dim}) where {dim} =
    Tens(_self_otimes(get_array(t)), get_basis(t), (get_var(t)..., get_var(t)...))

Tensors.otimes(t::TensOrthonormal{2, dim}) where {dim} =
    Tens(_self_otimes(get_array(t)), get_basis(t))


scontract(t1::AbstractArray, t2::AbstractArray) = _contract_trailing(t1, t2, Val(1))

scontract(t1::AbstractArray{T1, 1}, t2::AbstractArray{T2, 1}) where {T1, T2} =
    _full_contraction(t1, t2)

for TT1 in (Vec, SecondOrderTensor), TT2 in (Vec, SecondOrderTensor)
    @eval scontract(S1::$TT1, S2::$TT2) = dot(S1, S2)
end

"""
    dot(t1::AbstractTens{order1,dim}, t2::AbstractTens{order2,dim})

Define a contracted product between two tensors

`a ⋅ b = aⁱbⱼ`
"""
function LinearAlgebra.dot(
        t1::AbstractTens{order1, dim},
        t2::AbstractTens{order2, dim},
    ) where {order1, order2, dim}
    nt1, nt2 = same_basis(t1, t2)
    var = (invvar(get_var(nt1)[end]), get_var(nt2)[(begin + 1):end]...)
    nt2 = express(nt2, get_basis(nt2), var)
    data = scontract(get_array(nt1), get_array(nt2))
    var = (get_var(nt1)[begin:(end - 1)]..., get_var(nt2)[(begin + 1):end]...)
    return Tens(data, get_basis(nt1), var)
end

function LinearAlgebra.dot(
        t1::TensOrthonormal{order1, dim},
        t2::TensOrthonormal{order2, dim},
    ) where {order1, order2, dim}
    nt1, nt2 = same_basis(t1, t2)
    data = scontract(get_array(nt1), get_array(nt2))
    return Tens(data, get_basis(nt1))
end

function LinearAlgebra.dot(t1::AbstractTens{1, dim}, t2::AbstractTens{1, dim}) where {dim}
    nt1, nt2 = same_basis(t1, t2)
    var = (invvar(get_var(nt1)[end]), get_var(nt2)[(begin + 1):end]...)
    nt2 = express(nt2, get_basis(nt2), var)
    return scontract(get_array(nt1), get_array(nt2))
end

function LinearAlgebra.dot(
        t1::TensOrthonormal{1, dim},
        t2::TensOrthonormal{1, dim},
    ) where {dim}
    nt1, nt2 = same_basis(t1, t2)
    return scontract(get_array(nt1), get_array(nt2))
end

LinearAlgebra.dot(v::TensOrthonormal{1, dim}) where {dim} = dot(get_array(v))

LinearAlgebra.dot(t::TensOrthonormal{2, dim}) where {dim} = dot(get_array(t))

"""
    norm(t::AbstractTens)

Euclidean (Frobenius) norm of a tensor: the square root of its **full**
self-contraction, `√(uᵢuⁱ)` at order 1, `√(tᵢⱼtⁱʲ)` at order 2 and
`√(tᵢⱼₖₗtⁱʲᵏˡ)` at order 4 — a scalar in every case, and independent of the
basis in which the tensor is expressed.

At order 2 this is the double contraction `⊡`, not the single contraction
`⋅`: `t ⋅ t` is again a second-order tensor, so `√(dot(t, t))` would return a
matrix rather than a number (and would then break `isapprox`, which needs a
scalar distance).
"""
LinearAlgebra.norm(t::Union{AbstractTens{1}, AbstractTens{2}, AbstractTens{4}}) = √(sqnorm(t))

"""
    sqnorm(t::AbstractTens)

The square of `norm(t)`, the self-contraction `uᵢuⁱ`, `tᵢⱼtⁱʲ` or
`tᵢⱼₖₗtⁱʲᵏˡ`, computed without the square root: smooth at zero, where the
derivative of `norm` is not, and exact for symbolic components. On complex
components it is Hermitian, `t̄ ⋅ t`: the contractions do not conjugate (they
are contractions), so the self-contraction of a complex tensor is not its squared
norm. Real element types, duals and symbols included, keep the plain
self-contraction, with no `abs` in a symbolic result.
"""
sqnorm(u::AbstractTens{1}) = dot(u, u)
sqnorm(t::AbstractTens{2}) = dcontract(t, t)
sqnorm(t::AbstractTens{4}) = qcontract(t, t)
sqnorm(u::AbstractTens{1, dim, <:Complex}) where {dim} = real(dot(conj(u), u))
sqnorm(t::AbstractTens{2, dim, <:Complex}) where {dim} = real(dcontract(conj(t), t))
sqnorm(t::AbstractTens{4, dim, <:Complex}) where {dim} = real(qcontract(conj(t), t))

export sqnorm

"""
    isapprox(t1::AbstractTens, t2::AbstractTens; kwargs...)

Approximate equality of two tensors as mathematical objects: two tensors written
in different bases, or stored in different concrete types (a `TensISO`, a `TensTI`
and a plain `Tens` may denote the same tensor), compare equal whenever they are
the same tensor — which the generic `AbstractArray` fallback cannot guarantee,
since it subtracts the stored components.

The comparison is made exactly as `==` makes it, in one basis with one set of
variances (`same_basis_same_var`): two tensors already written in the same basis
with the same variances, orthonormal or not, are compared by their arrays with no
change at all; otherwise the second is written in the basis `choose_best_basis`
picks, never in the canonical basis by default (in a moving frame, a curvilinear
one above all, its components are the costly ones).

`kwargs` are forwarded to the array `isapprox` (`atol`, `rtol`, `nans`): the
tolerance applies to the components in that common basis.
"""
function Base.isapprox(
        t1::AbstractTens{order, dim}, t2::AbstractTens{order, dim}; kwargs...
    ) where {order, dim}
    a, b = same_basis_same_var(t1, t2)
    return isapprox(get_array(a), get_array(b); kwargs...)
end

"""
    contract(t::AbstractTens{order,dim}, i::Integer, j::Integer)

Calculate the tensor obtained after contraction with respect to the indices `i` and `j`
"""
function contract(t::AbstractTens{order, dim}, i::Integer, j::Integer) where {order, dim}
    var = ntuple(k -> k == j ? invvar(get_var(t, i)) : get_var(t, k), order)
    nt = express(t, get_basis(t), var)
    data = contract(get_array(nt), i, j)
    m = min(i, j)
    M = max(i, j)
    # `get_var(nt)` once, and the two contracted slots dropped by a predicate
    # rather than by splatting three ranges: the ranges are not compile-time
    # constants, so each `...` built a tuple of unknown length and the result
    # was type-unstable — `ntuple` with `Val(order - 2)` gives the compiler the
    # length it needs.
    v = get_var(nt)
    var = ntuple(k -> v[k < m ? k : (k < M - 1 ? k + 1 : k + 2)], Val(order - 2))
    return Tens(data, get_basis(nt), var)
end

function contract(t::AbstractTens{2, dim}, i::Integer, j::Integer) where {dim}
    var = ntuple(k -> k == j ? invvar(get_var(t, i)) : get_var(t, k), 2)
    nt = express(t, get_basis(t), var)
    return contract(get_array(nt), i, j)
end

contract(t::TensOrthonormal{order, dim}, i::Integer, j::Integer) where {order, dim} =
    Tens(contract(get_array(t), i, j), get_basis(t))

contract(t::TensOrthonormal{2, dim}, i::Integer, j::Integer) where {dim} =
    contract(get_array(t), i, j)

LinearAlgebra.tr(t::AbstractTens{2}) = contract(t, 1, 2)

"""
    dcontract(t1::AbstractTens{order1,dim}, t2::AbstractTens{order2,dim})

Define a double contracted product between two tensors

`𝛔 ⊡ 𝛆 = σⁱʲεᵢⱼ`
`𝛔 = ℂ ⊡ 𝛆`

# Examples
```julia
julia> 𝛆 = Tens(SymmetricTensor{2,3}((i, j) -> symbols("ε\$i\$j", real = true))) ;

julia> k, μ = symbols("k μ", real =true) ;

julia> ℂ = 3k * t𝕁() + 2μ * t𝕂() ;

julia> 𝛔 = ℂ ⊡ 𝛆
Tens{2, 3, Sym, Sym, SymmetricTensor{2, 3, Sym, 6}, CanonicalBasis{3, Sym}}
# data: 3×3 SymmetricTensor{2, 3, Sym, 6}:
 ε11*(k + 4*μ/3) + ε22*(k - 2*μ/3) + ε33*(k - 2*μ/3)                                              2⋅ε21⋅μ                                              2⋅ε31⋅μ
                                             2⋅ε21⋅μ  ε11*(k - 2*μ/3) + ε22*(k + 4*μ/3) + ε33*(k - 2*μ/3)                                              2⋅ε32⋅μ
                                             2⋅ε31⋅μ                                              2⋅ε32⋅μ  ε11*(k - 2*μ/3) + ε22*(k - 2*μ/3) + ε33*(k + 4*μ/3)
# var: (:cont, :cont)
# basis: 3×3 Tensor{2, 3, Sym, 9}:
 1  0  0
 0  1  0
 0  0  1
```
"""
function Tensors.dcontract(
        t1::AbstractTens{order1, dim},
        t2::AbstractTens{order2, dim},
    ) where {order1, order2, dim}
    nt1, nt2 = same_basis(t1, t2)
    var = (invvar(get_var(nt1)[end - 1]), invvar(get_var(nt1)[end]), get_var(nt2)[(begin + 2):end]...)
    nt2 = express(nt2, get_basis(nt2), var)
    data = Tensors.dcontract(get_array(nt1), get_array(nt2))
    var = (get_var(nt1)[begin:(end - 2)]..., get_var(nt2)[(begin + 2):end]...)
    return Tens(data, get_basis(nt1), var)
end

function Tensors.dcontract(
        t1::TensOrthonormal{order1, dim},
        t2::TensOrthonormal{order2, dim},
    ) where {order1, order2, dim}
    nt1, nt2 = same_basis(t1, t2)
    data = Tensors.dcontract(get_array(nt1), get_array(nt2))
    return Tens(data, get_basis(nt1))
end

"""
    dotdot(v1::AbstractTens{order1,dim}, S::AbstractTens{orderS,dim}, v2::AbstractTens{order2,dim})

Define a bilinear operator `𝐯₁⋅𝕊⋅𝐯₂`

# Examples
```julia
julia> n = Tens(Sym[0, 0, 1]) ;

julia> k, μ = symbols("k μ", real =true) ;

julia> ℂ = 3k * t𝕁() + 2μ * t𝕂() ;

julia> dotdot(n,ℂ,n) # Acoustic tensor
3×3 Tens{2, 3, Sym, Sym, Tensor{2, 3, Sym, 9}, CanonicalBasis{3, Sym}}:
 μ  0          0
 0  μ          0
 0  0  k + 4*μ/3
```
"""
function Tensors.dotdot(
        v1::AbstractTens{order1, dim},
        S::AbstractTens{orderS, dim},
        v2::AbstractTens{order2, dim},
    ) where {order1, orderS, order2, dim}
    nS, nv1, nv2 = _same_basis3(S, v1, v2)
    var = (invvar(get_var(nS)[begin]),)
    nv1 = express(nv1, get_basis(nv1), var)
    var = (invvar(get_var(nS)[end]),)
    nv2 = express(nv2, get_basis(nv2), var)
    data = dotdot(get_array(nv1), get_array(nS), get_array(nv2))
    var = (get_var(nS)[begin + 1], get_var(nS)[end - 1])
    return Tens(data, get_basis(nS), var)
end

function Tensors.dotdot(
        v1::TensOrthonormal{order1, dim},
        S::TensOrthonormal{orderS, dim},
        v2::TensOrthonormal{order2, dim},
    ) where {order1, orderS, order2, dim}
    nS, nv1, nv2 = _same_basis3(S, v1, v2)
    data = dotdot(get_array(nv1), get_array(nS), get_array(nv2))
    return Tens(data, get_basis(nS))
end

# The three operands in one basis. Two calls to `same_basis(S, ⋅)` could each
# move `S` to another basis, leaving `v1` in the first and `S` in the second.
# `get_array` rather than `.data`: the data of a structured tensor is its
# coefficient tuple, not its components.
function _same_basis3(S::AbstractTens, v1::AbstractTens, v2::AbstractTens)
    nS, nv1 = same_basis(S, v1)
    nS, nv2 = same_basis(nS, v2)
    return nS, express(nv1, get_basis(nS)), nv2
end

"""
    qcontract(t1::AbstractTens{order1,dim}, t2::AbstractTens{order2,dim})

Define a quadruple contracted product between two tensors

`𝔸 ⊙ 𝔹 = AᵢⱼₖₗBⁱʲᵏˡ`

# Examples
```julia
julia> 𝕀 = t𝕀(Sym) ; 𝕁 = t𝕁(Sym) ; 𝕂 = t𝕂(Sym) ;

julia> 𝕀 ⊙ 𝕀
6

julia> 𝕁 ⊙ 𝕀
1

julia> 𝕂 ⊙ 𝕀
5

julia> 𝕂 ⊙ 𝕁
0
```
"""
function qcontract(
        t1::AbstractTens{order1, dim},
        t2::AbstractTens{order2, dim},
    ) where {order1, order2, dim}
    nt1, nt2 = same_basis(t1, t2)
    var = (
        invvar(get_var(nt1)[end - 3]),
        invvar(get_var(nt1)[end - 2]),
        invvar(get_var(nt1)[end - 1]),
        invvar(get_var(nt1)[end]),
        get_var(nt2)[(begin + 4):end]...,
    )
    nt2 = express(nt2, get_basis(nt2), var)
    data = qcontract(get_array(nt1), get_array(nt2))
    var = (get_var(nt1)[begin:(end - 4)]..., get_var(nt2)[(begin + 4):end]...)
    return Tens(data, get_basis(nt1), var)
end

function qcontract(
        t1::TensOrthonormal{order1, dim},
        t2::TensOrthonormal{order2, dim},
    ) where {order1, order2, dim}
    nt1, nt2 = same_basis(t1, t2)
    data = qcontract(get_array(nt1), get_array(nt2))
    return Tens(data, get_basis(nt1))
end

function qcontract(t1::AbstractTens{4, dim}, t2::AbstractTens{4, dim}) where {dim}
    nt1, nt2 = same_basis(t1, t2)
    var = (
        invvar(get_var(nt1)[end - 3]),
        invvar(get_var(nt1)[end - 2]),
        invvar(get_var(nt1)[end - 1]),
        invvar(get_var(nt1)[end]),
        get_var(nt2)[(begin + 4):end]...,
    )
    nt2 = express(nt2, get_basis(nt2), var)
    return qcontract(get_array(nt1), get_array(nt2))
end

function qcontract(t1::TensOrthonormal{4, dim}, t2::TensOrthonormal{4, dim}) where {dim}
    nt1, nt2 = same_basis(t1, t2)
    return qcontract(get_array(nt1), get_array(nt2))
end

"""
    otimesu(t1::AbstractTens{order1,dim}, t2::AbstractTens{order2,dim})

Define a special tensor product between two tensors of at least second order

`(𝐚 ⊠ 𝐛) ⊡ 𝐩 = 𝐚⋅𝐩⋅𝐛 = aⁱᵏbʲˡpₖₗ eᵢ⊗eⱼ`
"""
function Tensors.otimesu(
        t1::AbstractTens{order1, dim},
        t2::AbstractTens{order2, dim},
    ) where {order1, order2, dim}
    nt1, nt2 = same_basis(t1, t2)
    data = otimesu(get_array(nt1), get_array(nt2))
    var = (
        get_var(nt1)[begin:(end - 1)]...,
        get_var(nt2)[begin],
        get_var(nt1)[end],
        get_var(nt2)[(begin + 1):end]...,
    )
    return Tens(data, get_basis(nt1), var)
end

function Tensors.otimesu(
        t1::TensOrthonormal{order1, dim},
        t2::TensOrthonormal{order2, dim},
    ) where {order1, order2, dim}
    nt1, nt2 = same_basis(t1, t2)
    data = otimesu(get_array(nt1), get_array(nt2))
    return Tens(data, get_basis(nt1))
end

function Tensors.otimesl(
        t1::AbstractTens{order1, dim},
        t2::AbstractTens{order2, dim},
    ) where {order1, order2, dim}
    nt1, nt2 = same_basis(t1, t2)
    data = otimesl(get_array(nt1), get_array(nt2))
    var = (
        get_var(nt1)[begin:(end - 1)]...,
        get_var(nt2)[begin + 1],
        get_var(nt1)[end],
        get_var(nt2)[begin],
        get_var(nt2)[(begin + 2):end]...,
    )
    return Tens(data, get_basis(nt1), var)
end

function Tensors.otimesl(
        t1::TensOrthonormal{order1, dim},
        t2::TensOrthonormal{order2, dim},
    ) where {order1, order2, dim}
    nt1, nt2 = same_basis(t1, t2)
    data = otimesl(get_array(nt1), get_array(nt2))
    return Tens(data, get_basis(nt1))
end

otimesul(S1::SecondOrderTensor{dim}, S2::SecondOrderTensor{dim}) where {dim} =
    symmetric(otimesu(S1, S2))

"""
    otimesul(t1::AbstractTens{order1,dim}, t2::AbstractTens{order2,dim})

Define a special tensor product between two tensors of at least second order

`(𝐚 ⊠ˢ 𝐛) ⊡ 𝐩 = (𝐚 ⊠ 𝐛) ⊡ (𝐩 + ᵗ𝐩)/2  = 1/2(aⁱᵏbʲˡ+aⁱˡbʲᵏ) pₖₗ eᵢ⊗eⱼ`
"""
function otimesul(
        t1::AbstractTens{order1, dim},
        t2::AbstractTens{order2, dim},
    ) where {order1, order2, dim}
    nt1, nt2 = same_basis(t1, t2)
    var = (get_var(nt1)[end - 1], get_var(nt1)[end], get_var(nt2)[(begin + 2):end]...)
    nt2 = express(nt2, get_basis(nt2), var)
    data = otimesul(get_array(nt1), get_array(nt2))
    var = (get_var(nt1)..., get_var(nt2)...)
    return Tens(data, get_basis(nt1), var)
end

function otimesul(
        t1::TensOrthonormal{order1, dim},
        t2::TensOrthonormal{order2, dim},
    ) where {order1, order2, dim}
    nt1, nt2 = same_basis(t1, t2)
    data = otimesul(get_array(nt1), get_array(nt2))
    return Tens(data, get_basis(nt1))
end


"""
    sotimes(t1::AbstractTens{order1,dim}, t2::AbstractTens{order2,dim})

Define a symmetric tensor product between two tensors

`(aⁱeᵢ) ⊗ˢ (bʲeⱼ) = 1/2(aⁱbʲ + aʲbⁱ) eᵢ⊗eⱼ`
"""
function sotimes(
        t1::AbstractTens{order1, dim},
        t2::AbstractTens{order2, dim},
    ) where {order1, order2, dim}
    nt1, nt2 = same_basis(t1, t2)
    var = (get_var(nt1)[end], get_var(nt2)[(begin + 1):end]...)
    nt2 = express(nt2, get_basis(nt2), var)
    data = sotimes(get_array(nt1), get_array(nt2))
    var = (get_var(nt1)..., get_var(nt2)...)
    return Tens(data, get_basis(nt1), var)
end

function sotimes(
        t1::AbstractTens{1, dim},
        t2::AbstractTens{1, dim},
    ) where {dim}
    if isequal(t1, t2)
        return otimes(t1)
    else
        nt1, nt2 = same_basis(t1, t2)
        var = (get_var(nt1)[end], get_var(nt2)[(begin + 1):end]...)
        nt2 = express(nt2, get_basis(nt2), var)
        data = sotimes(get_array(nt1), get_array(nt2))
        var = (get_var(nt1)..., get_var(nt2)...)
        return Tens(data, get_basis(nt1), var)
    end
end

function sotimes(
        t1::TensOrthonormal{order1, dim},
        t2::TensOrthonormal{order2, dim},
    ) where {order1, order2, dim}
    nt1, nt2 = same_basis(t1, t2)
    data = sotimes(get_array(nt1), get_array(nt2))
    return Tens(data, get_basis(nt1))
end

function sotimes(
        t1::TensOrthonormal{1, dim},
        t2::TensOrthonormal{1, dim},
    ) where {dim}
    if isequal(t1, t2)
        return otimes(t1)
    else
        nt1, nt2 = same_basis(t1, t2)
        data = sotimes(get_array(nt1), get_array(nt2))
        return Tens(data, get_basis(nt1))
    end
end

Base.transpose(t::TensArray{order, dim, T, <:SecondOrderTensor}) where {order, dim, T} =
    Tens(transpose(get_array(t)), get_basis(t), (get_var(t)[2], get_var(t)[1]))

Base.transpose(t::TensOrthonormal{order, dim, T, <:SecondOrderTensor}) where {order, dim, T} =
    Tens(transpose(get_array(t)), get_basis(t))

Base.transpose(t::TensArray{order, dim, T, <:FourthOrderTensor}) where {order, dim, T} = Tens(
    Tensors.transpose(get_array(t)),
    get_basis(t),
    (get_var(t)[2], get_var(t)[1], get_var(t)[4], get_var(t)[3]),
)

Base.transpose(t::TensOrthonormal{order, dim, T, <:FourthOrderTensor}) where {order, dim, T} =
    Tens(Tensors.transpose(get_array(t)), get_basis(t))

Tensors.majortranspose(t::TensArray{order, dim, T, <:FourthOrderTensor}) where {order, dim, T} =
    Tens(
    _majortranspose(get_array(t)),
    get_basis(t),
    (get_var(t)[3], get_var(t)[4], get_var(t)[1], get_var(t)[2]),
)

Tensors.majortranspose(
    t::TensOrthonormal{order, dim, T, <:FourthOrderTensor},
) where {order, dim, T} = Tens(_majortranspose(get_array(t)), get_basis(t))

Tensors.minortranspose(t::TensArray{order, dim, T, <:FourthOrderTensor}) where {order, dim, T} =
    Tens(
    minortranspose(get_array(t)),
    get_basis(t),
    (get_var(t)[2], get_var(t)[1], get_var(t)[4], get_var(t)[3]),
)

Tensors.minortranspose(
    t::TensOrthonormal{order, dim, T, <:FourthOrderTensor},
) where {order, dim, T} = Tens(minortranspose(get_array(t)), get_basis(t))
