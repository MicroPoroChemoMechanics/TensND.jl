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

@inline Base.:-(t::AbstractTens) = Tens(.-(get_array(t)), get_basis(t), get_var(t))
@inline Base.:*(α::Number, t::AbstractTens) = Tens(α * get_array(t), get_basis(t), get_var(t))
@inline Base.:*(t::AbstractTens, α::Number) = Tens(α * get_array(t), get_basis(t), get_var(t))
@inline Base.:/(t::AbstractTens, α::Number) = Tens(get_array(t) / α, get_basis(t), get_var(t))

@inline Base.inv(t::AbstractTens{2}) =
    Tens(inv(get_array(t)), get_basis(t), (invvar(get_var(t, 2)), invvar(get_var(t, 1))))
@inline Base.inv(t::AbstractTens{4}) = Tens(
    inv(get_array(t)),
    get_basis(t),
    (
        invvar(get_var(t, 3)),
        invvar(get_var(t, 4)),
        invvar(get_var(t, 1)),
        invvar(get_var(t, 2)),
    ),
)

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

Tensors.otimes(v::AbstractTens{1, dim}) where {dim} =
    Tens(otimes(get_array(v)), get_basis(v), (get_var(v)..., get_var(v)...))

Tensors.otimes(v::TensOrthonormal{1, dim}) where {dim} =
    Tens(otimes(get_array(v)), get_basis(v))

@inline function Tensors.otimes(S::SymmetricTensor{2, dim}) where {dim}
    return SymmetricTensor{4, dim}(
        @inline function (i, j, k, l)
            return @inbounds S[i, j] * S[k, l]
        end
    )
end

Tensors.otimes(t::AbstractTens{2, dim}) where {dim} =
    Tens(otimes(get_array(t)), get_basis(t), (get_var(t)..., get_var(t)...))

Tensors.otimes(t::TensOrthonormal{2, dim}) where {dim} =
    Tens(otimes(get_array(t)), get_basis(t))


scontract(t1::AbstractArray, t2::AbstractArray) = _contract_trailing(t1, t2, Val(1))

scontract(t1::AbstractArray{T1, 1}, t2::AbstractArray{T2, 1}) where {T1, T2} =
    dot(AbstractArray{T1}(t1), AbstractArray{T2}(t2))

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
    nt2 = change_tens(nt2, get_basis(nt2), var)
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
    nt2 = change_tens(nt2, get_basis(nt2), var)
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
LinearAlgebra.norm(u::AbstractTens{1, dim}) where {dim} = √(dot(u, u))

LinearAlgebra.norm(t::AbstractTens{2, dim}) where {dim} = √(dcontract(t, t))

LinearAlgebra.norm(t::AbstractTens{4, dim}) where {dim} = √(qcontract(t, t))

"""
    isapprox(t1::AbstractTens, t2::AbstractTens; kwargs...)

Approximate equality of two tensors, compared **component-wise in the
canonical basis**. Two tensors expressed in different bases, or stored in
different concrete types (a `TensISO`, a `TensTI` and a plain `Tens` may all
denote the same tensor), therefore compare equal whenever they are the same
tensor — which the generic `AbstractArray` fallback cannot guarantee, since
it subtracts the stored components.

`kwargs` are forwarded to the underlying `isapprox` (`atol`, `rtol`, `nans`).
"""
function Base.isapprox(
        t1::AbstractTens{order, dim}, t2::AbstractTens{order, dim}; kwargs...
    ) where {order, dim}
    return isapprox(
        Array(components_canon(t1)), Array(components_canon(t2)); kwargs...
    )
end

"""
    contract(t::AbstractTens{order,dim}, i::Integer, j::Integer)

Calculate the tensor obtained after contraction with respect to the indices `i` and `j`
"""
function contract(t::AbstractTens{order, dim}, i::Integer, j::Integer) where {order, dim}
    var = ntuple(k -> k == j ? invvar(get_var(t, i)) : get_var(t, k), order)
    nt = change_tens(t, get_basis(t), var)
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
    nt = change_tens(t, get_basis(t), var)
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
    nt2 = change_tens(nt2, get_basis(nt2), var)
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
    nS, nv1 = same_basis(S, v1)
    nS, nv2 = same_basis(S, v2)
    var = (invvar(get_var(nS)[begin]),)
    nv1 = change_tens(nv1, get_basis(nv1), var)
    var = (invvar(get_var(nS)[end]),)
    nv2 = change_tens(nv2, get_basis(nv2), var)
    data = dotdot(nv1.data, nS.data, nv2.data)
    var = (get_var(nS)[begin + 1], get_var(nS)[end - 1])
    return Tens(data, get_basis(nS), var)
end

function Tensors.dotdot(
        v1::TensOrthonormal{order1, dim},
        S::TensOrthonormal{orderS, dim},
        v2::TensOrthonormal{order2, dim},
    ) where {order1, orderS, order2, dim}
    nS, nv1 = same_basis(S, v1)
    nS, nv2 = same_basis(S, v2)
    data = dotdot(nv1.data, nS.data, nv2.data)
    return Tens(data, get_basis(nS))
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
    nt2 = change_tens(nt2, get_basis(nt2), var)
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
    nt2 = change_tens(nt2, get_basis(nt2), var)
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
    nt2 = change_tens(nt2, get_basis(nt2), var)
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
    nt2 = change_tens(nt2, get_basis(nt2), var)
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
        nt2 = change_tens(nt2, get_basis(nt2), var)
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
    majortranspose(get_array(t)),
    get_basis(t),
    (get_var(t)[3], get_var(t)[4], get_var(t)[1], get_var(t)[2]),
)

Tensors.majortranspose(
    t::TensOrthonormal{order, dim, T, <:FourthOrderTensor},
) where {order, dim, T} = Tens(majortranspose(get_array(t)), get_basis(t))

Tensors.minortranspose(t::TensArray{order, dim, T, <:FourthOrderTensor}) where {order, dim, T} =
    Tens(
    minortranspose(get_array(t)),
    get_basis(t),
    (get_var(t)[2], get_var(t)[1], get_var(t)[4], get_var(t)[3]),
)

Tensors.minortranspose(
    t::TensOrthonormal{order, dim, T, <:FourthOrderTensor},
) where {order, dim, T} = Tens(minortranspose(get_array(t)), get_basis(t))
