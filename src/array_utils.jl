δkron(T::Type{<:Number}, i::Integer, j::Integer) = i == j ? one(T) : zero(T)

struct Id2{dim, T <: Number} <: AbstractMatrix{T} end
@pure Base.size(::Id2{dim}) where {dim} = (dim, dim)
Base.getindex(::Id2{dim, T}, i::Integer, j::Integer) where {dim, T} = δkron(T, i, j)
function Base.replace_in_print_matrix(::Id2, i::Integer, j::Integer, s::AbstractString)
    return i == j ? s : Base.replace_with_centered_mark(s)
end

isidentity(a::AbstractMatrix{T}) where {T} = a ≈ I
isdiagonal(a::AbstractMatrix{T}) where {T} = norm(a - Diagonal(a)) <= eps(T)

# isapprox(x::Num, y::Num; kwargs...) =

for OP in (:(tsimplify), :(tfactor), :(tsubs), :(ttrigsimp), :(texpand_trig))
    @eval $OP(x, args...; kwargs...) = x
end

for OP in (:(simplify), :(factor), :(subs), :(diff))
    @eval $(Symbol("t", OP))(x::T, args...; kwargs...) where {T <: Sym} = SymPy.$OP(x, args...; kwargs...)
end
for OP in (:(trigsimp), :(expand_trig))
    @eval $(Symbol("t", OP))(x::T, args...; kwargs...) where {T <: Sym} = sympy.$OP(x, args...; kwargs...)
end
for OP in (:(tsimplify), :(tfactor), :(tsubs), :(tdiff), :(ttrigsimp), :(texpand_trig))
    @eval $OP(m::AbstractArray{T}, args...; kwargs...) where {T <: Sym} = $OP.(m, args...; kwargs...)
    @eval $OP(m::Array{T}, args...; kwargs...) where {T <: Sym} = $OP.(m, args...; kwargs...)
    @eval $OP(m::Symmetric{T}, args...; kwargs...) where {T <: Sym} = Symmetric($OP.(m, args...; kwargs...))
end

tsimplify(x::Num, args...; kwargs...) = Symbolics.simplify(x, args...; kwargs...)
tsubs(x::Num, d...) = substitute(x, Dict(d...))
function tdiff(y::Num, x...; kwargs...)
    for xᵢ in x
        y = Symbolics.Differential(xᵢ)(y)
    end
    return expand_derivatives(y)
end
for OP in (:(tsimplify), :(tsubs), :(tdiff))
    @eval $OP(m::AbstractArray{Num}, args...; kwargs...) = $OP.(m, args...; kwargs...)
    @eval $OP(m::Array{Num}, args...; kwargs...) = $OP.(m, args...; kwargs...)
    @eval $OP(m::Symmetric{Num}, args...; kwargs...) = Symmetric($OP.(m, args...; kwargs...))
end

const SymType = Union{Sym, Num}

isdiagonal(a::AbstractMatrix{T}) where {T <: SymType} = isdiag(a)
isidentity(a::AbstractMatrix{T}) where {T <: SymType} = isone(a)

@inline LinearAlgebra.issymmetric(t::Tensor{2, 2, T}) where {T <: Union{AbstractFloat, Complex{AbstractFloat}}} = @inbounds t[1, 2] ≈ t[2, 1]

@inline function LinearAlgebra.issymmetric(t::Tensor{2, 3, T}) where {T <: Union{AbstractFloat, Complex{AbstractFloat}}}
    return @inbounds t[1, 2] ≈ t[2, 1] && t[1, 3] ≈ t[3, 1] && t[2, 3] ≈ t[3, 2]
end

@inline LinearAlgebra.issymmetric(t::Tensor{2, 2, Num}) = @inbounds iszero(t[1, 2] - t[2, 1])

@inline function LinearAlgebra.issymmetric(t::Tensor{2, 3, Num})
    return @inbounds iszero(t[1, 2] - t[2, 1]) && iszero(t[1, 3] - t[3, 1]) && iszero(t[2, 3] - t[3, 2])
end

function Tensors.isminorsymmetric(t::Tensor{4, dim, T}) where {dim, T <: Union{AbstractFloat, Complex{AbstractFloat}}}
    @inbounds for l in 1:dim, k in l:dim, j in 1:dim, i in j:dim
        if !(t[i, j, k, l] ≈ t[j, i, k, l]) || !(t[i, j, k, l] ≈ t[i, j, l, k])
            return false
        end
    end
    return true
end

function Tensors.isminorsymmetric(t::Tensor{4, dim, Num}) where {dim}
    @inbounds for l in 1:dim, k in l:dim, j in 1:dim, i in j:dim
        if !iszero(t[i, j, k, l] - t[j, i, k, l]) || !iszero(t[i, j, k, l] - t[i, j, l, k])
            return false
        end
    end
    return true
end

function Tensors.ismajorsymmetric(t::FourthOrderTensor{dim, T}) where {dim, T <: Union{AbstractFloat, Complex{AbstractFloat}}}
    @inbounds for l in 1:dim, k in l:dim, j in 1:dim, i in j:dim
        if !(t[i, j, k, l] ≈ t[k, l, i, j])
            return false
        end
    end
    return true
end

function Tensors.ismajorsymmetric(t::FourthOrderTensor{dim, Num}) where {dim}
    @inbounds for l in 1:dim, k in l:dim, j in 1:dim, i in j:dim
        if !(t[i, j, k, l] - t[k, l, i, j] == zero(Num))
            return false
        end
    end
    return true
end

@inline function Tensors.majortranspose(S::SymmetricTensor{4, dim}) where {dim}
    return SymmetricTensor{4, dim}(
        @inline function (i, j, k, l)
            return @inbounds S[k, l, i, j]
        end
    )
end

function Tensors.otimes(
        t1::AbstractArray{T1, order1},
        t2::AbstractArray{T2, order2},
    ) where {T1, T2, order1, order2}
    ec1 = ntuple(i -> i, order1)
    ec2 = ntuple(i -> order1 + i, order2)
    ec3 = ntuple(i -> i, order1 + order2)
    return einsum(EinCode((ec1, ec2), ec3), (AbstractArray{T1}(t1), AbstractArray{T2}(t2)))
end

function contract(t::AbstractArray{T, order}, i::Integer, j::Integer) where {T, order}
    m = min(i, j)
    M = max(i, j)
    ec1 = ntuple(k -> k == j ? i : k, order)
    ec2 = (Tuple(1:(m - 1))..., Tuple((m + 1):(M - 1))..., Tuple((M + 1):order)...)
    return einsum(EinCode((ec1,), ec2), (AbstractArray{T}(t),))
end

contract(t::AbstractArray{T, 2}, ::Integer, ::Integer) where {T} = tr(t)

function Tensors.dcontract(
        t1::AbstractArray{T1, order1},
        t2::AbstractArray{T2, order2},
    ) where {T1, T2, order1, order2}
    newc = order1 + order2
    ec1 = (ntuple(i -> i, order1 - 2)..., newc, newc + 1)
    ec2 = (newc, newc + 1, ntuple(i -> order1 - 2 + i, order2 - 2)...)
    ec3 = ntuple(i -> i, order1 + order2 - 4)
    return einsum(EinCode((ec1, ec2), ec3), (AbstractArray{T1}(t1), AbstractArray{T2}(t2)))
end

Tensors.dcontract(t1::AbstractArray{T1, 2}, t2::AbstractArray{T2, 2}) where {T1, T2} =
    dot(AbstractArray{T1}(t1), AbstractArray{T2}(t2))

function Tensors.dotdot(
        v1::AbstractArray{T1, order1},
        S::AbstractArray{TS, orderS},
        v2::AbstractArray{T2, order2},
    ) where {T1, TS, T2, order1, orderS, order2}
    newc = order1 + orderS
    ec1 = (ntuple(i -> i, order1 - 1)..., newc)
    ecS = (newc, ntuple(i -> order1 - 1 + i, orderS - 1)...)
    ec3 = ntuple(i -> i, order1 + orderS - 2)
    v1S = einsum(EinCode((ec1, ecS), ec3), (AbstractArray{T1}(v1), AbstractArray{TS}(S)))
    newc += order2
    ecv1S = (ntuple(i -> i, order1 + orderS - 3)..., newc)
    ec2 = (newc, ntuple(i -> order1 + orderS - 3 + i, order2 - 1)...)
    ec3 = ntuple(i -> i, newc - 4)
    return einsum(EinCode((ecv1S, ec2), ec3), (v1S, AbstractArray{T2}(v2)))
end

function qcontract(
        t1::AbstractArray{T1, order1},
        t2::AbstractArray{T2, order2},
    ) where {T1, T2, order1, order2}
    newc = order1 + order2
    ec1 = (ntuple(i -> i, order1 - 4)..., newc, newc + 1, newc + 2, newc + 3)
    ec2 = (newc, newc + 1, newc + 2, newc + 3, ntuple(i -> order1 - 4 + i, order2 - 4)...)
    ec3 = ntuple(i -> i, order1 + order2 - 8)
    return einsum(EinCode((ec1, ec2), ec3), (AbstractArray{T1}(t1), AbstractArray{T2}(t2)))
end

qcontract(t1::AbstractArray{T1, 4}, t2::AbstractArray{T2, 4}) where {T1, T2} =
    dot(AbstractArray{T1}(t1), AbstractArray{T2}(t2))

function Tensors.otimesu(
        t1::AbstractArray{T1, order1},
        t2::AbstractArray{T2, order2},
    ) where {T1, T2, order1, order2}
    ec1 = (ntuple(i -> i, order1 - 1)..., order1 + 1)
    ec2 = (order1, ntuple(i -> order1 + 1 + i, order2 - 1)...)
    ec3 = ntuple(i -> i, order1 + order2)
    return einsum(EinCode((ec1, ec2), ec3), (AbstractArray{T1}(t1), AbstractArray{T2}(t2)))
end

function Tensors.otimesl(
        t1::AbstractArray{T1, order1},
        t2::AbstractArray{T2, order2},
    ) where {T1, T2, order1, order2}
    ec1 = (ntuple(i -> i, order1 - 1)..., order1 + 2)
    ec2 = (order1, order1 + 1, ntuple(i -> order1 + 2 + i, order2 - 2)...)
    ec3 = ntuple(i -> i, order1 + order2)
    return einsum(EinCode((ec1, ec2), ec3), (AbstractArray{T1}(t1), AbstractArray{T2}(t2)))
end

otimesul(t1::AbstractArray{T1}, t2::AbstractArray{T2}) where {T1, T2} =
    (otimesu(t1, t2) + otimesl(t1, t2)) / promote_type(T1, T2)(2)

function sotimes(
        t1::AbstractArray{T1, order1},
        t2::AbstractArray{T2, order2},
    ) where {T1, T2, order1, order2}
    ec1 = ntuple(i -> i, order1)
    ec2 = ntuple(i -> order1 + i, order2)
    ec3 = ntuple(i -> i, order1 + order2)
    t3 = einsum(EinCode((ec1, ec2), ec3), (AbstractArray{T1}(t1), AbstractArray{T2}(t2)))
    ec1 = (ntuple(i -> i, order1 - 1)..., order1 + 1)
    ec2 = (order1, ntuple(i -> order1 + 1 + i, order2 - 1)...)
    ec3 = ntuple(i -> i, order1 + order2)
    t4 = einsum(EinCode((ec1, ec2), ec3), (AbstractArray{T1}(t1), AbstractArray{T2}(t2)))
    return (t3 + t4) / promote_type(T1, T2)(2)
end

@inline function sotimes(S1::Vec{dim}, S2::Vec{dim}) where {dim}
    return SymmetricTensor{2, dim}(
        @inline function (i, j)
            return @inbounds (S1[i] * S2[j] + S1[j] * S2[i]) / 2
        end
    )
end

@inline function sotimes(S1::SecondOrderTensor{dim}, S2::SecondOrderTensor{dim}) where {dim}
    TensorType = Tensors.getreturntype(
        otimes,
        Tensors.get_base(typeof(S1)),
        Tensors.get_base(typeof(S2)),
    )
    return TensorType(
        @inline function (i, j, k, l)
            return @inbounds (S1[i, j] * S2[k, l] + S1[i, k] * S2[j, l]) / 2
        end
    )
end

Tensors.otimes(α::Number, t::AbstractArray) = α * t
Tensors.otimes(t::AbstractArray, α::Number) = α * t

sotimes(α::Number, t::AbstractArray) = α * t
sotimes(t::AbstractArray, α::Number) = α * t

const ⊙ = qcontract
const ⊠ = otimesu
const ⊠ˢ = otimesul
const ⊗ˢ = sotimes
const sboxtimes = otimesul

export isidentity, contract, qcontract, otimesu, otimesul, sboxtimes, sotimes, ⊙, ⊠, ⊠ˢ, ⊗ˢ
export tsimplify, tfactor, tsubs, tdiff, ttrigsimp, texpand_trig
export ⋅, ⊡, ⊗
