##############################################################################
# Lazy tensors: the zero tensor and the multiples of the identity.            #
#                                                                            #
# `TensZero{order, dim, T}` stores nothing and `TensIdentity{order, dim, T}`  #
# stores its factor λ alone: λ 𝟏 at the order 2, λ 𝕀 at the order 4, 𝕀 the  #
# symmetric identity (`𝕀ᵢⱼₖₗ = (δᵢₖδⱼₗ + δᵢₗδⱼₖ)/2`, `tens_Id4`). Both are     #
# tensors like the others (components, bases, Kelvin–Mandel, display); what   #
# they add is that an operation on them reduces to its result without        #
# forming any component: the zero tensor is neutral for the sums and         #
# absorbing for the products and contractions, λ 𝟏 contracts as λ times the  #
# identity map and λ 𝕀 ⊡ 𝐀 = λ 𝐀 for a tensor with the minor symmetries.    #
# An operation without a closed form goes through `TensISO`, which both are  #
# exactly, and the promotions of the structured tensors.                     #
##############################################################################

"""
    TensZero{order, dim, T}()

The zero tensor of order `order` in dimension `dim`, of element type `T`,
which stores nothing. It is neutral for `+` and `-`, and the products and
contractions that involve it return the zero tensor of the order of their
result (a zero of type `T` at the order 0) without forming any component.

# Examples
```jldoctest
julia> 𝟎 = TensZero{2, 3, Float64}();

julia> A = TensISO{3}(2.0);

julia> 𝟎 + A === A
true

julia> 𝟎 ⊗ A
TensZero{4, 3, Float64}
→ lazy: 𝟎
```

See also [`TensIdentity`](@ref).
"""
struct TensZero{order, dim, T <: Number} <: AbstractTens{order, dim, T} end

"""
    TensIdentity{order, dim}(λ = true)

The multiple ``\\lambda\\,\\boldsymbol{1}`` (order 2) or ``\\lambda\\,\\mathbb{I}``
(order 4) of the identity in dimension `dim`, which stores `λ` alone;
``\\mathbb{I}`` is the symmetric identity,
``\\mathbb{I}_{ijkl} = (\\delta_{ik}\\delta_{jl} + \\delta_{il}\\delta_{jk})/2``, that of
[`tens_Id4`](@ref). Its contractions reduce to the identity map:
``\\lambda\\boldsymbol{1}\\cdot\\mathbf{A} = \\lambda\\mathbf{A}``,
``\\lambda\\mathbb{I}:\\mathbf{A} = \\lambda\\mathbf{A}`` for a tensor with the
minor symmetries (its symmetric part otherwise), and its inverse is
``\\lambda^{-1}\\boldsymbol{1}`` or ``\\lambda^{-1}\\mathbb{I}``.

# Examples
```jldoctest
julia> 𝕀 = TensIdentity{4, 3}(2.0);

julia> C = TensISO{3}(3.0, 2.0);

julia> 𝕀 ⊡ C == 2C
true

julia> inv(𝕀)
TensIdentity{4, 3, Float64}
→ lazy: (0.5) 𝕀
```

See also [`TensZero`](@ref), [`TensISO`](@ref).
"""
struct TensIdentity{order, dim, T <: Number} <: AbstractTens{order, dim, T}
    λ::T
    function TensIdentity{order, dim}(λ::T = true) where {order, dim, T <: Number}
        order in (2, 4) || throw(
            ArgumentError(
                "TensIdentity: the identity is of order 2 (𝟏) or 4 (𝕀, symmetric), got $order"
            )
        )
        return new{order, dim, T}(λ)
    end
end

const _LazyTens{order, dim, T} = Union{TensZero{order, dim, T}, TensIdentity{order, dim, T}}

# ── Interface of AbstractTens ────────────────────────────────────────────────

Base.size(::_LazyTens{order, dim}) where {order, dim} = ntuple(_ -> dim, Val(order))
Base.getindex(::TensZero{order, dim, T}, ::Vararg{Integer, order}) where {order, dim, T} = zero(T)
Base.getindex(t::TensIdentity{2}, i::Integer, j::Integer) = t.λ * δkron(typeof(t.λ), i, j)
Base.getindex(t::TensIdentity{4}, i::Integer, j::Integer, k::Integer, l::Integer) =
    t.λ * (δkron(typeof(t.λ), i, k) * δkron(typeof(t.λ), j, l) + δkron(typeof(t.λ), i, l) * δkron(typeof(t.λ), j, k)) / 2

get_data(::TensZero) = ()
get_data(t::TensIdentity) = (t.λ,)
get_array(t::_LazyTens) = Array(t)
get_basis(::_LazyTens{order, dim, T}) where {order, dim, T} = CanonicalBasis{dim, T}()
get_var(::_LazyTens{order}) where {order} = ntuple(_ -> :cont, Val(order))
get_var(::_LazyTens, i::Integer) = :cont
components(t::_LazyTens) = get_array(t)
# both are the same in every orthonormal basis, whatever the variances
components(t::_LazyTens{order, dim}, ::OrthonormalBasis{dim}, ::NTuple{order, Symbol}) where {order, dim} = get_array(t)
components(t::_LazyTens{order}, ::NTuple{order, Symbol}) where {order} = get_array(t)
express(t::_LazyTens{order, dim}, ℬ::OrthonormalBasis{dim}) where {order, dim} = Tens(_static_tensor(t), ℬ)
express(t::_LazyTens{order, dim}, ℬ::OrthonormalBasis{dim}, ::NTuple{order, Symbol}) where {order, dim} =
    Tens(_static_tensor(t), ℬ)

@inline _static_tensor(t::_LazyTens{order, dim}) where {order, dim} =
    order in (2, 4) && dim in (2, 3) ?
    SymmetricTensor{order, dim, eltype(t)}((ij...) -> @inbounds t[ij...]) :
    _static_tensor(get_array(t), Val(dim), Val(order))

# The isotropic tensor each one is.
TensISO(t::TensIdentity{2, dim}) where {dim} = TensISO{dim}(t.λ)
TensISO(t::TensIdentity{4, dim}) where {dim} = TensISO{dim}(t.λ, t.λ)
TensISO(::TensZero{2, dim, T}) where {dim, T} = TensISO{dim}(zero(T))
TensISO(::TensZero{4, dim, T}) where {dim, T} = TensISO{dim}(zero(T), zero(T))

isisotropic(::_LazyTens) = true
istransverselyisotropic(::_LazyTens) = false
isorthotropic(::_LazyTens) = false
LinearAlgebra.issymmetric(::_LazyTens{2}) = true
LinearAlgebra.issymmetric(::_LazyTens{4}) = true
Tensors.isminorsymmetric(::_LazyTens{4}) = true
Tensors.ismajorsymmetric(::_LazyTens{4}) = true
Base.iszero(::TensZero) = true
Base.iszero(t::TensIdentity) = iszero(t.λ)

# ── Kelvin–Mandel ────────────────────────────────────────────────────────────

_mandel(::TensZero{2, dim, T}) where {dim, T} = zeros(T, dim * (dim + 1) ÷ 2)
_mandel(::TensZero{4, dim, T}) where {dim, T} = zeros(T, dim * (dim + 1) ÷ 2, dim * (dim + 1) ÷ 2)
_mandel(t::TensIdentity{2, dim}) where {dim} = [i ≤ dim ? t.λ : zero(t.λ) for i in 1:(dim * (dim + 1) ÷ 2)]
_mandel(t::TensIdentity{4, dim}) where {dim} =
    [i == j ? t.λ : zero(t.λ) for i in 1:(dim * (dim + 1) ÷ 2), j in 1:(dim * (dim + 1) ÷ 2)]

# ── Display ──────────────────────────────────────────────────────────────────

Base.show(io::IO, ::TensZero) = print(io, "𝟎")
Base.show(io::IO, t::TensIdentity{2}) = print(io, "(", t.λ, ") 𝟏")
Base.show(io::IO, t::TensIdentity{4}) = print(io, "(", t.λ, ") 𝕀")
function Base.show(io::IO, ::MIME"text/plain", t::_LazyTens)
    println(io, typeof(t))
    print(io, "→ lazy: ")
    return show(io, t)
end

# ── Scalars, zero and one ────────────────────────────────────────────────────

Base.:-(t::TensZero) = t
Base.:-(t::TensIdentity{order, dim}) where {order, dim} = TensIdentity{order, dim}(-t.λ)
Base.:*(::Number, t::TensZero{order, dim, T}) where {order, dim, T} = TensZero{order, dim, T}()
Base.:*(t::TensZero, α::Number) = α * t
Base.:/(t::TensZero{order, dim, T}, α::Number) where {order, dim, T} = TensZero{order, dim, typeof(zero(T) / α)}()
Base.:*(α::Number, t::TensIdentity{order, dim}) where {order, dim} = TensIdentity{order, dim}(α * t.λ)
Base.:*(t::TensIdentity, α::Number) = α * t
Base.:/(t::TensIdentity{order, dim}, α::Number) where {order, dim} = TensIdentity{order, dim}(t.λ / α)
# against the division of a real tensor by a symbol
Base.:/(t::TensZero{order, dim, <:Real}, α::Symbolics.Num) where {order, dim} = TensZero{order, dim, typeof(zero(eltype(t)) / α)}()
Base.:/(t::TensIdentity{order, dim, <:Real}, α::Symbolics.Num) where {order, dim} = TensIdentity{order, dim}(t.λ / α)

Base.zero(::_LazyTens{order, dim, T}) where {order, dim, T} = TensZero{order, dim, T}()
Base.one(::_LazyTens{order, dim, T}) where {order, dim, T} = TensIdentity{order, dim}(one(T))
Base.transpose(t::_LazyTens) = t
Base.adjoint(t::TensZero) = t
Base.adjoint(t::TensIdentity{order, dim}) where {order, dim} = TensIdentity{order, dim}(conj(t.λ))

Base.inv(t::TensIdentity{order, dim}) where {order, dim} = TensIdentity{order, dim}(inv(t.λ))
Base.inv(t::TensZero) = throw(DomainError(t, "the zero tensor has no inverse"))
Base.literal_pow(::typeof(^), t::TensIdentity{order, dim}, ::Val{p}) where {order, dim, p} =
    TensIdentity{order, dim}(t.λ^p)
# the power −1 of λ 𝟏 against that of the matrices
Base.literal_pow(::typeof(^), t::TensIdentity{2}, ::Val{-1}) = inv(t)
# the positive powers of the zero tensor are zero, its power 0 is the identity
Base.literal_pow(::typeof(^), t::TensZero, ::Val{p}) where {p} = _zero_power(t, Val(p > 0), Val(p == 0))
Base.literal_pow(::typeof(^), t::TensZero{2}, ::Val{-1}) = inv(t)
_zero_power(t::TensZero, ::Val{true}, ::Val{false}) = t
_zero_power(t::TensZero, ::Val{false}, ::Val{true}) = one(t)
_zero_power(t::TensZero, ::Val{false}, ::Val{false}) = inv(t)

LinearAlgebra.tr(::TensZero{2, dim, T}) where {dim, T} = zero(T)
LinearAlgebra.tr(t::TensIdentity{2, dim}) where {dim} = dim * t.λ
LinearAlgebra.norm(::TensZero{order, dim, T}) where {order, dim, T} = abs(zero(T))
LinearAlgebra.norm(t::TensIdentity{2, dim}) where {dim} = abs(t.λ) * sqrt(oftype(float(abs(t.λ)), dim))
LinearAlgebra.norm(t::TensIdentity{4, dim}) where {dim} = abs(t.λ) * sqrt(oftype(float(abs(t.λ)), dim * (dim + 1) ÷ 2))

Base.:(==)(::TensZero{order, dim}, ::TensZero{order, dim}) where {order, dim} = true
Base.:(==)(a::TensIdentity{order, dim}, b::TensIdentity{order, dim}) where {order, dim} = a.λ == b.λ

# ── Sums ─────────────────────────────────────────────────────────────────────

# The zero tensor is neutral and two identities add their factors; any other
# sum with an identity goes through the isotropic tensor it is.
_lazy_add(::TensZero, B::AbstractTens) = B
_lazy_sub(::TensZero, B::AbstractTens) = -B
for (F, OP) in ((:_lazy_add, :+), (:_lazy_sub, :-))
    @eval begin
        $F(A::AbstractTens, ::TensZero) = A
        $F(a::TensIdentity, B::AbstractTens) = $OP(TensISO(a), B)
        $F(A::AbstractTens, b::TensIdentity) = $OP(A, TensISO(b))
        # the operands of a sum share their order and dimension
        $F(a::TensZero, b::TensZero) = TensZero{get_order(a), get_dim(a), promote_type(eltype(a), eltype(b))}()
        $F(::TensZero, b::TensIdentity) = $OP(b)
        $F(a::TensIdentity, ::TensZero) = a
        $F(a::TensIdentity, b::TensIdentity) = TensIdentity{get_order(a), get_dim(a)}($OP(a.λ, b.λ))
    end
end
# `+B` of a tensor is the tensor itself.
Base.:+(t::_LazyTens) = t

# ── Products and contractions ────────────────────────────────────────────────

# The zero of the order of a result: a scalar at the order 0.
_lazy_zero(::Val{0}, ::Val{dim}, ::Type{T}) where {dim, T} = zero(T)
_lazy_zero(::Val{order}, ::Val{dim}, ::Type{T}) where {order, dim, T} = TensZero{order, dim, T}()

# The order of `A` contracted with `B` over `n` indices, checked.
function _contracted_order(o1::Int, o2::Int, n::Int)
    (o1 ≥ n && o2 ≥ n) || throw(
        DimensionMismatch("a contraction over $n indices of tensors of orders $o1 and $o2")
    )
    return o1 + o2 - 2n
end

# The product of a zero tensor with any tensor, over n contracted indices.
_zero_product(n::Int, A::AbstractTens{o1, dim}, B::AbstractTens{o2, dim}) where {o1, o2, dim} =
    _lazy_zero(Val(_contracted_order(o1, o2, n)), Val(dim), promote_type(eltype(A), eltype(B)))

# λ 𝟏 is the identity map for the simple contraction, on either side.
_identity_dot(a::TensIdentity{2}, B::AbstractTens) = a.λ * B
_identity_dot(A::AbstractTens, b::TensIdentity{2}) = b.λ * A
_identity_dot(a::TensIdentity{2}, b::TensIdentity{2}) = TensIdentity{2, get_dim(a)}(a.λ * b.λ)

# λ 𝟏 ⊡ 𝐀 contracts the first two indices of 𝐀, 𝐀 ⊡ λ 𝟏 its last two; λ 𝕀 ⊡ 𝐀
# is λ times the symmetric part of 𝐀 in its first two indices, 𝐀 itself for a
# tensor that has that symmetry by its type, and the mirror image on the right.
_identity_ddot(a::TensIdentity{2}, B::AbstractTens) = a.λ * contract(B, 1, 2)
_identity_ddot(A::AbstractTens{order}, b::TensIdentity{2}) where {order} = b.λ * contract(A, order - 1, order)
_identity_ddot(a::TensIdentity{4}, B::AbstractTens) = a.λ * _sym_pair(B, Val(:first))
_identity_ddot(A::AbstractTens, b::TensIdentity{4}) = b.λ * _sym_pair(A, Val(:last))
# two identities: the product of their factors on the identity of the order of
# the result, the trace of 𝟏 for 𝟏 ⊡ 𝟏
_identity_ddot(a::TensIdentity{4}, b::TensIdentity{4}) = TensIdentity{4, get_dim(a)}(a.λ * b.λ)
_identity_ddot(a::TensIdentity{4}, b::TensIdentity{2}) = TensIdentity{2, get_dim(a)}(a.λ * b.λ)
_identity_ddot(a::TensIdentity{2}, b::TensIdentity{4}) = TensIdentity{2, get_dim(a)}(a.λ * b.λ)
_identity_ddot(a::TensIdentity{2}, b::TensIdentity{2}) = get_dim(a) * a.λ * b.λ

# The tensors whose type carries the symmetry of a pair of indices: the
# structured ones and the lazy ones (both pairs), and those stored as
# Tensors.jl symmetric tensors.
const _PairSymmetric{dim} = Union{
    TensISO{2, dim}, TensISO{4, dim}, TensTI{2}, TensTI{4}, TensOrtho, TensCubic, _LazyTens{2, dim}, _LazyTens{4, dim},
    TensOrthonormal{2, dim, <:Any, <:SymmetricTensor{2, dim}}, TensOrthonormal{4, dim, <:Any, <:SymmetricTensor{4, dim}},
}
_sym_pair(A::_PairSymmetric, ::Val) = A
function _sym_pair(A::AbstractTens{order, dim}, ::Val{side}) where {order, dim, side}
    # the canonical components, where the variances do not matter
    ℬ = CanonicalBasis{dim, eltype(A)}()
    a = components(A, ℬ, ntuple(_ -> :cont, Val(order)))
    p = side === :first ? (2, 1, ntuple(k -> k + 2, Val(order - 2))...) :
        (ntuple(identity, Val(order - 2))..., order, order - 1)
    return Tens((a .+ permutedims(a, p)) ./ 2, ℬ)
end

# ── Methods ──────────────────────────────────────────────────────────────────
#
# Each operation dispatches to the functions above for a lazy operand on
# either side. The generic methods take any other tensor; those on the
# structured types and on a few orders exist already for every operation, and
# the methods that pair a lazy operand with each of them settle what would
# otherwise be ambiguous.

const _LAZY_PARTNERS = (
    :(TensISO{2}), :(TensISO{4}), :(TensTI{2}), :(TensTI{4}), :(TensOrtho{2}), :(TensOrtho{4}), :TensCubic,
    :(AbstractTens{1}), :(AbstractTens{4}),
)

for (OP, F) in ((:+, :_lazy_add), (:-, :_lazy_sub))
    @eval begin
        Base.$OP(A::_LazyTens{order, dim}, B::AbstractTens{order, dim}) where {order, dim} = $F(A, B)
        Base.$OP(A::AbstractTens{order, dim}, B::_LazyTens{order, dim}) where {order, dim} = $F(A, B)
        Base.$OP(A::_LazyTens{order, dim}, B::_LazyTens{order, dim}) where {order, dim} = $F(A, B)
        Base.$OP(A::_LazyTens{order, dim}, B::TensISO{order, dim}) where {order, dim} = $F(A, B)
        Base.$OP(A::TensISO{order, dim}, B::_LazyTens{order, dim}) where {order, dim} = $F(A, B)
    end
end

for (OP, n) in (
        (:(LinearAlgebra.dot), 1), (:(Tensors.dcontract), 2), (:qcontract, 4),
        (:(Tensors.otimes), 0), (:(Tensors.otimesu), 0), (:otimesul, 0), (:sotimes, 0),
    )
    @eval begin
        $OP(A::TensZero{o1, dim}, B::AbstractTens{o2, dim}) where {o1, o2, dim} = _zero_product($n, A, B)
        $OP(A::AbstractTens{o1, dim}, B::TensZero{o2, dim}) where {o1, o2, dim} = _zero_product($n, A, B)
        $OP(A::TensZero, B::TensZero) = _zero_product($n, A, B)
    end
    for P in _LAZY_PARTNERS
        @eval begin
            $OP(A::TensZero, B::$P) = _zero_product($n, A, B)
            $OP(A::$P, B::TensZero) = _zero_product($n, A, B)
        end
    end
end

for (OP, F) in ((:(LinearAlgebra.dot), :_identity_dot), (:(Tensors.dcontract), :_identity_ddot))
    orders = OP == :(LinearAlgebra.dot) ? (2,) : (2, 4)
    for o in orders
        @eval begin
            $OP(A::TensIdentity{$o, dim}, B::AbstractTens{order, dim}) where {order, dim} = $F(A, B)
            $OP(A::AbstractTens{order, dim}, B::TensIdentity{$o, dim}) where {order, dim} = $F(A, B)
            $OP(A::TensIdentity{$o}, B::TensZero) = _zero_product($(OP == :(LinearAlgebra.dot) ? 1 : 2), A, B)
            $OP(A::TensZero, B::TensIdentity{$o}) = _zero_product($(OP == :(LinearAlgebra.dot) ? 1 : 2), A, B)
        end
        for P in _LAZY_PARTNERS
            @eval begin
                $OP(A::TensIdentity{$o}, B::$P) = $F(A, B)
                $OP(A::$P, B::TensIdentity{$o}) = $F(A, B)
            end
        end
    end
end
LinearAlgebra.dot(A::TensIdentity{2}, B::TensIdentity{2}) = _identity_dot(A, B)
for o1 in (2, 4), o2 in (2, 4)
    @eval Tensors.dcontract(A::TensIdentity{$o1}, B::TensIdentity{$o2}) = _identity_ddot(A, B)
end

export TensZero, TensIdentity
