struct TensISO{order, dim, T, N} <: AbstractTens{order, dim, T}
    data::NTuple{N, T}
    TensISO{dim}(λ::T) where {dim, T} = new{2, dim, T, 1}((λ,))
    TensISO{dim}(α::T1, β::T2) where {dim, T1, T2} = new{4, dim, promote_type(T1, T2), 2}((α, β))
    TensISO{dim}(data::NTuple{N, T}) where {dim, N, T} = TensISO{dim}(data...)
    TensISO{order, dim, T}() where {order, dim, T} =
        new{order, dim, T, order ÷ 2}(ntuple(_ -> one(T), Val(order ÷ 2)))
end

@pure get_order(::TensISO{order}) where {order} = order
@pure get_dim(::TensISO{order, dim}) where {order, dim} = dim
@pure Base.eltype(::Type{TensISO{order, dim, T}}) where {order, dim, T} = T
@pure Base.length(::TensISO{order, dim, T, N}) where {order, dim, T, N} = dim^order
@pure datanumber(::TensISO{order, dim, T, N}) where {order, dim, T, N} = N
@pure Base.size(::TensISO{order, dim}) where {order, dim} = ntuple(_ -> dim, Val(order))
Base.getindex(t::TensISO{2}, i::Integer, j::Integer) = t.data[1] * I[i, j]
Base.getindex(
    t::TensISO{4, dim},
    i::Integer,
    j::Integer,
    k::Integer,
    l::Integer,
) where {dim} =
    (t.data[1] - t.data[2]) * I[i, j] * I[k, l] / dim +
    t.data[2] * (I[i, k] * I[j, l] + I[i, l] * I[j, k]) / 2
function Base.replace_in_print_matrix(
        ::TensISO{2},
        i::Integer,
        j::Integer,
        s::AbstractString,
    )
    return i == j ? s : Base.replace_with_centered_mark(s)
end

"""
    tens_Id2(::Val{dim}, ::Val{T}) where {dim,T<:Number}

Identity tensor of second order `𝟏ᵢⱼ = δᵢⱼ = 1 if i=j otherwise 0`

# Examples
```julia
julia> 𝟏 = t𝟏() ; KM(𝟏)
6-element Vector{Sym}:
 1
 1
 1
 0
 0
 0

julia> 𝟏.data
3×3 SymmetricTensor{2, 3, Sym, 6}:
 1  0  0
 0  1  0
 0  0  1
```
"""
tens_Id2(::Val{dim} = Val(3), ::Val{T} = Val(Sym)) where {dim, T <: Number} = TensISO{2, dim, T}()

"""
    tens_Id4(::Val{dim} = Val(3), ::Val{T} = Val(Sym))

Symmetric identity tensor of fourth order  `𝕀 = 𝟏 ⊠ˢ 𝟏` i.e. `(𝕀)ᵢⱼₖₗ = (δᵢₖδⱼₗ+δᵢₗδⱼₖ)/2`

# Examples
```julia
julia> 𝕀 = t𝕀() ; KM(𝕀)
6×6 Matrix{Sym}:
 1  0  0  0  0  0
 0  1  0  0  0  0
 0  0  1  0  0  0
 0  0  0  1  0  0
 0  0  0  0  1  0
 0  0  0  0  0  1
```
"""
tens_Id4(::Val{dim} = Val(3), ::Val{T} = Val(Sym)) where {dim, T <: Number} = TensISO{4, dim, T}()

"""
    tens_J4(::Val{dim} = Val(3), ::Val{T} = Val(Sym))

Spherical projector of fourth order  `𝕁 = (𝟏 ⊗ 𝟏) / dim` i.e. `(𝕁)ᵢⱼₖₗ = δᵢⱼδₖₗ/dim`

# Examples
```julia
julia> 𝕁 = t𝕁() ; KM(𝕁)
6×6 Matrix{Sym}:
 1/3  1/3  1/3  0  0  0
 1/3  1/3  1/3  0  0  0
 1/3  1/3  1/3  0  0  0
   0    0    0  0  0  0
   0    0    0  0  0  0
   0    0    0  0  0  0
```
"""
tens_J4(::Val{dim} = Val(3), ::Val{T} = Val(Sym)) where {dim, T <: Number} =
    TensISO{dim}(one(T), zero(T))

"""
    tens_K4(::Val{dim} = Val(3), ::Val{T} = Val(Sym))

Deviatoric projector of fourth order  `𝕂 = 𝕀 - 𝕁` i.e. `(𝕂)ᵢⱼₖₗ = (δᵢₖδⱼₗ+δᵢₗδⱼₖ)/2 - δᵢⱼδₖₗ/dim`

# Examples
```julia
julia> 𝕂 = t𝕂() ; KM(𝕂)
6×6 Matrix{Sym}:
  2/3  -1/3  -1/3  0  0  0
 -1/3   2/3  -1/3  0  0  0
 -1/3  -1/3   2/3  0  0  0
    0     0     0  1  0  0
    0     0     0  0  1  0
    0     0     0  0  0  1
```
"""
tens_K4(::Val{dim} = Val(3), ::Val{T} = Val(Sym)) where {dim, T <: Number} =
    TensISO{dim}(zero(T), one(T))


"""
    iso_projectors(::Val{dim} = Val(3), ::Val{T} = Val(Sym))

Return the three fourth-order isotropic tensors `(𝕀, 𝕁, 𝕂)` — the symmetric
identity, spherical projector, and deviatoric projector.  Any isotropic
4th-order tensor can be written as `α·𝕁 + β·𝕂`.

# Examples
```julia
julia> 𝕀, 𝕁, 𝕂 = iso_projectors();

julia> 𝕁 + 𝕂 == 𝕀
true
```

See also [`tens_Id4`](@ref), [`tens_J4`](@ref), [`tens_K4`](@ref).
"""
iso_projectors(::Val{dim} = Val(3), ::Val{T} = Val(Sym)) where {dim, T <: Number} =
    tens_Id4(Val(dim), Val(T)), tens_J4(Val(dim), Val(T)), tens_K4(Val(dim), Val(T))

"""
    ISO(args...)

Legacy alias of [`iso_projectors`](@ref) kept for backward compatibility.
"""
ISO(args...) = iso_projectors(Val.(args)...)

for FUNC in (:tens_Id2, :tens_Id4, :tens_J4, :tens_K4, :iso_projectors)
    @eval $FUNC(args...) = $FUNC(Val.(args)...)
end

get_data(t::TensISO) = t.data
get_array(t::TensISO) = Array(t)
get_basis(::TensISO{order, dim, T}) where {order, dim, T} = CanonicalBasis{dim, T}()
get_var(::TensISO{order}) where {order} = ntuple(_ -> :cont, Val(order))
get_var(::TensISO, i::Integer) = :cont
components(t::TensISO{order, dim, T}) where {order, dim, T} = get_array(t)
components(
    t::TensISO{order, dim, T},
    ::OrthonormalBasis{dim, T},
    ::NTuple{order, Symbol},
) where {order, dim, T} = get_array(t)
components(t::TensISO{order, dim, T}, ::NTuple{order, Symbol}) where {order, dim, T} =
    get_array(t)

change_tens(t::TensISO{order, dim, T}, ℬ::OrthonormalBasis{dim, T}) where {order, dim, T} =
    Tens(get_array(t), ℬ)
change_tens(
    t::TensISO{order, dim, T},
    ℬ::OrthonormalBasis{dim, T},
    ::NTuple{order, Symbol},
) where {order, dim, T} = Tens(get_array(t), ℬ)

# Scalar arithmetic (-, α*A, A*α, A/α) defined in structured_tens_ops.jl
for OP in (:+, :-, :*)
    @eval @inline Base.$OP(
        A1::TensISO{order, dim},
        A2::TensISO{order, dim},
    ) where {order, dim} = TensISO{dim}($OP.(get_data(A1), get_data(A2)))
    @eval @inline Base.$OP(
        A1::TensISO{order, dim, T, N},
        A2::UniformScaling,
    ) where {order, dim, T, N} = TensISO{dim}($OP.(get_data(A1), ntuple(_ -> A2.λ, N)))
    @eval @inline Base.$OP(
        A1::UniformScaling,
        A2::TensISO{order, dim, T, N},
    ) where {order, dim, T, N} = TensISO{dim}($OP.(ntuple(_ -> A1.λ, N), get_data(A2)))
    @eval @inline function Base.$OP(
            A1::TensISO{order, dim, T, N},
            A2::AbstractTens{order, dim},
        ) where {order, dim, T, N}
        m1 = components(A1, get_basis(A2), get_var(A2))
        return Tens($OP(m1, get_array(A2)), get_basis(A2), get_var(A2))
    end
    @eval @inline function Base.$OP(
            A1::AbstractTens{order, dim},
            A2::TensISO{order, dim, T, N},
        ) where {order, dim, T, N}
        m2 = components(A2, get_basis(A1), get_var(A1))
        return Tens($OP(get_array(A1), m2), get_basis(A1), get_var(A1))
    end
end
for OP in (:(==), :(<=), :(>=), :(<), :(>))
    @eval @inline Base.$OP(
        A1::TensISO{order, dim},
        A2::TensISO{order, dim},
    ) where {order, dim} = all($OP.(get_data(A1), get_data(A2)))
end
@inline Base.inv(A::TensISO{order, dim, T}) where {order, dim, T} =
    TensISO{dim}(one(T) ./ get_data(A))
@inline Base.one(A::TensISO{order, dim, T}) where {order, dim, T} =
    TensISO{dim}(one.(get_data(A)))
@inline Base.zero(A::TensISO{order, dim, T}) where {order, dim, T} =
    TensISO{dim}(zero.(get_data(A)))

for FUNC in (:one, :zero)
    @eval begin
        @inline Base.$FUNC(A::AbstractTens{4, dim, T}) where {dim, T} =
            TensISO{dim}($FUNC(T), $FUNC(T))
        @inline Base.$FUNC(A::AbstractTens{2, dim, T}) where {dim, T} = TensISO{dim}($FUNC(T))
    end
end

@inline Base.literal_pow(::typeof(^), A::TensISO, ::Val{-1}) = inv(A)
@inline Base.literal_pow(::typeof(^), A::TensISO, ::Val{0}) = one(A)
@inline Base.literal_pow(::typeof(^), A::TensISO, ::Val{1}) = A
@inline Base.literal_pow(
    ::typeof(^),
    A::TensISO{order, dim, T},
    ::Val{p},
) where {order, dim, T, p} = TensISO{dim}(get_data(A) .^ (p))

@inline Base.transpose(A::TensISO) = A
@inline Base.adjoint(A::TensISO) = A

# ── Display ───────────────────────────────────────────────────────────────────

function Base.show(io::IO, A::TensISO{4})
    return print(io, "(", get_data(A)[1], ") 𝕁 + (", get_data(A)[2], ") 𝕂")
end
function Base.show(io::IO, A::TensISO{2})
    return print(io, "(", get_data(A)[1], ") 𝟏")
end

for OP in (:show, :print, :display)
    @eval function Base.$OP(A::TensISO{4})
        $OP(typeof(A))
        print("→ decomposition: ")
        println("(", get_data(A)[1], ") 𝕁 + (", get_data(A)[2], ") 𝕂")
        print("→ KM: ")
        return $OP(KM(A))
    end
    @eval function Base.$OP(A::TensISO{2})
        $OP(typeof(A))
        print("→ decomposition: ")
        return println("(", get_data(A)[1], ") 𝟏")
    end
end

intrinsic(A::TensISO{4}) = println("(", get_data(A)[1], ") 𝕁 + (", get_data(A)[2], ") 𝕂")
intrinsic(A::TensISO{2}) = println("(", get_data(A)[1], ") 𝟏")

# ── Rebuild helper (used by symbolic ops) ─────────────────────────────────────

_rebuild(::TensISO{order, dim}, new_data) where {order, dim} = TensISO{dim}(new_data)

# Symbolic helpers (tsimplify, tsubs, …) defined in structured_tens_ops.jl


"""
    KM(v::AllIsotropic{dim}; kwargs...)

Kelvin-Mandel vector or matrix representation
"""
KM(A::TensISO{order, dim}) where {order, dim} = tomandel(SymmetricTensor{order, dim}(A))

Tensors.otimes(A::TensISO{2, dim}, B::TensISO{2, dim}) where {dim} =
    TensISO{dim}(dim * get_data(A)[1] * get_data(B)[1], zero(eltype(A)))

scontract(A::TensISO{2, dim}, B::TensISO{2, dim}) where {dim} =
    TensISO{dim}(get_data(A)[1] * get_data(B)[1])

Tensors.otimes(A::TensISO{2, dim}) where {dim} =
    TensISO{dim}(dim * get_data(A)[1]^2, zero(eltype(A)))

scontract(A::TensISO{2, dim}) where {dim} =
    TensISO{dim}(get_data(A)[1]^2)

scontract(A::TensISO{2, dim}, B::AbstractArray) where {dim} = get_data(A)[1] * B
scontract(A::AbstractArray, B::TensISO{2, dim}) where {dim} = A * get_data(B)[1]

LinearAlgebra.dot(A::TensISO{2, dim}, B::TensISO{2, dim}) where {dim} = scontract(A, B)
for T in (AbstractArray, AbstractTens)
    @eval LinearAlgebra.dot(A::TensISO{2, dim}, B::$T) where {dim} = scontract(A, B)
    @eval LinearAlgebra.dot(A::$T, B::TensISO{2, dim}) where {dim} = scontract(A, B)
end

Tensors.dcontract(A::TensISO{2, dim}, B::TensISO{2, dim}) where {dim} =
    dim * get_data(A)[1] * get_data(B)[1]
Tensors.dcontract(A::TensISO{4, dim}, B::TensISO{2, dim}) where {dim} =
    TensISO{dim}(get_data(A)[1] * get_data(B)[1])
Tensors.dcontract(A::TensISO{2, dim}, B::TensISO{4, dim}) where {dim} =
    TensISO{dim}(get_data(A)[1] * get_data(B)[1])
Tensors.dcontract(A::TensISO{4, dim}, B::TensISO{4, dim}) where {dim} =
    TensISO{dim}(get_data(A)[1] * get_data(B)[1], get_data(A)[2] * get_data(B)[2])

Tensors.dcontract(A::TensISO{2, dim}, B::AbstractTens{order, dim}) where {order, dim} = get_data(A)[1] * contract(B, 1, 2)
Tensors.dcontract(A::AbstractTens{order, dim}, B::TensISO{2, dim}) where {order, dim} = contract(A, order - 1, order) * get_data(B)[1]

Tensors.dcontract(A::TensISO{4, dim}, B::TensOrthonormal{2}) where {dim} =
    get_data(A)[2] * B + (get_data(A)[1] - get_data(A)[2]) * tr(B) * I / dim
Tensors.dcontract(A::TensOrthonormal{2}, B::TensISO{4, dim}) where {dim} =
    A * get_data(B)[2] + tr(A) * (get_data(B)[1] - get_data(B)[2]) * I / dim

function Tensors.dcontract(
        A::TensISO{4, dim, T},
        B::AllTensOrthogonal{order, dim},
    ) where {order, dim, T}
    nB = TensOrthonormal(B)
    m = get_array(nB)
    ec1 = ntuple(i -> i, order)
    ec2 = (2, 1, ntuple(i -> i + 2, order - 2)...)
    m2 = einsum(EinCode((ec1,), ec2), (m,))
    newm =
        get_data(A)[2] * (m + m2) / 2 +
        (get_data(A)[1] - get_data(A)[2]) * Id2{dim, T}() ⊗ contract(m, 1, 2) / dim
    return Tens(newm, get_basis(nB))
end

function Tensors.dcontract(
        A::AllTensOrthogonal{order, dim},
        B::TensISO{4, dim, T},
    ) where {order, dim, T}
    nA = TensOrthonormal(A)
    m = get_array(nA)
    ec1 = ntuple(i -> i, order)
    ec2 = (ntuple(i -> i, order - 2)..., order, order - 1)
    m2 = einsum(EinCode((ec1,), ec2), (m,))
    newm =
        (m + m2) * get_data(B)[2] / 2 +
        contract(m, order - 1, order) ⊗ Id2{dim, T}() * (get_data(B)[1] - get_data(B)[2]) / dim
    return Tens(newm, get_basis(nA))
end

for order in (2, 4)
    for OP in (:+, :-, :*)
        @eval @inline Base.$OP(
            A1::AbstractTensor{$order, dim, T},
            A2::UniformScaling{T},
        ) where {dim, T <: Number} = $OP(A1, A2.λ * one(A1))
        @eval @inline Base.$OP(
            A1::UniformScaling{T},
            A2::AbstractTensor{$order, dim, T},
        ) where {dim, T <: Number} = $OP(A1.λ * one(A2), A2)
        @eval @inline Base.$OP(
            A1::AbstractTensor{$order, dim, T},
            A2::UniformScaling{T},
        ) where {dim, T <: SymType} = $OP(A1, A2.λ * one(A1))
        @eval @inline Base.$OP(
            A1::UniformScaling{T},
            A2::AbstractTensor{$order, dim, T},
        ) where {dim, T <: SymType} = $OP(A1.λ * one(A2), A2)
    end
end

Tensors.dotdot(v1::AbstractTens{1}, S::TensISO{2, dim}, v2::AbstractTens{1}) where {dim} =
    get_data(S)[1] * v1 ⋅ v2
Tensors.dotdot(v1::AbstractTens{1}, S::TensISO{4, dim}, v2::AbstractTens{1}) where {dim} =
    (get_data(S)[1] - get_data(S)[2]) * (v1 ⊗ v2) / dim +
    get_data(S)[2] * (v2 ⊗ v1 + v1 ⋅ v2 * I) / 2

Tensors.dotdot(a1::AbstractTens{2}, S::TensISO{4, dim}, a2::AbstractTens{2}) where {dim} =
    (get_data(S)[1] - get_data(S)[2]) * tr(a1) * tr(a2) / dim + get_data(S)[2] * a1 ⊡ a2

qcontract(A::TensISO{4, 3}, B::TensISO{4, 3}) =
    get_data(A)[1] * get_data(B)[1] + 5 * get_data(A)[2] * get_data(B)[2]

qcontract(A::TensISO{4, 2}, B::TensISO{4, 2}) =
    get_data(A)[1] * get_data(B)[1] + 2 * get_data(A)[2] * get_data(B)[2]

function qcontract(A::TensISO{4, dim, T}, B::AllTensOrthogonal{order, dim}) where {order, dim, T}
    nB = TensOrthonormal(B)
    m = get_array(nB)
    newm =
        get_data(A)[2] *
        (contract(contract(m, 1, 3), 1, 2) + contract(contract(m, 1, 4), 1, 2)) / 2 +
        (get_data(A)[1] - get_data(A)[2]) * contract(contract(m, 1, 2), 1, 2) / dim
    return Tens(newm, get_basis(nB))
end

function qcontract(A::AllTensOrthogonal{order, dim}, B::TensISO{4, dim, T}) where {order, dim, T}
    nA = TensOrthonormal(A)
    m = get_array(nA)
    newm =
        (
        contract(contract(m, order - 2, order), order - 1, order) +
            contract(contract(m, order - 3, order), order - 1, order)
    ) * get_data(B)[2] / 2 +
        contract(contract(m, order - 1, order), order - 1, order) *
        (get_data(B)[1] - get_data(B)[2]) / dim
    return Tens(newm, get_basis(nA))
end

isotropify(A::AbstractArray{T, 2}) where {T} = TensISO{size(A)[1]}(tr(A) / size(A)[1])

function isotropify(A::AbstractArray{T, 4}) where {T}
    dim = size(A)[1]
    α = tens_J4(dim, T) ⊙ A
    # 5 = dim(deviatoric space) = dim(𝕂) for 3D; generalises to dim*(dim+1)/2 - 1
    dimK = dim * (dim + 1) / 2 - 1
    β = (tens_K4(dim, T) ⊙ A) / dimK
    return TensISO{dim}(α, β)
end

TensISO(A::AbstractArray) = isotropify(Tens(A))

function proj_tens(::Val{:ISO}, A::AbstractArray)
    norm = x -> tsimplify(√(sum(x .^ 2)))
    nA = norm(A)
    if nA == zero(eltype(A))
        return zero(A), nA, nA
    else
        B = isotropify(A)
        d = norm(B - A)
        return B, d, d / nA
    end
end

is_ISO(::TensISO) = true
# is_ISO(A::AbstractArray; ε = 1.0e-6) is defined in tens_projection.jl using
# the closed-form projection (O(1)) and a relative-residual tolerance.
is_TI(::TensISO) = false
is_ORTHO(::TensISO) = false

LinearAlgebra.issymmetric(::TensISO) = true
Tensors.isminorsymmetric(::TensISO{4}) = true
Tensors.ismajorsymmetric(::TensISO{4}) = true

export TensISO, tens_Id2, tens_Id4, tens_J4, tens_K4, ISO, iso_projectors, isotropify, is_ISO, is_TI, is_ORTHO
