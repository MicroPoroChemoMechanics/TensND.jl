δkron(T::Type{<:Number}, i::Integer, j::Integer) = i == j ? one(T) : zero(T)

struct Id2{dim, T <: Number} <: AbstractMatrix{T} end
@pure Base.size(::Id2{dim}) where {dim} = (dim, dim)
Base.getindex(::Id2{dim, T}, i::Integer, j::Integer) where {dim, T} = δkron(T, i, j)
function Base.replace_in_print_matrix(::Id2, i::Integer, j::Integer, s::AbstractString)
    return i == j ? s : Base.replace_with_centered_mark(s)
end

"""
    isidentity(a::AbstractMatrix) → Bool

Test whether `a` is the identity matrix, within the tolerance of `≈` for
floating-point element types and exactly for symbolic ones.

# Examples
```jldoctest
julia> isidentity([1.0 0.0; 0.0 1.0])
true

julia> isidentity([1.0 1e-20; 0.0 1.0])
true

julia> isidentity([1.0 0.5; 0.0 1.0])
false
```
"""
isidentity(a::AbstractMatrix{T}) where {T} = a ≈ I
isdiagonal(a::AbstractMatrix{T}) where {T} = norm(a - Diagonal(a)) <= eps(T)

# ── Symbolic passes ──────────────────────────────────────────────────────────
#
# Each of these applies a symbolic operation elementwise and is a **no-op on
# numeric types**, so generic code may call them unconditionally without paying
# anything on a `Float64` path. The fallback below is what provides that.

"""
    tsimplify(x, args...; kwargs...)

Simplify a scalar, an array or a tensor elementwise.

Dispatches to `SymPy.simplify` for `Sym`, to `Symbolics.simplify` for `Num`,
and is the **identity** on every other type — so a routine written to be
simplification-aware runs unchanged, and at full speed, on numeric input.

See also [`tfactor`](@ref), [`tsubs`](@ref), [`tdiff`](@ref),
[`ttrigsimp`](@ref), [`texpand_trig`](@ref).

# Examples
```jldoctest
julia> tsimplify(3.0)          # no-op on numbers
3.0

julia> tsimplify([1.0, 2.0])
2-element Vector{Float64}:
 1.0
 2.0
```
"""
function tsimplify end

"""
    tfactor(x, args...; kwargs...)

Factorize an expression elementwise (`SymPy.factor`); identity on non-symbolic
types. See [`tsimplify`](@ref).
"""
function tfactor end

"""
    tsubs(x, substitutions...)

Substitute into an expression elementwise (`SymPy.subs` / `Symbolics.substitute`);
identity on non-symbolic types.

Substitutions are given as `old => new` pairs. See [`tsimplify`](@ref).
"""
function tsubs end

"""
    tdiff(y, x...; kwargs...)

Differentiate an expression elementwise with respect to one or more symbols
(`SymPy.diff` / `Symbolics.Differential`); identity on non-symbolic types.

Repeated symbols give higher derivatives. See [`tsimplify`](@ref).
"""
function tdiff end

"""
    ttrigsimp(x, args...; kwargs...)

Trigonometric simplification elementwise (`sympy.trigsimp`); identity on
non-symbolic types. Often succeeds where [`tsimplify`](@ref) does not, on
expressions built from a rotated basis.
"""
function ttrigsimp end

"""
    texpand_trig(x, args...; kwargs...)

Trigonometric expansion elementwise (`sympy.expand_trig`); identity on
non-symbolic types. The counterpart of [`ttrigsimp`](@ref).
"""
function texpand_trig end

"""
    tlimit(x, s, v)
    tlimit(x, s, v, dir)

Limit of an expression elementwise as the symbol `s` tends to `v`
(`sympy.limit`), optionally one-sided with `dir = "+"` or `"-"`; identity on
non-symbolic types, a numeric value having no free symbol to send anywhere.

This is the pass a degenerate limit of a closed form needs — a vanishing
regularization, an incompressible phase (`k => oo`), a flat inclusion — and it
keeps the tensor's structure: applied to a `TensISO`, a `TensTI` or a
`TensOrtho` it takes the limit of the few canonical coefficients and rebuilds
the same type, never going through the `dim^order` components.

Unlike its siblings there is no `Symbolics.Num` method: Symbolics has no limit,
and returning the input unchanged would be a silent wrong answer, so a `Num`
raises instead.

See also [`tsimplify`](@ref), [`tsubs`](@ref).

# Examples
```jldoctest
julia> tlimit(3.0, :whatever, 0)      # no-op on numbers
3.0
```
"""
function tlimit end

for OP in (:(tsimplify), :(tfactor), :(tsubs), :(ttrigsimp), :(texpand_trig), :(tlimit))
    @eval $OP(x, args...; kwargs...) = x
end

# Tuples, elementwise.
#
# This is not cosmetic. `get_data` of every structured tensor returns an
# `NTuple` of coefficients, and `structured_tens_ops.jl` implements the whole
# family as `_rebuild(A, OP(get_data(A)))`. A `Tuple` is not an `AbstractArray`,
# so without these methods the generic identity fallback above caught it and
# `tsimplify(::TensISO{4,3,Sym})`, `tsubs(::TensTI{4,Sym,5}, …)` and their
# siblings were **silent no-ops** — they returned the tensor unchanged and
# reported success. Only the generic `Tens` path, whose data is an `Array`,
# ever simplified.
for OP in (
        :(tsimplify), :(tfactor), :(tsubs), :(tdiff),
        :(ttrigsimp), :(texpand_trig), :(tlimit),
    )
    @eval $OP(m::Tuple, args...; kwargs...) = map(x -> $OP(x, args...; kwargs...), m)
end

tlimit(x::T, args...; kwargs...) where {T <: Sym} = sympy.limit(x, args...; kwargs...)
tlimit(::Num, args...; kwargs...) = throw(
    ArgumentError(
        "tlimit: Symbolics.jl provides no limit. Use the SymPy backend " *
            "(`Sym`) for a limit passage, or substitute an explicit value with `tsubs`."
    )
)

for OP in (:(simplify), :(factor), :(subs), :(diff))
    @eval $(Symbol("t", OP))(x::T, args...; kwargs...) where {T <: Sym} = SymPy.$OP(x, args...; kwargs...)
end
for OP in (:(trigsimp), :(expand_trig))
    @eval $(Symbol("t", OP))(x::T, args...; kwargs...) where {T <: Sym} = sympy.$OP(x, args...; kwargs...)
end
for OP in (:(tsimplify), :(tfactor), :(tsubs), :(tdiff), :(ttrigsimp), :(texpand_trig), :(tlimit))
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

# Element types whose symmetry predicates should use a *tolerance* rather than
# exact equality: floating-point numbers, and `ForwardDiff.Dual` numbers built
# on them.
#
# `ForwardDiff.Dual <: Real` but **not** `<: AbstractFloat`, so before this
# union existed a `Dual`-valued tensor fell through to the exact fallback and a
# round-off of a few ulp was enough to report it as non-minor-symmetric. The
# consequence was not local: `_KM_of_array` then built a 9×9 matrix instead of
# a 6×6 one and every `proj_tens` call on a `Dual` input died with a
# `DimensionMismatch`, i.e. automatic differentiation through a projection was
# impossible. Symbolic types keep their own exact methods and are unaffected.
const ApproxType = Union{AbstractFloat, Complex{<:AbstractFloat}, ForwardDiff.Dual}

"""
    is_hard_numeric(::Type) -> Bool
    is_hard_numeric(x)      -> Bool

Whether comparisons on that scalar type yield an honest `Bool`, so a value of it
may drive `if`, `&&` or a tolerance test.

**`T <: Real` is not this predicate.** `Symbolics.Num <: Real`, yet `<`, `≈` and
`iszero` on a `Num` return *symbolic* expressions, and using one in a boolean
context throws `TypeError: non-boolean (Num) used in boolean context`.
`ForwardDiff.Dual` is also `<: Real` and *does* compare to a `Bool`, so the two
cannot be separated by `<: Real` either way. `SymPy.Sym` is not `<: Real` at
all, which is why the bug it guards against only ever showed up on the
Symbolics side.

The list is explicit and closed, and the default is the safe answer: an unknown
scalar type is assumed *not* comparable, which pushes callers onto their
structural (`isequal`) branch. Add a type here only after checking that its `<`
really returns a `Bool`.

This is the companion of [`ApproxType`](@ref): `ApproxType` says *use a
tolerance rather than exact equality*, `is_hard_numeric` says *a comparison is
allowed at all*. `Int` and `Rational` are hard numeric but not `ApproxType`.

```jldoctest
julia> TensND.is_hard_numeric(Float64), TensND.is_hard_numeric(Sym)
(true, false)
```
"""
is_hard_numeric(::Type{<:AbstractFloat}) = true
is_hard_numeric(::Type{<:Integer}) = true
is_hard_numeric(::Type{<:Rational}) = true
is_hard_numeric(::Type{ForwardDiff.Dual{T, V, N}}) where {T, V, N} = is_hard_numeric(V)
is_hard_numeric(::Type) = false
is_hard_numeric(x) = is_hard_numeric(typeof(x))

isdiagonal(a::AbstractMatrix{T}) where {T <: SymType} = isdiag(a)
isidentity(a::AbstractMatrix{T}) where {T <: SymType} = isone(a)

@inline LinearAlgebra.issymmetric(t::Tensor{2, 2, T}) where {T <: ApproxType} = @inbounds t[1, 2] ≈ t[2, 1]

@inline function LinearAlgebra.issymmetric(t::Tensor{2, 3, T}) where {T <: ApproxType}
    return @inbounds t[1, 2] ≈ t[2, 1] && t[1, 3] ≈ t[3, 1] && t[2, 3] ≈ t[3, 2]
end


"""
    _num_equal(a::Num, b::Num) -> Bool

Whether two `Symbolics.Num` are the same value, for the symmetry predicates.

**`iszero(a - b)` alone is not enough**, and assuming it was made every
symmetry predicate answer `false` on a perfectly symmetric symbolic tensor:
Symbolics folds `x - x` to zero only when `x` is an *atom*. For anything built
up, `(3k - 2μ)/3 - (3k - 2μ)/3` stays an unsimplified sum and `iszero` says
`false`. The consequence was not local — `_store_symmetric` then kept the full
9×9 form for a `Num` fourth-order tensor where a `Float64` one gets 6×6, which
is the same breakage the `ApproxType` union was introduced to cure for
`ForwardDiff.Dual`.

`isequal` is SymbolicUtils' *structural* equality: exact, cheap, and it settles
the case that matters here, two components built from the same expression.
`iszero` is kept as the second chance, for expressions that differ
syntactically but cancel. Neither calls `simplify`: this runs on every tensor
construction, and a predicate that is conservative (a false `false` costs
storage, never correctness) must not cost an algebraic simplification.
"""
@inline _num_equal(a, b) = isequal(a, b) || iszero(a - b)

@inline LinearAlgebra.issymmetric(t::Tensor{2, 2, Num}) = @inbounds _num_equal(t[1, 2], t[2, 1])

@inline function LinearAlgebra.issymmetric(t::Tensor{2, 3, Num})
    return @inbounds _num_equal(t[1, 2], t[2, 1]) && _num_equal(t[1, 3], t[3, 1]) &&
        _num_equal(t[2, 3], t[3, 2])
end

function Tensors.isminorsymmetric(t::Tensor{4, dim, T}) where {dim, T <: ApproxType}
    @inbounds for l in 1:dim, k in l:dim, j in 1:dim, i in j:dim
        if !(t[i, j, k, l] ≈ t[j, i, k, l]) || !(t[i, j, k, l] ≈ t[i, j, l, k])
            return false
        end
    end
    return true
end

function Tensors.isminorsymmetric(t::Tensor{4, dim, Num}) where {dim}
    @inbounds for l in 1:dim, k in l:dim, j in 1:dim, i in j:dim
        if !_num_equal(t[i, j, k, l], t[j, i, k, l]) ||
                !_num_equal(t[i, j, k, l], t[i, j, l, k])
            return false
        end
    end
    return true
end

function Tensors.ismajorsymmetric(t::FourthOrderTensor{dim, T}) where {dim, T <: ApproxType}
    @inbounds for l in 1:dim, k in l:dim, j in 1:dim, i in j:dim
        if !(t[i, j, k, l] ≈ t[k, l, i, j])
            return false
        end
    end
    return true
end

function Tensors.ismajorsymmetric(t::FourthOrderTensor{dim, Num}) where {dim}
    @inbounds for l in 1:dim, k in l:dim, j in 1:dim, i in j:dim
        # This one used to read `a - b == zero(Num)`, which on a `Num` builds a
        # symbolic *equation* rather than answering, so `!(…)` threw
        # `TypeError: non-boolean (Num) used in boolean context` — a hard
        # failure where its siblings merely returned a wrong `false`.
        if !_num_equal(t[i, j, k, l], t[k, l, i, j])
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

"""
    otimes(t1::AbstractArray, t2::AbstractArray)

Tensor (outer) product: `out[i…, j…] = t1[i…] * t2[j…]`.

No index is summed over, so there is nothing for a contraction engine to do
here. Going through `einsum` anyway cost 3.5 µs and 4.1 kB for two 3×3 arrays,
against 65 ns for the outer product written directly — the machinery, not the
arithmetic.

It is the identity case of `_interleaved_otimes`: `t1`'s indices go to
output positions `1…order1` and `t2`'s to `order1+1…order1+order2`, so each
operand is reshaped with singleton axes where the other's indices sit and the
product is one broadcast.

An earlier version wrote this as `vec(t1) .* transpose(vec(t2))`, which is the
same layout but **fails on `Tensors.Vec`**: `transpose` of a first-order
`Tensors` array is deliberately discontinued upstream, so a first-order operand
raised instead of multiplying. Singleton axes transpose nothing and have no such
restriction. Broadcasting also keeps this generic over `Dual` and symbolic
element types, where BLAS could not be used.
"""
function Tensors.otimes(
        t1::AbstractArray{T1, order1},
        t2::AbstractArray{T2, order2},
    ) where {T1, T2, order1, order2}
    return _interleaved_otimes(
        t1, t2, ntuple(i -> i, Val(order1)),
        ntuple(i -> order1 + i, Val(order2)), Val(order1 + order2),
    )
end

"""
    contract(t::AbstractArray, i::Integer, j::Integer)

Contract (trace over) indices `i` and `j` of a single tensor, lowering its order
by two. For an order-2 array this is the ordinary trace.

# Examples
```jldoctest
julia> contract([1.0 2.0; 3.0 4.0], 1, 2)
5.0
```
"""
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

"""
    qcontract(t1::AbstractArray, t2::AbstractArray)
    t1 ⊙ t2

Quadruple contraction: the last four indices of `t1` against the first four of
`t2`, lowering the total order by eight.

Between two order-4 tensors it is the **Frobenius scalar product**
`T ⊙ T' = Tᵢⱼₖₗ T'ᵢⱼₖₗ` — the inner product every projection in
[`proj_tens`](@ref) minimizes against.

# Examples
```jldoctest
julia> 𝕀, 𝕁, 𝕂 = ISO(Val(3), Val(Float64));

julia> (𝕁 ⊙ 𝕁, 𝕂 ⊙ 𝕂, 𝕁 ⊙ 𝕂)
(1.0, 5.0, 0.0)
```

`𝕁 ⊙ 𝕁 = 1` and `𝕂 ⊙ 𝕂 = 5` are the dimensions of the spherical and
deviatoric subspaces in 3-D.

See also `dcontract`, [`contract`](@ref).
"""
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

"""
    _interleaved_otimes(t1, t2, ec1, ec2, Val(n))

Outer product whose operand indices land at prescribed output positions:
`out[…] = t1[…] * t2[…]`, with `t1`'s `m`-th index at output position `ec1[m]`
and `t2`'s at `ec2[m]`.

`otimesu`, `otimesl` and the second term of `sotimes` are all of this shape —
nothing is summed, the two index lists merely interleave. Both lists are
**increasing**, so neither operand needs its axes permuted: giving each one
singleton axes where the other's indices go makes the whole thing a single
broadcast. One allocation, no contraction engine, no permutation pass, and
generic over `Dual` and symbolic element types.
"""
@inline function _otimes_shapes(
        t1::AbstractArray, t2::AbstractArray,
        ec1::NTuple{o1, Int}, ec2::NTuple{o2, Int}, ::Val{n},
    ) where {o1, o2, n}
    s1 = ntuple(
        k -> (m = findfirst(isequal(k), ec1); m === nothing ? 1 : size(t1, m)), Val(n)
    )
    s2 = ntuple(
        k -> (m = findfirst(isequal(k), ec2); m === nothing ? 1 : size(t2, m)), Val(n)
    )
    return s1, s2
end

@inline function _interleaved_otimes(
        t1::AbstractArray, t2::AbstractArray,
        ec1::NTuple{o1, Int}, ec2::NTuple{o2, Int}, v::Val{n},
    ) where {o1, o2, n}
    s1, s2 = _otimes_shapes(t1, t2, ec1, ec2, v)
    return reshape(t1, s1) .* reshape(t2, s2)
end

"""
    otimesu(t1::AbstractArray, t2::AbstractArray)
    t1 ⊠ t2

Modified tensor product, `(a ⊠ b)ᵢⱼₖₗ = aᵢₖ bⱼₗ`.

With the order-2 identity, `𝟙 = 𝟏 ⊠ 𝟏` is the identity of **all** order-2
tensors, where [`otimesul`](@ref) gives the identity of the *symmetric* ones.

Unlike its symmetrized counterpart it inverts termwise:
`(a ⊠ b)⁻¹ = a⁻¹ ⊠ b⁻¹`.

See [Tensor algebra](@ref th-tensor-algebra) for the full set of identities.
"""

function Tensors.otimesu(
        t1::AbstractArray{T1, order1},
        t2::AbstractArray{T2, order2},
    ) where {T1, T2, order1, order2}
    ec1 = (ntuple(i -> i, Val(order1 - 1))..., order1 + 1)
    ec2 = (order1, ntuple(i -> order1 + 1 + i, Val(order2 - 1))...)
    return _interleaved_otimes(t1, t2, ec1, ec2, Val(order1 + order2))
end

function Tensors.otimesl(
        t1::AbstractArray{T1, order1},
        t2::AbstractArray{T2, order2},
    ) where {T1, T2, order1, order2}
    ec1 = (ntuple(i -> i, Val(order1 - 1))..., order1 + 2)
    ec2 = (order1, order1 + 1, ntuple(i -> order1 + 2 + i, Val(order2 - 2))...)
    return _interleaved_otimes(t1, t2, ec1, ec2, Val(order1 + order2))
end

"""
    otimesul(t1::AbstractArray, t2::AbstractArray)
    sboxtimes(t1, t2)
    t1 ⊠ˢ t2

Symmetrized modified tensor product,
`(a ⊠ˢ b)ᵢⱼₖₗ = (aᵢₖ bⱼₗ + aᵢₗ bⱼₖ) / 2`, i.e. the half-sum of
[`otimesu`](@ref) and `otimesl`.

`𝕀 = 𝟏 ⊠ˢ 𝟏` is the identity of the **symmetric** order-2 tensors, which is
what [`tens_Id4`](@ref) returns.

!!! warning "It does not invert termwise"
    `(a ⊠ˢ b)⁻¹ ≠ a⁻¹ ⊠ˢ b⁻¹` unless `a` and `b` are **proportional**.
    Commuting is not sufficient, and taking `b = 𝟏` does not help. This is why
    inversion is implemented per symmetry class rather than by a single generic
    formula.
"""
otimesul(t1::AbstractArray{T1}, t2::AbstractArray{T2}) where {T1, T2} =
    (otimesu(t1, t2) + otimesl(t1, t2)) / promote_type(T1, T2)(2)

"""
    sotimes(t1::AbstractArray, t2::AbstractArray)
    t1 ⊗ˢ t2

Tensor product symmetrized over the **last** index of `t1` and the **first** of
`t2`. For two vectors, `u ⊗ˢ v = (u ⊗ v + v ⊗ u) / 2`.

It is the product used by [`SYMGRAD`](@ref), which is why that operator returns
a linearized strain tensor directly.

See also `otimes`, [`otimesul`](@ref).
"""
function sotimes(
        t1::AbstractArray{T1, order1},
        t2::AbstractArray{T2, order2},
    ) where {T1, T2, order1, order2}
    # The symmetrized product is the plain one averaged with its `otimesu`
    # partner. Both are outer products with singleton axes, so the average is
    # written as **one** fused broadcast rather than as two arrays plus their
    # sum: three allocations become one.
    n = Val(order1 + order2)
    a1, a2 = _otimes_shapes(
        t1, t2, ntuple(i -> i, Val(order1)),
        ntuple(i -> order1 + i, Val(order2)), n,
    )
    b1, b2 = _otimes_shapes(
        t1, t2, (ntuple(i -> i, Val(order1 - 1))..., order1 + 1),
        (order1, ntuple(i -> order1 + 1 + i, Val(order2 - 1))...), n,
    )
    # Kept as a division by 2, not a multiplication by 0.5: the two are
    # bit-identical in binary floating point but render differently for
    # symbolic element types, and this function is used symbolically.
    two = promote_type(T1, T2)(2)
    return (reshape(t1, a1) .* reshape(t2, a2) .+ reshape(t1, b1) .* reshape(t2, b2)) ./
        two
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
export tsimplify, tfactor, tsubs, tdiff, ttrigsimp, texpand_trig, tlimit
export is_hard_numeric
export ⋅, ⊡, ⊗
