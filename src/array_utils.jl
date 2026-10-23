δkron(T::Type{<:Number}, i::Integer, j::Integer) = i == j ? one(T) : zero(T)

struct Id2{dim, T <: Number} <: AbstractMatrix{T} end
Base.size(::Id2{dim}) where {dim} = (dim, dim)
Base.getindex(::Id2{dim, T}, i::Integer, j::Integer) where {dim, T} = δkron(T, i, j)
function Base.replace_in_print_matrix(::Id2, i::Integer, j::Integer, s::AbstractString)
    return i == j ? s : Base.replace_with_centered_mark(s)
end

"""
    isidentity(a::AbstractMatrix) → Bool

Test whether `a` is the identity matrix, within the tolerance of `≈` for
floating-point element types and exactly for symbolic ones. A matrix of
`ForwardDiff.Dual` numbers is the identity only if its derivatives vanish too.

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
isidentity(a::AbstractMatrix) = a ≈ I && !any(carries_partials, a)

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

"""
    ApproxType

Element types whose symmetry predicates use a *tolerance* rather than exact
equality: floating-point numbers, and `ForwardDiff.Dual` numbers built on them.

`ForwardDiff.Dual <: Real` but **not** `<: AbstractFloat`, so before this union
existed a `Dual`-valued tensor fell through to the exact fallback, where a
round-off of a few ulp was enough to report it as non-minor-symmetric. The
consequence was not local: `_KM_of_array` then built a 9×9 matrix instead of a
6×6 one, and every `proj_tens` call on a `Dual` input died with a
`DimensionMismatch` — automatic differentiation through a projection was
impossible. Symbolic types keep their own exact methods and are unaffected.

Companion of [`is_hard_numeric`](@ref), which answers a different question:
`ApproxType` says *use a tolerance rather than exact equality*,
`is_hard_numeric` says *a comparison is allowed at all*.
"""
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

# Exact element types (integers, rationals) compare exactly: `eps` is not
# defined on them. The tolerance is that of the floating-point type underneath.
isdiagonal(a::AbstractMatrix) = isdiag(a)
function isdiagonal(a::AbstractMatrix{T}) where {T <: ApproxType}
    off = a - Diagonal(a)
    return norm(off) <= _approx_eps(T) && !any(carries_partials, off)
end

"""
    carries_partials(x) -> Bool

Whether `x` carries a nonzero derivative, at any nesting level of
`ForwardDiff.Dual`. A plain number, or a symbolic one, carries none.

The structural predicates [`isidentity`](@ref) and `isdiagonal` decide whether
a basis collapses to a simpler one, and they compare **values** — `≈` and
`norm` see only the value of a `Dual`. A rotation by a dual angle at zero is the
identity in value, yet its derivative is the generator of the rotation:
collapsing it to the canonical basis returned zero derivatives. The predicates
therefore also require the derivatives to vanish.
"""
carries_partials(::Number) = false
carries_partials(x::ForwardDiff.Dual) =
    any(p -> !iszero(p) || carries_partials(p), ForwardDiff.partials(x)) ||
    carries_partials(ForwardDiff.value(x))

"""
    primal(x)

`x` stripped of every level of `ForwardDiff.Dual`: the value a choice that must not
depend on derivatives is made on (a quadrature order, a scaling, a bracket, a
pivot), so that the derivatives stay where the mathematics puts them. Any other
number, symbolic ones included, is returned as it is.

See also [`carries_partials`](@ref).
"""
primal(x::ForwardDiff.Dual) = primal(ForwardDiff.value(x))
primal(x) = x

"""
    primal_type(T) -> Type

The floating-point type under every level of `ForwardDiff.Dual` in `T`: the type
of `float(primal(x))` for `x::T`, to build a tolerance, a `typemax` or quadrature
weights that carry no derivative. An integer or rational type gives its float
type; any other (a complex or symbolic type) is returned as it is.
"""
primal_type(::Type{D}) where {D <: ForwardDiff.Dual} = primal_type(ForwardDiff.valtype(D))
primal_type(::Type{T}) where {T <: Real} = is_hard_numeric(T) ? float(T) : T
primal_type(::Type{T}) where {T} = T

"""
    dual_depth(T) -> Int

The number of levels of `ForwardDiff.Dual` nested in `T`: 0 for a plain number, 2
for a second derivative taken through nested `ForwardDiff` calls.
"""
dual_depth(::Type{<:ForwardDiff.Dual{<:Any, V}}) where {V} = 1 + dual_depth(V)
dual_depth(::Type) = 0

"""
    dual_channels(x) -> Tuple

The value of `x` followed by every partial, at every nesting level of
`ForwardDiff.Dual`, flattened into one tuple: `(x,)` for a plain number. A
convergence test that must hold for the derivatives as well as for the value
measures each channel separately.
"""
dual_channels(x::Number) = (x,)
dual_channels(d::ForwardDiff.Dual) = (
    dual_channels(ForwardDiff.value(d))...,
    foldl((a, b) -> (a..., b...), map(dual_channels, ForwardDiff.partials(d).values); init = ())...,
)

"""
    float_like(T) -> Type

`float(T)` for an integer type, `T` itself otherwise: a float, a dual or a
symbolic type stays as it is (`float` of a SymPy type answers `Float64`, which a
symbolic quantity cannot be converted to).
"""
float_like(::Type{T}) where {T <: Integer} = float(T)
float_like(::Type{T}) where {T} = T

isdiagonal(a::AbstractMatrix{T}) where {T <: SymType} = isdiag(a)
isidentity(a::AbstractMatrix{T}) where {T <: SymType} = isone(a)

# ── Index symmetries of a Tensors.jl array, as TensND needs them ─────────────
#
# For the storage decision (`_store_symmetric`) and nowhere else. They are
# TensND's own functions: Tensors.jl's `issymmetric`, `isminorsymmetric` and
# `ismajorsymmetric` keep their meaning on Tensors.jl's types (an exact
# comparison) for everyone, TensND loaded or not.

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

# Exact for exact element types (floats too: tolerance is `_store_symmetric`'s
# business), structural for a `Num`, whose `==` is a symbolic equation.
@inline _sym_eq(a, b) = a == b
@inline _sym_eq(a::Num, b::Num) = _num_equal(a, b)

_issym(::Tensors.SymmetricTensor) = true
_issym(t::Tensors.Tensor{2, dim}) where {dim} =
    all(_sym_eq(t[i, j], t[j, i]) for j in 1:dim for i in (j + 1):dim)
_issym(t::Tensors.Tensor{4}) = _isminorsym(t)

_isminorsym(::Tensors.SymmetricTensor{4}) = true
function _isminorsym(t::Tensors.Tensor{4, dim}) where {dim}
    return all(
        _sym_eq(t[i, j, k, l], t[j, i, k, l]) && _sym_eq(t[i, j, k, l], t[i, j, l, k])
            for l in 1:dim for k in l:dim for j in 1:dim for i in j:dim
    )
end
_ismajorsym(t::Tensors.FourthOrderTensor{dim}) where {dim} =
    all(_sym_eq(t[i, j, k, l], t[k, l, i, j]) for l in 1:dim for k in l:dim for j in 1:dim for i in j:dim)

# The major transpose of a `SymmetricTensor{4}` (Tensors.jl defines it on
# `Tensor{4}` only).
_majortranspose(S::Tensors.SymmetricTensor{4, dim}) where {dim} =
    Tensors.SymmetricTensor{4, dim}((i, j, k, l) -> @inbounds S[k, l, i, j])
_majortranspose(t) = Tensors.majortranspose(t)

"""
    otimes(t1::AbstractArray, t2::AbstractArray)

Tensor (outer) product: `out[i…, j…] = t1[i…] * t2[j…]`.

No index is summed over, so there is nothing to contract: the product is
written directly.

It is the identity case of `_interleaved_otimes`: `t1`'s indices go to
output positions `1…order1` and `t2`'s to `order1+1…order1+order2`, so each
operand is reshaped with singleton axes where the other's indices sit and the
product is one broadcast.

It is not written `vec(t1) .* transpose(vec(t2))`, the same layout: `transpose`
of a first-order `Tensors` array is deliberately discontinued upstream, so a
`Tensors.Vec` operand would raise instead of multiplying. Singleton axes
transpose nothing and have no such restriction. Broadcasting also keeps this generic over `Dual` and symbolic
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
    m, M = minmax(i, j)
    n = size(t, m)
    (m != M && size(t, M) == n) ||
        throw(DimensionMismatch("cannot trace index $i against index $j of an array of size $(size(t))"))
    # Position, in `t`, of each index of the result: all but `m` and `M`.
    kept = ntuple(q -> q < m ? q : (q < M - 1 ? q + 1 : q + 2), Val(order - 2))
    out = Array{T}(undef, map(q -> size(t, q), kept))
    @inbounds for I in CartesianIndices(out)
        s = zero(T)
        for x in 1:n
            idx = ntuple(q -> q == m || q == M ? x : I[q < m ? q : (q < M ? q - 1 : q - 2)], Val(order))
            s += t[idx...]
        end
        out[I] = s
    end
    return out
end

contract(t::AbstractArray{T, 2}, ::Integer, ::Integer) where {T} = tr(t)

Tensors.dcontract(t1::AbstractArray, t2::AbstractArray) = _contract_trailing(t1, t2, Val(2))

# `Σ aᵢ bᵢ` over all the components, as every contraction is. Not `dot`, which
# conjugates its first argument (the Hermitian product): it gave `Σ āᵢⱼ bᵢⱼ` on
# complex components and `conjugate(a)` on a `SymPy.Sym` not declared real.
# `vec` also keeps a TensND tensor, itself an `AbstractArray`, away from the
# tensor `dot`, which is a single contraction.
_full_contraction(t1::AbstractArray, t2::AbstractArray) = transpose(vec(t1)) * vec(t2)

Tensors.dcontract(t1::AbstractArray{T1, 2}, t2::AbstractArray{T2, 2}) where {T1, T2} =
    _full_contraction(t1, t2)

Tensors.dotdot(v1::AbstractArray, S::AbstractArray, v2::AbstractArray) =
    _contract_trailing(_contract_trailing(v1, S, Val(1)), v2, Val(1))

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
qcontract(t1::AbstractArray, t2::AbstractArray) = _contract_trailing(t1, t2, Val(4))

qcontract(t1::AbstractArray{T1, 4}, t2::AbstractArray{T2, 4}) where {T1, T2} =
    _full_contraction(t1, t2)

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

# A full `Tensor` even for two symmetric operands: symmetrizing the middle pair
# of indices breaks the minor symmetry `i ↔ j`. Same products, same order as
# the array method above, hence the same bits.
@inline function sotimes(S1::SecondOrderTensor{dim}, S2::SecondOrderTensor{dim}) where {dim}
    return Tensor{4, dim}(
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
"""
    isparallel(a, b; rtol = 1.0e-10) -> Bool

Whether the directions `a` and `b` (3-vectors, tuples or 1st-order tensors, read in
canonical components) are parallel, in the same or the opposite sense:
``(\\underline{a}\\cdot\\underline{b})^2 = |\\underline{a}|^2\\,|\\underline{b}|^2``.

For comparable numbers (`is_hard_numeric`), within the relative tolerance `rtol`;
a `ForwardDiff.Dual` is compared on its value. For symbolic ones the test is
structural: the difference is simplified and must vanish exactly, so a symbolic
answer is never a false positive, and is `false` when the backend cannot decide.
No `dot` and no `abs`: both conjugate a `SymPy.Sym`, and `Symbolics` does not
fold `abs` of a literal.
"""
function isparallel(a, b; rtol::Real = 1.0e-10)
    va, vb = _extract_vec(a), _extract_vec(b)
    d = sum(va .* vb)
    nn = sum(va .* va) * sum(vb .* vb)
    return _isparallel(Val(is_hard_numeric(promote_type(typeof(d), typeof(nn)))), d^2, nn, rtol)
end
_isparallel(::Val{true}, dd, nn, rtol) = !iszero(nn) && abs(dd - nn) <= rtol * nn
function _isparallel(::Val{false}, dd, nn, rtol)
    verdict = iszero(_simplified(dd - nn))
    return verdict isa Bool ? verdict : false
end
_simplified(x::SymType) = tsimplify(x)
_simplified(x) = x

export is_hard_numeric, carries_partials, primal, isparallel, primal_type, dual_depth, dual_channels, float_like
# `dotdot` belongs here too: its docstring offers `dotdot(n, ℂ, n)` as the
# acoustic tensor, which a `using TensND` alone could not reach.
export ⋅, ⊡, ⊗, dotdot
