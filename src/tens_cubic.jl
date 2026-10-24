# ─────────────────────────────────────────────────────────────────────────────
#  tens_cubic.jl — 4th-order tensors of CUBIC symmetry.
#
#  A hybrid of the two types that surround it, and it is worth saying which half
#  comes from where.
#
#  **The algebra is `TensISO`'s.** Under the octahedral group the space of
#  minor-symmetric 4th-order tensors splits as A1g + Eg + T2g, of dimensions
#  1 + 2 + 3 — three inequivalent irreducible representations, each of
#  multiplicity one. The commutant is therefore spanned by three **mutually
#  orthogonal projectors** summing to the identity, exactly as 𝕁 and 𝕂 span the
#  isotropic one. So sums, products and inverses are all componentwise:
#
#      ℂ = α𝕁 + β𝔼 + γ𝕋 ,   ℂ⁻¹ = α⁻¹𝕁 + β⁻¹𝔼 + γ⁻¹𝕋 ,
#      ℂ ⊡ ℂ′ = αα′𝕁 + ββ′𝔼 + γγ′𝕋 .
#
#  Two consequences follow that neither `TensTI` nor `TensOrtho` enjoys. The
#  double contraction of two cubic tensors is **cubic**, where TI widens from
#  N=5 to N=6 and orthotropic widens out of its class entirely. And every cubic
#  tensor is automatically **major-symmetric**, the three projectors being
#  symmetric: a localization tensor, which in general has no major symmetry,
#  recovers it here.
#
#  **The geometry is `TensOrtho`'s.** Unlike an isotropic tensor, a cubic one is
#  only cubic relative to its cube axes, so the type carries an orthonormal
#  frame. Where it differs from orthotropy is that the three projectors are
#  invariant under the *whole* octahedral group, so two frames related by a
#  signed permutation of the axes describe the same tensor with the same
#  coefficients — see `_same_cube_frame`. Permuting the axes of an orthotropic
#  tensor permutes C₁₁, C₂₂, C₃₃; permuting those of a cubic one changes
#  nothing.
#
#  Storage is the projector coefficients `(α, β, γ)`, following `TensISO`'s
#  `(3k, 2μ)` rather than the engineering constants, because that is what makes
#  the algebra above literal. `tens_cubic` / `arg_cubic` convert.
# ─────────────────────────────────────────────────────────────────────────────

"""
    TensCubic{T, B<:OrthonormalBasis{3}} <: AbstractTens{4,3,T}

Cubic 4th-order tensor with cube axes `frame` and **three** independent
coefficients `(α, β, γ)` on the three orthogonal projectors of the octahedral
group:

```math
\\mathbb C = \\alpha\\,\\mathbb J + \\beta\\,\\mathbb E + \\gamma\\,\\mathbb T ,
```

where, in Kelvin-Mandel form and in the cube frame,

```math
\\mathbb J = \\tfrac13 \\underline v\\,\\underline v^{\\mathsf T}, \\quad
\\underline v = (1,1,1,0,0,0)^{\\mathsf T}, \\qquad
\\mathbb E = \\begin{pmatrix} \\mathbb I_3 - \\tfrac13\\underline 1\\,\\underline 1^{\\mathsf T} & 0 \\\\ 0 & 0\\end{pmatrix},
\\qquad
\\mathbb T = \\begin{pmatrix} 0 & 0 \\\\ 0 & \\mathbb I_3\\end{pmatrix} .
```

They are mutually orthogonal, idempotent, and sum to the identity, with traces
`1, 2, 3`. The Kelvin-Mandel matrix in the cube frame is therefore

    [[C₁₁, C₁₂, C₁₂,  0,    0,    0  ],
     [C₁₂, C₁₁, C₁₂,  0,    0,    0  ],
     [C₁₂, C₁₂, C₁₁,  0,    0,    0  ],
     [ 0,   0,   0,  2C₄₄,  0,    0  ],
     [ 0,   0,   0,   0,   2C₄₄,  0  ],
     [ 0,   0,   0,   0,    0,   2C₄₄]]

with `α = C₁₁ + 2C₁₂ = 3K`, `β = C₁₁ - C₁₂`, `γ = 2C₄₄`. Use
[`tens_cubic`](@ref) to build one from `(C₁₁, C₁₂, C₄₄)` and [`arg_cubic`](@ref)
to read them back.

The isotropic class is the case `β = γ`: `TensISO{3}(α, β)` is
`α𝕁 + β(𝔼 + 𝕋)`, since `𝕂 = 𝔼 + 𝕋`. What a two-constant approximation of a
cubic tensor discards is measured by [`cubic_anisotropy`](@ref).

| Parametrization    | Role                | Invariance group | Stored coefficients                  |
|--------------------|---------------------|------------------|--------------------------------------|
| `TensCubic{T, B}`  | 4th-order cubic     | `Oₕ`             | `data = (α, β, γ)`, `frame::B`       |

One shape only, so the type carries neither order nor length: at order 2 a
cubic tensor is isotropic (a `TensISO{2}`), and the product of two cubic
tensors on a cube is cubic and major-symmetric, so there is no general shape
(`Oₕ` and its rotation group `O` have the same invariants at even order). The
inverse stays cubic, and a cubic tensor is its own transpose.

As for [`TensOrtho`](@ref), the frame's element type `B` is a free type
parameter independent of the data eltype `T`: differentiating with respect to
the constants does not require a `Dual`-typed geometric frame.

# Examples
```julia
julia> C = tens_cubic(10.0, 4.0, 2.0, CanonicalBasis{3,Float64}());

julia> arg_cubic(C)
(10.0, 4.0, 2.0)

julia> arg_cubic(inv(C) ⊡ C)          # the identity, exactly
(1.0, 0.0, 0.5)
```

See also [`tens_cubic`](@ref), [`arg_cubic`](@ref),
[`cubic_anisotropy`](@ref), [`iso_to_cubic`](@ref), [`cubic_to_ortho`](@ref).
"""
struct TensCubic{T, B <: OrthonormalBasis{3}} <: AbstractTens{4, 3, T}
    data::NTuple{3, T}     # (α, β, γ) on (𝕁, 𝔼, 𝕋)
    frame::B               # cube axes (e₁, e₂, e₃)
end

# Outer constructor preserving the `TensCubic{T}(data, frame)` call form used by
# `_rebuild` and the constructors below — infers `B` from the frame.
TensCubic{T}(data::NTuple{3}, frame::B) where {T, B <: OrthonormalBasis{3}} =
    TensCubic{T, B}(data, frame)

# ── Traits ───────────────────────────────────────────────────────────────────

Base.eltype(::Type{TensCubic{T}}) where {T} = T
Base.length(::TensCubic) = 81
Base.size(::TensCubic) = (3, 3, 3, 3)

get_basis(::TensCubic{T}) where {T} = CanonicalBasis{3, T}()
get_var(::TensCubic) = (:cont, :cont, :cont, :cont)
get_var(::TensCubic, ::Integer) = :cont
get_data(t::TensCubic) = t.data
frame(t::TensCubic) = t.frame

# ── Rebuild helper (used by the symbolic ops) ────────────────────────────────
_rebuild(t::TensCubic, new_data) = TensCubic{eltype(new_data)}(new_data, frame(t))

# ── Constructors ─────────────────────────────────────────────────────────────

"""
    TensCubic(α, β, γ, frame)

Cubic tensor from its three **projector** coefficients on `(𝕁, 𝔼, 𝕋)`. For the
engineering constants use [`tens_cubic`](@ref).
"""
function TensCubic(α, β, γ, frame::OrthonormalBasis{3})
    T = promote_type(typeof(α), typeof(β), typeof(γ), eltype(frame))
    return TensCubic{T}((T(α), T(β), T(γ)), frame)
end

"""
    tens_cubic(C₁₁, C₁₂, C₄₄, frame) -> TensCubic

Cubic tensor from the three elastic constants, `C₄₄ = C₂₃₂₃`:

```math
\\alpha = C_{11} + 2C_{12}, \\qquad \\beta = C_{11} - C_{12},
\\qquad \\gamma = 2C_{44}.
```

`C₁₁ = C₁₂ + 2C₄₄` gives back an isotropic tensor — the identity worth
remembering, since the departure from it is what [`cubic_anisotropy`](@ref)
measures.
"""
tens_cubic(C₁₁, C₁₂, C₄₄, frame::OrthonormalBasis{3}) =
    TensCubic(C₁₁ + 2C₁₂, C₁₁ - C₁₂, 2C₄₄, frame)

"""
    TensCubic(KMmat::AbstractMatrix, frame)

Build a `TensCubic` from a 6×6 Kelvin-Mandel matrix **expressed in the cube
frame**. Only `KM[1,1]`, `KM[1,2]` and `KM[4,4]` are read: the matrix is assumed
to have the cubic structure, and anything else in it is discarded silently. To
project a general tensor onto the class instead, and to learn how far it was,
use [`best_fit_cubic`](@ref) or `proj_tens(Val(:CUBIC), …)`.
"""
TensCubic(KMmat::AbstractMatrix, frame::OrthonormalBasis{3}) =
    tens_cubic(KMmat[1, 1], KMmat[1, 2], KMmat[4, 4] / 2, frame)

"""
    arg_cubic(t::TensCubic) -> (C₁₁, C₁₂, C₄₄)

The three elastic constants, the inverse of [`tens_cubic`](@ref).
"""
function arg_cubic(t::TensCubic)
    α, β, γ = get_data(t)
    return (α / 3 + 2β / 3, α / 3 - β / 3, γ / 2)
end

"""
    cubic_anisotropy(t::TensCubic) -> Real

``(C_{11} - C_{12} - 2C_{44})/C_{11} = (\\beta - \\gamma)/C_{11}``: a
Zener-type measure of how far a cubic tensor departs from an isotropic one.

Exactly zero when `β = γ`, which is exactly when the tensor is isotropic. It is
precisely the quantity a two-constant approximation of a cubic tensor throws
away.
"""
function cubic_anisotropy(t::TensCubic)
    α, β, γ = get_data(t)
    C₁₁ = α / 3 + 2β / 3
    return (β - γ) / C₁₁
end

# ── get_array and getindex ───────────────────────────────────────────────────

"""
    get_array(t::TensCubic{T}) -> Array{T,4}

The 3×3×3×3 component array in the canonical frame.

Evaluated through `_ortho_entry` with `C₁₁ = C₂₂ = C₃₃`, `C₁₂ = C₁₃ = C₂₃` and
`C₄₄ = C₅₅ = C₆₆`, which is what a cubic tensor *is* as an orthotropic one — no
second copy of that algebra.
"""
get_array(t::TensCubic) = Array(_static_tensor(t))

"""
    getindex(t::TensCubic, i, j, k, l)

Single component in closed form, so that the generic algorithms that walk a
tensor element by element — `Array`, `norm`, broadcasting, `show`, `tomandel` —
do not each pay for materializing all 81.
"""
Base.@propagate_inbounds function Base.getindex(
        t::TensCubic, i::Integer, j::Integer, k::Integer, l::Integer
    )
    C₁₁, C₁₂, C₄₄ = arg_cubic(t)
    return _ortho_entry(
        vecbasis(t.frame, :cov), C₁₁, C₁₁, C₁₁, C₁₂, C₁₂, C₁₂, C₄₄, C₄₄, C₄₄, i, j, k, l
    )
end

Base.IndexStyle(::Type{<:TensCubic}) = IndexCartesian()

# ── Kelvin-Mandel ────────────────────────────────────────────────────────────

"""
    mandel(t::TensCubic)

The 6×6 Kelvin-Mandel matrix in the **canonical** frame. Use
[`mandel_material`](@ref) for the sparse form in the cube frame.
"""
mandel(t::TensCubic) = Matrix(_km_in_frame(t))

"""
    mandel_material(t::TensCubic)

The 6×6 Kelvin-Mandel matrix in the cube frame, where a cubic tensor has only
three distinct entries.
"""
function mandel_material(t::TensCubic{T}) where {T}
    C₁₁, C₁₂, C₄₄ = arg_cubic(t)
    z = zero(T)
    return @SMatrix [
        C₁₁ C₁₂ C₁₂  z    z    z;
        C₁₂ C₁₁ C₁₂  z    z    z;
        C₁₂ C₁₂ C₁₁  z    z    z;
        z   z   z   2C₄₄  z    z;
        z   z   z    z   2C₄₄  z;
        z   z   z    z    z   2C₄₄
    ]
end

# ── Frames that describe the same cube ───────────────────────────────────────

"""
    _same_cube_frame(A::TensCubic, B::TensCubic) -> Bool

Whether two cube frames describe the **same** cube, which for this class is a
weaker and more useful question than whether they are equal.

The three projectors are invariant under the whole octahedral group, so two
frames related by a signed permutation of the axes carry the *same*
coefficients: `Fᴬᵀ Fᴮ` being a signed permutation matrix is exactly the
condition, and it costs nine comparisons rather than enumerating 24 rotations.
That is genuinely specific to cubic symmetry — permuting the axes of an
orthotropic tensor permutes `C₁₁, C₂₂, C₃₃`, so `TensOrtho` must and does
require frame equality.

On an element type where a comparison does not return an honest `Bool` — a
symbolic frame — the test is not attempted and plain equality is required
instead.
"""
function _same_cube_frame(A::TensCubic, B::TensCubic)
    FA, FB = vecbasis(frame(A), :cov), vecbasis(frame(B), :cov)
    is_hard_numeric(real(promote_type(eltype(FA), eltype(FB)))) ||
        return frame(A) == frame(B)
    M = transpose(FA) * FB
    for i in 1:3
        hits = 0
        for j in 1:3
            a = abs(M[i, j])
            if a > 0.5
                isapprox(a, one(a); atol = 1.0e-10) || return false
                hits += 1
            elseif a > 1.0e-10
                return false
            end
        end
        hits == 1 || return false
    end
    return true
end

@inline _check_same_reference(A::TensCubic, B::TensCubic) =
    _same_cube_frame(A, B) ||
    throw(ArgumentError("TensCubic operation requires frames describing the same cube"))

# ── Arithmetic ───────────────────────────────────────────────────────────────
# Scalar operations (-, α*A, A*α, A/α) come from structured_tens_ops.jl.

@inline function Base.:+(A::TensCubic, B::TensCubic)
    _check_same_reference(A, B)
    return _rebuild(A, get_data(A) .+ get_data(B))
end
@inline function Base.:-(A::TensCubic, B::TensCubic)
    _check_same_reference(A, B)
    return _rebuild(A, get_data(A) .- get_data(B))
end

# ── Inverse, identity, zero ──────────────────────────────────────────────────

"""
    inv(t::TensCubic) -> TensCubic

`α⁻¹𝕁 + β⁻¹𝔼 + γ⁻¹𝕋`. The three projectors being orthogonal and idempotent,
the inverse is componentwise — no adjugate, no factorization, and the result is
still cubic on the same frame.
"""
Base.inv(t::TensCubic) = _rebuild(t, map(inv, get_data(t)))

"""
    one(t::TensCubic)

The 4th-order symmetric identity `𝕀 = 𝕁 + 𝔼 + 𝕋`, cubic on the same frame.
"""
Base.one(t::TensCubic{T}) where {T} = _rebuild(t, (one(T), one(T), one(T)))

"""
    zero(t::TensCubic)

The zero tensor, kept **in the class and on the same frame** rather than
narrowed to a `TensISO`. Narrowing would silently change the type of an
accumulator initialized with `zero`, which has caught this package out before.
"""
Base.zero(t::TensCubic{T}) where {T} = _rebuild(t, (zero(T), zero(T), zero(T)))

@inline Base.literal_pow(::typeof(^), A::TensCubic, ::Val{-1}) = inv(A)
@inline Base.literal_pow(::typeof(^), A::TensCubic, ::Val{0}) = one(A)
@inline Base.literal_pow(::typeof(^), A::TensCubic, ::Val{1}) = A
@inline Base.literal_pow(::typeof(^), A::TensCubic, ::Val{n}) where {n} =
    _rebuild(A, map(x -> x^n, get_data(A)))

@inline Base.transpose(A::TensCubic) = A
@inline Base.adjoint(A::TensCubic) = A

# ── Double contraction ───────────────────────────────────────────────────────

"""
    dcontract(A::TensCubic, B::TensCubic) -> TensCubic

`(α_Aα_B, β_Aβ_B, γ_Aγ_B)` on the same frame.

**The product stays in the class**, which neither of the other structured types
manages: `dcontract(::TensTI{4}, ::TensTI{4})` widens from five Walpole
coefficients to six, and `dcontract(::TensOrtho, ::TensOrtho)` from nine
orthotropic constants to twelve, because the product of two symmetric 3×3
blocks is not symmetric unless they commute. Here the three projectors are one-dimensional in
the sense that matters — each irreducible representation appears once — so
every cubic tensor commutes with every other on the same cube, and the product
is cubic and still major-symmetric.

If the two frames do not describe the same cube the product is generally fully
anisotropic and the generic route is taken.
"""
function Tensors.dcontract(A::TensCubic, B::TensCubic)
    _same_cube_frame(A, B) ||
        return Tensors.dcontract(generic_tens(A), generic_tens(B))
    return _rebuild(A, get_data(A) .* get_data(B))
end

# ── Symmetry ─────────────────────────────────────────────────────────────────

LinearAlgebra.issymmetric(::TensCubic) = true
Tensors.isminorsymmetric(::TensCubic) = true
# Automatic, and not by construction: the three projectors are symmetric, so any
# combination of them is. It is why a cubic *localization* tensor is
# major-symmetric although a localization tensor in general is not.
Tensors.ismajorsymmetric(::TensCubic) = true

# ── Class predicates and the unified accessors ───────────────────────────────

isisotropic(::TensCubic) = false
istransverselyisotropic(::TensCubic) = false
isorthotropic(::TensCubic) = false

"""
    iscubic(t) -> Bool

Whether `t` is stored as a [`TensCubic`](@ref). A *type-level* query, like
[`isisotropic`](@ref) and [`istransverselyisotropic`](@ref): it says what the container guarantees,
not whether the components happen to satisfy a tighter symmetry.
"""
iscubic(::TensCubic) = true
iscubic(::Any) = false

symmetry(::TensCubic) = :CUBIC
reference(t::TensCubic) = frame(t)

# ── change_tens / components ─────────────────────────────────────────────────
#
# Same reasoning as for `TensTI{4}` and `TensOrtho`, and the same three methods.
# A `TensCubic` stores its coefficients with `get_basis` canonical and the cube
# frame kept separately, so re-expressing it in another orthonormal basis is a
# genuine multi-index rotation: wrapping the stored components in the new basis
# would relabel them and silently return a different physical tensor. The
# rotation is done in closed form, for `TensOrtho` as well: in an orthonormal
# basis ℬ the tensor has the same coefficients with its frame written in ℬ, and
# `_ortho_entry` gives its components there, exact for every element type,
# without the 81 canonical ones rotated index by index. A left-handed ℬ needs no
# care: each frame vector enters the formula an even number of times.
function change_tens(t::Union{TensOrtho{4, T}, TensCubic{T}}, ℬ::OrthonormalBasis{3, T}) where {T}
    F = vecbasis(frame(t), :cov)
    Tf = promote_type(T, eltype(F))
    E = transpose(SMatrix{3, 3, Tf}(vecbasis(ℬ, :cov))) * SMatrix{3, 3, Tf}(F)
    C = _ortho_coefficients(t)
    S = SymmetricTensor{4, 3}((i, j, k, l) -> @inbounds _ortho_entry(E, C..., i, j, k, l))
    return Tens(S, ℬ, get_var(t))
end

# The nine orthotropic coefficients, a cubic tensor being an orthotropic one
# with C₁₁ = C₂₂ = C₃₃, C₁₂ = C₁₃ = C₂₃ and C₄₄ = C₅₅ = C₆₆.
_ortho_coefficients(t::TensOrtho{4}) = get_data(t)
function _ortho_coefficients(t::TensCubic)
    C₁₁, C₁₂, C₄₄ = arg_cubic(t)
    return (C₁₁, C₁₁, C₁₁, C₁₂, C₁₂, C₁₂, C₄₄, C₄₄, C₄₄)
end
components(t::TensCubic{T}, ℬ::OrthonormalBasis{3, T}, v::NTuple{4, Symbol}) where {T} =
    components(_as_canonical_tens(t), ℬ, v)
components(t::TensCubic) = get_array(t)
components(t::TensCubic, ::NTuple{4, Symbol}) = get_array(t)

# ── Display ──────────────────────────────────────────────────────────────────

function Base.show(io::IO, A::TensCubic)
    α, β, γ = get_data(A)
    print(io, "(", α, ") 𝕁 + (", β, ") 𝔼 + (", γ, ") 𝕋")
    return print(io, "\n  cube frame: ", vecbasis(A.frame, :cov))
end

function pprint(A::TensCubic)
    α, β, γ = get_data(A)
    C₁₁, C₁₂, C₄₄ = arg_cubic(A)
    println("(", α, ") 𝕁 + (", β, ") 𝔼 + (", γ, ") 𝕋")
    println("  C₁₁ = ", C₁₁, ", C₁₂ = ", C₁₂, ", C₄₄ = ", C₄₄)
    return println("  cube frame: ", vecbasis(A.frame, :cov))
end

##############################################################################
# Exports
##############################################################################

export TensCubic, tens_cubic, arg_cubic, cubic_anisotropy

# Structured tensors: their decomposition (the compact `show`), then their
# Kelvin–Mandel matrix.
function Base.show(io::IO, ::MIME"text/plain", A::Union{TensISO, TensTI, TensOrtho, TensCubic})
    println(io, typeof(A))
    print(io, "→ decomposition: ")
    show(io, A)
    return _show_km(io, A)
end

_show_km(io::IO, A::AbstractTens{2}) = nothing
function _show_km(io::IO, A::AbstractTens{4})
    print(io, "\n→ KM: ")
    return show(io, MIME"text/plain"(), mandel(A))
end
function _show_km(io::IO, A::TensOrtho{4})
    print(io, "\n→ KM (material frame): ")
    show(io, MIME"text/plain"(), mandel_material(A))
    print(io, "\n→ KM (canonical frame): ")
    return show(io, MIME"text/plain"(), mandel(A))
end
