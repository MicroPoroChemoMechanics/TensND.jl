#########################################
# Random tensors: uniform orientations
#########################################

# Every draw below has a law that no orientation favors. A rotation is drawn
# from the Haar measure of SO(3), an axis uniformly on the unit sphere, and the
# coefficients of a tensor through its Kelvin–Mandel blocks, from laws that are
# invariant under the orthogonal changes of those blocks (see
# `theory/rotations.md` and `manual/random_tensors.md`).

"""
    random_rotation([rng], [T = Float64]) -> SMatrix{3, 3, T}

A rotation drawn from the Haar measure of SO(3): the uniform law on the
rotations, which no orientation favors. It is the rotation of a unit
quaternion drawn uniformly on the unit sphere of ℝ⁴, the normalization of four
independent standard normal numbers (see [Random rotations](@ref th-random-rotations)).

Drawing the three Euler angles uniformly does **not** give this law: the axis
`e₃` would then fall near the poles as often as near the equator.

See also [`random_axis`](@ref), [`random_basis`](@ref), [`rand`](@ref Base.rand(::Type{<:TensND._RandomClass})).
"""
random_rotation(rng::AbstractRNG, ::Type{T} = Float64) where {T <: AbstractFloat} =
    SMatrix{3, 3, T}(rand(rng, RotMatrix{3, T}))
random_rotation(::Type{T} = Float64) where {T <: AbstractFloat} = random_rotation(Random.default_rng(), T)

"""
    random_axis([rng], [T = Float64]) -> NTuple{3, T}

A unit vector drawn uniformly on the unit sphere: three independent standard
normal numbers, normalized. It is the axis `R⋅e₃` of a [`random_rotation`](@ref)
`R`, and the law of the axis of a random transversely isotropic tensor.
"""
random_axis(rng::AbstractRNG, ::Type{T} = Float64) where {T <: AbstractFloat} =
    Tuple(normalize(randn(rng, SVector{3, T})))
random_axis(::Type{T} = Float64) where {T <: AbstractFloat} = random_axis(Random.default_rng(), T)

"""
    random_basis([rng], [T = Float64]) -> RotatedBasis{3, T}

The orthonormal basis of the columns of a [`random_rotation`](@ref): a frame
of uniformly random orientation, the law of the frame of a random orthotropic
or cubic tensor.
"""
random_basis(rng::AbstractRNG, ::Type{T} = Float64) where {T <: AbstractFloat} =
    RotatedBasis(random_rotation(rng, T), Val(3), Val(:typed))
random_basis(::Type{T} = Float64) where {T <: AbstractFloat} = random_basis(Random.default_rng(), T)

# ── The laws of the coefficients ─────────────────────────────────────────────

# An eigenvalue: standard normal, or log-normal when the tensor is to be
# positive definite.
_random_eigenvalue(rng, ::Type{T}, definite::Bool) where {T} = definite ? exp(randn(rng, T)) : randn(rng, T)

# An orthogonal matrix from the Haar measure of O(n): the Q factor of a matrix
# of independent standard normal entries, the signs of its columns fixed by the
# diagonal of R. The law of the Gaussian matrix is invariant under a left
# orthogonal factor, and so is the Q of the decomposition made unique.
function _haar_orthogonal(rng, ::Type{T}, ::Val{n}) where {T, n}
    F = qr(randn(rng, SMatrix{n, n, T}))
    return F.Q * Diagonal(sign.(diag(F.R)))
end

# A symmetric n×n block: from the Gaussian orthogonal ensemble (standard normal
# diagonal, off-diagonal of variance 1/2), or, positive definite, with Haar
# eigenvectors and log-normal eigenvalues. Both laws are invariant under every
# orthogonal change of basis of the block.
function _random_symmetric(rng, ::Type{T}, ::Val{n}, definite::Bool) where {T, n}
    if definite
        Q = _haar_orthogonal(rng, T, Val(n))
        A = Q * Diagonal(exp.(randn(rng, SVector{n, T}))) * transpose(Q)
    else
        A = randn(rng, SMatrix{n, n, T})
    end
    return (A + transpose(A)) / 2
end

# A block without symmetry: independent standard normal entries (the Ginibre
# ensemble), invariant under orthogonal changes on either side.
_random_general(rng, ::Type{T}, ::Val{n}) where {T, n} = randn(rng, SMatrix{n, n, T})

_refuse_definite(definite::Bool, what) =
    definite && throw(
    ArgumentError(
        "definite = true asks for a positive definite tensor, which needs the major " *
            "symmetry; $what has none: drop `definite`, or draw the symmetric shape"
    )
)

# The classes `rand` draws. A class given without its parameters (`TensCubic`,
# `TensTI{4}`) is not a subtype of `AbstractTens`, whose element type is
# bounded, so the union of the classes is what the methods dispatch on.
const _RandomClass = Union{TensISO, TensTI, TensOrtho, TensCubic, Tens}

# ── rand ──────────────────────────────────────────────────────────────

"""
    rand([rng], class; definite = false, ...) -> AbstractTens
    rand([rng], class, dims...) -> Array

A tensor of `class` with random coefficients and, for a class with a
direction, a random orientation. `class` is the type of the tensor, given
either by its class with the default shape and `Float64` coefficients, or in
full for another shape or element type:

| `class` | returns | direction (keyword) |
|:--|:--|:--|
| `TensISO{4, dim}`, `TensISO{2, dim}` | isotropic | none |
| `TensTI{4}` (`N = 5`), `TensTI{4, T, N}` (`N = 5, 6, 8`) | transversely isotropic | `axis = random_axis(rng, T)`, or `frame` |
| `TensTI{2}` (`N = 2`), `TensTI{2, T, N}` (`N = 2, 3`) | transversely isotropic | `axis = random_axis(rng, T)`, or `frame` |
| `TensOrtho{4}` (`N = 9`), `TensOrtho{4, T, N}` (`N = 9, 12`) | orthotropic | `frame = random_basis(rng, T)` |
| `TensOrtho{2}`, `TensOrtho{2, T}` | orthotropic | `frame = random_basis(rng, T)` |
| `TensCubic`, `TensCubic{T}` | cubic | `frame = random_basis(rng, T)` |
| `Tens{order, dim}`, `Tens{order, dim, T}` (`order = 2, 4`) | generic | none: its law is invariant |

Passing `axis` or `frame` fixes the direction and leaves only the coefficients
random; a transversely isotropic tensor takes either, the axis being the third
vector of the frame. Without a generator `rng`, `Random.default_rng()` is used;
`rand(class, n)` draws an array of `n` tensors of the default law, as `rand`
does for numbers; with keywords, an array is a comprehension,
`[rand(rng, class; definite = true) for _ in 1:n]`. `definite = true` draws a positive definite tensor, for the shapes
with the major symmetry (the default ones); the others refuse it. A generic
tensor takes `symmetry`: `Val(:major)` (the default, minor and major
symmetries, the class of elastic moduli), `Val(:minor)` or `Val(:none)` at
order 4, `Val(:symmetric)` (the default) or `Val(:none)` at order 2.

**The laws.** The Kelvin–Mandel matrix of the tensor in its own frame has the
block structure of its class. Each symmetric block is drawn from the Gaussian
orthogonal ensemble, or, positive definite, with Haar eigenvectors and
log-normal eigenvalues; each block without symmetry has independent standard
normal entries; a scalar block is standard normal, or log-normal. Each law is
invariant under the orthogonal changes of basis of its block, so that a
generic tensor needs no random frame, its law being already invariant under
rotations. See [Random tensors](@ref man-random-tensors).

# Examples
```julia
julia> using Random; rng = Xoshiro(1);

julia> ℂ = rand(rng, TensCubic; definite = true);   # a cubic crystal, random cube

julia> ℂ isa TensCubic && minimum(eigvals(Symmetric(Matrix(mandel(ℂ))))) > 0
true

julia> rand(rng, TensTI{4}; axis = (0.0, 0.0, 1.0)) isa TensTI{4, Float64, 5}
true
```
"""
Random.rand(class::Type{<:_RandomClass}; kwargs...) = rand(Random.default_rng(), class; kwargs...)

# `rand(class, dims...)` and `rand!`: the draws of the default law, one per entry.
Random.rand(rng::AbstractRNG, ::Random.SamplerType{C}) where {C <: _RandomClass} = rand(rng, C)

# Isotropic.
Random.rand(rng::AbstractRNG, ::Type{TensISO{order, dim}}; kwargs...) where {order, dim} =
    rand(rng, TensISO{order, dim, Float64}; kwargs...)
Random.rand(rng::AbstractRNG, ::Type{TensISO{4, dim, T}}; definite::Bool = false) where {dim, T} =
    TensISO{dim}(_random_eigenvalue(rng, T, definite), _random_eigenvalue(rng, T, definite))
Random.rand(rng::AbstractRNG, ::Type{TensISO{2, dim, T}}; definite::Bool = false) where {dim, T} =
    TensISO{dim}(_random_eigenvalue(rng, T, definite))

# The axis of a transversely isotropic draw: the one given, the third vector of
# the frame given (the convention of `frame_from_axis`), or a random one.
_ti_axis(rng, ::Type{T}, ::Nothing, ::Nothing) where {T} = random_axis(rng, T)
_ti_axis(rng, ::Type, axis, ::Nothing) = axis
_ti_axis(rng, ::Type, ::Nothing, frame::OrthonormalBasis{3}) = Tuple(vecbasis(frame, :cov)[:, 3])
_ti_axis(rng, ::Type, axis, frame) =
    throw(ArgumentError("a transversely isotropic tensor takes its axis, or a frame whose third vector is the axis, not both"))

# Transversely isotropic: the 2×2 block coupling the axial and the transverse
# normal directions, and the two eigenvalues of multiplicity two.
Random.rand(rng::AbstractRNG, ::Type{TensTI{4}}; kwargs...) = rand(rng, TensTI{4, Float64, 5}; kwargs...)
function Random.rand(rng::AbstractRNG, ::Type{TensTI{4, T, 5}}; axis = nothing, frame = nothing, definite::Bool = false) where {T}
    B = _random_symmetric(rng, T, Val(2), definite)
    return TensTI{4}(B[1, 1], B[2, 2], B[1, 2], _random_eigenvalue(rng, T, definite), _random_eigenvalue(rng, T, definite), _ti_axis(rng, T, axis, frame))
end
function Random.rand(rng::AbstractRNG, ::Type{TensTI{4, T, 6}}; axis = nothing, frame = nothing, definite::Bool = false) where {T}
    _refuse_definite(definite, "a TensTI{4} of 6 coefficients")
    B = _random_general(rng, T, Val(2))
    return TensTI{4}(B[1, 1], B[2, 2], B[1, 2], B[2, 1], randn(rng, T), randn(rng, T), _ti_axis(rng, T, axis, frame))
end
function Random.rand(rng::AbstractRNG, ::Type{TensTI{4, T, 8}}; axis = nothing, frame = nothing, definite::Bool = false) where {T}
    _refuse_definite(definite, "a TensTI{4} of 8 coefficients")
    B = _random_general(rng, T, Val(2))
    return TensTI{4}(B[1, 1], B[2, 2], B[1, 2], B[2, 1], randn(rng, T), randn(rng, T), randn(rng, T), randn(rng, T), _ti_axis(rng, T, axis, frame))
end
Random.rand(rng::AbstractRNG, ::Type{TensTI{2}}; kwargs...) = rand(rng, TensTI{2, Float64, 2}; kwargs...)
Random.rand(rng::AbstractRNG, ::Type{TensTI{2, T, 2}}; axis = nothing, frame = nothing, definite::Bool = false) where {T} =
    TensTI{2}(_random_eigenvalue(rng, T, definite), _random_eigenvalue(rng, T, definite), _ti_axis(rng, T, axis, frame))
function Random.rand(rng::AbstractRNG, ::Type{TensTI{2, T, 3}}; axis = nothing, frame = nothing, definite::Bool = false) where {T}
    _refuse_definite(definite, "a TensTI{2} of 3 coefficients")
    return TensTI{2}(randn(rng, T), randn(rng, T), randn(rng, T), _ti_axis(rng, T, axis, frame))
end

# Orthotropic: the 3×3 normal block and the three shear eigenvalues, 2C₄₄,
# 2C₅₅, 2C₆₆ in the Kelvin–Mandel matrix.
Random.rand(rng::AbstractRNG, ::Type{TensOrtho{4}}; kwargs...) = rand(rng, TensOrtho{4, Float64, 9}; kwargs...)
function Random.rand(rng::AbstractRNG, ::Type{TensOrtho{4, T, 9}}; frame = random_basis(rng, T), definite::Bool = false) where {T}
    B = _random_symmetric(rng, T, Val(3), definite)
    s = ntuple(_ -> _random_eigenvalue(rng, T, definite) / 2, Val(3))
    return TensOrtho{4}(B[1, 1], B[2, 2], B[3, 3], B[1, 2], B[1, 3], B[2, 3], s..., frame)
end
function Random.rand(rng::AbstractRNG, ::Type{TensOrtho{4, T, 12}}; frame = random_basis(rng, T), definite::Bool = false) where {T}
    _refuse_definite(definite, "a TensOrtho{4} of 12 coefficients")
    B = _random_general(rng, T, Val(3))
    s = ntuple(_ -> randn(rng, T) / 2, Val(3))
    return TensOrtho{4}(B[1, 1], B[2, 2], B[3, 3], B[1, 2], B[1, 3], B[2, 3], s..., B[2, 1], B[3, 1], B[3, 2], frame)
end
Random.rand(rng::AbstractRNG, ::Type{TensOrtho{2}}; kwargs...) = rand(rng, TensOrtho{2, Float64}; kwargs...)
Random.rand(rng::AbstractRNG, ::Type{TensOrtho{2, T}}; frame = random_basis(rng, T), definite::Bool = false) where {T} =
    TensOrtho{2}(ntuple(_ -> _random_eigenvalue(rng, T, definite), Val(3))..., frame)

# Cubic: the three eigenvalues, on 𝕁, 𝔼 and 𝕋.
Random.rand(rng::AbstractRNG, ::Type{TensCubic}; kwargs...) = rand(rng, TensCubic{Float64}; kwargs...)
Random.rand(rng::AbstractRNG, ::Type{TensCubic{T}}; frame = random_basis(rng, T), definite::Bool = false) where {T} =
    TensCubic(ntuple(_ -> _random_eigenvalue(rng, T, definite), Val(3))..., frame)

# Generic: the whole Kelvin–Mandel matrix, or the whole array.
Random.rand(rng::AbstractRNG, ::Type{Tens{order, dim}}; kwargs...) where {order, dim} =
    rand(rng, Tens{order, dim, Float64}; kwargs...)
Random.rand(
    rng::AbstractRNG, ::Type{Tens{4, dim, T}}; symmetry::Val = Val(:major), definite::Bool = false
) where {dim, T} = _random_generic(rng, Val(4), Val(dim), T, symmetry, definite)
Random.rand(
    rng::AbstractRNG, ::Type{Tens{2, dim, T}}; symmetry::Val = Val(:symmetric), definite::Bool = false
) where {dim, T} = _random_generic(rng, Val(2), Val(dim), T, symmetry, definite)

_random_generic(rng, ::Val{4}, ::Val{dim}, ::Type{T}, ::Val{:major}, definite) where {dim, T} =
    from_mandel(_random_symmetric(rng, T, Val(dim * (dim + 1) ÷ 2), definite))
function _random_generic(rng, ::Val{4}, ::Val{dim}, ::Type{T}, ::Val{:minor}, definite) where {dim, T}
    _refuse_definite(definite, "a generic tensor of order 4 with the minor symmetries only")
    return from_mandel(_random_general(rng, T, Val(dim * (dim + 1) ÷ 2)))
end
function _random_generic(rng, ::Val{4}, ::Val{dim}, ::Type{T}, ::Val{:none}, definite) where {dim, T}
    _refuse_definite(definite, "a generic tensor of order 4 without symmetry")
    return Tens(Tensor{4, dim}(randn(rng, T, dim, dim, dim, dim)))
end
_random_generic(rng, ::Val{2}, ::Val{dim}, ::Type{T}, ::Val{:symmetric}, definite) where {dim, T} =
    Tens(SymmetricTensor{2, dim, T}(_random_symmetric(rng, T, Val(dim), definite)))
function _random_generic(rng, ::Val{2}, ::Val{dim}, ::Type{T}, ::Val{:none}, definite) where {dim, T}
    _refuse_definite(definite, "a generic tensor of order 2 without symmetry")
    return Tens(Tensor{2, dim}(randn(rng, T, dim, dim)))
end
_random_generic(rng, ::Val{order}, ::Val, ::Type, ::Val{S}, definite) where {order, S} =
    throw(
    ArgumentError(
        "symmetry = Val(:$S) is not a symmetry of a generic tensor of order $order: " *
            (order == 4 ? "use Val(:major), Val(:minor) or Val(:none)" : "use Val(:symmetric) or Val(:none)")
    )
)

export random_rotation, random_axis, random_basis
