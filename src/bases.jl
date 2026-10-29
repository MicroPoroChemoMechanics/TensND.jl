abstract type AbstractBasis{dim, T <: Number} <: AbstractMatrix{T} end

Base.size(::AbstractBasis{dim}) where {dim} = (dim, dim)

# `Val(:cov)` is spelled out rather than reached through `vecbasis(ℬ, :cov)`:
# the Symbol overload builds the `Val` from a runtime value, so the dispatch was
# dynamic on EVERY scalar access — 67 ns to read one entry.  Since
# `AbstractBasis <: AbstractMatrix`, that cost was paid by every generic
# traversal of a basis.
Base.getindex(ℬ::AbstractBasis, i::Integer, j::Integer) =
    getindex(vecbasis(ℬ, Val(:cov)), i, j)

# Equality compares the primal bases as a whole.  The `AbstractArray` fallback
# would do it entry by entry through the `getindex` above, i.e. `2·dim²`
# separate `vecbasis` calls: measured at 1117 ns for two 3×3 bases, which made
# it the dominant cost of every check of a shared frame — more expensive than
# the tensor algebra it was guarding.  Semantics are unchanged (still an
# elementwise comparison of the primal bases, so a `RotatedBasis` holding the
# identity still equals a `CanonicalBasis`); only the number of dispatches is.
Base.:(==)(A::AbstractBasis{dim}, B::AbstractBasis{dim}) where {dim} =
    vecbasis(A, Val(:cov)) == vecbasis(B, Val(:cov))
Base.eltype(::Type{AbstractBasis{dim, T}}) where {dim, T} = T
get_dim(::AbstractBasis{dim}) where {dim} = dim

# Build a Unicode sub/superscript string from an integer.
# `minus_char` is the sign character; `digit_map` converts a digit 0–9 to its
# Unicode sub/superscript codepoint.
function _scriptnumber(i::Integer, minus_char::Char, digit_map)
    c = i < 0 ? Char[minus_char] : Char[]
    for d in reverse(digits(abs(i)))
        push!(c, digit_map(d))
    end
    return join(c)
end

subscriptnumber(i::Integer) =
    _scriptnumber(i, Char(0x208B), d -> Char(0x2080 + d))

# Superscript digits 0–3 use non-contiguous codepoints; 4–9 are contiguous.
const _SUPERSCRIPT_DIGITS = (
    Char(0x2070), Char(0x00B9), Char(0x00B2), Char(0x00B3),
    Char(0x2074), Char(0x2075), Char(0x2076), Char(0x2077), Char(0x2078), Char(0x2079),
)
superscriptnumber(i::Integer) =
    _scriptnumber(i, Char(0x207B), d -> _SUPERSCRIPT_DIGITS[d + 1])

diago(M::AbstractMatrix) = [M[i, i] for i in 1:min(size(M)...)]

"""
    Basis(v::AbstractMatrix{T}, ::Val{:cov})
    Basis{dim, T<:Number}()
    Basis(θ::T<:Number, ϕ::T<:Number, ψ::T<:Number)

Basis built from a square matrix `v` where columns correspond either to
- primal vectors ie `eᵢ=v[:,i]` if `var=:cov` as by default
- dual vectors ie `eⁱ=v[:,i]` if `var=:cont`.

Basis without any argument refers to the canonical basis (`CanonicalBasis`) in `Rᵈⁱᵐ` (by default `dim=3` and `T=Sym`)

Basis can also be built from Euler angles (`RotatedBasis`) `θ` in 2D and `(θ, ϕ, ψ)` in 3D

The attributes of this object can be obtained by
- `vecbasis(ℬ, :cov)`: square matrix defining the primal basis `eᵢ=e[:,i]`
- `vecbasis(ℬ, :cont)`: square matrix defining the dual basis `eⁱ=E[:,i]`
- `metric(ℬ, :cov)`: square matrix defining the covariant components of the metric tensor `gᵢⱼ=eᵢ⋅eⱼ=g[i,j]`
- `metric(ℬ, :cont)`: square matrix defining the contravariant components of the metric tensor `gⁱʲ=eⁱ⋅eʲ=gⁱʲ[i,j]`

# Examples
```julia
julia> v = Sym[1 0 0; 0 1 0; 0 1 1] ; ℬ = Basis(v)
Basis{3, Sym}
# basis: 3×3 Tensor{2, 3, Sym, 9}:
 1  0  0
 0  1  0
 0  1  1
# dual basis: 3×3 Tensor{2, 3, Sym, 9}:
 1  0   0
 0  1  -1
 0  0   1
# covariant metric tensor: 3×3 SymmetricTensor{2, 3, Sym, 6}:
 1  0  0
 0  2  1
 0  1  1
# contravariant metric tensor: 3×3 SymmetricTensor{2, 3, Sym, 6}:
 1   0   0
 0   1  -1
 0  -1   2

julia> θ, ϕ, ψ = symbols("θ, ϕ, ψ", real = true) ; ℬʳ = Basis(θ, ϕ, ψ) ; display(vecbasis(ℬʳ, :cov))
3×3 Tensor{2, 3, Sym, 9}:
 -sin(ψ)⋅sin(ϕ) + cos(θ)⋅cos(ψ)⋅cos(ϕ)  -sin(ψ)⋅cos(θ)⋅cos(ϕ) - sin(ϕ)⋅cos(ψ)  sin(θ)⋅cos(ϕ)
  sin(ψ)⋅cos(ϕ) + sin(ϕ)⋅cos(θ)⋅cos(ψ)  -sin(ψ)⋅sin(ϕ)⋅cos(θ) + cos(ψ)⋅cos(ϕ)  sin(θ)⋅sin(ϕ)
                        -sin(θ)⋅cos(ψ)                          sin(θ)⋅sin(ψ)         cos(θ)
```
"""
struct Basis{dim, T} <: AbstractBasis{dim, T}
    eᵢ::Matrix{T} # Primal basis `eᵢ=eᵢ[:,i]`
    eⁱ::Matrix{T} # Dual basis `eⁱ=eⁱ[:,i]`
    gᵢⱼ::Symmetric{T, Matrix{T}} # Metric tensor `gᵢⱼ=eᵢ⋅eⱼ=gᵢⱼ[i,j]`
    gⁱʲ::Symmetric{T, Matrix{T}} # Inverse of the metric tensor `gⁱʲ=eⁱ⋅eʲ=gⁱʲ[i,j]`
    function Basis(
            eᵢ::AbstractMatrix{T},
            eⁱ::AbstractMatrix{T},
            gᵢⱼ::AbstractMatrix{T},
            gⁱʲ::AbstractMatrix{T},
        ) where {T}
        dim = size(eᵢ, 1)
        dim == size(eᵢ, 2) || throw(DimensionMismatch("a basis is built from a square matrix"))
        if isidentity(eᵢ)
            return CanonicalBasis{dim, T}()
        elseif isidentity(gᵢⱼ)
            return RotatedBasis(eᵢ)
        elseif isdiag(gᵢⱼ)
            χ = sqrt.(diago(gᵢⱼ))
            return OrthogonalBasis(RotatedBasis(eᵢ .* transpose(inv.(χ))), χ)
        else
            eᵢ = Matrix(eᵢ)
            gᵢⱼ = Symmetric(gᵢⱼ)
            gⁱʲ = Symmetric(gⁱʲ)
            eⁱ = Matrix(eⁱ)
            new{dim, T}(eᵢ, eⁱ, gᵢⱼ, gⁱʲ)
        end
    end
    function Basis(ℬ::AbstractBasis{dim, T}, χᵢ::V) where {dim, T, V}
        Χ = collect(χᵢ)
        invΧ = inv.(Χ)
        all(isone, Χ) && return ℬ
        eᵢ = vecbasis(ℬ, :cov) .* transpose(Χ)
        eⁱ = vecbasis(ℬ, :cont) .* transpose(invΧ)
        gᵢⱼ = Symmetric(Χ .* metric(ℬ, :cov) .* transpose(Χ))
        gⁱʲ = Symmetric(invΧ .* metric(ℬ, :cont) .* transpose(invΧ))
        return new{dim, promote_type(T, eltype(Χ))}(eᵢ, eⁱ, gᵢⱼ, gⁱʲ)
    end
    function Basis(eᵢ::AbstractMatrix{T}, ::Val{:cov}) where {T}
        dim = size(eᵢ, 1)
        dim == size(eᵢ, 2) || throw(DimensionMismatch("a basis is built from a square matrix"))
        if isidentity(eᵢ)
            return CanonicalBasis{dim, T}()
        else
            eᵢ = Matrix(eᵢ)
            gᵢⱼ = tsimplify(Symmetric(transpose(eᵢ) * eᵢ))
            if isidentity(gᵢⱼ)
                return RotatedBasis(eᵢ)
            elseif isdiagonal(gᵢⱼ)
                χ = sqrt.(diago(gᵢⱼ))
                return OrthogonalBasis(RotatedBasis(eᵢ .* transpose(inv.(χ))), χ)
            else
                if T <: SymType
                    gⁱʲ = tsimplify(Symmetric(inv(Matrix(gᵢⱼ))))
                else
                    gⁱʲ = tsimplify(inv(gᵢⱼ))
                end
                eⁱ = tsimplify(eᵢ * transpose(gⁱʲ))
                new{dim, T}(eᵢ, eⁱ, gᵢⱼ, gⁱʲ)
            end
        end
    end
    function Basis(eⁱ::AbstractMatrix{T}, ::Val{:cont}) where {T}
        dim = size(eⁱ, 1)
        dim == size(eⁱ, 2) || throw(DimensionMismatch("a basis is built from a square matrix"))
        if isidentity(eⁱ)
            return CanonicalBasis{dim, T}()
        else
            eⁱ = Matrix(eⁱ)
            gⁱʲ = tsimplify(Symmetric(transpose(eⁱ) * eⁱ))
            if isidentity(gⁱʲ)
                return RotatedBasis(eⁱ)
            elseif isdiagonal(gⁱʲ)
                uχ = inv.(sqrt.(diago(gⁱʲ)))
                return OrthogonalBasis(RotatedBasis(eⁱ .* transpose(uχ)), uχ)
            else
                if T <: SymType
                    gᵢⱼ = tsimplify(Symmetric(inv(Matrix(gⁱʲ))))
                else
                    gᵢⱼ = tsimplify(inv(gⁱʲ))
                end
                eᵢ = tsimplify(eⁱ * transpose(gᵢⱼ))
                new{dim, T}(eᵢ, eⁱ, gᵢⱼ, gⁱʲ)
            end
        end
    end
    Basis(v::AbstractMatrix{T}, var::Symbol) where {T} = _by_variance(Basis, var, v)
    # An `AbstractBasis` is an `AbstractMatrix`, so `Basis(ℬ, χᵢ)` above — which
    # rescales a basis by a tuple of factors — is weighed against the three
    # constructors that read a variance marker as their second argument. The
    # second argument tells them apart: a `Symbol` or a `Val` asks to reinterpret
    # the matrix, never to rescale a basis. These three say so.
    Basis(ℬ::AbstractBasis{dim, T}, var::Symbol) where {dim, T} =
        _by_variance(Basis, var, Matrix(ℬ))
    Basis(ℬ::AbstractBasis{dim, T}, ::Val{:cov}) where {dim, T} =
        Basis(Matrix(ℬ), Val(:cov))
    Basis(ℬ::AbstractBasis{dim, T}, ::Val{:cont}) where {dim, T} =
        Basis(Matrix(ℬ), Val(:cont))
    Basis(v::AbstractMatrix{T}) where {T} = Basis(v, :cov)
    Basis(θ::T1, ϕ::T2, ψ::T3 = 0) where {T1, T2, T3} = RotatedBasis(θ, ϕ, ψ)
    Basis(θ::T) where {T} = RotatedBasis(θ)
    Basis{dim, T}() where {dim, T} = CanonicalBasis{dim, T}()
    Basis() = CanonicalBasis()
end


"""
    CanonicalBasis{dim, T}

Canonical basis of dimension `dim` (default: 3) and type `T` (default: Sym)

The attributes of this object can be obtained by
- `vecbasis(ℬ, :cov)`: square matrix defining the primal basis `eᵢ=e[:,i]=δᵢⱼ`
- `vecbasis(ℬ, :cont)`: square matrix defining the dual basis `eⁱ=E[:,i]=δᵢⱼ`
- `metric(ℬ, :cov)`: square matrix defining the covariant components of the metric tensor `gᵢⱼ=eᵢ⋅eⱼ=g[i,j]=δᵢⱼ`
- `metric(ℬ, :cont)`: square matrix defining the contravariant components of the metric tensor `gⁱʲ=eⁱ⋅eʲ=gⁱʲ[i,j]=δᵢⱼ`

# Examples
```julia
julia> ℬ = CanonicalBasis()
CanonicalBasis{3, Sym}
# basis: 3×3 TensND.LazyIdentity{3, Sym}:
 1  0  0
 0  1  0
 0  0  1
# dual basis: 3×3 TensND.LazyIdentity{3, Sym}:
 1  0  0
 0  1  0
 0  0  1
# covariant metric tensor: 3×3 TensND.LazyIdentity{3, Sym}:
 1  0  0
 0  1  0
 0  0  1
# contravariant metric tensor: 3×3 TensND.LazyIdentity{3, Sym}:
 1  0  0
 0  1  0
 0  0  1

julia> ℬ₂ = CanonicalBasis{2, Float64}()
CanonicalBasis{2, Float64}
# basis: 2×2 TensND.LazyIdentity{2, Float64}:
 1.0  0.0
 0.0  1.0
# dual basis: 2×2 TensND.LazyIdentity{2, Float64}:
 1.0  0.0
 0.0  1.0
# covariant metric tensor: 2×2 TensND.LazyIdentity{2, Float64}:
 1.0  0.0
 0.0  1.0
# contravariant metric tensor: 2×2 TensND.LazyIdentity{2, Float64}:
 1.0  0.0
 0.0  1.0
```
"""
struct CanonicalBasis{dim, T} <: AbstractBasis{dim, T} end
CanonicalBasis() = CanonicalBasis{3, Sym}()


_check_orthogonal(R::AbstractMatrix{T}) where {T} = _check_orthogonal(R, Val(is_hard_numeric(T)))
_check_orthogonal(::AbstractMatrix, ::Val{false}) = nothing
function _check_orthogonal(R::AbstractMatrix, ::Val{true})
    transpose(R) * R ≈ I || throw(
        ArgumentError(
            "RotatedBasis needs an orthogonal matrix, got one with RᵀR = " *
                "$(transpose(R) * R); use `Basis(R)` for a general basis"
        )
    )
    return nothing
end

"""
    RotatedBasis(θ::T<:Number, ϕ::T<:Number, ψ::T<:Number)
    RotatedBasis(θ::T<:Number)

Orthonormal basis of dimension `dim` (default: 3) and type `T` (default: Sym) built from Euler angles `θ` in 2D and `(θ, ϕ, ψ)` in 3D

# Examples
```julia
julia> θ, ϕ, ψ = symbols("θ, ϕ, ψ", real = true) ; ℬʳ = RotatedBasis(θ, ϕ, ψ) ; display(vecbasis(ℬʳ, :cov))
3×3 Tensor{2, 3, Sym, 9}:
 -sin(ψ)⋅sin(ϕ) + cos(θ)⋅cos(ψ)⋅cos(ϕ)  -sin(ψ)⋅cos(θ)⋅cos(ϕ) - sin(ϕ)⋅cos(ψ)  sin(θ)⋅cos(ϕ)
  sin(ψ)⋅cos(ϕ) + sin(ϕ)⋅cos(θ)⋅cos(ψ)  -sin(ψ)⋅sin(ϕ)⋅cos(θ) + cos(ψ)⋅cos(ϕ)  sin(θ)⋅sin(ϕ)
                        -sin(θ)⋅cos(ψ)                          sin(θ)⋅sin(ψ)         cos(θ)
```
"""
struct RotatedBasis{dim, T} <: AbstractBasis{dim, T}
    eᵢ::Matrix{T} # Primal basis `eᵢ=e[:,i]`
    eⁱ::Matrix{T} # Dual basis `eⁱ=E[:,i]`
    # Read by `angles(ℬ)` only, never on a computation's path: a field of
    # concrete type would be copied with every tensor that carries the basis.
    angles::NamedTuple
    # The dimension as a type parameter, so that a caller who knows it (a basis
    # permuted or converted from another) gets a concrete type back. `check`
    # is false only where `R` is orthogonal by construction.
    function RotatedBasis(R::AbstractMatrix{T}, ::Val{dim}, check::Bool = true) where {T <: Number, dim}
        size(R) == (dim, dim) ||
            throw(DimensionMismatch("RotatedBasis: a $(dim)×$(dim) matrix is needed, got $(size(R))"))
        isidentity(R) && return CanonicalBasis{dim, T}()
        # The dual basis is taken to be `R` itself, which holds only for an
        # orthogonal `R`. Checked where a comparison means something; a
        # symbolic `R` is the caller's responsibility.
        check && _check_orthogonal(R)
        eᵢ = eⁱ = Matrix(R)
        return new{dim, T}(eᵢ, eⁱ, angles(R))
    end
    # A rotated basis even when `R` is the identity: for a frame built inside the
    # package, whose type must not depend on its value (a `RotatedBasis` holding
    # the identity equals the `CanonicalBasis`, see `==` below).
    function RotatedBasis(R::AbstractMatrix{T}, ::Val{dim}, ::Val{:typed}) where {T <: Number, dim}
        _check_orthogonal(R)
        eᵢ = eⁱ = Matrix(R)
        return new{dim, T}(eᵢ, eⁱ, angles(R))
    end
    function RotatedBasis(θ::T1, ϕ::T2, ψ::T3 = 0) where {T1 <: Number, T2 <: Number, T3 <: Number}
        Rz = RotZYZ(ϕ, θ, ψ)
        R = Matrix(Rz)
        isidentity(R) && return CanonicalBasis{3, promote_type(T1, T2, T3)}()
        # The type of the rotation, not of the angles: integer angles give
        # floating-point cosines.
        return new{3, eltype(R)}(R, R, angles(Rz))
    end
    function RotatedBasis(θ::T) where {T <: Number}
        cθ = cos(θ)
        sθ = sin(θ)
        eᵢ = [cθ -sθ; sθ cθ]
        isidentity(eᵢ) && return CanonicalBasis{2, T}()
        return new{2, eltype(eᵢ)}(eᵢ, eᵢ, angles(eᵢ))
    end
    function RotatedBasis(θ::T) where {T <: SymType}
        cθ = cos(θ)
        sθ = sin(θ)
        eᵢ = [cθ -sθ; sθ cθ]
        isidentity(eᵢ) && return CanonicalBasis{2, T}()
        return new{2, eltype(eᵢ)}(eᵢ, eᵢ, (θ = θ,))
    end
end


RotatedBasis(R::AbstractMatrix{<:Number}) = RotatedBasis(R, Val(size(R, 1)))

const OrthonormalBasis{dim, T} = Union{CanonicalBasis{dim, T}, RotatedBasis{dim, T}}

"""
    OrthogonalBasis{dim,T}

Orthogonal but **not** normalized basis: an orthonormal basis scaled by one
factor `χᵢ` per direction.

Its metric is diagonal, `gᵢⱼ = diag(χᵢ²)`, and its dual vectors are
`𝐞ⁱ = 𝐞ᵢ/χᵢ²`, so variance still matters but only diagonally. This is exactly
the situation of the **natural basis** of a curvilinear coordinate system, where
the `χᵢ` are the Lamé coefficients — see
[Curvilinear differential calculus](@ref th-curvilinear).

Built by `Basis(ℬ, χᵢ)` from an orthonormal basis and a tuple of scaling
factors. Normalizing it returns the orthonormal basis it came from.
"""
struct OrthogonalBasis{dim, T, P <: OrthonormalBasis{dim, T}} <: AbstractBasis{dim, T}
    parent::P
    λ::Vector{T}
    eᵢ::Matrix{T}
    eⁱ::Matrix{T}
    gᵢⱼ::Diagonal{T, Vector{T}}
    gⁱʲ::Diagonal{T, Vector{T}}
    function OrthogonalBasis(parent::OrthonormalBasis{dim, T0}, λ::Vector) where {dim, T0}
        # Scale factors of a wider type than the frame (duals on a `Float64`
        # frame) widen the basis rather than being converted to the frame's.
        T = promote_type(T0, eltype(λ))
        parent = _with_eltype(T, parent)
        λ = T.(λ)
        if all(isone.(λ))
            return parent
        else
            eᵢ = [λ[j] * parent[i, j] for i in 1:dim, j in 1:dim]
            eⁱ = [parent[i, j] / λ[j] for i in 1:dim, j in 1:dim]
            return new{dim, T, typeof(parent)}(parent, λ, eᵢ, eⁱ, Diagonal(λ .^ 2), Diagonal(inv.(λ) .^ 2))
        end
    end
end

# Scaling an orthonormal basis gives an orthogonal one.
function Basis(ℬ::OrthonormalBasis{dim}, χᵢ::Union{AbstractVector, Tuple}) where {dim}
    Χ = collect(χᵢ)
    return all(isone, Χ) ? ℬ : OrthogonalBasis(ℬ, Χ)
end

# The same orthonormal basis with entries of type `T`.
_with_eltype(::Type{T}, ::CanonicalBasis{dim}) where {dim, T} = CanonicalBasis{dim, T}()
_with_eltype(::Type{T}, ℬ::RotatedBasis{dim, T}) where {dim, T} = ℬ
_with_eltype(::Type{T}, ℬ::RotatedBasis{dim}) where {T, dim} =
    RotatedBasis(Matrix{T}(vecbasis(ℬ, :cov)), Val(dim), false)

relevant_OrthonormalBasis(ℬ::OrthogonalBasis) = ℬ.parent
relevant_OrthonormalBasis(ℬ::OrthonormalBasis) = ℬ
relevant_OrthonormalBasis(::Basis{dim, T}) where {dim, T} = CanonicalBasis{dim, T}()


"""
    CylindricalBasis(θ)

Local orthonormal basis `(𝐞ʳ, 𝐞ᶿ, 𝐞ᶻ)` of the cylindrical system at azimuth `θ`,
i.e. `RotatedBasis(0, θ, 0)`.

Convenience constructor for a one-off local frame; for differential calculus use
[`coorsys_cylindrical`](@ref), which also carries the Lamé coefficients and the
Christoffel symbols.
"""
@inline CylindricalBasis(θ) = RotatedBasis(0, θ, 0)

"""
    SphericalBasis(θ, ϕ)

Local orthonormal basis `(𝐞ᶿ, 𝐞ᵠ, 𝐞ʳ)` of the spherical system at polar angle
`θ` and azimuth `ϕ`, i.e. `RotatedBasis(θ, ϕ, 0)`.

Note the ordering: the radial vector comes **last**, so that `θ = ϕ = 0` returns
the canonical basis in the canonical order. See
[Rotations and Euler angles](@ref th-rotations).

Convenience constructor for a one-off local frame; for differential calculus use
[`coorsys_spherical`](@ref).
"""
@inline SphericalBasis(θ, ϕ) = RotatedBasis(θ, ϕ, 0)

"""
    AllOrthogonalBasis{dim,T}

Union of every basis whose metric is **diagonal**: the orthonormal ones
([`CanonicalBasis`](@ref), [`RotatedBasis`](@ref)) and the merely orthogonal
[`OrthogonalBasis`](@ref).

Used to dispatch the component conversions that need no dense metric.
"""
const AllOrthogonalBasis{dim, T} = Union{OrthonormalBasis{dim, T}, OrthogonalBasis{dim, T}}

angles(M::AbstractMatrix{T}, ::Val{2}) where {T} =
    (θ = atan(M[2, 1] - M[1, 2], M[1, 1] + M[2, 2]),)
function angles(M::AbstractMatrix{T}, ::Val{3}) where {T}
    R = RotZYZ(M)
    return (θ = R.theta2, ϕ = R.theta1, ψ = R.theta3)
end

"""
    angles(M::AbstractMatrix{T})

Determine the Euler angles corresponding to the input matrix supposed to be a rotation matrix or at least a similarity

# Examples
```julia
julia> θ, ϕ, ψ = symbols("θ, ϕ, ψ", real = true) ; ℬʳ = RotatedBasis(θ, ϕ, ψ) ; display(vecbasis(ℬʳ, :cov))
3×3 Tensor{2, 3, Sym, 9}:
 -sin(ψ)⋅sin(ϕ) + cos(θ)⋅cos(ψ)⋅cos(ϕ)  -sin(ψ)⋅cos(θ)⋅cos(ϕ) - sin(ϕ)⋅cos(ψ)  sin(θ)⋅cos(ϕ)
  sin(ψ)⋅cos(ϕ) + sin(ϕ)⋅cos(θ)⋅cos(ψ)  -sin(ψ)⋅sin(ϕ)⋅cos(θ) + cos(ψ)⋅cos(ϕ)  sin(θ)⋅sin(ϕ)
                        -sin(θ)⋅cos(ψ)                          sin(θ)⋅sin(ψ)         cos(θ)

julia> angles(ℬʳ)
(θ = θ, ϕ = ϕ, ψ = ψ)
```
"""
angles(M::AbstractMatrix{T}) where {T} = angles(M, Val(size(M, 1)))
angles(ℬ::RotatedBasis) = ℬ.angles

angles(v::AbstractVector{T}, ::Val{2}) where {T} = (θ = atan(v[2], v[1]),)
angles(v::AbstractVector{T}, ::Val{3}) where {T} =
    (θ = atan(√(v[1]^2 + v[2]^2), v[3]), ϕ = atan(v[2], v[1]))
angles(v::AbstractVector{T}) where {T} = angles(v, Val(size(v, 1)))


invvar(::Val{:cov}) = :cont
invvar(::Val{:cont}) = :cov
invvar(var::Symbol) = _by_variance(invvar, var)

# A variance is a run-time `Symbol`: branching on its two values keeps the call
# static, where `Val(var)` was dispatched at run time on every call.
@inline function _by_variance(f, var::Symbol, args...)
    var === :cov && return f(args..., Val(:cov))
    var === :cont && return f(args..., Val(:cont))
    throw(ArgumentError("a variance is :cov or :cont, not :$var"))
end


"""
    vecbasis(ℬ::AbstractBasis, var = :cov)

Return the primal (if `var = :cov`) or dual (if `var = :cont`) basis
"""
vecbasis(ℬ::AbstractBasis, ::Val{:cov}) = ℬ.eᵢ
vecbasis(ℬ::AbstractBasis, ::Val{:cont}) = ℬ.eⁱ
vecbasis(::CanonicalBasis{dim, T}, ::Val{:cov}) where {dim, T} = Id2{dim, T}()
vecbasis(::CanonicalBasis{dim, T}, ::Val{:cont}) where {dim, T} = Id2{dim, T}()

vecbasis(ℬ::AbstractBasis, var::Symbol) = _by_variance(vecbasis, var, ℬ)
vecbasis(ℬ::AbstractBasis) = vecbasis(ℬ, Val(:cov))
vecbasis(ℬ::AbstractBasis, i::Integer, j::Integer, var::Symbol = :cov) = vecbasis(ℬ, var)[i, j]
vecbasis(ℬ::AbstractBasis, i::Integer, var::Symbol = :cov) = vecbasis(ℬ, var)[:, i]

strvecbasis(::AbstractBasis, i::Integer, ::Val{:cov}; vec = "𝐞") = vec * subscriptnumber(i)
strvecbasis(::AbstractBasis, i::Integer, ::Val{:cont}; vec = "𝐞") = vec * superscriptnumber(i)
strvecbasis(ℬ::AbstractBasis, i::Integer, var::Symbol = :cov; vec = "𝐞") =
    _by_variance((ℬ, i, v) -> strvecbasis(ℬ, i, v; vec), var, ℬ, i)

const dsubscriptchar = Dict(
    "a" => "ₐ",
    "e" => "ₑ",
    "h" => "ₕ",
    "i" => "ᵢ",
    "j" => "ⱼ",
    "k" => "ₖ",
    "l" => "ₗ",
    "m" => "ₘ",
    "n" => "ₙ",
    "o" => "ₒ",
    "p" => "ₚ",
    "r" => "ᵣ",
    "s" => "ₛ",
    "t" => "ₜ",
    "u" => "ᵤ",
    "v" => "ᵥ",
    "x" => "ₓ",
    "β" => "ᵦ",
    "γ" => "ᵧ",
    "ρ" => "ᵨ",
    "ϕ" => "ᵩ",
    "χ" => "ᵪ",
)

const dsuperscriptchar = Dict(
    "a" => "ᵃ",
    "b" => "ᵇ",
    "c" => "ᶜ",
    "d" => "ᵈ",
    "e" => "ᵉ",
    "f" => "ᶠ",
    "g" => "ᵍ",
    "h" => "ʰ",
    "i" => "ⁱ",
    "j" => "ʲ",
    "k" => "ᵏ",
    "l" => "ˡ",
    "m" => "ᵐ",
    "n" => "ⁿ",
    "o" => "ᴼ",
    "p" => "ᵖ",
    "r" => "ʳ",
    "s" => "ˢ",
    "t" => "ᵗ",
    "u" => "ᶸ",
    "v" => "ᵛ",
    "w" => "ʷ",
    "x" => "ˣ",
    "y" => "ʸ",
    "z" => "ᶻ",
    "β" => "ᵝ",
    "γ" => "ᵞ",
    "ε" => "ᵋ",
    "θ" => "ᶿ",
    "ι" => "ᶥ",
    "ϕ" => "ᵠ",
    "χ" => "ᵡ",
)

subscriptchar(s::String) = s ∈ keys(dsubscriptchar) ? dsubscriptchar[s] : s ∈ keys(dsuperscriptchar) ? dsuperscriptchar[s] : s
superscriptchar(s::String) = s ∈ keys(dsuperscriptchar) ? dsuperscriptchar[s] : s ∈ keys(dsubscriptchar) ? dsubscriptchar[s] : s


strvecbasis(::AbstractBasis, i::AbstractString, ::Val{:cov}; vec = "𝐞") = vec * subscriptchar(i)
strvecbasis(::AbstractBasis, i::AbstractString, ::Val{:cont}; vec = "𝐞") = vec * superscriptchar(i)
strvecbasis(ℬ::AbstractBasis, i::AbstractString, var::Symbol = :cov; vec = "𝐞") =
    _by_variance((ℬ, i, v) -> strvecbasis(ℬ, i, v; vec), var, ℬ, i)


"""
    metric(ℬ::AbstractBasis, var = :cov)

Return the covariant (if `var = :cov`) or contravariant (if `var = :cont`) metric matrix
"""
metric(ℬ::AbstractBasis, ::Val{:cov}) = ℬ.gᵢⱼ
metric(ℬ::AbstractBasis, ::Val{:cont}) = ℬ.gⁱʲ
metric(::OrthonormalBasis{dim, T}, ::Val{:cov}) where {dim, T} = Id2{dim, T}()
metric(::OrthonormalBasis{dim, T}, ::Val{:cont}) where {dim, T} = Id2{dim, T}()

metric(ℬ::AbstractBasis, var::Symbol) = _by_variance(metric, var, ℬ)
metric(ℬ::AbstractBasis) = metric(ℬ, Val(:cov))
metric(ℬ::AbstractBasis, i::Integer, j::Integer, var::Symbol = :cov) = metric(ℬ, var)[i, j]

"""
    normalize(ℬ::AbstractBasis, var::Symbol = :cov)

Build a basis after normalization of column vectors of input matrix `v` where columns define either
- primal vectors ie `eᵢ=v[:,i]/norm(v[:,i])` if `var = :cov` as by default
- dual vector ie `eⁱ=v[:,i]/norm(v[:,i])` if `var = :cont`.
"""
function LinearAlgebra.normalize(ℬ::AbstractBasis, var::Symbol = :cov)
    # `var` is a `Symbol`: untyped, the method tied with LinearAlgebra's
    # `normalize(::AbstractArray, ::Real)` on a basis and a number.
    w = copy(vecbasis(ℬ, var))
    for i in 1:size(w, 2)
        w[:, i] /= norm(w[:, i])
    end
    return Basis(w, var)
end

"""
    isorthogonal(ℬ::AbstractBasis)

Check whether the basis `ℬ` is orthogonal
"""
isorthogonal(ℬ::AbstractBasis) = isdiagonal(metric(ℬ))

isorthogonal(::OrthonormalBasis) = true

isorthogonal(::OrthogonalBasis) = true

"""
    isorthonormal(ℬ::AbstractBasis)

Check whether the basis `ℬ` is orthonormal
"""
isorthonormal(ℬ::AbstractBasis) = isidentity(metric(ℬ))

isorthonormal(::OrthonormalBasis) = true

for OP in (:(tsimplify), :(tfactor), :(tsubs), :(tdiff), :(ttrigsimp), :(texpand_trig))
    @eval $OP(b::AbstractBasis{dim, <:Sym}, args...; kwargs...) where {dim} =
        Basis($OP(b.eᵢ, args...; kwargs...))
    @eval $OP(b::CanonicalBasis{dim, <:Sym}, args...; kwargs...) where {dim} = b
end
for OP in (:(tsimplify), :(tsubs), :(tdiff))
    @eval $OP(b::AbstractBasis{dim, Num}, args...; kwargs...) where {dim} =
        Basis($OP(b.eᵢ, args...; kwargs...))
    @eval $OP(b::CanonicalBasis{dim, Num}, args...; kwargs...) where {dim} = b
end


#####################
# Display Functions #
#####################
function Base.show(io::IO, ::MIME"text/plain", ℬ::AbstractBasis)
    println(io, typeof(ℬ))
    for (label, M) in (
            ("→ basis: ", vecbasis(ℬ, :cov)), ("\n→ dual basis: ", vecbasis(ℬ, :cont)),
            ("\n→ covariant metric tensor: ", metric(ℬ, :cov)),
            ("\n→ contravariant metric tensor: ", metric(ℬ, :cont)),
        )
        print(io, label)
        show(io, MIME"text/plain"(), M)
    end
    return nothing
end

"""
    euler_basis(angles::Tuple; T = Float64) -> AbstractBasis{3}

The orthonormal basis of the ZYZ Euler angles `angles`, zero to three of them:
missing trailing angles are zero, and the angles are promoted to one type —
`Int` or `Irrational` to a float, a dual staying a dual, a symbolic angle
keeping its type. When every angle is zero the basis is the canonical one, of
element type `Float64` for real numbers (duals included: its entries are exact
zeros and ones) and of the symbolic type otherwise, promoted with `T` — the
element type of the data the basis will carry.

```julia
julia> euler_basis(()) isa CanonicalBasis{3, Float64}
true

julia> euler_basis((π / 2,)) == RotatedBasis(π / 2, 0.0, 0.0)
true
```
"""
function euler_basis(angles::Tuple{Vararg{Number}}; T::Type = Float64)
    n = length(angles)
    n ≤ 3 || throw(ArgumentError("euler_basis: at most three ZYZ angles, got $n: $angles"))
    θ = _euler_triple(angles)
    return all(iszero, θ) ? CanonicalBasis{3, canonical_eltype(promote_type(T, eltype(θ)))}() :
        RotatedBasis(θ...)
end

# Three angles of one type: a float for numbers (`float` turns an `Int` or an
# `Irrational` into one and leaves a dual a dual), the symbolic type otherwise
# (`Symbolics.Num` subtypes `Real` without converting from `Float64`).
function _euler_triple(angles::Tuple)
    T0 = promote_type(map(typeof, angles)...)
    T = is_hard_numeric(T0) ? float(T0) : T0
    return ntuple(i -> i ≤ length(angles) ? convert(T, angles[i]) : zero(T), 3)
end
_euler_triple(::Tuple{}) = (0.0, 0.0, 0.0)

"""
    canonical_eltype(T) -> Type

The element type of the canonical basis that carries data of type `T`, the rule
[`euler_basis`](@ref) follows: `Float64` for any real number (duals included: the
entries of the canonical basis are exact zeros and ones), `T` itself for a
symbolic or complex type, whose exact `0` and `1` add no float to an exact result.

```julia
julia> canonical_eltype(Int)
Float64
```
"""
canonical_eltype(::Type{T}) where {T <: Real} = Float64
canonical_eltype(::Type{T}) where {T <: Number} = T

"""
    permute_basis(ℬ, σ; negate = dim) -> RotatedBasis

The orthonormal basis whose `k`-th vector is the `σ[k]`-th vector of `ℬ`, with
the orientation of `ℬ`: an odd permutation reverses the orientation, so the
vector of index `negate` (the last one by default) is reversed as well. A shape
invariant under the reflections of its own axes (an ellipsoid, a crack) is the
same shape in the permuted basis, its semi-axes permuted by `σ`. A basis that
is not orthonormal is refused, as by [`RotatedBasis`](@ref).

```julia
julia> permute_basis(RotatedBasis(0.3, 0.2, 0.1), (2, 1, 3)) ;   # e₂, e₁, -e₃
```
"""
function permute_basis(ℬ::AbstractBasis{dim, T}, σ::NTuple{dim, Int}; negate::Int = dim) where {dim, T}
    isperm(σ) || throw(ArgumentError("permute_basis: $σ is not a permutation of 1:$dim"))
    1 ≤ negate ≤ dim || throw(ArgumentError("permute_basis: no vector of index $negate in dimension $dim"))
    odd = isodd(count(σ[i] > σ[j] for i in 1:dim for j in (i + 1):dim))
    E = vecbasis(ℬ, :cov)
    s(k) = odd && k == negate ? -one(T) : one(T)
    # A typed comprehension rather than `hcat`, which Symbolics overloads. The
    # vectors of an orthonormal basis, permuted and signed, are orthonormal: the
    # check is left to the bases that are not orthonormal by type.
    return RotatedBasis(T[s(k) * E[i, σ[k]] for i in 1:dim, k in 1:dim], Val(dim), !isorthonormal(ℬ))
end

export Basis, CanonicalBasis, RotatedBasis, CylindricalBasis, SphericalBasis, OrthogonalBasis, AllOrthogonalBasis
export get_dim, vecbasis, metric, angles, isorthogonal, isorthonormal, permute_basis, euler_basis, canonical_eltype
