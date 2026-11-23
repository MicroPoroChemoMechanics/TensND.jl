"""
    levi_civita(T::Type{<:Number} = default_number_type())

Builds an Array{T,3} of Levi-Civita Symbol `ϵᵢⱼₖ = (i-j) (j-k) (k-i) / 2`

# Examples
```jldoctest
julia> ε = levi_civita(Sym)
3×3×3 Array{Sym{PythonCall.Py}, 3}:
[:, :, 1] =
 0   0  0
 0   0  1
 0  -1  0

[:, :, 2] =
 0  0  -1
 0  0   0
 1  0   0

[:, :, 3] =
  0  1  0
 -1  0  0
  0  0  0
``` 
"""
levi_civita(T::Type{<:Number} = default_number_type()) = [T(T((i - j) * (j - k) * (k - i)) / T(2)) for i in 1:3, j in 1:3, k in 1:3]

"""
    𝐞(i::Integer, dim::Int = 3, T::Type{<:Number} = default_number_type())

Vector of the canonical basis

# Examples
```jldoctest
julia> 𝐞(1)
TensND.TensCanonical{1, 3, Sym{PythonCall.Py}, Vec{3, Sym{PythonCall.Py}}}
→ array: 3-element Vec{3, Sym{PythonCall.Py}}:
 1
 0
 0
→ basis: 3×3 TensND.Id2{3, Sym{PythonCall.Py}}:
 1  ⋅  ⋅
 ⋅  1  ⋅
 ⋅  ⋅  1
``` 
"""
𝐞(::Val{i}, ::Val{dim} = Val(3), ::Val{T} = Val(default_number_type())) where {i, dim, T <: Number} =
    Tens(Vec{dim}(j -> j == i ? one(T) : zero(T)))

"""
    init_cartesian(coords = symbols("x y z", real = true))

Returns the coordinates, unit vectors and basis of the cartesian basis

# Examples
```jldoctest
julia> coords, vectors, ℬ = init_cartesian() ; x, y, z = coords ; 𝐞₁, 𝐞₂, 𝐞₃ = vectors ;
``` 
"""
init_cartesian(coords) = Tuple(coords),
    ntuple(i -> 𝐞(Val(i), Val(length(coords)), Val(eltype(coords))), length(coords)),
    CanonicalBasis{length(coords), eltype(coords)}()

init_cartesian(dim::Integer) = init_cartesian(Val(dim))
# The default coordinates are SymPy symbols, defined by the extension that
# SymPyPythonCall loads (`init_cartesian(::Val{2})`, `::Val{3}`).
init_cartesian(::Val{dim}) where {dim} = throw(
    ArgumentError(
        "init_cartesian($dim): the default coordinates are SymPy symbols, for dimensions 2 " *
            "and 3; load SymPyPythonCall (`using SymPyPythonCall`), or pass the coordinates"
    )
)


"""
    𝐞ᵖ(i::Integer, θ::T = zero(default_number_type()); canonical = false)

Vector of the polar basis

# Examples
```jldoctest
julia> θ = symbols("θ", real = true) ;

julia> 𝐞ᵖ(1, θ)
TensND.TensRotated{1, 2, Sym{PythonCall.Py}, Vec{2, Sym{PythonCall.Py}}, RotatedBasis{2, Sym{PythonCall.Py}}}
→ array: 2-element Vec{2, Sym{PythonCall.Py}}:
 1
 0
→ basis: 2×2 Matrix{Sym{PythonCall.Py}}:
 cos(θ)  -sin(θ)
 sin(θ)   cos(θ)
``` 
"""
𝐞ᵖ(::Val{1}, θ::T = zero(default_number_type()); canonical = false) where {T <: Number} =
    canonical ? Tens(Vec{2}([cos(θ), sin(θ)])) :
    Tens(Vec{2}([one(T), zero(T)]), Basis(θ))
𝐞ᵖ(::Val{2}, θ::T = zero(default_number_type()); canonical = false) where {T <: Number} =
    canonical ? Tens(Vec{2}([-sin(θ), cos(θ)])) :
    Tens(Vec{2}([zero(T), one(T)]), Basis(θ))

"""
    init_polar(coords = (symbols("r θ", real = true)); canonical = false)

Returns the coordinates, base vectors and basis of the polar basis

# Examples
```jldoctest
julia> coords, vectors, ℬᵖ = init_polar() ; r, θ = coords ; 𝐞ʳ, 𝐞ᶿ = vectors ;
``` 
"""
init_polar(coords; canonical = false) = Tuple(coords),
    ntuple(i -> 𝐞ᵖ(Val(i), coords[2]; canonical = canonical), 2),
    Basis(coords[2])

"""
    𝐞ᶜ(i::Integer, θ::T = zero(default_number_type()); canonical = false)

Vector of the cylindrical basis

# Examples
```jldoctest
julia> θ = symbols("θ", real = true) ;

julia> 𝐞ᶜ(1, θ)
TensND.TensRotated{1, 3, Sym{PythonCall.Py}, Vec{3, Sym{PythonCall.Py}}, RotatedBasis{3, Sym{PythonCall.Py}}}
→ array: 3-element Vec{3, Sym{PythonCall.Py}}:
 1
 0
 0
→ basis: 3×3 Matrix{Sym{PythonCall.Py}}:
 cos(θ)  -sin(θ)  0
 sin(θ)   cos(θ)  0
      0        0  1
``` 
"""
𝐞ᶜ(::Val{1}, θ::T = zero(default_number_type()); canonical = false) where {T <: Number} =
    canonical ? Tens(Vec{3}([cos(θ), sin(θ), zero(T)])) :
    Tens(Vec{3}([one(T), zero(T), zero(T)]), CylindricalBasis(θ))
𝐞ᶜ(::Val{2}, θ::T = zero(default_number_type()); canonical = false) where {T <: Number} =
    canonical ? Tens(Vec{3}([-sin(θ), cos(θ), zero(T)])) :
    Tens(Vec{3}([zero(T), one(T), zero(T)]), CylindricalBasis(θ))
𝐞ᶜ(::Val{3}, θ::T = zero(default_number_type()); canonical = false) where {T <: Number} =
    canonical ? Tens(Vec{3}([zero(T), zero(T), one(T)])) :
    Tens(Vec{3}([zero(T), zero(T), one(T)]), CylindricalBasis(θ))

"""
    init_cylindrical(coords = (symbols("r", positive = true), symbols("θ z", real = true)...); canonical = false)

Returns the coordinates, base vectors and basis of the cylindrical basis

# Examples
```jldoctest
julia> coords, vectors, ℬᶜ = init_cylindrical() ; r, θ, z = coords ; 𝐞ʳ, 𝐞ᶿ, 𝐞ᶻ = vectors ;
``` 
"""
init_cylindrical(coords; canonical = false) = Tuple(coords),
    ntuple(i -> 𝐞ᶜ(Val(i), coords[2]; canonical = canonical), 3),
    CylindricalBasis(coords[2])

"""
    𝐞ˢ(i::Integer, θ::T = zero(default_number_type()), ϕ::T = zero(default_number_type()), ψ::T = zero(default_number_type()); canonical = false)

Vector of the basis rotated with the 3 Euler angles `θ, ϕ, ψ` (spherical if `ψ=0`)

# Examples
```jldoctest
julia> θ, ϕ, ψ = symbols("θ, ϕ, ψ", real = true) ;

julia> 𝐞ˢ(1, θ, ϕ, ψ)
TensND.TensRotated{1, 3, Sym{PythonCall.Py}, Vec{3, Sym{PythonCall.Py}}, RotatedBasis{3, Sym{PythonCall.Py}}}
→ array: 3-element Vec{3, Sym{PythonCall.Py}}:
 1
 0
 0
→ basis: 3×3 Matrix{Sym{PythonCall.Py}}:
 -sin(ψ)⋅sin(ϕ) + cos(θ)⋅cos(ψ)⋅cos(ϕ)  …  sin(θ)⋅cos(ϕ)
  sin(ψ)⋅cos(ϕ) + sin(ϕ)⋅cos(θ)⋅cos(ψ)     sin(θ)⋅sin(ϕ)
                        -sin(θ)⋅cos(ψ)            cos(θ)
``` 
"""
function 𝐞ˢ(
        ::Val{1},
        θ::T1 = zero(default_number_type()),
        ϕ::T2 = zero(θ),
        ψ::T3 = zero(ϕ);
        canonical = false,
    ) where {T1 <: Number, T2 <: Number, T3 <: Number}
    if canonical
        return Tens(
            Vec{3}(
                [
                    -sin(ψ) * sin(ϕ) + cos(θ) * cos(ψ) * cos(ϕ),
                    sin(ψ) * cos(ϕ) + sin(ϕ) * cos(θ) * cos(ψ),
                    -sin(θ) * cos(ψ),
                ]
            ),
        )
    else
        T = promote_type(T1, T2, T3)
        return Tens(Vec{3}([one(T), zero(T), zero(T)]), Basis(θ, ϕ, ψ))
    end
end
function 𝐞ˢ(
        ::Val{2},
        θ::T1 = zero(default_number_type()),
        ϕ::T2 = zero(θ),
        ψ::T3 = zero(ϕ);
        canonical = false,
    ) where {T1 <: Number, T2 <: Number, T3 <: Number}
    if canonical
        return Tens(
            Vec{3}(
                [
                    -sin(ψ) * cos(θ) * cos(ϕ) - sin(ϕ) * cos(ψ),
                    -sin(ψ) * sin(ϕ) * cos(θ) + cos(ψ) * cos(ϕ),
                    sin(θ) * sin(ψ),
                ]
            ),
        )
    else
        T = promote_type(T1, T2, T3)
        return Tens(Vec{3}([zero(T), one(T), zero(T)]), Basis(θ, ϕ, ψ))
    end
end
function 𝐞ˢ(
        ::Val{3},
        θ::T1 = zero(default_number_type()),
        ϕ::T2 = zero(θ),
        ψ::T3 = zero(ϕ);
        canonical = false,
    ) where {T1 <: Number, T2 <: Number, T3 <: Number}
    if canonical
        return Tens(Vec{3}([sin(θ) * cos(ϕ), sin(θ) * sin(ϕ), cos(θ)]))
    else
        T = promote_type(T1, T2, T3)
        return Tens(Vec{3}([zero(T), zero(T), one(T)]), Basis(θ, ϕ, ψ))
    end
end

for eb in (:𝐞, :𝐞ᵖ, :𝐞ᶜ, :𝐞ˢ)
    @eval $eb(i::Integer, args...; kwargs...) = $eb(Val(i), args...; kwargs...)
end
# The documented form of `𝐞`, with a plain dimension and element type: the
# forwarding above passes them on unchanged, where `Val`s are expected.
𝐞(i::Integer, dim::Integer, ::Type{T} = default_number_type()) where {T <: Number} = 𝐞(Val(i), Val(dim), Val(T))

"""
    init_spherical(coords = (symbols("θ ϕ", real = true)..., symbols("r", positive = true)); canonical = false)

Return the coordinates, base vectors and basis of the spherical basis.
Take care that the order of the 3 vectors is `𝐞ᶿ, 𝐞ᵠ, 𝐞ʳ` so that
the basis coincides with the canonical one when the angles are null and in consistency
the coordinates are ordered as `θ, ϕ, r`.

# Examples
```jldoctest
julia> coords, vectors, ℬˢ = init_spherical() ; θ, ϕ, r = coords ; 𝐞ᶿ, 𝐞ᵠ, 𝐞ʳ  = vectors ;
``` 
"""
init_spherical(coords; canonical = false) = Tuple(coords),
    ntuple(i -> 𝐞ˢ(Val(i), coords[1:2]...; canonical = canonical), 3),
    SphericalBasis(coords[1:2]...)

"""
    init_rotated(coords = symbols("θ ϕ ψ", real = true); canonical = false)

Return the angles, base vectors and basis of the rotated basis.
Here the coordinates are angles and do not represent a valid parametrization of `ℝ³`

# Examples
```jldoctest
julia> angles, vectors, ℬʳ = init_rotated() ; θ, ϕ, ψ = angles ; 𝐞ᶿ, 𝐞ᵠ, 𝐞ʳ = vectors ;
```
"""
init_rotated(angles; canonical = false) = Tuple(angles),
    ntuple(i -> 𝐞ˢ(Val(i), angles...; canonical = canonical), 3),
    Basis(angles...)

"""
    rotation_matrix(θ, ϕ, ψ = 0)

The 3×3 rotation matrix of the Euler angles `θ, ϕ, ψ` (ZYZ: `θ` the inclination
from ``\\mathbf{e}_3``, `ϕ` the azimuth, `ψ` the spin about the rotated axis).
Its third column is the unit vector of inclination `θ` and azimuth `ϕ`. With a
single angle, `rotation_matrix(θ)` is the 2D
rotation.

# Examples
```jldoctest
julia> cθ, cϕ, cψ, sθ, sϕ, sψ = symbols("cθ cϕ cψ sθ sϕ sψ", real = true) ;

julia> d = Dict(cos(θ) => cθ, cos(ϕ) => cϕ, cos(ψ) => cψ, sin(θ) => sθ, sin(ϕ) => sϕ, sin(ψ) => sψ) ;
ERROR: UndefVarError: `θ` not defined in `Main`
Suggestion: check for spelling errors or missing imports.
Stacktrace:
 [1] top-level scope
   @ none:1

julia> subs.(rotation_matrix(θ, ϕ, ψ),d...)
ERROR: UndefVarError: `θ` not defined in `Main`
Suggestion: check for spelling errors or missing imports.
Stacktrace:
 [1] top-level scope
   @ none:1
```
"""
rotation_matrix(θ, ϕ, ψ = 0) = RotZYZ(ϕ, θ, ψ)

"""
    rotation_matrix(θ)

The 2D rotation matrix of angle `θ`, as a `Tensor{2, 2}`.

# Examples
```jldoctest
julia> rotation_matrix(θ)
ERROR: UndefVarError: `θ` not defined in `Main`
Suggestion: check for spelling errors or missing imports.
Stacktrace:
 [1] top-level scope
   @ none:1
```
"""
rotation_matrix(θ) = Tensor{2, 2}((cos(θ), sin(θ), -sin(θ), cos(θ)))


"""
    rotation_tensor4(θ, ϕ = 0, ψ = 0)

The fourth-order tensor ``\\mathbf{R} \\boxtimes^s \\mathbf{R}`` that rotates a
symmetric second-order tensor, ``\\mathbf{R}\\,\\boldsymbol{\\varepsilon}\\,\\mathbf{R}^T``,
`R = rotation_matrix(θ, ϕ, ψ)`; its Kelvin-Mandel matrix is the 6×6 rotation
[`mandel_rotation`](@ref)`(θ, ϕ, ψ)`.

# Examples
```jldoctest
julia> cθ, cϕ, cψ, sθ, sϕ, sψ = symbols("cθ cϕ cψ sθ sϕ sψ", real = true) ;

julia> d = Dict(cos(θ) => cθ, cos(ϕ) => cϕ, cos(ψ) => cψ, sin(θ) => sθ, sin(ϕ) => sϕ, sin(ψ) => sψ) ;
ERROR: UndefVarError: `θ` not defined in `Main`
Suggestion: check for spelling errors or missing imports.
Stacktrace:
 [1] top-level scope
   @ none:1

julia> R = Tens(subs.(rotation_matrix(θ, ϕ, ψ),d...))
ERROR: UndefVarError: `θ` not defined in `Main`
Suggestion: check for spelling errors or missing imports.
Stacktrace:
 [1] top-level scope
   @ none:1

julia> RR = R ⊠ˢ R
ERROR: UndefVarError: `R` not defined in `Main`
Suggestion: add an appropriate import or assignment. This global was declared but not assigned.
Stacktrace:
 [1] top-level scope
   @ none:1

julia> R6 = from_mandel(subs.(mandel(rotation_tensor4(θ, ϕ, ψ)),d...))
ERROR: UndefVarError: `θ` not defined in `Main`
Suggestion: check for spelling errors or missing imports.
Stacktrace:
 [1] top-level scope
   @ none:1

julia> R6 == RR
ERROR: UndefVarError: `R6` not defined in `Main`
Suggestion: add an appropriate import or assignment. This global was declared but not assigned.
Stacktrace:
 [1] top-level scope
   @ none:1
```
"""
function rotation_tensor4(θ, ϕ = 0, ψ = 0)
    R = TensCanonical(rotation_matrix(θ, ϕ, ψ))
    return sboxtimes(R, R)
end

export levi_civita
export 𝐞, 𝐞ᵖ, 𝐞ᶜ, 𝐞ˢ
export init_cartesian, init_polar, init_cylindrical, init_spherical, init_rotated
export rotation_matrix, rotation_tensor4
