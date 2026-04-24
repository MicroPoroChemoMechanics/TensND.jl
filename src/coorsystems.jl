abstract type AbstractCoorSystem{dim, T <: Number} <: Any end


"""
    CoorSystemSym(OM::AbstractTens{1,dim,Sym},coords::NTuple{dim,Sym},bnorm::AbstractBasis{dim,Sym},χᵢ::NTuple{dim},
                  tmp_coords::NTuple = (),params::NTuple = ();rules::Dict = Dict(),tmp_var::Dict = Dict(),to_coords::Dict = Dict()) where {dim}
    CoorSystemSym(OM::AbstractTens{1,dim,Sym},coords::NTuple{dim,Sym},
                  tmp_coords::NTuple = (),params::NTuple = ();rules::Dict = Dict(),tmp_var::Dict = Dict(),to_coords::Dict = Dict()) where {dim}

Define a new coordinate system either from
1. the position vector `OM`, the coordinates `coords`, the basis of unit vectors (`𝐞ᵢ`) `bnorm` and the Lamé coefficients `χᵢ`

    In this case the natural basis is formed by the vectors `𝐚ᵢ = χᵢ 𝐞ᵢ` directly calculated from the input data.

1. or the position vector `OM` and the coordinates `coords`

    In this case the natural basis is formed by the vectors `𝐚ᵢ = ∂ᵢOM` i.e. by the derivative of the position vector with respect to the `iᵗʰ` coordinate

Optional parameters can be provided:
- `tmp_coords` contains temporary variables depending on coordinates (in order to allow symbolic simplifications)
- `params` contains possible parameters involved in `OM`
- `rules` contains a `Dict` with substitution rules to facilitate the simplification of formulas
- `tmp_var` contains a `Dict` with substitution of coordinates by temporary variables
- `to_coords` indicates how to eliminate the temporary variables to come back to the actual coordinates before derivation for Examples

# Examples
```julia
julia> ϕ, p = symbols("ϕ p", real = true) ;

julia> p̄, q, q̄, c = symbols("p̄ q q̄ c", positive = true) ;

julia> coords = (ϕ, p, q) ; tmp_coords = (p̄, q̄) ; params = (c,) ;

julia> OM = Tens(c * [p̄ * q̄ * cos(ϕ), p̄ * q̄ * sin(ϕ), p * q]) ;

julia> Spheroidal = CoorSystemSym(OM, coords, tmp_coords, params; tmp_var = Dict(1-p^2 => p̄^2, q^2-1 => q̄^2), to_coords = Dict(p̄ => √(1-p^2), q̄ => √(q^2-1))) ;
```
"""
# Build the canonical basis vectors (eᵢ), natural covariant (aᵢ) and
# contravariant (aⁱ) vectors from the normalized basis and Lamé coefficients.
function _build_basis_vectors(normalized_basis::AbstractBasis{dim, T}, χᵢ::NTuple{dim}) where {dim, T}
    eᵢ = ntuple(
        i -> Tens(Vec{dim}(j -> j == i ? one(T) : zero(T)), normalized_basis, (:cov,)),
        dim,
    )
    aᵢ = ntuple(
        i -> Tens(Vec{dim}(j -> j == i ? χᵢ[i] : zero(T)), normalized_basis, (:cov,)),
        dim,
    )
    aⁱ = ntuple(
        i -> Tens(Vec{dim}(j -> j == i ? inv(χᵢ[i]) : zero(T)), normalized_basis, (:cont,)),
        dim,
    )
    return eᵢ, aᵢ, aⁱ
end

struct CoorSystemSym{dim, T <: Number, VEC, BNORM, BNAT} <: AbstractCoorSystem{dim, T}
    OM::VEC
    coords::NTuple{dim, T}
    normalized_basis::BNORM
    natural_basis::BNAT
    aᵢ::NTuple{dim}
    χᵢ::NTuple{dim}
    aⁱ::NTuple{dim}
    eᵢ::NTuple{dim}
    Γ::Array{T, 3}
    tmp_coords::NTuple
    params::NTuple
    rules::Dict
    tmp_var::Dict
    to_coords::Dict
    function CoorSystemSym(
            OM::VEC,
            coords::NTuple{dim, T},
            normalized_basis::AbstractBasis{dim, T},
            χᵢ::NTuple{dim},
            tmp_coords::NTuple = (),
            params::NTuple = ();
            rules::Dict = Dict(),
            tmp_var::Dict = Dict(),
            to_coords::Dict = Dict()
        ) where {dim, T, VEC}
        simp(t) = length(rules) > 0 ? tsimplify(tsubs(tsimplify(t), rules...)) : tsimplify(t)
        eᵢ, aᵢ, aⁱ = _build_basis_vectors(normalized_basis, χᵢ)
        Γ = simp(
            compute_Christoffel(
                coords,
                χᵢ,
                metric(normalized_basis, :cov),
                metric(normalized_basis, :cont),
            )
        )
        natural_basis = Basis(normalized_basis, χᵢ)
        return new{dim, T, VEC, typeof(normalized_basis), typeof(natural_basis)}(
            OM,
            coords,
            normalized_basis,
            natural_basis,
            aᵢ,
            χᵢ,
            aⁱ,
            eᵢ,
            Γ,
            tmp_coords,
            params,
            rules,
            tmp_var,
            to_coords,
        )
    end
    function CoorSystemSym(
            OM::VEC,
            coords::NTuple{dim, T},
            tmp_coords::NTuple = (),
            params::NTuple = ();
            rules::Dict = Dict(),
            tmp_var::Dict = Dict(),
            to_coords::Dict = Dict()
        ) where {dim, T, VEC}
        simp(t) = length(rules) > 0 ? tsimplify(tsubs(t, rules...)) : tsimplify(t)
        chvar(t, d) = length(d) > 0 ? tsubs(t, d...) : t
        OMc = chvar(OM, to_coords)
        aᵢ = ntuple(i -> simp(chvar(∂(OMc, coords[i]), tmp_var)), dim)
        χᵢ = ntuple(i -> simp(norm(aᵢ[i])), dim)
        eᵢ = ntuple(i -> simp(aᵢ[i] / χᵢ[i]), dim)
        χᵢ = ntuple(i -> simp(chvar(χᵢ[i], to_coords)), dim)
        eᵢ = ntuple(i -> simp(chvar(eᵢ[i], to_coords)), dim)
        normalized_basis = Basis(tsimplify(hcat(components_canon.(eᵢ)...)))
        eᵢ, aᵢ, aⁱ = _build_basis_vectors(normalized_basis, χᵢ)
        Γ = compute_Christoffel(
            coords,
            χᵢ,
            metric(normalized_basis, :cov),
            metric(normalized_basis, :cont),
        )
        natural_basis = Basis(normalized_basis, χᵢ)
        return new{dim, T, VEC, typeof(normalized_basis), typeof(natural_basis)}(
            OMc,
            coords,
            normalized_basis,
            natural_basis,
            aᵢ,
            χᵢ,
            aⁱ,
            eᵢ,
            Γ,
            tmp_coords,
            params,
            rules,
            tmp_var,
            to_coords,
        )
    end
end

with_tmp_var(CS::CoorSystemSym, t) = length(CS.tmp_var) > 0 ? tsubs(t, CS.tmp_var...) : t
only_coords(CS::CoorSystemSym, t) = length(CS.to_coords) > 0 ? tsubs(t, CS.to_coords...) : t

getcoords(CS::CoorSystemSym) = CS.coords
getcoords(CS::CoorSystemSym, i::Integer) = getcoords(CS)[i]

@pure get_dim(::AbstractCoorSystem{dim}) where {dim} = dim

getOM(CS::AbstractCoorSystem) = CS.OM

normalized_basis(CS::AbstractCoorSystem) = CS.normalized_basis
natural_basis(CS::AbstractCoorSystem) = CS.natural_basis

Lame(CS::AbstractCoorSystem) = CS.χᵢ
Christoffel(CS::AbstractCoorSystem) = CS.Γ
# simprules(t, CS::AbstractCoorSystem) = length(CS.rules) > 0 ? tsimplify(tsubs(tsimplify(t), CS.rules...)) : tsimplify(t)
simprules(t, CS::AbstractCoorSystem) = length(CS.rules) > 0 ? tsubs(t, CS.rules...) : t

natvec(CS::AbstractCoorSystem, ::Val{:cov}) = CS.aᵢ
natvec(CS::AbstractCoorSystem, ::Val{:cont}) = CS.aⁱ
natvec(CS::AbstractCoorSystem, var = :cov) = natvec(CS, Val(var))
natvec(CS::AbstractCoorSystem, i::Integer, var = :cov) = natvec(CS, var)[i]

unitvec(CS::AbstractCoorSystem) = CS.eᵢ
unitvec(CS::AbstractCoorSystem, i::Integer) = unitvec(CS)[i]

function compute_Christoffel(coords, χ, γ, invγ)
    dim = length(coords)
    gᵢⱼ = [γ[i, j] * χ[i] * χ[j] for i in 1:dim, j in 1:dim]
    gⁱʲ = [invγ[i, j] / (χ[i] * χ[j]) for i in 1:dim, j in 1:dim]
    ∂g = [tdiff(gᵢⱼ[i, j], coords[k]) for i in 1:dim, j in 1:dim, k in 1:dim]
    Γᵢⱼₖ =
        [(∂g[i, k, j] + ∂g[j, k, i] - ∂g[i, j, k]) / 2 for i in 1:dim, j in 1:dim, k in 1:dim]
    return ein"ijl,lk->ijk"(Γᵢⱼₖ, gⁱʲ)
end

"""
    ∂(t::AbstractTens{order,dim,T,A},xᵢ::T)

Return the derivative of the tensor `t` with respect to the variable `x_i`

# Examples
```julia

julia> (θ, ϕ, r), (𝐞ᶿ, 𝐞ᵠ, 𝐞ʳ), ℬˢ = init_spherical() ;

julia> ∂(𝐞ʳ, ϕ) == sin(θ) * 𝐞ᵠ
true

julia> ∂(𝐞ʳ ⊗ 𝐞ʳ,θ)
Tens.TensRotated{2, 3, Sym, SymmetricTensor{2, 3, Sym, 6}}
# data: 3×3 SymmetricTensor{2, 3, Sym, 6}:
 0  0  1
 0  0  0
 1  0  0
# basis: 3×3 Tensor{2, 3, Sym, 9}:
 cos(θ)⋅cos(ϕ)  -sin(ϕ)  sin(θ)⋅cos(ϕ)
 sin(ϕ)⋅cos(θ)   cos(ϕ)  sin(θ)⋅sin(ϕ)
       -sin(θ)        0         cos(θ)
# var: (:cont, :cont)
```
"""
∂(t::AbstractTens{order, dim, T}, xᵢ...) where {order, dim, T <: SymType} =
    change_tens(Tens(tdiff(components_canon(t), xᵢ...)), get_basis(t), get_var(t))

∂(t::SymType, xᵢ...) = tdiff(t, xᵢ...)

function ∂(
        t::AbstractTens{order, dim, T},
        i::Integer,
        CS::CoorSystemSym{dim},
    ) where {order, dim, T <: SymType}
    t = only_coords(CS, t)
    ℬ = natural_basis(CS)
    var = ntuple(_ -> :cont, order)
    t = Array(components(t, ℬ, var))
    Γ = Christoffel(CS)
    data = tdiff(t, getcoords(CS, i))
    for o in 1:order
        ec1 = ntuple(j -> j == o ? order + 1 : j, order)
        ec2 = (order + 1, o)
        ec3 = ntuple(j -> j, order)
        data += einsum(EinCode((ec1, ec2), ec3), (t, view(Γ, i, :, :)))
    end
    return change_tens(Tens(simprules(data, CS), ℬ, var), normalized_basis(CS), var)
end

∂(t::T, i::Integer, CS::CoorSystemSym{dim, T}) where {dim, T <: SymType} =
    tdiff(only_coords(CS, t), getcoords(CS, i))

function ∂(
        t::AbstractTens{order, dim, T},
        x::T,
        CS::CoorSystemSym{dim, T},
    ) where {order, dim, T <: SymType}
    ind = findfirst(i -> i == x, getcoords(CS))
    return ind === nothing ? zero(t) : ∂(t, ind, CS)
end

function ∂(
        t::T,
        x::T,
        CS::CoorSystemSym{dim, T},
    ) where {dim, T <: SymType}
    ind = findfirst(i -> i == x, getcoords(CS))
    return ind === nothing ? zero(t) : ∂(t, ind, CS)
end

"""
    GRAD(t::Union{t,AbstractTens{order,dim,T}},CS::CoorSystemSym{dim}) where {order,dim,T<:Number}

Calculate the gradient of `t` with respect to the coordinate system `CS`
"""
GRAD(t::Union{T, AbstractTens{order, dim, T}}, CS::CoorSystemSym{dim, T}) where {order, dim, T <: SymType} =
    sum([∂(t, i, CS) ⊗ natvec(CS, i, :cont) for i in 1:dim])


"""
    SYMGRAD(t::Union{T,AbstractTens{order,dim,T}},CS::CoorSystemSym{dim}) where {order,dim,T<:Number}

Calculate the symmetrized gradient of `T` with respect to the coordinate system `CS`
"""
SYMGRAD(
    t::Union{T, AbstractTens{order, dim, T}},
    CS::CoorSystemSym{dim, T},
) where {order, dim, T <: SymType} = sum([∂(t, i, CS) ⊗ˢ natvec(CS, i, :cont) for i in 1:dim])


"""
    DIV(t::AbstractTens{order,dim,Sym},CS::CoorSystemSym{dim}) where {order,dim,T<:Number}

Calculate the divergence  of `T` with respect to the coordinate system `CS`
"""
DIV(t::AbstractTens{order, dim, T}, CS::CoorSystemSym{dim, T}) where {order, dim, T <: SymType} =
    sum([∂(t, i, CS) ⋅ natvec(CS, i, :cont) for i in 1:dim])


"""
    LAPLACE(t::Union{Sym,AbstractTens{order,dim,Sym}},CS::CoorSystemSym{dim}) where {order,dim,T<:Number}

Calculate the Laplace operator of `T` with respect to the coordinate system `CS`
"""
LAPLACE(
    t::Union{T, AbstractTens{order, dim, T}},
    CS::CoorSystemSym{dim, T},
) where {order, dim, T <: SymType} = DIV(GRAD(t, CS), CS)

"""
    HESS(t::Union{Sym,AbstractTens{order,dim,Sym}},CS::CoorSystemSym{dim}) where {order,dim,T<:Number}

Calculate the Hessian of `T` with respect to the coordinate system `CS`
"""
HESS(t::Union{T, AbstractTens{order, dim, T}}, CS::CoorSystemSym{dim, T}) where {order, dim, T <: SymType} =
    GRAD(GRAD(t, CS), CS)

"""
    coorsys_cartesian(coords = symbols("x y z", real = true))

Return the cartesian coordinate system

# Examples
```julia
julia> Cartesian = coorsys_cartesian() ; 𝐗 = getcoords(Cartesian) ; 𝐄 = unitvec(Cartesian) ; ℬ = get_basis(Cartesian)

julia> 𝛔 = Tens(SymmetricTensor{2,3}((i, j) -> SymFunction("σ\$i\$j", real = true)(𝐗...))) ;

julia> DIV(𝛔, CScar)
Tens.TensCanonical{1, 3, Sym, Vec{3, Sym}}
# data: 3-element Vec{3, Sym}:
 Derivative(σ11(x, y, z), x) + Derivative(σ21(x, y, z), y) + Derivative(σ31(x, y, z), z)
 Derivative(σ21(x, y, z), x) + Derivative(σ22(x, y, z), y) + Derivative(σ32(x, y, z), z)
 Derivative(σ31(x, y, z), x) + Derivative(σ32(x, y, z), y) + Derivative(σ33(x, y, z), z)
# basis: 3×3 Tens.LazyIdentity{3, Sym}:
 1  0  0
 0  1  0
 0  0  1
# var: (:cont,)
``` 
"""
function coorsys_cartesian(coords::NTuple{dim, T} = symbols("x y z", real = true)) where {dim, T <: SymType}
    𝐗, 𝐄, ℬ = init_cartesian(coords)
    OM = sum([𝐗[i] * 𝐄[i] for i in 1:dim])
    χᵢ = ntuple(_ -> one(eltype(coords)), dim)
    return CoorSystemSym(OM, coords, ℬ, χᵢ)
end
coorsys_cartesian(::Val{<:Sym}, coords = symbols("x y z", real = true)) = coorsys_cartesian(coords)
coorsys_cartesian(::Val{Num}, coords = Tuple(@variables x y z)) = coorsys_cartesian(coords)


"""
    coorsys_polar(coords = (symbols("r", positive = true), symbols("θ", real = true)); canonical = false)

Return the polar coordinate system

# Examples
```julia
julia> Polar = coorsys_polar() ; r, θ = getcoords(Polar) ; 𝐞ʳ, 𝐞ᶿ = unitvec(Polar) ; ℬᵖ = get_basis(Polar)

julia> f = SymFunction("f", real = true)(r, θ) ;

julia> LAPLACE(f, Polar)
                               2
                              ∂
                             ───(f(r, θ))
                               2
               ∂             ∂θ
  2            ──(f(r, θ)) + ────────────
 ∂             ∂r                 r
───(f(r, θ)) + ──────────────────────────
  2                        r
∂r
``` 
"""
function coorsys_polar(
        coords::NTuple{2, T} = (symbols("r", positive = true), symbols("θ", real = true));
        canonical = false
    ) where {T <: SymType}
    (r, θ), (𝐞ʳ, 𝐞ᶿ), ℬᵖ = init_polar(coords, canonical = canonical)
    OM = r * 𝐞ʳ
    return CoorSystemSym(OM, coords, ℬᵖ, (one(eltype(coords)), r))
end
coorsys_polar(
    ::Val{<:Sym}, coords = (
        symbols("r", positive = true),
        symbols("θ", real = true),
    ); canonical = false
) = coorsys_polar(coords; canonical = canonical)
coorsys_polar(::Val{Num}, coords = Tuple(@variables r θ); canonical = false) = coorsys_polar(coords; canonical = canonical)


"""
    coorsys_cylindrical(coords = (symbols("r", positive = true), symbols("θ", real = true), symbols("z", real = true)); canonical = false)

Return the cylindrical coordinate system

# Examples
```julia
julia> Cylindrical = coorsys_cylindrical() ; rθz = getcoords(Cylindrical) ; 𝐞ʳ, 𝐞ᶿ, 𝐞ᶻ = unitvec(Cylindrical) ; ℬᶜ = get_basis(Cylindrical)

julia> 𝐯 = Tens(Vec{3}(i -> SymFunction("v\$(rθz[i])", real = true)(rθz...)), ℬᶜ) ;

julia> DIV(𝐯, Cylindrical)
                                                  ∂
                                    vr(r, θ, z) + ──(vθ(r, θ, z))
∂                 ∂                               ∂θ
──(vr(r, θ, z)) + ──(vz(r, θ, z)) + ─────────────────────────────
∂r                ∂z                              r
``` 
"""
function coorsys_cylindrical(
        coords::NTuple{3, T} = (
            symbols("r", positive = true),
            symbols("θ", real = true),
            symbols("z", real = true),
        );
        canonical = false
    ) where {T <: SymType}
    (r, θ, z), (𝐞ʳ, 𝐞ᶿ, 𝐞ᶻ), ℬᶜ = init_cylindrical(coords, canonical = canonical)
    OM = r * 𝐞ʳ + z * 𝐞ᶻ

    return CoorSystemSym(OM, coords, ℬᶜ, (one(eltype(coords)), r, one(eltype(coords))))
end
coorsys_cylindrical(
    ::Val{<:Sym}, coords = (
        symbols("r", positive = true),
        symbols("θ", real = true), symbols("z", real = true),
    ); canonical = false
) = coorsys_cylindrical(coords; canonical = canonical)
coorsys_cylindrical(::Val{Num}, coords = Tuple(@variables r θ z); canonical = false) = coorsys_cylindrical(coords; canonical = canonical)


"""
    coorsys_spherical(coords = (symbols("θ", real = true), symbols("ϕ", real = true), symbols("r", positive = true)); canonical = false)

Return the spherical coordinate system

# Examples
```julia
julia> Spherical = coorsys_spherical() ; θ, ϕ, r = getcoords(Spherical) ; 𝐞ᶿ, 𝐞ᵠ, 𝐞ʳ = unitvec(Spherical) ; ℬˢ = get_basis(Spherical)

julia> for σⁱʲ ∈ ("σʳʳ", "σᶿᶿ", "σᵠᵠ") @eval \$(Symbol(σⁱʲ)) = SymFunction(\$σⁱʲ, real = true)(\$r) end ;

julia> 𝛔 = σʳʳ * 𝐞ʳ ⊗ 𝐞ʳ + σᶿᶿ * 𝐞ᶿ ⊗ 𝐞ᶿ + σᵠᵠ * 𝐞ᵠ ⊗ 𝐞ᵠ ;

julia> div𝛔 = DIV(𝛔, Spherical)
Tens.TensRotated{1, 3, Sym, Vec{3, Sym}}
# data: 3-element Vec{3, Sym}:
                              (-σᵠᵠ(r) + σᶿᶿ(r))*cos(θ)/(r*sin(θ))
                                                                 0
 Derivative(σʳʳ(r), r) + (σʳʳ(r) - σᵠᵠ(r))/r + (σʳʳ(r) - σᶿᶿ(r))/r
# basis: 3×3 Tensor{2, 3, Sym, 9}:
 cos(θ)⋅cos(ϕ)  -sin(ϕ)  sin(θ)⋅cos(ϕ)
 sin(ϕ)⋅cos(θ)   cos(ϕ)  sin(θ)⋅sin(ϕ)
       -sin(θ)        0         cos(θ)
# var: (:cont,)

julia> div𝛔 ⋅ 𝐞ʳ
d            σʳʳ(r) - σᵠᵠ(r)   σʳʳ(r) - σᶿᶿ(r)
──(σʳʳ(r)) + ─────────────── + ───────────────
dr                  r                 r
``` 
"""
function coorsys_spherical(
        coords::NTuple{3, T} = (
            symbols("θ", real = true),
            symbols("ϕ", real = true),
            symbols("r", positive = true),
        );
        canonical = false
    ) where {T <: SymType}
    (θ, ϕ, r), (𝐞ᶿ, 𝐞ᵠ, 𝐞ʳ), ℬˢ = init_spherical(coords, canonical = canonical)
    OM = r * 𝐞ʳ
    rules = Dict(abs(sin(θ)) => sin(θ), transpose(tan(θ)) => tan(θ), 1 // 1 => 1, 2 // 1 => 2)
    return CoorSystemSym(OM, coords, ℬˢ, (r, r * sin(θ), one(eltype(coords))); rules = rules)
end
coorsys_spherical(
    ::Val{<:Sym}, coords = (
        symbols("θ", real = true),
        symbols("ϕ", real = true), symbols("r", positive = true),
    ); canonical = false
) = coorsys_spherical(coords; canonical = canonical)
coorsys_spherical(::Val{Num}, coords = Tuple(@variables θ ϕ r); canonical = false) = coorsys_spherical(coords; canonical = canonical)

"""
    coorsys_spheroidal(coords = (symbols("ϕ", real = true),symbols("p", real = true),symbols("q", positive = true),),
                            c = symbols("c", positive = true),tmp_coords = (symbols("p̄ q̄", positive = true)...,),)

Return the spheroidal coordinate system

# Examples
```julia
julia> Spheroidal = coorsys_spheroidal() ; OM = getOM(Spheroidal)
Tens.TensCanonical{1, 3, Sym, Vec{3, Sym}}
# data: 3-element Vec{3, Sym}:
 c⋅p̄⋅q̄⋅cos(ϕ)
 c⋅p̄⋅q̄⋅sin(ϕ)
          c⋅p⋅q
# basis: 3×3 Tens.LazyIdentity{3, Sym}:
 1  0  0
 0  1  0
 0  0  1
# var: (:cont,)

julia> LAPLACE(OM[1]^2, Spheroidal)
2
``` 
"""
function coorsys_spheroidal(
        coords::NTuple{3, T} = (
            symbols("ϕ", real = true),
            symbols("p", real = true),
            symbols("q", positive = true),
        ),
        c = symbols("c", positive = true),
        tmp_coords = symbols("p̄ q̄", positive = true),
    ) where {T <: SymType}
    ϕ, p, q = coords
    params = (c,)
    p̄, q̄ = tmp_coords
    OM = Tens(c * [p̄ * q̄ * cos(ϕ), p̄ * q̄ * sin(ϕ), p * q])
    # OM = Tens(c * [√(1 - p^2) * √(q^2 - 1) * cos(ϕ), √(1 - p^2) * √(q^2 - 1) * sin(ϕ), p * q])
    ℬ = RotatedBasis(
        T[
            -sin(ϕ) -p * sqrt(q^2 - 1) * cos(ϕ) / sqrt(q^2 - p^2) q * sqrt(1 - p^2) * cos(ϕ) / sqrt(q^2 - p^2)
            cos(ϕ) -p * sqrt(q^2 - 1) * sin(ϕ) / sqrt(q^2 - p^2) q * sqrt(1 - p^2) * sin(ϕ) / sqrt(q^2 - p^2)
            0 q * sqrt(1 - p^2) / sqrt(q^2 - p^2) p * sqrt(q^2 - 1) / sqrt(q^2 - p^2)
        ],
    )
    χᵢ = (
        c * sqrt(1 - p^2) * sqrt(q^2 - 1),
        c * sqrt(q^2 - p^2) / sqrt(1 - p^2),
        c * sqrt(q^2 - p^2) / sqrt(q^2 - 1),
    )
    tmp_var = Dict(1 - p^2 => p̄^2, q^2 - 1 => q̄^2)
    to_coords = Dict(p̄ => √(1 - p^2), q̄ => √(q^2 - 1))
    return CoorSystemSym(
        OM,
        coords,
        ℬ,
        χᵢ,
        tmp_coords,
        params;
        tmp_var = tmp_var,
        to_coords = to_coords
    )
end
coorsys_spheroidal(
    ::Val{<:Sym}, coords = (
        symbols("ϕ", real = true),
        symbols("p", real = true),
        symbols("q", positive = true),
    ),
    c = symbols("c", positive = true),
    tmp_coords = symbols("p̄ q̄", positive = true)
) = coorsys_spheroidal(coords, c, tmp_coords)
coorsys_spheroidal(
    ::Val{Num}, coords = Tuple(@variables ϕ p q),
    c = (@variables c)[1], tmp_coords = Tuple(@variables p̄ q̄)
) = coorsys_spheroidal(coords, c, tmp_coords)


for cs in (:coorsys_cartesian, :coorsys_polar, :coorsys_cylindrical, :coorsys_spherical, :coorsys_spheroidal)
    @eval $cs(T::Type, args...; kwargs...) = $cs(Val(T), args...; kwargs...)
end


"""
    @set_coorsys CS
    @set_coorsys(CS)

Set a coordinate system in order to avoid precising it in differential operators

# Examples
```julia
julia> Spherical = coorsys_spherical() ; θ, ϕ, r = getcoords(Spherical) ; 𝐞ᶿ, 𝐞ᵠ, 𝐞ʳ = unitvec(Spherical) ; vec = ("𝐞ᶿ", "𝐞ᵠ", "𝐞ʳ") ;

julia> @set_coorsys Spherical

julia> intrinsic(GRAD(𝐞ʳ),vec)
(1/r)𝐞ᶿ⊗𝐞ᶿ + (1/r)𝐞ᵠ⊗𝐞ᵠ

julia> intrinsic(DIV(𝐞ʳ ⊗ 𝐞ʳ),vec)
(2/r)𝐞ʳ

julia> LAPLACE(1/r)
0
``` 
"""
macro set_coorsys(CS = coorsys_cartesian(), vec = '𝐞', coords = nothing)
    m = @__MODULE__
    return quote
        $m.∂(t::AbstractTens{order, dim, T}, i::Integer) where {order, dim, T <: SymType} = $m.∂(t, i, $(esc(CS)))
        $m.∂(t::AbstractTens{order, dim, T}, x::Sym) where {order, dim, T <: SymType} = $m.∂(t, x, $(esc(CS)))
        $m.∂(t::SymType, i::Integer) = $m.∂(t, i, $(esc(CS)))
        $m.∂(t::SymType, x::Sym) = $m.∂(t, x, $(esc(CS)))
        $m.GRAD(t::Union{T, AbstractTens}) where {T <: SymType} = $m.GRAD(t, $(esc(CS)))
        $m.SYMGRAD(t::Union{T, AbstractTens}) where {T <: SymType} = $m.SYMGRAD(t, $(esc(CS)))
        $m.DIV(t::AbstractTens) = $m.DIV(t, $(esc(CS)))
        $m.LAPLACE(t::Union{T, AbstractTens}) where {T <: SymType} = $m.LAPLACE(t, $(esc(CS)))
        $m.HESS(t::Union{T, AbstractTens}) where {T <: SymType} = $m.HESS(t, $(esc(CS)))

        if $(esc(coords)) === nothing
            coords = string.(getcoords($(esc(CS))))
        end
        dim = get_dim($(esc(CS)))
        if length(coords) == dim - 1
            coords = (coords..., dim)
        end
        ℬ = normalized_basis($(esc(CS)))
        $m.intrinsic(t::AbstractTens{order, dim, T}) where {order, dim, T <: SymType} = intrinsic(change_tens(t, ℬ); vec = $(esc(vec)), coords = coords)

        # Base.show(t::AbstractTens{order,dim,T}) where {order,dim,T<:SymType} = intrinsic(change_tens(t, ℬ); vec = $(esc(vec)), coords = coords)
        # Base.print(t::AbstractTens{order,dim,T}) where {order,dim,T<:SymType} = intrinsic(change_tens(t, ℬ); vec = $(esc(vec)), coords = coords)
        # Base.display(t::AbstractTens{order,dim,T}) where {order,dim,T<:SymType} = intrinsic(change_tens(t, ℬ); vec = $(esc(vec)), coords = coords)

    end
end

function intrinsic(t::AbstractTens{order, dim, T}, CS::AbstractCoorSystem; vec = '𝐞') where {order, dim, T}
    coords = string.(getcoords(CS))
    ℬ = normalized_basis(CS)
    return intrinsic(change_tens(t, ℬ); vec = vec, coords = coords)
end

export ∂, CoorSystemSym, Lame, Christoffel
export GRAD, SYMGRAD, DIV, LAPLACE, HESS
export normalized_basis, natural_basis, natvec, unitvec, getcoords, getOM
export coorsys_cartesian, coorsys_polar, coorsys_cylindrical, coorsys_spherical, coorsys_spheroidal
export @set_coorsys
