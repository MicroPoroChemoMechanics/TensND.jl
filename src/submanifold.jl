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
struct SubManifoldSym{dim, VEC, BNORM, BNAT, TENSA, TENSB} <: AbstractCoorSystem{dim, Sym}
    OM::VEC
    coords::NTuple
    normalized_basis::BNORM
    natural_basis::BNAT
    aᵢ::NTuple{dim}
    χᵢ::NTuple{dim}
    aⁱ::NTuple{dim}
    eᵢ::NTuple{dim}
    a::TENSA
    b::TENSB
    Γ::Array{Sym, 3}
    tmp_coords::NTuple
    params::NTuple
    rules::Dict
    tmp_var::Dict
    to_coords::Dict
    function SubManifoldSym(
            OM::VEC,
            coords::NTuple{dimm1, Sym},
            tmp_coords::NTuple = (),
            params::NTuple = ();
            rules::Dict = Dict(),
            tmp_var::Dict = Dict(),
            to_coords::Dict = Dict(),
        ) where {VEC, dimm1}
        dim = dimm1 + 1
        simp(t) = length(rules) > 0 ? tsimplify(tsubs(tsimplify(t), rules...)) : tsimplify(t)
        chvar(t, d) = length(d) > 0 ? tsubs(t, d...) : t
        OMc = chvar(OM, to_coords)
        aᵢ = ntuple(i -> simp(chvar(∂(OMc, coords[i]), tmp_var)), dimm1)
        χᵢ = ntuple(i -> simp(norm(aᵢ[i])), dimm1)
        eᵢ = ntuple(i -> simp(aᵢ[i] / χᵢ[i]), dimm1)
        χᵢ = (ntuple(i -> simp(chvar(χᵢ[i], to_coords)), dimm1)..., one(Sym))
        eᵢ = ntuple(i -> simp(chvar(eᵢ[i], to_coords)), dimm1)
        A₀ = tsimplify(hcat(components_canon.(eᵢ)...))
        n = [tsimplify(det(hcat(A₀, [j == i ? one(Sym) : zero(Sym) for j in 1:dim]))) for i in 1:dim]
        n = n / tsimplify(norm(n))
        A = hcat(A₀, n)
        normalized_basis = Basis(A)
        eᵢ = ntuple(
            i -> Tens(
                Vec{dim}(j -> j == i ? one(Sym) : zero(Sym)),
                normalized_basis,
                (:cov,),
            ),
            dim,
        )
        aᵢ = ntuple(
            i -> Tens(Vec{dim}(j -> j == i ? χᵢ[i] : zero(Sym)), normalized_basis, (:cov,)),
            dim,
        )
        aⁱ = ntuple(
            i -> Tens(
                Vec{dim}(j -> j == i ? inv(χᵢ[i]) : zero(Sym)),
                normalized_basis,
                (:cont,),
            ),
            dim,
        )
        natural_basis = Basis(normalized_basis, χᵢ)
        a₀ = metric(natural_basis, :cov)
        a = Tens(SymmetricTensor{2, dim, Sym}((i, j) -> i < dim && j < dim ? a₀[i, j] : zero(Sym)), natural_basis, (:cov, :cov))
        b = Tens(SymmetricTensor{2, dim, Sym}((i, j) -> i < dim && j < dim ? aᵢ[dim] ⋅ simp(chvar(∂(chvar(aᵢ[j], to_coords), coords[i]), tmp_var)) : zero(Sym)), natural_basis, (:cov, :cov))
        Γ₀ = compute_Christoffel(
            coords,
            χᵢ,
            metric(normalized_basis, :cov),
            metric(normalized_basis, :cont),
        )
        Γ₁ = cat(Γ₀, b[1:(dim - 1), 1:(dim - 1)], dims = 3)
        bc = change_tens(b, (:cov, :cont))
        Γ = cat(Γ₁, reshape(-bc[1:(dim - 1), 1:dim], dim - 1, 1, dim), dims = 2)
        return new{dim, typeof(OM), typeof(normalized_basis), typeof(natural_basis), typeof(a), typeof(b)}(
            OMc,
            coords,
            normalized_basis,
            natural_basis,
            aᵢ,
            χᵢ,
            aⁱ,
            eᵢ,
            a,
            b,
            Γ,
            tmp_coords,
            params,
            rules,
            tmp_var,
            to_coords,
        )
    end
end

normal(SM::SubManifoldSym{dim}) where {dim} = natvec(SM, :cov)[dim]

submetric(SM::SubManifoldSym) = SM.a

curvature(SM::SubManifoldSym) = SM.b

Riemann(SM::SubManifoldSym{dim}) where {dim} = SM.Γ[1:(dim - 1), 1:(dim - 1), 1:(dim - 1)]

function ∂(
        t::AbstractTens{order, dim, Sym},
        i::Integer,
        SM::SubManifoldSym{dim},
    ) where {order, dim}
    t = only_coords(SM, t)
    ℬ = natural_basis(SM)
    var = ntuple(_ -> :cont, order)
    t = Array(components(t, ℬ, var))
    Γ = Christoffel(SM)
    data = tdiff(t, getcoords(SM, i))
    for o in 1:order
        ec1 = ntuple(j -> j == o ? order + 1 : j, order)
        ec2 = (order + 1, o)
        ec3 = ntuple(j -> j, order)
        data += einsum(EinCode((ec1, ec2), ec3), (t, view(Γ, i, :, :)))
    end
    return change_tens(Tens(simprules(data, SM), ℬ, var), normalized_basis(SM), var)
end

∂(t::Sym, i::Integer, SM::SubManifoldSym{dim}) where {dim} =
    tdiff(only_coords(SM, t), getcoords(SM, i))

function ∂(
        t::AbstractTens{order, dim, Sym},
        x::Sym,
        SM::SubManifoldSym{dim},
    ) where {order, dim}
    ind = findfirst(i -> i == x, getcoords(SM))
    return isnothing(ind) ? zero(t) : ∂(t, ind, SM)
end

function ∂(
        t::Sym,
        x::Sym,
        SM::SubManifoldSym{dim},
    ) where {dim}
    ind = findfirst(i -> i == x, getcoords(SM))
    return isnothing(ind) ? zero(t) : ∂(t, ind, SM)
end

"""
    GRAD(T::Union{Sym,AbstractTens{order,dim,Sym}},SM::SubManifoldSym{dim}) where {order,dim}

Calculate the gradient of `T` with respect to the coordinate system `SM`
"""
GRAD(T::Union{Sym, AbstractTens{order, dim, Sym}}, SM::SubManifoldSym{dim}) where {order, dim} =
    sum([∂(T, i, SM) ⊗ natvec(SM, i, :cont) for i in 1:(dim - 1)])


"""
    SYMGRAD(T::Union{Sym,AbstractTens{order,dim,Sym}},SM::SubManifoldSym{dim}) where {order,dim}

Calculate the symmetrized gradient of `T` with respect to the coordinate system `SM`
"""
SYMGRAD(
    T::Union{Sym, AbstractTens{order, dim, Sym}},
    SM::SubManifoldSym{dim},
) where {order, dim} = sum([∂(T, i, SM) ⊗ˢ natvec(SM, i, :cont) for i in 1:(dim - 1)])

"""
    DIV(T::AbstractTens{order,dim,Sym},SM::SubManifoldSym{dim}) where {order,dim}

Calculate the divergence  of `T` with respect to the coordinate system `SM`
"""
DIV(T::AbstractTens{order, dim, Sym}, SM::SubManifoldSym{dim}) where {order, dim} =
    sum([∂(T, i, SM) ⋅ natvec(SM, i, :cont) for i in 1:(dim - 1)])

"""
    LAPLACE(T::Union{Sym,AbstractTens{order,dim,Sym}},SM::SubManifoldSym{dim}) where {order,dim}

Calculate the Laplace operator of `T` with respect to the coordinate system `SM`
"""
LAPLACE(
    T::Union{Sym, AbstractTens{order, dim, Sym}},
    SM::SubManifoldSym{dim},
) where {order, dim} = DIV(GRAD(T, SM), SM)

"""
    HESS(T::Union{Sym,AbstractTens{order,dim,Sym}},SM::SubManifoldSym{dim}) where {order,dim}

Calculate the Hessian of `T` with respect to the coordinate system `SM`
"""
HESS(T::Union{Sym, AbstractTens{order, dim, Sym}}, SM::SubManifoldSym{dim}) where {order, dim} =
    GRAD(GRAD(T, SM), SM)

export SubManifoldSym
export normal, submetric, curvature, Riemann
