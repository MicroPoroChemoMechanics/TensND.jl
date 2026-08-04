"""
    SubManifoldSym(OM::AbstractTens{1,dim,Sym}, coords::NTuple{dim-1,Sym},
                   tmp_coords::NTuple = (), params::NTuple = ();
                   rules::Dict = Dict(), tmp_var::Dict = Dict(), to_coords::Dict = Dict())

Hypersurface of dimension `dim-1` embedded in `Rᵈⁱᵐ`, parametrized by the
position vector `OM` of the `dim-1` surface coordinates `coords`.

Unlike [`CoorSystemSym`](@ref), which describes a chart of the whole space, a
`SubManifoldSym` describes an embedded surface: it carries one coordinate fewer
than the ambient dimension, and it completes the tangent frame with the **unit
normal**, so that the stored bases are still `dim`-dimensional.

# Construction

The `dim-1` tangent vectors are the derivatives of the position vector,

    𝐚ᵢ = ∂ᵢOM ,  χᵢ = ‖𝐚ᵢ‖ ,  𝐞ᵢ = 𝐚ᵢ / χᵢ    (i = 1 … dim-1)

and the frame is closed by the unit normal `𝐧`, obtained as the generalized
cross product of the tangent vectors. `𝐧` is stored last, with Lamé coefficient
`χ_dim = 1`; `normal(SM)` returns it.

# Fundamental forms

Two order-2 tensors describe the surface, both stored with a vanishing
last row and column so that they live in the tangent block:

- the **first fundamental form** `𝐚` (induced metric), returned by
  [`submetric`](@ref) — `aᵢⱼ = 𝐚ᵢ ⋅ 𝐚ⱼ`;
- the **second fundamental form** `𝐛` (curvature tensor), returned by
  [`curvature`](@ref) — `bᵢⱼ = 𝐧 ⋅ ∂ᵢ𝐚ⱼ`.

The stored connection array `Γ` holds the Gauss–Weingarten equations of the
embedded frame rather than a plain Christoffel array: `Γ[:,:,dim]` is `𝐛`
(Gauss formula, the normal component of `∂ᵢ𝐚ⱼ`) and `Γ[:,dim,:]` is `-𝐛` with
one index raised (Weingarten formula, the tangential derivative of `𝐧`). The
purely tangential block, i.e. the intrinsic Christoffel symbols of the surface
metric, is what [`Riemann`](@ref) returns.

The optional `tmp_coords`, `params`, `rules`, `tmp_var` and `to_coords`
arguments have exactly the meaning they have for [`CoorSystemSym`](@ref): they
drive the symbolic simplification of intermediate expressions.

# Examples

```julia
julia> θ, ϕ = symbols("θ ϕ", real = true) ; R = symbols("R", positive = true) ;

julia> OM = Tens(R * [sin(θ)cos(ϕ), sin(θ)sin(ϕ), cos(θ)]) ;

julia> Sphere = SubManifoldSym(OM, (θ, ϕ), (), (R,)) ;

julia> normal(Sphere)          # outward unit normal 𝐞ʳ

julia> curvature(Sphere)       # 𝐛 = -𝐚/R for a sphere of radius R
```

See also [`CoorSystemSym`](@ref), [`normal`](@ref), [`submetric`](@ref),
[`curvature`](@ref), [`Riemann`](@ref).
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
        # Plain componentwise derivative of the position vector.
        #
        # This must NOT go through the two-argument `∂`: `@set_coorsys` adds a
        # method `∂(::AbstractTens, ::Sym)` bound to whatever system was made
        # the default, and that method is more specific than the variadic
        # fallback, so it wins. It then looks `coords[i]` up among *its own*
        # coordinates, fails to find it, and returns `zero(t)` — leaving every
        # tangent vector null, the frame matrix singular, and the constructor
        # throwing `NonInvertibleMatrixError`. Calling `∂` with an explicit
        # `Val` disambiguator would work too; going straight to `tdiff` is
        # simpler and makes the intent obvious.
        ∂coord(t, x) = change_tens(Tens(tdiff(components_canon(t), x)), get_basis(t), get_var(t))
        aᵢ = ntuple(i -> simp(chvar(∂coord(OMc, coords[i]), tmp_var)), dimm1)
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
        # Same reason as above: `∂coord`, not the two-argument `∂`.
        b = Tens(SymmetricTensor{2, dim, Sym}((i, j) -> i < dim && j < dim ? aᵢ[dim] ⋅ simp(chvar(∂coord(chvar(aᵢ[j], to_coords), coords[i]), tmp_var)) : zero(Sym)), natural_basis, (:cov, :cov))
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

"""
    normal(SM::SubManifoldSym{dim}) → AbstractTens{1,dim,Sym}

Unit normal `𝐧` of the hypersurface, i.e. the last vector of the natural basis.
Its orientation is that of the generalized cross product of the tangent vectors
`∂ᵢOM` taken in the order of `coords`, so reversing two coordinates reverses
`𝐧` — and with it the sign of [`curvature`](@ref).

# Examples
```julia
julia> θ, ϕ = symbols("θ ϕ", real = true) ; R = symbols("R", positive = true) ;

julia> Sphere = SubManifoldSym(Tens(R * [sin(θ)cos(ϕ), sin(θ)sin(ϕ), cos(θ)]), (θ, ϕ), (), (R,)) ;

julia> normal(Sphere)   # the outward radial vector 𝐞ʳ
```

See also [`SubManifoldSym`](@ref), [`curvature`](@ref).
"""
normal(SM::SubManifoldSym{dim}) where {dim} = natvec(SM, :cov)[dim]

"""
    submetric(SM::SubManifoldSym) → AbstractTens{2,dim,Sym}

First fundamental form `𝐚` of the hypersurface (the metric induced by the
embedding), `aᵢⱼ = 𝐚ᵢ ⋅ 𝐚ⱼ`.

It is returned as a `dim`-dimensional order-2 tensor whose last row and column
vanish, so that it can be combined directly with [`curvature`](@ref) and with
tensors expressed in the full embedded frame.

For a sphere of radius `R` parametrized by `(θ, ϕ)`,
`𝐚 = diag(R², R² sin²θ, 0)`.

See also [`SubManifoldSym`](@ref), [`curvature`](@ref).
"""
submetric(SM::SubManifoldSym) = SM.a

"""
    curvature(SM::SubManifoldSym) → AbstractTens{2,dim,Sym}

Second fundamental form `𝐛` of the hypersurface, `bᵢⱼ = 𝐧 ⋅ ∂ᵢ𝐚ⱼ`, with `𝐧`
the unit normal returned by [`normal`](@ref).

Like [`submetric`](@ref) it is stored as a `dim`-dimensional order-2 tensor with
a vanishing last row and column. **Its sign follows the orientation of `𝐧`**:
with the outward normal, a sphere of radius `R` gives `𝐛 = -𝐚/R`, hence the
principal curvatures `-1/R`.

The mixed form `𝐛` with one index raised is the Weingarten (shape) operator; its
trace is the mean curvature and its determinant, restricted to the tangent
block, the Gaussian curvature.

See also [`SubManifoldSym`](@ref), [`submetric`](@ref), [`normal`](@ref).
"""
curvature(SM::SubManifoldSym) = SM.b

"""
    Riemann(SM::SubManifoldSym{dim}) → Array{Sym,3}

Intrinsic connection coefficients of the hypersurface: the purely tangential
block `Γ[1:dim-1, 1:dim-1, 1:dim-1]` of the stored Gauss–Weingarten array, with
the convention `Γ[i,j,k] = Γᵏᵢⱼ`.

!!! warning "The name is misleading"
    Despite its name this function returns the **Christoffel symbols** of the
    induced metric, *not* the Riemann curvature tensor. For a sphere of radius
    `R` parametrized by `(θ, ϕ)` it returns `Γᶿ_ϕϕ = -sinθ cosθ` and
    `Γᵠ_θϕ = Γᵠ_ϕθ = cotθ`, which are connection coefficients — they are not
    tensor components and they vanish in a suitable chart, unlike a curvature.
    The intrinsic curvature is recovered from these by the Gauss equation, or
    directly from [`curvature`](@ref) and [`submetric`](@ref). The name is kept
    for backward compatibility.

See also [`SubManifoldSym`](@ref), [`curvature`](@ref), [`Christoffel`](@ref).
"""
Riemann(SM::SubManifoldSym{dim}) where {dim} = SM.Γ[1:(dim - 1), 1:(dim - 1), 1:(dim - 1)]

function ∂(
        t::AbstractTens{order, dim, T},
        i::Integer,
        SM::SubManifoldSym{dim},
    ) where {order, dim, T <: SymType}
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

∂(t::SymType, i::Integer, SM::SubManifoldSym{dim}) where {dim} =
    tdiff(only_coords(SM, t), getcoords(SM, i))

function ∂(
        t::AbstractTens{order, dim, T},
        x::SymType,
        SM::SubManifoldSym{dim},
    ) where {order, dim, T <: SymType}
    ind = findfirst(i -> i == x, getcoords(SM))
    return isnothing(ind) ? zero(t) : ∂(t, ind, SM)
end

function ∂(
        t::SymType,
        x::SymType,
        SM::SubManifoldSym{dim},
    ) where {dim}
    ind = findfirst(i -> i == x, getcoords(SM))
    return isnothing(ind) ? zero(t) : ∂(t, ind, SM)
end

"""
    GRAD(T::Union{Sym,AbstractTens{order,dim,Sym}},SM::SubManifoldSym{dim}) where {order,dim}

Calculate the gradient of `T` with respect to the coordinate system `SM`
"""
GRAD(
    t::Union{T, AbstractTens{order, dim, T}},
    SM::SubManifoldSym{dim},
) where {order, dim, T <: SymType} =
    sum([∂(t, i, SM) ⊗ natvec(SM, i, :cont) for i in 1:(dim - 1)])


"""
    SYMGRAD(T::Union{Sym,AbstractTens{order,dim,Sym}},SM::SubManifoldSym{dim}) where {order,dim}

Calculate the symmetrized gradient of `T` with respect to the coordinate system `SM`
"""
SYMGRAD(
    t::Union{T, AbstractTens{order, dim, T}},
    SM::SubManifoldSym{dim},
) where {order, dim, T <: SymType} =
    sum([∂(t, i, SM) ⊗ˢ natvec(SM, i, :cont) for i in 1:(dim - 1)])

"""
    DIV(T::AbstractTens{order,dim,Sym},SM::SubManifoldSym{dim}) where {order,dim}

Calculate the divergence  of `T` with respect to the coordinate system `SM`
"""
DIV(
    t::AbstractTens{order, dim, T},
    SM::SubManifoldSym{dim},
) where {order, dim, T <: SymType} =
    sum([∂(t, i, SM) ⋅ natvec(SM, i, :cont) for i in 1:(dim - 1)])

"""
    LAPLACE(T::Union{Sym,AbstractTens{order,dim,Sym}},SM::SubManifoldSym{dim}) where {order,dim}

Calculate the Laplace operator of `T` with respect to the coordinate system `SM`
"""
LAPLACE(
    t::Union{T, AbstractTens{order, dim, T}},
    SM::SubManifoldSym{dim},
) where {order, dim, T <: SymType} = DIV(GRAD(t, SM), SM)

"""
    HESS(T::Union{Sym,AbstractTens{order,dim,Sym}},SM::SubManifoldSym{dim}) where {order,dim}

Calculate the Hessian of `T` with respect to the coordinate system `SM`
"""
HESS(
    t::Union{T, AbstractTens{order, dim, T}},
    SM::SubManifoldSym{dim},
) where {order, dim, T <: SymType} = GRAD(GRAD(t, SM), SM)

export SubManifoldSym
export normal, submetric, curvature, Riemann
