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

Surface indices run from `1` to `dim-1` and are written `α, β, γ`; ambient
indices run to `dim` and are written `i, j, k`.

The `dim-1` tangent vectors are the derivatives of the position vector,

    𝐚_α = ∂_α OM ,  χ_α = ‖𝐚_α‖ ,  𝐞_α = 𝐚_α / χ_α    (α = 1 … dim-1)

and the frame is closed by the unit normal `𝐧`, obtained as the generalized
cross product of the tangent vectors. `𝐧` is stored last, with Lamé coefficient
`χ_dim = 1`; `normal(SM)` returns it.

# Fundamental forms

Two order-2 tensors describe the surface, both stored with a vanishing
last row and column so that they live in the tangent block:

- the **first fundamental form** `𝐚` (induced metric), returned by
  [`submetric`](@ref) — `a_αβ = 𝐚_α ⋅ 𝐚_β`;
- the **second fundamental form** `𝐛` (curvature tensor), returned by
  [`curvature`](@ref) — `b_αβ = 𝐧 ⋅ ∂_α 𝐚_β`.

The stored array `Γ` holds the Gauss–Weingarten equations of the embedded frame
rather than a plain Christoffel array: `Γ[:,:,dim]` is `𝐛` (Gauss formula, the
normal component of `∂_α 𝐚_β`) and `Γ[:,dim,:]` is `-𝐛` with one index raised
(Weingarten formula, the tangential derivative of `𝐧`). Its purely tangential
block — the **intrinsic connection coefficients** `Γᵞ_αβ` of the induced
metric — is what [`connection`](@ref) returns.

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
[`curvature`](@ref), [`connection`](@ref).
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
        # Shared with `CoorSystemSym`: builds 𝐞ᵢ, 𝐚ᵢ and 𝐚ⁱ with the *correct*
        # variances. Duplicating it here is what let the covariant/contravariant
        # swap survive in the submanifold path after it was fixed for charts.
        eᵢ, aᵢ, aⁱ = _build_basis_vectors(normalized_basis, χᵢ)
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
`∂_α OM` taken in the order of `coords`, so reversing two coordinates reverses
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
embedding), `a_αβ = 𝐚_α ⋅ 𝐚_β`, with Greek surface indices running to `dim-1`.

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

Second fundamental form `𝐛` of the hypersurface, `b_αβ = 𝐧 ⋅ ∂_α 𝐚_β`, with `𝐧`
the unit normal returned by [`normal`](@ref).

Like [`submetric`](@ref) it is stored as a `dim`-dimensional order-2 tensor with
a vanishing last row and column. **Its sign follows the orientation of `𝐧`**:
with the outward normal, a sphere of radius `R` gives `𝐛 = -𝐚/R`, hence the
principal curvatures `-1/R`.

The mixed form `b_α{}^β` is the Weingarten (shape) operator; its
trace is the mean curvature and its determinant, restricted to the tangent
block, the Gaussian curvature.

See also [`SubManifoldSym`](@ref), [`submetric`](@ref), [`normal`](@ref).
"""
curvature(SM::SubManifoldSym) = SM.b

"""
    connection(SM::SubManifoldSym{dim}) → Array{Sym,3}

**Intrinsic connection coefficients** of the hypersurface — the Christoffel
symbols `Γᵞ_αβ` of the induced metric — indexed `Γ[α,β,γ]`, the contravariant
index last.

Surface indices run from `1` to `dim-1` and are written with Greek letters
`α, β, γ` throughout, to distinguish them from ambient indices `i, j, k` which
run to `dim`. This is the purely tangential block of the Gauss–Weingarten array
returned by [`Christoffel`](@ref).

For a sphere of radius `R` parametrized by `(θ, φ)`:

    Γᶿ_φφ = −sinθ cosθ ,   Γᵠ_θφ = Γᵠ_φθ = cotθ

These are **connection coefficients, not curvature**: they are not the
components of a tensor, they transform inhomogeneously under a change of chart,
and they can be made to vanish at any single point. The intrinsic curvature
follows from them by the Gauss equation, or directly from [`curvature`](@ref)
and [`submetric`](@ref).

!!! note "Renamed from `Riemann`"
    This function was called `Riemann`, which named the wrong object entirely —
    it never returned a Riemann curvature tensor. `Riemann` still works and
    forwards here, with a deprecation warning.

See also [`SubManifoldSym`](@ref), [`Christoffel`](@ref), [`curvature`](@ref).
"""
connection(SM::SubManifoldSym{dim}) where {dim} = SM.Γ[1:(dim - 1), 1:(dim - 1), 1:(dim - 1)]

@deprecate Riemann(SM::SubManifoldSym) connection(SM)

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
export normal, submetric, curvature, connection, Riemann
