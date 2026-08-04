abstract type AbstractCoorSystem{dim, T <: Number} <: Any end


# Build the canonical basis vectors (eᵢ), natural covariant (aᵢ) and
# contravariant (aⁱ) vectors from the normalized basis and Lamé coefficients.
function _build_basis_vectors(normalized_basis::AbstractBasis{dim, T}, χᵢ::NTuple{dim}) where {dim, T}
    # Variances matter here, and only on a NON-orthogonal chart: components
    # `(0,…,1,…,0)` with variance `:cont` expand on the basis vectors 𝐞ᵢ, with
    # `:cov` on the *dual* vectors 𝐞ⁱ. The two coincide when the normalized
    # basis is orthonormal — which every predefined system is — so a swap here
    # is invisible until someone builds a skew chart, and then it silently
    # returns the dual frame and every differential operator is wrong.
    eᵢ = ntuple(
        i -> Tens(Vec{dim}(j -> j == i ? one(T) : zero(T)), normalized_basis, (:cont,)),
        dim,
    )
    aᵢ = ntuple(
        i -> Tens(Vec{dim}(j -> j == i ? χᵢ[i] : zero(T)), normalized_basis, (:cont,)),
        dim,
    )
    aⁱ = ntuple(
        i -> Tens(Vec{dim}(j -> j == i ? inv(χᵢ[i]) : zero(T)), normalized_basis, (:cov,)),
        dim,
    )
    return eᵢ, aᵢ, aⁱ
end

"""
    ChartCore{dim,T,VEC,BNORM,BNAT}

The geometry shared by every symbolic coordinate system: the position vector,
the coordinates, the two bases, the natural and dual frames, the Lamé
coefficients, the connection array, and the symbolic-simplification settings.

Held by both [`CoorSystemSym`](@ref) and [`SubManifoldSym`](@ref) — which
differ only in what they *add* to it (a submanifold adds the two fundamental
forms) and in how many directions they differentiate along ([`nderiv`](@ref)).
Declaring the block once is what keeps the two in step: it was previously
duplicated field for field, and every fix to one had to be repeated in the
other.

Reached through [`core`](@ref); user code should go through the accessors
([`getcoords`](@ref), [`Lame`](@ref), [`Christoffel`](@ref), …) instead.
"""
struct ChartCore{dim, T, VEC, BNORM, BNAT}
    OM::VEC
    coords::NTuple
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
end

"""
    core(CS::AbstractCoorSystem) → ChartCore

The shared geometry block of `CS`. Internal: prefer the named accessors.
"""
core(CS::AbstractCoorSystem) = CS.core

"""
    CoorSystemSym{dim,T,VEC,BNORM,BNAT}

Curvilinear coordinate system with **exact symbolic** derivatives.

Stores the position vector, the coordinate symbols, the natural and normalized
bases, the Lamé coefficients and the Christoffel symbols as expressions, so the
differential operators [`GRAD`](@ref), [`SYMGRAD`](@ref), [`DIV`](@ref),
[`LAPLACE`](@ref) and [`HESS`](@ref) return closed-form results valid
everywhere, not values at a point.

Predefined systems: [`coorsys_cartesian`](@ref), [`coorsys_polar`](@ref),
[`coorsys_cylindrical`](@ref), [`coorsys_spherical`](@ref),
[`coorsys_spheroidal`](@ref). [`@set_coorsys`](@ref) makes one the default so
the operators take a single argument.

For pointwise evaluation by automatic differentiation, use
[`CoorSystemNum`](@ref) instead; the operators have the same names and meaning.

# Construction

    CoorSystemSym(OM, coords, tmp_coords = (), params = ();
                  rules = Dict(), tmp_var = Dict(), to_coords = Dict())
    CoorSystemSym(OM, coords, bnorm, χᵢ, tmp_coords = (), params = (); ...)

The first form derives the natural basis from `𝐚ᵢ = ∂ᵢOM`; the second takes the
normalized basis and the Lamé coefficients directly, which is much faster when
they are known in closed form.

The optional arguments drive symbolic simplification, without which a
non-trivial chart produces correct but unusable nested radicals:

| Argument | Role |
|:--|:--|
| `tmp_coords` | auxiliary symbols standing for compound expressions |
| `params` | constants appearing in `OM` |
| `tmp_var` | substitutions replacing expressions by those symbols |
| `to_coords` | how to eliminate them again before differentiating |
| `rules` | rewrite rules applied after each simplification |

See [Adding a coordinate system](@ref dev-adding-coorsystem) for the recipe and
[`coorsys_spheroidal`](@ref) for a worked example.
"""
struct CoorSystemSym{dim, T <: Number, VEC, BNORM, BNAT} <: AbstractCoorSystem{dim, T}
    core::ChartCore{dim, T, VEC, BNORM, BNAT}
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
            ChartCore{dim, T, VEC, typeof(normalized_basis), typeof(natural_basis)}(
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
        # Plain componentwise derivative of the position vector — deliberately
        # NOT the two-argument `∂`.
        #
        # `@set_coorsys` installs a method `∂(::AbstractTens, ::Sym)` bound to
        # whatever system was made the default, and it is more specific than the
        # variadic fallback, so it wins here. It then looks `coords[i]` up among
        # *its own* coordinates, does not find it, and returns `zero(t)`: every
        # natural vector comes out null, the frame matrix is singular, and the
        # constructor dies with `NonInvertibleMatrixError`. In other words,
        # building any new coordinate system after a `@set_coorsys` used to
        # fail. `SubManifoldSym` carries the same fix.
        ∂coord(t, xᵢ) = change_tens(Tens(tdiff(components_canon(t), xᵢ)), get_basis(t), get_var(t))
        aᵢ = ntuple(i -> simp(chvar(∂coord(OMc, coords[i]), tmp_var)), dim)
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
            ChartCore{dim, T, VEC, typeof(normalized_basis), typeof(natural_basis)}(
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
        )
    end
end

# These four were restricted to `CoorSystemSym` while every neighboring
# accessor below dispatches on `AbstractCoorSystem`.  `SubManifoldSym` is also
# an `AbstractCoorSystem` and carries the very same fields, so the restriction
# made `∂`/`GRAD`/`DIV`/`LAPLACE`/`HESS` on a submanifold fail with a
# `MethodError` on `only_coords`.
with_tmp_var(CS::AbstractCoorSystem, t) = length(core(CS).tmp_var) > 0 ? tsubs(t, core(CS).tmp_var...) : t
only_coords(CS::AbstractCoorSystem, t) = length(core(CS).to_coords) > 0 ? tsubs(t, core(CS).to_coords...) : t

"""
    getcoords(CS::AbstractCoorSystem) → NTuple
    getcoords(CS::AbstractCoorSystem, i::Integer)

Coordinate symbols of `CS`, or the `i`-th one.

Note the ordering of the spherical system, `(θ, ϕ, r)` and not `(r, θ, ϕ)`, so
that `θ = ϕ = 0` reproduces the canonical basis in the canonical order.
"""
getcoords(CS::AbstractCoorSystem) = core(CS).coords
getcoords(CS::AbstractCoorSystem, i::Integer) = getcoords(CS)[i]

@pure get_dim(::AbstractCoorSystem{dim}) where {dim} = dim

"""
    nderiv(CS::AbstractCoorSystem) → Int

Number of independent differentiation directions of `CS`.

`dim` for a chart of the whole space, `dim-1` for a [`SubManifoldSym`](@ref),
which differentiates along its surface coordinates only. Every differential
operator loops to `nderiv(CS)`, which is what lets a single implementation serve
both.
"""
@pure nderiv(::AbstractCoorSystem{dim}) where {dim} = dim

"""
    getOM(CS::AbstractCoorSystem) → AbstractTens{1}

Position vector of `CS` as a function of its coordinates — the map the natural
basis, the Lamé coefficients and the Christoffel symbols are all derived from.
"""
getOM(CS::AbstractCoorSystem) = core(CS).OM

normalized_basis(CS::AbstractCoorSystem) = core(CS).normalized_basis
natural_basis(CS::AbstractCoorSystem) = core(CS).natural_basis

"""
    Lame(CS::AbstractCoorSystem) → NTuple
    Lame(CS::CoorSystemNum, x₀::AbstractVector) → Vector

Lamé coefficients `χᵢ = ‖𝐚ᵢ‖`, the norms of the natural basis vectors.

They relate the natural and normalized bases, `𝐞ᵢ = 𝐚ᵢ/χᵢ`, and for an
orthogonal system give the line element `ds² = Σᵢ χᵢ² (dqⁱ)²`. Symbolic systems
return expressions; [`CoorSystemNum`](@ref) evaluates them at a point.

`Lame(coorsys_spherical())` is `(r, r sin(θ), 1)`.

See also [Curvilinear differential calculus](@ref th-curvilinear).
"""
Lame(CS::AbstractCoorSystem) = core(CS).χᵢ
"""
    Christoffel(CS::AbstractCoorSystem) → Array{T,3}
    Christoffel(CS::CoorSystemNum, x₀::AbstractVector) → Array{T,3}

Christoffel symbols of `CS`, `Γᵏᵢⱼ = ∂ᵢ𝐚ⱼ ⋅ 𝐚ᵏ`, symmetric in `(i,j)`.

!!! note "Storage convention"
    The array is indexed `Γ[i,j,k]` `= Γᵏᵢⱼ` — the **contravariant index last**.
    The same convention is used by the `Γ_func` closure of
    [`CoorSystemNum`](@ref).

They are what distinguishes a derivative on a curvilinear chart from a plain
partial derivative; every one of them vanishes exactly for a Cartesian chart.

See also [`Lame`](@ref), [Curvilinear differential calculus](@ref th-curvilinear).
"""
Christoffel(CS::AbstractCoorSystem) = core(CS).Γ
# simprules(t, CS::AbstractCoorSystem) = length(CS.rules) > 0 ? tsimplify(tsubs(tsimplify(t), CS.rules...)) : tsimplify(t)
simprules(t, CS::AbstractCoorSystem) = length(core(CS).rules) > 0 ? tsubs(t, core(CS).rules...) : t

natvec(CS::AbstractCoorSystem, ::Val{:cov}) = core(CS).aᵢ
natvec(CS::AbstractCoorSystem, ::Val{:cont}) = core(CS).aⁱ
natvec(CS::AbstractCoorSystem, var = :cov) = natvec(CS, Val(var))
natvec(CS::AbstractCoorSystem, i::Integer, var = :cov) = natvec(CS, var)[i]

unitvec(CS::AbstractCoorSystem) = core(CS).eᵢ
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
        CS::AbstractCoorSystem{dim},
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

∂(t::T, i::Integer, CS::AbstractCoorSystem{dim}) where {dim, T <: SymType} =
    tdiff(only_coords(CS, t), getcoords(CS, i))

# `T` stays coupled to the tensor here, unlike the other operators: loosening it
# to `AbstractTens{order,dim}` makes this method ambiguous with the variadic
# plain-derivative fallback `∂(t::AbstractTens{order,dim,T}, xᵢ...)`, which can
# absorb `(x, CS)` as its trailing arguments.
function ∂(
        t::AbstractTens{order, dim, T},
        x::T,
        CS::AbstractCoorSystem{dim},
    ) where {order, dim, T <: SymType}
    ind = findfirst(i -> i == x, getcoords(CS))
    return ind === nothing ? zero(t) : ∂(t, ind, CS)
end

function ∂(
        t::T,
        x::T,
        CS::AbstractCoorSystem{dim},
    ) where {dim, T <: SymType}
    ind = findfirst(i -> i == x, getcoords(CS))
    return ind === nothing ? zero(t) : ∂(t, ind, CS)
end

"""
    GRAD(t::Union{T,AbstractTens{order,dim,T}}, CS::CoorSystemSym{dim,T}) where {order,dim,T<:SymType}

Gradient of the scalar or tensor field `t` with respect to the coordinate
system `CS`, raising the order by one:

    GRAD(t) = Σᵢ ∂ᵢt ⊗ 𝐚ⁱ

!!! note "The derivative index comes last"
    The dual natural vector `𝐚ⁱ` is appended **on the right**, so for a vector
    field `𝐯` the components of `GRAD(𝐯)` are `(∇𝐯)ᵢⱼ = ∂ⱼvᵢ`. This is the
    convention that makes [`DIV`](@ref) the contraction of the *last* index and
    `LAPLACE = DIV ∘ GRAD` come out right; a library using the opposite
    convention differs by a transpose.

If `CS` has been made the default with [`@set_coorsys`](@ref), the second
argument may be omitted.

# Examples
```julia
julia> Spherical = coorsys_spherical() ; θ, ϕ, r = getcoords(Spherical) ;

julia> @set_coorsys Spherical

julia> GRAD(r)          # = 𝐞ʳ

julia> 𝐞ᶿ, 𝐞ᵠ, 𝐞ʳ = unitvec(Spherical) ; GRAD(𝐞ʳ)   # = (𝐞ᶿ⊗𝐞ᶿ + 𝐞ᵠ⊗𝐞ᵠ)/r
```

See also [`SYMGRAD`](@ref), [`DIV`](@ref), [`LAPLACE`](@ref), [`HESS`](@ref).
"""
GRAD(t::Union{T, AbstractTens{order, dim}}, CS::AbstractCoorSystem{dim}) where {order, dim, T <: SymType} =
    sum([∂(t, i, CS) ⊗ natvec(CS, i, :cont) for i in 1:nderiv(CS)])


"""
    SYMGRAD(t::Union{T,AbstractTens{order,dim,T}}, CS::CoorSystemSym{dim,T}) where {order,dim,T<:SymType}

Symmetrized gradient of `t` with respect to the coordinate system `CS`:

    SYMGRAD(t) = Σᵢ ∂ᵢt ⊗ˢ 𝐚ⁱ

Applied to a displacement field this is the linearized strain tensor,
`𝛆 = (∇𝛏 + ᵗ∇𝛏)/2`, which is why it is a primitive rather than a composition of
[`GRAD`](@ref) and a transpose.

If `CS` has been made the default with [`@set_coorsys`](@ref), the second
argument may be omitted.

# Examples
```julia
julia> Cylindrical = coorsys_cylindrical() ; r, θ, z = getcoords(Cylindrical) ;

julia> 𝐞ʳ, 𝐞ᶿ, 𝐞ᶻ = unitvec(Cylindrical) ; @set_coorsys Cylindrical

julia> 𝛏 = SymFunction("ξʳ", real = true)(r, z) * 𝐞ʳ + SymFunction("ξᶻ", real = true)(r, z) * 𝐞ᶻ ;

julia> SYMGRAD(𝛏)       # axisymmetric strain tensor, with εᶿᶿ = ξʳ/r
```

See also [`GRAD`](@ref), [`DIV`](@ref).
"""
SYMGRAD(
    t::Union{T, AbstractTens{order, dim}},
    CS::AbstractCoorSystem{dim},
) where {order, dim, T <: SymType} = sum([∂(t, i, CS) ⊗ˢ natvec(CS, i, :cont) for i in 1:nderiv(CS)])


"""
    DIV(t::AbstractTens{order,dim,T}, CS::CoorSystemSym{dim,T}) where {order,dim,T<:SymType}

Divergence of the tensor field `t` (of order ≥ 1) with respect to the coordinate
system `CS`, lowering the order by one:

    DIV(t) = Σᵢ ∂ᵢt ⋅ 𝐚ⁱ

!!! note "The contracted index is the last one"
    Consistently with [`GRAD`](@ref), the contraction acts on the **last** index
    of `t`: for an order-2 field, `DIV(𝛔)ᵢ = ∂ⱼσᵢⱼ`. For a symmetric field such
    as a stress tensor the distinction is immaterial, but it is not for a general
    order-2 field.

If `CS` has been made the default with [`@set_coorsys`](@ref), the second
argument may be omitted.

# Examples
```julia
julia> Spherical = coorsys_spherical() ; θ, ϕ, r = getcoords(Spherical) ;

julia> 𝐞ᶿ, 𝐞ᵠ, 𝐞ʳ = unitvec(Spherical) ; @set_coorsys Spherical

julia> DIV(𝐞ʳ)          # = 2/r

julia> 𝛔 = SymFunction("σʳʳ", real = true)(r) * 𝐞ʳ ⊗ 𝐞ʳ ;

julia> simplify(DIV(𝛔)) # radial equilibrium operator
```

See also [`GRAD`](@ref), [`LAPLACE`](@ref).
"""
DIV(t::AbstractTens{order, dim, T}, CS::AbstractCoorSystem{dim}) where {order, dim, T <: SymType} =
    sum([∂(t, i, CS) ⋅ natvec(CS, i, :cont) for i in 1:nderiv(CS)])


"""
    LAPLACE(t::Union{T,AbstractTens{order,dim,T}}, CS::CoorSystemSym{dim,T}) where {order,dim,T<:SymType}

Laplacian of the scalar or tensor field `t`, defined as the composition

    LAPLACE(t) = DIV(GRAD(t))

and therefore preserving the order of `t`.

If `CS` has been made the default with [`@set_coorsys`](@ref), the second
argument may be omitted.

# Examples
```julia
julia> Polar = coorsys_polar() ; r, θ = getcoords(Polar) ; @set_coorsys Polar

julia> LAPLACE(SymFunction("f", real = true)(r, θ))   # the polar Laplacian

julia> n = symbols("n", integer = true) ; simplify(LAPLACE(r^n * cos(n*θ)))   # = 0
```

See also [`GRAD`](@ref), [`DIV`](@ref), [`HESS`](@ref).
"""
LAPLACE(
    t::Union{T, AbstractTens{order, dim}},
    CS::AbstractCoorSystem{dim},
) where {order, dim, T <: SymType} = DIV(GRAD(t, CS), CS)

"""
    HESS(t::Union{T,AbstractTens{order,dim,T}}, CS::CoorSystemSym{dim,T}) where {order,dim,T<:SymType}

Hessian of `t`, defined as the second gradient

    HESS(t) = GRAD(GRAD(t))

and therefore raising the order of `t` by two. For a scalar field the trace of
`HESS(t)` is [`LAPLACE`](@ref)`(t)`.

If `CS` has been made the default with [`@set_coorsys`](@ref), the second
argument may be omitted.

# Examples
```julia
julia> Spherical = coorsys_spherical() ; θ, ϕ, r = getcoords(Spherical) ;

julia> @set_coorsys Spherical

julia> simplify(HESS(1/r))    # the kernel of the 3-D Laplace equation
```

See also [`GRAD`](@ref), [`LAPLACE`](@ref).
"""
HESS(t::Union{T, AbstractTens{order, dim}}, CS::AbstractCoorSystem{dim}) where {order, dim, T <: SymType} =
    GRAD(GRAD(t, CS), CS)

"""
    coorsys_cartesian(coords = symbols("x y z", real = true))

Return the cartesian coordinate system, in which the natural and normalized
bases both coincide with the canonical one and all Christoffel symbols vanish.

# Examples

The divergence of a general symmetric order-2 field, which in cartesian
coordinates reduces to the plain sum of partial derivatives:

```julia
julia> Cartesian = coorsys_cartesian() ; 𝐗 = getcoords(Cartesian) ;

julia> ℬ = normalized_basis(Cartesian) ;

julia> 𝛔 = Tens(SymmetricTensor{2,3}((i, j) -> SymFunction("σ\$i\$j", real = true)(𝐗...))) ;

julia> get_array(DIV(𝛔, Cartesian))
 Derivative(σ11(x, y, z), x) + Derivative(σ21(x, y, z), y) + Derivative(σ31(x, y, z), z)
 Derivative(σ21(x, y, z), x) + Derivative(σ22(x, y, z), y) + Derivative(σ32(x, y, z), z)
 Derivative(σ31(x, y, z), x) + Derivative(σ32(x, y, z), y) + Derivative(σ33(x, y, z), z)
```

which is `DIV(𝛔)ᵢ = ∂ⱼσᵢⱼ`.

See also [`coorsys_polar`](@ref), [`coorsys_cylindrical`](@ref),
[`coorsys_spherical`](@ref), [`@set_coorsys`](@ref).
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

julia> print_tensor(GRAD(𝐞ʳ),vec)
(1/r)𝐞ᶿ⊗𝐞ᶿ + (1/r)𝐞ᵠ⊗𝐞ᵠ

julia> print_tensor(DIV(𝐞ʳ ⊗ 𝐞ʳ),vec)
(2/r)𝐞ʳ

julia> LAPLACE(1/r)
0
``` 
"""
# ── The default coordinate system ────────────────────────────────────────────
#
# `@set_coorsys` used to `@eval` single-argument methods for `∂`, `GRAD`, `DIV`,
# … straight into this module. That is global mutation of the method table from
# a macro, and it caused two concrete failures:
#
#   * `∂(t, x)` already means "plain derivative with respect to the symbol `x`"
#     (see the `∂` docstring). The generated method was *more specific* than
#     that fallback, so after any `@set_coorsys` the same call silently became a
#     covariant derivative in the default chart — which then failed to find `x`
#     among its own coordinates and returned `zero(t)`. Every constructor that
#     differentiates a position vector was broken by it.
#   * it invalidates compiled code and is invisible at the call site.
#
# The default is now a plain `Ref`, and the single-argument methods are defined
# once, below. `∂` deliberately gets **no** default-chart method: it keeps its
# one unambiguous meaning.

const _DEFAULT_COORSYS = Ref{Any}(nothing)
const _DEFAULT_COORSYS_VEC = Ref{Char}('𝐞')
const _DEFAULT_COORSYS_COORDS = Ref{Any}(nothing)

"""
    set_coorsys!(CS; vec = '𝐞', coords = nothing) → CS

Make `CS` the default coordinate system, so that [`GRAD`](@ref),
[`SYMGRAD`](@ref), [`DIV`](@ref), [`LAPLACE`](@ref), [`HESS`](@ref) and
[`print_tensor`](@ref) may be called with a single argument.

`vec` and `coords` control how [`print_tensor`](@ref) names the basis vectors.
The macro [`@set_coorsys`](@ref) is a thin wrapper over this function.

!!! note "`∂` is deliberately not affected"
    `∂(t, x)` always means the plain derivative of `t` with respect to the
    symbol `x`, whether or not a default chart is set. For the covariant
    derivative, pass the system explicitly: `∂(t, x, CS)`.

See also [`default_coorsys`](@ref), [`unset_coorsys!`](@ref).
"""
function set_coorsys!(CS::AbstractCoorSystem; vec::Char = '𝐞', coords = nothing)
    _DEFAULT_COORSYS[] = CS
    _DEFAULT_COORSYS_VEC[] = vec
    _DEFAULT_COORSYS_COORDS[] = coords
    return CS
end

"""
    default_coorsys() → AbstractCoorSystem

The coordinate system installed by [`set_coorsys!`](@ref) or
[`@set_coorsys`](@ref). Throws if none has been set.
"""
function default_coorsys()
    CS = _DEFAULT_COORSYS[]
    CS === nothing && error(
        "no default coordinate system is set: call `set_coorsys!(CS)` " *
            "(or `@set_coorsys CS`), or pass the system explicitly as the last argument."
    )
    return CS
end

"""
    unset_coorsys!()

Forget the default coordinate system installed by [`set_coorsys!`](@ref).
"""
function unset_coorsys!()
    _DEFAULT_COORSYS[] = nothing
    _DEFAULT_COORSYS_COORDS[] = nothing
    return nothing
end

# Single-argument forms, defined once and for all.
GRAD(t::Union{T, AbstractTens}) where {T <: SymType} = GRAD(t, default_coorsys())
SYMGRAD(t::Union{T, AbstractTens}) where {T <: SymType} = SYMGRAD(t, default_coorsys())
DIV(t::AbstractTens) = DIV(t, default_coorsys())
LAPLACE(t::Union{T, AbstractTens}) where {T <: SymType} = LAPLACE(t, default_coorsys())
HESS(t::Union{T, AbstractTens}) where {T <: SymType} = HESS(t, default_coorsys())

"""
    @set_coorsys CS [vec] [coords]

Make `CS` the default coordinate system for the differential operators, so they
can be called with one argument.

Equivalent to [`set_coorsys!`](@ref)`(CS; vec, coords)`; kept as a macro for
familiarity. Unlike earlier versions it **does not define any method** — it only
stores the system — and it leaves the meaning of `∂(t, x)` untouched.

# Examples
```julia
julia> Spherical = coorsys_spherical() ; θ, ϕ, r = getcoords(Spherical) ;

julia> @set_coorsys Spherical

julia> LAPLACE(1 / r)      # one argument
0
```

See also [`set_coorsys!`](@ref), [`default_coorsys`](@ref),
[`unset_coorsys!`](@ref).
"""
macro set_coorsys(CS = coorsys_cartesian(), vec = '𝐞', coords = nothing)
    m = @__MODULE__
    return quote
        $m.set_coorsys!($(esc(CS)); vec = $(esc(vec)), coords = $(esc(coords)))
    end
end

# One-argument `print_tensor` uses the default chart when one is set, so that a
# result prints with the chart's own coordinate names (`𝐞ʳ`, `𝐞ᶿ`, …) rather
# than with numbered basis vectors. Without a default it falls through to the
# generic method in `tens.jl`.
function print_tensor(t::AbstractTens{order, dim, T}) where {order, dim, T <: SymType}
    CS = _DEFAULT_COORSYS[]
    (CS === nothing || get_dim(CS) != dim) &&
        return print_tensor(t; vec = _DEFAULT_COORSYS_VEC[], coords = ntuple(i -> i, dim))
    coords = _DEFAULT_COORSYS_COORDS[]
    coords === nothing && (coords = string.(getcoords(CS)))
    length(coords) == dim - 1 && (coords = (coords..., dim))
    return print_tensor(change_tens(t, normalized_basis(CS)); vec = _DEFAULT_COORSYS_VEC[], coords = coords)
end

function print_tensor(t::AbstractTens{order, dim, T}, CS::AbstractCoorSystem; vec = '𝐞') where {order, dim, T}
    coords = string.(getcoords(CS))
    ℬ = normalized_basis(CS)
    return print_tensor(change_tens(t, ℬ); vec = vec, coords = coords)
end

export ∂, CoorSystemSym, Lame, Christoffel, nderiv
export GRAD, SYMGRAD, DIV, LAPLACE, HESS
export normalized_basis, natural_basis, natvec, unitvec, getcoords, getOM
export coorsys_cartesian, coorsys_polar, coorsys_cylindrical, coorsys_spherical, coorsys_spheroidal
export @set_coorsys, set_coorsys!, default_coorsys, unset_coorsys!
