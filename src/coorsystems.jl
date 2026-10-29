abstract type AbstractCoorSystem{dim, T <: Number} <: Any end

"""
    AbstractGeometry

What a [`CoordinateSystem`](@ref) knows of its chart, and so how its differential
operators work: [`ExactGeometry`](@ref) holds closed-form expressions of the
coordinate symbols, [`PointwiseGeometry`](@ref) functions of the point.
"""
abstract type AbstractGeometry end


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
    ExactGeometry{dim,T,VEC,BNORM,BNAT} <: AbstractGeometry

The geometry of a chart in closed form, as expressions of its coordinate symbols
(SymPy's `Sym` or Symbolics' `Num`): the position vector, the coordinates, the
two bases, the natural and dual frames, the Lamé coefficients, the connection
array, and the symbolic-simplification settings. Differential operators on it
take and return expressions.

Held by both a [`CoordinateSystem`](@ref) and a [`Submanifold`](@ref) — which
differ only in what they *add* to it (a submanifold adds the two fundamental
forms) and in how many directions they differentiate along ([`nderiv`](@ref)).
Declaring the block once keeps the two in step: a fix to one is a fix to
the other.

Reached through `core` (internal); user code should go through the accessors
([`getcoords`](@ref), [`lame_coefficients`](@ref), [`christoffel`](@ref), …) instead.
"""
struct ExactGeometry{dim, T, VEC, BNORM, BNAT} <: AbstractGeometry
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
    core(CS::AbstractCoorSystem) → ExactGeometry

The closed-form geometry of `CS`. Internal: prefer the named accessors.
"""
core(CS::AbstractCoorSystem) = CS.geometry

"""
    CoordinateSystem{dim,T,G<:AbstractGeometry} <: AbstractCoorSystem{dim,T}

Curvilinear coordinate system of the `dim`-dimensional space, of one of two
geometries:

| Geometry | Fields | Operators |
|:--|:--|:--|
| [`ExactGeometry`](@ref) | expressions of the coordinate symbols (`Sym` or `Num`) | take and return expressions, closed-form results valid everywhere |
| [`PointwiseGeometry`](@ref) | functions of the point | take a function of the point, return one, by automatic differentiation |

The differential operators [`grad`](@ref), [`symgrad`](@ref), [`divergence`](@ref),
[`laplacian`](@ref) and [`hess`](@ref) have the same names and meaning on both.

Predefined systems: [`coorsys_cartesian`](@ref), [`coorsys_polar`](@ref),
[`coorsys_cylindrical`](@ref), [`coorsys_spherical`](@ref),
[`coorsys_spheroidal`](@ref) — exact by default, pointwise with
`PointwiseGeometry` as first argument (e.g. `coorsys_polar(PointwiseGeometry)`).
[`@set_coorsys`](@ref) makes an exact one the default so the operators take a
single argument.

# Construction

    CoordinateSystem(OM, coords, tmp_coords = (), params = ();
                     rules = Dict(), tmp_var = Dict(), to_coords = Dict())
    CoordinateSystem(OM, coords, bnorm, χᵢ, tmp_coords = (), params = (); ...)
    CoordinateSystem(OM_func::Function, dim::Integer, T::Type = Float64)

The first two build an exact system from the position vector `OM`, an
expression of the coordinate symbols `coords`. The first derives the natural
basis from `𝐚ᵢ = ∂ᵢOM`; the second takes the normalized basis and the Lamé
coefficients directly, which is much faster when they are known in closed form.
The third builds a pointwise system from `OM_func`, a function of the point
returning its Cartesian components, by automatic differentiation.

The optional arguments of the exact forms drive symbolic simplification,
without which a non-trivial chart produces correct but unusable nested radicals:

| Argument | Role |
|:--|:--|
| `tmp_coords` | auxiliary symbols standing for compound expressions |
| `params` | constants appearing in `OM` |
| `tmp_var` | substitutions replacing expressions by those symbols |
| `to_coords` | how to eliminate them again before differentiating |
| `rules` | rewrite rules applied after each simplification |

See [Adding a coordinate system](@ref dev-adding-coorsystem) for the recipe and
[`coorsys_spheroidal`](@ref) for a worked example.

!!! note "Renamed in 0.8"
    `CoorSystemSym` and `CoorSystemNum` are the two geometries of this one type;
    the former names still work, with a deprecation warning.
"""
struct CoordinateSystem{dim, T <: Number, G <: AbstractGeometry} <: AbstractCoorSystem{dim, T}
    geometry::G
end

# The exact systems: a closed-form geometry, built by the two constructors below.
const _ExactSystem{dim, T} = CoordinateSystem{dim, T, <:ExactGeometry}
_exact_system(g::ExactGeometry{dim, T}) where {dim, T} = CoordinateSystem{dim, T, typeof(g)}(g)

function CoordinateSystem(
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
    return _exact_system(
        ExactGeometry{dim, T, VEC, typeof(normalized_basis), typeof(natural_basis)}(
            OM, coords, normalized_basis, natural_basis, aᵢ, χᵢ, aⁱ, eᵢ, Γ,
            tmp_coords, params, rules, tmp_var, to_coords,
        )
    )
end

# `Tuple{T, Vararg{T, M}}` rather than `NTuple{dim, T}`: a coordinate system
# with no coordinate is not a thing, and the empty tuple the latter also
# accepts leaves `T` bound by nothing.
function CoordinateSystem(
        OM::VEC,
        coords::Tuple{T, Vararg{T, M}},
        tmp_coords::NTuple = (),
        params::NTuple = ();
        rules::Dict = Dict(),
        tmp_var::Dict = Dict(),
        to_coords::Dict = Dict()
    ) where {T, M, VEC}
    dim = M + 1
    simp(t) = length(rules) > 0 ? tsimplify(tsubs(t, rules...)) : tsimplify(t)
    chvar(t, d) = length(d) > 0 ? tsubs(t, d...) : t
    OMc = chvar(OM, to_coords)
    # Plain componentwise derivative of the position vector — deliberately
    # NOT the two-argument `∂`.
    #
    # The natural vectors are the derivatives of the canonical components
    # of the position vector, taken here directly: what they are must not
    # depend on any method a chart adds to `∂`. `Submanifold` does the same.
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
    return _exact_system(
        ExactGeometry{dim, T, VEC, typeof(normalized_basis), typeof(natural_basis)}(
            OMc, coords, normalized_basis, natural_basis, aᵢ, χᵢ, aⁱ, eᵢ, Γ,
            tmp_coords, params, rules, tmp_var, to_coords,
        )
    )
end

# On `AbstractCoorSystem`, like every accessor below: `Submanifold` carries
# the same fields, and `∂`/`grad`/`divergence`/`laplacian`/`hess` on a submanifold reach
# `only_coords`.
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

get_dim(::AbstractCoorSystem{dim}) where {dim} = dim

"""
    nderiv(CS::AbstractCoorSystem) → Int

Number of independent differentiation directions of `CS`.

`dim` for a chart of the whole space, `dim-1` for a [`Submanifold`](@ref),
which differentiates along its surface coordinates only. Every differential
operator loops to `nderiv(CS)`, which is what lets a single implementation serve
both.
"""
nderiv(::AbstractCoorSystem{dim}) where {dim} = dim

"""
    position_vector(CS::AbstractCoorSystem) → AbstractTens{1}

Position vector of `CS` as a function of its coordinates — the map the natural
basis, the Lamé coefficients and the Christoffel symbols are all derived from.
"""
position_vector(CS::AbstractCoorSystem) = core(CS).OM

normalized_basis(CS::AbstractCoorSystem) = core(CS).normalized_basis
natural_basis(CS::AbstractCoorSystem) = core(CS).natural_basis

"""
    lame_coefficients(CS::AbstractCoorSystem) → NTuple
    lame_coefficients(CS::CoordinateSystem{dim,T,<:PointwiseGeometry}, x₀::AbstractVector) → Vector

Lamé coefficients `χᵢ = ‖𝐚ᵢ‖`, the norms of the natural basis vectors.

They relate the natural and normalized bases, `𝐞ᵢ = 𝐚ᵢ/χᵢ`, and for an
orthogonal system give the line element `ds² = Σᵢ χᵢ² (dqⁱ)²`. Symbolic systems
return expressions; a [`PointwiseGeometry`](@ref) evaluates them at a point.

`lame_coefficients(coorsys_spherical())` is `(r, r sin(θ), 1)`.

See also [Curvilinear differential calculus](@ref th-curvilinear).
"""
lame_coefficients(CS::AbstractCoorSystem) = core(CS).χᵢ
"""
    christoffel(CS::AbstractCoorSystem) → Array{T,3}
    christoffel(CS::CoordinateSystem{dim,T,<:PointwiseGeometry}, x₀::AbstractVector) → Array{T,3}

Christoffel symbols of `CS`, `Γᵏᵢⱼ = ∂ᵢ𝐚ⱼ ⋅ 𝐚ᵏ`, symmetric in `(i,j)`.

!!! note "Storage convention"
    The array is indexed `Γ[i,j,k]` `= Γᵏᵢⱼ` — the **contravariant index last**.
    The same convention is used by the `Γ_func` closure of
    [`PointwiseGeometry`](@ref).

They are what distinguishes a derivative on a curvilinear chart from a plain
partial derivative; every one of them vanishes exactly for a Cartesian chart.

See also [`lame_coefficients`](@ref), [Curvilinear differential calculus](@ref th-curvilinear).
"""
christoffel(CS::AbstractCoorSystem) = core(CS).Γ
# simprules(t, CS::AbstractCoorSystem) = length(CS.rules) > 0 ? tsimplify(tsubs(tsimplify(t), CS.rules...)) : tsimplify(t)
simprules(t, CS::AbstractCoorSystem) = length(core(CS).rules) > 0 ? tsubs(t, core(CS).rules...) : t

natvec(CS::AbstractCoorSystem, ::Val{:cov}) = core(CS).aᵢ
natvec(CS::AbstractCoorSystem, ::Val{:cont}) = core(CS).aⁱ
natvec(CS::AbstractCoorSystem, var::Symbol = :cov) = _by_variance(natvec, var, CS)
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
    # Γᵢⱼᵏ = Γᵢⱼₗ gˡᵏ: the last index raised by the inverse metric.
    return _contract_trailing(Γᵢⱼₖ, gⁱʲ, Val(1))
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
    Γ = christoffel(CS)
    data = tdiff(t, getcoords(CS, i))
    for o in 1:order
        data += _contract_index(t, view(Γ, i, :, :), o)
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

# Two methods, not one union. `Union{SymType, AbstractTens{order, dim}}` leaves
# `order` and `dim` bound by nothing when the argument is a `Num` or a `Sym`,
# and spelling the union the other way round, `Union{T, AbstractTens}`, leaves
# `T` unbound when it is a tensor. Split, each signature binds what it names.
"""
    grad(t::Union{T,AbstractTens{order,dim,T}}, CS::CoordinateSystem{dim,T,<:ExactGeometry}) where {order,dim,T<:SymType}

Gradient of the scalar or tensor field `t` with respect to the coordinate
system `CS`, raising the order by one:

    grad(t) = Σᵢ ∂ᵢt ⊗ 𝐚ⁱ

!!! note "The derivative index comes last"
    The dual natural vector `𝐚ⁱ` is appended **on the right**, so for a vector
    field `𝐯` the components of `grad(𝐯)` are `(∇𝐯)ᵢⱼ = ∂ⱼvᵢ`. This is the
    convention that makes [`divergence`](@ref) the contraction of the *last* index and
    `laplacian = divergence ∘ grad` come out right; a library using the opposite
    convention differs by a transpose.

If `CS` has been made the default with [`@set_coorsys`](@ref), the second
argument may be omitted.

# Examples
```julia
julia> Spherical = coorsys_spherical() ; θ, ϕ, r = getcoords(Spherical) ;

julia> @set_coorsys Spherical

julia> grad(r)          # = 𝐞ʳ

julia> 𝐞ᶿ, 𝐞ᵠ, 𝐞ʳ = unitvec(Spherical) ; grad(𝐞ʳ)   # = (𝐞ᶿ⊗𝐞ᶿ + 𝐞ᵠ⊗𝐞ᵠ)/r
```

See also [`symgrad`](@ref), [`divergence`](@ref), [`laplacian`](@ref), [`hess`](@ref).
"""
grad(t::AbstractTens{order, dim}, CS::AbstractCoorSystem{dim}) where {order, dim} =
    sum([∂(t, i, CS) ⊗ natvec(CS, i, :cont) for i in 1:nderiv(CS)])
grad(t::SymType, CS::AbstractCoorSystem) =
    sum([∂(t, i, CS) ⊗ natvec(CS, i, :cont) for i in 1:nderiv(CS)])


# Two methods, not one union. `Union{SymType, AbstractTens{order, dim}}` leaves
# `order` and `dim` bound by nothing when the argument is a `Num` or a `Sym`,
# and spelling the union the other way round, `Union{T, AbstractTens}`, leaves
# `T` unbound when it is a tensor. Split, each signature binds what it names.
"""
    symgrad(t::Union{T,AbstractTens{order,dim,T}}, CS::CoordinateSystem{dim,T,<:ExactGeometry}) where {order,dim,T<:SymType}

Symmetrized gradient of `t` with respect to the coordinate system `CS`:

    symgrad(t) = Σᵢ ∂ᵢt ⊗ˢ 𝐚ⁱ

Applied to a displacement field this is the linearized strain tensor,
`𝛆 = (∇𝛏 + ᵗ∇𝛏)/2`, which is why it is a primitive rather than a composition of
[`grad`](@ref) and a transpose.

If `CS` has been made the default with [`@set_coorsys`](@ref), the second
argument may be omitted.

# Examples
```julia
julia> Cylindrical = coorsys_cylindrical() ; r, θ, z = getcoords(Cylindrical) ;

julia> 𝐞ʳ, 𝐞ᶿ, 𝐞ᶻ = unitvec(Cylindrical) ; @set_coorsys Cylindrical

julia> 𝛏 = SymFunction("ξʳ", real = true)(r, z) * 𝐞ʳ + SymFunction("ξᶻ", real = true)(r, z) * 𝐞ᶻ ;

julia> symgrad(𝛏)       # axisymmetric strain tensor, with εᶿᶿ = ξʳ/r
```

See also [`grad`](@ref), [`divergence`](@ref).
"""
symgrad(t::AbstractTens{order, dim}, CS::AbstractCoorSystem{dim}) where {order, dim} =
    sum([∂(t, i, CS) ⊗ˢ natvec(CS, i, :cont) for i in 1:nderiv(CS)])
symgrad(t::SymType, CS::AbstractCoorSystem) =
    sum([∂(t, i, CS) ⊗ˢ natvec(CS, i, :cont) for i in 1:nderiv(CS)])


"""
    divergence(t::AbstractTens{order,dim,T}, CS::CoordinateSystem{dim,T,<:ExactGeometry}) where {order,dim,T<:SymType}

Divergence of the tensor field `t` (of order ≥ 1) with respect to the coordinate
system `CS`, lowering the order by one:

    divergence(t) = Σᵢ ∂ᵢt ⋅ 𝐚ⁱ

!!! note "The contracted index is the last one"
    Consistently with [`grad`](@ref), the contraction acts on the **last** index
    of `t`: for an order-2 field, `divergence(𝛔)ᵢ = ∂ⱼσᵢⱼ`. For a symmetric field such
    as a stress tensor the distinction is immaterial, but it is not for a general
    order-2 field.

If `CS` has been made the default with [`@set_coorsys`](@ref), the second
argument may be omitted.

# Examples
```julia
julia> Spherical = coorsys_spherical() ; θ, ϕ, r = getcoords(Spherical) ;

julia> 𝐞ᶿ, 𝐞ᵠ, 𝐞ʳ = unitvec(Spherical) ; @set_coorsys Spherical

julia> divergence(𝐞ʳ)          # = 2/r

julia> 𝛔 = SymFunction("σʳʳ", real = true)(r) * 𝐞ʳ ⊗ 𝐞ʳ ;

julia> simplify(divergence(𝛔)) # radial equilibrium operator
```

See also [`grad`](@ref), [`laplacian`](@ref).
"""
divergence(t::AbstractTens{order, dim, T}, CS::AbstractCoorSystem{dim}) where {order, dim, T <: SymType} =
    sum([∂(t, i, CS) ⋅ natvec(CS, i, :cont) for i in 1:nderiv(CS)])


# Two methods, not one union. `Union{SymType, AbstractTens{order, dim}}` leaves
# `order` and `dim` bound by nothing when the argument is a `Num` or a `Sym`,
# and spelling the union the other way round, `Union{T, AbstractTens}`, leaves
# `T` unbound when it is a tensor. Split, each signature binds what it names.
"""
    laplacian(t::Union{T,AbstractTens{order,dim,T}}, CS::CoordinateSystem{dim,T,<:ExactGeometry}) where {order,dim,T<:SymType}

Laplacian of the scalar or tensor field `t`, defined as the composition

    laplacian(t) = divergence(grad(t))

and therefore preserving the order of `t`.

If `CS` has been made the default with [`@set_coorsys`](@ref), the second
argument may be omitted.

# Examples
```julia
julia> Polar = coorsys_polar() ; r, θ = getcoords(Polar) ; @set_coorsys Polar

julia> laplacian(SymFunction("f", real = true)(r, θ))   # the polar Laplacian

julia> n = symbols("n", integer = true) ; simplify(laplacian(r^n * cos(n*θ)))   # = 0
```

See also [`grad`](@ref), [`divergence`](@ref), [`hess`](@ref).
"""
laplacian(t::AbstractTens{order, dim}, CS::AbstractCoorSystem{dim}) where {order, dim} =
    divergence(grad(t, CS), CS)
laplacian(t::SymType, CS::AbstractCoorSystem) = divergence(grad(t, CS), CS)

# Two methods, not one union. `Union{SymType, AbstractTens{order, dim}}` leaves
# `order` and `dim` bound by nothing when the argument is a `Num` or a `Sym`,
# and spelling the union the other way round, `Union{T, AbstractTens}`, leaves
# `T` unbound when it is a tensor. Split, each signature binds what it names.
"""
    hess(t::Union{T,AbstractTens{order,dim,T}}, CS::CoordinateSystem{dim,T,<:ExactGeometry}) where {order,dim,T<:SymType}

Hessian of `t`, defined as the second gradient

    hess(t) = grad(grad(t))

and therefore raising the order of `t` by two. For a scalar field the trace of
`hess(t)` is [`laplacian`](@ref)`(t)`.

If `CS` has been made the default with [`@set_coorsys`](@ref), the second
argument may be omitted.

# Examples
```julia
julia> Spherical = coorsys_spherical() ; θ, ϕ, r = getcoords(Spherical) ;

julia> @set_coorsys Spherical

julia> simplify(hess(1/r))    # the kernel of the 3-D Laplace equation
```

See also [`grad`](@ref), [`laplacian`](@ref).
"""
hess(t::AbstractTens{order, dim}, CS::AbstractCoorSystem{dim}) where {order, dim} =
    grad(grad(t, CS), CS)
hess(t::SymType, CS::AbstractCoorSystem) = grad(grad(t, CS), CS)

# `Tuple{T, Vararg{T, M}}`: at least one coordinate, which also binds `T` —
# `NTuple{dim, T}` accepts the empty tuple and leaves it bound by nothing.
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

julia> get_array(divergence(𝛔, Cartesian))
 Derivative(σ11(x, y, z), x) + Derivative(σ21(x, y, z), y) + Derivative(σ31(x, y, z), z)
 Derivative(σ21(x, y, z), x) + Derivative(σ22(x, y, z), y) + Derivative(σ32(x, y, z), z)
 Derivative(σ31(x, y, z), x) + Derivative(σ32(x, y, z), y) + Derivative(σ33(x, y, z), z)
```

which is `divergence(𝛔)ᵢ = ∂ⱼσᵢⱼ`.

See also [`coorsys_polar`](@ref), [`coorsys_cylindrical`](@ref),
[`coorsys_spherical`](@ref), [`@set_coorsys`](@ref).
"""
function coorsys_cartesian(coords::Tuple{T, Vararg{T, M}}) where {T <: SymType, M}
    dim = M + 1
    𝐗, 𝐄, ℬ = init_cartesian(coords)
    OM = sum([𝐗[i] * 𝐄[i] for i in 1:dim])
    χᵢ = ntuple(_ -> one(eltype(coords)), dim)
    return CoordinateSystem(OM, coords, ℬ, χᵢ)
end
coorsys_cartesian(::Val{<:Sym}, coords) = coorsys_cartesian(coords)
coorsys_cartesian(::Val{Num}, coords = Tuple(@variables x y z)) = coorsys_cartesian(coords)


"""
    coorsys_polar(coords = (symbols("r", positive = true), symbols("θ", real = true)); canonical = false)

Return the polar coordinate system

# Examples
```julia
julia> Polar = coorsys_polar() ; r, θ = getcoords(Polar) ; 𝐞ʳ, 𝐞ᶿ = unitvec(Polar) ; ℬᵖ = get_basis(Polar)

julia> f = SymFunction("f", real = true)(r, θ) ;

julia> laplacian(f, Polar)
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
function coorsys_polar(coords::NTuple{2, T}; canonical = false) where {T <: SymType}
    (r, θ), (𝐞ʳ, 𝐞ᶿ), ℬᵖ = init_polar(coords, canonical = canonical)
    OM = r * 𝐞ʳ
    return CoordinateSystem(OM, coords, ℬᵖ, (one(eltype(coords)), r))
end
coorsys_polar(::Val{<:Sym}, coords; canonical = false) = coorsys_polar(coords; canonical = canonical)
coorsys_polar(::Val{Num}, coords = Tuple(@variables r θ); canonical = false) = coorsys_polar(coords; canonical = canonical)


"""
    coorsys_cylindrical(coords = (symbols("r", positive = true), symbols("θ", real = true), symbols("z", real = true)); canonical = false)

Return the cylindrical coordinate system

# Examples
```julia
julia> Cylindrical = coorsys_cylindrical() ; rθz = getcoords(Cylindrical) ; 𝐞ʳ, 𝐞ᶿ, 𝐞ᶻ = unitvec(Cylindrical) ; ℬᶜ = get_basis(Cylindrical)

julia> 𝐯 = Tens(Vec{3}(i -> SymFunction("v\$(rθz[i])", real = true)(rθz...)), ℬᶜ) ;

julia> divergence(𝐯, Cylindrical)
                                                  ∂
                                    vr(r, θ, z) + ──(vθ(r, θ, z))
∂                 ∂                               ∂θ
──(vr(r, θ, z)) + ──(vz(r, θ, z)) + ─────────────────────────────
∂r                ∂z                              r
``` 
"""
function coorsys_cylindrical(coords::NTuple{3, T}; canonical = false) where {T <: SymType}
    (r, θ, z), (𝐞ʳ, 𝐞ᶿ, 𝐞ᶻ), ℬᶜ = init_cylindrical(coords, canonical = canonical)
    OM = r * 𝐞ʳ + z * 𝐞ᶻ

    return CoordinateSystem(OM, coords, ℬᶜ, (one(eltype(coords)), r, one(eltype(coords))))
end
coorsys_cylindrical(::Val{<:Sym}, coords; canonical = false) = coorsys_cylindrical(coords; canonical = canonical)
coorsys_cylindrical(::Val{Num}, coords = Tuple(@variables r θ z); canonical = false) = coorsys_cylindrical(coords; canonical = canonical)


"""
    coorsys_spherical(coords = (symbols("θ", real = true), symbols("ϕ", real = true), symbols("r", positive = true)); canonical = false)

Return the spherical coordinate system

# Examples
```julia
julia> Spherical = coorsys_spherical() ; θ, ϕ, r = getcoords(Spherical) ; 𝐞ᶿ, 𝐞ᵠ, 𝐞ʳ = unitvec(Spherical) ; ℬˢ = get_basis(Spherical)

julia> for σⁱʲ ∈ ("σʳʳ", "σᶿᶿ", "σᵠᵠ") @eval \$(Symbol(σⁱʲ)) = SymFunction(\$σⁱʲ, real = true)(\$r) end ;

julia> 𝛔 = σʳʳ * 𝐞ʳ ⊗ 𝐞ʳ + σᶿᶿ * 𝐞ᶿ ⊗ 𝐞ᶿ + σᵠᵠ * 𝐞ᵠ ⊗ 𝐞ᵠ ;

julia> div𝛔 = divergence(𝛔, Spherical)
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
function coorsys_spherical(coords::NTuple{3, T}; canonical = false) where {T <: SymType}
    (θ, ϕ, r), (𝐞ᶿ, 𝐞ᵠ, 𝐞ʳ), ℬˢ = init_spherical(coords, canonical = canonical)
    OM = r * 𝐞ʳ
    rules = Dict(abs(sin(θ)) => sin(θ), transpose(tan(θ)) => tan(θ), 1 // 1 => 1, 2 // 1 => 2)
    return CoordinateSystem(OM, coords, ℬˢ, (r, r * sin(θ), one(eltype(coords))); rules = rules)
end
coorsys_spherical(::Val{<:Sym}, coords; canonical = false) = coorsys_spherical(coords; canonical = canonical)
coorsys_spherical(::Val{Num}, coords = Tuple(@variables θ ϕ r); canonical = false) = coorsys_spherical(coords; canonical = canonical)

"""
    coorsys_spheroidal(coords = (symbols("ϕ", real = true),symbols("p", real = true),symbols("q", positive = true),),
                            c = symbols("c", positive = true),tmp_coords = (symbols("p̄ q̄", positive = true)...,),)

Return the spheroidal coordinate system

# Examples
```julia
julia> Spheroidal = coorsys_spheroidal() ; OM = position_vector(Spheroidal)
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

julia> laplacian(OM[1]^2, Spheroidal)
2
``` 
"""
function coorsys_spheroidal(coords::NTuple{3, T}, c, tmp_coords) where {T <: SymType}
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
    return CoordinateSystem(
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
coorsys_spheroidal(::Val{<:Sym}, coords, c, tmp_coords) = coorsys_spheroidal(coords, c, tmp_coords)
coorsys_spheroidal(
    ::Val{Num}, coords = Tuple(@variables ϕ p q),
    c = (@variables c)[1], tmp_coords = Tuple(@variables p̄ q̄)
) = coorsys_spheroidal(coords, c, tmp_coords)


for cs in (:coorsys_cartesian, :coorsys_polar, :coorsys_cylindrical, :coorsys_spherical, :coorsys_spheroidal)
    @eval $cs(T::Type, args...; kwargs...) = $cs(Val(T), args...; kwargs...)
end


# ── The default coordinate system ────────────────────────────────────────────
#
# The default chart is data held in a reference, and the single-argument forms
# of `grad`, `divergence`, … are ordinary methods defined once, below. Methods
# evaluated into the module by `@set_coorsys` would mutate the method table
# globally — invalidating compiled code, invisibly at the call site — and a
# default-chart `∂(t, x)` would be more specific than the plain derivative with
# respect to the symbol `x` and silently turn it into a covariant derivative.
# `∂` therefore has **no** default-chart method: it keeps its one meaning.

# One reference holds the whole default — the chart, the letter of its vectors
# and the names of its coordinates — so that a reader never sees parts of two
# different `set_coorsys!` calls.
const _DEFAULT_COORSYS = Ref{Any}(nothing)      # `nothing` or `(; CS, vec, coords)`

"""
    set_coorsys!(CS; vec = '𝐞', coords = nothing) → CS

Make `CS` the default coordinate system, so that [`grad`](@ref),
[`symgrad`](@ref), [`divergence`](@ref), [`laplacian`](@ref), [`hess`](@ref) and
[`pprint`](@ref) may be called with a single argument.

`vec` and `coords` control how [`pprint`](@ref) names the basis vectors.
The macro [`@set_coorsys`](@ref) is a thin wrapper over this function.

!!! note "`∂` is deliberately not affected"
    `∂(t, x)` always means the plain derivative of `t` with respect to the
    symbol `x`, whether or not a default chart is set. For the covariant
    derivative, pass the system explicitly: `∂(t, x, CS)`.

See also [`default_coorsys`](@ref), [`unset_coorsys!`](@ref).
"""
function set_coorsys!(CS::AbstractCoorSystem; vec::Char = '𝐞', coords = nothing)
    _DEFAULT_COORSYS[] = (; CS, vec, coords)
    return CS
end

"""
    default_coorsys() → AbstractCoorSystem

The coordinate system installed by [`set_coorsys!`](@ref) or
[`@set_coorsys`](@ref). Throws if none has been set.
"""
function default_coorsys()
    d = _DEFAULT_COORSYS[]
    d === nothing && error(
        "no default coordinate system is set: call `set_coorsys!(CS)` " *
            "(or `@set_coorsys CS`), or pass the system explicitly as the last argument."
    )
    return d.CS
end

"""
    unset_coorsys!()

Forget the default coordinate system installed by [`set_coorsys!`](@ref), with
the letter and the coordinate names it gave [`pprint`](@ref).
"""
function unset_coorsys!()
    _DEFAULT_COORSYS[] = nothing
    return nothing
end

# Single-argument forms, defined once and for all.
#
# One method per kind of argument: a single `Union{T, AbstractTens} where
# {T <: SymType}` would bind `T` to nothing for an `AbstractTens`, a signature
# Julia accepts and can never dispatch on.
for OP in (:grad, :symgrad, :laplacian, :hess)
    @eval $OP(t::AbstractTens) = $OP(t, default_coorsys())
    @eval $OP(t::T) where {T <: SymType} = $OP(t, default_coorsys())
end
divergence(t::AbstractTens) = divergence(t, default_coorsys())

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

julia> laplacian(1 / r)      # one argument
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

# One-argument `pprint` uses the default chart when one is set, so that a
# result prints with the chart's own coordinate names (`𝐞ʳ`, `𝐞ᶿ`, …) rather
# than with numbered basis vectors. Without a default it falls through to the
# generic method in `tens.jl`.
function pprint(t::AbstractTens{order, dim, T}) where {order, dim, T <: SymType}
    d = _DEFAULT_COORSYS[]
    d === nothing && return pprint(t; vec = '𝐞', coords = ntuple(i -> i, dim))
    get_dim(d.CS) != dim && return pprint(t; vec = d.vec, coords = ntuple(i -> i, dim))
    coords = d.coords === nothing ? string.(getcoords(d.CS)) : d.coords
    length(coords) == dim - 1 && (coords = (coords..., dim))
    return pprint(change_tens(t, normalized_basis(d.CS)); vec = d.vec, coords = coords)
end

function pprint(t::AbstractTens{order, dim, T}, CS::AbstractCoorSystem; vec = '𝐞') where {order, dim, T}
    coords = string.(getcoords(CS))
    ℬ = normalized_basis(CS)
    return pprint(change_tens(t, ℬ); vec = vec, coords = coords)
end

export ∂, AbstractGeometry, ExactGeometry, PointwiseGeometry, CoordinateSystem, lame_coefficients, christoffel, nderiv
export grad, symgrad, divergence, laplacian, hess
export normalized_basis, natural_basis, natvec, unitvec, getcoords, position_vector
export coorsys_cartesian, coorsys_polar, coorsys_cylindrical, coorsys_spherical, coorsys_spheroidal
export @set_coorsys, set_coorsys!, default_coorsys, unset_coorsys!
