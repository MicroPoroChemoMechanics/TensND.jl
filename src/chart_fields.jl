# ── Tensor fields of constant components in the local frame of a chart ───────
#
# A material deposited layer by layer on a sphere, wound around an axis or
# grown along the radius has the same tensor of properties at every point *in
# the frame that turns with the point*: the unit vectors of a curvilinear
# system. `ChartTens` holds such a field: a pattern, whose canonical components
# are those of the field in the local frame, and the chart whose unit vectors
# make that frame. At a point it is an ordinary tensor (`tensor_at`), with
# numbers, dual numbers or symbols; in an exact coordinate system it is the
# tensor field on which the differential operators act.
#
# The frame is that of `unit_vectors(chart, ξ)`, in the order of the coordinates
# of the chart: (e_r, e_θ, e_z) for the cylindrical one, (e_θ, e_ϕ, e_r) for the
# spherical one.

"""
    ChartTens(pattern::AbstractTens, chart::AbstractChart)

The tensor field whose components in the local frame of `chart` (its unit
vectors, in the order of its coordinates) are, at every point, the components
of `pattern` in the canonical basis: a tensor that turns with the point. A
`pattern` transversely isotropic about ``\\underline{e}_3`` with a
[`SphericalChart`](@ref) is transversely isotropic about the radial direction
at every point; an orthotropic one with a [`CylindricalChart`](@ref) is the
cylindrical orthotropy of a wound material.

[`tensor_at`](@ref) gives the tensor at a point, or the field in an exact
coordinate system. Fields of the same chart add, and a field is multiplied by a
number and inverted, as their patterns are.

# Examples
```jldoctest
julia> F = ChartTens(TensTI{2}(1.0, 5.0, (0.0, 0.0, 1.0)), SphericalChart());   # 5 along the radius, 1 across

julia> K = tensor_at(F, (0.0, 2.0, 0.0));                                        # at a point of the axis 𝐞₂

julia> canonical_components(K) ≈ [1 0 0; 0 5 0; 0 0 1]
true
```

See also [`radial_ti`](@ref), [`radial_axis`](@ref).
"""
struct ChartTens{order, dim, T, P <: AbstractTens{order, dim, T}, C <: AbstractChart{dim}}
    pattern::P
    chart::C
    ChartTens(pattern::P, chart::C) where {order, dim, T, P <: AbstractTens{order, dim, T}, C <: AbstractChart{dim}} =
        new{order, dim, T, P, C}(pattern, chart)
end

Base.eltype(::Type{<:ChartTens{order, dim, T}}) where {order, dim, T} = T
Base.eltype(F::ChartTens) = eltype(typeof(F))
Base.:(==)(F::ChartTens, G::ChartTens) = F.chart == G.chart && F.pattern == G.pattern

function _same_chart(F::ChartTens, G::ChartTens)
    F.chart == G.chart || throw(
        ArgumentError("the two fields are written in the frames of different charts: $(F.chart) and $(G.chart)")
    )
    return F.chart
end
Base.:+(F::ChartTens, G::ChartTens) = ChartTens(F.pattern + G.pattern, _same_chart(F, G))
Base.:-(F::ChartTens, G::ChartTens) = ChartTens(F.pattern - G.pattern, _same_chart(F, G))
Base.:-(F::ChartTens) = ChartTens(-F.pattern, F.chart)
Base.:*(x::Number, F::ChartTens) = ChartTens(x * F.pattern, F.chart)
Base.:*(F::ChartTens, x::Number) = x * F
Base.:/(F::ChartTens, x::Number) = ChartTens(F.pattern / x, F.chart)
Base.inv(F::ChartTens) = ChartTens(inv(F.pattern), F.chart)

Base.show(io::IO, F::ChartTens) = print(io, "ChartTens(", F.pattern, ", ", F.chart, ")")

"""
    tensor_at(F::ChartTens, x) -> AbstractTens
    tensor_at(F::ChartTens, CS::AbstractCoorSystem) -> Tens

The tensor of the field `F` at the point of Cartesian components `x`: the
pattern, read in the frame of the unit vectors of the chart at that point. An
isotropic pattern is returned as it is, and a transversely isotropic one keeps
its class, its axis being carried by the frame. `x` may hold dual numbers; it is
not on the singular set of the chart (its axis, its center), where the frame is
not defined.

With a coordinate system `CS`, the field itself, as a tensor of constant
components in the normalized basis of `CS`, on which `symgrad`, `divergence`, …
act. `CS` must be the system of the chart of `F`
(`coorsys_spherical()` for a [`SphericalChart`](@ref), …): this is not checked.
"""
tensor_at(F::ChartTens, x) = _carried(F.pattern, unit_vectors(F.chart, cartesian_to_coords(F.chart, x)))
tensor_at(F::ChartTens, CS::AbstractCoorSystem) = Tens(canonical_components(F.pattern), normalized_basis(CS))

# The pattern carried into the frame whose vectors are the columns of `Q`.
_carried(t::TensISO, Q) = t
function _carried(t::TensTI{order, T, N}, Q) where {order, T, N}
    S = promote_type(T, eltype(Q))
    n = Q * SVector{3}(axis(t))
    return TensTI{order, S, N}(map(S, get_data(t)), (S(n[1]), S(n[2]), S(n[3])))
end
_carried(t::AbstractTens, Q) = Tens(canonical_components(t), RotatedBasis(Q))

"""
    radial_axis(chart::AbstractChart) -> NTuple

The components, in the local frame of `chart`, of its unit radial vector:
`(0, 0, 1)` for a [`SphericalChart`](@ref), whose coordinates are `(θ, ϕ, r)`,
and `(1, 0, 0)` for a [`CylindricalChart`](@ref), `(r, θ, z)`. It is the axis to
give to a pattern that is transversely isotropic about the radius.
"""
radial_axis(::SphericalChart) = (0, 0, 1)
radial_axis(::CylindricalChart) = (1, 0, 0)

"""
    radial_ti(chart, ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆) -> ChartTens
    radial_ti(chart, c₁, c₂, c₃, c₄, c₅; form, role = :compliance) -> ChartTens
    radial_ti(chart; n, l, k, m, p) -> ChartTens
    radial_ti(chart, aₜ, aᵣ) -> ChartTens
    radial_ti(chart, t::TensISO) -> ChartTens

A field transversely isotropic about the radial direction of `chart`
([`radial_axis`](@ref)), at every point.

At order 4, from its coefficients on the Walpole basis about the radius, or
from Hill's moduli about it: with ``r`` the radial direction and ``t, t'`` two
orthogonal directions across it, ``n = C_{rrrr}``, ``l = C_{rrtt}``,
``k = \\tfrac12(C_{tttt} + C_{ttt't'})``, ``m = \\tfrac12(C_{tttt} - C_{ttt't'})``
and ``p = C_{rtrt}``, that is ``\\ell_1 = n``, ``\\ell_2 = 2k``,
``\\ell_3 = \\sqrt2\\,l``, ``\\ell_5 = 2m``, ``\\ell_6 = 2p``. At order 2, from
the coefficients `aₜ` across the radius and `aᵣ` along it. An isotropic tensor
is written as the field it is in that class.

With `form`, the five constants are read as [`ti_from_constants`](@ref) reads
them, the radial direction being its axis 3 and the directions across the
radius its axes 1 and 2: `form = :components` takes
``(C_{tttt}, C_{ttt't'}, C_{ttrr}, C_{rrrr}, C_{trtr})``, `:hill` Hill's moduli
in the order ``(k, l, m, n, p)`` (the keywords of the third method, by
position), `:engineering` the
engineering constants ``(E_t, E_r, \\nu_{tt'}, \\nu_{rt}, G_{rt})`` and `:hoenig`
Hoenig's parameters, the last two as a compliance or, with `role = :stiffness`,
as the stiffness of the same material. `form = :walpole`, the default, is the
first method. [`hill_moduli`](@ref) reads Hill's moduli of the field back, and
`ti_constants(F.pattern; form, role)` its constants in another form.
"""
function radial_ti end

Base.@constprop :aggressive function radial_ti(
        chart::AbstractChart{3}, c₁, c₂, c₃, c₄, c₅; form::Symbol = :walpole, role::Symbol = :compliance
    )
    return ChartTens(ti_from_constants(c₁, c₂, c₃, c₄, c₅, radial_axis(chart); form, role), chart)
end
radial_ti(chart::AbstractChart{3}; n, l, k, m, p) = ChartTens(_ti_from_hill(k, l, m, n, p, radial_axis(chart)), chart)
radial_ti(chart::AbstractChart{3}, aₜ, aᵣ) = ChartTens(TensTI{2}(aₜ, aᵣ, radial_axis(chart)), chart)
radial_ti(chart::AbstractChart{3}, t::TensISO{4, 3}) = ChartTens(t + TensTI{4}(0, 0, 0, 0, 0, radial_axis(chart)), chart)
radial_ti(chart::AbstractChart{3}, t::TensISO{2, 3}) = (a = first(get_data(t)); radial_ti(chart, a, a))

hill_moduli(F::ChartTens{4}) = hill_moduli(F.pattern)

export ChartTens, tensor_at, radial_axis, radial_ti
