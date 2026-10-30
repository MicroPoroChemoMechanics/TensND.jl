# ── Isotropic elastic moduli ───────────────────────────────────────────────────
#
# A 3-D isotropic stiffness is ℂ = 3k 𝕁 + 2μ 𝕂, and its compliance
# 𝕊 = 𝕁/(3k) + 𝕂/(2μ): the two coefficients on the projectors, `α` and `β`,
# are (3k, 2μ) for a stiffness and (1/(3k), 1/(2μ)) for a compliance.
#
# Six moduli describe the material: the bulk modulus k, Young's modulus E,
# Lamé's first parameter λ, the shear modulus μ, Poisson's ratio ν and the
# P-wave (oedometric) modulus M = λ + 2μ. Any two of them give the others: the
# stiffness or the compliance is built from any pair
# (`iso_stiffness_young_lame(E, λ)`, `iso_compliance_young_poisson(E, ν)`, …),
# and any modulus or pair is read from either (`lame_lambda(ℂ)`,
# `young_shear(𝕊; role = :compliance)`, …), each by its own formula on the two
# coefficients, never through an inverse.
#
# In 2-D, `TensISO{4,2}` is the plane-strain stiffness of a 3-D material,
# ℂ = 2k 𝕁 + 2μ 𝕂, or its compliance: `k` is the in-plane (area) bulk modulus
# λ + μ, and E, λ, μ, ν and M are those of the 3-D material. The builders take
# the dimension as a leading `Val(2)` (3 by default).
#
# No formula branches on a value, so that symbolic moduli go through, except
# for the two pairs whose conversion solves a quadratic, which take the root of
# the physical branch (see their docstrings).

# The role of an isotropic tensor, read from the keyword: the readers take
# `role = :stiffness` (the default) or `role = :compliance`. Propagated as a
# constant, so that the role written at the call is resolved at compile time.
Base.@constprop :aggressive function _role(role::Symbol)
    role === :stiffness && return Val(:stiffness)
    role === :compliance && return Val(:compliance)
    throw(ArgumentError("unknown role `:$role`; use `:stiffness` or `:compliance`"))
end

"""
    iso_stiffness(k, μ) → TensISO{4,3}
    iso_stiffness(Val(dim), k, μ) → TensISO{4,dim}

Isotropic stiffness ``\\mathbb{C} = d\\,k\\,\\mathbb{J} + 2\\mu\\,\\mathbb{K}`` from the bulk
modulus `k` and the shear modulus `μ`, in dimension ``d`` (3 by default; in 2-D,
`k` is the in-plane bulk modulus ``\\lambda + \\mu`` of the plane-strain stiffness); the
reciprocal of [`bulk_shear`](@ref). The compliance is [`iso_compliance`](@ref)`(k, μ)`.

`iso_stiffness_bulk_shear`, `iso_stiffness_k_mu` and `iso_stiffness_k_μ` are the
same function. The stiffness from any other pair of moduli is
`iso_stiffness_<pair>`: see [`iso_stiffness_young_poisson`](@ref).
"""
iso_stiffness(k, μ) = TensISO{3}(3 * k, 2 * μ)
iso_stiffness(::Val{3}, k, μ) = iso_stiffness(k, μ)
iso_stiffness(::Val{2}, k, μ) = TensISO{2}(2 * k, 2 * μ)

"""
    iso_compliance(k, μ) → TensISO{4,3}
    iso_compliance(Val(dim), k, μ) → TensISO{4,dim}

Isotropic compliance ``\\mathbb{S} = \\mathbb{J}/(d\\,k) + \\mathbb{K}/(2\\mu)`` from the
bulk modulus `k` (in 2-D, the in-plane ``\\lambda + \\mu``) and the shear modulus `μ`,
written directly, without an inverse: the inverse of [`iso_stiffness`](@ref)`(k, μ)`.

`iso_compliance_bulk_shear`, `iso_compliance_k_mu` and `iso_compliance_k_μ` are
the same function. The compliance from any other pair of moduli is
`iso_compliance_<pair>`: see [`iso_stiffness_young_poisson`](@ref) for the table.
"""
iso_compliance(k, μ) = iso_compliance(Val(3), k, μ)
iso_compliance(::Val{dim}, k, μ) where {dim} = TensISO{dim}(inv(dim * k), inv(2 * μ))

# The stiffness and the compliance of the 3-D material (λ, μ), or of its
# plane-strain restriction.
_iso_lame_shear(::Val{dim}, λ, μ) where {dim} = TensISO{dim}(dim * λ + 2 * μ, 2 * μ)
_iso_compliance_lame_shear(::Val{dim}, λ, μ) where {dim} = TensISO{dim}(inv(dim * λ + 2 * μ), inv(2 * μ))

"""
    iso_stiffness_young_poisson(E, ν) → TensISO{4,3}
    iso_stiffness_young_poisson(Val(dim), E, ν) → TensISO{4,dim}

Isotropic stiffness from Young's modulus and Poisson's ratio,
``k = E/(3(1 - 2\\nu))`` and ``\\mu = E/(2(1 + \\nu))``; the reciprocal of
[`young_poisson`](@ref). With `Val(2)`, the plane-strain stiffness of that
material.

The stiffness is built in the same way from any pair of the six moduli, named
in the order bulk `k`, Young `E`, Lamé `λ`, shear `μ`, Poisson `ν`, P-wave `M`,
and the compliance likewise by `iso_compliance_<pair>`, written directly
(``\\mathbb{S} = (1 - 2\\nu)/E\\;\\mathbb{J} + (1 + \\nu)/E\\;\\mathbb{K}`` from ``(E, \\nu)``):

| pair | stiffness | compliance | aliases (`iso_stiffness_` or `iso_compliance_` then) |
| :--- | :--- | :--- | :--- |
| ``(k, E)`` | `iso_stiffness_bulk_young` | `iso_compliance_bulk_young` | `k_E` |
| ``(k, \\lambda)`` | `iso_stiffness_bulk_lame` | `iso_compliance_bulk_lame` | `k_lambda`, `k_λ` |
| ``(k, \\mu)`` | [`iso_stiffness`](@ref), `iso_stiffness_bulk_shear` | [`iso_compliance`](@ref), `iso_compliance_bulk_shear` | `k_mu`, `k_μ` |
| ``(k, \\nu)`` | `iso_stiffness_bulk_poisson` | `iso_compliance_bulk_poisson` | `k_nu`, `k_ν` |
| ``(k, M)`` | `iso_stiffness_bulk_pwave` | `iso_compliance_bulk_pwave` | |
| ``(E, \\lambda)`` | `iso_stiffness_young_lame` | `iso_compliance_young_lame` | `E_lambda`, `E_λ` |
| ``(E, \\mu)`` | `iso_stiffness_young_shear` | `iso_compliance_young_shear` | `E_mu`, `E_μ` |
| ``(E, \\nu)`` | `iso_stiffness_young_poisson` | `iso_compliance_young_poisson` | `E_nu`, `E_ν` |
| ``(E, M)`` | `iso_stiffness_young_pwave` | `iso_compliance_young_pwave` | |
| ``(\\lambda, \\mu)`` | `iso_stiffness_lame_shear` | `iso_compliance_lame_shear` | `lambda_mu`, `λ_μ` |
| ``(\\lambda, \\nu)`` | `iso_stiffness_lame_poisson` | `iso_compliance_lame_poisson` | `lambda_nu`, `λ_ν` |
| ``(\\lambda, M)`` | `iso_stiffness_lame_pwave` | `iso_compliance_lame_pwave` | |
| ``(\\mu, \\nu)`` | `iso_stiffness_shear_poisson` | `iso_compliance_shear_poisson` | `mu_nu`, `μ_ν` |
| ``(\\mu, M)`` | `iso_stiffness_shear_pwave` | `iso_compliance_shear_pwave` | |
| ``(\\nu, M)`` | `iso_stiffness_poisson_pwave` | `iso_compliance_poisson_pwave` | |

The pair ``(\\lambda, \\nu)`` does not determine ``\\mu`` at ``\\nu = 0``, where
``\\lambda = 0`` too. Each takes the dimension as an optional leading `Val(2)` or
`Val(3)`. The pairs with the P-wave modulus have no short alias: `k_M` or `μ_M`
already name the moduli of a matrix in MeanFieldHomogenization.

Two conversions solve a quadratic. From ``(E, \\lambda)``,
``\\mu = (E - 3\\lambda + R)/4`` with ``R = \\sqrt{E^2 + 9\\lambda^2 + 2E\\lambda}``, the
only root with ``-1 < \\nu < 1/2``. From ``(E, M)``, two materials share the pair,
one with ``\\nu \\geq 0`` and one with ``\\nu \\leq 0``; the one kept is ``\\nu \\geq 0``,
``\\lambda = (M - E + S)/4`` and ``\\mu = (3M + E - S)/8`` with
``S = \\sqrt{E^2 + 9M^2 - 10EM}``. In 2-D, ``(k, E)`` gives
``\\mu = (3k - \\sqrt{9k^2 - 4Ek})/2``, the root of ``\\nu \\geq -1/4``: the other,
``-1/2 - \\nu``, has the same pair of moduli.
"""
iso_stiffness_young_poisson(E, ν) = iso_stiffness(E / (3 * (1 - 2 * ν)), E / (2 * (1 + ν)))
iso_stiffness_young_poisson(::Val{3}, E, ν) = iso_stiffness_young_poisson(E, ν)
iso_stiffness_young_poisson(v::Val{2}, E, ν) = _iso_lame_shear(v, _lame_shear(Val((:young, :poisson)), E, ν)...)

iso_compliance_young_poisson(E, ν) = iso_compliance_young_poisson(Val(3), E, ν)
iso_compliance_young_poisson(::Val{3}, E, ν) = TensISO{3}((1 - 2 * ν) / E, (1 + ν) / E)
iso_compliance_young_poisson(v::Val{2}, E, ν) = _iso_compliance_lame_shear(v, _lame_shear(Val((:young, :poisson)), E, ν)...)

# (λ, μ) of the 3-D material from a pair of moduli without `k`.
function _lame_shear(::Val{(:young, :lame)}, E, λ)
    R = sqrt(E^2 + 9 * λ^2 + 2 * E * λ)
    return λ, (E - 3 * λ + R) / 4
end
_lame_shear(::Val{(:young, :shear)}, E, μ) = μ * (E - 2 * μ) / (3 * μ - E), μ
_lame_shear(::Val{(:young, :poisson)}, E, ν) = E * ν / ((1 + ν) * (1 - 2 * ν)), E / (2 * (1 + ν))
function _lame_shear(::Val{(:young, :pwave)}, E, M)
    S = sqrt(E^2 + 9 * M^2 - 10 * E * M)
    return (M - E + S) / 4, (3 * M + E - S) / 8
end
_lame_shear(::Val{(:lame, :shear)}, λ, μ) = λ, μ
_lame_shear(::Val{(:lame, :poisson)}, λ, ν) = λ, λ * (1 - 2 * ν) / (2 * ν)
_lame_shear(::Val{(:lame, :pwave)}, λ, M) = λ, (M - λ) / 2
_lame_shear(::Val{(:shear, :poisson)}, μ, ν) = 2 * μ * ν / (1 - 2 * ν), μ
_lame_shear(::Val{(:shear, :pwave)}, μ, M) = M - 2 * μ, μ
_lame_shear(::Val{(:poisson, :pwave)}, ν, M) = M * ν / (1 - ν), M * (1 - 2 * ν) / (2 * (1 - ν))

# μ from the bulk modulus `k` of the dimension (3-D, or in-plane λ + μ in 2-D)
# and a second modulus.
_shear_from_bulk(::Val{3}, ::Val{:young}, k, E) = 3 * k * E / (9 * k - E)
_shear_from_bulk(::Val{2}, ::Val{:young}, k, E) = (3 * k - sqrt(9 * k^2 - 4 * E * k)) / 2
_shear_from_bulk(::Val{3}, ::Val{:lame}, k, λ) = 3 * (k - λ) / 2
_shear_from_bulk(::Val{2}, ::Val{:lame}, k, λ) = k - λ
_shear_from_bulk(::Val{3}, ::Val{:poisson}, k, ν) = 3 * k * (1 - 2 * ν) / (2 * (1 + ν))
_shear_from_bulk(::Val{2}, ::Val{:poisson}, k, ν) = k * (1 - 2 * ν)
_shear_from_bulk(::Val{3}, ::Val{:pwave}, k, M) = 3 * (M - k) / 4
_shear_from_bulk(::Val{2}, ::Val{:pwave}, k, M) = M - k

# The six moduli, in the order of the pairs, with their symbols in ASCII and
# in Unicode; the P-wave modulus has none, see above.
const _ISO_MODULI = (
    (:bulk, (:k, :k)), (:young, (:E, :E)), (:lame, (:lambda, :λ)),
    (:shear, (:mu, :μ)), (:poisson, (:nu, :ν)), (:pwave, nothing),
)

# The builders of the pairs that are not written above, stiffness and compliance.
for i in eachindex(_ISO_MODULI), j in (i + 1):length(_ISO_MODULI)
    p, q = first(_ISO_MODULI[i]), first(_ISO_MODULI[j])
    for (kind, base, fromλμ) in ((:stiffness, :iso_stiffness, :_iso_lame_shear), (:compliance, :iso_compliance, :_iso_compliance_lame_shear))
        (p, q) in ((:bulk, :shear), (:young, :poisson)) && continue
        f = Symbol(base, :_, p, :_, q)
        doc = "    $f(a, b) → TensISO{4,3}\n    $f(Val(dim), a, b) → TensISO{4,dim}\n\n" *
            "Isotropic $kind from its $(p) and $(q) moduli (Poisson's ratio for `poisson`, " *
            "Lamé's first parameter for `lame`), in 3-D or, with `Val(2)`, in plane strain. " *
            "See [`iso_stiffness_young_poisson`](@ref) for the table of the pairs and their conventions."
        if p === :bulk
            @eval begin
                @doc $doc $f(a, b) = $f(Val(3), a, b)
                $f(v::Val, k, b) = $base(v, k, _shear_from_bulk(v, Val($(QuoteNode(q))), k, b))
            end
        else
            @eval begin
                @doc $doc $f(a, b) = $f(Val(3), a, b)
                $f(v::Val, a, b) = $fromλμ(v, _lame_shear(Val(($(QuoteNode(p)), $(QuoteNode(q)))), a, b)...)
            end
        end
    end
end
const iso_stiffness_bulk_shear = iso_stiffness
const iso_compliance_bulk_shear = iso_compliance

# ── The moduli read from a stiffness or a compliance ───────────────────────────
#
# One kernel per modulus, role and dimension, on the coefficients `(α, β)` of
# 𝕁 and 𝕂: on a stiffness, `α = d k`, `β = 2μ`; on a compliance, `α = 1/(d k)`,
# `β = 1/(2μ)`. The formulas on a compliance are written directly, never through
# its inverse; those on a stiffness are those of 0.11.0, bit for bit.

_bulk(::Val{:stiffness}, ::Val{d}, α, β) where {d} = α / d
_bulk(::Val{:compliance}, ::Val{d}, α, β) where {d} = inv(d * α)
_shear(::Val{:stiffness}, ::Val, α, β) = β / 2
_shear(::Val{:compliance}, ::Val, α, β) = inv(2 * β)
_lame(::Val{:stiffness}, ::Val{d}, α, β) where {d} = (α - β) / d
_lame(::Val{:compliance}, ::Val{d}, α, β) where {d} = (β - α) / (d * α * β)
_pwave(r::Val, v::Val, α, β) = _lame(r, v, α, β) + 2 * _shear(r, v, α, β)
_pwave(::Val{:compliance}, ::Val{d}, α, β) where {d} = (β + (d - 1) * α) / (d * α * β)
_young(r::Val{:stiffness}, v::Val{3}, α, β) = (k = _bulk(r, v, α, β); μ = _shear(r, v, α, β); 9 * k * μ / (3 * k + μ))
_young(r::Val{:stiffness}, v::Val{2}, α, β) = 2 * _shear(r, v, α, β) * (1 + _poisson(r, v, α, β))
_young(::Val{:compliance}, ::Val{3}, α, β) = 3 / (α + 2 * β)
_young(::Val{:compliance}, ::Val{2}, α, β) = (3 * β - α) / (2 * β^2)
_poisson(r::Val{:stiffness}, v::Val{3}, α, β) = poisson_ratio(_bulk(r, v, α, β), _shear(r, v, α, β))
_poisson(::Val{:stiffness}, ::Val{2}, α, β) = (α - β) / (2 * α)
_poisson(::Val{:compliance}, ::Val{3}, α, β) = (β - α) / (α + 2 * β)
_poisson(::Val{:compliance}, ::Val{2}, α, β) = (β - α) / (2 * β)

"""
    poisson_ratio(k, μ) → ν

Poisson's ratio of a 3-D isotropic material from its bulk and shear moduli,
``\\nu = (3k - 2\\mu)/(2(3k + \\mu))``. On a tensor, see [`lame_lambda`](@ref).
"""
poisson_ratio(k, μ) = (3 * k - 2 * μ) / (2 * (3 * k + μ))

"""
    lame_lambda(t::TensISO{4,dim}; role = :stiffness)

Lamé's first parameter ``\\lambda`` of an isotropic stiffness
``\\mathbb{C} = (d\\lambda + 2\\mu)\\,\\mathbb{J} + 2\\mu\\,\\mathbb{K}``, or, with
`role = :compliance`, of the material whose compliance is `t`; in 2-D, that of
the 3-D material whose plane-strain stiffness (or compliance) is `t`.

Each modulus has its reader, which computes that modulus alone, from the two
coefficients of `t` on ``\\mathbb{J}`` and ``\\mathbb{K}``, by a formula of its own for
each role (never through an inverse): [`bulk_modulus`](@ref) (in 2-D, the
in-plane ``\\lambda + \\mu``), [`young_modulus`](@ref), `lame_lambda`,
[`shear_modulus`](@ref), [`poisson_ratio`](@ref), [`pwave_modulus`](@ref)
(``M = \\lambda + 2\\mu``). On a 3-D compliance ``\\mathbb{S} = a\\,\\mathbb{J} + b\\,\\mathbb{K}``,
for instance, ``E = 3/(a + 2b)`` and ``\\nu = (b - a)/(a + 2b)``.

Each pair has one too, named as its builder without `iso_stiffness_` (`bulk_shear`,
`young_poisson`, `lame_shear`, `young_pwave`, …), with the same aliases (`k_mu`,
`E_nu`, `λ_μ`, …), the same `role` keyword, and the two moduli of its readers:
the conversion of a pair of moduli to another is
`young_shear(iso_stiffness_lame_poisson(λ, ν))`.
"""
function lame_lambda end

"The bulk modulus of an isotropic stiffness, in 2-D the in-plane one ``\\lambda + \\mu``, or of a compliance with `role = :compliance`; see [`lame_lambda`](@ref)."
function bulk_modulus end
"The shear modulus of an isotropic stiffness, or of a compliance with `role = :compliance`; see [`lame_lambda`](@ref)."
function shear_modulus end
"Young's modulus of an isotropic stiffness, or of a compliance with `role = :compliance`; see [`lame_lambda`](@ref)."
function young_modulus end
"The P-wave (oedometric) modulus ``M = \\lambda + 2\\mu`` of an isotropic stiffness, or of a compliance with `role = :compliance`; see [`lame_lambda`](@ref)."
function pwave_modulus end

for (reader, kernel) in (
        (:bulk_modulus, :_bulk), (:shear_modulus, :_shear), (:lame_lambda, :_lame),
        (:young_modulus, :_young), (:poisson_ratio, :_poisson), (:pwave_modulus, :_pwave),
    )
    @eval Base.@constprop :aggressive function $reader(t::TensISO{4, dim}; role::Symbol = :stiffness) where {dim}
        α, β = get_data(t)
        return $kernel(_role(role), Val(dim), α, β)
    end
end

"""
    bulk_shear(t::TensISO{4,dim}; role = :stiffness) → (k, μ)

Bulk and shear moduli of an isotropic **stiffness** ``\\mathbb{C} = d\\,k\\,\\mathbb{J} +
2\\mu\\,\\mathbb{K}`` in dimension ``d``, or of a **compliance** with `role = :compliance`:
``k`` relates the trace of the stress to the trace of the strain,
``\\mathrm{tr}\\,\\boldsymbol{\\sigma} = d\\,k\\,\\mathrm{tr}\\,\\boldsymbol{\\varepsilon}``. In 3-D it
is the bulk modulus; in 2-D the in-plane (area) bulk modulus, ``\\lambda + \\mu`` for the
plane-strain restriction of a 3-D material. For a tensor that is not a `TensISO`,
take its isotropic part first: [`isotropify`](@ref) (the exact rotation average)
or [`best_fit_iso`](@ref) (the closest isotropic tensor).

`k_mu` and `k_μ` are the same function, under the symbols of the moduli.
"""
function bulk_shear end

"""
    young_poisson(t::TensISO{4,dim}; role = :stiffness) → (E, ν)

Young's modulus and Poisson's ratio of an isotropic **stiffness**, or of a
**compliance** with `role = :compliance`. In 3-D, ``E = 9k\\mu/(3k + \\mu)`` and
``\\nu`` = [`poisson_ratio`](@ref)`(k, μ)` on a stiffness, ``E = 3/(a + 2b)`` and
``\\nu = (b - a)/(a + 2b)`` on a compliance ``a\\,\\mathbb{J} + b\\,\\mathbb{K}``. In 2-D,
those of the 3-D material whose plane-strain stiffness (or compliance) is `t`:
on ``\\mathbb{C} = \\alpha\\,\\mathbb{J} + \\beta\\,\\mathbb{K}``, ``\\nu = (\\alpha - \\beta)/(2\\alpha)``
and ``E = 2\\mu(1 + \\nu)``, ``\\mu = \\beta/2``.

See also [`bulk_shear`](@ref), [`shear_poisson`](@ref), [`iso_stiffness_young_poisson`](@ref),
and [`lame_lambda`](@ref) for every other modulus and pair.

`E_nu` and `E_ν` are the same function, under the symbols of the moduli.
"""
function young_poisson end

"""
    shear_poisson(t::TensISO{4,dim}; role = :stiffness) → (μ, ν)

Shear modulus and Poisson's ratio of an isotropic **stiffness**, or of a
**compliance** with `role = :compliance`, the pair the Kelvin and plane-strain
Green functions are written with; in 2-D, ``\\nu`` is that of the 3-D material
whose plane-strain stiffness is `t` (see [`young_poisson`](@ref)).

`mu_nu` and `μ_ν` are the same function, under the symbols of the moduli.
"""
function shear_poisson end

# The readers of the pairs: the two moduli of their single readers.
const _ISO_READERS = Dict(
    :bulk => :bulk_modulus, :young => :young_modulus, :lame => :lame_lambda,
    :shear => :shear_modulus, :poisson => :poisson_ratio, :pwave => :pwave_modulus,
)
for i in eachindex(_ISO_MODULI), j in (i + 1):length(_ISO_MODULI)
    p, q = first(_ISO_MODULI[i]), first(_ISO_MODULI[j])
    f, rp, rq = Symbol(p, :_, q), _ISO_READERS[p], _ISO_READERS[q]
    @eval Base.@constprop :aggressive $f(t::TensISO{4}; role::Symbol = :stiffness) = ($rp(t; role), $rq(t; role))
    (p, q) in ((:bulk, :shear), (:young, :poisson), (:shear, :poisson)) && continue
    doc = "    $f(t::TensISO{4}; role = :stiffness) → ($(p), $(q))\n\nThe $(p) and $(q) moduli of an " *
        "isotropic stiffness, or of a compliance with `role = :compliance` (`$rp(t; role)`, " *
        "`$rq(t; role)`); see [`lame_lambda`](@ref) for the conventions."
    @eval @doc $doc $f(::TensISO{4})
end


# The symbols: exact aliases of the functions named in words, shorter and closer
# to the formulas, in ASCII and with Greek letters.
for i in eachindex(_ISO_MODULI), j in (i + 1):length(_ISO_MODULI)
    (p, ps), (q, qs) = _ISO_MODULI[i], _ISO_MODULI[j]
    symbols = ps === nothing || qs === nothing ? () : unique(((ps[1], qs[1]), (ps[2], qs[2])))
    for (a, b) in symbols
        @eval begin
            const $(Symbol(a, :_, b)) = $(Symbol(p, :_, q))
            const $(Symbol(:iso_stiffness_, a, :_, b)) = $(Symbol(:iso_stiffness_, p, :_, q))
            const $(Symbol(:iso_compliance_, a, :_, b)) = $(Symbol(:iso_compliance_, p, :_, q))
            export $(Symbol(a, :_, b)), $(Symbol(:iso_stiffness_, a, :_, b)), $(Symbol(:iso_compliance_, a, :_, b))
        end
    end
    @eval export $(Symbol(p, :_, q)), $(Symbol(:iso_stiffness_, p, :_, q)), $(Symbol(:iso_compliance_, p, :_, q))
end

export iso_stiffness, iso_compliance, poisson_ratio, lame_lambda, bulk_modulus, shear_modulus, young_modulus, pwave_modulus

# ── Plane strain and plane stress ─────────────────────────────────────────────

"""
    plane_strain_stiffness(ℂ::TensISO{4,3}) → TensISO{4,2}

The in-plane stiffness of an isotropic material under **plane strain**
(``\\varepsilon_{i3} = 0``): the in-plane block of ``\\mathbb{C}``, with the same ``\\lambda``
and ``\\mu``, ``2(\\lambda + \\mu)\\,\\mathbb{J} + 2\\mu\\,\\mathbb{K}`` in two dimensions. It is
the convention of every `TensISO{4,2}` here: its [`young_poisson`](@ref) are those
of `ℂ`.

See also [`plane_stress_stiffness`](@ref).
"""
plane_strain_stiffness(ℂ::TensISO{4, 3}) = _iso_lame_shear(Val(2), lame_lambda(ℂ), shear_modulus(ℂ))

"""
    plane_stress_stiffness(ℂ::TensISO{4,3}) → TensISO{4,2}

The in-plane stiffness of an isotropic material under **plane stress**
(``\\sigma_{i3} = 0``): ``\\varepsilon_{33}`` condensed out, which replaces ``\\lambda``
by ``\\lambda^* = 2\\lambda\\mu/(\\lambda + 2\\mu) = E\\nu/(1 - \\nu^2)`` and keeps ``\\mu``.
Read as a plane-strain stiffness, as every `TensISO{4,2}` is, its moduli are
[`plane_stress_to_plane_strain`](@ref)`(E, ν)`.

See also [`plane_strain_stiffness`](@ref).
"""
function plane_stress_stiffness(ℂ::TensISO{4, 3})
    λ, μ = lame_shear(ℂ)
    return _iso_lame_shear(Val(2), 2 * λ * μ / (λ + 2 * μ), μ)
end

"""
    plane_stress_to_plane_strain(E, ν) → (E′, ν′)

The moduli with which a plane-strain computation reproduces a plane-stress
state of the material ``(E, \\nu)``: ``E' = E(1 + 2\\nu)/(1 + \\nu)^2``,
``\\nu' = \\nu/(1 + \\nu)``. The plane-strain stiffness of ``(E', \\nu')`` is the
plane-stress stiffness of ``(E, \\nu)``; [`plane_strain_to_plane_stress`](@ref) is
the reciprocal.
"""
plane_stress_to_plane_strain(E, ν) = (E * (1 + 2 * ν) / (1 + ν)^2, ν / (1 + ν))

"""
    plane_strain_to_plane_stress(E, ν) → (E′, ν′)

The moduli with which a plane-stress computation reproduces a plane-strain
state of the material ``(E, \\nu)``: ``E' = E/(1 - \\nu^2)``, ``\\nu' = \\nu/(1 - \\nu)``;
the reciprocal of [`plane_stress_to_plane_strain`](@ref).
"""
plane_strain_to_plane_stress(E, ν) = (E / (1 - ν^2), ν / (1 - ν))

export plane_strain_stiffness, plane_stress_stiffness, plane_stress_to_plane_strain, plane_strain_to_plane_stress
