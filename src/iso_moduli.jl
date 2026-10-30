# ── Isotropic elastic moduli ───────────────────────────────────────────────────
#
# A 3-D isotropic stiffness is ℂ = 3k 𝕁 + 2μ 𝕂: its coefficients on the two
# projectors are 3k and 2μ. A compliance is the inverse, 𝕊 = 𝕁/(3k) + 𝕂/(2μ),
# and `inv` is exact on `TensISO`, so the moduli of a compliance are those of
# `inv(𝕊)`: no second family of functions.
#
# Six moduli describe the material: the bulk modulus k, Young's modulus E,
# Lamé's first parameter λ, the shear modulus μ, Poisson's ratio ν and the
# P-wave (oedometric) modulus M = λ + 2μ. Any two of them give the others: the
# stiffness is built from any pair (`iso_stiffness_young_lame(E, λ)`, …), and
# any modulus or pair is read from a stiffness (`lame_lambda(ℂ)`,
# `young_shear(ℂ)`, …), so that a conversion is a composition of the two.
#
# In 2-D, `TensISO{4,2}` is the plane-strain stiffness of a 3-D material,
# ℂ = 2k 𝕁 + 2μ 𝕂: `k` is the in-plane (area) bulk modulus λ + μ, as
# `bulk_shear` reads it, and E, λ, μ, ν and M are those of the 3-D material.
# The builders take the dimension as a leading `Val(2)` (3 by default).
#
# No formula branches on a value, so that symbolic moduli go through, except
# for the two pairs whose conversion solves a quadratic, which take the root of
# the physical branch (see their docstrings).

"""
    iso_stiffness(k, μ) → TensISO{4,3}
    iso_stiffness(Val(dim), k, μ) → TensISO{4,dim}

Isotropic stiffness ``\\mathbb{C} = d\\,k\\,\\mathbb{J} + 2\\mu\\,\\mathbb{K}`` from the bulk
modulus `k` and the shear modulus `μ`, in dimension ``d`` (3 by default; in 2-D,
`k` is the in-plane bulk modulus ``\\lambda + \\mu`` of the plane-strain stiffness); the
reciprocal of [`bulk_shear`](@ref). The compliance is `inv(iso_stiffness(k, μ))`.

`iso_stiffness_bulk_shear`, `iso_stiffness_k_mu` and `iso_stiffness_k_μ` are the
same function. The stiffness from any other pair of moduli is
`iso_stiffness_<pair>`: see [`iso_stiffness_young_poisson`](@ref).
"""
iso_stiffness(k, μ) = TensISO{3}(3 * k, 2 * μ)
iso_stiffness(::Val{3}, k, μ) = iso_stiffness(k, μ)
iso_stiffness(::Val{2}, k, μ) = TensISO{2}(2 * k, 2 * μ)

# The stiffness of the 3-D material (λ, μ), or its plane-strain restriction.
_iso_lame_shear(::Val{dim}, λ, μ) where {dim} = TensISO{dim}(dim * λ + 2 * μ, 2 * μ)

"""
    bulk_shear(ℂ::TensISO{4,dim}) → (k, μ)

Bulk and shear moduli of an isotropic **stiffness** ``\\mathbb{C} = d\\,k\\,\\mathbb{J} +
2\\mu\\,\\mathbb{K}`` in dimension ``d``: ``k`` relates the trace of the stress to the
trace of the strain, ``\\mathrm{tr}\\,\\boldsymbol{\\sigma} = d\\,k\\,\\mathrm{tr}\\,\\boldsymbol{\\varepsilon}``.
In 3-D it is the bulk modulus; in 2-D the in-plane (area) bulk modulus, ``\\lambda + \\mu``
for the plane-strain restriction of a 3-D material. For a compliance
``\\mathbb{S}``, use `bulk_shear(inv(𝕊))`. For a tensor that is not a `TensISO`, take its
isotropic part first: [`isotropify`](@ref) (the exact rotation average) or
[`best_fit_iso`](@ref) (the closest isotropic tensor).

`k_mu` and `k_μ` are the same function, under the symbols of the moduli.
"""
function bulk_shear(ℂ::TensISO{4, dim}) where {dim}
    α, β = get_data(ℂ)
    return α / dim, β / 2
end

"""
    poisson_ratio(k, μ) → ν
    poisson_ratio(ℂ::TensISO{4}) → ν

Poisson's ratio of a 3-D isotropic material from its bulk and shear moduli,
``\\nu = (3k - 2\\mu)/(2(3k + \\mu))``, or of an isotropic stiffness (in 2-D, that of
the 3-D material whose plane-strain stiffness it is; see [`young_poisson`](@ref)).
"""
poisson_ratio(k, μ) = (3 * k - 2 * μ) / (2 * (3 * k + μ))

"""
    young_poisson(ℂ::TensISO{4,3}) → (E, ν)
    young_poisson(ℂ::TensISO{4,2}) → (E, ν)

Young's modulus and Poisson's ratio of an isotropic **stiffness**. In 3-D,
``E = 9k\\mu/(3k + \\mu)`` and ``\\nu`` = [`poisson_ratio`](@ref)`(k, μ)`. In 2-D, those
of the 3-D material whose plane-strain stiffness is ``\\mathbb{C}``: with
``\\mathbb{C} = \\alpha\\,\\mathbb{J} + \\beta\\,\\mathbb{K}``, ``\\nu = (\\alpha - \\beta)/(2\\alpha)``
and ``E = 2\\mu(1 + \\nu)``, ``\\mu = \\beta/2``. For a compliance, use `young_poisson(inv(𝕊))`.

See also [`bulk_shear`](@ref), [`shear_poisson`](@ref), [`iso_stiffness_young_poisson`](@ref),
and [`lame_lambda`](@ref) for every other modulus and pair.

`E_nu` and `E_ν` are the same function, under the symbols of the moduli.
"""
function young_poisson(ℂ::TensISO{4, 3})
    k, μ = bulk_shear(ℂ)
    return 9 * k * μ / (3 * k + μ), poisson_ratio(k, μ)
end
function young_poisson(ℂ::TensISO{4, 2})
    μ, ν = shear_poisson(ℂ)
    return 2 * μ * (1 + ν), ν
end

"""
    shear_poisson(ℂ::TensISO{4,dim}) → (μ, ν)

Shear modulus and Poisson's ratio of an isotropic **stiffness**, the pair the
Kelvin and plane-strain Green functions are written with; in 2-D, ``\\nu`` is
that of the 3-D material whose plane-strain stiffness is ``\\mathbb{C}`` (see
[`young_poisson`](@ref)).

`mu_nu` and `μ_ν` are the same function, under the symbols of the moduli.
"""
shear_poisson(ℂ::TensISO{4, 3}) = (bulk_shear(ℂ)[2], young_poisson(ℂ)[2])
function shear_poisson(ℂ::TensISO{4, 2})
    α, β = get_data(ℂ)
    return β / 2, (α - β) / (2 * α)
end

"""
    iso_stiffness_young_poisson(E, ν) → TensISO{4,3}
    iso_stiffness_young_poisson(Val(dim), E, ν) → TensISO{4,dim}

Isotropic stiffness from Young's modulus and Poisson's ratio,
``k = E/(3(1 - 2\\nu))`` and ``\\mu = E/(2(1 + \\nu))``; the reciprocal of
[`young_poisson`](@ref). With `Val(2)`, the plane-strain stiffness of that
material.

The stiffness is built in the same way from any pair of the six moduli, named
in the order bulk `k`, Young `E`, Lamé `λ`, shear `μ`, Poisson `ν`, P-wave `M`:

| pair | builder | aliases |
| :--- | :--- | :--- |
| ``(k, E)`` | `iso_stiffness_bulk_young` | `iso_stiffness_k_E` |
| ``(k, \\lambda)`` | `iso_stiffness_bulk_lame` | `iso_stiffness_k_lambda`, `iso_stiffness_k_λ` |
| ``(k, \\mu)`` | [`iso_stiffness`](@ref), `iso_stiffness_bulk_shear` | `iso_stiffness_k_mu`, `iso_stiffness_k_μ` |
| ``(k, \\nu)`` | `iso_stiffness_bulk_poisson` | `iso_stiffness_k_nu`, `iso_stiffness_k_ν` |
| ``(k, M)`` | `iso_stiffness_bulk_pwave` | |
| ``(E, \\lambda)`` | `iso_stiffness_young_lame` | `iso_stiffness_E_lambda`, `iso_stiffness_E_λ` |
| ``(E, \\mu)`` | `iso_stiffness_young_shear` | `iso_stiffness_E_mu`, `iso_stiffness_E_μ` |
| ``(E, \\nu)`` | `iso_stiffness_young_poisson` | `iso_stiffness_E_nu`, `iso_stiffness_E_ν` |
| ``(E, M)`` | `iso_stiffness_young_pwave` | |
| ``(\\lambda, \\mu)`` | `iso_stiffness_lame_shear` | `iso_stiffness_lambda_mu`, `iso_stiffness_λ_μ` |
| ``(\\lambda, \\nu)`` | `iso_stiffness_lame_poisson` | `iso_stiffness_lambda_nu`, `iso_stiffness_λ_ν` |
| ``(\\lambda, M)`` | `iso_stiffness_lame_pwave` | |
| ``(\\mu, \\nu)`` | `iso_stiffness_shear_poisson` | `iso_stiffness_mu_nu`, `iso_stiffness_μ_ν` |
| ``(\\mu, M)`` | `iso_stiffness_shear_pwave` | |
| ``(\\nu, M)`` | `iso_stiffness_poisson_pwave` | |

The pair ``(\\lambda, \\nu)`` does not determine ``\\mu`` at ``\\nu = 0``, where
``\\lambda = 0`` too. Each takes the dimension as an optional leading `Val(2)` or `Val(3)`. The
pairs with the P-wave modulus have no short alias: `k_M` or `μ_M` already
name the moduli of a matrix in MeanFieldHomogenization.

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

# The builders of the pairs that are not written above.
for i in eachindex(_ISO_MODULI), j in (i + 1):length(_ISO_MODULI)
    p, q = first(_ISO_MODULI[i]), first(_ISO_MODULI[j])
    (p, q) in ((:bulk, :shear), (:young, :poisson)) && continue
    f = Symbol(:iso_stiffness_, p, :_, q)
    doc = "    $f(a, b) → TensISO{4,3}\n    $f(Val(dim), a, b) → TensISO{4,dim}\n\n" *
        "Isotropic stiffness from its $(p) and $(q) moduli (Poisson's ratio for `poisson`, " *
        "Lamé's first parameter for `lame`), in 3-D or, with `Val(2)`, in plane strain. " *
        "See [`iso_stiffness_young_poisson`](@ref) for the table of the pairs and their conventions."
    if p === :bulk
        @eval begin
            @doc $doc $f(a, b) = $f(Val(3), a, b)
            $f(v::Val, k, b) = iso_stiffness(v, k, _shear_from_bulk(v, Val($(QuoteNode(q))), k, b))
        end
    else
        @eval begin
            @doc $doc $f(a, b) = $f(Val(3), a, b)
            $f(v::Val, a, b) = _iso_lame_shear(v, _lame_shear(Val(($(QuoteNode(p)), $(QuoteNode(q)))), a, b)...)
        end
    end
end
const iso_stiffness_bulk_shear = iso_stiffness

"""
    lame_lambda(ℂ::TensISO{4,dim})

Lamé's first parameter ``\\lambda`` of an isotropic **stiffness**,
``\\mathbb{C} = (d\\lambda + 2\\mu)\\,\\mathbb{J} + 2\\mu\\,\\mathbb{K}``: in 2-D, that of the 3-D
material whose plane-strain stiffness is ``\\mathbb{C}``.

Each modulus has its reader: [`bulk_modulus`](@ref) (in 2-D, the in-plane
``\\lambda + \\mu``), [`young_modulus`](@ref), `lame_lambda`, [`shear_modulus`](@ref),
[`poisson_ratio`](@ref), [`pwave_modulus`](@ref) (``M = \\lambda + 2\\mu``). Each pair
has one too, named as its builder without `iso_stiffness_` (`bulk_shear`,
`young_poisson`, `lame_shear`, `young_pwave`, …), with the same aliases (`k_mu`,
`E_nu`, `λ_μ`, …): the conversion of a pair of moduli to another is
`young_shear(iso_stiffness_lame_poisson(λ, ν))`.
"""
function lame_lambda(ℂ::TensISO{4, dim}) where {dim}
    α, β = get_data(ℂ)
    return (α - β) / dim
end

"The bulk modulus of an isotropic stiffness, in 2-D the in-plane one ``\\lambda + \\mu``; see [`bulk_shear`](@ref)."
bulk_modulus(ℂ::TensISO{4}) = first(bulk_shear(ℂ))

"The shear modulus of an isotropic stiffness; see [`bulk_shear`](@ref)."
shear_modulus(ℂ::TensISO{4}) = last(bulk_shear(ℂ))

"Young's modulus of an isotropic stiffness; see [`young_poisson`](@ref)."
young_modulus(ℂ::TensISO{4}) = first(young_poisson(ℂ))

poisson_ratio(ℂ::TensISO{4}) = last(young_poisson(ℂ))

"The P-wave (oedometric) modulus ``M = \\lambda + 2\\mu`` of an isotropic stiffness; see [`lame_lambda`](@ref)."
pwave_modulus(ℂ::TensISO{4}) = lame_lambda(ℂ) + 2 * shear_modulus(ℂ)

# The readers of the pairs that are not written above.
for i in eachindex(_ISO_MODULI), j in (i + 1):length(_ISO_MODULI)
    p, q = first(_ISO_MODULI[i]), first(_ISO_MODULI[j])
    (p, q) in ((:bulk, :shear), (:young, :poisson), (:shear, :poisson)) && continue
    f = Symbol(p, :_, q)
    rp, rq = (Symbol(m === :lame ? :lame_lambda : m === :poisson ? :poisson_ratio : Symbol(m, :_modulus)) for m in (p, q))
    doc = "    $f(ℂ::TensISO{4}) → ($(p), $(q))\n\nThe $(p) and $(q) moduli of an isotropic stiffness " *
        "(`$rp(ℂ)`, `$rq(ℂ)`); see [`lame_lambda`](@ref) for the conventions."
    @eval @doc $doc $f(ℂ::TensISO{4}) = ($rp(ℂ), $rq(ℂ))
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
            export $(Symbol(a, :_, b)), $(Symbol(:iso_stiffness_, a, :_, b))
        end
    end
    @eval export $(Symbol(p, :_, q)), $(Symbol(:iso_stiffness_, p, :_, q))
end

export iso_stiffness, poisson_ratio, lame_lambda, bulk_modulus, shear_modulus, young_modulus, pwave_modulus

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
