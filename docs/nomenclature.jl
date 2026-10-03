# The nomenclature of the documentation, read from `nomenclature.toml`: the page
# `nomenclature.md` renders it, and `make.jl` writes it as the JSON the hints
# over the equations read (`src/.vitepress/nomenclature.json`, not versioned).
# Ported from MeanFieldHomogenization, itself ported from ChemistryLab.

using JSON
using Markdown
using TOML

const NOMENCLATURE_FILE = joinpath(@__DIR__, "nomenclature.toml")

const NOMENCLATURE_GROUPS = [
    "Identities and projectors", "Bases and coordinates", "Products and operators",
    "Structured tensors", "Elastic moduli", "Differential geometry",
    "Rotations and projection",
]

"""
    nomenclature_entries() -> Vector{NamedTuple}

Every entry of `nomenclature.toml`, with an `id`.
"""
function nomenclature_entries()
    raw = TOML.parsefile(NOMENCLATURE_FILE)["symbol"]
    return map(enumerate(raw)) do (k, e)
        (;
            id = "nomen-$k", group = e["group"], tex = e["tex"], label = e["label"],
            name = e["name"], unit = get(e, "unit", ""),
            pages = String.(get(e, "pages", String[])), script = get(e, "script", false),
        )
    end
end

"""
    write_nomenclature_json(path)

The entries, as the JSON the plugin that typesets the formulas and the hints of
the theme read.
"""
function write_nomenclature_json(path)
    entries = [
        Dict(
            "id" => e.id, "tex" => e.tex, "label" => e.label, "name" => e.name,
            "unit" => e.unit, "pages" => e.pages, "script" => e.script,
        ) for e in nomenclature_entries()
    ]
    mkpath(dirname(path))
    open(io -> JSON.print(io, entries, 1), path, "w")
    return path
end

# The HTML the hints show (<sub>, <sup>, <b>, <u>), as the text of a Markdown
# table: Unicode subscripts and superscripts where every character has one,
# `_x` and `^x` otherwise.
const _SUBSCRIPTS = Dict(zip("0123456789+-−=()aehijklmnoprstuvx", "₀₁₂₃₄₅₆₇₈₉₊₋₋₌₍₎ₐₑₕᵢⱼₖₗₘₙₒₚᵣₛₜᵤᵥₓ"))
const _SUPERSCRIPTS = Dict(zip("0123456789+-−=()ni", "⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁻⁼⁽⁾ⁿⁱ"))
function _plain(s::AbstractString)
    shift(t, table, mark) = all(c -> haskey(table, c), t) ? join(table[c] for c in t) : mark * t
    s = replace(s, r"<sub>(.*?)</sub>" => m -> shift(match(r"<sub>(.*?)</sub>", m)[1], _SUBSCRIPTS, "_"))
    s = replace(s, r"<sup>(.*?)</sup>" => m -> shift(match(r"<sup>(.*?)</sup>", m)[1], _SUPERSCRIPTS, "^"))
    return replace(s, r"</?[bu]>" => "")
end

# A page of the documentation as a link from `nomenclature.md`.
_nomenclature_page_link(p) = endswith(p, "/") ? "`$p`" : "[`$(first(splitext(p)))`]($p)"

"""
    nomenclature_markdown() -> Markdown.MD

The tables of the nomenclature page, one per group: the symbol, its meaning, its
unit, and the pages a meaning is restricted to.
"""
function nomenclature_markdown()
    entries = nomenclature_entries()
    io = IOBuffer()
    for g in NOMENCLATURE_GROUPS
        rows = [e for e in entries if e.group == g]
        isempty(rows) && continue
        println(io, "## ", g, "\n")
        println(io, "| Symbol | Meaning | Unit | Where |")
        println(io, "|:--|:--|:--|:--|")
        for e in rows
            where = isempty(e.pages) ? "" : join(_nomenclature_page_link.(e.pages), ", ")
            println(io, "| ", "``", e.tex, "``", " | ", _plain(e.name), " | ", _plain(e.unit), " | ", where, " |")
        end
        println(io)
    end
    return Markdown.parse(String(take!(io)))
end
