# The data a page reads, shown where it reads them.
#
# A page that loads a file of `data/literature/` shows its content folded under
# the block that loads it, as a page shows the scripts it includes
# (`source_links.jl`): what was transcribed, from where, and how it was checked,
# laid out as a card and tables rather than as the JSON itself. Included by
# `make.jl`, after `source_links.jl`, whose `source_url` it reuses.
#
# Each function returns a `Base.HTML`, emitted by an `@example` block whose only
# line ends in `# hide`, so that the page shows the panel and no code.

using JSON

"The HTML text of `x`, with the five characters that HTML reads escaped."
_esc(x) = replace(string(x), '&' => "&amp;", '<' => "&lt;", '>' => "&gt;", '"' => "&quot;", '\'' => "&#39;")

"A cell: a missing value is a dash, a number is set in figures, anything else is text."
_cell(::Nothing) = """<td class="na">—</td>"""
_cell(x::Real) = """<td class="num">$(_esc(x))</td>"""
_cell(x) = """<td>$(_esc(x))</td>"""

"An HTML table of `header` (strings) over `rows` (vectors of already rendered cells)."
function _html_table(header, rows; class = "")
    head = join(("<th>$(_esc(h))</th>" for h in header))
    body = join(("<tr>" * join(r) * "</tr>" for r in rows), "\n")
    return """<div class="data-table-wrap"><table class="data-table $class"><thead><tr>$head</tr></thead><tbody>
    $body
    </tbody></table></div>"""
end

# Every link of a panel opens a new tab, as VitePress does for the external links
# of a page. It must: the sites of the MPCM packages share one origin, and
# VitePress's router takes a link of the same origin without a `target` for a
# page of the current site, which it cannot find (a 404 until the page is
# reloaded).
"A link to the DOI `d`, or nothing when there is none."
_doi_link(::Nothing) = ""
_doi_link(d) = """ <a target="_blank" rel="noreferrer" href="https://doi.org/$(_esc(d))">doi:$(_esc(d))</a>"""

"The folded panel: `summary` in the bar, `body` once unfolded."
_panel(summary, body) = Base.HTML(
    """<details class="data-file"><summary>$summary</summary><div class="data-file-body">
    $body
    </div></details>"""
)

"The unit of a quantity as it reads in a table: a pure number has none."
_unit(u) = u === nothing || u == "1" ? "" : """ <span class="unit">$(_esc(u))</span>"""

const _KIND_BADGE = Dict(
    "published" => "published", "fitted" => "fitted by the authors", "derived" => "derived",
)
# Most values are published as they are; a badge marks the others only.
_kind(k) = k === nothing || k == "published" ? "" :
    """<span class="badge kind-$(_esc(k))">$(_esc(get(_KIND_BADGE, k, k)))</span>"""

"""
    _sparkline(rows) -> String

An inline SVG of a two-column numerical series, its points joined in order and
scaled to the box: a glance at the shape of the curve the table lists.
"""
function _sparkline(rows; w = 160, h = 40, pad = 3)
    pts = [(float(r[1]), float(r[2])) for r in rows if r[1] isa Real && r[2] isa Real]
    length(pts) < 2 && return ""
    ## An abscissa spread over decades (a time since loading) is read on a
    ## logarithmic scale, shifted by its smallest positive value so that 0 stays.
    xpos = filter(>(0), first.(pts))
    if all(≥(0), first.(pts)) && !isempty(xpos) && maximum(xpos) > 100 * minimum(xpos)
        x₀ = minimum(xpos)
        pts = [(log1p(x / x₀), y) for (x, y) in pts]
    end
    (x0, x1), (y0, y1) = extrema(first.(pts)), extrema(last.(pts))
    sx = x1 > x0 ? (w - 2pad) / (x1 - x0) : 0.0
    sy = y1 > y0 ? (h - 2pad) / (y1 - y0) : 0.0
    path = join(("$(round(pad + sx * (x - x0); digits = 1)),$(round(h - pad - sy * (y - y0); digits = 1))" for (x, y) in pts), " ")
    return """<svg class="sparkline" viewBox="0 0 $w $h" width="$w" height="$h" aria-hidden="true"><polyline points="$path"/></svg>"""
end

"A series `{columns, rows}` as a folded table with its sparkline."
function _series(name, s)
    cols, rows = s["columns"], s["rows"]
    tbl = _html_table(cols, [[_cell(x) for x in r] for r in rows]; class = "series")
    return """<details class="data-series"><summary>$(_esc(name)): $(length(rows)) points $(_sparkline(rows))</summary>
    $tbl
    </details>"""
end

_is_series(x) = x isa AbstractDict && haskey(x, "columns") && haskey(x, "rows")

"""
A table of records (dictionaries with the same keys): the scalar fields as
columns, each nested series folded under the table, named by the record's first
field.
"""
function _records(name, recs::AbstractVector)
    isempty(recs) && return ""
    keys_ = collect(keys(first(recs)))
    scalar = filter(k -> !_is_series(first(recs)[k]), keys_)
    nested = filter(k -> _is_series(first(recs)[k]), keys_)
    tbl = _html_table(scalar, [[_cell(get(r, k, nothing)) for k in scalar] for r in recs])
    series = join(
        (_series("$(r[first(scalar)]), $k", r[k]) for r in recs for k in nested if _is_series(get(r, k, nothing))),
        "\n"
    )
    return """<h4>Table <code>$(_esc(name))</code>, $(length(recs)) records</h4>
    $tbl
    $series"""
end

"""
    show_literature(key) -> Base.HTML

The file `data/literature/<key>.json`, folded: a card for its source (citation,
DOI, where in the source the values are read, how they were transcribed and
checked), the table of its quantities, its tables and its notes, and a link to
the file at the commit the site is built from.
"""
function show_literature(key::AbstractString)
    rel = "data/literature/$key.json"
    file = joinpath(REPO_ROOT, rel)
    isfile(file) || throw(ArgumentError("show_literature: no file $file"))
    d = JSON.parsefile(file)
    get(d, "schema", nothing) == "tensnd-literature/1" ||
        throw(ArgumentError("show_literature: $rel is not of schema tensnd-literature/1"))
    src, tr = d["source"], d["transcription"]
    doi = _doi_link(src["doi"])
    checked = get(tr, "checked_against_source", false) ?
        """<span class="badge ok">checked against the source</span>""" :
        """<span class="badge warn">not checked against the source</span>"""
    card = """<div class="data-source">
    <p class="citation">$(_esc(src["citation"]))$doi</p>
    <p><span class="label">Read in</span> $(_esc(src["location"])) $checked</p>
    <p><span class="label">Transcribed from</span> $(_esc(tr["from"]))</p>
    $(get(tr, "note", nothing) === nothing ? "" : "<p class=\"note\">$(_esc(tr["note"]))</p>")
    </div>"""
    qs = d["quantities"]
    quantities = isempty(qs) ? "" : "<h4>Quantities</h4>\n" * _html_table(
            ["name", "value", "description", "where"],
            [
                [
                    "<td><code>$(_esc(k))</code></td>",
                    """<td class="num">$(_esc(q["value"]))$(_unit(q["unit"]))</td>""",
                    "<td>$(_esc(q["description"])) $(_kind(get(q, "kind", nothing)))</td>",
                    """<td class="where">$(_esc(q["location"]))</td>""",
                ] for (k, q) in qs
            ]
        )
    tables = join((_records(k, v) for (k, v) in d["tables"]), "\n")
    notes = isempty(d["notes"]) ? "" :
        "<h4>Notes</h4>\n<ul>" * join(("<li>$(_esc(n))</li>" for n in d["notes"])) * "</ul>"
    link = """<p class="data-link">The file <a target="_blank" rel="noreferrer" href="$(source_url(rel))"><code>$rel</code></a>
    on GitHub, at the commit this page was built from.</p>"""
    nq = length(qs)
    nt = sum(length, values(d["tables"]); init = 0)
    count = join(filter(!isempty, [nq > 0 ? "$nq quantities" : "", nt > 0 ? "$nt records" : ""]), ", ")
    summary = """<span class="data-icon" aria-hidden="true"></span>Data read above: <code>$rel</code>
    <span class="data-count">$count</span>"""
    return _panel(summary, join([card, quantities, tables, notes, link], "\n"))
end
