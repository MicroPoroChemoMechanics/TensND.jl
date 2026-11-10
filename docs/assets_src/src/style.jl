# =============================================================================
#  style.jl — the shared look of every drawn figure of the documentation: one
#  palette, one panel, one way of writing a label. Same conventions as the
#  figures of MeanFieldHomogenization.jl, so that the two packages look alike.
#
#  Drawing: Luxor. Labels: Typst, through Typstry — mathematics in New Computer
#  Modern Math with the notation of the documentation (underlined vectors,
#  bold italic second-order tensors, blackboard fourth-order tensors), words in
#  DejaVu Sans. Every label is compiled to SVG by Typst and placed as vector
#  paths, so the output SVG is self-contained.
#
#  Colors are named by ROLE (matrix, inclusion, crack, coating…), never used as
#  raw hex codes in the figure files, so that one change here recolors every
#  figure consistently.
# =============================================================================

using Luxor
using Typstry

# ── Palette ───────────────────────────────────────────────────────────────────
# The Julia colors (Julia colors), with tints for large areas. Every figure
# sits on its own white panel: an SVG served as an image follows the operating
# system's color scheme, not the site's toggle, so it must not rely on the page
# background.

const INK = "#24292f"          # text and outlines
const MUTED = "#57606a"        # secondary text, guides
const GUIDE = "#8c959f"        # dashed guides, zoom cones
const PANEL = "#ffffff"
const PANEL_EDGE = "#d0d7de"

const BLUE = "#4063D8"
const GREEN = "#389826"
const PURPLE = "#9558B2"
const RED = "#CB3C33"
const ORANGE = "#E07B39"
const TEAL = "#1F9E89"
const GOLD = "#D4A017"
const SLATE = "#4A5568"

"Mix `hex` with white: `t = 0` keeps the color, `t = 1` gives white."
function tint(hex::AbstractString, t::Real)
    c = parse(Luxor.Colors.RGB, hex)
    f(x) = x + (1 - x) * t
    return Luxor.Colors.RGB(f(c.r), f(c.g), f(c.b))
end

# Roles: one color per coordinate, the same in every figure of the family.
const COORD_Q = BLUE           # q, τ, ρ: the ellipses and ellipsoids
const COORD_P = ORANGE         # p, ν: the hyperbolas of two sheets
const COORD_MU = GREEN         # μ: the hyperboloids of one sheet
const VECTOR = "#1a5fb4"
const CONVERGES = tint(GREEN, 0.82)
const DIVERGES = tint(RED, 0.55)

# ── Sheet ─────────────────────────────────────────────────────────────────────

"Where the sheets go: `docs/src/assets/`, under the relative path each one names."
const OUT = Ref(normpath(joinpath(@__DIR__, "..", "..", "src", "assets")))

"""
    sheet(f, relpath, w, h; panel = true)

Draw `f()` into `OUT[]/relpath.svg`, a `w × h` drawing with the origin at the
top-left corner, on a white rounded panel.
"""
function sheet(f, relpath::AbstractString, w, h; panel::Bool = true)
    path = joinpath(OUT[], relpath * ".svg")
    mkpath(dirname(path))
    empty!(_GLYPHS)
    empty!(AVOID)
    empty!(_PLACED)
    Drawing(w, h, path)
    origin()
    translate(-w / 2, -h / 2)
    if panel
        setcolor(PANEL)
        box(Point(w / 2, h / 2), w - 2, h - 2, 12, :fill)
        setcolor(PANEL_EDGE)
        setline(1.2)
        box(Point(w / 2, h / 2), w - 2, h - 2, 12, :stroke)
    end
    setline(1.6)
    setcolor(INK)
    f()
    finish()
    place_labels(path)
    compact_svg(path)
    println("  ", rpad(relpath * ".svg", 44), w, "×", h, "  ", round(Int, filesize(path) / 1024), " kB")
    return path
end

_short(x::AbstractString) = (v = round(parse(Float64, x); digits = 2); v == 0 ? "0" : rstrip(rstrip(string(v), '0'), '.'))

"""
    compact_svg(path)

Rewrite Cairo's SVG output without its redundancy: every number rounded to the
hundredth of a drawing unit (Typst glyph outlines come out with six decimals),
and the attributes that only restate SVG defaults dropped. Nothing visible
changes; the files shrink by about half.
"""
function compact_svg(path::AbstractString)
    s = read(path, String)
    s = replace(s, r"-?\d+\.\d{3,}" => _short)
    for a in (
            " fill-opacity=\"1\"", " stroke-opacity=\"1\"", " stroke-linecap=\"butt\"",
            " stroke-linejoin=\"miter\"", " fill-rule=\"nonzero\"",
        )
        s = replace(s, a => "")
    end
    write(path, s)
    return nothing
end

# ── Labels (Typst) ────────────────────────────────────────────────────────────

const TEXT_FONT = "DejaVu Sans"
const REF_PT = 20.0          # labels are compiled at this size, then scaled

"A compiled label: its size as Luxor measures it, its body, its glyph definitions."
struct Compiled
    width::Float64
    height::Float64
    body::String
    defs::Vector{Pair{String, String}}
end

const _LABELS = Dict{Tuple{String, String, String}, Compiled}()
const _GLYPHS = Dict{String, String}()                       # id => definition, for one sheet
const _PLACED = Tuple{String, Float64, Float64, Float64}[]  # (body, x, y, scale), for one sheet

_hex(c::AbstractString) = c
_hex(c) = "#" * Luxor.Colors.hex(Luxor.Colors.RGB(c))

const _DEF = r"<(symbol|clipPath|linearGradient|radialGradient|pattern)\b[^>]*\bid=\"([^\"]+)\".*?</\1>"s

function _compile(src::AbstractString, color, font::AbstractString)
    key = (String(src), _hex(color), String(font))
    return get!(_LABELS, key) do
        # The Typst compiler shipped with Typstry, called directly: a page that
        # shrinks to its content, no margin, no background.
        doc = """
        #set page(width: auto, height: auto, margin: 0pt, fill: none)
        #set text(size: $(REF_PT)pt, fill: rgb("$(_hex(color))"), font: "$(font)")
        #show math.equation: set text(font: "New Computer Modern Math")
        $(src)
        """
        io, err = IOBuffer(), IOBuffer()
        try
            Typstry.run_typst(
                c -> pipeline(c; stdin = IOBuffer(doc), stdout = io, stderr = err),
                TypstCommand(["compile", "--format", "svg", "-", "-"]),
            )
        catch
            error("Typst rejected the label `$(src)`:\n", String(take!(err)))
        end
        svg = String(take!(io))
        # Luxor's measure of the label sets the layout, as it always has; the
        # label itself is Typst's own SVG, whose glyphs are symbols named by a
        # hash of their outline, so that a sheet defines each glyph once.
        img = readsvg(svg)
        defs = [m[2] => m.match for m in eachmatch(_DEF, svg)]
        inner = svg[(findfirst('>', svg) + 1):(findlast("</svg>", svg)[1] - 1)]
        body = replace(inner, r"<defs\b.*?</defs>"s => "")
        Compiled(img.width, img.height, body, defs)
    end
end

"""
    label(src, pos; size = 18, color = INK, halign = :center, valign = :middle, bg = nothing)

Typeset the Typst source `src` (words, or mathematics between `\$…\$`) at
`pos`, with a font size of `size` drawing units, on a plate of color `bg` if
given. Returns `(width, height)`. The label is written into the sheet when it is
finished, above the drawing.
"""
function label(
        src::AbstractString, pos::Point;
        size::Real = 18, color = INK, halign::Symbol = :center, valign::Symbol = :middle,
        font::AbstractString = TEXT_FONT, bg = nothing,
    )
    lab = _compile(src, color, font)
    k = size / REF_PT
    w, h = lab.width * k, lab.height * k
    x = halign === :left ? pos.x : halign === :right ? pos.x - w : pos.x - w / 2
    y = valign === :top ? pos.y : valign === :bottom ? pos.y - h : pos.y - h / 2
    if bg !== nothing   # a backing plate, for a label over a busy drawing
        @layer (setcolor(bg); box(Point(x + w / 2, y + h / 2), w + 8, h + 6, 4, :fill))
    end
    # the user point (x, y) in the sheet's own coordinates
    M = getmatrix()
    push!(_PLACED, (lab.body, M[1] * x + M[3] * y + M[5], M[2] * x + M[4] * y + M[6], k))
    foreach(((id, d),) -> get!(_GLYPHS, id, d), lab.defs)
    return (w, h)
end

"Write the labels of a finished sheet into its SVG: the glyphs once, then each label."
function place_labels(path::AbstractString)
    isempty(_PLACED) && return nothing
    s = read(path, String)
    root = findfirst("<svg", s)[1]
    occursin("xmlns:xlink", s[root:findnext('>', s, root)]) ||
        (s = replace(s, "<svg " => "<svg xmlns:xlink=\"http://www.w3.org/1999/xlink\" "; count = 1))
    io = IOBuffer()
    print(io, "<defs>")
    foreach(id -> print(io, _GLYPHS[id]), sort!(collect(keys(_GLYPHS))))
    print(io, "</defs>")
    for (body, x, y, k) in _PLACED
        print(io, "<g transform=\"matrix(", k, " 0 0 ", k, " ", x, " ", y, ")\">", body, "</g>")
    end
    i = findlast("</svg>", s)[1]
    write(path, s[1:(i - 1)] * String(take!(io)) * s[i:end])
    return nothing
end

"Mathematics: `m(\"bb(C)^I\")` is the Typst source `\$bb(C)^I\$`."
m(s::AbstractString) = "\$" * s * "\$"

# ── Primitives ────────────────────────────────────────────────────────────────

"An ellipse of semi-axes `(a, b)` centered at `c`, its `a` axis at angle `θ`."
function ell(c::Point, a, b, θ = 0.0; fill = nothing, stroke = INK, lw = 1.4, dash = nothing)
    @layer begin
        translate(c)
        rotate(θ)
        if fill !== nothing
            setcolor(fill)
            ellipse(O, 2a, 2b, :fill)
        end
        if stroke !== nothing
            setcolor(stroke)
            setline(lw)
            dash === nothing ? setdash("solid") : setdash(dash)
            ellipse(O, 2a, 2b, :stroke)
            setdash("solid")
        end
    end
    return nothing
end

"Points of the drawing that a label must keep clear of (see `smart_label`)."
const AVOID = Point[]

"Record the segment `a`–`b` in `AVOID`, a point every 3 drawing units at most."
function avoid_segment!(a::Point, b::Point)
    n = max(1, ceil(Int, distance(a, b) / 3))
    foreach(t -> push!(AVOID, a + (b - a) * (t / n)), 0:n)
    return nothing
end

"A straight arrow from `a` to `b`."
function arrow_ab(a::Point, b::Point; color = INK, lw = 1.6, head = 9)
    avoid_segment!(a, b)
    @layer begin
        setcolor(color)
        setline(lw)
        arrow(a, b; arrowheadlength = head, arrowheadangle = π / 7, linewidth = lw)
    end
    return nothing
end

"A polyline through `pts`."
function curve(pts::AbstractVector{Point}; color = INK, lw = 1.4, dash = nothing, close = false)
    foreach(i -> avoid_segment!(pts[i], pts[i + 1]), 1:(length(pts) - 1))
    @layer begin
        setcolor(color)
        setline(lw)
        dash === nothing ? setdash("solid") : setdash(dash)
        poly(pts, :stroke; close = close)
    end
    return nothing
end

"A filled dot of radius `r` at `p`."
dot(p::Point; r = 3.2, color = INK) = (push!(AVOID, p); @layer (setcolor(color); circle(p, r, :fill)))

"Words typeset verbatim: no Typst markup is interpreted."
txt(s::AbstractString) = "#\"" * replace(s, "\\" => "\\\\", "\"" => "\\\"") * "\""

"Record the outline of the box `x0..x0+w`, `y0..y0+h` in `AVOID`."
function avoid_box!(x0, y0, w, h)
    c = [Point(x0, y0), Point(x0 + w, y0), Point(x0 + w, y0 + h), Point(x0, y0 + h), Point(x0, y0)]
    foreach(i -> avoid_segment!(c[i], c[i + 1]), 1:4)
    return nothing
end


"""
    smart_label(src, anchor; win, size = 16, color = INK, radii = (4, 7, 11, 16, 22))

Typeset `src` next to `anchor`, at the offset (distance and direction) whose
box keeps the largest clearance from everything already drawn (`AVOID`); the box
stays inside `win = (center, w, h)`. Its outline is then recorded, so that the
next label avoids it.
"""
function smart_label(
        src::AbstractString, anchor::Point; win, size = 16, color = INK,
        radii = (4, 7, 11, 16, 22), angles = 0:20:340, leader::Bool = false,
    )
    lab = _compile(src, color, TEXT_FONT)
    k = size / REF_PT
    w, h = lab.width * k, lab.height * k
    wc, ww, wh = win
    best, bestc, bestu, bestr = -Inf, anchor, Point(1, 0), 0.0
    for r in radii, θ in angles
        u = Point(cosd(θ), -sind(θ))
        c = anchor + r * u + Point(u.x * w / 2, u.y * h / 2)
        (abs(c.x - wc.x) + w / 2 > ww / 2 - 4 || abs(c.y - wc.y) + h / 2 > wh / 2 - 4) && continue
        samples = [c + Point(i * w / 2, j * h / 2) for i in -1:0.5:1, j in -1:0.5:1]
        cl = minimum(minimum(distance(q, a) for a in AVOID) for q in samples)
        score = min(cl, 5.0) - 0.03 * r
        if score > best
            best, bestc, bestu, bestr = score, c, u, r
        end
    end
    if leader && bestr >= 11   # far from its object: a thin leader
        @layer (setcolor(MUTED); setline(0.8); line(anchor + 4 * bestu, anchor + (bestr + 1) * bestu, :stroke))
    end
    label(src, bestc; size = size, color = color)
    corners = [bestc + Point(i * w / 2, j * h / 2) for (i, j) in ((-1, -1), (1, -1), (1, 1), (-1, 1), (-1, -1))]
    foreach(i -> avoid_segment!(corners[i], corners[i + 1]), 1:4)
    return bestc
end
