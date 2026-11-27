# =============================================================================
#  symmetries.jl — the symmetry classes of the tensors of the plane.
#
#  One sheet, `symmetries/plane_classes`: a shape that has the symmetries of
#  each class (a disc, a square, a rectangle, a shape without symmetry), its
#  mirror lines, the frame the class carries, and its number of constants.
#  The shapes of the square and of the rectangle are turned by the same angle,
#  that of their frame.
# =============================================================================

const MIRROR = ORANGE          # a mirror line
const SHAPE = tint(BLUE, 0.80)

# A line through `c` at the angle `φ` (counted from the horizontal, upward).
mirror_line(c::Point, φ, half) =
    curve([c + half * Point(cos(φ), -sin(φ)), c - half * Point(cos(φ), -sin(φ))]; color = MIRROR, lw = 1.5, dash = [6, 4])

# The polygon of the corners `(x, y)`, turned by `θ` about `c`.
function turned_polygon(c::Point, corners, θ; fill = SHAPE)
    pts = [c + Point(x * cos(θ) - y * sin(θ), -(x * sin(θ) + y * cos(θ))) for (x, y) in corners]
    @layer (setcolor(fill); poly(pts, :fill; close = true))
    curve(pts; color = INK, lw = 1.5, close = true)
    return pts
end

# The frame (e₁, e₂) at the angle `θ`.
function plane_frame(c::Point, θ; len = 46)
    e₁, e₂ = Point(cos(θ), -sin(θ)), Point(-sin(θ), -cos(θ))
    arrow_ab(c, c + len * e₁; color = VECTOR, lw = 2.0)
    arrow_ab(c, c + len * e₂; color = VECTOR, lw = 2.0)
    # beside the arrows, clear of the mirror lines they lie on
    label(m("underline(e)_1"), c + (len - 4) * e₁ - 15 * e₂; size = 16, color = VECTOR)
    label(m("underline(e)_2"), c + (len - 4) * e₂ - 17 * e₁; size = 16, color = VECTOR)
    return nothing
end

function plane_classes()
    return sheet("symmetries/plane_classes", 860, 300) do
        θ = π / 9
        xs = (115.0, 325.0, 535.0, 745.0)
        y = 118.0
        # isotropic: every line through the center is a mirror line
        c = Point(xs[1], y)
        @layer (setcolor(SHAPE); circle(c, 58, :fill))
        @layer (setcolor(INK); setline(1.5); circle(c, 58, :stroke))
        foreach(φ -> mirror_line(c, φ, 76), range(0, π; length = 7)[1:6])
        # square: the two axes of the frame and the two diagonals
        c = Point(xs[2], y)
        turned_polygon(c, ((52, 52), (-52, 52), (-52, -52), (52, -52)), θ)
        foreach(φ -> mirror_line(c, θ + φ, 84), (0, π / 4, π / 2, 3π / 4))
        plane_frame(c, θ)
        # orthotropic: the two axes of the frame
        c = Point(xs[3], y)
        turned_polygon(c, ((66, 36), (-66, 36), (-66, -36), (66, -36)), θ)
        foreach(φ -> mirror_line(c, θ + φ, 84), (0, π / 2))
        plane_frame(c, θ)
        # no symmetry
        c = Point(xs[4], y)
        turned_polygon(c, ((62, 14), (20, 50), (-58, 30), (-40, -44), (30, -36)), 0.0)
        for (x, name, n) in zip(xs, ("isotropic", "square", "orthotropic", "general"), (2, 3, 4, 6))
            label(txt(name), Point(x, 222); size = 17)
            label(txt("$n constants"), Point(x, 246); size = 15, color = MUTED)
        end
        legend_item(Point(xs[2] - 60, 278), MIRROR, txt("mirror line"))
    end
end

symmetries_figures() = (plane_classes(); nothing)
