# =============================================================================
#  coordinates.jl — the confocal coordinate systems of the theory pages:
#  the spheroidal meridian plane, the three sections of the ellipsoidal chart,
#  the domain of convergence of the re-expansion of exterior harmonics.
#  Docs: theory/spheroidal.md, theory/ellipsoidal.md,
#        theory/spheroidal_harmonics.md.
# =============================================================================

"The frame of a panel: center `c`, scale `s` pixels per unit, y axis up."
struct Frame
    c::Point
    s::Float64
end
(f::Frame)(u, v) = f.c + Point(f.s * u, -f.s * v)

"Run `body()` with drawing clipped to the window of size `w × h` centered at `c`."
function windowed(body, c::Point, w, h)
    @layer begin
        box(c, w, h, :clip)
        body()
    end
    @layer (setcolor(PANEL_EDGE); setline(1.0); box(c, w, h, :stroke))
    return nothing
end

"Both axes of a frame, with arrowheads, and their names."
function axes_of(f::Frame, w, h, horiz::AbstractString, vert::AbstractString)
    arrow_ab(f(-w / 2 / f.s + 0.04, 0), f(w / 2 / f.s - 0.02, 0); color = MUTED, lw = 1.1, head = 7)
    arrow_ab(f(0, -h / 2 / f.s + 0.04), f(0, h / 2 / f.s - 0.02); color = MUTED, lw = 1.1, head = 7)
    p1 = f(w / 2 / f.s - 0.02, 0) + Point(-4, 12)
    wh = label(m(horiz), p1; size = 16, halign = :right, valign = :top, color = MUTED)
    avoid_box!(p1.x - wh[1], p1.y, wh[1], wh[2])
    p2 = f(0, h / 2 / f.s - 0.02) + Point(12, 4)
    wv = label(m(vert), p2; size = 16, halign = :left, valign = :top, color = MUTED)
    avoid_box!(p2.x, p2.y, wv[1], wv[2])
    return nothing
end

"The pair of unit vectors at `M` (screen), right-angle mark and names."
function frame_vectors(M::Point, a, b, na, nb, ca, cb; win, len = 38)
    # a, b: unit vectors in drawing axes (u to the right, v up)
    sa, sb = Point(a[1], -a[2]), Point(b[1], -b[2])
    k = 9.0
    @layer begin
        setcolor(MUTED)
        setline(1.0)
        poly([M + k * sa, M + k * (sa + sb), M + k * sb], :stroke)
    end
    arrow_ab(M, M + len * sa; color = ca, lw = 2.2, head = 9)
    arrow_ab(M, M + len * sb; color = cb, lw = 2.2, head = 9)
    dot(M; r = 3.4)
    smart_label(m(na), M + len * sa; win = win, size = 17, color = ca, radii = (3, 6, 9))
    smart_label(m(nb), M + len * sb; win = win, size = 17, color = cb, radii = (3, 6, 9))
    return nothing
end

"A legend line sample followed by its text."
function legend_item(p::Point, color, src; lw = 1.8)
    @layer (setcolor(color); setline(lw); line(p, p + Point(26, 0), :stroke))
    label(src, p + Point(34, 0); size = 16, halign = :left)
    return nothing
end

# ── 1. The meridian plane of the two spheroidal systems ──────────────────────

"One panel: confocal conics, foci, the point M and its two unit vectors."
function meridian_panel(center::Point, W, H; prolate::Bool)
    f = Frame(center, prolate ? 64.0 : 62.0)
    c = 1.0
    ρ(p, q) = c * sqrt(1 - p^2) * (prolate ? sqrt(q^2 - 1) : sqrt(q^2 + 1))
    z(p, q) = c * p * q
    # the drawn values
    qs = prolate ? (1.6, 2.0, 2.3) : (0.5, 0.9, 1.5, 2.2)
    ps = prolate ? (0.25, 0.55, 0.75) : (0.5, 0.75, 0.92)
    q0, p0 = prolate ? (2.0, 0.55) : (0.9, 0.5)
    @assert q0 in qs && p0 in ps
    # M and the unit vectors, from the formulas of the page
    Mρ, Mz = ρ(p0, q0), z(p0, q0)
    d = prolate ? sqrt(q0^2 - p0^2) : sqrt(q0^2 + p0^2)
    sq = prolate ? sqrt(q0^2 - 1) : sqrt(q0^2 + 1)
    ep = (-p0 * sq / d, q0 * sqrt(1 - p0^2) / d)
    eq = prolate ? (q0 * sqrt(1 - p0^2) / d, p0 * sq / d) : (q0 * sqrt(1 - p0^2) / d, p0 * sq / d)
    @assert isapprox(ep[1]^2 + ep[2]^2, 1; atol = 1.0e-12) && isapprox(eq[1]^2 + eq[2]^2, 1; atol = 1.0e-12)
    @assert isapprox(ep[1] * eq[1] + ep[2] * eq[2], 0; atol = 1.0e-12)
    # M lies on the drawn ellipse q0 and on the drawn hyperbola p0
    a, b = prolate ? (sqrt(q0^2 - 1), q0) : (sqrt(q0^2 + 1), q0)
    @assert isapprox((Mρ / a)^2 + (Mz / b)^2, 1; atol = 1.0e-12)
    σ = prolate ? 1 : -1   # prolate: (z/p)² − ϱ²/(1−p²) = 1; oblate: the opposite sign
    @assert isapprox(σ * ((Mz / p0)^2 - (Mρ / sqrt(1 - p0^2))^2), 1; atol = 1.0e-12)
    # e^q is normal to the ellipse q = q0 and e^p to the hyperbola p = p0: parallel to the gradients of their equations
    gq = (2Mρ / a^2, 2Mz / b^2)
    @assert abs(gq[1] * eq[2] - gq[2] * eq[1]) < 1.0e-12
    gp = σ .* (-2Mρ / (1 - p0^2), 2Mz / p0^2)
    @assert abs(gp[1] * ep[2] - gp[2] * ep[1]) < 1.0e-12

    windowed(center, W, H) do
        # ellipses
        for q in qs
            θ = range(0, 2π; length = 241)
            pts = prolate ? [f(sqrt(q^2 - 1) * cos(t), q * sin(t)) for t in θ] :
                [f(sqrt(q^2 + 1) * cos(t), q * sin(t)) for t in θ]
            curve(pts; color = COORD_Q, lw = 1.5)
        end
        # hyperbolas, both signs of p, both sides of the axis
        for p in ps, sg in (1, -1)
            pp = sg * p
            if prolate
                t = range(-2.4, 2.4; length = 241)
                pts = [f(sqrt(1 - p^2) * sinh(u), pp * cosh(u)) for u in t]
                curve(pts; color = COORD_P, lw = 1.5)
            else
                u = range(0, 2.4; length = 161)
                for sd in (1, -1)
                    curve([f(sd * sqrt(1 - p^2) * cosh(v), pp * sinh(v)) for v in u]; color = COORD_P, lw = 1.5)
                end
            end
        end
        axes_of(f, W, H, "rho.alt", "x_3")
        # foci and c
        if prolate
            curve([f(0, 0), f(0, c)]; color = INK, lw = 2.8)
            dot(f(0, c)); dot(f(0, -c))
        else
            curve([f(0, 0), f(c, 0)]; color = INK, lw = 2.8)
            dot(f(c, 0)); dot(f(-c, 0))
        end
    end
    return (f, (Mρ, Mz), ep, eq, c)
end

function spheroidal_meridian()
    return sheet("coordinates/spheroidal_meridian", 740, 440) do
        W, H = 340.0, 350.0
        for (k, prolate) in enumerate((true, false))
            center = Point(10 + W / 2 + (k - 1) * (W + 20), 48 + H / 2)
            label(txt(prolate ? "prolate" : "oblate"), Point(center.x, 24); size = 18)
            f, (Mρ, Mz), ep, eq, c = meridian_panel(center, W, H; prolate = prolate)
            M = f(Mρ, Mz)
            win = (center, W, H)
            if prolate
                Fp, Fm = f(0, c), f(0, -c)
                @layer (setcolor(MUTED); setline(1.3); setdash([5, 3]); line(M, Fp, :stroke); line(M, Fm, :stroke))
                avoid_segment!(M, Fp); avoid_segment!(M, Fm)
                frame_vectors(M, ep, eq, "underline(e)^p", "underline(e)^q", COORD_P, COORD_Q; win = win)
                smart_label(m("r_(+)"), Fp + 0.28 * (M - Fp); win = win, size = 17, color = MUTED, radii = (2, 5, 9, 14, 20))
                smart_label(m("r_(-)"), M + 0.5 * (Fm - M); win = win, size = 17, color = MUTED, radii = (2, 5, 9, 14, 20))
                smart_label(m("F_(+)"), Fp; win = win, size = 16, radii = (4, 7))
                smart_label(m("F_(-)"), Fm; win = win, size = 16, radii = (4, 7))
                smart_label(m("c"), f(0, c / 2); win = win, size = 16, radii = (3, 6, 9))
            else
                frame_vectors(M, ep, eq, "underline(e)^p", "underline(e)^tau", COORD_P, COORD_Q; win = win)
                smart_label(m("c"), f(c / 2, 0); win = win, size = 16, radii = (3, 6, 9))
                label(txt("focal circle"), f(0, 0) + Point(0, 14); size = 12, color = MUTED, valign = :top)
                for sg in (1, -1)
                    @layer (setcolor(MUTED); setline(0.9); line(f(sg * 0.2c, 0) + Point(0, 12), f(sg * c, 0) + Point(-sg * 3, 3), :stroke))
                end
            end
            smart_label(m("M"), M; win = win, size = 17, radii = (4, 7, 11, 16, 22, 30), angles = 0:15:345)
            # legend below the window
            y = 48 + H + 20
            x0 = center.x - 112
            legend_item(Point(x0, y), COORD_Q, m(prolate ? "q = \"const\"" : "tau = \"const\""))
            legend_item(Point(x0 + 130, y), COORD_P, m("p = \"const\""))
        end
    end
end

# ── 2. The three sections of the ellipsoidal chart ───────────────────────────

const H2, H3 = 1.2, 0.7
const H1 = sqrt(H2^2 - H3^2)

"x = a cosh t, y = b sinh t, both branches, for |t| ≤ T."
hyper_pts(f::Frame, a, b, T; swap = false) =
    [swap ? f(b * sinh(t), a * cosh(t)) : f(a * cosh(t), b * sinh(t)) for t in range(-T, T; length = 241)]

function ellipsoidal_sections()
    return sheet("coordinates/ellipsoidal_sections", 1060, 446) do
        W, H = 330.0, 330.0
        s = 72.0
        ρs = (1.5, 1.8, 2.1)
        μs = (0.8, 0.95)
        νs = (0.3, 0.5)
        T = 3.0
        titles = ("x_3 = 0", "x_2 = 0", "x_1 = 0")
        axnames = (("x_1", "x_2"), ("x_1", "x_3"), ("x_2", "x_3"))
        for k in 1:3
            center = Point(15 + W / 2 + (k - 1) * (W + 20), 46 + H / 2)
            f = Frame(center, s)
            label(m(titles[k]), Point(center.x, 24); size = 18)
            windowed(center, W, H) do
                if k == 1      # x3 = 0: x1²/κ + x2²/(κ−h3²) = 1
                    for r in ρs
                        curve([f(r * cos(t), sqrt(r^2 - H3^2) * sin(t)) for t in range(0, 2π; length = 241)]; color = COORD_Q)
                    end
                    for u in μs
                        curve([f(u * cos(t), sqrt(u^2 - H3^2) * sin(t)) for t in range(0, 2π; length = 241)]; color = COORD_MU)
                    end
                    for n in νs, sg in (1, -1)
                        curve([f(sg * n * cosh(t), sqrt(H3^2 - n^2) * sinh(t)) for t in range(-T, T; length = 241)]; color = COORD_P)
                    end
                    # focal ellipse x1²/h2² + x2²/h1² = 1
                    curve([f(H2 * cos(t), H1 * sin(t)) for t in range(0, 2π; length = 241)]; color = INK, lw = 2.8)
                elseif k == 2  # x2 = 0: x1²/κ − x3²/(h2²−κ) = 1, κ ∈ (0, h2²); ρ ellipses
                    for r in ρs
                        curve([f(r * cos(t), sqrt(r^2 - H2^2) * sin(t)) for t in range(0, 2π; length = 241)]; color = COORD_Q)
                    end
                    for (vals, col) in ((μs, COORD_MU), (νs, COORD_P)), v in vals, sg in (1, -1)
                        curve([f(sg * v * cosh(t), sqrt(H2^2 - v^2) * sinh(t)) for t in range(-T, T; length = 241)]; color = col)
                    end
                    # focal hyperbola x1²/h3² − x3²/h1² = 1
                    for sg in (1, -1)
                        curve([f(sg * H3 * cosh(t), H1 * sinh(t)) for t in range(-T, T; length = 241)]; color = INK, lw = 2.8)
                    end
                else           # x1 = 0: ρ ellipses and μ hyperbolas, no ν
                    for r in ρs
                        curve([f(sqrt(r^2 - H3^2) * cos(t), sqrt(r^2 - H2^2) * sin(t)) for t in range(0, 2π; length = 241)]; color = COORD_Q)
                    end
                    for u in μs, sg in (1, -1)
                        curve([f(sg * sqrt(u^2 - H3^2) * cosh(t), sqrt(H2^2 - u^2) * sinh(t)) for t in range(-T, T; length = 241)]; color = COORD_MU)
                    end
                end
                axes_of(f, W, H, axnames[k][1], axnames[k][2])
            end
            # lengths and foci
            marks = k == 1 ? ((H3, 0, "h_3", :b), (H2, 0, "h_2", :b), (0, H1, "h_1", :r)) :
                k == 2 ? ((H3, 0, "h_3", :b), (H2, 0, "h_2", :b)) : ((H1, 0, "h_1", :b),)
            k == 1 && (dot(f(-H3, 0)); dot(f(0, -H1)); dot(f(-H2, 0)))
            k == 2 && (dot(f(-H3, 0)); dot(f(-H2, 0)))
            k == 3 && dot(f(-H1, 0))
            # each length is read on a thin dashed guide that runs to the border of the panel,
            # where its name stands, clear of the curves
            foreach(((u, v, name, side),) -> dot(f(u, v)), marks)
            for (u, v, name, side) in marks
                a = f(u, v)
                if side === :b
                    e = Point(a.x, center.y + H / 2)
                    @layer (setcolor(GUIDE); setline(0.9); setdash([3, 3]); line(a, e, :stroke))
                    label(m(name), e + Point(0, 11); size = 16)
                else
                    e = Point(center.x + W / 2, a.y)
                    @layer (setcolor(GUIDE); setline(0.9); setdash([3, 3]); line(a, e, :stroke))
                    label(m(name), e + Point(-6, -3); size = 16, halign = :right, valign = :bottom)
                end
            end
        end
        # legend
        y = 46 + H + 40
        x0 = 1060 / 2 - 250
        legend_item(Point(x0, y), COORD_Q, m("rho = \"const\""))
        legend_item(Point(x0 + 150, y), COORD_MU, m("mu = \"const\""))
        legend_item(Point(x0 + 300, y), COORD_P, m("nu = \"const\""))
        legend_item(Point(x0 + 450, y), INK, txt("focal curve"); lw = 3)
    end
end

# ── 3. Domain of convergence of the re-expansion ─────────────────────────────

"Diagonal hatching of the ellipse of semi-axes `(a, b)` centered at `c`, spacing `d`."
function hatch_ellipse(c::Point, a, b; d = 9.0, color = DIVERGES)
    @layer begin
        translate(c)
        ellipse(O, 2a, 2b, :clip)
        setcolor(color)
        setline(1.2)
        for x in -(a + b):d:(a + b)
            line(Point(x - b, b), Point(x + b, -b), :stroke)
        end
    end
    return nothing
end

function reexpansion_domain()
    return sheet("coordinates/reexpansion_domain", 740, 440) do
        W, H = 340.0, 360.0
        for k in 1:2
            center = Point(10 + W / 2 + (k - 1) * (W + 20), 50 + H / 2)
            label(m(k == 1 ? "c_1 < c_2" : "c_1 > c_2"), Point(center.x, 26); size = 18)
            s, c1, c2 = k == 1 ? (112.0, 0.5, 1.0) : (125.0, 1.0, 0.6)
            f = Frame(center, s)
            @layer (setcolor(CONVERGES); box(center, W, H, :fill))
            if k == 2
                q2 = c1 / c2
                @assert isapprox(sqrt(c1^2 - c2^2), c2 * sqrt(q2^2 - 1); atol = 1.0e-12)  # same spheroid
                a, b = c1 * s, sqrt(c1^2 - c2^2) * s
                @layer (setcolor(PANEL); ell(center, b, a; fill = PANEL, stroke = nothing))
                hatch_ellipse(center, b, a)
                ell(center, b, a; stroke = INK, lw = 2.0, dash = [7, 4])
            end
            windowed(center, W, H) do
                axes_of(f, W, H, "rho.alt", "x_3")
            end
            # focal segments: source (blue, right of the axis), target (orange, left)
            dx = 3.5
            curve([f(0, -c1) + Point(dx, 0), f(0, c1) + Point(dx, 0)]; color = COORD_Q, lw = 3.2)
            curve([f(0, -c2) + Point(-dx, 0), f(0, c2) + Point(-dx, 0)]; color = COORD_P, lw = 4.2)
            for v in (-c1, c1)
                curve([f(0, v) + Point(dx, 0), f(0, v) + Point(dx + 7, 0)]; color = COORD_Q, lw = 2)
            end
            for v in (-c2, c2)
                curve([f(0, v) + Point(-dx, 0), f(0, v) + Point(-dx - 7, 0)]; color = COORD_P, lw = 2)
            end
            k == 1 ? label(m("c_1"), f(0, c1) + Point(dx + 11, 0); size = 17, halign = :left, color = COORD_Q) :
                label(m("c_1"), f(0, c1) + Point(dx + 4, -9); size = 17, halign = :left, valign = :bottom, color = COORD_Q)
            label(m("c_2"), f(0, c2) + Point(-dx - 11, 0); size = 17, halign = :right, color = COORD_P, bg = k == 2 ? PANEL : nothing)
            if k == 1
                label(txt("converges"), center + Point(-W / 2 + 14, -H / 2 + 20); size = 15, halign = :left, color = tint(GREEN, 0.0))
            else
                label(txt("converges"), center + Point(-W / 2 + 14, -H / 2 + 20); size = 15, halign = :left, color = GREEN)
                label(txt("diverges"), center + Point(14, 72); size = 15, halign = :left, color = RED, bg = PANEL)
                a, b = c1 * s, sqrt(c1^2 - c2^2) * s
                # equatorial semi-axis
                @layer (setcolor(INK); setline(1.2); setdash([3, 3]); line(center, center + Point(b, 0), :stroke))
                dot(center + Point(b, 0); r = 2.8)
                label(m("sqrt(c_1^2 - c_2^2)"), center + Point(b / 2 + 8, 8); size = 14, valign = :top, bg = PANEL)
                label(m("q_2 = c_1 / c_2"), center + Point(b + 12, -a * 0.62); size = 17, halign = :left, bg = CONVERGES)
            end
        end
    end
end

coordinates_figures() = (spheroidal_meridian(); ellipsoidal_sections(); reexpansion_domain(); nothing)
