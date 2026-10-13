// ── The figures drawn by FrameFigure.vue ────────────────────────────────────
// One entry per figure: its sliders, and how it is drawn from their values.
// The formulas are read off the source of TensND:
//
//   euler        `_rot3_raw` (src/tens_projection.jl): R = Rz(φ) Ry(θ) Rz(ψ);
//   polar,
//   cylindrical,
//   spherical,
//   spheroidal   `coorsys_*` (src/coorsystems.jl): the position <u>OM</u> and
//                the normalized basis, in the order of the coordinates.
//
// Colors follow one rule: the k-th coordinate (or angle) owns `pal.q[k]`, and
// so do its slider, the angle or length it measures, the surface on which it is
// constant and the unit vector of the normalized basis, which crosses that
// surface at right angles. Construction lines are drawn in `pal.guide`.

import {
  type Annotation, type Bounds, type Palette, type Trace, type Vec2, type Vec3,
  add, arcPoints, arrow2, arrow3, fill2, grid, label2, label3, layout2, layout3,
  line2, line3, point2, point3, range, scale, sector3, surface3, vec,
} from './kit'

export interface Slider {
  /** The symbol, as Plotly/HTML markup (`<i>θ</i>`). */
  symbol: string
  min: number
  max: number
  step: number
  init: number
  /**
   * An angle: `min`, `max`, `step` and `init` then count steps of π/24
   * (`ANGLE_DIVISIONS` in ./kit), the slider shows a fraction of π and the
   * formulas receive radians.
   */
  angle: boolean
}

export interface FrameSpec {
  sliders: Slider[]
  /** Label of the checkbox showing the surfaces (lines in 2-D) of constant coordinate; none if absent. */
  iso?: string
  height: number
  build(q: number[], pal: Palette, iso: boolean, uid: string): { data: Trace[]; layout: Record<string, unknown> }
}

const E1: Vec3 = [1, 0, 0]
const E2: Vec3 = [0, 1, 0]
const E3: Vec3 = [0, 0, 1]
const O: Vec3 = [0, 0, 0]
const TWO_PI = 2 * Math.PI

// The canonical basis at the origin, with its labels.
function canonical3(pal: Palette, L = 1, data: Trace[], ann: Annotation[]): void {
  ;[E1, E2, E3].forEach((e, k) => {
    data.push(...arrow3(O, scale(L, e), pal.canonical, 0.1, 5))
    ann.push(label3(scale(L + 0.13, e), vec('e', `<sub>${k + 1}</sub>`), pal.canonical))
  })
  ann.push(label3(O, '<i>O</i>', pal.ink, { xshift: -12, yshift: -10 }))
}

// The normalized basis at M: one arrow per coordinate, in the color of that coordinate.
function localBasis3(
  M: Vec3, basis: Vec3[], labels: string[], pal: Palette, L: number, data: Trace[], ann: Annotation[],
): void {
  basis.forEach((e, k) => {
    data.push(...arrow3(M, scale(L, e), pal.q[k], 0.1, 6))
    ann.push(label3(add(M, scale(L + 0.14, e)), labels[k], pal.q[k]))
  })
  data.push(point3(M, pal.ink, 4))
  ann.push(label3(M, '<i>M</i>', pal.ink, { xshift: -12, yshift: 12 }))
}

// A label at the middle of an angular sector, just outside its arc.
function sectorLabel(c: Vec3, u: Vec3, v: Vec3, r: number, a: number, text: string, color: string): Annotation {
  const m = a / 2
  return label3(add(c, add(scale(r * Math.cos(m), u), scale(r * Math.sin(m), v))), text, color)
}

const sym = (s: string) => `<i>${s}</i>`

// The coordinate systems are seen from a little higher than the Euler angles,
// so that the surfaces of constant coordinate are not seen edge on.
const EYE: Vec3 = [1.3, 0.62, 0.8]

// ── Euler angles ────────────────────────────────────────────────────────────
// R = Rz(φ) Ry(θ) Rz(ψ). Read from the left, as rotations about the axes they
// produce: φ about e3 turns e2 into the line of nodes e^φ; θ about it tilts e3
// into n = R·e3; ψ about n turns the pair (e^θ, e^φ) in the plane tangent to
// the unit sphere at n into (e'1, e'2). (e^θ, e^φ, n) is the frame
// R(θ, φ, 0)·(e1, e2, e3), the normalized basis of the spherical coordinates.

const euler: FrameSpec = {
  sliders: [
    { symbol: sym('θ'), min: 0, max: 24, step: 1, init: 8, angle: true },
    { symbol: sym('φ'), min: -24, max: 24, step: 1, init: 10, angle: true },
    { symbol: sym('ψ'), min: -24, max: 24, step: 1, init: 6, angle: true },
  ],
  height: 520,
  build([t, f, s], pal, _iso, uid) {
    const ct = Math.cos(t), st = Math.sin(t), cf = Math.cos(f), sf = Math.sin(f)
    const cs = Math.cos(s), ss = Math.sin(s)
    const n: Vec3 = [st * cf, st * sf, ct]
    const uf: Vec3 = [cf, sf, 0]                  // direction of the meridian plane
    const eth: Vec3 = [ct * cf, ct * sf, -st]     // R(θ, φ, 0)·e1
    const eph: Vec3 = [-sf, cf, 0]                // R(θ, φ, 0)·e2, the line of nodes
    const e1p = add(scale(cs, eth), scale(ss, eph))
    const e2p = add(scale(-ss, eth), scale(cs, eph))
    const data: Trace[] = []
    const ann: Annotation[] = []

    // The unit sphere, its equator and the meridian through n.
    data.push(surface3(grid((a, b) => [Math.sin(b) * Math.cos(a), Math.sin(b) * Math.sin(a), Math.cos(b)],
      range(0, TWO_PI, 60), range(0, Math.PI, 30)), pal.body, 0.12))
    data.push(line3(arcPoints(O, E1, E2, 1, 0, TWO_PI, 96), pal.guide, 2))
    data.push(line3(arcPoints(O, E3, uf, 1, -Math.PI, Math.PI, 96), pal.guide, 2, 'dash'))
    // The axis n, and its projection on the equatorial plane.
    data.push(line3([O, n], pal.guide, 3, 'dash'))
    data.push(line3([n, [st * cf, st * sf, 0]], pal.guide, 2, 'dot'))
    data.push(line3([O, uf], pal.guide, 2, 'dash'))
    // The line of nodes, about which θ tilts e3.
    data.push(line3([scale(-1.25, eph), scale(1.25, eph)], pal.q[0], 2, 'dot'))

    canonical3(pal, 1.3, data, ann)

    // φ in the equatorial plane, θ in the meridian plane.
    data.push(...sector3(O, E1, E2, 0.4, 0, f, pal.q[1]))
    ann.push(sectorLabel(O, E1, E2, 0.53, f, sym('φ'), pal.q[1]))
    data.push(...sector3(O, E3, uf, 0.55, 0, t, pal.q[0]))
    ann.push(sectorLabel(O, E3, uf, 0.68, t, sym('θ'), pal.q[0]))

    // The tangent plane at n, the reference pair (e^θ, e^φ) and ψ between e^θ and e'1.
    const h = 0.55
    data.push(surface3(grid((a, b) => add(n, add(scale(a, eth), scale(b, eph))),
      [-h, h], [-h, h]), pal.local, 0.16))
    data.push(line3([n, add(n, scale(0.62, eth))], pal.guide, 3, 'dash'))
    data.push(line3([n, add(n, scale(0.62, eph))], pal.guide, 3, 'dash'))
    ann.push(label3(add(n, scale(0.72, eth)), vec('e', `<sup>${sym('θ')}</sup>`), pal.guide, { size: 15 }))
    ann.push(label3(add(n, scale(0.72, eph)), vec('e', `<sup>${sym('φ')}</sup>`), pal.guide, { size: 15 }))
    // Lifted off the plane by a hair, so that the two do not fight for the same pixels.
    const lift = add(n, scale(0.006, n))
    data.push(...sector3(lift, eth, eph, 0.3, 0, s, pal.q[2]))
    ann.push(sectorLabel(lift, eth, eph, 0.42, s, sym('ψ'), pal.q[2]))

    // The rotated basis, attached at n.
    const L = 0.5
    ;[e1p, e2p, n].forEach((e, k) => {
      data.push(...arrow3(n, scale(L, e), pal.local, 0.09, 6))
      ann.push(label3(add(n, scale(L + 0.13, e)), vec('e', `′<sub>${k + 1}</sub>`), pal.local))
    })
    data.push(point3(n, pal.local, 4))

    const b: Bounds = [[-1.55, 1.55], [-1.55, 1.55], [-1.55, 1.6]]
    return { data, layout: layout3(pal, b, ann, uid) }
  },
}

// ── Polar coordinates (r, θ) ────────────────────────────────────────────────
// OM = (r cos θ, r sin θ); normalized basis (e^r, e^θ).

const polar: FrameSpec = {
  sliders: [
    { symbol: sym('r'), min: 0.3, max: 1.2, step: 0.01, init: 1.0, angle: false },
    { symbol: sym('θ'), min: -24, max: 24, step: 1, init: 6, angle: true },
  ],
  iso: 'Lines of constant coordinate',
  height: 460,
  build([r, t], pal, iso, uid) {
    const c = Math.cos(t), s = Math.sin(t)
    const er: Vec2 = [c, s], et: Vec2 = [-s, c]
    const M: Vec2 = [r * c, r * s]
    const data: Trace[] = []
    const ann: Annotation[] = []
    const arc = (rad: number, a0: number, a1: number, n = 72): Vec2[] =>
      range(a0, a1, n).map((a) => [rad * Math.cos(a), rad * Math.sin(a)] as Vec2)

    if (iso) {
      data.push(line2(arc(r, 0, TWO_PI, 144), pal.q[0], 2.5))
      data.push(line2([[0, 0], [1.95 * c, 1.95 * s]], pal.q[1], 2.5))
      // Opposite M for the circle, past the end of the ray for the ray: where
      // neither the other line nor the basis at M can be.
      // Anchored on the side of the label facing away from O, so that it never
      // runs back over its line.
      const outward = (x: number, y: number) => ({
        size: 14,
        xanchor: x > 0.3 ? 'left' : x < -0.3 ? 'right' : 'center',
        yanchor: y > 0.3 ? 'bottom' : y < -0.3 ? 'top' : 'middle',
      })
      ann.push(label2([-(r + 0.06) * c, -(r + 0.06) * s], `${sym('r')} = const`, pal.q[0], outward(-c, -s)))
      ann.push(label2([2.02 * c, 2.02 * s], `${sym('θ')} = const`, pal.q[1], outward(c, s)))
    } else {
      data.push(line2([[0, 0], M], pal.guide, 1.5, 'dot'))
    }

    // θ from e1, r along OM.
    const sec = arc(0.32, 0, t)
    data.push(fill2([[0, 0], ...sec, [0, 0]], pal.q[1], 0.3))
    data.push(line2(sec, pal.q[1], 3))
    const m = t / 2
    ann.push(label2([0.44 * Math.cos(m), 0.44 * Math.sin(m)], sym('θ'), pal.q[1]))
    // r as a dimension line beside OM, on the side away from e^θ.
    const off = (k: number, p: Vec2): Vec2 => [p[0] - k * et[0], p[1] - k * et[1]]
    data.push(line2([off(0.14, [0, 0]), off(0.14, M)], pal.q[0], 2.5))
    data.push(line2([off(0.09, [0, 0]), off(0.19, [0, 0])], pal.q[0], 2.5))
    data.push(line2([off(0.09, M), off(0.19, M)], pal.q[0], 2.5))
    ann.push(label2(off(0.28, [0.5 * M[0], 0.5 * M[1]]), sym('r'), pal.q[0]))

    // Canonical basis at O, normalized basis at M.
    ann.push(arrow2([0, 0], [0.9, 0], pal.canonical), arrow2([0, 0], [0, 0.9], pal.canonical))
    ann.push(label2([1.02, -0.07], vec('e', '<sub>1</sub>'), pal.canonical))
    ann.push(label2([-0.1, 0.98], vec('e', '<sub>2</sub>'), pal.canonical))
    const L = 0.5
    ann.push(arrow2(M, [L * er[0], L * er[1]], pal.q[0], 3), arrow2(M, [L * et[0], L * et[1]], pal.q[1], 3))
    ann.push(label2([M[0] + (L + 0.12) * er[0], M[1] + (L + 0.12) * er[1]], vec('e', `<sup>${sym('r')}</sup>`), pal.q[0]))
    ann.push(label2([M[0] + (L + 0.12) * et[0], M[1] + (L + 0.12) * et[1]], vec('e', `<sup>${sym('θ')}</sup>`), pal.q[1]))

    data.push(point2([0, 0], pal.ink, 7), point2(M, pal.ink, 7))
    ann.push(label2([0, 0], '<i>O</i>', pal.ink, { xshift: -12, yshift: -12 }))
    ann.push(label2(M, '<i>M</i>', pal.ink, { xshift: -14, yshift: 12 }))
    return { data, layout: layout2(pal, [-2.3, 2.3], [-2.25, 2.25], ann, uid) }
  },
}

// ── Cylindrical coordinates (r, θ, z) ───────────────────────────────────────
// OM = (r cos θ, r sin θ, z); normalized basis (e^r, e^θ, e^z).

const cylindrical: FrameSpec = {
  sliders: [
    { symbol: sym('r'), min: 0.3, max: 1.5, step: 0.01, init: 1.0, angle: false },
    { symbol: sym('θ'), min: -24, max: 24, step: 1, init: -4, angle: true },
    { symbol: sym('z'), min: -0.9, max: 1.2, step: 0.01, init: 0.7, angle: false },
  ],
  iso: 'Surfaces of constant coordinate',
  height: 560,
  build([r, t, z], pal, iso, uid) {
    const c = Math.cos(t), s = Math.sin(t)
    const er: Vec3 = [c, s, 0], et: Vec3 = [-s, c, 0]
    const P: Vec3 = [r * c, r * s, 0]
    const M: Vec3 = [r * c, r * s, z]
    const [z0, z1, R] = [-1.2, 1.5, 1.8]
    const data: Trace[] = []
    const ann: Annotation[] = []

    if (iso) {
      data.push(surface3(grid((a, h) => [r * Math.cos(a), r * Math.sin(a), h], range(0, TWO_PI, 72), [z0, z1]), pal.q[0], 0.16))
      data.push(surface3(grid((rho, h) => [rho * c, rho * s, h], [0, R], [z0, z1]), pal.q[1], 0.18))
      data.push(surface3(grid((a, rho) => [rho * Math.cos(a), rho * Math.sin(a), z], range(0, TWO_PI, 72), [0, R]), pal.q[2], 0.15))
    }

    // The lines through M along which one coordinate varies.
    data.push(line3([[0, 0, z], add([0, 0, z], scale(R, er))], pal.guide, 2, 'dot'))
    data.push(line3(arcPoints([0, 0, z], E1, E2, r, 0, TWO_PI, 96), pal.guide, 2, 'dot'))
    data.push(line3([[P[0], P[1], z0], [P[0], P[1], z1]], pal.guide, 2, 'dot'))

    canonical3(pal, 1, data, ann)

    // θ in the plane z = 0, r from the axis, z above that plane.
    data.push(...sector3(O, E1, E2, 0.4, 0, t, pal.q[1]))
    ann.push(sectorLabel(O, E1, E2, 0.53, t, sym('θ'), pal.q[1]))
    data.push(line3([O, P], pal.q[0], 5))
    ann.push(label3(add(scale(0.5, P), scale(0.12, et)), sym('r'), pal.q[0]))
    data.push(line3([P, M], pal.q[2], 5))
    ann.push(label3(add(P, [0.1 * c, 0.1 * s, z / 2]), sym('z'), pal.q[2]))

    localBasis3(M, [er, et, E3], [vec('e', `<sup>${sym('r')}</sup>`), vec('e', `<sup>${sym('θ')}</sup>`),
      vec('e', `<sup>${sym('z')}</sup>`)], pal, 0.55, data, ann)
    const b: Bounds = [[-R, R], [-R, R], [z0 - 0.1, z1 + 0.25]]
    return { data, layout: layout3(pal, b, ann, uid, EYE) }
  },
}

// ── Spherical coordinates (θ, φ, r) ─────────────────────────────────────────
// OM = (r sin θ cos φ, r sin θ sin φ, r cos θ); normalized basis (e^θ, e^φ, e^r),
// in that order, which is R(θ, φ, 0)·(e1, e2, e3).

const spherical: FrameSpec = {
  sliders: [
    { symbol: sym('θ'), min: 0, max: 24, step: 1, init: 8, angle: true },
    { symbol: sym('φ'), min: -24, max: 24, step: 1, init: -4, angle: true },
    { symbol: sym('r'), min: 0.3, max: 1.5, step: 0.01, init: 1.1, angle: false },
  ],
  iso: 'Surfaces of constant coordinate',
  height: 560,
  build([t, f, r], pal, iso, uid) {
    const ct = Math.cos(t), st = Math.sin(t), cf = Math.cos(f), sf = Math.sin(f)
    const eth: Vec3 = [ct * cf, ct * sf, -st]
    const eph: Vec3 = [-sf, cf, 0]
    const er: Vec3 = [st * cf, st * sf, ct]
    const uf: Vec3 = [cf, sf, 0]
    const M = scale(r, er)
    const P: Vec3 = [M[0], M[1], 0]
    const R = 1.75
    const data: Trace[] = []
    const ann: Annotation[] = []

    if (iso) {
      data.push(surface3(grid((a, b) => [Math.sin(b) * Math.cos(a), Math.sin(b) * Math.sin(a), Math.cos(b)].map((x) => r * x) as Vec3,
        range(0, TWO_PI, 60), range(0, Math.PI, 30)), pal.q[2], 0.14))
      data.push(surface3(grid((a, l) => [l * st * Math.cos(a), l * st * Math.sin(a), l * ct], range(0, TWO_PI, 72), [0, 1.6]), pal.q[0], 0.16))
      data.push(surface3(grid((rho, h) => [rho * cf, rho * sf, h], [0, 1.6], [-1.45, 1.45]), pal.q[1], 0.15))
    }

    // The lines through M along which one coordinate varies.
    data.push(line3(arcPoints(O, E3, uf, r, 0, Math.PI, 72), pal.guide, 2, 'dot'))
    data.push(line3(arcPoints([0, 0, r * ct], E1, E2, r * st, 0, TWO_PI, 96), pal.guide, 2, 'dot'))
    data.push(line3([M, scale(R, er)], pal.guide, 2, 'dot'))
    data.push(line3([M, P], pal.guide, 2, 'dot'))
    data.push(line3([O, P], pal.guide, 2, 'dot'))

    canonical3(pal, 1, data, ann)

    // θ from e3 in the meridian plane, φ from e1 in the equatorial one, r along OM.
    data.push(...sector3(O, E3, uf, 0.45, 0, t, pal.q[0]))
    ann.push(sectorLabel(O, E3, uf, 0.58, t, sym('θ'), pal.q[0]))
    data.push(...sector3(O, E1, E2, 0.38, 0, f, pal.q[1]))
    ann.push(sectorLabel(O, E1, E2, 0.51, f, sym('φ'), pal.q[1]))
    data.push(line3([O, M], pal.q[2], 5))
    ann.push(label3(add(scale(0.5, M), scale(-0.14, eth)), sym('r'), pal.q[2]))

    localBasis3(M, [eth, eph, er], [vec('e', `<sup>${sym('θ')}</sup>`), vec('e', `<sup>${sym('φ')}</sup>`),
      vec('e', `<sup>${sym('r')}</sup>`)], pal, 0.55, data, ann)
    const b: Bounds = [[-R, R], [-R, R], [-R, R]]
    return { data, layout: layout3(pal, b, ann, uid, EYE) }
  },
}

// ── Prolate spheroidal coordinates (φ, p, q) ────────────────────────────────
// OM = c (√(1−p²) √(q²−1) cos φ, √(1−p²) √(q²−1) sin φ, p q), −1 ≤ p ≤ 1,
// q ≥ 1, drawn for c = 1: q = const is a prolate spheroid of foci (0, 0, ±c),
// p = const one sheet of a two-sheeted hyperboloid of the same foci, φ = const
// a half-plane bounded by the axis. Normalized basis (e^φ, e^p, e^q).

const spheroidal: FrameSpec = {
  sliders: [
    { symbol: sym('φ'), min: -24, max: 24, step: 1, init: -4, angle: true },
    { symbol: sym('p'), min: -0.95, max: 0.95, step: 0.01, init: 0.7, angle: false },
    { symbol: sym('q'), min: 1.05, max: 2, step: 0.01, init: 1.5, angle: false },
  ],
  iso: 'Surfaces of constant coordinate',
  height: 580,
  build([f, p, q], pal, iso, uid) {
    const cf = Math.cos(f), sf = Math.sin(f)
    const sp = Math.sqrt(1 - p * p), sq = Math.sqrt(q * q - 1), D = Math.sqrt(q * q - p * p)
    const at = (pp: number, qq: number, a: number): Vec3 => {
      const rho = Math.sqrt(1 - pp * pp) * Math.sqrt(qq * qq - 1)
      return [rho * Math.cos(a), rho * Math.sin(a), pp * qq]
    }
    const M = at(p, q, f)
    const eph: Vec3 = [-sf, cf, 0]
    const ep: Vec3 = [-p * sq * cf / D, -p * sq * sf / D, q * sp / D]
    const eq: Vec3 = [q * sp * cf / D, q * sp * sf / D, p * sq / D]
    const qmax = 2.3
    const data: Trace[] = []
    const ann: Annotation[] = []

    if (iso) {
      // Sampled through p = cos β and q = cosh η, which spreads the points evenly.
      data.push(surface3(grid((a, b) => at(Math.cos(b), q, a), range(0, TWO_PI, 60), range(0, Math.PI, 36)), pal.q[2], 0.14))
      data.push(surface3(grid((a, e) => at(p, Math.cosh(e), a), range(0, TWO_PI, 60), range(0, Math.acosh(qmax), 24)), pal.q[1], 0.18))
      data.push(surface3(grid((rho, h) => [rho * cf, rho * sf, h], [0, 1.9], [-2, 2]), pal.q[0], 0.14))
    }

    // The lines through M along which one coordinate varies.
    data.push(line3(arcPoints([0, 0, M[2]], E1, E2, Math.hypot(M[0], M[1]), 0, TWO_PI, 96), pal.guide, 2, 'dot'))
    data.push(line3(range(0, Math.PI, 72).map((b) => at(Math.cos(b), q, f)), pal.guide, 2, 'dot'))
    data.push(line3(range(0, Math.acosh(qmax), 48).map((e) => at(p, Math.cosh(e), f)), pal.guide, 2, 'dot'))
    data.push(line3([[0, 0, M[2]], M], pal.guide, 2, 'dot'))

    // The foci, and φ from e1 in the equatorial plane.
    data.push(line3([[0, 0, -1], [0, 0, 1]], pal.ink, 4))
    data.push(point3([0, 0, 1], pal.ink, 4), point3([0, 0, -1], pal.ink, 4))
    ann.push(label3([0, 0, 1], '<i>c</i>', pal.ink, { xshift: 14 }))
    ann.push(label3([0, 0, -1], '−<i>c</i>', pal.ink, { xshift: 16 }))

    canonical3(pal, 1, data, ann)
    const uf: Vec3 = [cf, sf, 0]
    data.push(line3([O, scale(1.1, uf)], pal.guide, 2, 'dot'))
    data.push(...sector3(O, E1, E2, 0.42, 0, f, pal.q[0]))
    ann.push(sectorLabel(O, E1, E2, 0.55, f, sym('φ'), pal.q[0]))

    localBasis3(M, [eph, ep, eq], [vec('e', `<sup>${sym('φ')}</sup>`), vec('e', `<sup>${sym('p')}</sup>`),
      vec('e', `<sup>${sym('q')}</sup>`)], pal, 0.55, data, ann)
    const b: Bounds = [[-2.1, 2.1], [-2.1, 2.1], [-qmax, qmax]]
    return { data, layout: layout3(pal, b, ann, uid, EYE) }
  },
}

export const FRAMES: Record<string, FrameSpec> = { euler, polar, cylindrical, spherical, spheroidal }

