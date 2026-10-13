// ── The interactive figures of the theory pages ─────────────────────────────
// What the figure components (EulerAnglesFigure, CoordinateSystemFigure) share:
// Plotly, loaded once and only on a page that draws one; the colors of the
// current scheme; and the few shapes they are built from — arrows, arcs, angular
// sectors, surfaces and labels.
//
// Everything is computed in the browser, because the sliders move the picture
// continuously: there is no grid of states precomputed in Julia to choose from.
// The formulas are those of the pages, which read them off the source of
// TensND (`_rot3_raw`, `coorsys_*`).

export type Vec3 = [number, number, number]
export type Trace = Record<string, unknown>
export type Annotation = Record<string, unknown>

// Plotly is a bundled dependency (`plotly.js-gl3d-dist-min`, as in
// MeanFieldHomogenization), imported lazily and once whatever the number of
// figures on the page. The gl3d bundle also carries the 2-D `scatter` trace.
let plotlyModule: Promise<any> | null = null

export async function loadPlotly(): Promise<any> {
  if (!plotlyModule) plotlyModule = import('plotly.js-gl3d-dist-min')
  // CommonJS bundle: depending on Vite's interop the object is the default
  // export or the namespace itself.
  const mod = (await plotlyModule) as any
  return mod.default ?? mod
}

export const PLOTLY_CONFIG = { responsive: true, displaylogo: false }

// ── Colors ──────────────────────────────────────────────────────────────────
// By role, from the `--tensnd-fig-*` tokens of `theme/overrides.css`, which are
// defined once per color scheme. Plotly cannot read a CSS variable (most of a 3D
// scene is WebGL), so the figure reads their current values each time it is
// drawn, and is drawn again when the reader toggles the scheme.
//
// `q` holds one color per coordinate, in the order of the chart: the slider of
// the coordinate, the angle or length it measures, the surface on which it is
// constant and the unit vector of the normalized basis that crosses that
// surface all share it.

export interface Palette {
  ink: string
  guide: string
  canonical: string
  local: string
  body: string
  q: [string, string, string]
}

export function palette(): Palette {
  const style = getComputedStyle(document.documentElement)
  const v = (name: string) => style.getPropertyValue(`--tensnd-fig-${name}`).trim()
  return {
    ink: v('ink'),
    guide: v('guide'),
    canonical: v('canonical'),
    local: v('local'),
    body: v('body'),
    q: [v('q1'), v('q2'), v('q3')],
  }
}

// VitePress switches the scheme by toggling the `dark` class on `<html>`.
export function onSchemeChange(callback: () => void): () => void {
  const root = document.documentElement
  let dark = root.classList.contains('dark')
  const observer = new MutationObserver(() => {
    if (root.classList.contains('dark') === dark) return
    dark = !dark
    callback()
  })
  observer.observe(root, { attributes: true, attributeFilter: ['class'] })
  return () => observer.disconnect()
}

// ── Vectors ─────────────────────────────────────────────────────────────────

// The sliders move an angle by steps of π/ANGLE_DIVISIONS and show it as an exact
// fraction of π. With 24 divisions every value falls on a usual denominator
// (2, 3, 4, 6, 8, 12 or 24), in steps small enough for the picture to move smoothly.
export const ANGLE_DIVISIONS = 24

const gcd = (a: number, b: number): number => (b === 0 ? a : gcd(b, a % b))

/** `k π / ANGLE_DIVISIONS` as a reduced fraction of π: `0`, `π/4`, `−5π/12`, `π`. */
export function piFraction(k: number): string {
  if (k === 0) return '0'
  const g = gcd(Math.abs(k), ANGLE_DIVISIONS)
  const num = Math.abs(k) / g, den = ANGLE_DIVISIONS / g
  const top = num === 1 ? 'π' : `${num}π`
  return `${k < 0 ? '−' : ''}${top}${den === 1 ? '' : `/${den}`}`
}

export const add = (a: Vec3, b: Vec3): Vec3 => [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
export const sub = (a: Vec3, b: Vec3): Vec3 => [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
export const scale = (k: number, a: Vec3): Vec3 => [k * a[0], k * a[1], k * a[2]]
export const dot = (a: Vec3, b: Vec3): number => a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
export const cross = (a: Vec3, b: Vec3): Vec3 => [
  a[1] * b[2] - a[2] * b[1],
  a[2] * b[0] - a[0] * b[2],
  a[0] * b[1] - a[1] * b[0],
]
export const norm = (a: Vec3): number => Math.hypot(a[0], a[1], a[2])
export const unit = (a: Vec3): Vec3 => scale(1 / norm(a), a)

/** `c + r (cos t u + sin t v)` for `t` from `a0` to `a1`: an arc of the plane (u, v). */
export function arcPoints(c: Vec3, u: Vec3, v: Vec3, r: number, a0: number, a1: number, n = 48): Vec3[] {
  const pts: Vec3[] = []
  for (let k = 0; k <= n; k++) {
    const t = a0 + ((a1 - a0) * k) / n
    pts.push(add(c, add(scale(r * Math.cos(t), u), scale(r * Math.sin(t), v))))
  }
  return pts
}

/** The grid `f(u, v)` as the three matrices of a `surface` trace. */
export function grid(f: (u: number, v: number) => Vec3, us: number[], vs: number[]) {
  const X: number[][] = [], Y: number[][] = [], Z: number[][] = []
  for (const v of vs) {
    const x: number[] = [], y: number[] = [], z: number[] = []
    for (const u of us) {
      const p = f(u, v)
      x.push(p[0]); y.push(p[1]); z.push(p[2])
    }
    X.push(x); Y.push(y); Z.push(z)
  }
  return { X, Y, Z }
}

export const range = (a: number, b: number, n: number): number[] =>
  Array.from({ length: n + 1 }, (_, k) => a + ((b - a) * k) / n)

// ── 3-D traces ──────────────────────────────────────────────────────────────

export function line3(pts: Vec3[], color: string, width = 4, dash = 'solid'): Trace {
  return {
    type: 'scatter3d', mode: 'lines',
    x: pts.map((p) => p[0]), y: pts.map((p) => p[1]), z: pts.map((p) => p[2]),
    line: { color, width, dash }, hoverinfo: 'skip', showlegend: false,
  }
}

export function point3(p: Vec3, color: string, size = 5): Trace {
  return {
    type: 'scatter3d', mode: 'markers', x: [p[0]], y: [p[1]], z: [p[2]],
    marker: { color, size }, hoverinfo: 'skip', showlegend: false,
  }
}

/**
 * An arrow from `o` to `o + d`. The head is a cone meshed explicitly, so that it
 * keeps its length whatever the scene: a Plotly `cone` trace rescales itself
 * from the data instead.
 */
export function arrow3(o: Vec3, d: Vec3, color: string, head = 0.11, width = 6): Trace[] {
  const L = norm(d)
  const u = scale(1 / L, d)
  const tip = add(o, d)
  const h = Math.min(head, 0.6 * L)
  const r = 0.36 * h
  const a: Vec3 = Math.abs(u[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0]
  const e1 = unit(sub(a, scale(dot(a, u), u)))
  const e2 = cross(u, e1)
  const base = sub(tip, scale(h, u))
  const n = 18
  const V: Vec3[] = [tip, base]
  for (let k = 0; k < n; k++) {
    const t = (2 * Math.PI * k) / n
    V.push(add(base, add(scale(r * Math.cos(t), e1), scale(r * Math.sin(t), e2))))
  }
  const I: number[] = [], J: number[] = [], K: number[] = []
  for (let k = 0; k < n; k++) {
    const i = k + 2, j = ((k + 1) % n) + 2
    I.push(0, 1); J.push(i, j); K.push(j, i)
  }
  return [
    line3([o, sub(tip, scale(0.9 * h, u))], color, width),
    {
      type: 'mesh3d', x: V.map((p) => p[0]), y: V.map((p) => p[1]), z: V.map((p) => p[2]),
      i: I, j: J, k: K, color, flatshading: true,
      lighting: { ambient: 0.85, diffuse: 0.3, specular: 0.05 },
      hoverinfo: 'skip', showlegend: false,
    },
  ]
}

/** A filled angular sector, center `c`, in the plane (u, v), with its arc. */
export function sector3(
  c: Vec3, u: Vec3, v: Vec3, r: number, a0: number, a1: number, color: string, opacity = 0.32,
): Trace[] {
  const pts = arcPoints(c, u, v, r, a0, a1)
  const V = [c, ...pts]
  const I: number[] = [], J: number[] = [], K: number[] = []
  for (let k = 1; k < V.length - 1; k++) { I.push(0); J.push(k); K.push(k + 1) }
  return [
    {
      type: 'mesh3d', x: V.map((p) => p[0]), y: V.map((p) => p[1]), z: V.map((p) => p[2]),
      i: I, j: J, k: K, color, opacity, flatshading: true,
      lighting: { ambient: 1, diffuse: 0, specular: 0 }, hoverinfo: 'skip', showlegend: false,
    },
    line3(pts, color, 5),
  ]
}

/**
 * A one-color surface. Plotly has no flat-color surface, so the color is given
 * as a two-stop colorscale.
 */
export function surface3(
  g: { X: number[][]; Y: number[][]; Z: number[][] }, color: string, opacity = 0.22,
): Trace {
  return {
    type: 'surface', x: g.X, y: g.Y, z: g.Z, opacity, showscale: false,
    colorscale: [[0, color], [1, color]],
    lighting: { ambient: 0.8, diffuse: 0.35, specular: 0.05, roughness: 1, fresnel: 0 },
    contours: { x: { highlight: false }, y: { highlight: false }, z: { highlight: false } },
    hoverinfo: 'skip', showlegend: false,
  }
}

// ── Labels ──────────────────────────────────────────────────────────────────
// Scene annotations, drawn in the SVG layer above the WebGL scene: the only
// text Plotly underlines, and the vectors of this documentation are underlined
// (Notation page). They stay upright and the same size however the scene turns.

const LABEL_FONT = 'STIX Two Text, Cambria, Times New Roman, serif'

export interface LabelOptions { size?: number; xshift?: number; yshift?: number; xanchor?: string; yanchor?: string }

export function label3(p: Vec3, text: string, color: string, o: LabelOptions = {}): Annotation {
  return {
    x: p[0], y: p[1], z: p[2], text, showarrow: false,
    font: { family: LABEL_FONT, color, size: o.size ?? 17 },
    xanchor: o.xanchor ?? 'center', yanchor: o.yanchor ?? 'middle',
    xshift: o.xshift ?? 0, yshift: o.yshift ?? 0,
  }
}

/** Vector `name` underlined, with an index set as `<sub>` or `<sup>` markup. */
export const vec = (name: string, index = ''): string => `<u><i>${name}</i></u>${index}`

// ── Layouts ─────────────────────────────────────────────────────────────────
// A transparent figure without axes: a picture of a frame has no coordinate to
// read off a grid, the canonical basis is drawn instead. The ranges are fixed,
// so that the picture does not rescale while a slider moves, and `uirevision`
// keeps the camera the reader has chosen across redraws.

export type Bounds = [[number, number], [number, number], [number, number]]

export function layout3(
  pal: Palette, bounds: Bounds, annotations: Annotation[], uid: string,
  eye: Vec3 = [1.22, 0.57, 0.52],
): Record<string, unknown> {
  const span = bounds.map(([a, b]) => b - a)
  const g = Math.cbrt(span[0] * span[1] * span[2])
  const axis = (k: number) => ({ visible: false, range: bounds[k], showspikes: false })
  return {
    paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
    font: { color: pal.ink }, margin: { l: 0, r: 0, t: 0, b: 0 }, showlegend: false,
    modebar: { bgcolor: 'rgba(0,0,0,0)', color: pal.guide, activecolor: pal.ink },
    uirevision: uid,
    scene: {
      bgcolor: 'rgba(0,0,0,0)', hovermode: false,
      aspectmode: 'manual',
      aspectratio: { x: span[0] / g, y: span[1] / g, z: span[2] / g },
      camera: { eye: { x: eye[0], y: eye[1], z: eye[2] }, up: { x: 0, y: 0, z: 1 } },
      xaxis: axis(0), yaxis: axis(1), zaxis: axis(2),
      annotations,
    },
  }
}

// ── 2-D ─────────────────────────────────────────────────────────────────────

export type Vec2 = [number, number]

export function line2(pts: Vec2[], color: string, width = 2.5, dash = 'solid'): Trace {
  return {
    type: 'scatter', mode: 'lines', x: pts.map((p) => p[0]), y: pts.map((p) => p[1]),
    line: { color, width, dash }, hoverinfo: 'skip', showlegend: false,
  }
}

export function fill2(pts: Vec2[], color: string, opacity = 0.3): Trace {
  return {
    type: 'scatter', mode: 'lines', x: pts.map((p) => p[0]), y: pts.map((p) => p[1]),
    fill: 'toself', fillcolor: color, opacity, line: { width: 0, color },
    hoverinfo: 'skip', showlegend: false,
  }
}

export function point2(p: Vec2, color: string, size = 8): Trace {
  return {
    type: 'scatter', mode: 'markers', x: [p[0]], y: [p[1]], marker: { color, size },
    hoverinfo: 'skip', showlegend: false,
  }
}

/** An arrow from `o` to `o + d`, as an annotation: Plotly draws its head. */
export function arrow2(o: Vec2, d: Vec2, color: string, width = 2.5): Annotation {
  return {
    x: o[0] + d[0], y: o[1] + d[1], ax: o[0], ay: o[1],
    xref: 'x', yref: 'y', axref: 'x', ayref: 'y', text: '', showarrow: true,
    arrowhead: 2, arrowsize: 1.1, arrowwidth: width, arrowcolor: color,
  }
}

export function label2(p: Vec2, text: string, color: string, o: LabelOptions = {}): Annotation {
  return {
    x: p[0], y: p[1], xref: 'x', yref: 'y', text, showarrow: false,
    font: { family: LABEL_FONT, color, size: o.size ?? 17 },
    xanchor: o.xanchor ?? 'center', yanchor: o.yanchor ?? 'middle',
    xshift: o.xshift ?? 0, yshift: o.yshift ?? 0,
  }
}

export function layout2(
  pal: Palette, xr: [number, number], yr: [number, number], annotations: Annotation[], uid: string,
): Record<string, unknown> {
  const axis = (r: [number, number]) => ({
    visible: false, range: r, showgrid: false, zeroline: false, fixedrange: false,
  })
  return {
    paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
    font: { color: pal.ink }, margin: { l: 0, r: 0, t: 0, b: 0 }, showlegend: false,
    modebar: { bgcolor: 'rgba(0,0,0,0)', color: pal.guide, activecolor: pal.ink },
    uirevision: uid, hovermode: false, dragmode: 'pan',
    xaxis: axis(xr), yaxis: { ...axis(yr), scaleanchor: 'x', scaleratio: 1 },
    annotations,
  }
}
