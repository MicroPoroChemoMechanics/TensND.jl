// ── Interactive 3D figures ──────────────────────────────────────────────────
// A tutorial emits, for each figure, a `<div>` plus an inline `<script>`
// calling `window.tensndPlotly`, which is `drawPlotly` below; the writer
// wraps that output in `<ClientOnly v-exec-scripts>`, so the call happens on
// mount and again after every client-side navigation.
//
// Plotly is a bundled dependency (`plotly.js-gl3d-dist-min`), imported lazily so
// that only pages carrying a figure pay for it, and awaited through a single
// shared promise so that a page with fifteen scenes still loads it once.
let plotlyModule: Promise<any> | null = null

async function loadPlotly(): Promise<any> {
  if (!plotlyModule) plotlyModule = import('plotly.js-gl3d-dist-min')
  // The gl3d bundle is CommonJS; depending on how Vite interops it, the object
  // is either the module's default export or the module namespace itself.
  const mod = (await plotlyModule) as any
  return mod.default ?? mod
}

// A figure may give any color as `var(--token)`, so that a sketch on a
// transparent background stays legible in both color schemes. Plotly cannot
// read a CSS variable (most of a 3D scene is drawn in WebGL), so each token is
// replaced by its current value before drawing, and the figure is drawn again
// with the new values whenever the reader toggles the scheme. Only figures that
// use a token are kept for that second pass.
const TOKEN = /^var\((--[\w-]+)\)$/
const themed = new Map<string, { data: unknown[]; layout: unknown }>()

function resolveTokens(x: unknown, style: CSSStyleDeclaration, seen: { token: boolean }): any {
  if (typeof x === 'string') {
    const m = TOKEN.exec(x)
    if (!m) return x
    seen.token = true
    return style.getPropertyValue(m[1]).trim() || x
  }
  if (Array.isArray(x)) return x.map((v) => resolveTokens(v, style, seen))
  if (x !== null && typeof x === 'object') {
    const out: Record<string, unknown> = {}
    for (const [k, v] of Object.entries(x)) out[k] = resolveTokens(v, style, seen)
    return out
  }
  return x
}

let themeObserver: MutationObserver | null = null

// VitePress switches the scheme by toggling the `dark` class on `<html>`.
function watchColorScheme(): void {
  if (themeObserver) return
  const root = document.documentElement
  let dark = root.classList.contains('dark')
  themeObserver = new MutationObserver(() => {
    if (root.classList.contains('dark') === dark) return
    dark = !dark
    void redrawThemed()
  })
  themeObserver.observe(root, { attributes: true, attributeFilter: ['class'] })
}

// Writes into `live`, the figure as it now stands, the current value of every
// token of `spec`, the figure as it was emitted. Everything else in `live` is
// left alone: the slider position, the camera and the zoom the reader has set.
function overlayTokens(spec: unknown, live: any, style: CSSStyleDeclaration): void {
  if (spec === null || typeof spec !== 'object' || live === null || typeof live !== 'object') return
  for (const [k, v] of Object.entries(spec as Record<string, unknown>)) {
    if (typeof v === 'string') {
      const m = TOKEN.exec(v)
      if (m) live[k] = style.getPropertyValue(m[1]).trim() || v
    } else if (k in live) {
      overlayTokens(v, live[k], style)
    }
  }
}

async function redrawThemed(): Promise<void> {
  const Plotly = await loadPlotly()
  const style = getComputedStyle(document.documentElement)
  for (const [uid, spec] of themed) {
    const gd = document.getElementById(uid) as any
    if (!gd || !gd.data) {
      themed.delete(uid) // left behind by a client-side navigation
      continue
    }
    // `react` wants new objects to diff against the ones it drew.
    const data = JSON.parse(JSON.stringify(gd.data))
    const layout = JSON.parse(JSON.stringify(gd.layout))
    overlayTokens(spec.data, data, style)
    overlayTokens(spec.layout, layout, style)
    await Plotly.react(gd, data, layout)
  }
}

export async function drawPlotly(uid: string, data: unknown[], layout: unknown): Promise<void> {
  if (!document.getElementById(uid)) return // navigated away before the import resolved
  const Plotly = await loadPlotly()
  // The element may have been torn down while the bundle was in flight.
  if (!document.getElementById(uid)) return
  const seen = { token: false }
  const style = getComputedStyle(document.documentElement)
  const resolvedData = resolveTokens(data, style, seen)
  const resolvedLayout = resolveTokens(layout, style, seen)
  if (seen.token) {
    themed.set(uid, { data, layout })
    watchColorScheme()
  }
  await Plotly.newPlot(uid, resolvedData, resolvedLayout, { responsive: true })
}
