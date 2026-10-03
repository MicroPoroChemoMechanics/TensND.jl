// Which symbols of the nomenclature a formula holds, read from its TeX source.
//
// The entries come from `nomenclature.json`, which make.jl writes from
// docs/nomenclature.toml. The plugin that typesets the formulas calls
// `symbolsOf(tex, page)` on each one and stores the ids it returns on the
// rendered equation; the theme shows those entries when the equation is hovered.
// Plain JavaScript, so that node can test it without the build.
//
// Ported from MeanFieldHomogenization, which took it from ChemistryLab with three
// changes these formulas need as well: an upright label is a token of its own
// (`\mathrm{GRAD}` is not the letters G, R, A, D), a subscript and a superscript
// are told apart (`\underline{e}_i` is not `\underline{e}^i`), and an underline
// makes one symbol of its argument, as a vector is written `\underline{n}`.

// Commands whose argument is skipped: text that is not a symbol.
const DROPPED = new Set(['label', 'tag', 'mbox', 'textbf', 'textit'])
// Commands whose argument is a label when it is a single word: `\mathrm{GRAD}`,
// `\mathrm{tr}`. Several words are prose, and are skipped.
const LABELS = new Set(['mathrm', 'text', 'textrm', 'operatorname'])
// Commands that make one symbol of their argument: \dot{B}, \mathbb{C}.
const ACCENTS = new Set([
  'mathcal', 'mathbf', 'boldsymbol', 'mathit', 'mathsf', 'mathbb', 'dot', 'ddot', 'bar',
  'hat', 'tilde', 'vec', 'overline', 'underline', 'mathring', 'widehat', 'widetilde',
])

// The next atom of `s` from `i`: a braced group (its content), a command, or a
// character. Returns [text, next index].
function atom(s, i) {
  while (i < s.length && /\s/.test(s[i])) i++
  if (i >= s.length) return ['', i]
  if (s[i] === '{') {
    let depth = 0
    let j = i
    for (; j < s.length; j++) {
      if (s[j] === '{') depth++
      else if (s[j] === '}') { depth--; if (depth === 0) break }
    }
    return [s.slice(i + 1, j), j + 1]
  }
  if (s[i] === '\\') {
    const m = /^\\([A-Za-z]+|.)/.exec(s.slice(i))
    return [m[0], i + m[0].length]
  }
  return [s[i], i + 1]
}

/**
 * The symbols of a TeX string, in order: `{ t, script }`, `t` a command
 * (`\sigma`), an accented symbol (`\mathbb{C}`), a label (`\mathrm{tr}`) or a
 * character, `script` false outside any script and 'sub' or 'sup' inside one.
 */
export function tokenize(tex, inScript = false, out = []) {
  let i = 0
  while (i < tex.length) {
    const c = tex[i]
    if (/\s/.test(c) || c === '{' || c === '}' || c === '&' || c === "'") { i++; continue }
    if (c === '_' || c === '^') {
      const [a, j] = atom(tex, i + 1)
      tokenize(a, c === '_' ? 'sub' : 'sup', out)
      i = j
      continue
    }
    if (c === '\\') {
      const m = /^\\([A-Za-z]+)/.exec(tex.slice(i))
      if (!m) { i += 2; continue } // \, \; \! \{ and the like
      const name = m[1]
      i += m[0].length
      if (DROPPED.has(name)) { const [, j] = atom(tex, i); i = j; continue }
      if (LABELS.has(name)) {
        const [a, j] = atom(tex, i)
        const word = a.trim()
        if (/^[A-Za-z0-9]+$/.test(word)) out.push({ t: `\\mathrm{${word}}`, script: inScript })
        i = j
        continue
      }
      if (ACCENTS.has(name)) {
        const [a, j] = atom(tex, i)
        out.push({ t: `\\${name}{${a.replace(/\s+/g, '')}}`, script: inScript })
        i = j
        continue
      }
      out.push({ t: `\\${name}`, script: inScript })
      continue
    }
    out.push({ t: c, script: inScript })
    i++
  }
  return out
}

function holdsOn(entry, page) {
  return entry.pages.some((p) => (p.endsWith('/') ? page.startsWith(p) : page === p))
}

/**
 * The entries meant on `page` (a path under docs/src, such as
 * "theory/hill_tensors.md"): for each TeX form, the entries scoped to the page if
 * there are any, the unscoped one otherwise. A subscript and a symbol written
 * alike are two forms: the index `i` is not the imaginary unit.
 */
export function entriesFor(entries, page) {
  const byTex = new Map()
  for (const e of entries) {
    const key = `${e.script ? 'script:' : ''}${e.tex}`
    if (!byTex.has(key)) byTex.set(key, [])
    byTex.get(key).push(e)
  }
  const chosen = []
  for (const group of byTex.values()) {
    const scoped = group.filter((e) => e.pages.length > 0 && holdsOn(e, page))
    chosen.push(...(scoped.length > 0 ? scoped : group.filter((e) => e.pages.length === 0)))
  }
  return chosen
}

// An index: `\underline{e}_i` and `\underline{e}_j` are the same entry, and an
// index the entry does not name does not hide the rest of the symbol, so
// `\Gamma^k_{ji}` is the entry written `\Gamma^k_{ij}`.
const isIndex = (tok) => tok.script === 'sub' && /^[ijkIJK]$/.test(tok.t)

/**
 * The tokens of `toks` from `p` that the entry sequence `seq` matches, or null.
 * Index subscripts match one another, and one that the entry does not name is
 * skipped once the symbol has started.
 */
function matchAt(toks, taken, p, seq, script) {
  const used = []
  let q = p
  for (let k = 0; k < seq.length; k++) {
    for (;;) {
      if (q >= toks.length || taken[q]) return null
      const tok = toks[q]
      const sameScript = k === 0 && script ? tok.script !== false : tok.script === seq[k].script
      if (sameScript && (tok.t === seq[k].t || (isIndex(tok) && isIndex(seq[k])))) break
      if (k > 0 && isIndex(tok)) { used.push(q); q++; continue }
      return null
    }
    used.push(q)
    q++
  }
  return used
}

/**
 * The ids of the entries a formula holds, in the order they first appear. A
 * longer symbol is matched first and its tokens are then taken: the `\Gamma` of
 * `\Gamma^k_{ij}` is not also the array of connection coefficients. An entry
 * marked `script` matches only inside a subscript or a superscript, any other
 * one only outside.
 */
export function symbolsOf(tex, page, entries) {
  const toks = tokenize(tex)
  const taken = new Array(toks.length).fill(false)
  const found = []
  const candidates = entriesFor(entries, page)
    .map((e) => ({ e, seq: tokenize(e.tex) }))
    .filter(({ seq }) => seq.length > 0)
    .sort((a, b) => b.seq.length - a.seq.length)
  for (const { e, seq } of candidates) {
    for (let p = 0; p < toks.length; p++) {
      const used = matchAt(toks, taken, p, seq, e.script)
      if (!used) continue
      for (const q of used) taken[q] = true
      found.push({ id: e.id, at: p })
    }
  }
  const seen = new Set()
  return found
    .sort((a, b) => a.at - b.at)
    .map((f) => f.id)
    .filter((id) => (seen.has(id) ? false : (seen.add(id), true)))
}
