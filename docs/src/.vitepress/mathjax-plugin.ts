// adapter from https://github.com/orgs/vuepress-theme-hope/discussions/5178#discussioncomment-15642629
// mathjax-plugin.ts
// @ts-ignore
import MathJax from '@mathjax/src'
import type { Plugin as VitePlugin } from 'vite'
import type MarkdownIt from 'markdown-it'
import { tex as mdTex } from '@mdit/plugin-tex'
// The symbols of the nomenclature, written by make.jl from docs/nomenclature.toml.
import nomenclature from './nomenclature.json'
// @ts-ignore
import { symbolsOf } from './nomenclature-match.mjs'

const mathjaxStyleModuleID = 'virtual:mathjax-styles.css'

// Width, in px, of the box MathJax lays math out in at build time: MathJax's
// own default (80 ex of 8 px), stated so that the warning below can name it.
// With `displayOverflow: 'linebreak'`, a wider display equation is cut into
// lines baked into the static SVG, whatever the reader's screen.
const BUILD_WIDTH_PX = 640

interface MathJaxOptions {
  font?: string
}

async function initializeMathJax(options: MathJaxOptions = {}) {
  const font = options.font || 'mathjax-newcm'

  const config: any = {
    loader: {
      load: [
        'input/tex',
        'output/svg',
        '[tex]/boldsymbol',
        '[tex]/braket',
        '[tex]/mathtools',
      ],
      paths: { mathjax: '@mathjax/src/bundle' },
    },
    tex: {
      tags: 'ams',
      packages: {
        '[+]': ['boldsymbol', 'braket', 'mathtools'],
      },
    },
    output: {
      font,
      displayOverflow: 'linebreak',
      mtextInheritFont: true,
    },
    svg: {
      fontCache: 'none', // critical: avoids async font loading
    },
  }

  await MathJax.init(config)
  await MathJax.startup.document.outputJax.font.loadDynamicFiles()
}

export function mathjaxPlugin(options: MathJaxOptions = {}) {
  let adaptor: any
  let initialized = false
  // The page being rendered, a path under docs/src such as
  // "theory/hill_tensors.md", from the environment VitePress passes to
  // `md.render`: a symbol's meaning can depend on the page.
  let currentPage = ''

  async function ensureInitialized() {
    if (!initialized) {
      await initializeMathJax(options)
      adaptor = MathJax.startup.adaptor
      initialized = true
    }
  }

  // Each equation is reported once per process: the dev server renders a page
  // again on every reload, and a build may render it once per bundle.
  const reportedBreaks = new Set<string>()

  // The equation's lines were chosen by MathJax, not by its author: it is wider
  // than the build box, or carries a top-level `\\` (the pages have none; their
  // multi-line equations use `aligned`, `gathered` or `array`). The break usually
  // falls on a relation and leaves an orphan `=` line.
  function warnLineBroken(content: string) {
    const key = `${currentPage}\0${content}`
    if (reportedBreaks.has(key)) return
    reportedBreaks.add(key)
    const tex = content.replace(/\s+/g, ' ').trim()
    console.warn(
      `[mathjax] ${currentPage}: a display equation was broken into lines ` +
        `(wider than the ${BUILD_WIDTH_PX} px build width, or a top-level \\\\); ` +
        `lay it out with aligned or gathered:\n    ` +
        (tex.length > 160 ? `${tex.slice(0, 157)}...` : tex),
    )
  }

  function renderMath(content: string, displayMode: boolean): string {
    if (!initialized) {
      throw new Error('MathJax not initialized')
    }

    const node = MathJax.tex2svg(content, {
      display: displayMode,
      containerWidth: BUILD_WIDTH_PX,
    })

    // Prevent Vue from touching MathJax output
    adaptor.setAttribute(node, 'v-pre', '')

    // The symbols of the nomenclature this formula holds, with their meaning on
    // the page being rendered; the theme shows them when the formula is hovered.
    // The nomenclature page itself lists them already.
    if (currentPage !== 'nomenclature.md') {
      const ids = symbolsOf(content, currentPage, nomenclature)
      if (ids.length > 0) adaptor.setAttribute(node, 'data-nomen', ids.join(' '))
    }

    let html = adaptor.outerHTML(node)

    if (displayMode && html.includes('data-mjx-linestack')) warnLineBroken(content)

    // Preserve spaces inside mjx-break (SVG only)
    html = html.replace(
      /<mjx-break(.*?)>(.*?)<\/mjx-break>/g,
      (_: string, attr: string, inner: string) =>
        `<mjx-break${attr}>${inner.replace(/ /g, '&nbsp;')}</mjx-break>`,
    )
    
    // Wrap only display equations (not inline math)
    html = html.replace(
      /(<mjx-container[^>]*display="true"[^>]*>)([\s\S]*?)(<\/mjx-container>)/,
      '<div class="mjx-scroll-wrapper">$1$2$3</div>'
    )

    return html
  }

  function getMathJaxStyles(): string {
    return initialized
      ? adaptor.textContent(MathJax.svgStylesheet()) || ''
      : ''
  }

  function resetMathJax(): void {
    if (!initialized) return
    MathJax.texReset()
    MathJax.typesetClear()
  }

  function viteMathJax(): VitePlugin {
    const virtualModuleID = '\0' + mathjaxStyleModuleID

    return {
      name: 'mathjax-styles',

      resolveId(id) {
        if (id === mathjaxStyleModuleID) {
          return virtualModuleID
        }
      },

      async load(id) {
        if (id === virtualModuleID) {
          await ensureInitialized()
          return getMathJaxStyles()
        }
      },
    }
  }

  function mdMathJax(md: MarkdownIt): void {
    mdTex(md, {
      render: renderMath,
    })

    const orig = md.render
    md.render = function (...args) {
      const env: any = args[1]
      currentPage = (env && env.relativePath) || ''
      resetMathJax()
      return orig.apply(this, args)
    }
  }

  const init = ensureInitialized()

  return {
    vitePlugin: viteMathJax(),
    markdownConfig: mdMathJax,
    styleModuleID: mathjaxStyleModuleID,
    init,
  }
}