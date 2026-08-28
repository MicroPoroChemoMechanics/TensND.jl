import { defineConfig } from 'vitepress'
import { withMermaid } from 'vitepress-plugin-mermaid'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import { mathjaxPlugin } from './mathjax-plugin'
import { juliaReplTransformer } from './julia-repl-transformer'
import footnote from "markdown-it-footnote";
import path from 'path'

const mathjax = mathjaxPlugin()

function getBaseRepository(base: string): string {
  if (!base || base === '/') return '/';
  const parts = base.split('/').filter(Boolean);
  return parts.length > 0 ? `/${parts[0]}/` : '/';
}

const baseTemp = {
  base: 'REPLACE_ME_DOCUMENTER_VITEPRESS',// TODO: replace this in makedocs!
}

const navTemp = {
  nav: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
}

const sidebarTemp = {
  sidebar: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
}

// DocumenterVitepress mirrors the whole `pages` tree into the navbar, which for
// this manual means ten dropdowns holding 158 entries — a duplicate of the
// sidebar that fills the bar edge to edge and pushes the GitHub link off screen.
// Keep the main reading path as dropdowns, fold the reference material into one
// "More" menu, and shorten the two labels that were doing the most damage.
const RENAME: Record<string, string> = {
  'Finite-element coupling': 'FE coupling',
  'Tools and migration': 'Tools',
}
const MORE = ['Developer', 'API', 'References']

// Order is the one declared in `pages` — the navbar must not tell a different
// story from the sidebar. Only two things are done to it: the reference
// material is folded into one menu, and two long labels are shortened.
function curateNav(items: any[]): any[] {
  const more = items.filter((i) => MORE.includes(i.text))
  const out = items
    .filter((i) => !more.includes(i))
    .map((i) => ({ ...i, text: RENAME[i.text] ?? i.text }))
  if (more.length) out.push({ text: 'Reference', items: more })
  return out
}

// VitePress renders a sidebar group expanded unless it says otherwise, and
// DocumenterVitepress hard-codes `collapsed: false`. With six levels of nesting
// the panel opens several screens tall, so every group starts closed instead;
// the group holding the current page still opens on its own.
function collapseGroups(node: any): any {
  if (Array.isArray(node)) return node.map(collapseGroups)
  if (node && typeof node === 'object') {
    const out: any = { ...node }
    if (Array.isArray(out.items)) {
      out.items = out.items.map(collapseGroups)
      out.collapsed = true
    }
    return out
  }
  return node
}

const nav = [
  ...curateNav(navTemp.nav as unknown as any[]),
  {
    component: 'VersionPicker'
  }
]

const sidebar = collapseGroups(sidebarTemp.sidebar as unknown as any)

// https://vitepress.dev/reference/site-config
// `withMermaid` wraps the config so that ```mermaid fences are rendered by
// mermaid itself; it replaces DocumenterMermaid, whose renderer was tied to
// Documenter's own HTML output.
export default withMermaid(defineConfig({
  base: 'REPLACE_ME_DOCUMENTER_VITEPRESS',// TODO: replace this in makedocs!
  title: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  description: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  lastUpdated: true,
  cleanUrls: true,
  outDir: 'REPLACE_ME_DOCUMENTER_VITEPRESS', // This is required for MarkdownVitepress to work correctly...
  head: [
    ['link', { rel: 'icon', href: 'REPLACE_ME_DOCUMENTER_VITEPRESS_FAVICON' }],
    ['script', {src: `${getBaseRepository(baseTemp.base)}versions.js`}],
    // ['script', {src: '/versions.js'], for custom domains, I guess if deploy_url is available.
    ['script', {src: `${baseTemp.base}siteinfo.js`}],
    // REPLACE_ME_DOCUMENTER_VITEPRESS_NOINDEX
  ],
  
  markdown: {
    codeTransformers: [juliaReplTransformer()],
    config(md) {
      md.use(tabsMarkdownPlugin);
      md.use(footnote);
      mathjax.markdownConfig(md);
    },
    theme: {
      light: "github-light",
      dark: "github-dark"
    },
  },
  vite: {
    plugins: [
      mathjax.vitePlugin,
    ],
    define: {
      __DEPLOY_ABSPATH__: JSON.stringify('REPLACE_ME_DOCUMENTER_VITEPRESS_DEPLOY_ABSPATH'),
    },
    resolve: {
      alias: {
        '@': path.resolve(__dirname, '../components')
      }
    },
    optimizeDeps: {
      exclude: [ 
        '@nolebase/vitepress-plugin-enhanced-readabilities/client',
        'vitepress',
        '@nolebase/ui',
      ], 
    }, 
    ssr: { 
      noExternal: [ 
        // If there are other packages that need to be processed by Vite, you can add them here.
        '@nolebase/vitepress-plugin-enhanced-readabilities',
        '@nolebase/ui',
      ], 
    },
  },
  mermaid: {
    // Let the diagram keep its natural size and be constrained by CSS instead;
    // mermaid's own max-width shrinks the font until wide flowcharts are unreadable.
    flowchart: { useMaxWidth: false },
    sequence: { useMaxWidth: false },
  },
  themeConfig: {
    outline: 'deep',
    logo: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
    search: {
      provider: 'local',
      options: {
        detailedView: true
      }
    },
    nav,
    sidebar,
    sidebarDrawer: 'REPLACE_ME_DOCUMENTER_VITEPRESS_SIDEBAR_DRAWER',
    editLink: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
    socialLinks: [
      { icon: 'github', link: 'REPLACE_ME_DOCUMENTER_VITEPRESS' }
    ],
    footer: {
      message: 'Made with <a href="https://luxdl.github.io/DocumenterVitepress.jl/dev/" target="_blank"><strong>DocumenterVitepress.jl</strong></a><br>',
      copyright: `© Copyright ${new Date().getUTCFullYear()}.`
    }
  }
}))
