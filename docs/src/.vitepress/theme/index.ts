// .vitepress/theme/index.ts
import { h } from 'vue'
import DefaultTheme from 'vitepress/theme'
import type { Theme as ThemeConfig } from 'vitepress'
import 'virtual:mathjax-styles.css';

import { 
  NolebaseEnhancedReadabilitiesMenu, 
  NolebaseEnhancedReadabilitiesScreenMenu, 
} from '@nolebase/vitepress-plugin-enhanced-readabilities/client'

import VersionPicker from "@/VersionPicker.vue"
import AuthorBadge from '@/AuthorBadge.vue'
import Authors from '@/Authors.vue'
import SidebarDrawerToggle from '@/SidebarDrawerToggle.vue'
// __DV_PLUGIN_COMPONENT_IMPORTS__

import { enhanceAppWithTabs } from 'vitepress-plugin-tabs/client'

import '@nolebase/vitepress-plugin-enhanced-readabilities/client/style.css'
import './style.css' // You could setup your own, or else a default will be copied.
import './docstrings.css' // You could setup your own, or else a default will be copied.
import './overrides.css' // You could setup your own, or else a default will be copied.

// `v-exec-scripts` runs the <script> tags inside a `v-html`'d block: innerHTML never executes
// scripts, so we re-create each one. `src` scripts are awaited so order holds (bundle before
// its callers). Used on interactive text/html output (WGLMakie/Bonito, Plotly) which the writer
// wraps in <ClientOnly> + this directive.
async function activateScripts(container: Element): Promise<void> {
  for (const old of Array.from(container.querySelectorAll('script'))) {
    const fresh = document.createElement('script')
    for (const attr of Array.from(old.attributes)) fresh.setAttribute(attr.name, attr.value)
    fresh.textContent = old.textContent
    const hasSrc = old.hasAttribute('src')
    const ran = hasSrc
      ? new Promise<void>((resolve) => {
          fresh.addEventListener('load', () => resolve(), { once: true })
          fresh.addEventListener('error', () => resolve(), { once: true })
        })
      : Promise.resolve()
    old.replaceWith(fresh) // inserting is what runs it; inline scripts run synchronously
    await ran
  }
}

// ── Citation hints ──────────────────────────────────────────────────────────
// A numeric citation reads as `[42]`, which tells the reader nothing without a
// trip to the References page. Hovering one pops the full entry instead.
//
// The entries are read from the References page itself, fetched once on the
// first hover and indexed by the anchor DocumenterCitations puts on each item,
// so there is nothing to keep in sync: whatever the bibliography renders is
// what the hint shows.
let bibIndex: Promise<Map<string, string>> | null = null

async function loadBibliography(base: string): Promise<Map<string, string>> {
  const map = new Map<string, string>()
  try {
    const html = await (await fetch(`${base}references`)).text()
    const doc = new DOMParser().parseFromString(html, 'text/html')
    for (const a of Array.from(doc.querySelectorAll('a[id]'))) {
      const item = a.closest('li')
      if (item) map.set(a.id, (item.textContent || '').trim().replace(/\s+/g, ' '))
    }
  } catch {
    /* offline or missing page: hints simply stay silent */
  }
  return map
}

function citationKey(a: HTMLAnchorElement): string | null {
  const href = a.getAttribute('href') || ''
  const i = href.indexOf('/references#')
  return i === -1 ? null : decodeURIComponent(href.slice(i + '/references#'.length))
}

function installCitationHints(base: string): void {
  const tip = document.createElement('div')
  tip.className = 'tensnd-cite-hint'
  tip.setAttribute('role', 'tooltip')
  document.body.appendChild(tip)
  let current: HTMLAnchorElement | null = null

  const hide = () => { tip.classList.remove('is-visible'); current = null }

  document.addEventListener('mouseover', async (e) => {
    const a = (e.target as HTMLElement | null)?.closest?.('a[href*="/references#"]') as HTMLAnchorElement | null
    if (!a || a === current) return
    const key = citationKey(a)
    if (!key) return
    current = a
    if (!bibIndex) bibIndex = loadBibliography(base)
    const text = (await bibIndex).get(key)
    if (!text || current !== a) return
    tip.textContent = text
    tip.classList.add('is-visible')
    // Place it under the citation, nudged back inside the viewport if needed.
    const r = a.getBoundingClientRect()
    const w = Math.min(460, window.innerWidth - 24)
    tip.style.width = `${w}px`
    tip.style.left = `${Math.max(12, Math.min(r.left, window.innerWidth - w - 12))}px`
    tip.style.top = `${r.bottom + window.scrollY + 8}px`
  })

  document.addEventListener('mouseout', (e) => {
    const a = (e.target as HTMLElement | null)?.closest?.('a[href*="/references#"]')
    if (a && a === current) hide()
  })
  window.addEventListener('scroll', hide, { passive: true })
}

export const Theme: ThemeConfig = {
  extends: DefaultTheme,
  Layout() {
    return h(DefaultTheme.Layout, null, {
      'nav-bar-content-after': () => [
        h(NolebaseEnhancedReadabilitiesMenu), // Enhanced Readabilities menu
      ],
      // A enhanced readabilities menu for narrower screens (usually smaller than iPad Mini)
      'nav-screen-content-after': () => h(NolebaseEnhancedReadabilitiesScreenMenu),
      // Sidebar drawer toggle button (to the left of search bar)
      'nav-bar-content-before': () => h(SidebarDrawerToggle),
    })
  },
  enhanceApp({ app, router, siteData }) {
    enhanceAppWithTabs(app);
    if (typeof window !== "undefined") installCitationHints(siteData.value.base)
    app.component('VersionPicker', VersionPicker);
    app.component('AuthorBadge', AuthorBadge)
    app.component('Authors', Authors)

    // Execute the scripts inside interactive `text/html` outputs (WGLMakie/Bonito, Plotly, …)
    // once their `<ClientOnly>` wrapper has mounted on the client. `mounted` fires on the
    // initial client render and again whenever the page component is remounted by a
    // client-side navigation, so figures initialise instead of staying blank.
    app.directive('exec-scripts', {
      mounted(el: HTMLElement) {
        activateScripts(el)
      },
    })
    // __DV_PLUGIN_COMPONENT_REGISTRATIONS__
  }
}
export default Theme