<script setup lang="ts">
// An interactive figure of a frame: Euler angles or a predefined coordinate
// system (`kind`, a key of FRAMES in ./figures/frames.ts), with one slider per
// angle or coordinate and, for a coordinate system, a checkbox showing the
// surfaces on which each coordinate is constant.
//
// Used from a page as
//
//     ```@raw html
//     <ClientOnly><FrameFigure kind="spherical" /></ClientOnly>
//     ```
//
// `ClientOnly` because the figure is WebGL and is drawn after mounting; Plotly
// is imported then, and only on a page that has a figure.
import { computed, onBeforeUnmount, onMounted, reactive, ref, watch } from 'vue'
import { FRAMES } from './figures/frames'
import { ANGLE_DIVISIONS, PLOTLY_CONFIG, loadPlotly, onSchemeChange, palette, piFraction } from './figures/kit'

const props = defineProps<{ kind: string }>()

const spec = FRAMES[props.kind]
if (!spec) throw new Error(`FrameFigure: unknown kind "${props.kind}"`)

const values = reactive(spec.sliders.map((s) => s.init))
const iso = ref(false)
const plot = ref<HTMLDivElement | null>(null)
// Unique on the page: the id Plotly keeps the camera under across redraws.
const uid = `tensnd-frame-${props.kind}-${Math.random().toString(36).slice(2, 8)}`

let Plotly: any = null
let stopWatchingScheme: (() => void) | null = null
let queued = false

function format(k: number): string {
  const s = spec.sliders[k]
  return s.angle ? piFraction(values[k]) : values[k].toFixed(2)
}

async function draw(): Promise<void> {
  if (!Plotly || !plot.value) return
  const q = values.map((v, k) => (spec.sliders[k].angle ? (v * Math.PI) / ANGLE_DIVISIONS : v))
  const { data, layout } = spec.build(q, palette(), iso.value, uid)
  await Plotly.react(plot.value, data, layout, PLOTLY_CONFIG)
}

// A slider fires many events per frame; draw at most once per frame.
function schedule(): void {
  if (queued) return
  queued = true
  requestAnimationFrame(() => {
    queued = false
    void draw()
  })
}

watch([values, iso], schedule)

onMounted(async () => {
  Plotly = await loadPlotly()
  await draw()
  stopWatchingScheme = onSchemeChange(() => void draw())
})

onBeforeUnmount(() => {
  stopWatchingScheme?.()
  if (Plotly && plot.value) Plotly.purge(plot.value)
})

// Plotly fits a 3D scene to the height of its box, so on a narrow screen a
// fixed height would crop the sides: the box never gets much taller than wide.
const height = computed(() => `min(${spec.height}px, 92vw)`)
</script>

<template>
  <figure class="tensnd-frame-figure">
    <div ref="plot" class="tensnd-frame-plot" :style="{ height }"></div>
    <div class="tensnd-frame-controls">
      <label v-for="(s, k) in spec.sliders" :key="k" class="tensnd-frame-slider">
        <span class="tensnd-frame-symbol" :style="{ color: `var(--tensnd-fig-q${k + 1})` }" v-html="s.symbol"></span>
        <input
          v-model.number="values[k]"
          type="range"
          :min="s.min"
          :max="s.max"
          :step="s.step"
          :style="{ accentColor: `var(--tensnd-fig-q${k + 1})` }"
        />
        <output>{{ format(k) }}</output>
      </label>
      <label v-if="spec.iso" class="tensnd-frame-check">
        <input v-model="iso" type="checkbox" />
        {{ spec.iso }}
      </label>
    </div>
  </figure>
</template>
