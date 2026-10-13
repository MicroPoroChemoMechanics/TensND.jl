<script setup lang="ts">
// An interactive figure of a frame: Euler angles or a predefined coordinate
// system (`kind`, a key of FRAMES in ./figures/frames.ts), with one slider per
// angle, coordinate or parameter, radio buttons when the figure has several
// forms (prolate or oblate spheroidal coordinates) and, for a coordinate system,
// a checkbox showing the surfaces on which each coordinate is constant.
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
const variant = ref(0)
// The sliders of the current form: a form may move the bounds of some of them.
const sliders = computed(() =>
  spec.sliders.map((s, k) => ({ ...s, ...(spec.variants?.[variant.value]?.sliders?.[k] ?? {}) })),
)
const sliderColor = (k: number) => sliders.value[k].color ?? `var(--tensnd-fig-q${k + 1})`
const plot = ref<HTMLDivElement | null>(null)
// Unique on the page: the id Plotly keeps the camera under across redraws.
const uid = `tensnd-frame-${props.kind}-${Math.random().toString(36).slice(2, 8)}`

let Plotly: any = null
let stopWatchingScheme: (() => void) | null = null
let queued = false

function format(k: number): string {
  const s = sliders.value[k]
  // A true minus sign, as in the fractions of π.
  return s.angle ? piFraction(values[k]) : values[k].toFixed(2).replace('-', '−')
}

async function draw(): Promise<void> {
  if (!Plotly || !plot.value) return
  const q = values.map((v, k) => (sliders.value[k].angle ? (v * Math.PI) / ANGLE_DIVISIONS : v))
  const { data, layout } = spec.build(q, palette(), iso.value, uid, variant.value)
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

// A value the new form does not allow is brought back to its nearest bound.
watch(variant, () => {
  sliders.value.forEach((s, k) => {
    values[k] = Math.min(s.max, Math.max(s.min, values[k]))
  })
})

watch([values, iso, variant], schedule)

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
      <label v-for="(s, k) in sliders" :key="k" class="tensnd-frame-slider">
        <span class="tensnd-frame-symbol" :style="{ color: sliderColor(k) }" v-html="s.symbol"></span>
        <input
          v-model.number="values[k]"
          type="range"
          :min="s.min"
          :max="s.max"
          :step="s.step"
          :style="{ accentColor: sliderColor(k) }"
        />
        <output>{{ format(k) }}</output>
      </label>
      <span v-if="spec.variants" class="tensnd-frame-variants" role="radiogroup">
        <label v-for="(v, i) in spec.variants" :key="i" class="tensnd-frame-check">
          <input v-model="variant" type="radio" :name="uid" :value="i" />
          {{ v.label }}
        </label>
      </span>
      <label v-if="spec.iso" class="tensnd-frame-check">
        <input v-model="iso" type="checkbox" />
        {{ spec.iso }}
      </label>
    </div>
  </figure>
</template>
