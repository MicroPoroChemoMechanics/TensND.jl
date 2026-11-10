# `docs/assets_src` — the generators of the drawn figures

Every figure that is *drawn* rather than plotted from data lives here as the
program that produces it, and its output is committed under `docs/src/assets/`.
The documentation build never runs these scripts, so `Luxor` and `Typstry` stay
out of `docs/Project.toml`. The conventions are those of the figures of
MeanFieldHomogenization.jl, so that the two packages look alike.

- Drawing: [Luxor](https://github.com/JuliaGraphics/Luxor.jl).
- Labels: [Typst](https://typst.app), through
  [Typstry](https://github.com/jakobjpeters/Typstry.jl) — mathematics in New
  Computer Modern Math with the notation of the documentation (underlined
  vectors, italic coordinates), words in DejaVu Sans. Typst's own SVG is written
  into each sheet, every glyph defined once, so the files are self-contained.

| file | output under `docs/src/assets/` |
| :--- | :--- |
| `src/style.jl` | — the palette by role, the sheet and its white panel, `label` and `smart_label`, the primitives |
| `src/coordinates.jl` | `coordinates/`: the meridian plane of the prolate and oblate systems, the three sections of the ellipsoidal chart, the domain of convergence of the re-expansion of exterior spheroidal harmonics |

## Regenerating

```bash
julia --project=docs/assets_src -e 'using Pkg; Pkg.instantiate()'
julia --project=docs/assets_src docs/assets_src/make.jl                 # every family
julia --project=docs/assets_src docs/assets_src/make.jl coordinates     # one family
```

Families: `coordinates`.

Then **look at the result**: `rsvg-convert -z 2 -o out.png in.svg` renders a
sheet at readable size. Overlapping text is the failure mode of drawings like
these, and nothing in the build reports it. `smart_label` puts a label where
its box keeps the largest clearance from what is already drawn; it does not
replace the look.

## Rules the sheets follow

1. **One palette, by role.** Colors are named in `style.jl` by what they stand
   for (`COORD_Q` for `q`, `τ`, `ρ`; `COORD_P` for `p`, `ν`; `COORD_MU` for `μ`)
   and never written as hex codes in the figure files.
2. **Every sheet carries its own white panel.** An SVG served as an image
   follows the operating system's color scheme, not the site's theme toggle, so
   it must not rely on the page background.
3. **The symbols are the ones on the page.** A drawing carries the notation of
   the formulas it illustrates, and the figure checks them: the generator asserts
   that the point `M` lies on the drawn curves and that its two unit vectors are
   orthonormal and normal to them.
