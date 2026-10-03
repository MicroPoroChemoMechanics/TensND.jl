# [Nomenclature](@id nomenclature)

This page lists the symbols of the formulas of the documentation, grouped by
subject. A symbol that means different things in different chapters is listed
once for each meaning, with the pages where that meaning holds; a meaning listed
without pages holds everywhere else. Hovering an equation on any page shows the
symbols it holds, with their meaning on that page.

Only the symbols that recur are listed. Following the rule of
[Notation and conventions](@ref th-notation), a symbol used on one page only is
defined on that page, and the typography — the typeface that carries the order
of a tensor, the position of an index, the products and their Julia operators —
is set out there.

```@eval
using TensND
include(joinpath(pkgdir(TensND), "docs", "nomenclature.jl"))
nomenclature_markdown()
```
