```@meta
EditURL = "../../../../scripts/31_performance.jl"
```

# Performance of the structured types

`TensND` is symbolic *and* numerical. The structured types
([`TensISO`](@ref), [`TensTI`](@ref), [`TensOrtho`](@ref)) exist so that
products and inverses are computed on a handful of stored coefficients rather
than on ``3^4=81`` components. This tutorial measures what that buys, and
where it stops applying.

Theory: [Isotropic tensors](@ref th-isotropic),
[The Walpole basis](@ref th-walpole), [Orthotropy](@ref th-orthotropy).

````@example performance
using TensND
using LinearAlgebra
using BenchmarkTools
using Printf

arr(x) = get_array(x)
````

## What each type stores

The whole point in one table: the number of scalars held, against the 81
components of a dense order-4 tensor in 3-D.

````@example performance
ℬ = CanonicalBasis{3, Float64}()
n = (0.0, 0.0, 1.0)

iso = TensISO{3}(3 * 20.0, 2 * 8.0)
ti = TensTI{4}(12.0, 13.0, 3.5, 7.0, 4.0, n)
ort = TensOrtho(10.0, 8.0, 9.0, 3.0, 2.0, 4.0, 2.5, 3.0, 1.5, ℬ)
dense = Tens(arr(ort))

for (name, t) in ("TensISO" => iso, "TensTI{4,5}" => ti, "TensOrtho" => ort)
    @printf "  %-12s stores %2d scalars   (dense: 81 components)\n" name length(get_data(t))
end
````

## Double contraction

Each structured type has a closed-form `⊡` that stays in the stored
representation; the dense route expands all 81 components through an
`einsum`.

````@example performance
function bench_dcontract()
    rows = Any[]
    for (name, a, b) in (
            ("TensISO ⊡ TensISO", iso, TensISO{3}(3 * 30.0, 2 * 12.0)),
            ("TensTI ⊡ TensTI", ti, TensTI{4}(10.0, 9.0, 2.0, 5.0, 3.0, n)),
            ("TensOrtho ⊡ TensOrtho", ort, TensOrtho(6.0, 5.0, 7.0, 2.0, 1.0, 3.0, 1.2, 2.2, 0.8, ℬ)),
        )
        ts = @belapsed $a ⊡ $b
        aa, bb = arr(a), arr(b)
        td = @belapsed $aa ⊡ $bb
        push!(rows, (name, ts, td))
    end
    return rows
end

println("  operation                     structured        dense        speed-up")
for (name, ts, td) in bench_dcontract()
    @printf "  %-26s %8.1f ns   %10.1f ns   ×%6.1f\n" name 1.0e9 * ts 1.0e9 * td td / ts
end
````

The isotropic case is two scalar products; the Walpole case a ``2\times2``
matrix product plus two scalar products; the orthotropic case one
``3\times3`` product plus three scalar products.

## Inversion

The same story, and a starker one: inverting a dense minor-symmetric order-4
tensor means inverting a ``6\times6`` matrix.

````@example performance
function bench_inv()
    rows = Any[]
    for (name, t) in ("TensISO" => iso, "TensTI{4,5}" => ti, "TensOrtho" => ort)
        ts = @belapsed inv($t)
        km = Matrix(KM(t))
        td = @belapsed inv($km)
        push!(rows, (name, ts, td))
    end
    return rows
end

println("  type            structured        6×6 inverse     speed-up")
for (name, ts, td) in bench_inv()
    @printf "  %-14s %8.1f ns   %12.1f ns   ×%6.1f\n" name 1.0e9 * ts 1.0e9 * td td / ts
end
````

## Allocations

A structured operation allocates a single small object — the tuple of stored
coefficients — where the dense route allocates intermediate arrays of 81
components. The figure below is the *whole* cost of the result, and it scales
with the number of stored scalars, not with ``3^4``.

````@example performance
for (name, f) in (
        "TensISO ⊡ TensISO" => () -> iso ⊡ iso,
        "inv(TensISO)" => () -> inv(iso),
        "TensTI ⊡ TensTI" => () -> ti ⊡ ti,
        "inv(TensTI)" => () -> inv(ti),
        "inv(TensOrtho)" => () -> inv(ort),
    )
    b = @benchmark $f()
    @printf "  %-22s %6d bytes in %d allocations\n" name b.memory b.allocs
end
````

## Element access

A structured type has no component array to index into, so `getindex`
reconstructs the component from the stored coefficients. That is cheap, but
**repeated** indexing is not the way to use these types: if a whole array is
needed, ask for it once with [`get_array`](@ref).

````@example performance
@printf "  getindex on TensOrtho : %6.1f ns\n" 1.0e9 * @belapsed $ort[1, 2, 1, 2]
@printf "  getindex on a dense   : %6.1f ns\n" 1.0e9 * @belapsed $(arr(ort))[1, 2, 1, 2]
@printf "  get_array(TensOrtho)  : %6.1f ns  (once, then index freely)\n" 1.0e9 * @belapsed get_array($ort)
````

## Where the structure stops paying

The closed forms apply only while the result **stays in the class**. Two cases
leave it, and both fall back to the generic route:

- a product of two `TensOrtho` expressed in *different* material frames is
  generally fully anisotropic;
- a product of two orthotropic tensors in the *same* frame is orthotropic but
  **not major-symmetric** — twelve constants where `TensOrtho` stores nine —
  so the result is returned as a generic tensor
  ([Orthotropy](@ref th-orthotropy)).

````@example performance
ort2 = TensOrtho(6.0, 5.0, 7.0, 2.0, 1.0, 3.0, 1.2, 2.2, 0.8, ℬ)
ort_rot = TensOrtho(6.0, 5.0, 7.0, 2.0, 1.0, 3.0, 1.2, 2.2, 0.8, Basis(0.3, 0.7, 0.2))

println("same frame      : ", typeof(ort ⊡ ort2))
println("different frames: ", typeof(ort ⊡ ort_rot))
````

Both agree with the dense computation to machine precision — the fallback is
exact, not approximate:

````@example performance
(
    norm(arr(ort ⊡ ort2) - arr(ort) ⊡ arr(ort2)),
    norm(arr(ort ⊡ ort_rot) - arr(ort) ⊡ arr(ort_rot)),
)
````

## Mixed-class arithmetic

Operations between different classes promote to the least constrained of the
two, in closed form where possible
([`iso_to_ortho`](@ref), [`walpole_to_ortho`](@ref)):

````@example performance
for (name, x, y) in (
        ("TensISO ⊡ TensTI", iso, ti),
        ("TensISO ⊡ TensOrtho", iso, ort),
        ("TensTI ⊡ TensISO", ti, iso),
    )
    @printf "  %-22s → %s\n" name string(typeof(x ⊡ y))
end
````

## Summary

| | stored scalars | `⊡` | `inv` |
|:--|:--|:--|:--|
| [`TensISO`](@ref) | 2 | two scalar products | two reciprocals |
| [`TensTI`](@ref) `{4,5}` | 5 | ``2\times2`` product + 2 scalars | ``2\times2`` inverse + 2 reciprocals |
| [`TensOrtho`](@ref) | 9 | ``3\times3`` product + 3 scalars | ``3\times3`` inverse + 3 reciprocals |
| dense | 81 | `einsum` over 81 components | ``6\times6`` inverse |

Use a structured type whenever the physics guarantees the symmetry. Convert to
a dense array once, with [`get_array`](@ref), when component-by-component
access is genuinely needed.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

