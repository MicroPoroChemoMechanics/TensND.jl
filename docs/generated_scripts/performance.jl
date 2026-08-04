import Pkg
Pkg.activate(joinpath(@__DIR__, "..", "docs"); io = devnull)
# The `docs` environment declares every dependency the tutorials use
# (SymPy, Symbolics, NLopt, BenchmarkTools, ...) and resolves TensND
# itself through `[sources]`. These four lines are stripped from the
# generated page and notebook, which already run inside that project.

using TensND
using LinearAlgebra
using BenchmarkTools
using Printf

arr(x) = get_array(x)

ℬ = CanonicalBasis{3, Float64}()
n = (0.0, 0.0, 1.0)

iso = TensISO{3}(3 * 20.0, 2 * 8.0)
ti = TensTI{4}(12.0, 13.0, 3.5, 7.0, 4.0, n)
ort = TensOrtho(10.0, 8.0, 9.0, 3.0, 2.0, 4.0, 2.5, 3.0, 1.5, ℬ)
dense = Tens(arr(ort))

for (name, t) in ("TensISO" => iso, "TensTI{4,5}" => ti, "TensOrtho" => ort)
    @printf "  %-12s stores %2d scalars   (dense: 81 components)\n" name length(get_data(t))
end

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

@printf "  getindex on TensOrtho : %6.1f ns\n" 1.0e9 * @belapsed $ort[1, 2, 1, 2]
@printf "  getindex on a dense   : %6.1f ns\n" 1.0e9 * @belapsed $(arr(ort))[1, 2, 1, 2]
@printf "  get_array(TensOrtho)  : %6.1f ns  (once, then index freely)\n" 1.0e9 * @belapsed get_array($ort)

ort2 = TensOrtho(6.0, 5.0, 7.0, 2.0, 1.0, 3.0, 1.2, 2.2, 0.8, ℬ)
ort_rot = TensOrtho(6.0, 5.0, 7.0, 2.0, 1.0, 3.0, 1.2, 2.2, 0.8, Basis(0.3, 0.7, 0.2))

println("same frame      : ", typeof(ort ⊡ ort2))
println("different frames: ", typeof(ort ⊡ ort_rot))

(
    norm(arr(ort ⊡ ort2) - arr(ort) ⊡ arr(ort2)),
    norm(arr(ort ⊡ ort_rot) - arr(ort) ⊡ arr(ort_rot)),
)

for (name, x, y) in (
        ("TensISO ⊡ TensTI", iso, ti),
        ("TensISO ⊡ TensOrtho", iso, ort),
        ("TensTI ⊡ TensISO", ti, iso),
    )
    @printf "  %-22s → %s\n" name string(typeof(x ⊡ y))
end

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
