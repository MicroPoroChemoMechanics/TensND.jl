# Lock of the structured and generic tensor algebra: every value of the battery
# in `golden/cases.jl` against `golden/values.txt`, recorded in hex.
#
# On the machine that recorded the values every case is equal to the bit, and
# `TENSND_GOLDEN_BITWISE=true` checks exactly that: a refactoring that reorders
# floating-point operations moves values in their last bits and fails on
# purpose. Elsewhere a BLAS reduction — the 81-term `dot` behind `⊙`, hence
# behind `best_fit_iso` — sums in an order set by the width of the processor's
# vector registers, so the same code differs by an ulp or two; the comparison
# then allows 8 ulps of the scale of the case, its largest component. An error
# recorded or raised, and a change of length, are always compared exactly. A
# case that moves is reported with its largest relative difference, and the
# record is regenerated (`golden/generate.jl`) by a commit that says so.
include(joinpath(@__DIR__, "golden", "cases.jl"))

const GOLDEN_BITWISE = get(ENV, "TENSND_GOLDEN_BITWISE", "false") == "true"

_golden_values(s) = [reinterpret(Float64, parse(UInt64, h; base = 16)) for h in split(s, ',')]

function _golden_close(now, then)
    now == then && return true
    (GOLDEN_BITWISE || startswith(now, "ERROR") || startswith(then, "ERROR")) && return false
    a, b = _golden_values(now), _golden_values(then)
    return length(a) == length(b) && maximum(abs.(a .- b)) <= 8 * eps(maximum(abs, b))
end

@testsection "Bitwise lock" begin
    recorded = Dict(
        first(split(l, ' '; limit = 2)) => last(split(l, ' '; limit = 2))
            for l in eachline(joinpath(@__DIR__, "golden", "values.txt"))
    )
    cases = golden_cases()
    @test Set(first.(cases)) == Set(keys(recorded))
    for (name, f) in cases
        now = golden_record(f)
        then = get(recorded, name, "")
        close = _golden_close(now, then)
        if !close && !startswith(now, "ERROR") && !startswith(then, "ERROR")
            a, b = _golden_values(now), _golden_values(then)
            gap = length(a) == length(b) ? maximum(abs.(a .- b) ./ max.(abs.(b), 1.0e-300)) : Inf
            @info "golden case $name moved" max_relative_difference = gap
        end
        @test close
    end
end
