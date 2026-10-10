# Bitwise lock of the structured and generic tensor algebra: every value of the
# battery in `golden/cases.jl` equal to the bit to `golden/values.txt`. A
# change that reorders floating-point operations moves values in their last
# bits and fails here on purpose; the largest relative difference is reported,
# and the record is regenerated (`golden/generate.jl`) by a commit that says so.
include(joinpath(@__DIR__, "golden", "cases.jl"))

@testsection "Bitwise lock" begin
    recorded = Dict(
        first(split(l, ' '; limit = 2)) => last(split(l, ' '; limit = 2))
            for l in eachline(joinpath(@__DIR__, "golden", "values.txt"))
    )
    parse_hex(s) = [reinterpret(Float64, parse(UInt64, h; base = 16)) for h in split(s, ',')]
    cases = golden_cases()
    @test Set(first.(cases)) == Set(keys(recorded))
    for (name, f) in cases
        now = golden_record(f)
        then = get(recorded, name, "")
        if now != then && !startswith(now, "ERROR") && !startswith(then, "ERROR")
            a, b = parse_hex(now), parse_hex(then)
            gap = length(a) == length(b) ? maximum(abs.(a .- b) ./ max.(abs.(b), 1.0e-300)) : Inf
            @info "golden case $name moved" max_relative_difference = gap
        end
        @test now == then
    end
end
