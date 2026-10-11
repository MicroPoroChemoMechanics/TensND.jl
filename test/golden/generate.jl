# Writes `values.txt`, the record the bitwise lock compares against:
#     julia --project=<env with TensND> test/golden/generate.jl
# Regenerated only on purpose, by a commit that says which values moved and why.
include(joinpath(@__DIR__, "cases.jl"))
open(joinpath(@__DIR__, "values.txt"), "w") do io
    for (name, f) in golden_cases()
        println(io, name, " ", golden_record(f))
    end
end
