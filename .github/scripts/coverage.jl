# Build `lcov.info` for Codecov.
#
# Why not `julia-actions/julia-processcoverage`?  That action calls
# CoverageTools' `amend_coverage_from_src!`, which marks as "executable" every
# line the *source parser* believes to be code — including lines Julia's
# instrumentation never emits a count for.  Methods generated inside `@eval`
# loops are the main offender, and both TensND and MeanFieldHom use that style
# heavily.  Such lines can never be covered by any test, yet they sit in the
# denominator forever: TensND measures 1352/1418 = 95 % of the lines Julia
# actually instruments, but the amended ratio is 1352/2278 = 59 %.
#
# So we mix the two:
#   • a file that was loaded      → raw `.cov` counts, no amendment (faithful);
#   • a file never loaded at all  → amended, i.e. every executable line counted
#     as missed, so genuine gaps (`src/submanifold.jl`, `ext/…Ext.jl`) stay
#     visible at 0 % instead of silently vanishing from the report.

# CoverageTools is a dependency of Coverage
using Coverage
using Coverage.CoverageTools

const DIRS = ["src", "ext"]
# `.cov` files are also dropped next to the test sources; wipe them too so a
# local run leaves the working tree clean.
const CLEAN_DIRS = ["src", "ext", "test", "docs"]

collect_dirs() = reduce(vcat, (process_folder(d) for d in DIRS if isdir(d)); init = FileCoverage[])

"""
Raw counts for loaded files, amended (all-missed) counts for never-loaded ones.
"""
function merged_coverage()
    raw = withenv("DISABLE_AMEND_COVERAGE_FROM_SRC" => "yes") do
        collect_dirs()
    end
    amended = Dict(fc.filename => fc for fc in collect_dirs())
    out = FileCoverage[]
    for fc in raw
        # `process_cov` returns an all-`nothing` vector when no `.cov` file
        # exists, i.e. when nothing in the file was ever executed.
        push!(out, any(!isnothing, fc.coverage) ? fc : amended[fc.filename])
    end
    return out
end

function report(fcs)
    tot = 0
    cov = 0
    for fc in sort(fcs, by = x -> x.filename)
        lines = count(!isnothing, fc.coverage)
        lines == 0 && continue
        hit = count(x -> !isnothing(x) && x > 0, fc.coverage)
        tot += lines
        cov += hit
        println(rpad(fc.filename, 46), lpad(hit, 5), "/", lpad(lines, 5), "  ", round(100 * hit / lines, digits = 1), "%")
    end
    println("TOTAL ", cov, "/", tot, " = ", round(100 * cov / tot, digits = 2), "%")
    return (cov, tot)
end

fcs = merged_coverage()
report(fcs)
LCOV.writefile("lcov.info", fcs)
clean_folder.(filter(isdir, CLEAN_DIRS))
