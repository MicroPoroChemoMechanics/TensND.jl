# =============================================================================
#  make.jl — the drawn figures of the documentation, written into
#  docs/src/assets/ (committed; the documentation build never runs this).
#
#      julia --project=docs/assets_src docs/assets_src/make.jl                 # all
#      julia --project=docs/assets_src docs/assets_src/make.jl coordinates     # one family
#
#  Families: coordinates, symmetries.
# =============================================================================

for f in ("style.jl", "coordinates.jl", "symmetries.jl")
    include(joinpath(@__DIR__, "src", f))
end

const FAMILIES = Dict{String, Function}(
    "coordinates" => coordinates_figures,
    "symmetries" => symmetries_figures,
)

function main(args)
    names = isempty(args) ? sort(collect(keys(FAMILIES))) : args
    for n in names
        haskey(FAMILIES, n) || error("unknown family $n; expected $(sort(collect(keys(FAMILIES))))")
        println("── ", n)
        FAMILIES[n]()
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
