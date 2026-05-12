using TensND
using Documenter
using DocumenterCitations
using SymPy

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "references.bib");
    style = :numeric,
)

DocMeta.setdocmeta!(
    TensND,
    :DocTestSetup,
    :(using TensND, LinearAlgebra, SymPy, Tensors, OMEinsum, Rotations);
    recursive = true,
)

use_KaTeX = false

# Use MathJax syntax to define macros:
# - macro names should not be preceded by "\\"
macros = Dict(
    # "cellindices" => "\\mathcal{P}",
    # "conj" => "\\operatorname{conj}",
    # "D" => "\\mathrm d",
    # "dbldot" => "\\mathbin{\\mathord{:}}",
    # "dft" => "\\operatorname{DFT}",
    # "element" => "\\mathrm{e}",
    # "fftfreq" => "\\operatorname{Z}",
    # "I" => "\\mathrm i",
    # "integers" => "\\mathbb Z",
    # "naturals" => "\\mathbb N",
    # "PI" => "\\mathrm{\\pi}",
    # "reals" => "\\mathbb R",
    # "sinc" => "\\operatorname{sinc}",
    # "symotimes" => "\\stackrel{\\mathrm{s}}{\\otimes}",
    # "strain" => "\\boldsymbol{\\varepsilon}",
    # "strains" => "\\mathcal E",
    # "stress" => "\\boldsymbol{\\sigma}",
    # "stresses" => "\\mathcal S",
    # #"tens" => "\\@ifstar{\\boldsymbol}{\\mathbf}",
    # "tens" => "\\boldsymbol",
    # "tensors" => "\\mathcal T",
    # "tr" => "\\operatorname{tr}",
    # "tuple" => "\\mathsf",
    # "vec" => "\\boldsymbol",
)

if use_KaTeX
    macros2 = Dict()
    for (key, val) ∈ pairs(macros)
        macros2["\\"*key] = val
    end
    macros = macros2
end

mathengine =
    use_KaTeX ?
    KaTeX(
        Dict{Symbol,Any}(
            :delimiters => Dict{Any,Any}[
                Dict(:left => "\$", :right => "\$", display => false),
                Dict(:left => "\$\$", :right => "\$\$", display => true),
                Dict(:left => "\\[", :right => "\\]", display => true),
            ],
            :macros => macros,
        ),
    ) :
    MathJax3(
        Dict{Symbol,Any}(
            :tex => Dict{String,Any}(
                "macros" => macros,
                #"packages" => ["base", "ams", "autoload", "bm"],
                "inlineMath" => [["\$", "\$"], ["\\(", "\\)"]],
                "tags" => "ams",
            ),
            :options => Dict(
                "ignoreHtmlClass" => "tex2jax_ignore",
                "processHtmlClass" => "tex2jax_process",
            ),
        ),
        true,
    )

# format = Documenter.LaTeX(platform = "none")

makedocs(
    modules = [TensND],
    authors = "Jean-François Barthélémy and contributors",
    remotes = nothing,
    sitename = "TensND.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://MicroPoroChemoMechanics.codeberg.page/TensND.jl",
        repolink = "https://codeberg.org/MicroPoroChemoMechanics/TensND.jl",
        mathengine = mathengine,
        size_threshold_warn = 1_000_000,
        size_threshold = 1_000_000,
    ),
    plugins = [bib],
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "man/getting_started.md",
            "man/bases.md",
            "man/tensors.md",
            "man/coorsystems.md",
        ],
        "Tutorials" => [
            "tuto/nlayersphere.md",
            "tuto/coorsystems_num.md",
            "tuto/green_function.md",
            "tuto/projection.md",
        ],
        "API" => "api.md",
        "References" => "references.md",
    ],
)

# deploydocs(; repo = "codeberg.org/MicroPoroChemoMechanics/TensND.jl", devbranch = "main")
println("Build complet ! Consultez docs/build/")
