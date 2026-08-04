# Build the TensND.jl documentation.
#
#     julia --project=docs docs/make.jl
#
# Before paying for a full build, run the cheap guard that catches the one
# mistake able to break every page at once — a published script whose
# `Pkg.activate` lost its `#jl` marker (see `docs/literate.jl`):
#
#     grep -Ln '#jl' $(grep -l 'Pkg.activate' scripts/*.jl)   # prints nothing
#
# `docs/literate.jl` re-checks this at build time and errors out early.

using Documenter
using DocumenterCitations
# Renders the ```mermaid blocks (type hierarchy, operator dependency graphs).
# The diagram source is the page itself — no regenerated asset to keep in sync.
using DocumenterMermaid
using TensND
using SymPy

# Generates the tutorial pages (+ companion notebooks and clean scripts) from
# the curated `scripts/` demos. Must run before `makedocs`, so the generated
# markdown exists when the `pages` list below references it.
include("literate.jl")

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

makedocs(
    modules = [TensND],
    authors = "Jean-François Barthélémy and contributors",
    remotes = nothing,
    sitename = "TensND.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://MicroPoroChemoMechanics.github.io/TensND.jl",
        repolink = "https://github.com/MicroPoroChemoMechanics/TensND.jl",
        edit_link = "main",
        collapselevel = 1,
        mathengine = Documenter.MathJax3(
            Dict{Symbol, Any}(
                :tex => Dict{String, Any}(
                    "inlineMath" => [["\$", "\$"], ["\\(", "\\)"]],
                    "tags" => "ams",
                ),
                :options => Dict(
                    "ignoreHtmlClass" => "tex2jax_ignore",
                    "processHtmlClass" => "tex2jax_process",
                ),
            ),
            true,
        ),
        # Symbolic pages emit large SymPy expressions and the Walpole pages
        # embed several 6×6 matrices per section; the 200 KiB default is far
        # too tight for either.
        size_threshold = 1_500_000,
        size_threshold_warn = 800_000,
        example_size_threshold = 500_000,
        assets = ["assets/favicon.ico", "assets/custom.css"],
    ),
    plugins = [bib],
    pages = [
        "Home" => "index.md",
        # Ordered as a reading path with two pillars. First the *algebra*:
        # conventions, the products, then what a basis is (variance is the one
        # notion the Echoes manual deliberately does without), then the storage
        # convention everything else is expressed in, then the symmetry classes
        # from the most constrained to the least, and finally the projection
        # that maps an arbitrary tensor onto one of them. Then the *analysis*:
        # differential calculus on a curvilinear chart, and its restriction to
        # an embedded surface.
        "Theory" => [
            "theory/index.md",
            "theory/notation.md",
            "theory/tensor_algebra.md",
            "theory/bases_variance.md",
            "theory/kelvin_mandel.md",
            "theory/rotations.md",
            "theory/isotropic.md",
            "theory/walpole.md",
            "theory/walpole_extended.md",
            "theory/ti_parametrizations.md",
            "theory/orthotropy.md",
            "theory/projection.md",
            "theory/curvilinear.md",
            "theory/submanifolds.md",
        ],
        # Task-oriented. Every mathematical claim links back to Theory rather
        # than restating it.
        "Manual" => [
            "manual/installation.md",
            "manual/getting_started.md",
            "manual/bases.md",
            "manual/tensors.md",
            "manual/structured_tensors.md",
            "manual/parametrizations.md",
            "manual/projection.md",
            "manual/coorsystems.md",
            "manual/coorsystems_num.md",
            "manual/submanifolds.md",
            "manual/symbolic_and_numeric.md",
        ],
        # Every page under `generated/` is produced from `scripts/` by Literate
        # (see `docs/literate.jl`). That is an implementation detail the reader
        # has no reason to care about, so the grouping is thematic.
        "Tutorials" => [
            "tutorials/index.md",
            "Tensor algebra" => [
                "tutorials/generated/bases_variance.md",
                "tutorials/generated/tensor_algebra.md",
                "tutorials/generated/walpole.md",
                "tutorials/generated/walpole_extended.md",
                "tutorials/generated/projection.md",
            ],
            "Differential calculus" => [
                "tutorials/generated/symbolic_operators.md",
                "tutorials/generated/christoffel.md",
                "tutorials/generated/numerical_operators.md",
                "tutorials/generated/submanifolds.md",
                "tutorials/generated/operator_identities.md",
            ],
            "Applications in mechanics" => [
                "tutorials/generated/green_function.md",
                "tutorials/generated/cluster.md",
                "tutorials/generated/sphere_problems.md",
                "tutorials/generated/hill_tensors.md",
            ],
            "Types, differentiation and performance" => [
                "tutorials/generated/ad_interop.md",
                "tutorials/generated/performance.md",
            ],
        ],
        "Developer" => [
            "developer/architecture.md",
            "developer/adding_a_coorsystem.md",
            "developer/testing_conventions.md",
        ],
        "API" => [
            "api/index.md",
            "api/bases.md",
            "api/tensors.md",
            "api/structured.md",
            "api/walpole.md",
            "api/projection.md",
            "api/special.md",
            "api/coorsystems.md",
            "api/coorsystems_num.md",
            "api/submanifold.md",
            "api/symbolic.md",
            "api/full_index.md",
        ],
        "References" => "references.md",
    ],
    # Reports any *exported* symbol whose docstring is not included anywhere in
    # the documentation. This is what keeps the curated API pages honest: add an
    # export without placing it on a themed page and the build says so, instead
    # of the symbol quietly vanishing from the docs. Preferred over a catch-all
    # `@autodocs` page, which would document everything a second time and emit a
    # "duplicate docs" warning per symbol.
    checkdocs = :exports,
    warnonly = true,
)

deploydocs(;
    repo = "github.com/MicroPoroChemoMechanics/TensND.jl.git",
    devbranch = "main",
    push_preview = false,
)
