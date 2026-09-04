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
# VitePress renders the site from the Markdown Documenter emits, and typesets
# every formula at build time into static SVG: no MathJax bundle ever reaches
# the reader's browser.
using DocumenterVitepress
# Renders the ```mermaid blocks (type hierarchy, operator dependency graphs).
# The diagram source is the page itself — no regenerated asset to keep in sync.
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

# ── Stopgap: CitationSiteNode in the Markdown writer ─────────────────────────
#
# GUARDED, and the guard is what makes one file serve both versions.
#
# `CitationSiteNode` exists only in DocumenterCitations 1.5, while
# DocumenterVitepress declares a weak dependency `DocumenterCitations = "1 - 1.4"`
# — an upstream cap, so the docs environment resolves to 1.4 and the name is not
# there. Referring to it unconditionally is an `UndefVarError` raised before the
# first page is built, which is exactly how this was found. On 1.4 the citations
# need no help at all: the node is what 1.5 introduced.
#
# When DocumenterVitepress raises its cap, widening the bound in
# `docs/Project.toml` is the only edit needed — this method starts applying again
# on its own.
if isdefined(DocumenterCitations, :CitationSiteNode)
# DocumenterCitations 1.5 wraps every expanded citation in a `CitationSiteNode`,
# whose only purpose is to give the citation an HTML anchor so the bibliography
# can link back to it. Its own docstring calls it "transparent in any output
# format other than HTML", and both the LaTeX writer and MDFlatten implement it
# as "render my children".
#
# DocumenterVitepress 0.3.5 ships a DocumenterCitations extension, but it covers
# only `BibliographyNode`. With no method for `CitationSiteNode`, the writer
# falls through to its generic branch, which prints `Markdown.plain(element)` —
# so every one of this manual's citations came out as the literal text
# `DocumenterCitations.CitationSiteNode("kachanov1992-cite-1")`.
#
# The same one-line treatment as the other non-HTML writers. Remove this once
# DocumenterVitepress covers the node upstream.
    function DocumenterVitepress.render(
            io::IO,
            mime::MIME"text/plain",
            node::Documenter.MarkdownAST.Node,
            ::DocumenterCitations.CitationSiteNode,
            page,
            doc;
            kwargs...,
        )
        return DocumenterVitepress.render(
            io, mime, node, node.children, page, doc; kwargs...,
        )
    end
end

# ── Stopgap: heading anchors that contain LaTeX ──────────────────────────────
# DocumenterVitepress builds each heading as `## <text> {#<slug>}`, where the
# slug is Documenter's anchor label passed through its own
# `sanitized_anchor_label` — whose comment says "vitepress doesn't like special
# markdown characters in the id slug", but which only strips `[ ] ( ) *`.
#
# A heading such as `## The COD tensor ``\boldsymbol{B}``` yields the slug
# `The-COD-tensor-\boldsymbol{B}`. VitePress's `{#...}` parser rejects the
# backslash and the braces, so it treats the whole suffix as *text*: the heading
# renders as "The COD tensor {#The-COD-tensor-\boldsymbol{B}}", the formula is
# dropped, and the same garbage lands in the "On this page" outline. Twenty-three
# headings across seven pages were affected.
#
# Stripping those characters from the slug is safe here: nothing links to those
# anchors (checked against every `](…#…)` destination in the generated Markdown),
# and this narrows to headings only, leaving docstring anchors — which legitimately
# carry braces, are emitted as raw `<a id=…>`, and *are* linked to — untouched.
#
# Remove once `sanitized_anchor_label` covers these characters upstream.
function DocumenterVitepress.render(
        io::IO,
        mime::MIME"text/plain",
        node::Documenter.MarkdownAST.Node,
        header::Documenter.AnchoredHeader,
        page,
        doc;
        kwargs...,
    )
    anchor = header.anchor
    label = DocumenterVitepress.sanitized_anchor_label(anchor)
    id = replace(replace(label, r"[\\{}]" => ""), " " => "-")
    heading = first(node.children)
    println(io)
    print(io, "#"^(heading.element.level), " ")
    heading_iob = IOBuffer()
    DocumenterVitepress.render(heading_iob, mime, node, heading.children, page, doc; kwargs...)
    print(io, rstrip(String(take!(heading_iob))))
    print(io, " {#$(id)}")
    if haskey(kwargs, :inventory)
        item = DocumenterVitepress.InventoryItem(
            name = anchor.id,
            domain = "std",
            role = "label",
            dispname = DocumenterVitepress._get_inventory_dispname(
                anchor.id, Documenter.MDFlatten.mdflatten(anchor.node)
            ),
            priority = -1,
            uri = DocumenterVitepress._get_inventory_uri(doc, page, id),
        )
        push!(kwargs[:inventory], item)
    end
    println(io)
    return nothing
end

# ── Stopgap: ordered lists start at 2, and swallow their first item ──────────
# DocumenterVitepress numbers ordered-list items with `bullet(i) = "$(i+1). "`,
# but `enumerate` is already 1-based: every ordered list in the manual came out
# numbered from 2. It also emits no blank line before the list.
#
# Together those two produce the damage seen on the References page. A list whose
# first marker is `2.` cannot interrupt a paragraph — CommonMark allows that only
# for a list starting at `1.` — so, with no blank line to separate them, the first
# entry was absorbed into the preceding prose as plain text and the list began at
# `3.`. Eighteen ordered lists across the manual were affected; not one of them
# started at 1.
#
# Remove once the numbering is fixed upstream.
function DocumenterVitepress.render(
        io::IO,
        mime::MIME"text/plain",
        node::Documenter.MarkdownAST.Node,
        list::Documenter.MarkdownAST.List,
        page,
        doc;
        kwargs...,
    )
    bullet(i) = list.type === :ordered ? "$(i). " : "- "
    println(io)
    iob = IOBuffer()
    for (i, item) in enumerate(node.children)
        DocumenterVitepress.render(
            iob, mime, item, item.children, page, doc; prenewline = false, kwargs...
        )
        eachline = split(String(take!(iob)), '\n')
        # Continuation lines must line up with the text, i.e. under the marker's
        # full width. Upstream hard-codes two spaces, which fits `- ` but not
        # `1. `: a display equation inside an ordered item fell out of the list,
        # splitting it in two and restarting the numbering.
        pad = " "^length(bullet(i))
        eachline[2:end] .= pad .* eachline[2:end]
        final_string = join(eachline, '\n')
        endswith(final_string, '\n') || (final_string *= "\n")
        print(io, bullet(i))
        print(io, final_string)
    end
    return nothing
end

makedocs(
    modules = [TensND],
    authors = "Jean-François Barthélémy and contributors",
    remotes = nothing,
    sitename = "TensND.jl",
    # The favicon and the logo are picked up from `docs/src/assets`, the sidebar
    # is derived from `pages`, and there is no page-size ceiling to raise: the
    # large SymPy outputs travel through Vite rather than Documenter's guard.
    format = DocumenterVitepress.MarkdownVitepress(;
        repo = "https://github.com/MicroPoroChemoMechanics/TensND.jl",
        devbranch = "main",
        devurl = "dev",
        deploy_url = "https://MicroPoroChemoMechanics.github.io/TensND.jl",
        description = "Tensor algebra and differential geometry on any coordinate system, in Julia",
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
    # No blanket `warnonly = true`: it cancelled the very check the comment
    # above asks for, so an unplaced export warned and passed. Measured before
    # removing it — every exported, documented name is on a page.
)

# DocumenterVitepress writes a real directory per version rather than the
# symlinks Documenter used, so it needs its own `deploydocs`.
DocumenterVitepress.deploydocs(;
    repo = "github.com/MicroPoroChemoMechanics/TensND.jl.git",
    target = joinpath(@__DIR__, "build"),
    branch = "gh-pages",
    devbranch = "main",
    push_preview = false,
)
