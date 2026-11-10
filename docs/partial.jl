# ── Partial builds: `TENSND_DOCS_ONLY` ───────────────────────────────────────
#
# The full site is long to build, nearly all of the time being spent inside
# `@example` blocks, so a change touching one page pays for every other. This
# builds only the pages whose path contains one of a comma-separated list of
# patterns:
#
#     TENSND_DOCS_ONLY=theory/walpole julia --project=docs docs/make.jl
#     TENSND_DOCS_ONLY=manual/,theory/     julia --project=docs docs/make.jl
#
# FILTERING `pages` IS NOT ENOUGH, and this is the whole reason the file exists.
# Documenter's `SetupBuildDirectory` walks the **source directory** and adds
# every `.md` it finds; `pages` only decides the navigation. Pruning the page
# tree alone leaves every `@example` on the site still running. What does work is
# giving Documenter a different source tree, in which a kept page is linked
# rather than copied (Documenter walks with `follow_symlinks = true`), and a
# page tree filtered to match.
#
# WHAT A PARTIAL BUILD PROVES: that the `@example` blocks of the pages it kept
# run, and that those pages render.
#
# WHAT IT CANNOT PROVE, and says so: anything about links or about docstring
# coverage. Every `@ref` into a pruned page points at nothing, and every
# docstring whose API page was pruned is "missing", so `cross_references` and
# `missing_docs` are demoted to warnings for that build only, and deployment is
# refused. The full build, which tolerates no warning, is the only one that
# deploys.
#
# `index.md` and `references.md` are always kept: the site root the navigation
# is rendered around, and the page every citation on a kept page resolves into.
#
# THE THIRD PIECE IS IN `docs/src/.vitepress/config.mts`. Demoting
# `cross_references` makes Documenter emit the unresolved references literally,
# and VitePress refuses to build over dead links. Its `ignoreDeadLinks` is
# therefore driven by this same environment variable, so that the Julia half and
# the Node half of one decision cannot disagree.
#
# Ported from MeanFieldHomogenization's `docs/partial.jl`.

const ALWAYS_KEPT = ("index.md", "references.md")

"""
    page_leaves(node) -> Vector{String}

Every page path under `node`, in build order. A leaf is a bare path or the value
of a `"Title" => "path"` pair; a `"Title" => [...]` pair is a section.
"""
page_leaves(node::AbstractString) = [node]
page_leaves(node::Pair) = page_leaves(node.second)
page_leaves(node::AbstractVector) = reduce(vcat, page_leaves.(node); init = String[])

"""
    is_kept(page, patterns) -> Bool

Whether a source-relative page path survives the filter.
"""
is_kept(page, pats) = page in ALWAYS_KEPT || any(p -> occursin(p, page), pats)

"""
    keep_pages(node, patterns) -> node or `nothing`

The page tree with only the leaves that `is_kept` admits. Returns `nothing`
when nothing under `node` survives, which prunes a section header whose whole
contents were filtered out.
"""
keep_pages(node::AbstractString, pats) = is_kept(node, pats) ? node : nothing
function keep_pages(node::Pair, pats)
    kept = keep_pages(node.second, pats)
    return kept === nothing ? nothing : (node.first => kept)
end
function keep_pages(node::AbstractVector, pats)
    kept = Any[]
    for child in node
        k = keep_pages(child, pats)
        k === nothing || push!(kept, k)
    end
    return isempty(kept) ? nothing : kept
end

"""
    pruned_source(srcdir, patterns) -> String

A temporary directory mirroring `srcdir`, with the Markdown pages the filter
rejects left out. Everything that is not a page -- images, `references.bib`,
`.vitepress` -- is carried through, because a kept page may reference any of it.

Pages are linked and everything else is COPIED. Julia's `cp` reproduces a
symlink as a symlink, so a linked `.vitepress/config.mts` would let
DocumenterVitepress, which fills that file's `REPLACE_ME_DOCUMENTER_VITEPRESS`
markers in what it believes is the build copy, write straight through the link
and overwrite the template in the repository. That happened once in
ChemistryLab.
"""
function pruned_source(srcdir, pats)
    tmp = mktempdir(; prefix = "tensnd_docs_")
    for (root, _, files) in walkdir(srcdir)
        rel = relpath(root, srcdir)
        mkpath(normpath(joinpath(tmp, rel)))
        for f in files
            src = abspath(joinpath(root, f))
            dst = normpath(joinpath(tmp, rel, f))
            if endswith(f, ".md")
                is_kept(replace(normpath(joinpath(rel, f)), '\\' => '/'), pats) || continue
                # A symlink where the platform allows one, a copy where it does
                # not (Windows without Developer Mode); a page is only ever read.
                try
                    symlink(src, dst)
                catch
                    cp(src, dst; force = true)
                end
            else
                cp(src, dst; force = true)
            end
        end
    end
    return tmp
end

const DOCS_ONLY = strip(get(ENV, "TENSND_DOCS_ONLY", ""))
const PARTIAL_BUILD = !isempty(DOCS_ONLY)

const DOCS_SOURCE = if PARTIAL_BUILD
    patterns = filter(!isempty, strip.(split(DOCS_ONLY, ',')))

    # Tested against the leaves BEFORE filtering. The two pages above are kept
    # unconditionally, so the filtered tree is never empty, and a typo in the
    # variable would otherwise build a two-page site and report success.
    matched = filter(pg -> any(p -> occursin(p, pg), patterns), page_leaves(pages))
    isempty(matched) && error(
        "TENSND_DOCS_ONLY=\"$DOCS_ONLY\" matches no page. A pattern is matched " *
            "against the path as it appears in docs/pages.jl, for example " *
            "`theory/walpole` or `manual/`.",
    )

    global pages = keep_pages(pages, patterns)
    @warn """
    PARTIAL DOCUMENTATION BUILD -- not deployable, and it checks neither links
    nor docstring coverage. `cross_references` and `missing_docs` are demoted to
    warnings, because a reference into a pruned page, and a docstring whose API
    page was pruned, would otherwise be reported as broken. Run the full build
    before merging.
    Patterns: $(join(patterns, ", "))
    Pages   : $(join(matched, ", "))"""

    pruned_source(joinpath(@__DIR__, "src"), patterns)
else
    "src"
end
