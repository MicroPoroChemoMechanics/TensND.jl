# The source of the files the documentation runs, shown where they run.
#
# A page that `include`s a script of the repository shows the script itself,
# folded under the block that runs it, so that a reader sees the code without
# leaving the page, and links to it on GitHub at the commit the site was built
# from, so that the link and the code shown are the same file. Included by
# `make.jl`.

using Markdown

const REPO_URL = "https://github.com/MicroPoroChemoMechanics/TensND.jl"
const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

# The commit the site is built from: the one the CI checked out, else the head of
# the working tree, else the default branch when there is no git at all.
function source_ref(root = REPO_ROOT)
    sha = get(ENV, "GITHUB_SHA", "")
    isempty(sha) || return sha
    out = IOBuffer()
    found = success(pipeline(`git -C $root rev-parse HEAD`; stdout = out, stderr = devnull))
    return found ? strip(String(take!(out))) : "main"
end
const SOURCE_REF = source_ref()

"The GitHub address of `relpath` (relative to the repository root) at `SOURCE_REF`."
source_url(relpath) = "$(REPO_URL)/blob/$(SOURCE_REF)/$(relpath)"

"""
    show_source(relpath) -> Markdown.MD

The file `relpath` of the repository as a highlighted Julia block, preceded by its
link on GitHub. Called from an `@eval` block placed between the two `@raw html`
blocks that open and close a `<details>` element, so that it stays folded until
the reader asks for it.
"""
function show_source(relpath::AbstractString)
    file = joinpath(REPO_ROOT, relpath)
    isfile(file) || throw(ArgumentError("show_source: no file $(file); the path is relative to the repository root"))
    link = Markdown.Paragraph(
        Any[
            "The file ", Markdown.Link(Any[Markdown.Code(relpath)], source_url(relpath)),
            " on GitHub, at the commit this page was built from:",
        ]
    )
    return Markdown.MD(Any[link, Markdown.Code("julia", rstrip(read(file, String)))])
end
