#!/usr/bin/env julia
# Build and deploy docs to Codeberg Pages (pages branch).
# Usage:  julia --project=docs/ docs/deploy_docs.jl              # auto-detect branch or tag from git
#         julia --project=docs/ docs/deploy_docs.jl v0.1.8       # force stable deploy for a given tag
#
# Requires SSH alias "codeberg-docs" in ~/.ssh/config pointing to a deploy key
# with write access on the repository.
#
# Branch  →  deploys to dev/
# Tag     →  deploys to stable/ + vX.Y.Z/

using Pkg

docsdir = @__DIR__
Pkg.activate(docsdir)
println("✅ Activated docs/")

# Windows: mktempdir() uses TEMP which may have 8.3 short paths
# that MSYS2 git cannot handle correctly.
if Sys.iswindows()
    mkpath("C:\\tmp")
    ENV["TMP"]  = "C:/tmp"
    ENV["TEMP"] = "C:/tmp"
end

builddir   = joinpath(docsdir, "build")
repo_remote = "git@codeberg-docs:MicroPoroChemoMechanics/TensND.jl.git"

# ── Detect branch or tag ───────────────────────────────────────────────────────
tag = if !isempty(ARGS) && startswith(ARGS[1], "v") && tryparse(VersionNumber, ARGS[1][2:end]) !== nothing
    ARGS[1]
else
    try
        strip(readchomp(`git describe --exact-match --tags HEAD`))
    catch
        ""
    end
end
is_tag = !isempty(tag) && startswith(tag, "v") && tryparse(VersionNumber, tag[2:end]) !== nothing
ref = is_tag ? tag : strip(readchomp(`git rev-parse --abbrev-ref HEAD`))
println("$(is_tag ? "📦 TAG" : "🐛 BRANCH") $ref")

# ── Build docs (CI=false → Documenter skips its own deploydocs) ───────────────
ENV["CI"] = "true"
println("🔨 Building...")
include(joinpath(docsdir, "make.jl"))
println("✅ Build complete")

# ── Clone pages branch (or create orphan if it doesn't exist yet) ─────────────
pagesdir = mktempdir()
println("📥 Cloning pages branch into $pagesdir ...")

pages_exists = success(`git ls-remote --exit-code --heads $repo_remote pages`)
if pages_exists
    run(`git clone --branch pages --single-branch --depth 1 $repo_remote $pagesdir`)
else
    println("ℹ️  pages branch not found — creating orphan branch...")
    run(`git clone --single-branch --depth 1 $repo_remote $pagesdir`)
    cd(pagesdir) do
        run(`git checkout --orphan pages`)
        run(`git rm -rf .`)
    end
end

# ── Remove only the directories that will be replaced ─────────────────────────
# (keeps other versioned directories intact across deploys)
if is_tag
    rm(joinpath(pagesdir, "stable"); recursive = true, force = true)
    rm(joinpath(pagesdir, tag);      recursive = true, force = true)
else
    rm(joinpath(pagesdir, "dev");    recursive = true, force = true)
end

# ── Copy build into target subdirectory ───────────────────────────────────────
if is_tag
    cp(builddir, joinpath(pagesdir, "stable"))
    cp(builddir, joinpath(pagesdir, tag))
    println("$tag → stable/ + $tag/")
else
    cp(builddir, joinpath(pagesdir, "dev"))
    println("$ref → dev/")
end

# ── Generate versions.js and root index.html ──────────────────────────────────
version_dirs = sort(
    filter(d -> isdir(joinpath(pagesdir, d)) && startswith(d, "v") &&
                tryparse(VersionNumber, d[2:end]) !== nothing,
           readdir(pagesdir)),
    by = d -> VersionNumber(d[2:end]), rev = true,
)
versions = String[]
isdir(joinpath(pagesdir, "stable")) && push!(versions, "stable")
append!(versions, version_dirs)
isdir(joinpath(pagesdir, "dev")) && push!(versions, "dev")
newest = isempty(version_dirs) ? "dev" : version_dirs[1]
versions_list = join(["\"$v\"" for v in versions], ",")
write(joinpath(pagesdir, "versions.js"),
    "var DOC_VERSIONS = [$versions_list];\n" *
    "var DOCUMENTER_NEWEST = \"$newest\";\n" *
    "var DOCUMENTER_STABLE = \"stable\";\n")
println("✅ versions.js: $(join(versions, ", "))")
redirect = isdir(joinpath(pagesdir, "stable")) ? "stable" : "dev"
write(joinpath(pagesdir, "index.html"),
    "<!DOCTYPE html>\n<html>\n" *
    "<head><meta http-equiv=\"refresh\" content=\"0; url=$redirect/\"/></head>\n" *
    "<body><p>Redirecting to <a href=\"$redirect/\">documentation</a>.</p></body>\n" *
    "</html>\n")
println("✅ index.html → $redirect/")

# ── Commit and push ────────────────────────────────────────────────────────────
println("🚀 Pushing pages branch...")
cd(pagesdir) do
    run(`git config user.name "TensND Documenter"`)
    run(`git config user.email "docs@tensnd.codeberg.page"`)
    run(`git add -A`)
    if !isempty(strip(readchomp(`git status --porcelain`)))
        msg = is_tag ? "Deploy $tag" : "Deploy dev ($ref)"
        run(`git commit -m $msg`)
        run(`git push -u origin pages`)
        dest = is_tag ? tag[2:end] : "dev"
        println("✅ Deployed! → $(dest)/")
    else
        println("ℹ️  No changes to deploy.")
    end
end

println("🎉 Done!")
