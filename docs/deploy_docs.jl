#!/usr/bin/env julia
# Build and deploy docs to Codeberg Pages (pages branch).
# Usage:  julia --project=docs/ docs/deploy_docs.jl
#
# Requires SSH alias "codeberg-docs" in ~/.ssh/config pointing to a deploy key
# with write access on the repository.
#
# Branch push  →  deploys to docs/dev/
# Tag push     →  deploys to docs/stable/ + docs/vX.Y.Z/

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
tag = try
    strip(readchomp(`git describe --exact-match --tags HEAD`))
catch
    ""
end
is_tag = !isempty(tag) && startswith(tag, "v") && tryparse(VersionNumber, tag[2:end]) !== nothing
ref = is_tag ? tag : strip(readchomp(`git rev-parse --abbrev-ref HEAD`))
println("$(is_tag ? "📦 TAG" : "🐛 BRANCH") $ref")

# ── Build docs (CI=false → Documenter skips its own deploydocs) ───────────────
ENV["CI"] = "false"
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
    rm(joinpath(pagesdir, "stable");      recursive = true, force = true)
    rm(joinpath(pagesdir, tag[2:end]);    recursive = true, force = true)
else
    rm(joinpath(pagesdir, "dev");         recursive = true, force = true)
end

# ── Copy build into target subdirectory ───────────────────────────────────────
if is_tag
    cp(builddir, joinpath(pagesdir, "stable"))
    cp(builddir, joinpath(pagesdir, tag[2:end]))
    println("$tag → stable/ + $(tag[2:end])/")
else
    cp(builddir, joinpath(pagesdir, "dev"))
    println("$ref → dev/")
end

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
