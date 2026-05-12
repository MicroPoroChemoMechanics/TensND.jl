#!/usr/bin/env julia
# VSCode: Ctrl+F5 → creates pages branch if missing!

using Pkg

# Auto-activate
docsdir = normpath(joinpath(@__DIR__, "..", "docs"))
Pkg.activate(docsdir)
println("✅ Activated docs/ ($(basename(docsdir)))")

using LibGit2
using FilePathsBase
using Dates

builddir = joinpath(docsdir, "build")
repo = "MicroPoroChemoMechanics/TensND.jl"

println("📁 Build: $builddir")
println("🔑 SSH: codeberg-docs")

# ── Build ─────────────────────────────────────────────────────────────────────
ENV["CI"] = "false"
println("🔨 Building...")
include(joinpath(docsdir, "make.jl"))
println("✅ Build complete")

# ── Detect ref ────────────────────────────────────────────────────────────────
tag = success(`git describe --exact-match --tags HEAD`) ?
    strip(readchomp(`git describe --exact-match --tags HEAD`)) : ""
is_tag = !isempty(tag) && startswith(tag, "v")
ref = is_tag ? tag : strip(readchomp(`git rev-parse --abbrev-ref HEAD`))
println("$(is_tag ? "📦 TAG" : "🐛 BRANCH") $ref")

# ── Clone/Create pages branch ─────────────────────────────────────────────────
pagesdir = mktempdir()
repo_url_docs = "git@codeberg-docs:$repo.git"

println("📥 Cloning/creating pages branch...")
try
    run(`git clone --branch pages --single-branch --depth 1 $repo_url_docs $pagesdir`)
catch e
    println("ℹ️  pages branch not found, creating from main...")
    run(`git clone --single-branch --depth 1 $repo_url_docs $pagesdir`)
    cd(pagesdir) do
        run(`git checkout --orphan pages`)
        run(`git rm -rf .`)
    end
end

# Clean
for f in readdir(pagesdir; join=true)
    basename(f) ∉ [".git", ".domains"] && rm(f; recursive=true, force=true)
end

# ── Deploy structure ──────────────────────────────────────────────────────────
cp(builddir, pagesdir; force=true)

if is_tag
    cp(builddir, joinpath(pagesdir, "stable"); force=true)
    cp(builddir, joinpath(pagesdir, tag[2:end]); force=true)
    println("$tag → stable/ + $(tag[2:end])/")
else
    cp(builddir, joinpath(pagesdir, "dev"); force=true)
    println("$ref → dev/")
end

# ── Push ──────────────────────────────────────────────────────────────────────
println("🚀 Pushing pages branch...")

run(`git -C $pagesdir config user.name "TensND Documenter"`)
run(`git -C $pagesdir config user.email "docs@tensnd.codeberg.page"`)
run(`git -C $pagesdir add -A`)

if !isempty(readchomp(`git -C $pagesdir status --porcelain`))
    msg = "Initial deploy $ref ($(Dates.format(now(), "yyyy-MM-dd HH:MM")))"
    run(`git -C $pagesdir commit -m $msg`)
    run(`git -C $pagesdir push -u origin pages`)  # -u sets upstream
    println("✅ Pages branch created! https://tensnd.codeberg.page")
else
    println("ℹ️ No changes")
end

println("🎉 First deploy complete!")
