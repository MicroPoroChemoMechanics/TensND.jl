#!/usr/bin/env julia
# Build and deploy docs to Codeberg Pages (pages branch).
#
# Requires SSH alias "codeberg-docs" in ~/.ssh/config pointing to a deploy key
# with write access on the repository.
#
# ── How to choose the deployment target ───────────────────────────────────────
# Priority: command-line argument > DEPLOY_TARGET below > git auto-detect
#
# Command-line:
#   julia --project=docs/ docs/deploy_docs.jl           # auto-detect
#   julia --project=docs/ docs/deploy_docs.jl v0.1.8   # force tag
#
# In-file (set DEPLOY_TARGET below, takes effect when no argument is passed):
#   nothing     — auto-detect from git: branch → dev/,  tag → stable/ + vX.Y.Z/
#   "dev"       — force deploy to dev/
#   "project"   — read version from Project.toml → stable/ + vX.Y.Z/

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

builddir    = joinpath(docsdir, "build")
repo_remote = "git@codeberg-docs:MicroPoroChemoMechanics/TensND.jl.git"

# ── Deployment target ─────────────────────────────────────────────────────────
# DEPLOY_TARGET = nothing       # auto-detect from git
# DEPLOY_TARGET = "dev"       # force dev/ deploy
DEPLOY_TARGET = "project"   # stable/ + vX.Y.Z/ using version from Project.toml

function read_project_version()
    for line in eachline(joinpath(dirname(docsdir), "Project.toml"))
        m = match(r"^version\s*=\s*\"(.+)\"", line)
        m !== nothing && return "v" * m.captures[1]
    end
    error("Could not read version from Project.toml")
end

# ── Detect branch or tag ───────────────────────────────────────────────────────
tag = if !isempty(ARGS) && startswith(ARGS[1], "v") && tryparse(VersionNumber, ARGS[1][2:end]) !== nothing
    ARGS[1]                          # command-line argument wins
elseif DEPLOY_TARGET == "project"
    read_project_version()
elseif DEPLOY_TARGET == "dev"
    ""                               # empty string → dev mode below
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

# ── Copy build into target subdirectory and write siteinfo.js ─────────────────
function write_siteinfo(dir, version)
    write(joinpath(dir, "siteinfo.js"), "var DOCUMENTER_CURRENT_VERSION = \"$version\";\n")
end

if is_tag
    cp(builddir, joinpath(pagesdir, "stable"))
    write_siteinfo(joinpath(pagesdir, "stable"), "stable")
    cp(builddir, joinpath(pagesdir, tag))
    write_siteinfo(joinpath(pagesdir, tag), tag)
    println("$tag → stable/ + $tag/")
else
    cp(builddir, joinpath(pagesdir, "dev"))
    write_siteinfo(joinpath(pagesdir, "dev"), "dev")
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
    run(`git config core.autocrlf false`)
    run(`git add -A`)
    if !isempty(strip(readchomp(`git status --porcelain`)))
        msg = is_tag ? "Deploy $tag" : "Deploy dev ($ref)"
        run(`git commit -m $msg`)
        run(`git push -u origin pages`)
        dest = is_tag ? tag : "dev"
        println("✅ Deployed! → $(dest)/")
    else
        println("ℹ️  No changes to deploy.")
    end
end

println("🎉 Done!")
