# Cross-platform documentation build and deployment to Codeberg Pages.
#
# Usage:  julia --project=docs/ deploy_docs.jl
#
# Behaviour:
#   HEAD on a branch  →  deploys to docs/dev/
#   HEAD on a tag     →  deploys to docs/vX.Y.Z/ and updates docs/stable/
#
# Authentication — set one as a persistent environment variable:
#   PROJECT_ACCESS_TOKEN  Codeberg Personal Access Token with "repository" scope.
#                         Codeberg → Settings → Applications → Generate Token.
#                         Recommended on Windows (avoids SSH path issues).
#   DOCUMENTER_KEY        Base64-encoded SSH private key matching a deploy key
#                         (write access) added in the Codeberg repository settings.
#                         Linux/macOS: export DOCUMENTER_KEY=$(base64 -w0 ~/.ssh/documenter_codeberg)
#                         Windows:     see WORKFLOW.md

# ── Tag auto-detection ─────────────────────────────────────────────────────────
tag = try
    strip(readchomp(`git describe --exact-match --tags HEAD`))
catch
    ""
end
on_tag = !isempty(tag) && startswith(tag, "v") && tryparse(VersionNumber, tag[2:end]) !== nothing

if on_tag
    @info "Tagged release: $tag  →  stable/ + $tag/"
else
    @info "Branch push  →  dev/"
end

# ── Authentication ─────────────────────────────────────────────────────────────
use_https = haskey(ENV, "PROJECT_ACCESS_TOKEN")

if !use_https && !haskey(ENV, "DOCUMENTER_KEY")
    error("""
No authentication credential found. Set one of:

  PROJECT_ACCESS_TOKEN  (recommended — HTTPS, no SSH path issues on Windows)
    Codeberg → Settings → Applications → Generate Token  (scope: repository)
    Linux/macOS:  export PROJECT_ACCESS_TOKEN=<token>
    Windows:      [Environment]::SetEnvironmentVariable("PROJECT_ACCESS_TOKEN", "<token>", "User")

  DOCUMENTER_KEY  (SSH — base64-encoded private key)
    The public key must be added as a deploy key (write access) in the
    Codeberg repository settings. See WORKFLOW.md for one-time setup.
""")
end

@info use_https ? "Auth: HTTPS (PROJECT_ACCESS_TOKEN)" : "Auth: SSH (DOCUMENTER_KEY)"

# ── Woodpecker CI environment variables ────────────────────────────────────────
# Documenter.jl auto-detects Woodpecker CI when CI="woodpecker" and deploys to
# the "pages" branch by default — which is what Codeberg Pages serves.
ENV["CI"]                = "woodpecker"
ENV["CI_SYSTEM_VERSION"] = "2.0.0"
ENV["CI_FORGE_URL"]      = "https://codeberg.org"
ENV["CI_REPO"]           = "MicroPoroChemoMechanics/TensND.jl"

if on_tag
    ENV["CI_PIPELINE_EVENT"] = "tag"
    ENV["CI_COMMIT_REF"]     = "refs/tags/$tag"
    ENV["CI_COMMIT_TAG"]     = tag
else
    ENV["CI_PIPELINE_EVENT"] = "push"
    ENV["CI_COMMIT_REF"]     = "refs/heads/main"
    ENV["CI_COMMIT_TAG"]     = ""
end

# HTTPS mode: Woodpecker uses SSH when DOCUMENTER_KEY is present, HTTPS otherwise.
use_https && delete!(ENV, "DOCUMENTER_KEY")

# ── Windows: fix temp path for MSYS2 SSH compatibility ────────────────────────
# Julia's bundled git (MINGW) invokes SSH via sh.exe, which treats backslashes
# as escape characters. A TMP/TEMP pointing to a forward-slash path ensures the
# SSH config file created by Documenter has a path that sh.exe passes intact to
# the "ssh -F <path>" invocation.
if Sys.iswindows()
    mkpath("C:\\tmp")
    ENV["TMP"]  = "C:/tmp"
    ENV["TEMP"] = "C:/tmp"
end

# ── Build and deploy ───────────────────────────────────────────────────────────
include(joinpath(@__DIR__, "docs", "make.jl"))
