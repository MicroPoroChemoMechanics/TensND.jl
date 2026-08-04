# docs/literate.jl
#
# Runs Literate.jl over the `scripts/` demos before `makedocs`. Each published
# script produces three artifacts from a single source:
#
#   - a Documenter-ready markdown page  -> docs/src/tutorials/generated/
#     (executed by Documenter itself through `@example`; Literate's
#     `markdown(...; documenter = true)` leaves `execute = false`)
#   - a pre-run Jupyter notebook        -> docs/generated_notebooks/
#   - a cleaned standalone .jl script   -> docs/generated_scripts/
#
# `docs/make.jl` includes this file *before* `makedocs`, so the generated
# markdown exists by the time the `pages` list references it.
#
# The scripts stay runnable on their own (`julia scripts/NN_name.jl`); the
# Literate markers are the only thing that distinguishes them from ordinary
# demos. `PUBLISHED_SCRIPTS` maps each script to its **page name**: scripts
# carry a numeric prefix encoding a running order, pages carry thematic names,
# so inserting a tutorial never forces a renumbering.
#
# ┌─────────────────────────────────────────────────────────────────────────┐
# │ The `#jl` contract, which is not optional                               │
# │                                                                          │
# │ Every `import Pkg` / `Pkg.activate(...)` line in a published script MUST │
# │ end with `#jl`. That marker keeps the line in the standalone script and  │
# │ in the generated clean script, and strips it from the markdown and the   │
# │ notebook — which is what we want, because those already run inside the   │
# │ `docs` environment.                                                      │
# │                                                                          │
# │ Forgetting it does **not** break only that page. `Pkg.activate` mutates  │
# │ global process state, so Documenter switches project mid-build and every │
# │ later `@example` block, on every page, fails with                        │
# │ `Package TensND not found in current path`. Run the guard in the module  │
# │ docstring of `docs/make.jl` before paying for a full build.              │
# └─────────────────────────────────────────────────────────────────────────┘

using Literate

const SCRIPTS_DIR = joinpath(@__DIR__, "..", "scripts")
const TUTORIAL_MD_DIR = joinpath(@__DIR__, "src", "tutorials", "generated")
const NOTEBOOK_DIR = joinpath(@__DIR__, "generated_notebooks")
const CLEAN_SCRIPT_DIR = joinpath(@__DIR__, "generated_scripts")

# script file => generated page name (without extension)
const PUBLISHED_SCRIPTS = [
    # ── Tensor algebra ───────────────────────────────────────────────────────
    "01_bases_variance.jl" => "bases_variance",
    "02_tensor_algebra.jl" => "tensor_algebra",
    "03_walpole.jl" => "walpole",
    "04_walpole_extended.jl" => "walpole_extended",
    "05_projection.jl" => "projection",
    # ── Differential calculus ────────────────────────────────────────────────
    "10_symbolic_operators.jl" => "symbolic_operators",
    "11_christoffel.jl" => "christoffel",
    "12_numerical_operators.jl" => "numerical_operators",
    "13_submanifolds.jl" => "submanifolds",
    # ── Applications in mechanics ────────────────────────────────────────────
    "20_green_function.jl" => "green_function",
    "21_cluster.jl" => "cluster",
    "22_sphere_problems.jl" => "sphere_problems",
    "24_hill_tensors.jl" => "hill_tensors",
    # ── Types, differentiation and performance ───────────────────────────────
    "30_ad_interop.jl" => "ad_interop",
    "31_performance.jl" => "performance",
]

# Retired sources, absorbed by the scripts above:
#
#   Cylinder.jl, numerical_operators.jl  -> 10_symbolic_operators, 12_numerical_operators
#   Green.jl                             -> 20_green_function
#   cluster.jl                           -> 21_cluster
#   sphere.jl + docs/src/tuto/nlayersphere.md -> 22_sphere_problems
#   sif.jl, sifTI.jl                     -> 24_hill_tensors (sifTI.jl re-implemented
#                                           the Walpole basis locally; it now uses
#                                           the library)
#   isotropization.jl                    -> 02_tensor_algebra
#   misctest.jl (Revise sandbox)         -> 11_christoffel + 13_submanifolds

"""
    check_pkg_markers()

Fail loudly if a published script carries a bare `Pkg.activate` — the one
mistake that breaks the whole build rather than a single page. Cheap enough to
run on every build.
"""
function check_pkg_markers()
    offenders = String[]
    for (script, _) in PUBLISHED_SCRIPTS
        path = joinpath(SCRIPTS_DIR, script)
        isfile(path) || continue
        for line in eachline(path)
            occursin(r"^\s*(import\s+Pkg|Pkg\.activate)", line) &&
                !occursin("#jl", line) &&
                push!(offenders, "$script: $(strip(line))")
        end
    end
    isempty(offenders) || error(
        "Published scripts carry a bare `Pkg.activate` (missing the `#jl` marker).\n" *
            "This would switch the active project mid-build and break every @example\n" *
            "block in the documentation, not just these pages:\n  " *
            join(offenders, "\n  ")
    )
    return nothing
end

function build_tutorial_pages()
    check_pkg_markers()
    mkpath(TUTORIAL_MD_DIR)
    mkpath(NOTEBOOK_DIR)
    mkpath(CLEAN_SCRIPT_DIR)
    for (script, page) in PUBLISHED_SCRIPTS
        src = joinpath(SCRIPTS_DIR, script)
        isfile(src) || error("docs/literate.jl: published script not found: $src")
        Literate.markdown(src, TUTORIAL_MD_DIR; documenter = true, name = page)
        Literate.notebook(src, NOTEBOOK_DIR; name = page)
        Literate.script(src, CLEAN_SCRIPT_DIR; name = page)
    end
    return nothing
end

build_tutorial_pages()
