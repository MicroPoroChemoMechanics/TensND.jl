# The page tree. `make.jl` includes this file; a partial build
# (`TENSND_DOCS_ONLY`, see `partial.jl`) prunes it.
#
# Each chapter answers one question.
#
#   **Theory**    -- what are the objects and why do the formulas hold?
#                   Definitions and derivations, readable without running
#                   anything.
#   **Manual**    -- how is the package called? One page per kind of object,
#                   syntax first, on the smallest example that shows it.
#   **Tutorials** -- how is a calculation carried through? Worked cases
#                   generated from `scripts/` by Literate.jl and executed by the
#                   build.
#   **Developer** -- how is the source organized, and how is it extended?
#   **API**       -- the docstrings, by theme.
#
# Theory, Manual, Tutorials and API share the chapter names "Tensors and bases",
# "Material symmetries", "Coordinate systems" and "Surfaces", in this order.
#
# The Theory pages follow a reading path with two pillars: the algebra
# (conventions, products, bases and variance, the storage convention, the
# symmetry classes from the most to the least constrained, the projection onto
# one of them), then the analysis (differential calculus on a curvilinear chart
# and its restriction to an embedded surface).
#
# Every page under `tutorials/generated/` is produced from `scripts/` by
# Literate (see `literate.jl`), which the reader has no reason to care about, so
# the tutorials are grouped by subject.
#
# `pages` is a plain global rather than a `const`: `partial.jl` replaces it with
# the pruned tree of a partial build.

pages = [
    "Home" => "index.md",
    "Theory" => [
        "theory/index.md",
        "Tensors and bases" => [
            "theory/notation.md",
            "theory/tensor_algebra.md",
            "theory/bases_variance.md",
            "theory/kelvin_mandel.md",
            "theory/rotations.md",
        ],
        "Material symmetries" => [
            "theory/isotropic.md",
            "theory/walpole.md",
            "theory/walpole_extended.md",
            "theory/ti_parametrizations.md",
            "theory/orthotropy.md",
            "theory/cubic.md",
            "theory/projection.md",
        ],
        "Coordinate systems" => [
            "theory/curvilinear.md",
            "theory/classical_systems.md",
            "theory/spheroidal.md",
            "theory/ellipsoidal.md",
            "theory/coordinate_changes.md",
            "theory/spheroidal_harmonics.md",
        ],
        "Surfaces" => [
            "theory/submanifolds.md",
        ],
    ],
    "Manual" => [
        "manual/index.md",
        "First steps" => [
            "manual/installation.md",
            "manual/getting_started.md",
        ],
        "Tensors and bases" => [
            "manual/bases.md",
            "manual/tensors.md",
            "manual/random_tensors.md",
        ],
        "Material symmetries" => [
            "manual/structured_tensors.md",
            "manual/parametrizations.md",
            "manual/elastic_moduli.md",
            "manual/projection.md",
        ],
        "Coordinate systems" => [
            "manual/coorsystems.md",
            "manual/coorsystems_num.md",
            "manual/charts.md",
            "manual/spheroidal_harmonics.md",
        ],
        "Surfaces" => [
            "manual/submanifolds.md",
        ],
        "Symbolic and numeric" => [
            "manual/symbolic_and_numeric.md",
        ],
    ],
    "Tutorials" => [
        "tutorials/index.md",
        "Tensors and bases" => [
            "tutorials/generated/bases_variance.md",
            "tutorials/generated/tensor_algebra.md",
        ],
        "Material symmetries" => [
            "tutorials/generated/walpole.md",
            "tutorials/generated/walpole_extended.md",
            "tutorials/generated/projection.md",
        ],
        "Coordinate systems" => [
            "tutorials/generated/symbolic_operators.md",
            "tutorials/generated/christoffel.md",
            "tutorials/generated/numerical_operators.md",
            "tutorials/generated/operator_identities.md",
            "tutorials/generated/operators_catalog.md",
            "tutorials/generated/spheroidal_appendix.md",
            "tutorials/generated/ellipsoidal.md",
            "tutorials/generated/spheroidal_reexpansion.md",
        ],
        "Surfaces" => [
            "tutorials/generated/submanifolds.md",
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
        "Tensors and bases" => [
            "api/bases.md",
            "api/tensors.md",
            "api/special.md",
        ],
        "Material symmetries" => [
            "api/structured.md",
            "api/walpole.md",
            "api/cubic.md",
            "api/projection.md",
        ],
        "Coordinate systems" => [
            "api/coorsystems.md",
            "api/coorsystems_num.md",
            "api/charts.md",
            "api/spheroidal_harmonics.md",
        ],
        "Surfaces" => [
            "api/submanifold.md",
        ],
        "Symbolic" => [
            "api/symbolic.md",
        ],
        "api/full_index.md",
    ],
    "Nomenclature" => "nomenclature.md",
    "References" => "references.md",
]
