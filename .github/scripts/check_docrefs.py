#!/usr/bin/env python3
"""Catch the one class of Documenter cross-reference failure that only CI sees.

    python3 .github/scripts/check_docrefs.py [src_dir]

## The rule this models

Documenter resolves an `@ref` inside a docstring **in the module where that
docstring is defined**. So a docstring in `MeanFieldHomogenization.Superspheres`
that writes

    [`RVE`](@ref)

makes Documenter look for the binding `Superspheres.RVE`, which does not exist:
`RVE` lives in `Schemes`, and `Superspheres` does not bring it into scope. The
build then stops with

    Cannot resolve @ref for md"[`RVE`](@ref)"
    - No docstring found in doc for binding `...Superspheres.RVE`.
    - Fallback resolution in Main for `RVE` -> `...Schemes.RVE` is only
      allowed for fully qualified names

and the fix is in the last line: write the name out,

    [`RVE`](@ref MeanFieldHomogenization.Schemes.RVE)

which is what `src/Assemblies/as_rve.jl` and `src/Interactions/api.jl` already
do. This script flags every unqualified `@ref` whose target is not reachable
from the module the docstring lives in.

## What it does *not* model, and why it is still worth running

It is a lexical approximation, not Documenter. It reads definitions and
`export` lists with regular expressions, treats `using ..X` as making all of
`X` reachable, and assumes `Core` is reachable everywhere. It does not run
Julia, expand macros, or know which docstrings are actually listed on an API
page — so it can miss a case and it can, in principle, complain about one that
would resolve.

The doc build remains the authority. What this buys is finding the failure in a
second instead of in a CI run several minutes long, and finding it for
docstrings that no API page lists *yet* — those are latent, and break the build
the day someone lists them, far from the change that introduced them.
"""

import collections
import glob
import os
import re
import sys

DEF_PATTERNS = [
    re.compile(r"^\s*function\s+([A-Za-z_][\w!]*)\s*\("),
    re.compile(r"^\s*function\s+([A-Za-z_][\w!]*)\s+end\s*$"),
    re.compile(r"^\s*@[\w.]+\s+function\s+([A-Za-z_][\w!]*)\s*\("),
    re.compile(r"^\s*(?:mutable\s+)?struct\s+([A-Za-z_]\w*)"),
    re.compile(r"^\s*Base\.@kwdef\s+struct\s+([A-Za-z_]\w*)"),
    re.compile(r"^\s*abstract\s+type\s+([A-Za-z_]\w*)"),
    re.compile(r"^\s*const\s+([A-Za-z_]\w*)"),
    re.compile(r"^\s*macro\s+([A-Za-z_]\w*)"),
    re.compile(r"^\s*@inline\s+([A-Za-z_][\w!]*)\s*\("),
    re.compile(r"^\s*([A-Za-z_][\w!]*)\s*\([^)]*\)\s*(?:where[^=]*)?="),
]

REF = re.compile(r"\[`?([A-Za-z_][\w!.]*)`?\]\(@ref\)")


def defined_names(files):
    names = set()
    for path in files:
        text = open(path, encoding="utf-8").read()
        lines = text.splitlines()
        for line in lines:
            for pat in DEF_PATTERNS:
                m = pat.match(line)
                if m:
                    names.add(m.group(1))
        # `export` lists wrap across lines with a trailing comma.
        i = 0
        while i < len(lines):
            m = re.match(r"^\s*export\s+(.*)$", lines[i])
            if m:
                chunk = m.group(1)
                while chunk.rstrip().endswith(",") and i + 1 < len(lines):
                    i += 1
                    chunk += " " + lines[i].strip()
                for name in re.split(r"[,\s]+", chunk):
                    name = name.strip()
                    if name and not name.startswith("#"):
                        names.add(name)
            i += 1
    return names


MODULE = re.compile(r"^module\s+(\w+)\s*$", re.M)


def declared_module(files):
    """The module a directory declares, or `None` if it declares none.

    Two layouts have to be told apart, and getting it wrong makes the check
    useless in opposite ways.

    `MeanFieldHomogenization/src/Superspheres/` declares `module Superspheres`,
    so its files really are a separate namespace and a reference out of it needs
    qualifying. `ChemistryLab/src/equilibrium/` declares nothing: its files are
    `include`d straight into the one top-level module, so every name in the
    package is in reach and there is nothing to qualify. Treating the second
    layout as if it were the first reports every cross-directory reference as
    broken -- 37 of them, all false.

    Note it is the *module* name that matters, not the directory's: a
    `using ..X` names the module.
    """
    for path in files:
        m = MODULE.search(open(path, encoding="utf-8").read())
        if m:
            return m.group(1)
    return None


def main(root="src"):
    groups = collections.defaultdict(list)
    groups["<top level>"] = sorted(glob.glob(os.path.join(root, "*.jl")))
    for entry in sorted(os.listdir(root)):
        path = os.path.join(root, entry)
        if not os.path.isdir(path):
            continue
        files = sorted(glob.glob(path + "/**/*.jl", recursive=True))
        if not files:
            continue
        name = declared_module(files)
        if name is None:
            groups["<top level>"].extend(files)
        else:
            groups[name].extend(files)
    groups = dict(groups)

    defs = {name: defined_names(files) for name, files in groups.items()}

    reachable = collections.defaultdict(set)
    for name, files in groups.items():
        for path in files:
            text = open(path, encoding="utf-8").read()
            for m in re.finditer(r"^\s*(?:using|import)\s+\.\.(\w+)", text, re.M):
                reachable[name].add(m.group(1))

    bad = []
    for name, files in groups.items():
        reach = set(defs[name]) | set(defs.get("Core", set()))
        for other in reachable[name]:
            reach |= defs.get(other, set())
        for path in files:
            for lineno, line in enumerate(
                open(path, encoding="utf-8").read().splitlines(), 1
            ):
                for m in REF.finditer(line):
                    target = m.group(1)
                    if "." in target or target in reach:
                        continue
                    bad.append((path, lineno, name, target))

    if not bad:
        print(f"OK: every unqualified @ref under {root}/ is reachable from its module")
        return 0
    print(f"FAIL: {len(bad)} unqualified @ref out of reach of its own module")
    for path, lineno, mod, target in bad:
        print(f"  {path}:{lineno}  in {mod}: `{target}`")
    print("\nQualify them, e.g. [`RVE`](@ref MeanFieldHomogenization.Schemes.RVE),")
    print("or make them plain code spans if the target is private.")
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "src"))
