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

## The second rule: an `@ref` needs a *docstring*, not just a binding

Reachability is necessary and not sufficient. Documenter links an `@ref` to a
**docstring**, so a name that exists and is even exported, but that carries no
docstring at all, fails the same way:

    Cannot resolve @ref for md"[`FunctionRetention`](@ref)"
    - No docstring found in doc for binding `ChemistryLab.FunctionRetention`.

`checkdocs` does not catch this. It verifies that the docstrings a package
*has* are all published; a binding with none has nothing to publish and nothing
to complain about. Only the reference fails, and only at build time.

So this script also collects which names carry a docstring, by taking the first
code line after each triple-quoted block, and flags an `@ref` whose target is
defined but undocumented. That sweeps in a second trap for free: Julia detaches
a docstring silently when a **comment** sits between the closing quotes and the
definition, so an orphaned docstring reads here as an undocumented name.
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

# A name on a line of its own also takes a docstring: a doc block then `foo`,
# which is how a const or a function with no method here is documented.
BARE_NAME = re.compile(r"^\s*([A-Za-z_][\w!]*)\s*$")

# What the line after a doc block documents. Deliberately looser than
# `DEF_PATTERNS`, which is written to enumerate definitions: here a miss turns
# into a FAIL on a name that is plainly documented, so over-collecting is the
# safe direction. `DEF_PATTERNS` cannot be reused as it stands — its
# assignment-form pattern reads the signature as `[^)]*`, which stops at the
# first `)` and so misses every signature with a nested call in it, such as
# `tens_Id4(::Val{dim} = Val(3), ::Val{T} = Val(Sym)) where {dim, T} = ...`.
# That was 14 false positives on TensND alone.
DOC_TARGETS = [
    re.compile(
        r"^\s*(?:@\w[\w.]*\s+)?(?:mutable\s+)?"
        r"(?:function|struct|abstract\s+type|primitive\s+type|macro|const|"
        r"module|baremodule)\s+([A-Za-z_][\w!]*)"
    ),
    re.compile(r"^\s*(?:@\w[\w.]*\s+)?(?:Base\.)?([A-Za-z_][\w!]*)\s*[({=]"),
    BARE_NAME,
]


def documented_names(files):
    """Names carrying a docstring, paired block by block.

    Three shapes have to be told apart, and getting the third one wrong throws
    the pairing off for the rest of the file — every later closing line then
    reads as an opener:

        \"\"\"docstring\"\"\"                  the target is the next code line
        f(x) = ...

        @doc \"\"\"                        the target is on the CLOSING line,
        docstring                            which is how a type alias or a
        \"\"\" PhaseQuantities                binding with no definition of its
                                             own is documented

        \"\"\"docstring\"\"\" name             both on one line

    A **comment** between the closing quotes and the definition detaches the
    docstring in Julia, silently, so this stops at one and leaves the name out —
    which is the behavior wanted, since Documenter will not find it either.
    """
    q = '"' * 3
    names = set()
    for path in files:
        lines = open(path, encoding="utf-8").read().splitlines()
        i = 0
        while i < len(lines):
            head = lines[i].strip()
            if head.startswith("@doc"):
                head = head[4:].lstrip()
            if head.startswith("raw" + q):
                head = head[3:]
            if not head.startswith(q):
                i += 1
                continue
            rest = head[3:]
            if q in rest:                       # opens and closes on one line
                close, tail = i, rest.split(q, 1)[1]
            else:
                close = i + 1
                while close < len(lines) and q not in lines[close]:
                    close += 1
                tail = lines[close].split(q, 1)[1] if close < len(lines) else ""
            named = BARE_NAME.match(tail)       # `\"\"\" PhaseQuantities`
            if named:
                names.add(named.group(1))
                i = close + 1
                continue
            code = close + 1
            while code < len(lines) and lines[code].strip() == "":
                code += 1
            if code < len(lines):
                for pat in DOC_TARGETS:
                    m = pat.match(lines[code])
                    if m:
                        names.add(m.group(1))
                        break
            i = code + 1
    return names


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
    documented = set()
    for files in groups.values():
        documented |= documented_names(files)

    reachable = collections.defaultdict(set)
    for name, files in groups.items():
        for path in files:
            text = open(path, encoding="utf-8").read()
            for m in re.finditer(r"^\s*(?:using|import)\s+\.\.(\w+)", text, re.M):
                reachable[name].add(m.group(1))

    bad = []
    undocumented = []
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
                    if "." in target:
                        continue
                    if target not in reach:
                        bad.append((path, lineno, name, target))
                    elif target not in documented:
                        undocumented.append((path, lineno, name, target))

    if not bad and not undocumented:
        print(
            f"OK: every unqualified @ref under {root}/ is reachable from its "
            "module and has a docstring"
        )
        return 0
    if bad:
        print(f"FAIL: {len(bad)} unqualified @ref out of reach of its own module")
        for path, lineno, mod, target in bad:
            print(f"  {path}:{lineno}  in {mod}: `{target}`")
        print("\nQualify them, e.g. [`RVE`](@ref MeanFieldHomogenization.Schemes.RVE),")
        print("or make them plain code spans if the target is private.")
    if undocumented:
        print(f"FAIL: {len(undocumented)} @ref to a name that carries no docstring")
        for path, lineno, mod, target in undocumented:
            print(f"  {path}:{lineno}  in {mod}: `{target}`")
        print("\nDocumenter links an @ref to a docstring, not to a binding, so it")
        print("cannot resolve these. Give the target a docstring, or make the")
        print("reference a plain code span. A comment between the closing quotes")
        print("and the definition detaches a docstring silently and reads here as")
        print("a missing one.")
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "src"))
