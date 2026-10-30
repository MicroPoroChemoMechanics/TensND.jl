#!/usr/bin/env python3
"""Write, or check, the "See also" blocks that link the Theory, the Manual and
the Tutorials of the documentation to one another.

The pages that treat the same subject are listed in `docs/crosslinks.toml`, by
their `@id`. Every page named there gets, right under its title, a block that
links it to the other pages of its subjects, grouped by part and labeled with
their titles, which are read from the pages themselves. A tutorial's block is
written into its script under `scripts/`, the source Literate generates the
page from, and into the generated page.

    python3 .github/scripts/crosslinks.py           # rewrite the blocks
    python3 .github/scripts/crosslinks.py --check   # fail if one is stale

The check runs in the pre-push hook and in CI, so that a page added to, or
moved within, a subject cannot leave its links behind; Documenter fails the
build on an `@id` that no longer exists.
"""

import glob
import os
import re
import sys
import tomllib

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BEGIN = "<!-- crosslinks: written by .github/scripts/crosslinks.py from docs/crosslinks.toml; do not edit -->"
END = "<!-- /crosslinks -->"
PARTS = (("theory", "Theory"), ("manual", "Manual"), ("tutorials", "Tutorials"))
TITLE = re.compile(r"^(?:# )?# \[(?P<title>[^\]]+)\]\(@id (?P<id>[\w-]+)\)[ \t]*$", re.M)


def pages():
    """Every page that carries an `@id` title: its id → (title, files to write)."""
    found = {}
    sources = glob.glob(os.path.join(ROOT, "docs", "src", "**", "*.md"), recursive=True)
    sources += glob.glob(os.path.join(ROOT, "scripts", "*.jl"))
    for path in sources:
        text = open(path, encoding="utf-8").read()
        m = TITLE.search(text)
        if m is None:
            continue
        title, files = found.setdefault(m["id"], (m["title"], []))
        files.append(path)
    return found


def block(lines, script):
    """The block, as Markdown, or as the comments Literate turns into Markdown."""
    body = [BEGIN, '!!! tip "See also"', *(("    " + line) if line else "" for line in lines), END]
    return "\n".join(("# " + line if line else "#") if script else line for line in body)


def links_of(page, subjects, titles):
    """The pages that share a subject with `page`, part by part, in table order."""
    lines = []
    for key, label in PARTS:
        seen = []
        for subject in subjects:
            ids = [i for part, _ in PARTS for i in subject.get(part, [])]
            if page not in ids:
                continue
            for other in subject.get(key, []):
                if other != page and other not in seen:
                    seen.append(other)
        if seen:
            lines.append(f"**{label}**: " + " · ".join(f"[{titles[i]}](@ref {i})" for i in seen))
            lines.append("")
    return lines[:-1]


def rewrite(path, new_block):
    """The file with its block replaced, or inserted under its title, one blank
    line (a bare `#` in a script) on each side."""
    text = open(path, encoding="utf-8").read()
    script = path.endswith(".jl")
    begin = ("# " + BEGIN) if script else BEGIN
    end = ("# " + END) if script else END
    if begin in text:
        before, rest = text[: text.index(begin)], text[text.index(end) + len(end):]
    else:
        m = TITLE.search(text)
        before, rest = text[: m.end()], text[m.end():]
    before = before.rstrip("\n")
    rest = rest.lstrip("\n")
    if script:
        before = before.removesuffix("\n#")
        rest = rest[2:] if rest.startswith("#\n") else rest
        return before + "\n#\n" + new_block + "\n#\n" + rest
    return before + "\n\n" + new_block + "\n\n" + rest


def main(check):
    table = tomllib.load(open(os.path.join(ROOT, "docs", "crosslinks.toml"), "rb"))
    subjects = table["subject"]
    found = pages()
    named = {i for s in subjects for part, _ in PARTS for i in s.get(part, [])}
    missing = sorted(named - found.keys())
    if missing:
        sys.exit(f"crosslinks: no page has the @id {', '.join(missing)}")
    titles = {i: title for i, (title, _) in found.items()}
    stale = []
    for page in sorted(named):
        lines = links_of(page, subjects, titles)
        for path in found[page][1]:
            new = rewrite(path, block(lines, path.endswith(".jl")))
            if new != open(path, encoding="utf-8").read():
                stale.append(os.path.relpath(path, ROOT))
                if not check:
                    open(path, "w", encoding="utf-8").write(new)
    if check and stale:
        sys.exit("crosslinks: stale \"See also\" blocks, run .github/scripts/crosslinks.py:\n  " + "\n  ".join(stale))
    print(f"crosslinks: {len(named)} pages {'checked' if check else 'written'}"
          + (f", {len(stale)} updated" if stale and not check else ""))


if __name__ == "__main__":
    main("--check" in sys.argv[1:])
