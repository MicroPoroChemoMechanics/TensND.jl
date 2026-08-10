#!/usr/bin/env python3
"""US/UK English spelling scanner and converter for Julia package prose.

Scans docstrings, Julia comments, and Markdown docs -- never Julia code
tokens, identifiers, or kwargs -- for US/UK spelling variants, and can
convert prose to a chosen convention. Stdlib only, no dependencies.

Usage:
    spelling_tool.py scan    <path>... [--json] [--verbose] [--ignore-file PATH]
    spelling_tool.py convert <path>... --to {us,uk} [--dry-run]
                              [--include-context-sensitive] [--ignore-file PATH]
    spelling_tool.py check   <path>... [--convention {us,uk}] [--ignore-file PATH]
    spelling_tool.py init    <path> --convention {us,uk}

See SKILL.md in the same directory for the full design rationale.
"""
from __future__ import annotations

import argparse
import fnmatch
import json
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent

# What a whole-tree scan covers. This list is the enforcement boundary: a file
# outside it is never seen by `check .`, hence never by CI. It used to hold only
# `src`, `docs/src` and `scripts`, which left README.md, CHANGELOG.md, the test
# suite, weak-dependency extensions and `tools/` unchecked — the README being
# the most-read file in the repository. It also made the pre-commit hook
# stricter than CI, since the hook passes the staged files themselves.
DEFAULT_INCLUDE_GLOBS = [
    "*.md",                 # README, CHANGELOG, CONTRIBUTING at the root
    "src/**/*.jl",
    "ext/**/*.jl",
    "test/**/*.jl",
    "docs/**/*.md",
    "docs/**/*.jl",
    "scripts/**/*.jl",
    "tools/**/*.jl",
    "tools/**/*.py",
    "tools/**/*.md",
]
DEFAULT_EXCLUDE_GLOBS = [
    "**/generated/**",
    "docs/build/**",
    "**/.git/**",
]

DOCSTRING_KEYWORD_RE = re.compile(
    r"^\s*(?:@\w+(?:\([^)]*\))?\s+)*"
    r"(function|struct|mutable\s+struct|module|baremodule|macro|"
    r"abstract\s+type|primitive\s+type|const)\b"
)

# -ize/-ise and -yze/-yse family suffixes (analyze/analyse, paralyze/paralyse);
# order doesn't matter here, _ALL_SUFFIXES below sorts by length for correct
# alternation matching.
SUFFIX_PAIRS = [
    ("ization", "isation"),
    ("izer", "iser"),
    ("izing", "ising"),
    ("izes", "ises"),
    ("ized", "ised"),
    ("ize", "ise"),
    ("yzer", "yser"),
    ("yzing", "ysing"),
    ("yzes", "yses"),
    ("yzed", "ysed"),
    ("yze", "yse"),
]
US_TO_UK_SUFFIX = dict(SUFFIX_PAIRS)
UK_TO_US_SUFFIX = {v: k for k, v in SUFFIX_PAIRS}
_ALL_SUFFIXES = sorted(
    list(US_TO_UK_SUFFIX) + list(UK_TO_US_SUFFIX), key=len, reverse=True
)
IZE_ISE_RE = re.compile(
    r"\b([A-Za-z]{3,}?)(" + "|".join(_ALL_SUFFIXES) + r")\b", re.IGNORECASE
)


# --------------------------------------------------------------------------
# Data loading
# --------------------------------------------------------------------------

def load_json(name: str) -> dict:
    with open(SCRIPT_DIR / name, encoding="utf-8") as f:
        return json.load(f)


def silent_e_inflect(base: str) -> list[str]:
    """advise -> [advise, advised, advises, advising]"""
    assert base.endswith("e")
    stem = base[:-1]
    return [base, base + "d", base + "s", stem + "ing"]


def build_always_ise() -> tuple[set[str], set[str]]:
    data = load_json("always_ise.json")
    always_ise_forms = set()
    for base in data["always_ise_base"]:
        for form in silent_e_inflect(base):
            always_ise_forms.add(form.lower())
    invariant_words = {w.lower() for w in data["invariant_words"]}
    return always_ise_forms, invariant_words


def build_pair_index() -> tuple[dict, dict]:
    """Returns (simple_index, context_index): word.lower() -> metadata dict."""
    data = load_json("spelling_pairs.json")
    simple_index: dict[str, dict] = {}
    context_index: dict[str, dict] = {}

    def add(entries: list[dict], target: dict[str, dict]) -> None:
        for entry in entries:
            us_forms, uk_forms = entry["us"], entry["uk"]
            if len(us_forms) != len(uk_forms):
                raise ValueError(
                    f"Misaligned us/uk forms in spelling_pairs.json: {entry}"
                )
            for us_word, uk_word in zip(us_forms, uk_forms):
                us_l, uk_l = us_word.lower(), uk_word.lower()
                if us_l == uk_l:
                    continue  # identical in both dialects, not a variant
                target[us_l] = {"dialect": "us", "us": us_l, "uk": uk_l}
                target[uk_l] = {"dialect": "uk", "us": us_l, "uk": uk_l}

    add(data.get("simple", []), simple_index)
    add(data.get("context_sensitive", []), context_index)
    return simple_index, context_index


# --------------------------------------------------------------------------
# Length-preserving masking utilities (offsets stay valid against original text)
# --------------------------------------------------------------------------

def _blank(text: str, start: int, end: int) -> str:
    return text[:start] + "".join(" " if c != "\n" else "\n" for c in text[start:end]) + text[end:]


def mask_frontmatter(text: str) -> str:
    if not text.startswith("---\n") and not text.startswith("---\r\n"):
        return text
    m = re.search(r"\n---[ \t]*\r?\n", text)
    if not m:
        return text
    return _blank(text, 0, m.end())


def mask_fenced_code(text: str) -> str:
    lines = text.splitlines(keepends=True)
    out = []
    in_fence = False
    fence_marker = None
    offset = 0
    for line in lines:
        stripped = line.strip()
        is_fence_line = stripped.startswith("```") or stripped.startswith("~~~")
        if not in_fence and is_fence_line:
            in_fence = True
            fence_marker = stripped[:3]
            out.append(_blank(line, 0, len(line)))
        elif in_fence and stripped.startswith(fence_marker):
            in_fence = False
            fence_marker = None
            out.append(_blank(line, 0, len(line)))
        elif in_fence:
            out.append(_blank(line, 0, len(line)))
        else:
            out.append(line)
        offset += len(line)
    return "".join(out)


def mask_inline_code(text: str) -> str:
    def repl(m: re.Match) -> str:
        return " " * (m.end() - m.start())

    return re.sub(r"`[^`\n]+?`", repl, text)


def mask_link_urls(text: str) -> str:
    def repl(m: re.Match) -> str:
        return m.group(0)[: m.start(1) - m.start(0)] + " " * (
            m.end(1) - m.start(1)
        ) + m.group(0)[m.end(1) - m.start(0):]

    return re.sub(r"\]\(([^)]*)\)", repl, text)


def mask_literals(text: str, literals: list[str]) -> str:
    for lit in literals:
        if not lit:
            continue
        start = 0
        while True:
            idx = text.find(lit, start)
            if idx == -1:
                break
            text = _blank(text, idx, idx + len(lit))
            start = idx + len(lit)
    return text


def mask_markdown(text: str, literals: list[str]) -> str:
    text = mask_frontmatter(text)
    text = mask_fenced_code(text)
    text = mask_inline_code(text)
    text = mask_link_urls(text)
    text = mask_literals(text, literals)
    return text


# --------------------------------------------------------------------------
# Julia prose extraction: docstrings + comments, everything else left alone
# --------------------------------------------------------------------------

def _line_has_only_comment_or_blank(line: str) -> bool:
    stripped = line.strip()
    return not stripped or stripped.startswith("#")


def _next_meaningful_line(text: str, from_offset: int) -> str | None:
    """First non-blank line after from_offset. Deliberately does NOT skip
    over comment lines: a real docstring is immediately followed by the
    definition it documents (blank lines tolerated, but not a comment) --
    treating an interleaved comment as disqualifying keeps the classifier
    conservative and avoids misattributing an unrelated later definition to
    an unrelated earlier string literal."""
    for line in text[from_offset:].splitlines():
        if not line.strip():
            continue
        return line
    return None


def _is_first_statement(text: str, block_start: int) -> bool:
    before = text[:block_start]
    for line in before.splitlines():
        if line.startswith("#!"):
            continue
        if _line_has_only_comment_or_blank(line):
            continue
        return False
    return True


def find_triple_quote_blocks(text: str) -> list[tuple[int, int, int, int]]:
    """(block_start, block_end, interior_start, interior_end) for each \"\"\"...\"\"\"."""
    blocks = []
    for m in re.finditer(r'"""(.*?)"""', text, re.DOTALL):
        blocks.append((m.start(), m.end(), m.start() + 3, m.end() - 3))
    return blocks


def extract_jl_prose_spans(text: str, literals: list[str]) -> list[tuple[int, int, str]]:
    """Returns (start, end, masked_text_same_length) spans of prose to scan."""
    spans = []
    triple_blocks = find_triple_quote_blocks(text)

    for block_start, block_end, interior_start, interior_end in triple_blocks:
        next_line = _next_meaningful_line(text, block_end)
        is_docstring = _is_first_statement(text, block_start) or (
            next_line is not None and bool(DOCSTRING_KEYWORD_RE.match(next_line))
        )
        if not is_docstring:
            continue
        interior = text[interior_start:interior_end]
        masked = mask_markdown(interior, literals)
        spans.append((interior_start, interior_end, masked))

    # Blank out every triple-quoted block (docstring or not) before scanning
    # for comments, so a '#' inside string content is never treated as one.
    working = text
    for block_start, block_end, _, _ in triple_blocks:
        working = _blank(working, block_start, block_end)

    comment_spans = _find_comment_spans(working)
    for start, end in comment_spans:
        comment_text = text[start:end]
        masked = mask_inline_code(comment_text)
        masked = mask_literals(masked, literals)
        spans.append((start, end, masked))

    return spans


def _find_comment_spans(working_text: str) -> list[tuple[int, int]]:
    spans = []
    offset = 0
    for idx, line in enumerate(working_text.splitlines(keepends=True)):
        content = line.rstrip("\n").rstrip("\r")
        if idx == 0 and content.startswith("#!"):
            offset += len(line)
            continue
        in_str = False
        i = 0
        n = len(content)
        while i < n:
            c = content[i]
            if in_str and c == "\\":
                i += 2
                continue
            if c == '"':
                in_str = not in_str
                i += 1
                continue
            if c == "#" and not in_str:
                spans.append((offset + i, offset + n))
                break
            i += 1
        offset += len(line)
    return spans


# --------------------------------------------------------------------------
# Matching
# --------------------------------------------------------------------------

@dataclass
class Occurrence:
    start: int
    end: int
    word: str
    dialect: str  # "us" or "uk"
    family: str  # "us_form/uk_form"
    tier: str  # "simple" or "context_sensitive"


def _build_pair_regex(index: dict) -> re.Pattern | None:
    if not index:
        return None
    keys = sorted(index.keys(), key=len, reverse=True)
    return re.compile(r"\b(" + "|".join(re.escape(k) for k in keys) + r")\b", re.IGNORECASE)


def find_pair_occurrences(masked: str, base_offset: int, index: dict, regex: re.Pattern | None) -> list[Occurrence]:
    if regex is None:
        return []
    out = []
    for m in regex.finditer(masked):
        key = m.group(0).lower()
        info = index.get(key)
        if info is None:
            continue
        out.append(
            Occurrence(
                start=base_offset + m.start(),
                end=base_offset + m.end(),
                word=m.group(0),
                dialect=info["dialect"],
                family=f"{info['us']}/{info['uk']}",
                tier=info.get("tier", "simple"),
            )
        )
    return out


def find_ize_ise_occurrences(
    masked: str, base_offset: int, always_ise_forms: set[str], invariant_words: set[str]
) -> list[Occurrence]:
    out = []
    for m in IZE_ISE_RE.finditer(masked):
        stem, suffix = m.group(1), m.group(2).lower()
        word = m.group(0)
        word_l = word.lower()
        if word_l in invariant_words or word_l in always_ise_forms:
            continue
        if suffix in UK_TO_US_SUFFIX and stem.lower().endswith("w"):
            continue  # -wise family: elementwise, componentwise, otherwise...
        if suffix in US_TO_UK_SUFFIX:
            dialect = "us"
            us_suffix, uk_suffix = suffix, US_TO_UK_SUFFIX[suffix]
        elif suffix in UK_TO_US_SUFFIX:
            dialect = "uk"
            uk_suffix, us_suffix = suffix, UK_TO_US_SUFFIX[suffix]
        else:
            continue
        stem_l = stem.lower()
        out.append(
            Occurrence(
                start=base_offset + m.start(),
                end=base_offset + m.end(),
                word=word,
                dialect=dialect,
                family=f"{stem_l}{us_suffix}/{stem_l}{uk_suffix}",
                tier="simple",
            )
        )
    return out


def preserve_case(original: str, replacement: str) -> str:
    if original.isupper() and len(original) > 1:
        return replacement.upper()
    if original[:1].isupper():
        return replacement[:1].upper() + replacement[1:]
    return replacement


def occurrence_replacement(occ: Occurrence, target: str) -> str | None:
    """The replacement text to convert this occurrence to `target` dialect, or
    None if it is already in that dialect."""
    if occ.dialect == target:
        return None
    us_form, uk_form = occ.family.split("/")
    target_form = us_form if target == "us" else uk_form
    return preserve_case(occ.word, target_form)


# --------------------------------------------------------------------------
# File scanning
# --------------------------------------------------------------------------

def scan_text(path: Path, text: str, literals: list[str], simple_regex, context_regex, simple_index, context_index, always_ise_forms, invariant_words) -> list[Occurrence]:
    occurrences: list[Occurrence] = []
    if path.suffix == ".md":
        masked = mask_markdown(text, literals)
        spans = [(0, len(text), masked)]
    elif path.suffix == ".jl":
        spans = extract_jl_prose_spans(text, literals)
    else:
        return []

    for start, _end, masked in spans:
        occurrences.extend(find_pair_occurrences(masked, start, simple_index, simple_regex))
        occurrences.extend(find_pair_occurrences(masked, start, context_index, context_regex))
        occurrences.extend(find_ize_ise_occurrences(masked, start, always_ise_forms, invariant_words))
    occurrences.sort(key=lambda o: o.start)
    return occurrences


# --------------------------------------------------------------------------
# File discovery + ignore handling
# --------------------------------------------------------------------------

def load_ignore(root: Path, extra_ignore_file: str | None) -> tuple[list[str], list[str]]:
    globs: list[str] = []
    literals: list[str] = []
    candidates = []
    default = root / ".spelling-ignore"
    if default.exists():
        candidates.append(default)
    if extra_ignore_file:
        candidates.append(Path(extra_ignore_file))
    for p in candidates:
        for line in p.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("literal:"):
                literals.append(line[len("literal:"):].strip())
            else:
                globs.append(line)
    return globs, literals


def discover_files(root: Path, extra_ignore_file: str | None) -> tuple[list[Path], list[str]]:
    if root.is_file():
        _, literals = load_ignore(root.parent, extra_ignore_file)
        return [root], literals

    ignore_globs, literals = load_ignore(root, extra_ignore_file)
    exclude_globs = DEFAULT_EXCLUDE_GLOBS + ignore_globs

    found: dict[Path, None] = {}
    for pattern in DEFAULT_INCLUDE_GLOBS:
        for f in root.glob(pattern):
            if f.is_file():
                found[f] = None

    result = []
    for f in found:
        rel = f.relative_to(root).as_posix()
        if any(fnmatch.fnmatch(rel, g) for g in exclude_globs):
            continue
        result.append(f)
    return sorted(result), literals


# --------------------------------------------------------------------------
# Config (.spelling.json)
# --------------------------------------------------------------------------

def config_path(root: Path) -> Path:
    return (root if root.is_dir() else root.parent) / ".spelling.json"


def read_convention(root: Path) -> str | None:
    cfg = config_path(root)
    if not cfg.exists():
        return None
    try:
        data = json.loads(cfg.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None
    return data.get("convention")


def write_convention(root: Path, convention: str) -> Path:
    cfg = config_path(root)
    cfg.write_text(
        json.dumps({"convention": convention, "generated_by": "spelling_tool.py init", "version": 1}, indent=2)
        + "\n",
        encoding="utf-8",
    )
    return cfg


# --------------------------------------------------------------------------
# Commands
# --------------------------------------------------------------------------

def _load_engines():
    simple_index, context_index = build_pair_index()
    always_ise_forms, invariant_words = build_always_ise()
    simple_regex = _build_pair_regex(simple_index)
    context_regex = _build_pair_regex(context_index)
    return simple_index, context_index, simple_regex, context_regex, always_ise_forms, invariant_words


def _scan_root(root: Path, ignore_file: str | None, engines) -> dict[Path, list[Occurrence]]:
    simple_index, context_index, simple_regex, context_regex, always_ise_forms, invariant_words = engines
    files, literals = discover_files(root, ignore_file)
    results = {}
    for f in files:
        text = f.read_text(encoding="utf-8", errors="replace")
        occs = scan_text(
            f, text, literals, simple_regex, context_regex, simple_index, context_index,
            always_ise_forms, invariant_words,
        )
        if occs:
            results[f] = occs
    return results


def cmd_scan(args: argparse.Namespace) -> int:
    engines = _load_engines()
    report = {}
    for path_str in args.paths:
        root = Path(path_str).resolve()
        per_file = _scan_root(root, args.ignore_file, engines)
        if not per_file:
            continue
        families: dict[str, dict[str, int]] = defaultdict(lambda: {"us": 0, "uk": 0, "tier": "simple"})
        per_file_json = {}
        for f, occs in per_file.items():
            rel = str(f.relative_to(root)) if f.is_relative_to(root) else str(f)
            file_families: dict[str, dict[str, int]] = defaultdict(lambda: {"us": 0, "uk": 0})
            for o in occs:
                families[o.family]["tier"] = o.tier
                families[o.family][o.dialect] += 1
                file_families[o.family][o.dialect] += 1
            per_file_json[rel] = {
                fam: counts for fam, counts in file_families.items()
            }
            if args.verbose:
                per_file_json[rel]["_lines"] = [
                    {"line": f.read_text(encoding="utf-8", errors="replace").count("\n", 0, o.start) + 1,
                     "word": o.word, "family": o.family, "dialect": o.dialect}
                    for o in occs
                ]
        report[str(root)] = {"files": per_file_json, "families": families}

    if args.json:
        print(json.dumps(report, indent=2))
        return 0

    if not report:
        print("No US/UK spelling mixing detected.")
        return 0

    for root_str, data in report.items():
        n_files = len(data["files"])
        n_families = len(data["families"])
        print(f"\n{root_str}  ({n_files} file(s) with mixed spelling, {n_families} word famil{'y' if n_families == 1 else 'ies'})")
        net_us = net_uk = 0
        for fam, counts in sorted(data["families"].items()):
            us_n, uk_n = counts["us"], counts["uk"]
            net_us += us_n
            net_uk += uk_n
            tag = " (context-sensitive, manual review)" if counts["tier"] == "context_sensitive" else ""
            print(f"  {fam:35s} {us_n:4d} US / {uk_n:4d} UK{tag}")
        dominant = "US" if net_us >= net_uk else "UK"
        print(f"  DOMINANT: {dominant}  (net {abs(net_us - net_uk)} occurrences across {n_families} families)")
        if args.verbose:
            for rel, families in data["files"].items():
                lines = families.get("_lines")
                if not lines:
                    continue
                print(f"  -- {rel}")
                for entry in lines:
                    print(f"     L{entry['line']:<5d} {entry['dialect'].upper():2s}  {entry['word']!r}  ({entry['family']})")
    return 0


def cmd_convert(args: argparse.Namespace) -> int:
    engines = _load_engines()
    total_files = 0
    total_replacements = 0
    total_context_skipped = 0
    for path_str in args.paths:
        root = Path(path_str).resolve()
        per_file = _scan_root(root, args.ignore_file, engines)
        for f, occs in per_file.items():
            to_convert = []
            skipped_context = 0
            for o in occs:
                if o.tier == "context_sensitive" and not args.include_context_sensitive:
                    if o.dialect != args.to:
                        skipped_context += 1
                    continue
                repl = occurrence_replacement(o, args.to)
                if repl is not None:
                    to_convert.append((o, repl))
            if not to_convert and skipped_context == 0:
                continue
            total_context_skipped += skipped_context
            if not to_convert:
                continue
            text = f.read_text(encoding="utf-8")
            new_text_parts = []
            cursor = 0
            for o, repl in sorted(to_convert, key=lambda pair: pair[0].start):
                new_text_parts.append(text[cursor:o.start])
                new_text_parts.append(repl)
                cursor = o.end
            new_text_parts.append(text[cursor:])
            new_text = "".join(new_text_parts)

            rel = str(f.relative_to(root)) if f.is_relative_to(root) else str(f)
            print(f"  {rel}   {len(to_convert)} replacement(s)")
            total_files += 1
            total_replacements += len(to_convert)
            if not args.dry_run:
                f.write_text(new_text, encoding="utf-8")

    prefix = "DRY RUN -- no files written. " if args.dry_run else ""
    print(f"\n{prefix}{total_files} file(s) changed, {total_replacements} replacement(s) total.")
    if total_replacements:
        print("Review with: git diff")
    if total_context_skipped:
        print(
            f"{total_context_skipped} context-sensitive occurrence(s) skipped "
            "(license/practice) -- rerun with --include-context-sensitive to force (not recommended)."
        )
    return 0


def cmd_check(args: argparse.Namespace) -> int:
    engines = _load_engines()
    overall_ok = True
    for path_str in args.paths:
        root = Path(path_str).resolve()
        convention = args.convention or read_convention(root)
        if not convention:
            print(f"error: no convention configured for {root} (run 'init' or pass --convention)", file=sys.stderr)
            return 2
        other = "uk" if convention == "us" else "us"
        per_file = _scan_root(root, args.ignore_file, engines)
        violations = {}
        for f, occs in per_file.items():
            bad = [o for o in occs if o.tier == "simple" and o.dialect == other]
            if bad:
                violations[f] = bad
        if violations:
            overall_ok = False
            print(f"FAIL: {root} has {sum(len(v) for v in violations.values())} occurrence(s) of {other.upper()} spelling (convention is {convention.upper()})")
            for f, occs in sorted(violations.items()):
                rel = str(f.relative_to(root)) if f.is_relative_to(root) else str(f)
                words = ", ".join(sorted({o.word for o in occs}))
                print(f"  {rel}: {words}")
        else:
            n_files = sum(1 for f in discover_files(root, args.ignore_file)[0])
            print(f"OK: {root} is consistently {convention.upper()} English ({n_files} files checked)")
    return 0 if overall_ok else 1


def cmd_init(args: argparse.Namespace) -> int:
    root = Path(args.path).resolve()
    cfg = write_convention(root, args.convention)
    print(f"Wrote {cfg} (convention={args.convention})")
    return 0


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    p_scan = sub.add_parser("scan", help="Report US/UK spelling mixing (read-only).")
    p_scan.add_argument("paths", nargs="+")
    p_scan.add_argument("--json", action="store_true")
    p_scan.add_argument("--verbose", action="store_true")
    p_scan.add_argument("--ignore-file")
    p_scan.set_defaults(func=cmd_scan)

    p_convert = sub.add_parser("convert", help="Rewrite prose to the chosen convention.")
    p_convert.add_argument("paths", nargs="+")
    p_convert.add_argument("--to", choices=["us", "uk"], required=True)
    p_convert.add_argument("--dry-run", action="store_true")
    p_convert.add_argument("--include-context-sensitive", action="store_true")
    p_convert.add_argument("--ignore-file")
    p_convert.set_defaults(func=cmd_convert)

    p_check = sub.add_parser("check", help="Exit non-zero if the non-configured convention is found (for CI).")
    p_check.add_argument("paths", nargs="+")
    p_check.add_argument("--convention", choices=["us", "uk"])
    p_check.add_argument("--ignore-file")
    p_check.set_defaults(func=cmd_check)

    p_init = sub.add_parser("init", help="Write .spelling.json at the repo root.")
    p_init.add_argument("path")
    p_init.add_argument("--convention", choices=["us", "uk"], required=True)
    p_init.set_defaults(func=cmd_init)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
