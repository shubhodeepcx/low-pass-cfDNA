#!/usr/bin/env python3
"""Guard the prose: no long dashes anywhere in the writing.

Long dashes read as machine written and this project does not use them. The
check runs in CI so the rule holds without anyone having to remember it.

Two things are deliberately not flagged. Markdown table rules and the dashed
underlines of numpydoc sections are structure rather than punctuation, so a
run of three or more dashes only counts inside LaTeX, where it really does
typeset an em dash. And the IEEE document class emits its own "Abstract" and
"Index Terms" rules, which are part of the template rather than anything
written here.

The banned characters are built with ``chr`` so that this file does not trip
its own check.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SUFFIXES = {".md", ".py", ".tex", ".bib", ".txt", ".yml", ".yaml", ".toml"}
SKIP_DIRS = {".git", "__pycache__", ".pytest_cache", ".ruff_cache", "build",
             "dist", ".venv", "venv"}
SKIP_NAMES = {"ieee.csl", "main.bbl", "check_prose.py"}

BANNED_CHARS = {
    chr(0x2014): "em dash",
    chr(0x2013): "en dash",
    chr(0x2015): "horizontal bar",
}

#: In LaTeX, three dashes typeset an em dash. A line that is only dashes is a
#: rule or a separator, not punctuation.
LATEX_EM_DASH = re.compile(r"(?<!-)---(?!-)")
ONLY_DASHES = re.compile(r"^[\s|:-]*$")


def source_files() -> list[Path]:
    found = []
    for path in sorted(ROOT.rglob("*")):
        if not path.is_file() or path.suffix not in SUFFIXES:
            continue
        if path.name in SKIP_NAMES:
            continue
        if any(part in SKIP_DIRS for part in path.parts):
            continue
        found.append(path)
    return found


def problems_in(path: Path) -> list[str]:
    try:
        text = path.read_text(encoding="utf-8")
    except (UnicodeDecodeError, OSError):
        return []

    relative = path.relative_to(ROOT)
    found: list[str] = []
    for number, line in enumerate(text.splitlines(), start=1):
        for char, name in BANNED_CHARS.items():
            if char in line:
                found.append(f"{relative}:{number}: {name}: {line.strip()[:80]}")
        if path.suffix == ".tex" and not ONLY_DASHES.match(line):
            if LATEX_EM_DASH.search(line):
                found.append(f"{relative}:{number}: LaTeX em dash: "
                             f"{line.strip()[:80]}")
    return found


def main() -> int:
    files = source_files()
    problems = [problem for path in files for problem in problems_in(path)]
    if problems:
        print(f"{len(problems)} prose problems found:")
        for problem in problems:
            print(f"  {problem}")
        return 1
    print(f"prose check clean across {len(files)} files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
