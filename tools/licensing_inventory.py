# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT
#
# REUSE-IgnoreStart -- below this line SPDX tags appear as data
# (policy constants, evidence strings), not as declarations for this
# file. See docs/licensing.md.

"""Regenerate the per-file licensing inventory.

    python tools/licensing_inventory.py

Writes ``docs/licensing-inventory.csv``: one row per in-scope Python
file, recording what the audit concluded and what is in the file now.
The ``evidence`` column says *why* a file is categorised as it is, so a
reviewer can re-check the reasoning instead of trusting the label.
"""

from __future__ import annotations

import ast
import csv
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from license_headers import (  # noqa: E402
    LICENSE_EXCEPTIONS,
    REPO_ROOT,
    check_text,
    docstring_notices,
    parse_header,
    python_files,
    read_text,
)

OUTPUT = os.path.join(REPO_ROOT, "docs", "licensing-inventory.csv")

FIELDS = (
    "path",
    "category",
    "authorship",
    "license",
    "copyright",
    "evidence",
    "action",
    "status",
)

FIRST_PARTY_EVIDENCE = (
    "LICENSE and pyproject.toml declare MIT for this repository; no "
    "upstream notice, licence block or copied source found in the file"
)
RETIRED_EVIDENCE = (
    "carried a personal '(c) YEAR, holder' notice in its docstring "
    "before the 2026-09 audit; the maintainer confirmed copyright is "
    "institutional, so the claim was retired in favour of the holder in "
    "LICENSE and the name kept as an authorship line"
)

# `Author: Name (YEAR)` -- attribution, deliberately not a copyright
# claim. See docs/licensing.md.
_AUTHOR_RE = re.compile(r"^Author:\s*(.+?)\s*$", re.MULTILINE)


def authorship(text: str) -> list[str]:
    try:
        doc = ast.get_docstring(ast.parse(text))
    except SyntaxError:
        return []
    return _AUTHOR_RE.findall(doc or "")


def rows(repo_root: str = REPO_ROOT):
    for rel in python_files(repo_root):
        text = read_text(os.path.join(repo_root, rel))
        head = parse_header(text)
        authors = authorship(text)
        claims = docstring_notices(text)
        problems = check_text(rel, text)
        if rel in LICENSE_EXCEPTIONS:
            expression, evidence = LICENSE_EXCEPTIONS[rel]
            category = "third-party"
        else:
            expression = head.licenses[0] if head.licenses else ""
            evidence = RETIRED_EVIDENCE if authors else FIRST_PARTY_EVIDENCE
            category = "first-party"
        yield {
            "path": rel,
            "category": category,
            "authorship": "; ".join(authors) or "none",
            "license": expression,
            "copyright": " | ".join(head.copyrights),
            "evidence": evidence,
            "action": (
                "SPDX header added; personal copyright claim retired, "
                "authorship kept"
                if authors
                else "SPDX header added"
            ),
            "status": (
                "compliant"
                if not problems and not claims
                else "; ".join(problems)
            ),
        }


def write(output: str = OUTPUT, repo_root: str = REPO_ROOT) -> int:
    data = list(rows(repo_root))
    os.makedirs(os.path.dirname(output), exist_ok=True)
    with open(output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(data)

    def tally(field: str, value: str, equal: bool = True) -> int:
        return sum((row[field] == value) is equal for row in data)

    compliant = tally("status", "compliant")
    print(f"{len(data)} files -> {output}")
    print(f"  first-party : {tally('category', 'first-party')}")
    print(f"  third-party : {tally('category', 'third-party')}")
    print(f"  authored    : {tally('authorship', 'none', equal=False)}")
    print(f"  compliant   : {compliant}/{len(data)}")
    return 0 if compliant == len(data) else 1


if __name__ == "__main__":
    raise SystemExit(write())

# REUSE-IgnoreEnd
