# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Regenerate the per-file documentation inventory.

    python tools/docs_inventory.py

Writes `docs/documentation-inventory.csv`: one row per tracked Python
file with its role, how much of it is documented, and what still fails
a check. Complements `interrogate`, which reports one number for the
package and a per-file percentage but not the role or the reason.

Coverage here is counted the same way `interrogate` counts it, so the
`documented`/`documentable` columns sum to the figure the gate prints.
"""

from __future__ import annotations

import ast
import csv
import os
import subprocess

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTPUT = os.path.join(REPO, "docs", "documentation-inventory.csv")

FIELDS = (
    "path",
    "role",
    "loc",
    "documentable",
    "documented",
    "coverage_pct",
    "module_docstring",
    "ruff_findings",
)

# Longest prefix wins, so the research entries must precede the general
# src/potpourri one.
ROLES = (
    ("src/potpourri/research/", "research (exploratory study code)"),
    ("src/potpourri/models_multi_period/", "library: multi-period model"),
    ("src/potpourri/models/", "library: single-period model"),
    ("src/potpourri/technologies/", "library: flexible-device mix-in"),
    ("src/potpourri/benchmarks/", "library: benchmark loader"),
    ("src/potpourri/net_augmentation/", "library: network preprocessing"),
    ("src/potpourri/plotting/", "library: plotting"),
    ("src/potpourri/", "library: package root"),
    ("scripts/research/", "runnable study script"),
    ("scripts/", "runnable example script"),
    ("tests/installation_with_pip/", "test: installed package"),
    ("tests/", "test: unit and regression"),
    ("tools/", "development tooling"),
)


def role_of(path: str) -> str:
    """Classify a file by its directory.

    Args:
        path: Repo-relative path.

    Returns:
        A short human-readable role, or `"other"` when no prefix
        matches.
    """
    for prefix, role in ROLES:
        if path.startswith(prefix):
            return role
    return "other"


def tracked_python() -> list[str]:
    """Repo-relative Python files in scope, sorted.

    Tracked files plus untracked ones git does not ignore, so a newly
    added module appears in the inventory before it is recorded.

    Returns:
        The paths `git ls-files` reports for Python sources.

    Raises:
        RuntimeError: If git fails, so an empty inventory is never
            mistaken for a clean one.
    """
    patterns = ["*.py", "*.pyi", "*.pyw"]
    seen = set()
    for extra in ([], ["--others", "--exclude-standard"]):
        proc = subprocess.run(
            ["git", "ls-files", "-z", *extra, "--", *patterns],
            cwd=REPO,
            capture_output=True,
        )
        if proc.returncode != 0:
            raise RuntimeError(proc.stderr.decode().strip())
        seen.update(p for p in proc.stdout.decode().split("\0") if p)
    paths = sorted(seen)
    if not paths:
        raise RuntimeError("no Python files found; refusing to write")
    return paths


def counts(path: str) -> tuple[int, int, bool]:
    """Count documentable and documented definitions in one file.

    Counts the module itself plus every class and function, nested ones
    included, and skips `__init__` -- matching `[tool.interrogate]`.

    Args:
        path: Absolute path to the file.

    Returns:
        A `(documentable, documented, has_module_docstring)` triple.
    """
    tree = ast.parse(open(path, encoding="utf-8").read())
    has_module_doc = ast.get_docstring(tree) is not None
    total, done = 1, int(has_module_doc)
    for node in ast.walk(tree):
        if isinstance(
            node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)
        ):
            if node.name == "__init__":
                continue
            total += 1
            done += ast.get_docstring(node) is not None
    return total, done, has_module_doc


def ruff_findings() -> dict[str, int]:
    """Count the ruff findings still open per file.

    Returns:
        Mapping of repo-relative path to finding count. Files with none
        are absent.
    """
    proc = subprocess.run(
        ["ruff", "check", "--output-format", "concise", "."],
        cwd=REPO,
        capture_output=True,
        text=True,
    )
    out: dict[str, int] = {}
    for line in proc.stdout.splitlines():
        if ":" not in line:
            continue
        path = line.split(":", 1)[0]
        if path.endswith(".py"):
            out[path] = out.get(path, 0) + 1
    return out


def write(output: str = OUTPUT) -> int:
    """Write the inventory CSV and print a per-role summary.

    Args:
        output: Destination CSV path.

    Returns:
        0 always; this is a reporting tool, not a gate.
    """
    findings = ruff_findings()
    rows = []
    for rel in tracked_python():
        full = os.path.join(REPO, rel)
        total, done, has_doc = counts(full)
        rows.append(
            {
                "path": rel,
                "role": role_of(rel),
                "loc": sum(1 for _ in open(full, encoding="utf-8")),
                "documentable": total,
                "documented": done,
                "coverage_pct": round(100.0 * done / total, 1),
                "module_docstring": "yes" if has_doc else "no",
                "ruff_findings": findings.get(rel, 0),
            }
        )
    os.makedirs(os.path.dirname(output), exist_ok=True)
    with open(output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    print(f"{len(rows)} files -> {output}\n")
    print(f"{'role':38} {'files':>5} {'defs':>6} {'documented':>11} {'%':>7}")
    for role in sorted({r["role"] for r in rows}):
        group = [r for r in rows if r["role"] == role]
        total = sum(r["documentable"] for r in group)
        done = sum(r["documented"] for r in group)
        print(
            f"{role:38} {len(group):5} {total:6} {done:11} "
            f"{100.0 * done / total:6.1f}%"
        )
    total = sum(r["documentable"] for r in rows)
    done = sum(r["documented"] for r in rows)
    print(
        f"{'ALL TRACKED PYTHON':38} {len(rows):5} {total:6} {done:11} "
        f"{100.0 * done / total:6.1f}%"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(write())
