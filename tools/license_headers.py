# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT
#
# REUSE-IgnoreStart -- below this line SPDX tags appear as data
# (policy constants, test fixtures, help text), not as declarations
# for this file. See docs/licensing.md.

"""Licensing-header policy for potpourri: discovery, checking, fixing.

Single source of truth for the repository's Python licensing headers.
Consumed by :mod:`tools.check_license_headers` (read-only gate, run by
pre-commit and CI) and :mod:`tools.fix_license_headers` (writer, which
must be invoked deliberately and previews by default).

The policy is deliberately narrow.  It enforces the *mechanical* shape
of the header -- that a machine-readable notice physically exists in the
Python comment header -- and the first-party holder/licence wording that
``LICENSE`` already declares.  It does **not** decide who owns a file.
Provenance is an audit question; see ``docs/licensing.md``.

Two properties matter and are covered by tests:

* An ``SPDX-`` string inside a docstring, a string literal or a
  documentation example does **not** satisfy the check.  Only the
  physical leading comment block counts.
* Applying the fixer twice changes nothing the second time.
"""

from __future__ import annotations

import ast
import os
import re
import subprocess
from dataclasses import dataclass, field

# --------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Canonical first-party notice.  The holder wording is copied verbatim
# from LICENSE; do not abbreviate it.  The year range spans the
# repository's own history (earliest revision 2023-05-18) rather than
# the single year in LICENSE or the current year.
FIRST_PARTY_COPYRIGHT = (
    "2023-2026 Institute for High Voltage Equipment and Grids, "
    "Digitalization and Energy Economics (IAEW), RWTH Aachen University"
)

# SPDX expression for first-party files, matching LICENSE and the
# ``license`` field in pyproject.toml.
FIRST_PARTY_LICENSE = "MIT"

# Files whose licensing differs from the first-party default, as
# repo-relative path -> (SPDX expression, evidence).  Empty: the 2026-09
# audit found no third-party or mixed-origin Python source here.  Adding
# an entry is a statement about provenance and needs its evidence
# recorded in docs/licensing.md, not just a path.
LICENSE_EXCEPTIONS: dict[str, tuple[str, str]] = {}

# Paths excluded from the header requirement, as repo-relative path ->
# reason.  Kept empty on purpose: a blanket exclusion is how unresolved
# licensing questions get hidden behind a green check.
HEADER_EXEMPTIONS: dict[str, str] = {}

COPYRIGHT_TAG = "SPDX-FileCopyrightText:"
LICENSE_TAG = "SPDX-License-Identifier:"

# PEP 263 encoding declaration, valid only on the first or second line.
_CODING_RE = re.compile(r"^[ \t\f]*#.*?coding[:=][ \t]*([-_.a-zA-Z0-9]+)")

# An ``(c) YEAR, Holder`` notice as written in this repository's existing
# script docstrings.  Mirrored into SPDX metadata rather than invented.
_LEGACY_NOTICE_RE = re.compile(r"^\(c\)\s*(\d{4})\s*,\s*(.+?)\s*$")

BOM = "﻿"


class DiscoveryError(RuntimeError):
    """Raised when file discovery fails or yields nothing.

    An empty scan must never be reported as a pass.
    """


# --------------------------------------------------------------------
# Discovery
# --------------------------------------------------------------------


def _git(args: list[str], repo_root: str) -> str:
    try:
        proc = subprocess.run(
            ["git", *args], cwd=repo_root, capture_output=True
        )
    except OSError as exc:
        # Unusable working directory or no git binary: a tool failure,
        # which must not be mistaken for "nothing to check".
        raise DiscoveryError(f"cannot run git in {repo_root}: {exc}") from exc
    if proc.returncode != 0:
        raise DiscoveryError(
            "git {} failed in {}: {}".format(
                " ".join(args), repo_root, proc.stderr.decode().strip()
            )
        )
    return proc.stdout.decode()


def python_files(repo_root: str = REPO_ROOT) -> list[str]:
    """Repo-relative Python files in scope, sorted.

    Covers tracked files (including newly staged ones) and untracked
    files git does not ignore, so a new module is caught before it is
    recorded.  Files inside submodules are gitlinks to git and are
    therefore never returned.
    """
    patterns = ["*.py", "*.pyi", "*.pyw"]
    tracked = _git(["ls-files", "-z", "--", *patterns], repo_root)
    fresh = _git(
        ["ls-files", "-z", "--others", "--exclude-standard", "--", *patterns],
        repo_root,
    )
    found = {p for p in (tracked + fresh).split("\0") if p}
    if not found:
        raise DiscoveryError(
            f"no Python files discovered under {repo_root}; refusing to "
            "report an empty scan as success"
        )
    return sorted(found)


# --------------------------------------------------------------------
# Header parsing
# --------------------------------------------------------------------


@dataclass
class Header:
    """The physical comment header at the top of a Python file."""

    #: Number of leading lines that must stay first: BOM-bearing
    #: shebang, shebang, PEP 263 encoding declaration.
    prelude: int = 0
    #: Index one past the last line of the leading comment/blank block.
    end: int = 0
    #: Text after ``SPDX-FileCopyrightText:`` for each header line.
    copyrights: list[str] = field(default_factory=list)
    #: Text after ``SPDX-License-Identifier:`` for each header line.
    licenses: list[str] = field(default_factory=list)


def _is_comment(line: str) -> bool:
    return line.lstrip().startswith("#")


def split_bom(text: str) -> tuple[str, str]:
    """Separate a leading byte-order mark from the rest of ``text``.

    The BOM must stay the very first thing in the file, ahead of even a
    shebang, so it is peeled off before any line arithmetic and put back
    afterwards.
    """
    if text.startswith(BOM):
        return BOM, text[len(BOM) :]
    return "", text


def _tag_value(comment_body: str, tag: str) -> str:
    return comment_body[len(tag) :].strip()


def parse_header(text: str) -> Header:
    """Parse the leading comment block of ``text``.

    The header ends at the first line that is neither a comment nor
    blank, so a module docstring -- and anything inside it -- is outside
    the header by construction.
    """
    _, text = split_bom(text)
    lines = text.split("\n")
    head = Header()

    i = 0
    if lines and lines[0].startswith("#!"):
        i = 1
    # PEP 263: the encoding declaration is honoured only on line 1 or 2.
    if i < len(lines) and i < 2 and _CODING_RE.match(lines[i]):
        i += 1
    head.prelude = i

    j = i
    while j < len(lines) and (lines[j].strip() == "" or _is_comment(lines[j])):
        j += 1
    # Trailing blank lines separate the header from the body.
    while j > i and lines[j - 1].strip() == "":
        j -= 1
    head.end = j

    for line in lines[i:j]:
        stripped = line.lstrip()
        if not stripped.startswith("#"):
            continue
        body = stripped.lstrip("#").strip()
        if body.startswith(COPYRIGHT_TAG):
            head.copyrights.append(_tag_value(body, COPYRIGHT_TAG))
        elif body.startswith(LICENSE_TAG):
            head.licenses.append(_tag_value(body, LICENSE_TAG))
    return head


def canonical_expression(expression: str) -> str:
    """Canonical SPDX form of ``expression``.

    Validated with :mod:`packaging.licenses` (PEP 639), which ships with
    every modern setuptools/pytest environment, so the check needs no
    extra dependency and no solver.
    """
    from packaging.licenses import canonicalize_license_expression

    return str(canonicalize_license_expression(expression))


def docstring_notices(text: str) -> list[str]:
    """``(c) YEAR, Holder`` notices already present in the docstring.

    These are pre-existing claims in this repository.  The fixer mirrors
    them into SPDX metadata so a compact header does not silently drop a
    notice; it never invents one.
    """
    try:
        doc = ast.get_docstring(ast.parse(split_bom(text)[1]))
    except SyntaxError:
        return []
    if not doc:
        return []
    out = []
    for line in doc.split("\n"):
        match = _LEGACY_NOTICE_RE.match(line.strip())
        if match:
            out.append(f"{match.group(1)} {match.group(2)}")
    return out


# --------------------------------------------------------------------
# Checking
# --------------------------------------------------------------------


def expected_license(path: str) -> str:
    return LICENSE_EXCEPTIONS.get(path, (FIRST_PARTY_LICENSE, ""))[0]


def check_text(path: str, text: str) -> list[str]:
    """Policy problems in ``text``; empty means compliant."""
    if path in HEADER_EXEMPTIONS:
        return []
    problems: list[str] = []
    head = parse_header(text)
    wanted = expected_license(path)

    if not head.copyrights:
        problems.append(
            f"no '# {COPYRIGHT_TAG} ...' line in the Python comment header"
            + (
                " (found one only inside a string/docstring, which does"
                " not count)"
                if COPYRIGHT_TAG in text
                else ""
            )
        )
    elif any(not c for c in head.copyrights):
        problems.append(f"empty '# {COPYRIGHT_TAG}' line")

    if not head.licenses:
        problems.append(
            f"no '# {LICENSE_TAG} ...' line in the Python comment header"
            + (
                " (found one only inside a string/docstring, which does"
                " not count)"
                if LICENSE_TAG in text
                else ""
            )
        )
    elif len(set(head.licenses)) > 1:
        problems.append(
            "contradictory licence declarations in the header: "
            + ", ".join(sorted(set(head.licenses)))
        )
    elif len(head.licenses) > 1:
        problems.append(
            f"duplicated '# {LICENSE_TAG} {head.licenses[0]}' line "
            f"({len(head.licenses)} occurrences)"
        )
    else:
        declared = head.licenses[0]
        try:
            canonical = canonical_expression(declared)
        except Exception as exc:
            problems.append(f"invalid SPDX expression {declared!r}: {exc}")
        else:
            if canonical != canonical_expression(wanted):
                problems.append(
                    f"declares {canonical!r} but policy for this file is "
                    f"{wanted!r}"
                )
            elif canonical != declared:
                problems.append(
                    f"non-canonical SPDX expression {declared!r}; write "
                    f"{canonical!r}"
                )

    if head.copyrights and path not in LICENSE_EXCEPTIONS:
        if FIRST_PARTY_COPYRIGHT not in head.copyrights:
            problems.append(
                "first-party file is missing the canonical notice "
                f"'# {COPYRIGHT_TAG} {FIRST_PARTY_COPYRIGHT}'"
            )
    if len(head.copyrights) != len(set(head.copyrights)):
        problems.append("duplicated SPDX-FileCopyrightText line")

    for notice in docstring_notices(text):
        if notice not in head.copyrights:
            problems.append(
                f"docstring claims '(c) {notice}' but no matching "
                f"'# {COPYRIGHT_TAG} {notice}' line is in the header"
            )
    return problems


def read_text(full_path: str) -> str:
    with open(full_path, encoding="utf-8") as handle:
        return handle.read()


def check_repository(repo_root: str = REPO_ROOT) -> dict[str, list[str]]:
    """Map of repo-relative path -> problems, for every in-scope file."""
    failures = {}
    for rel in python_files(repo_root):
        problems = check_text(rel, read_text(os.path.join(repo_root, rel)))
        if problems:
            failures[rel] = problems
    return failures


# --------------------------------------------------------------------
# Fixing
# --------------------------------------------------------------------


def render_header(copyrights: list[str], license_id: str) -> list[str]:
    lines = [f"# {COPYRIGHT_TAG} {c}" for c in copyrights]
    lines.append("#")
    lines.append(f"# {LICENSE_TAG} {license_id}")
    return lines


def fix_text(path: str, text: str) -> tuple[str, str]:
    """Return ``(new_text, note)``.

    ``note`` is empty when nothing changed.  Refuses any file whose
    licensing the policy cannot state on the evidence, so ambiguity is
    never resolved by the writer.
    """
    if path in HEADER_EXEMPTIONS:
        return text, ""
    if path in LICENSE_EXCEPTIONS:
        return text, (
            "refused: third-party/mixed-origin licensing is recorded by "
            "hand, not written by the fixer"
        )

    bom, text = split_bom(text)
    head = parse_header(text)
    if head.licenses and set(head.licenses) != {FIRST_PARTY_LICENSE}:
        return bom + text, (
            "refused: header already declares "
            f"{sorted(set(head.licenses))}, which the first-party policy "
            "cannot overwrite"
        )
    foreign = [
        c
        for c in head.copyrights
        if c != FIRST_PARTY_COPYRIGHT and c not in docstring_notices(text)
    ]
    if foreign:
        return bom + text, (
            f"refused: header carries copyright notices {foreign} that the "
            "policy cannot verify or rewrite"
        )

    wanted = [FIRST_PARTY_COPYRIGHT]
    for notice in docstring_notices(text):
        if notice not in wanted:
            wanted.append(notice)
    for existing in head.copyrights:
        if existing not in wanted:
            wanted.append(existing)

    lines = text.split("\n")
    body = lines[head.end :]
    # Drop the SPDX lines we are replacing, keeping every other comment
    # (authorship notes, lint directives, tool pragmas) exactly as it is.
    kept = []
    for line in lines[head.prelude : head.end]:
        stripped = line.lstrip("﻿").lstrip().lstrip("#").strip()
        if stripped.startswith(COPYRIGHT_TAG) or stripped.startswith(
            LICENSE_TAG
        ):
            continue
        kept.append(line)
    while kept and kept[0].strip() in ("", "#"):
        kept.pop(0)
    while kept and kept[-1].strip() in ("", "#"):
        kept.pop()

    header = render_header(wanted, FIRST_PARTY_LICENSE)
    if kept:
        header = header + ["#"] + kept
    while body and body[0].strip() == "":
        body.pop(0)

    rebuilt = lines[: head.prelude] + header
    if body and not (len(body) == 1 and body[0] == ""):
        rebuilt += [""]
    rebuilt += body
    out = bom + "\n".join(rebuilt)
    # A previously empty file has no body to carry the final newline.
    if out and not out.endswith("\n"):
        out += "\n"
    # Idempotence by construction: the writer reports a change only when
    # the rebuilt text actually differs from what is on disk.
    if out == bom + text:
        return out, ""
    action = "added" if not head.licenses else "normalised"
    return out, action


def fix_repository(
    repo_root: str = REPO_ROOT, dry_run: bool = True
) -> dict[str, str]:
    """Apply :func:`fix_text` to every in-scope file.

    With ``dry_run`` (the default) nothing is written.
    """
    notes = {}
    for rel in python_files(repo_root):
        full = os.path.join(repo_root, rel)
        original = read_text(full)
        updated, note = fix_text(rel, original)
        if not note:
            continue
        notes[rel] = note
        if not dry_run and updated != original:
            with open(full, "w", encoding="utf-8", newline="") as handle:
                handle.write(updated)
    return notes


# REUSE-IgnoreEnd
