# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Refuse to publish a tag that disagrees with the version in the tree.

Three files carry the release version and two different jobs read them: the
PyPI job builds from ``pyproject.toml`` while the Zenodo job takes its
metadata from ``CITATION.cff``. Before this check only the first was compared
with the tag, and a mismatch made the workflow *skip* while still reporting
success -- so a stale ``CITATION.cff`` archived 0.6.0 on Zenodo under the
previous version, with nothing failing.

Pushing a tag is a deliberate act. If it disagrees with the tree, that is a
mistake, so this exits non-zero and the release stops.

Configuration comes from the environment, since this runs in CI:

* ``RELEASE_TAG`` -- the tag being released, e.g. ``v0.6.0``. Defaults to
  ``GITHUB_REF_NAME``.

Run it locally the same way::

    RELEASE_TAG=v0.6.0 python .github/scripts/check_release_version.py
"""

import os
import pathlib
import re
import sys
import tomllib

import yaml

# --- configuration -------------------------------------------------------
PYPROJECT = pathlib.Path("pyproject.toml")
CITATION = pathlib.Path("CITATION.cff")
CHANGELOG = pathlib.Path("CHANGELOG.md")
# -------------------------------------------------------------------------


def tag_version(tag: str) -> str:
    """The version a tag names: ``v0.6.0`` -> ``0.6.0``."""
    return tag[1:] if tag.startswith("v") else tag


def pyproject_version() -> str:
    """The version setuptools will build and PyPI will serve."""
    with PYPROJECT.open("rb") as handle:
        return tomllib.load(handle)["project"]["version"]


def citation_version() -> str:
    """The version the Zenodo job stamps on the archive."""
    return str(yaml.safe_load(CITATION.read_text(encoding="utf-8"))["version"])


def changelog_has_section(version: str) -> bool:
    """A released version must have its own CHANGELOG heading."""
    pattern = rf"^## \[{re.escape(version)}\]"
    return any(
        re.match(pattern, line)
        for line in CHANGELOG.read_text(encoding="utf-8").splitlines()
    )


def main() -> int:
    """Compare the tag with every file that carries the version."""
    tag = os.environ.get("RELEASE_TAG") or os.environ.get(
        "GITHUB_REF_NAME", ""
    )
    if not tag:
        print(
            "RELEASE_TAG and GITHUB_REF_NAME are both unset", file=sys.stderr
        )
        return 2

    want = tag_version(tag)
    problems = []

    for name, path, found in (
        ("pyproject.toml", PYPROJECT, pyproject_version()),
        ("CITATION.cff", CITATION, citation_version()),
    ):
        status = "ok" if found == want else "MISMATCH"
        print(f"  {name:<16} {found:<12} {status}")
        if found != want:
            problems.append(f"{path} says {found!r}, tag {tag} names {want!r}")

    if changelog_has_section(want):
        print(f"  {'CHANGELOG.md':<16} [{want}]{'':<4} ok")
    else:
        print(f"  {'CHANGELOG.md':<16} {'-':<12} MISSING SECTION")
        problems.append(f"{CHANGELOG} has no '## [{want}]' section")

    if problems:
        print(file=sys.stderr)
        print(f"Refusing to publish {tag}:", file=sys.stderr)
        for p in problems:
            print(f"  - {p}", file=sys.stderr)
        print(
            "\nBump the version everywhere and move the tag, rather than "
            "releasing a tree that disagrees with the tag.",
            file=sys.stderr,
        )
        return 1

    print(f"\n{tag} agrees with the tree; publishing {want}.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
