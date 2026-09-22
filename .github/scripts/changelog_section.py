# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Print one version's CHANGELOG section, for use as GitHub release notes.

The notes for a release are already written by the time it is tagged -- in
``CHANGELOG.md``. Extracting them here keeps the release page and the
changelog from drifting apart, and means a release needs no separate
hand-written body.

Configuration comes from the environment, since this runs in CI:

* ``RELEASE_TAG`` -- the tag being released, e.g. ``v0.6.0``. Defaults to
  ``GITHUB_REF_NAME``.

Run it locally the same way::

    RELEASE_TAG=v0.6.0 python .github/scripts/changelog_section.py
"""

import os
import pathlib
import re
import sys

# --- configuration -------------------------------------------------------
CHANGELOG = pathlib.Path("CHANGELOG.md")
REPO_URL = "https://github.com/RWTH-IAEW/opf-potpourri"
# -------------------------------------------------------------------------


def section(version: str, text: str) -> str | None:
    """The body under ``## [version] ...``, up to the next ``## `` heading."""
    lines = text.splitlines()
    start = None
    for i, line in enumerate(lines):
        if re.match(rf"^## \[{re.escape(version)}\]", line):
            start = i + 1
            break
    if start is None:
        return None
    end = len(lines)
    for i in range(start, len(lines)):
        if lines[i].startswith("## "):
            end = i
            break
    return "\n".join(lines[start:end]).strip("\n")


def main() -> int:
    """Print the tagged version's notes, ready for the release body."""
    tag = os.environ.get("RELEASE_TAG") or os.environ.get(
        "GITHUB_REF_NAME", ""
    )
    if not tag:
        print(
            "RELEASE_TAG and GITHUB_REF_NAME are both unset", file=sys.stderr
        )
        return 2
    version = tag[1:] if tag.startswith("v") else tag

    body = section(version, CHANGELOG.read_text(encoding="utf-8"))
    if body is None:
        print(f"no '## [{version}]' section in {CHANGELOG}", file=sys.stderr)
        return 1

    print(body)
    print()
    print(f"**Full changelog:** {REPO_URL}/blob/{tag}/CHANGELOG.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
