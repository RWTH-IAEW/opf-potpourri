# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Write first-party licensing headers.  Preview by default.

    python tools/fix_license_headers.py

prints what it *would* change and writes nothing.  To apply, set
``APPLY = True`` below and run it again.  This is never wired into CI:
the gate is read-only, and a machine must not repair a licence notice
on its own.

It only ever writes the first-party policy, and refuses any file whose
header already declares a different licence or an unverifiable
copyright holder.  Such a file is a provenance question -- resolve it by
hand against docs/licensing.md.
"""

from __future__ import annotations

import os
import sys

# ----------------------------- configuration -----------------------------
# False previews the change set; True writes it to disk.
APPLY = False
# ------------------------------------------------------------------------

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from license_headers import (  # noqa: E402
    REPO_ROOT,
    DiscoveryError,
    fix_repository,
)


def main(repo_root: str = REPO_ROOT, apply: bool = APPLY) -> int:
    """Write, or preview, the first-party headers.

    Args:
        repo_root: Repository to act on. Defaults to this checkout.
        apply: False previews the change set without writing.

    Returns:
        Process exit status: 1 if any file was refused, else 0.
    """
    try:
        notes = fix_repository(repo_root, dry_run=not apply)
    except DiscoveryError as exc:
        print(f"license-headers: DISCOVERY FAILED: {exc}", file=sys.stderr)
        return 2

    refused = {p: n for p, n in notes.items() if n.startswith("refused")}
    changed = {p: n for p, n in notes.items() if not n.startswith("refused")}

    verb = "wrote" if apply else "would write"
    for path in sorted(changed):
        print(f"{verb} {changed[path]} header: {path}")
    for path in sorted(refused):
        print(f"REFUSED {path}: {refused[path]}", file=sys.stderr)

    if not notes:
        print("license-headers: nothing to do, all headers already match.")
    elif not apply:
        print(
            f"\n{len(changed)} file(s) would change. Nothing was written. "
            "Set APPLY = True in this file to apply."
        )
    else:
        print(f"\n{len(changed)} file(s) written.")
    return 1 if refused else 0


if __name__ == "__main__":
    raise SystemExit(main())
