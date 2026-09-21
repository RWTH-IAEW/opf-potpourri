# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT
#
# REUSE-IgnoreStart -- below this line SPDX tags appear as data
# (policy constants, test fixtures, help text), not as declarations
# for this file. See docs/licensing.md.

"""Read-only licensing-header gate for every in-scope Python file.

Run by pre-commit and by CI.  It never writes: use
``tools/fix_license_headers.py`` for that.

    python tools/check_license_headers.py

Exit codes: 0 compliant, 1 policy violations, 2 discovery/tool error.
A discovery error is deliberately *not* a pass -- an empty or broken
scan must never look like success.
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from license_headers import (  # noqa: E402
    FIRST_PARTY_COPYRIGHT,
    FIRST_PARTY_LICENSE,
    REPO_ROOT,
    DiscoveryError,
    check_repository,
    python_files,
)


def main(repo_root: str = REPO_ROOT) -> int:
    """Check every in-scope Python file and report the failures.

    Args:
        repo_root: Repository to scan. Defaults to this checkout.

    Returns:
        Process exit status: 0 compliant, 1 policy violations, 2 a discovery or
            tool error. A broken scan is never reported as a pass.
    """
    try:
        scanned = python_files(repo_root)
        failures = check_repository(repo_root)
    except DiscoveryError as exc:
        print(f"license-headers: DISCOVERY FAILED: {exc}", file=sys.stderr)
        return 2
    except Exception as exc:  # tooling fault, not a clean bill of health
        print(f"license-headers: CHECK ERROR: {exc!r}", file=sys.stderr)
        return 2

    if not failures:
        print(f"license-headers: {len(scanned)} Python files OK")
        return 0

    for path in sorted(failures):
        for problem in failures[path]:
            print(f"{path}: {problem}")
    print(
        f"\nlicense-headers: {len(failures)} of {len(scanned)} Python files "
        "have a licensing-header problem.",
        file=sys.stderr,
    )
    print(
        "\nFirst-party files need this at the very top of the file, as real\n"
        "comments (a notice inside a docstring does not count):\n\n"
        f"    # SPDX-FileCopyrightText: {FIRST_PARTY_COPYRIGHT}\n"
        "    #\n"
        f"    # SPDX-License-Identifier: {FIRST_PARTY_LICENSE}\n\n"
        "Preview the fix:  python tools/fix_license_headers.py\n"
        "Apply it:         set APPLY = True in that file, then re-run.\n"
        "Copied or adapted third-party code must NOT be fixed this way;\n"
        "see docs/licensing.md.",
        file=sys.stderr,
    )
    return 1


if __name__ == "__main__":
    raise SystemExit(main())

# REUSE-IgnoreEnd
