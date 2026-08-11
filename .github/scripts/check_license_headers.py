#!/usr/bin/env python3
"""Check that all source files have the correct Microsoft license header.

This script verifies that Python source files contain the appropriate
copyright and license information.
"""

# --------------------------------------------------------------------------------------------
# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT License. See LICENSE in the project root for license information.
# --------------------------------------------------------------------------------------------

import argparse
import re
import sys
from pathlib import Path

PYTHON_LICENSE_PATTERNS = [
    # Hash-style comment header directly after module docstring or at top
    re.compile(
        r"\A(?:#!/usr/bin/env python3?\s*\n)?"  # Optional shebang
        r'(?:""".*?"""\s*\n)?'  # Optional module docstring
        r"(?:[ \t]*\n)*"  # Optional blank lines only (no code allowed)
        r"#\s*-+\s*\n"  # Opening dashes
        r"#\s*Copyright \(c\) Microsoft\s+Corporation\.\s+All\s+rights\s+reserved\.\s*\n"
        r"#\s*Licensed\s+under\s+the\s+MIT\s+License\.\s+See\s+LICENSE\s+in\s+the\s+project\s+root\s+for\s+license\s+information\.\s*\n"
        r"#\s*-+",  # Closing dashes
        re.MULTILINE | re.DOTALL,
    ),
]

PYTHON_EXTENSIONS = {".py"}


def check_file_header(filepath: Path) -> tuple[bool, str]:
    """Check if a file has the correct license header.

    Args:
        filepath: Path to the file to check

    Returns:
        Tuple of (success: bool, message: str)
    """
    try:
        content = filepath.read_text(encoding="utf-8")
    except Exception as e:
        return False, f"Error reading file: {e}"

    # Skip empty files
    if not content.strip():
        return True, "Empty file, skipping"

    suffix = filepath.suffix.lower()
    if suffix not in PYTHON_EXTENSIONS:
        return True, f"Unknown file type '{suffix}', skipping"

    # Check if any pattern matches in the beginning of the file
    header_section = content[:50000]

    for pattern in PYTHON_LICENSE_PATTERNS:
        if pattern.search(header_section):
            return True, "Valid Python license header found"

    return False, "Missing or incorrect Python license header"


def main(argv: list[str] | None = None) -> int:
    """Main entry point for the license header checker.

    Args:
        argv: Command line arguments

    Returns:
        Exit code (0 for success, 1 for failure)
    """
    parser = argparse.ArgumentParser(
        description="Check source files for correct license headers"
    )
    parser.add_argument(
        "files",
        nargs="*",
        type=Path,
        help="Files to check (if not provided, checks all source files)",
    )
    parser.add_argument(
        "--fix",
        action="store_true",
        help="Attempt to add missing license headers (not implemented)",
    )

    args = parser.parse_args(argv)

    if not args.files:
        print("No files provided to check")
        return 0

    failed_files: list[Path] = []
    checked_count = 0

    for filepath in args.files:
        # Skip files in build/ directories
        if "/build/" in str(filepath):
            continue

        # Only check Python files
        suffix = filepath.suffix.lower()
        if suffix not in PYTHON_EXTENSIONS:
            continue

        checked_count += 1
        success, _message = check_file_header(filepath)

        if not success:
            failed_files.append(filepath)

    if failed_files:
        print(f"\n{len(failed_files)} file(s) missing license headers:\n")

        print(f"Python files ({len(failed_files)}):")
        for filepath in failed_files:
            print(f"  - {filepath}")
        print("\nExpected header format:")
        print('  """Module description."""')
        print(
            "  # --------------------------------------------------------------------------------------------"
        )
        print("  # Copyright (c) Microsoft Corporation. All rights reserved.")
        print(
            "  # Licensed under the MIT License. See LICENSE in the project root for license information."
        )
        print(
            "  # --------------------------------------------------------------------------------------------"
        )
        print(
            "\nNote: License header must come directly after the module docstring (or at top if no docstring)."
        )

        return 1

    print(f"✅ All {checked_count} checked file(s) have valid license headers")
    return 0


if __name__ == "__main__":
    sys.exit(main())
