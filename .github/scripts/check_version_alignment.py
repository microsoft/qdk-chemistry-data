#!/usr/bin/env python3
"""Check that all version strings in the codebase are aligned.

This script verifies that version strings across multiple files are consistent:
- __version__ in qdk_chemistry_data/__init__.py
- version in pyproject.toml

Exit codes:
    0: All versions are aligned
    1: Versions are misaligned (with details printed)
"""

# --------------------------------------------------------------------------------------------
# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT License. See LICENSE in the project root for license information.
# --------------------------------------------------------------------------------------------

import re
import sys
from pathlib import Path
from typing import NamedTuple, TypedDict


class VersionLocation(NamedTuple):
    """Location of a version string in the codebase."""

    file_path: Path
    pattern: str
    description: str


class VersionInfo(TypedDict):
    """Information about a version found in the codebase."""

    version: str
    location: VersionLocation


def extract_version(file_path: Path, pattern: str) -> str | None:
    """Extract version string from a file using a regex pattern.

    Args:
        file_path: Path to the file to search
        pattern: Regex pattern with a capture group for the version

    Returns:
        The extracted version string, or None if not found
    """
    if not file_path.exists():
        return None

    try:
        content = file_path.read_text()
        match = re.search(pattern, content)
        if match:
            return match.group(1)
    except (OSError, UnicodeDecodeError) as e:
        print(f"Error reading {file_path}: {e}", file=sys.stderr)

    return None


def check_versions() -> int:
    """Check that all version strings are aligned.

    Returns:
        0 if all versions match, 1 if there are mismatches
    """
    # Define the root directory (repository root)
    repo_root = Path(__file__).parent.parent.parent

    # Define version locations
    version_locations = [
        VersionLocation(
            file_path=repo_root / "qdk_chemistry_data/__init__.py",
            pattern=r'__version__\s*=\s*["\']([^"\']+)["\']',
            description="Python __version__",
        ),
        VersionLocation(
            file_path=repo_root / "pyproject.toml",
            pattern=r'^version\s*=\s*["\']([^"\']+)["\']',
            description="pyproject.toml version",
        ),
    ]

    # Extract all versions
    versions: dict[str, VersionInfo] = {}
    missing: list[VersionLocation] = []

    for loc in version_locations:
        version = extract_version(loc.file_path, loc.pattern)
        if version is None:
            missing.append(loc)
        else:
            versions[loc.description] = {
                "version": version,
                "location": loc,
            }

    # Report missing versions
    if missing:
        print("❌ Version check failed: Missing version strings", file=sys.stderr)
        print(file=sys.stderr)
        for loc in missing:
            print(f"  ⚠️  Missing: {loc.description}", file=sys.stderr)
            print(f"      File: {loc.file_path}", file=sys.stderr)
        return 1

    # Get the canonical version (from pyproject.toml)
    canonical = versions.get("pyproject.toml version")
    if not canonical:
        print("❌ Version check failed: Canonical version not found", file=sys.stderr)
        return 1

    canonical_version = canonical["version"]

    # Check alignment
    misaligned = []
    for desc, ver_info in versions.items():
        if desc == "pyproject.toml version":
            continue  # Skip the canonical version itself

        actual = ver_info["version"]

        if actual != canonical_version:
            misaligned.append(
                {
                    "description": desc,
                    "expected": canonical_version,
                    "actual": actual,
                    "file": ver_info["location"].file_path,
                }
            )

    # Report results
    if misaligned:
        print(
            "❌ Version check failed: Version strings are not aligned", file=sys.stderr
        )
        print(file=sys.stderr)
        print(f"  ✓  Canonical version: {canonical_version}", file=sys.stderr)
        print(
            f"      Location: {canonical['location'].file_path}",
            file=sys.stderr,
        )
        print(file=sys.stderr)

        for item in misaligned:
            print(f"  ✗  {item['description']}: {item['actual']}", file=sys.stderr)
            print(f"      Expected: {item['expected']}", file=sys.stderr)
            print(f"      File: {item['file']}", file=sys.stderr)
            print(file=sys.stderr)

        print(
            "Please update the version strings to match the canonical version.",
            file=sys.stderr,
        )
        return 1

    # All versions aligned
    print(f"✓ All version strings are aligned: {canonical_version}")
    return 0


if __name__ == "__main__":
    sys.exit(check_versions())
