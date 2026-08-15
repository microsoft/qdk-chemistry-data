#!/usr/bin/env python3
"""Check that the committed JSON matches what regenerate.py --collate-only produces.

This ensures the JSON stays in sync with the raw_output/ logs and the
collation logic in regenerate.py. Exits 0 if they match, 1 if they differ.
"""

# --------------------------------------------------------------------------------------------
# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT License. See LICENSE in the project root for license information.
# --------------------------------------------------------------------------------------------

import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

# Keys that regenerate.py adds during collation but are not stored in the
# committed JSON (they are intermediate artefacts).
TRANSIENT_KEYS = {"output_file_text", "output_file_name", "workflow_timing"}


def normalize(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Remove transient keys and sort by name for stable comparison."""
    cleaned = []
    for mol in records:
        cleaned.append({k: v for k, v in mol.items() if k not in TRANSIENT_KEYS})
    return sorted(cleaned, key=lambda m: m["name"])


def project_to_schema(
    record: dict[str, Any], reference: dict[str, Any]
) -> dict[str, Any]:
    """Project a record down to only the keys (and sub-keys) in the reference.

    The regenerate script may parse additional fields from the logs that were
    intentionally excluded from the committed JSON. This function strips those
    extra fields so the comparison is fair.
    """
    result: dict[str, Any] = {}
    for k, ref_val in reference.items():
        if k not in record:
            continue
        rec_val = record[k]
        if isinstance(ref_val, dict) and isinstance(rec_val, dict):
            result[k] = project_to_schema(rec_val, ref_val)
        elif isinstance(ref_val, list) and isinstance(rec_val, list):
            # For lists of dicts, project each element
            if (
                ref_val
                and isinstance(ref_val[0], dict)
                and rec_val
                and isinstance(rec_val[0], dict)
            ):
                result[k] = [
                    project_to_schema(r, ref_val[min(i, len(ref_val) - 1)])
                    for i, r in enumerate(rec_val)
                ]
            else:
                result[k] = rec_val
        else:
            result[k] = rec_val
    return result


def diff_records(
    committed: list[dict[str, Any]], regenerated: list[dict[str, Any]]
) -> list[str]:
    """Return human-readable differences between two molecule lists."""
    errors: list[str] = []

    committed_by_name = {m["name"]: m for m in committed}
    regenerated_by_name = {m["name"]: m for m in regenerated}

    only_committed = set(committed_by_name) - set(regenerated_by_name)
    only_regenerated = set(regenerated_by_name) - set(committed_by_name)

    for name in sorted(only_committed):
        errors.append(
            f"{name}: present in committed JSON but not in regenerated output"
        )
    for name in sorted(only_regenerated):
        errors.append(
            f"{name}: present in regenerated output but not in committed JSON"
        )

    for name in sorted(set(committed_by_name) & set(regenerated_by_name)):
        c = committed_by_name[name]
        r = project_to_schema(regenerated_by_name[name], c)
        if c != r:
            # Find which keys differ
            all_keys = sorted(set(c.keys()) | set(r.keys()))
            for key in all_keys:
                if c.get(key) != r.get(key):
                    errors.append(f"{name}: field '{key}' differs")

    return errors


def main() -> int:
    """Scan all datasets and verify each JSON matches its regenerated output."""
    repo_root = Path(__file__).resolve().parents[2]
    datasets_dir = repo_root / "data" / "molecules"

    all_errors: list[str] = []

    for dataset_dir in sorted(datasets_dir.iterdir()):
        if not dataset_dir.is_dir():
            continue

        regenerate_script = dataset_dir / "regenerate.py"
        json_files = list(dataset_dir.glob("*.json"))

        if not regenerate_script.exists() or not json_files:
            continue

        for json_file in json_files:
            print(f"Checking {dataset_dir.name}/{json_file.name}...")

            # Load committed JSON
            with open(json_file, encoding="utf-8") as f:
                committed = json.load(f)

            # Run regenerate.py --collate-only in a disposable dataset copy.
            with tempfile.TemporaryDirectory() as temp_dir:
                working_dir = Path(temp_dir) / dataset_dir.name
                shutil.copytree(dataset_dir, working_dir)
                working_script = working_dir / regenerate_script.name
                working_json = working_dir / json_file.name

                result = subprocess.run(
                    [sys.executable, str(working_script), "--collate-only"],
                    cwd=str(working_dir),
                    capture_output=True,
                    text=True,
                )

                if result.returncode != 0:
                    all_errors.append(
                        f"{dataset_dir.name}: regenerate.py --collate-only failed:\n"
                        f"  stdout: {result.stdout[-500:]}\n"
                        f"  stderr: {result.stderr[-500:]}"
                    )
                    continue

                with open(working_json, encoding="utf-8") as f:
                    regenerated = json.load(f)

            # Normalize both for comparison
            committed_norm = normalize(committed)
            regenerated_norm = normalize(regenerated)

            if committed_norm != regenerated_norm:
                diffs = diff_records(committed_norm, regenerated_norm)
                for d in diffs:
                    all_errors.append(f"{dataset_dir.name}/{json_file.name}: {d}")

    if all_errors:
        print(f"\nFAILED: {len(all_errors)} difference(s) found:")
        for err in all_errors:
            print(f"  ✗ {err}")
        print("\nRun 'regenerate.py --collate-only' and commit the updated JSON.")
        return 1

    print(
        "\nAll regeneration checks passed — committed JSON matches regenerated output."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
