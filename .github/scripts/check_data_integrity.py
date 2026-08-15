#!/usr/bin/env python3
"""Validate dataset integrity: JSON structure, cross-references, and consistency.

This script checks that:
- Each dataset JSON file is valid and contains the expected schema
- All molecules referenced in JSON have corresponding XYZ, image, and log files
- XYZ content embedded in JSON matches the standalone XYZ files
- Numeric fields are within physically reasonable ranges
"""

# --------------------------------------------------------------------------------------------
# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT License. See LICENSE in the project root for license information.
# --------------------------------------------------------------------------------------------

import json
import sys
from pathlib import Path

# Required top-level keys for each molecule record
REQUIRED_MOLECULE_KEYS = {
    "name",
    "xyz",
    "structure",
    "scf_energy_hartree",
    "orbitals_summary",
    "casci_energies_hartree",
    "initial_casci_energy_hartree",
    "autocas_energy_hartree",
    "autocas_selected_indices",
    "hamiltonian_summaries",
    "sparse_ci_finder",
}

REQUIRED_STRUCTURE_KEYS = {
    "num_atoms",
    "composition",
    "total_mass_amu",
    "nuclear_repulsion_energy_eh",
}

REQUIRED_SPARSE_CI_KEYS = {
    "n_dets",
    "energy_hartree",
    "delta_e_mhartree",
    "determinants",
}

NON_DATA_JSON_FILENAMES = {"cgmanifest.json"}


def check_dataset(dataset_dir: Path) -> list[str]:
    """Validate a single dataset directory. Returns a list of error messages."""
    errors: list[str] = []
    json_files = [
        path
        for path in dataset_dir.glob("*.json")
        if path.name not in NON_DATA_JSON_FILENAMES
    ]

    if not json_files:
        errors.append(f"{dataset_dir.name}: No JSON data file found")
        return errors

    for json_file in json_files:
        try:
            with open(json_file, encoding="utf-8") as f:
                data = json.load(f)
        except json.JSONDecodeError as e:
            errors.append(f"{json_file.name}: Invalid JSON - {e}")
            continue

        if not isinstance(data, list):
            errors.append(f"{json_file.name}: Expected a JSON array at top level")
            continue

        if len(data) == 0:
            errors.append(f"{json_file.name}: Empty dataset")
            continue

        names_seen: set[str] = set()

        for i, mol in enumerate(data):
            fallback_label = f"index-{i}"
            if not isinstance(mol, dict):
                errors.append(
                    f"{json_file.name}/{fallback_label}: Expected a JSON object"
                )
                continue

            name = mol.get("name")
            if isinstance(name, str) and name:
                label = name
                has_valid_name = True
            else:
                label = fallback_label
                has_valid_name = False
            if "name" in mol and not has_valid_name:
                errors.append(
                    f"{json_file.name}/{label}: Invalid molecule name={name!r}"
                )

            # Check for duplicate names
            if has_valid_name and label in names_seen:
                errors.append(f"{json_file.name}: Duplicate molecule name '{label}'")
            if has_valid_name:
                names_seen.add(label)

            # Check required keys
            missing = REQUIRED_MOLECULE_KEYS - set(mol.keys())
            if missing:
                errors.append(f"{json_file.name}/{label}: Missing keys: {missing}")

            # Check structure sub-keys
            structure = mol.get("structure")
            if "structure" in mol and not isinstance(structure, dict):
                errors.append(
                    f"{json_file.name}/{label}: Expected structure to be a JSON object"
                )
            elif isinstance(structure, dict):
                struct_missing = REQUIRED_STRUCTURE_KEYS - set(structure.keys())
                if struct_missing:
                    errors.append(
                        f"{json_file.name}/{label}: "
                        f"Missing structure keys: {struct_missing}"
                    )

            # Check sparse_ci_finder sub-keys
            sparse_ci = mol.get("sparse_ci_finder")
            if "sparse_ci_finder" in mol and not isinstance(sparse_ci, dict):
                errors.append(
                    f"{json_file.name}/{label}: "
                    "Expected sparse_ci_finder to be a JSON object"
                )
            elif isinstance(sparse_ci, dict):
                sci_missing = REQUIRED_SPARSE_CI_KEYS - set(sparse_ci.keys())
                if sci_missing:
                    errors.append(
                        f"{json_file.name}/{label}: "
                        f"Missing sparse_ci_finder keys: {sci_missing}"
                    )

            # Check XYZ file exists and coordinates match
            xyz_file = dataset_dir / "xyz" / f"{label}.xyz"
            xyz = mol.get("xyz")
            if "xyz" in mol and not isinstance(xyz, str):
                errors.append(f"{json_file.name}/{label}: Expected xyz to be a string")
            elif has_valid_name and not xyz_file.exists():
                errors.append(
                    f"{json_file.name}/{label}: Missing XYZ file: {xyz_file.name}"
                )
            elif has_valid_name and isinstance(xyz, str):
                # Compare coordinate lines (skip line 2 which is a comment line
                # that may differ between the standalone file and JSON)
                file_lines = xyz_file.read_text(encoding="utf-8").strip().splitlines()
                json_lines = xyz.strip().splitlines()
                if len(file_lines) < 3 or len(json_lines) < 3:
                    errors.append(
                        f"{json_file.name}/{label}: XYZ data must contain an atom "
                        "count, comment line, and coordinates"
                    )
                elif file_lines[0] != json_lines[0]:
                    errors.append(
                        f"{json_file.name}/{label}: Atom count mismatch between "
                        f"JSON ({json_lines[0]}) and {xyz_file.name} ({file_lines[0]})"
                    )
                elif file_lines[2:] != json_lines[2:]:
                    errors.append(
                        f"{json_file.name}/{label}: Coordinate mismatch between "
                        f"JSON and {xyz_file.name}"
                    )

            # Check image file exists
            img_file = dataset_dir / "images" / f"{label}.png"
            if has_valid_name and not img_file.exists():
                errors.append(
                    f"{json_file.name}/{label}: Missing image: {img_file.name}"
                )

            # Check raw output file exists
            out_file = dataset_dir / "raw_output" / f"{label}.out"
            if has_valid_name and not out_file.exists():
                errors.append(
                    f"{json_file.name}/{label}: Missing raw output: {out_file.name}"
                )

            # Basic numeric sanity checks
            if isinstance(sparse_ci, dict):
                n_dets = sparse_ci.get("n_dets")
                if "n_dets" in sparse_ci and (
                    not isinstance(n_dets, int)
                    or isinstance(n_dets, bool)
                    or n_dets < 1
                ):
                    errors.append(
                        f"{json_file.name}/{label}: "
                        f"Invalid n_dets={n_dets} (expected integer >= 1)"
                    )

                delta_e = sparse_ci.get("delta_e_mhartree")
                if "delta_e_mhartree" in sparse_ci and (
                    not isinstance(delta_e, (int, float))
                    or isinstance(delta_e, bool)
                    or delta_e < 0
                ):
                    errors.append(
                        f"{json_file.name}/{label}: "
                        f"Invalid delta_e_mhartree={delta_e} "
                        "(expected non-negative number)"
                    )

            if isinstance(structure, dict):
                num_atoms = structure.get("num_atoms")
                if "num_atoms" in structure and (
                    not isinstance(num_atoms, int)
                    or isinstance(num_atoms, bool)
                    or num_atoms < 1
                ):
                    errors.append(
                        f"{json_file.name}/{label}: Invalid num_atoms={num_atoms}"
                    )

                mass = structure.get("total_mass_amu")
                if "total_mass_amu" in structure and (
                    not isinstance(mass, (int, float))
                    or isinstance(mass, bool)
                    or mass <= 0
                ):
                    errors.append(f"{json_file.name}/{label}: Invalid mass={mass}")

    return errors


def main() -> int:
    """Scan all datasets under data/molecules/ and report errors."""
    repo_root = Path(__file__).resolve().parents[2]
    molecules_dir = repo_root / "data" / "molecules"

    if not molecules_dir.exists():
        print("ERROR: data/molecules/ directory not found")
        return 1

    all_errors: list[str] = []

    for dataset_dir in sorted(molecules_dir.iterdir()):
        if not dataset_dir.is_dir():
            continue
        # Only check directories that contain a JSON file (i.e. actual datasets)
        if not any(
            path.name not in NON_DATA_JSON_FILENAMES
            for path in dataset_dir.glob("*.json")
        ):
            continue
        print(f"Checking {dataset_dir.name}...")
        errors = check_dataset(dataset_dir)
        all_errors.extend(errors)

    if all_errors:
        print(f"\nFAILED: {len(all_errors)} error(s) found:")
        for err in all_errors:
            print(f"  ✗ {err}")
        return 1

    print("\nAll dataset integrity checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
