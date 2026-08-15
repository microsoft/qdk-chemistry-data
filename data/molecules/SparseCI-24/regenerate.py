#!/usr/bin/env python3
"""Regenerate SparseCI-24 dataset from XYZ structures.

Usage::

  python regenerate.py                            # Run full pipeline (requires qdk-chemistry)
  python regenerate.py --collate-only             # Just re-parse existing raw_output/ logs
  python regenerate.py --images-only              # Just regenerate images (needs rdkit)

Requirements for full regeneration:
  - qdk-chemistry installed: pip install qdk-chemistry
    (or from source: https://github.com/microsoft/qdk-chemistry @ aab310b)
  - Dependencies: pyscf, macis

Images require: rdkit
"""

# --------------------------------------------------------------------------------------------
# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT License. See LICENSE in the project root for license information.
# --------------------------------------------------------------------------------------------

import argparse
import json
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

QDK_CHEMISTRY_COMMIT = "aab310b108342456bf6b1017d01ccb8e31d2d52d"
WORKFLOW_ARGS = [
    "--initial-active-space-solver",
    "macis_asci",
    "--autocas",
    "--max-determinants",
    "100",
]

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Draw, rdDetermineBonds

    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False


# =============================================================================
# Log parsing - extract structured data from workflow stdout
# =============================================================================


def parse_log(text: str) -> dict[str, Any]:
    """Parse workflow log text into structured data."""
    data: dict[str, Any] = {}

    # Structure Summary
    if m := re.search(r"Number of atoms: (\d+)", text):
        data.setdefault("structure", {})["num_atoms"] = int(m.group(1))
    if m := re.search(r"Composition: (.+)", text):
        data.setdefault("structure", {})["composition"] = m.group(1).strip()
    if m := re.search(r"Total mass: ([\d.]+) AMU", text):
        data.setdefault("structure", {})["total_mass_amu"] = float(m.group(1))
    if m := re.search(r"Nuclear Repulsion Energy: ([\d.]+) Eh", text):
        data.setdefault("structure", {})["nuclear_repulsion_energy_eh"] = float(
            m.group(1)
        )

    # SCF energy
    if m := re.search(r"SCF energy = ([-\d.]+) Hartree", text):
        data["scf_energy_hartree"] = float(m.group(1))

    # Orbitals Summary
    orb: dict[str, Any] = {}
    if m := re.search(r"AOs: (\d+)", text):
        orb["aos"] = int(m.group(1))
    if m := re.search(r"MOs: (\d+)", text):
        orb["mos"] = int(m.group(1))
    if m := re.search(r"Type: (\w+)", text):
        orb["type"] = m.group(1)
    if m := re.search(r"Active Orbitals: α=(\d+), β=(\d+)", text):
        orb["active_orbitals"] = {"alpha": int(m.group(1)), "beta": int(m.group(2))}
    if m := re.search(r"Inactive Orbitals: α=(\d+), β=(\d+)", text):
        orb["inactive_orbitals"] = {"alpha": int(m.group(1)), "beta": int(m.group(2))}
    if orb:
        data["orbitals_summary"] = orb

    # CASCI energies
    casci = re.findall(r"CASCI energy = ([-\d.]+) Hartree", text)
    if casci:
        data["casci_energies_hartree"] = [float(e) for e in casci]
        data["initial_casci_energy_hartree"] = float(casci[0])

    # AutoCAS
    if m := re.search(r"AutoCAS energy = ([-\d.]+) Hartree", text):
        data["autocas_energy_hartree"] = float(m.group(1))
    if m := re.search(r"AutoCAS.*selected active space.*indices: \[([\d, ]+)\]", text):
        data["autocas_selected_indices"] = [int(x) for x in m.group(1).split(",")]

    # Hamiltonian summaries
    hams = []
    for m in re.finditer(
        r"Active Orbitals: (\d+)\s+Total Orbitals: (\d+)\s+Core Energy: ([-\d.]+)", text
    ):
        hams.append(
            {
                "active_orbitals": int(m.group(1)),
                "total_orbitals": int(m.group(2)),
                "core_energy": float(m.group(3)),
            }
        )
    if hams:
        data["hamiltonian_summaries"] = hams

    # Sparse CI finder
    if m := re.search(
        r"Sparse CI finder \((\d+) dets\) = ([-\d.]+) Hartree \(ΔE = ([\d.]+) mHartree\)",
        text,
    ):
        sci: dict[str, Any] = {
            "n_dets": int(m.group(1)),
            "energy_hartree": float(m.group(2)),
            "delta_e_mhartree": float(m.group(3)),
        }
        # Parse determinants
        if dm := re.search(
            r"Stored wavefunction with \d+ determinants\s*\nDeterminants:\s*\n((?:  .+\n)+)",
            text,
        ):
            dets = []
            for line in dm.group(1).strip().split("\n"):
                if ":" in line:
                    det, coeff = line.strip().split(":")
                    dets.append({"det": det.strip(), "coeff": float(coeff.strip())})
            sci["determinants"] = dets
        data["sparse_ci_finder"] = sci

    # Timing
    if m := re.search(r"Start: ([\d\-T:.]+)", text):
        data.setdefault("workflow_timing", {})["start"] = m.group(1)
    if m := re.search(r"End: ([\d\-T:.]+)", text):
        data.setdefault("workflow_timing", {})["end"] = m.group(1)
    if m := re.search(r"Elapsed seconds: ([\d.]+)", text):
        data.setdefault("workflow_timing", {})["elapsed_seconds"] = float(m.group(1))

    return data


# =============================================================================
# Workflow execution
# =============================================================================


def find_workflow_script():
    """Find sample_sci_workflow.py from installed qdk-chemistry."""
    try:
        import qdk_chemistry

        pkg_dir = Path(qdk_chemistry.__file__).parent
        candidates = [
            pkg_dir / "sample_sci_workflow.py",
            pkg_dir.parent / "sample_sci_workflow.py",
            pkg_dir.parent / "samples" / "sample_sci_workflow.py",
        ]
        for p in candidates:
            if p.exists():
                return p
        raise FileNotFoundError(f"sample_sci_workflow.py not found near {pkg_dir}")
    except ImportError:
        sys.exit(
            "qdk-chemistry not installed.\n"
            "Install with: pip install qdk-chemistry\n"
            "Or from source: https://github.com/microsoft/qdk-chemistry"
        )


def run_workflow(xyz, log_path, script):
    cmd = [sys.executable, str(script), "--xyz", str(xyz.resolve()), *WORKFLOW_ARGS]

    print(f"  {xyz.stem}...", end=" ", flush=True)
    t0 = datetime.now()

    with open(log_path, "w") as f:
        f.write(
            f"=== WORKFLOW START ===\nStart: {t0.isoformat()}\nCommand: {' '.join(str(c) for c in cmd)}\n\n"
        )
        r = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
        f.write(
            f"\n\n=== WORKFLOW END ===\nEnd: {datetime.now().isoformat()}\nElapsed seconds: {(datetime.now() - t0).total_seconds():.2f}\n"
        )

    print(
        f"{'OK' if r.returncode == 0 else 'FAILED'} ({(datetime.now() - t0).total_seconds():.0f}s)"
    )
    return r.returncode == 0


# =============================================================================
# Image generation
# =============================================================================


def xyz_to_image(xyz_path, img_path):
    if not HAS_RDKIT:
        return False
    try:
        mol = Chem.MolFromXYZBlock(open(xyz_path).read())
        if not mol:
            return False
        rdDetermineBonds.DetermineBonds(Chem.Mol(mol), charge=0)
        mol2d = Chem.RemoveHs(mol)
        mol2d.RemoveAllConformers()
        AllChem.Compute2DCoords(mol2d)
        d = Draw.MolDraw2DCairo(300, 300)
        d.drawOptions().clearBackground = False
        d.drawOptions().bondLineWidth = 2
        d.DrawMolecule(mol2d)
        d.FinishDrawing()
        open(img_path, "wb").write(d.GetDrawingText())
        return True
    except Exception:
        return False


# =============================================================================
# Collation
# =============================================================================


def collate(xyz_dir, raw_dir, out_json, img_dir):
    img_dir.mkdir(exist_ok=True)
    molecules = []

    for xyz in sorted(xyz_dir.glob("*.xyz")):
        name = xyz.stem
        log = raw_dir / f"{name}.out"
        if not log.exists():
            print(f"  {name}: SKIPPED (no log)")
            continue

        print(f"  {name}...", end=" ", flush=True)

        log_text = open(log).read()
        data = parse_log(log_text)
        data["name"] = name
        data["output_file_name"] = f"{name}.out"
        data["output_file_text"] = log_text

        # Read xyz, clear comment line
        lines = open(xyz).readlines()
        lines[1] = "\n"
        data["xyz"] = "".join(lines)

        # Generate image
        if HAS_RDKIT:
            xyz_to_image(xyz, img_dir / f"{name}.png")

        molecules.append(data)
        print("OK")

    molecules.sort(key=lambda m: m["name"])
    json.dump(molecules, open(out_json, "w"), indent=2)
    print(f"\nWrote {len(molecules)} molecules to {out_json}")


# =============================================================================
# Main
# =============================================================================


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s (QDK/Chemistry commit {QDK_CHEMISTRY_COMMIT})",
    )
    p.add_argument(
        "--collate-only",
        action="store_true",
        help="Just re-parse existing raw_output/ logs",
    )
    p.add_argument("--images-only", action="store_true")
    p.add_argument("--molecules", nargs="+")
    args = p.parse_args()

    here = Path(__file__).parent.resolve()
    xyz_dir, raw_dir = here / "xyz", here / "raw_output"

    if args.images_only:
        if not HAS_RDKIT:
            sys.exit("rdkit required for --images-only")
        for xyz in sorted(xyz_dir.glob("*.xyz")):
            print(f"  {xyz.stem}...", end=" ")
            print(
                "OK"
                if xyz_to_image(xyz, here / "images" / f"{xyz.stem}.png")
                else "FAILED"
            )
        return

    if args.collate_only:
        collate(xyz_dir, raw_dir, here / "SparseCI-24.json", here / "images")
        return

    script = find_workflow_script()
    xyzs = sorted(xyz_dir.glob("*.xyz"))
    if args.molecules:
        xyzs = [x for x in xyzs if x.stem in args.molecules]

    print(f"Processing {len(xyzs)} molecules with {script}")
    raw_dir.mkdir(exist_ok=True)

    for xyz in xyzs:
        run_workflow(xyz, raw_dir / f"{xyz.stem}.out", script)

    collate(xyz_dir, raw_dir, here / "SparseCI-24.json", here / "images")


if __name__ == "__main__":
    main()
