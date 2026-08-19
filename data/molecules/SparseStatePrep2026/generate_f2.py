"""Regenerate the F2 wavefunction from its geometry."""

# --------------------------------------------------------------------------------------------
# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT License. See LICENSE in the project root for license information.
# --------------------------------------------------------------------------------------------

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
from qdk_chemistry.algorithms import create
from qdk_chemistry.data import Orbitals, Structure
from qdk_chemistry.utils import compute_valence_space_parameters

_F2_NUM_DETERMINANTS = 14
_F2_ACTIVE_ELECTRONS = 10
_F2_ACTIVE_ORBITALS = 8


def _determinant_bitstring(determinant: Any, num_orbitals: int) -> str:
    """Serialize a determinant as MSB-first beta-alpha blocks."""
    alpha, beta = determinant.to_binary_strings(num_orbitals)
    return beta[::-1] + alpha[::-1]


def generate_f2_wavefunction(
    xyz_path: Path,
    basis: str = "def2-svp",
) -> dict[str, Any]:
    """Regenerate the F2 wavefunction from xyz."""
    structure = Structure.from_xyz_file(xyz_path)
    hf_energy, hf_wavefunction = create("scf_solver").run(
        structure,
        0,
        1,
        basis_or_guess=basis,
    )

    assert compute_valence_space_parameters(hf_wavefunction, 0) == (14, 8), (
        f"Expected a (14e, 8o) valence space, "
        f"got {compute_valence_space_parameters(hf_wavefunction, 0)}"
    )

    active_orbitals = _F2_ACTIVE_ORBITALS
    selector = create("active_space_selector", "qdk_valence")
    selector.settings().set("num_active_electrons", _F2_ACTIVE_ELECTRONS)
    selector.settings().set("num_active_orbitals", active_orbitals)
    active_wavefunction = selector.run(hf_wavefunction)
    orbitals = active_wavefunction.get_orbitals()
    active_indices = list(orbitals.get_active_space_indices()[0])
    inactive_indices = list(orbitals.get_inactive_space_indices()[0])
    assert active_indices == list(range(4, 12)) and inactive_indices == [0, 1, 2, 3], (
        f"Expected active indices 4-11 and inactive indices [0, 1, 2, 3], "
        f"got {active_indices} and {inactive_indices}"
    )

    coefficients = np.array(orbitals.get_coefficients()[0], copy=True)
    energies = np.array(orbitals.get_energies()[0], copy=True)

    compatible_orbitals = Orbitals(
        coefficients,
        energies,
        orbitals.get_overlap_matrix(),
        orbitals.get_basis_set(),
        (active_indices, inactive_indices),
    )
    hamiltonian = create("hamiltonian_constructor", "qdk").run(compatible_orbitals)

    active_alpha = active_beta = _F2_ACTIVE_ELECTRONS // 2
    casci_energy, casci_wavefunction = create(
        "multi_configuration_calculator", "macis_cas"
    ).run(hamiltonian, active_alpha, active_beta)
    assert casci_energy <= hf_energy, (
        f"CASCI energy {casci_energy} is above the RHF energy {hf_energy}; "
        "the active space is inconsistent with the electron count"
    )

    coefficient_determinant_pairs = list(
        zip(
            casci_wavefunction.get_coefficients(),
            casci_wavefunction.get_active_determinants(),
            strict=True,
        )
    )
    ranked = sorted(
        coefficient_determinant_pairs,
        key=lambda pair: -abs(pair[0]),
    )
    selected_pairs = ranked[:_F2_NUM_DETERMINANTS]
    top_determinants = [determinant for _, determinant in selected_pairs]
    selected_order = [
        _determinant_bitstring(determinant, active_orbitals)
        for determinant in top_determinants
    ]
    projected_energy, projected_wavefunction = create(
        "projected_multi_configuration_calculator", "macis_pmc"
    ).run(hamiltonian, top_determinants)

    determinants = projected_wavefunction.get_active_determinants()
    coefficients = [
        complex(value) for value in projected_wavefunction.get_coefficients()
    ]
    max_imaginary = max((abs(value.imag) for value in coefficients), default=0.0)
    if max_imaginary > 1e-10:
        raise ValueError(
            "F2 generation produced complex coefficients "
            f"(max |imag| = {max_imaginary:.2e})."
        )
    projected_by_bitstring = {
        _determinant_bitstring(determinant, active_orbitals): (
            determinant,
            coefficient.real,
        )
        for determinant, coefficient in zip(determinants, coefficients, strict=True)
    }
    phase = 1.0 if projected_by_bitstring[selected_order[0]][1] > 0.0 else -1.0
    total_orbitals = len(hf_wavefunction.get_orbitals().get_energies_alpha())

    return {
        "structure": {
            "num_atoms": structure.get_num_atoms(),
            "composition": "F2",
            "total_mass_amu": structure.get_total_mass(),
            "nuclear_repulsion_energy_eh": (
                structure.calculate_nuclear_repulsion_energy()
            ),
        },
        "scf_energy_hartree": hf_energy,
        "orbitals_summary": {
            "aos": compatible_orbitals.get_num_atomic_orbitals(),
            "mos": total_orbitals,
            "active_orbitals": {
                "alpha": active_orbitals,
                "beta": active_orbitals,
            },
            "inactive_orbitals": {
                "alpha": len(inactive_indices),
                "beta": len(inactive_indices),
            },
        },
        "casci_energies_hartree": [casci_energy],
        "initial_casci_energy_hartree": casci_energy,
        "hamiltonian_summaries": [
            {
                "active_orbitals": active_orbitals,
                "total_orbitals": total_orbitals,
                "core_energy": hamiltonian.get_core_energy(),
            }
        ],
        "sparse_ci_finder": {
            "n_dets": len(determinants),
            "energy_hartree": projected_energy,
            "delta_e_mhartree": 1000 * (projected_energy - casci_energy),
            "determinants": [
                {
                    "det": str(projected_by_bitstring[bitstring][0]),
                    "coeff": phase * projected_by_bitstring[bitstring][1],
                }
                for bitstring in selected_order
            ],
        },
        "name": "F2",
        "xyz": xyz_path.read_text(),
    }


def main() -> None:
    """Parse command-line arguments and regenerate the F2 record."""
    parser = argparse.ArgumentParser(
        description="Regenerate the paper's historical F2 benchmark state."
    )
    parser.add_argument(
        "--xyz",
        default=Path(__file__).parent / "data" / "xyz" / "f2.xyz",
        type=Path,
        help="F2 geometry. Default: %(default)s",
    )
    parser.add_argument(
        "--output",
        default=Path(__file__).parent / "data" / "f2.json",
        type=Path,
        help="Generated F2 record. Default: %(default)s",
    )
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps([generate_f2_wavefunction(args.xyz)], indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
