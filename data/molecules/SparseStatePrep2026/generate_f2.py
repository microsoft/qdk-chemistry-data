"""Regenerate the paper's historical F2 benchmark state from its geometry."""

# --------------------------------------------------------------------------------------------
# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT License. See LICENSE.txt in the project root for license information.
# --------------------------------------------------------------------------------------------

import argparse
import json
import math
import os
from pathlib import Path
from typing import Any

import numpy as np
from qdk_chemistry.algorithms import create
from qdk_chemistry.data import Orbitals, Structure
from qdk_chemistry.utils import compute_valence_space_parameters

_F2_GENERATION_QDK_REVISION = "aab310b108342456bf6b1017d01ccb8e31d2d52d"
_F2_DETERMINANT_TIE_TOLERANCE = 5.0e-7
_F2_COMPATIBILITY_ROTATION_DEGREES = 67.0846143033431
_PAPER_F2_SUPPORT_ORDER = (
    "0001111100011111",
    "0101011101010111",
    "0011101100111011",
    "0111001101110011",
    "0011011101011011",
    "0101101100110111",
    "1000111110001111",
    "0011011100111011",
    "0101011101011011",
    "0101011100110111",
    "0101101100111011",
    "0011011101010111",
    "0011101101011011",
    "0011101100110111",
)


def _determinant_bitstring(determinant: Any, num_orbitals: int) -> str:
    """Serialize a determinant as the benchmark's MSB-first beta-alpha blocks."""
    alpha, beta = determinant.to_binary_strings(num_orbitals)
    return beta[::-1] + alpha[::-1]


def generate_f2_wavefunction(
    xyz_path: Path,
    basis: str = "def2-svp",
) -> dict[str, Any]:
    """Regenerate the paper's historical F2 benchmark state from its geometry."""
    required_threads = {
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
    }
    invalid_threads = {
        name: os.environ.get(name)
        for name, expected in required_threads.items()
        if os.environ.get(name) != expected
    }
    if invalid_threads:
        assignments = " ".join(
            f"{name}={value}" for name, value in required_threads.items()
        )
        raise RuntimeError(
            "F2 regeneration requires single-threaded native solvers. "
            f"Run with {assignments}; received {invalid_threads}."
        )

    structure = Structure.from_xyz_file(xyz_path)
    hf_energy, hf_wavefunction = create("scf_solver").run(
        structure,
        0,
        1,
        basis_or_guess=basis,
    )

    selected_electrons, active_orbitals = compute_valence_space_parameters(
        hf_wavefunction, 0
    )
    if (selected_electrons, active_orbitals) != (14, 8):
        raise RuntimeError(
            "The F2 compatibility gauge requires the historical 14-electron, "
            f"8-orbital selection; received {(selected_electrons, active_orbitals)}."
        )

    selector = create("active_space_selector", "qdk_valence")
    selector.settings().set("num_active_electrons", selected_electrons)
    selector.settings().set("num_active_orbitals", active_orbitals)
    active_wavefunction = selector.run(hf_wavefunction)
    orbitals = active_wavefunction.get_orbitals()
    active_indices = list(orbitals.get_active_space_indices()[0])
    inactive_indices = list(orbitals.get_inactive_space_indices()[0])
    if active_indices != list(range(2, 10)) or inactive_indices != [0, 1]:
        raise RuntimeError(
            "The F2 compatibility gauge requires active full-MO indices 2-9 "
            f"and inactive indices 0-1; received {active_indices} and "
            f"{inactive_indices}."
        )

    coefficients = np.array(orbitals.get_coefficients()[0], copy=True)
    energies = np.array(orbitals.get_energies()[0], copy=True)

    # Recover the paper's arbitrary basis inside the exactly degenerate pi pair.
    coefficients[:, [4, 5]] = coefficients[:, [5, 4]]
    coefficients[:, 3] *= -1.0
    coefficients[:, 6] *= -1.0
    theta = math.radians(_F2_COMPATIBILITY_ROTATION_DEGREES)
    rotation = np.array(
        [
            [math.cos(theta), -math.sin(theta)],
            [math.sin(theta), math.cos(theta)],
        ]
    )
    coefficients[:, [4, 5]] = coefficients[:, [4, 5]] @ rotation
    compatible_orbitals = Orbitals(
        coefficients,
        energies,
        orbitals.get_overlap_matrix(),
        orbitals.get_basis_set(),
        (active_indices, inactive_indices),
    )
    hamiltonian = create("hamiltonian_constructor", "qdk").run(compatible_orbitals)

    active_alpha = active_beta = 5
    casci_energy, casci_wavefunction = create(
        "multi_configuration_calculator", "macis_cas"
    ).run(hamiltonian, active_alpha, active_beta)

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
    num_determinants = len(_PAPER_F2_SUPPORT_ORDER)
    cutoff = abs(ranked[num_determinants - 1][0])
    above_cutoff = [
        pair for pair in ranked if abs(pair[0]) > cutoff + _F2_DETERMINANT_TIE_TOLERANCE
    ]
    cutoff_ties = [
        pair
        for pair in ranked
        if abs(abs(pair[0]) - cutoff) <= _F2_DETERMINANT_TIE_TOLERANCE
    ]
    cutoff_ties.sort(key=lambda pair: _determinant_bitstring(pair[1], active_orbitals))
    selected_pairs = above_cutoff + cutoff_ties[: num_determinants - len(above_cutoff)]
    top_determinants = [determinant for _, determinant in selected_pairs]
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
    if set(projected_by_bitstring) != set(_PAPER_F2_SUPPORT_ORDER):
        missing = sorted(set(_PAPER_F2_SUPPORT_ORDER) - set(projected_by_bitstring))
        unexpected = sorted(set(projected_by_bitstring) - set(_PAPER_F2_SUPPORT_ORDER))
        raise RuntimeError(
            "Deterministic F2 selection did not recover the paper support: "
            f"missing={missing}, unexpected={unexpected}."
        )

    phase = -1.0 if projected_by_bitstring[_PAPER_F2_SUPPORT_ORDER[0]][1] > 0.0 else 1.0
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
        "generation": {
            "qdk_chemistry_revision": _F2_GENERATION_QDK_REVISION,
            "thread_environment": required_threads,
            "selected_active_electrons": selected_electrons,
            "cas_electrons": {
                "alpha": active_alpha,
                "beta": active_beta,
            },
            "determinant_tie_tolerance": _F2_DETERMINANT_TIE_TOLERANCE,
            "compatibility_gauge": {
                "indexing": "zero-based full MO columns",
                "swap": [4, 5],
                "phase_flips": [3, 6],
                "rotation_pair": [4, 5],
                "rotation_degrees": _F2_COMPATIBILITY_ROTATION_DEGREES,
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
                for bitstring in _PAPER_F2_SUPPORT_ORDER
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
        default=Path(__file__).parent / "data" / "structures" / "f2.xyz",
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
