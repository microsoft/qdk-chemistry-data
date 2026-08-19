"""Benchmark different sparse state preparation methods on randomly generated
matrices and chemical wavefunctions.

This script:

1. Generates sparse isometry matrices across a grid of ``(num_qubits, num_configs)``
   using half-filling (``n_electrons = n_orbitals = num_qubits / 2``),
   Jordan-Wigner encoding (via ``generate_determinants_matrix``).
2. Loads chemical wavefunctions from ``data/input_wavefunctions.json``
    (excluding F2) and estimates resources for each.
3. Runs resource estimation for each method on each generated matrix or wavefunction.
4. Collects resource metrics (logical qubits, non-Clifford count, Clifford
   count) and saves a checkpoint JSON and a scaled line plot into ``output/``.
"""

# --------------------------------------------------------------------------------------------
# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT License. See LICENSE in the project root for license information.
# --------------------------------------------------------------------------------------------

import argparse
import json
import math
from collections.abc import Callable
from dataclasses import asdict
from math import comb
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from generate_random_matrix import generate_determinants_matrix
from qdk_chemistry.utils import Logger
from state_preparation_methods import (
    BenchmarkResult,
    Ramacciotti2024,
    ResourceEstimateData,
    Rupprecht2026,
    benchmark_environment,
    bitstring_from_qubit_occupations,
    gf2x,
    gf2x_binary_encoding,
)


def _save_checkpoint(
    data: list[BenchmarkResult],
    output_dir: Path,
    seed: int,
    qubits_list: list[int],
    num_configs_ratio: float,
) -> None:
    """Save current results to JSON atomically.

    Args:
        data: List of :class:`BenchmarkResult` accumulated so far.
        output_dir: Directory to write ``random_matrix_results.json`` into.
        seed: RNG seed recorded in the metadata.
        qubits_list: Qubit grid recorded in the metadata.
        num_configs_ratio: Configs-per-qubit ratio recorded in the metadata.
    """
    results = {
        "metadata": {
            "seed": seed,
            "qubits_list": qubits_list,
            "num_configs_ratio": num_configs_ratio,
            "environment": benchmark_environment(),
        },
        "data": [asdict(d) for d in data],
    }
    json_path = output_dir / "random_matrix_results.json"
    temp_path = json_path.with_suffix(".tmp")
    try:
        with open(temp_path, "w") as f:
            json.dump(results, f, indent=2, default=str)
            f.write("\n")
        temp_path.replace(json_path)
    except Exception as exc:
        Logger.warn(f"Failed to save checkpoint to {json_path}: {exc}")


def _estimate_method(
    method_name: str,
    bitstrings: list[str],
    coeffs: list[float],
    gf2x_max_qubits: int = 25,
) -> tuple[ResourceEstimateData, ResourceEstimateData] | None:
    """Run a single method and return its sparse and dense resource estimates.

    Args:
        method_name: One of ``"gf2x"``, ``"gf2x_binary_encoding"``,
            ``"Rupprecht2026"``, ``"Ramacciotti2024"``.
        bitstrings: MSB-first occupation bitstrings.
        coeffs: Normalised coefficients aligned with *bitstrings*.
        gf2x_max_qubits: Upper qubit limit for the ``"gf2x"`` method.

    Returns:
        ``(sparse_est, dense_est)`` pair, or ``None`` if the method was
        skipped or raised an exception.
    """
    num_qubits = len(bitstrings[0])
    num_configs = len(bitstrings)
    if method_name == "gf2x" and num_qubits > gf2x_max_qubits:
        Logger.info(
            f"Skipping {method_name} for q={num_qubits} (cap={gf2x_max_qubits})"
        )
        return None
    methods: dict[
        str, Callable[..., tuple[ResourceEstimateData, ResourceEstimateData]]
    ] = {
        "gf2x": gf2x,
        "gf2x_binary_encoding": gf2x_binary_encoding,
        "Rupprecht2026": Rupprecht2026,
        "Ramacciotti2024": Ramacciotti2024,
    }
    try:
        fn = methods[method_name]
        return fn(bitstrings, coeffs)

    except Exception as exc:
        Logger.warn(
            f"Method {method_name} failed (q={num_qubits}, configs={num_configs}): {exc}"
        )
    return None


def run_benchmark(
    wfn_json: Path,
    output_dir: Path,
    seed: int = 1000,
    qubits_list: list[int] | None = None,
    num_configs_ratio: float = 1,
    methods: list[str] | None = None,
    gf2x_max_qubits: int = 25,
) -> list[BenchmarkResult]:
    """Run the full benchmark scan over random matrices and chemical wavefunctions.

    Args:
        wfn_json: Path to the JSON file containing chemical wavefunctions.
            Molecules named ``"F2"`` are skipped.
        output_dir: Directory for checkpoint JSON files and cached matrices.
        seed: RNG seed for coefficient generation and matrix sampling.
        qubits_list: Qubit counts to benchmark. Defaults to
            ``list(range(4, 20, 4)) + list(range(20, 61, 10))``.
        num_configs_ratio: ``num_configs = ceil(ratio * num_qubits)``.
        methods: Method names to run. Defaults to all four methods.
        gf2x_max_qubits: Skip ``"gf2x"`` for systems larger than this.

    Returns:
        List of :class:`BenchmarkResult`, one per (method, system) pair.
    """
    if qubits_list is None:
        qubits_list = list(range(4, 20, 4)) + list(range(20, 61, 10))
    if methods is None:
        methods = METHOD_ORDER
    coeff_rng = np.random.default_rng(seed=seed)

    data: list[BenchmarkResult] = []

    # 1. Random Matrix Benchmark
    for num_qubits in qubits_list:
        n_orbitals = num_qubits // 2
        n_electrons = n_orbitals
        n_alpha = n_electrons // 2
        n_beta = n_electrons - n_alpha
        max_possible = comb(n_orbitals, n_alpha) * comb(n_orbitals, n_beta)
        num_configs = min(
            max(2, math.ceil(num_configs_ratio * num_qubits)), max_possible
        )

        Logger.info(
            f"Generating matrix: num_qubits={num_qubits}, n_orbitals={n_orbitals}, "
            f"n_electrons={n_electrons}, num_configs={num_configs}"
        )
        full_matrix = generate_determinants_matrix(
            n_electrons=n_electrons,
            n_orbitals=n_orbitals,
            n_dets=num_configs,
            seed=seed,
        )

        # Generate real coefficients
        coeffs_raw = coeff_rng.random(num_configs) - 0.5
        coeffs = coeffs_raw / np.linalg.norm(coeffs_raw)

        # Convert binary matrix rows to bitstrings and build a QDK Wavefunction
        bitstrings = [
            bitstring_from_qubit_occupations(row) for row in full_matrix[:num_configs]
        ]
        Logger.info(f"  q={num_qubits} configs={num_configs}")

        for method_name in methods:
            Logger.info(f"    Running {method_name} ...")
            est = _estimate_method(
                method_name, bitstrings, list(coeffs), gf2x_max_qubits
            )
            if est is None:
                continue

            sparse_est, dense_est = est
            data.append(
                BenchmarkResult(
                    method=method_name,
                    num_qubits=num_qubits,
                    num_dets=num_configs,
                    sparse=sparse_est,
                    dense=dense_est,
                    source="random",
                )
            )
            _save_checkpoint(
                data,
                output_dir,
                seed,
                qubits_list,
                num_configs_ratio,
            )

    # 2. Chemical Data Benchmark (from data/input_wavefunctions.json)
    if wfn_json.exists():
        with open(wfn_json, "r") as f:
            all_wfns = json.load(f)

        for mol_name, mol_data in all_wfns.items():
            # Skip F2 (handled by estimate_f2.py)
            if mol_name == "F2":
                continue

            bitstrings = mol_data["bitstrings"]
            coeffs = mol_data["coeffs"]
            num_qubits = mol_data["n_qubits"]
            num_configs = len(bitstrings)

            # Normalize coefficients
            norm = np.linalg.norm(coeffs)
            if norm > 0:
                coeffs = [c / norm for c in coeffs]

            Logger.info(f"Chemical: {mol_name}  q={num_qubits} configs={num_configs}")

            for method_name in methods:
                Logger.info(f"    Running {method_name} ...")
                est = _estimate_method(
                    method_name, bitstrings, list(coeffs), gf2x_max_qubits
                )
                if est is None:
                    continue

                sparse_est, dense_est = est
                data.append(
                    BenchmarkResult(
                        method=method_name,
                        num_qubits=num_qubits,
                        num_dets=num_configs,
                        sparse=sparse_est,
                        dense=dense_est,
                        source="chemical",
                        molecule=mol_name,
                    )
                )
                _save_checkpoint(
                    data,
                    output_dir,
                    seed,
                    qubits_list,
                    num_configs_ratio,
                )
    else:
        Logger.warn(f"Chemical data not found at {wfn_json} — skipping")

    return data


# Plotting Helpers
METHOD_LABELS = {
    "gf2x": "Ours w/o binary encoding",
    "gf2x_binary_encoding": "Ours w/ binary encoding",
    "Rupprecht2026": "Rupprecht2026",
    "Ramacciotti2024": "Ramacciotti2024",
}
METHOD_COLORS = {
    "gf2x": "#D55E00",
    "gf2x_binary_encoding": "#0072B2",
    "Rupprecht2026": "#CC79A7",
    "Ramacciotti2024": "#E69F00",
}
METHOD_MARKERS = {
    "gf2x": "D",
    "gf2x_binary_encoding": "^",
    "Rupprecht2026": "s",
    "Ramacciotti2024": "o",
}
# Fixed ordering so legend always matches visual curve order
METHOD_ORDER = [
    "gf2x",
    "Ramacciotti2024",
    "gf2x_binary_encoding",
    "Rupprecht2026",
]
# Markers for chemical data points, assigned by molecule name
MOLECULE_MARKERS = [
    "v",
    "<",
    ">",
    "P",
    "*",
    "X",
    "p",
    "h",
    "H",
    "d",
    "8",
    "s",
    "D",
    "1",
    "2",
    "3",
    "4",
    "+",
    "x",
]


def plot_scaled(data: list[BenchmarkResult], output_path: Path) -> None:
    """Plot example fits of resource counts vs num_qubits for each method."""
    subset = [d for d in data if d.source == "random"]
    subset_chem = [
        d
        for d in data
        if d.source == "chemical" and not (d.molecule or "").startswith("mol")
    ]

    if not subset:
        Logger.warn("No data found, skipping scaled plot.")
        return

    methods = [m for m in METHOD_ORDER if m in set(d.method for d in subset)]
    chem_molecules = sorted(set(d.molecule or "" for d in subset_chem))
    mol_to_marker = {
        m: MOLECULE_MARKERS[i % len(MOLECULE_MARKERS)]
        for i, m in enumerate(chem_molecules)
    }

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # (ylabel, is_symlog, value function, line alpha) per panel
    panels = [
        (
            "Ancillary Qubits",
            True,
            lambda d: d.combined.logical_qubits - d.num_qubits,
            0.9,
        ),
        ("Non-Clifford Count", False, lambda d: d.combined.non_clifford_count, 0.7),
        ("Clifford Count", False, lambda d: d.combined.clifford_count, 0.7),
    ]

    for ax, (ylabel, is_symlog, value_fn, line_alpha) in zip(axes, panels):
        for method in methods:
            mkr = METHOD_MARKERS.get(method, "o")
            clr = METHOD_COLORS[method]
            pts = sorted(
                (d for d in subset if d.method == method),
                key=lambda d: d.num_qubits,
            )
            xs = [d.num_qubits for d in pts]
            ys = [value_fn(d) for d in pts]

            ax.plot(
                xs,
                ys,
                linestyle="--",
                marker=mkr,
                label=METHOD_LABELS.get(method, method),
                color=clr,
                alpha=line_alpha,
            )

            # Chemical data — per molecule
            chem_pts = [d for d in subset_chem if d.method == method]
            for mol in chem_molecules:
                mol_pts = [d for d in chem_pts if d.molecule == mol]
                if mol_pts:
                    xs_chem = [d.num_qubits for d in mol_pts]
                    ys_chem = [value_fn(d) for d in mol_pts]
                    ax.scatter(
                        xs_chem,
                        ys_chem,
                        marker=mol_to_marker[mol],
                        color=clr,
                        edgecolors="black",
                        linewidths=1.0,
                        s=80,
                        alpha=0.85,
                        zorder=5,
                    )

        ax.set_xlabel("Number of Qubits", fontsize=16)
        ax.set_ylabel(ylabel, fontsize=16)
        ax.tick_params(axis="both", which="major", labelsize=14)
        if is_symlog:
            ax.set_yscale("symlog", linthresh=1, linscale=0.05)
            ax.set_ylim(bottom=0)
        else:
            ax.set_yscale("log")
        ax.grid(True, alpha=0.3)

    axes[2].legend(fontsize=14)

    # Molecule shape legend (center top) — black edge distinguishes chemical data
    if chem_molecules:
        mol_handles = [
            plt.Line2D(
                [0],
                [0],
                marker=mol_to_marker[m],
                color="w",
                markerfacecolor="gray",
                markeredgecolor="black",
                markeredgewidth=1.0,
                markersize=10,
                label=m,
            )
            for m in chem_molecules
        ]
        fig.legend(
            handles=mol_handles,
            title="Data Source",
            loc="lower center",
            bbox_to_anchor=(0.5, 1.0),
            ncol=max(1, (len(chem_molecules) + 1) // 2),
            framealpha=0.5,
            fontsize=12,
            title_fontsize=12,
        )

    plt.tight_layout()
    plt.savefig(output_path, dpi=600, bbox_inches="tight")
    plt.close()
    Logger.info(f"Saved scaled plot: {output_path}")


def main() -> None:
    """Parse command-line arguments and run the benchmark.

    Args are taken from ``sys.argv``.  Run with ``--help`` for usage.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Benchmark sparse state preparation methods on random and chemical "
            "wavefunctions."
        )
    )
    parser.add_argument(
        "--wfn_path",
        type=Path,
        default=Path(__file__).parent / "data" / "input_wavefunctions.json",
        metavar="PATH",
        help="Path to the molecule wavefunction JSON file. Default: %(default)s",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).parent / "output",
        metavar="DIR",
        help="Directory for JSON results and cached matrices. Default: %(default)s",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        metavar="N",
        help="RNG seed for coefficient generation. Default: %(default)s",
    )
    parser.add_argument(
        "--qubits",
        type=int,
        nargs="+",
        default=list(range(4, 20, 4)) + list(range(20, 61, 10)),
        metavar="N",
        help="Space-separated list of qubit counts to benchmark. Default: %(default)s",
    )
    parser.add_argument(
        "--configs-ratio",
        type=float,
        default=1.0,
        metavar="K",
        help="num_configs = ceil(K * num_qubits). Default: %(default)s",
    )
    parser.add_argument(
        "--max-gf2x-qubits",
        type=int,
        default=25,
        metavar="N",
        help="Skip 'gf2x' for systems larger than N qubits. Default: %(default)s",
    )
    args = parser.parse_args()
    output_dir: Path = args.output_dir
    figures_dir: Path = output_dir / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    Logger.info("Starting benchmark scan ...")
    data = run_benchmark(
        wfn_json=args.wfn_path,
        output_dir=output_dir,
        seed=args.seed,
        qubits_list=args.qubits,
        num_configs_ratio=args.configs_ratio,
        methods=METHOD_ORDER,
        gf2x_max_qubits=args.max_gf2x_qubits,
    )

    if data:
        plot_scaled(data, figures_dir / "random_matrix_results.png")
    else:
        Logger.warn("No data collected — skipping plots")

    Logger.info(f"Done. {len(data)} data points collected.")


if __name__ == "__main__":
    """Run the benchmark when this script is executed directly.

    The command-line to create the full benchmark and plot is:
    python estimate_random_matrix.py
    """
    Logger.set_global_level("info")
    main()
