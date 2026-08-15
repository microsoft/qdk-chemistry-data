#!/usr/bin/env python3
"""Generate visualization of sparse CI wavefunction analysis results.

Creates a scatter plot showing:
- X-axis: Number of determinants
- Y-axis: Multi-reference character (100 - |c₀|² × 100)
- Marker color: Energy error (ΔE)

All metrics are derived directly from SparseCI-24.json.

Usage::

    python plot_analysis.py

Requirements::

    pip install matplotlib numpy pandas
"""

# --------------------------------------------------------------------------------------------
# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT License. See LICENSE in the project root for license information.
# --------------------------------------------------------------------------------------------

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.style.use("seaborn-v0_8-whitegrid")


# =============================================================================
# Data loading
# =============================================================================


def load_data(json_path: Path) -> pd.DataFrame:
    """Load SparseCI-24.json and derive analysis columns."""
    with open(json_path) as f:
        raw = json.load(f)

    rows = []
    for entry in raw:
        sci = entry["sparse_ci_finder"]
        dets = sci["determinants"]
        c0 = max(abs(d["coeff"]) for d in dets)
        single_ref = c0**2 * 100
        multi_ref = 100 - single_ref

        rows.append(
            {
                "name": entry["name"],
                "n_dets": sci["n_dets"],
                "delta_e_mhartree": sci["delta_e_mhartree"],
                "energy_hartree": sci["energy_hartree"],
                "single_ref_character": single_ref,
                "multi_ref_character": multi_ref,
                "num_atoms": entry.get("structure", {}).get("num_atoms"),
                "composition": entry.get("structure", {}).get("composition"),
            }
        )

    return pd.DataFrame(rows)


# =============================================================================
# Scatter plot
# =============================================================================


def create_main_plot(df: pd.DataFrame, output_path: Path) -> None:
    """Scatter plot: x = n_dets, y = multi-ref character, color = ΔE."""
    df = df.reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(12, 8))

    # Small jitter so overlapping points separate
    rng = np.random.default_rng(42)
    jitter = rng.uniform(-0.15, 0.15, len(df))

    vmin, vmax = 0, df["delta_e_mhartree"].max() * 1.1
    cmap = "YlOrRd"

    scatter = ax.scatter(
        df["n_dets"] + jitter,
        df["multi_ref_character"],
        c=df["delta_e_mhartree"],
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        marker="o",
        s=150,
        alpha=0.85,
        edgecolors="black",
        linewidth=0.8,
        zorder=3,
    )

    cbar = plt.colorbar(scatter, ax=ax, pad=0.02, aspect=25)
    cbar.set_label("ΔE: Sparse CI − CASCI (mHartree)", fontsize=11)

    # Label notable multi-ref systems
    for idx, row in df.iterrows():
        if row["multi_ref_character"] > 10:
            name = row["name"].replace("_diradical", "")
            if "ozone" in name:
                name = "O₃"
            elif len(name) > 10:
                name = name[:8] + ".."
            ax.annotate(
                name,
                (row["n_dets"] + jitter[idx], row["multi_ref_character"]),
                xytext=(6, 2),
                textcoords="offset points",
                fontsize=8,
                alpha=0.9,
            )

    ax.set_xlabel("Number of Determinants", fontsize=12)
    ax.set_ylabel("Multi-Reference Character  (100 − |c₀|² × 100)  %", fontsize=12)
    ax.set_title("SparseCI-24 — Sparse CI Wavefunction Analysis", fontsize=14)

    ax.set_xlim(1.4, 6.8)
    ax.set_xticks([2, 3, 4, 5, 6])
    ax.set_ylim(-2, 55)

    plt.tight_layout()

    plt.savefig(output_path, dpi=150, bbox_inches="tight", facecolor="white")
    print(f"Plot saved to: {output_path}")

    plt.close()


# =============================================================================
# Main
# =============================================================================


def main():
    script_dir = Path(__file__).parent
    json_path = script_dir / "SparseCI-24.json"

    if not json_path.exists():
        print(f"Error: {json_path} not found")
        return

    df = load_data(json_path)
    print(f"Loaded data for {len(df)} molecules")
    create_main_plot(df, script_dir / "sparse_ci_analysis.png")
    print("Done.")


if __name__ == "__main__":
    main()
