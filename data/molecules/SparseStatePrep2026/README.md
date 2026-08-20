# Sparse State Preparation 2026 Scripts and Data

Scripts and figures for the paper "Clifford-efficient sparse state
preparation for molecular wavefunctions".

The results compare four sparse state preparation methods on random
determinant matrices and molecular wavefunctions.

## Figures

![Random and molecular resource estimates](output/figures/random_matrix_results.png)

![F2 resource estimates](output/figures/f2_matrix_results.png)

![F2 sparse and dense resource breakdown](output/figures/f2_matrix_results_stacked.png)

## Reproducibility

Use Python 3.12 in the qdk-chemistry dev container. Install the pinned main
dependencies:

```bash
python -m pip install --requirement requirements-lock.txt
```

`requirements-lock.txt` pins the direct dependencies only; their transitive
dependencies are resolved at install time and may drift.

The Rupprecht and Wolk reference implementation is distributed under Apache
License 2.0 through immutable
[Zenodo record 18234600](https://doi.org/10.5281/zenodo.18234600).

Run the resource estimates from this directory:

```bash
python estimate_random_matrix.py
python estimate_f2.py
```

## Methodology

The eight wavefunctions in `data/input_wavefunctions.json` were extracted from
records in the [SparseCI-24 dataset](../SparseCI-24/).
The paper's F2 benchmark is regenerated from `data/xyz/f2.xyz` with
def2-SVP restricted HF in `generate_f2.py`.

The random benchmark uses a fixed seed of 42. It constructs a half-filled
system with the number of configurations equal to the number of qubits and
samples excitations from the Hartree-Fock determinant.

The compared methods are:

- `gf2x`: QDK/Chemistry GF2+X sparse isometry.
- `gf2x_binary_encoding`: GF2+X with binary encoding.
- `Rupprecht2026`: batched sparse isometry from Rupprecht and Wolk 2026.
- `Ramacciotti2024`: permutation-based sparse state preparation from
  Ramacciotti et al. 2024.

## File Layout

```text
SparseStatePrep2026/
├── README.md
├── estimate_f2.py                    # Run detailed benchmark for F2 molecule
├── estimate_random_matrix.py         # Run random and molecular benchmarks
├── generate_f2.py                    # Regenerate the paper F2 wavefunction
├── generate_random_matrix.py         # Generate random determinant matrices
├── requirements-lock.txt             # Pinned main dependencies
├── state_preparation_methods.py      # Resource estimators for four methods
├── cgmanifest.json                   # Third-party component declarations
├── data/
│   ├── f2.json                       # F2 molecular record and wavefunction
│   ├── input_wavefunctions.json      # Eight SparseCI-24 wavefunctions
│   └── xyz/
│       └── f2.xyz                    # F2 molecular geometry
└── output/
    ├── random_matrix_results.json    # Random and molecular benchmark results
    ├── f2_matrix_results.json        # F2 benchmark results
    └── figures/                      # Three generated PNG figures
```

## Citation

- Rupprecht and Wolk, [Sparse quantum state preparation with improved Toffoli
  cost](https://arxiv.org/abs/2601.09388) (2026).
- Ramacciotti et al., [A simple quantum algorithm to efficiently prepare sparse
  states](https://arxiv.org/abs/2310.19309) (2024).
- [QDK/Chemistry](https://github.com/microsoft/qdk-chemistry).
- [Qualtran](https://github.com/quantumlib/Qualtran).

## License

See the repository [LICENSE](../../../LICENSE).
