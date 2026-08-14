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

Install the script dependencies:

```bash
pip install qdk-chemistry[qiskit] qualtran numpy matplotlib qiskit
```

The Rupprecht and Wolk reference implementation is distributed under Apache
License 2.0 through [Zenodo record 18234600](https://zenodo.org/records/18234600).
Install it as the importable `sparse_state_preparation` package before running
the resource estimation benchmark.

Run the full resource estimates from this directory:

```bash
python estimate_random_matrix.py
python estimate_f2.py
```

## Methodology

The wavefunctions in `data/input_wavefunctions.json` come from two sources.
Eight entries were extracted from records in the
[SparseCI-24 dataset](../SparseCI-24/). The fluorine entry is the
14-configuration truncated SCI wavefunction with the accompanying xyz geometry.

The random benchmark uses a fixed seed of 42. It constructs a half-filled
system with the number of configurations equal to the number of qubits and
samples excitations from the Hartree-Fock determinant.

Benchmark bitstrings are MSB-first. Backend adapters convert them to native
ordering, including QDK's `q[0]`-first configuration strings.

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
├── generate_random_matrix.py         # Generate random determinant matrices
├── state_preparation_methods.py      # Resource estimators for four methods
├── data/
│   ├── input_wavefunctions.json      # Nine molecular wavefunctions
│   └── structures/
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
