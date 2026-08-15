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

Wavefunction generation and resource estimation use separate qdk-chemistry
revisions:

- `aab310b108342456bf6b1017d01ccb8e31d2d52d` generates `data/f2.json`.
- `39ea175d191324e700a775caa9b933fe93d885d8` estimates resources with binary
  encoding from qdk-chemistry PR 435.

Use Python 3.12 in the qdk-chemistry development container. Install the
resource-estimation environment:

```bash
python -m pip install --requirement requirements-lock.txt
```

The Rupprecht and Wolk reference implementation is distributed under Apache
License 2.0 through immutable
[Zenodo record 18234600](https://doi.org/10.5281/zenodo.18234600).

Run the full resource estimates from this directory:

```bash
python estimate_random_matrix.py
python estimate_f2.py
```

Regenerate the F2 wavefunction in a separate environment:

```bash
python -m pip install --requirement requirements-lock.txt
CMAKE_BUILD_PARALLEL_LEVEL=1 \
  python -m pip install --no-deps --force-reinstall \
    --requirement requirements-f2-generation-lock.txt
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  python generate_f2.py
```

This generation was validated on Ubuntu 24.04 with Python 3.12.3, GCC 13.3.0,
CMake 3.28.3, NumPy 2.3.5, SciPy 1.18.0, PySCF 2.12.0, OpenBLAS 0.3.26, and
libgomp 14.2.0. The no-dependency reinstall deliberately retains QDK 1.30.0
from the resource environment while replacing only qdk-chemistry. Install the
Rupprecht and Wolk wheel from the Zenodo archive before importing the benchmark
scripts. Other BLAS implementations are not validated for byte-level
wavefunction reproduction.

## Methodology

The eight wavefunctions in `data/input_wavefunctions.json` were extracted from
records in the [SparseCI-24 dataset](../SparseCI-24/).

The paper's F2 benchmark is regenerated from `data/structures/f2.xyz` with
def2-SVP restricted HF. The historical pipeline selected a 14-electron,
8-orbital valence space but solved the 5-alpha/5-beta MACIS sector. The script
preserves that behavior, applies the recovered basis transformation within an
exactly degenerate orbital pair, and resolves the determinant-cutoff tie
lexicographically. The final determinants are ordered as in the paper because
the resource benchmark scans prefix subsets of that ordering. The exact gauge,
electron sectors, and energies are recorded in `data/f2.json`. The historical
electron-sector mismatch is retained for artifact reproduction; it is not a
consistent neutral-F2 active-space calculation.

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
├── generate_f2.py                    # Regenerate the paper F2 wavefunction
├── generate_random_matrix.py         # Generate random determinant matrices
├── requirements-f2-generation-lock.txt # F2 wavefunction environment
├── requirements-lock.txt             # Exact resource-estimation environment
├── state_preparation_methods.py      # Resource estimators for four methods
├── data/
│   ├── f2.json                       # F2 molecular record and wavefunction
│   ├── input_wavefunctions.json      # Eight SparseCI-24 wavefunctions
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
