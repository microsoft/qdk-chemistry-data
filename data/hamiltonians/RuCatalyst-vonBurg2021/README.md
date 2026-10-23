# Ruthenium catalyst XVIII FCIDUMP from von Burg et al. (2021)

Integral file for intermediate XVIII of the ruthenium-catalyzed conversion of carbon dioxide to methanol studied in [*Quantum computing enhanced computational catalysis*](https://doi.org/10.1103/PhysRevResearch.3.033055).

## File and provenance

| File | Active space (e, o) | `MS2` | Constant energy (Ha) |
| --- | ---: | ---: | ---: |
| [XVIII-cas4-fb-64e-56o.fcidump](./XVIII-cas4-fb-64e-56o.fcidump) | 64, 56 | 0 | -7180.18924584 |

Downloaded from [`aka.ms/fcidump/XVIII-cas4-fb-64e-56o`](https://aka.ms/fcidump/XVIII-cas4-fb-64e-56o), which resolves to the [Microsoft QDK v1.1.1 release asset](https://github.com/microsoft/qdk/releases/download/v1.1.1/XVIII-cas4-fb-64e-56o.fcidump).
The filename and contents are unmodified.
The supplement lists the corresponding model as `XVIII-cas4-fb-64e56o` (without the extra hyphen in the release filename).

SHA-256 (63,800,318 bytes):

```text
bdd81245ad09f93031e1fa7f0e7f5328b7f220b9b1c27a8688b7941837584be8  XVIII-cas4-fb-64e-56o.fcidump
```

## Model and conventions

- `XVIII` identifies the catalyst intermediate, not a transition state.
- `cas4` identifies the CAS(4,4)SCF calculation used to generate the molecular orbitals, not the active space stored in this file.
  The Hamiltonian contains 64 electrons in 56 spatial orbitals.
- `fb` denotes the full atomic orbital basis: ANO-RCC-VTZP for light elements and ANO-RCC-VQZP for ruthenium, with integrals generated using OpenMolcas.
- The absence of `highCD` denotes the default integral-generation Cholesky decomposition threshold of `1e-4`, not the `1e-8` high-accuracy variant.
  This is separate from subsequent Hamiltonian truncation in quantum resource estimation.
- The large integral active space adds the six phenyl groups' 18 pi and 18 pi-antibonding orbitals to the intermediate active space.

The file uses standard one-based FCIDUMP indices and chemists' two-electron ordering `(ij|kl)`; do not apply a LIQUiD reordering.
The header has `NORB=56`, `NELEC=64`, and `MS2=0` (32 alpha and 32 beta electrons).
Include the constant `0 0 0 0` record when calculating absolute energies; it is the scalar offset of the active-space Hamiltonian, not just nuclear repulsion.

## Citation

When using this file, please cite:

> V. von Burg et al., "Quantum computing enhanced computational catalysis," *Phys. Rev. Research* **3**, 033055 (2021). [doi:10.1103/PhysRevResearch.3.033055](https://doi.org/10.1103/PhysRevResearch.3.033055)

See also [arXiv:2007.14460v2](https://arxiv.org/abs/2007.14460v2), which includes the supplementary information: Section III and Table V (basis sets and orbital generation), Section III C (active-space selection), and Section IV, Table XVI (this integral model).
