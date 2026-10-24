# Ruthenium catalyst XVIII FCIDUMP from von Burg et al. (2021)

Hamiltonian for intermediate XVIII in ruthenium-catalyzed conversion of carbon dioxide to methanol.

## File

| File | Active space (e, o) | `MS2` | Core energy (Ha) |
| --- | ---: | ---: | ---: |
| [XVIII-cas4-fb-64e-56o.fcidump](./XVIII-cas4-fb-64e-56o.fcidump) | 64, 56 | 0 | -7180.18924584 |

The [aka.ms short link](https://aka.ms/fcidump/XVIII-cas4-fb-64e-56o) will point to this repository.
The file is unchanged from the [Microsoft QDK v1.1.1 release asset](https://github.com/microsoft/qdk/releases/download/v1.1.1/XVIII-cas4-fb-64e-56o.fcidump).

SHA-256:

```text
bdd81245ad09f93031e1fa7f0e7f5328b7f220b9b1c27a8688b7941837584be8  XVIII-cas4-fb-64e-56o.fcidump
```

## Model and conventions

- `cas4`: orbitals optimized by CASSCF with 4 electrons in 4 spatial orbitals, CAS(4,4).
  The stored integrals use a larger selection from that orbital basis: **64 electrons in 56 spatial orbitals**, not a CASSCF(64,56) calculation.
- `fb`: ANO-RCC-VTZP for light elements and ANO-RCC-VQZP for ruthenium.
- Generated with OpenMolcas using the default Cholesky decomposition threshold of `1e-4` (not the `highCD` variant).

Standard one-based FCIDUMP indices and chemists' ordering `(ij|kl)`.
Include the core energy (`0 0 0 0`) when computing total energies.

## Citation

> V. von Burg et al., "Quantum computing enhanced computational catalysis," *Phys. Rev. Research* **3**, 033055 (2021). [doi:10.1103/PhysRevResearch.3.033055](https://doi.org/10.1103/PhysRevResearch.3.033055)

See the preprint supplement: [Section III](https://arxiv.org/pdf/2007.14460v2#page=26) for basis sets and integral generation, and [Table XVI](https://arxiv.org/pdf/2007.14460v2#page=39) for this model (listed as `XVIII-cas4-fb-64e56o`).
