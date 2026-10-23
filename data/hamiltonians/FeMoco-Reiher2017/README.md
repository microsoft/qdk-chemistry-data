# FeMoco FCIDUMPs from Reiher et al. (2017)

Original integral files for the two FeMoco models in [*Elucidating reaction mechanisms on quantum computers*](https://arxiv.org/abs/1605.03590) (arXiv, May 2016; PNAS, July 2017).

We obtained these files directly from Markus Reiher and publish them with his permission.
Their contents are unmodified.
The filenames keep the descriptors from the names used when the files were first shared: spin, charge, orbital-generation calculation, and target active space.
We replaced the `FCIDUMP_` prefix with `FeMoco_` and added the `.fcidump` extension.
Some redistributed copies replaced the original names with names based on Reiher's surname, losing that information.

## Correction for the `nitrogenase-54e-54o` download

[`aka.ms/fcidump/nitrogenase-54e-54o`](https://aka.ms/fcidump/nitrogenase-54e-54o) previously pointed to [`nitrogenase-54e-54o.fcidump`](https://github.com/microsoft/qdk/releases/download/v1.1.1/nitrogenase-54e-54o.fcidump), an asset on the [Microsoft QDK v1.1.1 release](https://github.com/microsoft/qdk/releases/tag/v1.1.1).
That file uses LIQUiD two-electron index ordering under an FCIDUMP header.
Read as a standard FCIDUMP, it gives the wrong Hamiltonian.

The link now points to this repository.
Other corrupt or incorrectly converted copies remain in circulation; use these originals.

## Files

| File | Spin `S` | Charge | Active space (e, o) | Core energy (Ha) |
| --- | ---: | ---: | ---: | ---: |
| [FeMoco_S0_charge+3_orbsRAS24.16_targetCAS54.54_test3.fcidump](./FeMoco_S0_charge+3_orbsRAS24.16_targetCAS54.54_test3.fcidump) | 0 | +3 | 54, 54 | -13212.970326 |
| [FeMoco_S1.5_charge0_orbsRAS21.12_targetCAS65.57_test7_b.fcidump](./FeMoco_S1.5_charge0_orbsRAS21.12_targetCAS65.57_test7_b.fcidump) | 3/2 | 0 | 65, 57 | -13164.048069 |

- Spin `S`: total spin of the calculation that generated the orbitals (filename `S`).
  The FCIDUMP header has `MS2` = 2S, so the 65-electron file has 34 alpha and 31 beta electrons.
- Active space: electrons and orbitals in the Hamiltonian stored in the file (filename `targetCAS`).
- Core energy: the constant `0 0 0 0` record (nuclear repulsion plus the energy of the frozen orbitals).

The files use standard one-based FCIDUMP indices and chemists' ordering `(ij|kl)`; do not apply a LIQUiD reordering.
Add the core energy to obtain absolute energies.

SHA-256:

```text
84f3a34d61cb1552797cc2574159ec4b6efa641c380631fb4f66be654a34940e  FeMoco_S0_charge+3_orbsRAS24.16_targetCAS54.54_test3.fcidump
c3759d0bdc71801dc54f5f6536dd8b2d5fec16c6e0e57ae3dc22e4ef7a996073  FeMoco_S1.5_charge0_orbsRAS21.12_targetCAS65.57_test7_b.fcidump
```

## How Reiher et al. generated the files

Reiher et al. optimized two FeMoco structures with broken-symmetry DFT, ran a small CASSCF calculation on each, selected the target active space from its natural orbitals, and computed the integrals with Molcas.
The two files therefore use different structures and orbitals.
For details, see Appendix J of the arXiv preprint or Section 11 of the PNAS supporting information.

The spin and charge states were chosen to test quantum algorithms and do not necessarily match the resting state of nitrogenase.
Treat these files as benchmark Hamiltonians, not as models of the FeMoco ground state.

## Citation

When using these files, please cite:

> M. Reiher, N. Wiebe, K. M. Svore, D. Wecker, and M. Troyer, "Elucidating reaction mechanisms on quantum computers," *PNAS* **114**, 7555-7560 (2017). [doi:10.1073/pnas.1619152114](https://doi.org/10.1073/pnas.1619152114)

The [May 2016 arXiv preprint](https://arxiv.org/abs/1605.03590v1) describes both models in Appendix J, "Computational Methodology."
See also the [PNAS supporting information](https://pmc.ncbi.nlm.nih.gov/articles/instance/5530650/bin/pnas.1619152114.sapp.pdf).

Related later work: [von Burg et al., *Phys. Rev. Research* **3**, 033055 (2021)](https://doi.org/10.1103/PhysRevResearch.3.033055).
