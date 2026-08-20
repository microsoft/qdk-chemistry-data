"""Generate a determinant matrix mimicking chemical molecule orbitals.

Starting from the Hartree-Fock (HF) reference state, this module generates
excited Slater determinants by randomly selecting excitations (single, double,
etc.). Each determinant is represented as a binary occupation vector over
spin-orbitals.
"""

# --------------------------------------------------------------------------------------------
# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT License. See LICENSE in the project root for license information.
# --------------------------------------------------------------------------------------------

from typing import Optional
import numpy as np
from math import comb


def _hf_determinant(n_alpha: int, n_beta: int, n_orbitals: int) -> np.ndarray:
    """Build the Hartree-Fock reference determinant.

    Layout: [alpha_0, alpha_1, ..., alpha_{n_orbitals-1},
             beta_0,  beta_1,  ..., beta_{n_orbitals-1}]

    The HF state fills the lowest orbitals for each spin channel.

    Args:
        n_alpha (int): Number of alpha electrons.
        n_beta (int): Number of beta electrons.
        n_orbitals (int): Number of spatial orbitals.

    Returns:
        np.ndarray: Binary occupation vector of shape ``(2 * n_orbitals,)``
            with dtype ``int8``, where the first ``n_orbitals`` entries are
            alpha occupancies and the last ``n_orbitals`` entries are beta
            occupancies.
    """
    det = np.zeros(2 * n_orbitals, dtype=np.int8)
    det[:n_alpha] = 1  # alpha electrons in lowest orbitals
    det[n_orbitals : n_orbitals + n_beta] = 1  # beta electrons in lowest orbitals
    return det


def _random_excitation(
    det: np.ndarray,
    n_orbitals: int,
    rng: np.random.Generator,
    max_excitation_order: Optional[int] = None,
) -> Optional[np.ndarray]:
    """Generate a single random excitation from a reference determinant.

    Independently applies excitations in the alpha and beta spin channels.
    The excitation order is randomly chosen between 0 and max_excitation_order
    (inclusive); 0 means no excitation in that channel.
    Args:
        det (np.ndarray): Binary occupation vector of shape ``(2 * n_orbitals,)``
            representing the reference determinant (alpha block then beta block).
        n_orbitals (int): Number of spatial orbitals.
        rng (np.random.Generator): NumPy random generator instance used for
            all random choices.
        max_excitation_order (Optional[int]): Maximum number of
            occupied-to-virtual replacements per spin channel.  ``None`` means
            no upper bound beyond the physical limit.

    Returns:
        Optional[np.ndarray]: A new binary occupation vector of the same shape
            as ``det`` representing the excited determinant, or ``None`` if the
            randomly chosen excitation is the identity (no change).
    """
    new_det = det.copy()

    # Process alpha and beta channels independently
    for channel_start in (0, n_orbitals):
        channel = det[channel_start : channel_start + n_orbitals]
        occupied = np.where(channel == 1)[0]
        virtual = np.where(channel == 0)[0]

        if len(occupied) == 0 or len(virtual) == 0:
            continue

        max_exc = min(len(occupied), len(virtual))
        if max_excitation_order is not None:
            max_exc = min(max_exc, max_excitation_order)

        # Zero means no excitation in this channel.
        order = rng.integers(0, max_exc + 1)
        if order == 0:
            continue

        occ_indices = rng.choice(occupied, size=order, replace=False)
        vir_indices = rng.choice(virtual, size=order, replace=False)

        new_det[channel_start + occ_indices] = 0
        new_det[channel_start + vir_indices] = 1

    if np.array_equal(new_det, det):
        return None

    return new_det


def generate_determinants_matrix(
    n_electrons: int,
    n_orbitals: int,
    n_dets: int,
    seed: int = 0,
    max_excitation_order: Optional[int] = None,
    n_alpha: Optional[int] = None,
    n_beta: Optional[int] = None,
    include_hf: bool = True,
) -> np.ndarray:
    """Generate a matrix of Slater determinants by random excitation from HF.

    Each row is a determinant represented as a binary occupation vector over
    2*n_orbitals spin-orbitals (alpha block followed by beta block).

    This function does NOT enumerate the full determinant space; instead it
    randomly samples excitations, making it suitable for large systems where
    C(n_orbitals, n_electrons) is intractable.

    Args:
        n_electrons: Total number of electrons.
        n_orbitals:  Number of spatial orbitals (spin-orbitals = 2 * n_orbitals).
        n_dets:      Number of determinants to generate, including HF if requested.
        seed:        Random seed for reproducibility.
        max_excitation_order: Maximum excitation rank per spin channel.
                              None means up to the maximum possible.
        n_alpha:     Number of alpha electrons. Defaults to n_electrons // 2.
        n_beta:      Number of beta electrons. Defaults to n_electrons - n_alpha.
        include_hf:  Whether to always include the HF determinant as the first row.

    Returns:
        np.ndarray of shape (n_dets, 2 * n_orbitals) with entries 0 or 1.

    Raises:
        ValueError: If inputs are inconsistent.
    """
    if n_alpha is None:
        n_alpha = n_electrons // 2
    if n_beta is None:
        n_beta = n_electrons - n_alpha

    if n_alpha + n_beta != n_electrons:
        raise ValueError(
            f"n_alpha ({n_alpha}) + n_beta ({n_beta}) != n_electrons ({n_electrons})"
        )
    if n_alpha > n_orbitals or n_beta > n_orbitals:
        raise ValueError(
            f"Cannot place {n_alpha} alpha or {n_beta} beta electrons in "
            f"{n_orbitals} orbitals"
        )
    if n_dets < 1:
        raise ValueError("n_dets must be at least 1")

    # Upper bound on the number of unique determinants
    max_possible = comb(n_orbitals, n_alpha) * comb(n_orbitals, n_beta)
    if n_dets > max_possible:
        raise ValueError(
            f"Requested {n_dets} determinants but the total space has only "
            f"{max_possible}"
        )

    rng = np.random.default_rng(seed)
    hf = _hf_determinant(n_alpha, n_beta, n_orbitals)

    # Use bytes keys for fast duplicate checking (faster than tuple for large vectors)
    seen: set[bytes] = set()
    dets: list[np.ndarray] = []

    if include_hf:
        dets.append(hf)
        seen.add(hf.tobytes())

    # We attempt random excitations from the HF reference.
    # If the space is very sparse relative to n_dets, this may take many tries,
    # but we cap attempts to avoid infinite loops.
    max_attempts = n_dets * 200
    attempts = 0

    while len(dets) < n_dets and attempts < max_attempts:
        attempts += 1
        new_det = _random_excitation(hf, n_orbitals, rng, max_excitation_order)
        if new_det is None:
            continue
        key = new_det.tobytes()
        if key not in seen:
            seen.add(key)
            dets.append(new_det)

    if len(dets) < n_dets:
        raise RuntimeError(
            f"Only generated {len(dets)}/{n_dets} unique determinants after "
            f"{max_attempts} attempts. Try increasing max_excitation_order or "
            f"reducing n_dets."
        )

    return np.array(dets, dtype=np.int8)
