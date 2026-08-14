"""State preparation methods for sparse quantum wavefunctions.

All four methods accept MSB-first bitstrings and coefficients:

  - ``gf2x`` — GF2+X elimination-based sparse isometry via qdk_chemistry.
  - ``gf2x_binary_encoding`` — GF2+X with binary encoding via qdk_chemistry.
  - ``Rupprecht2026`` — Batched isometry from Rupprecht & Wolk (2026) via
    Qualtran bloqs.
  - ``Ramacciotti2024`` — Permutation-based sparse state preparation from
    Ramacciotti et al. (2024) via Qualtran.

Also provides helpers shared across methods: ``estimate_bloq`` and
``estimate_qdk_circuit``, and the shared ``ResourceEstimate`` result type.
"""

# --------------------------------------------------------------------------------------------
# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT License. See LICENSE in the project root for license information.
# --------------------------------------------------------------------------------------------

from collections import Counter
from dataclasses import dataclass
from typing import Any

import numpy as np
from qiskit.circuit import QuantumCircuit
from qiskit.circuit.controlflow import ControlFlowOp, IfElseOp
from qiskit.compiler import transpile

# Dependency checks
try:
    from qdk_chemistry.algorithms import create
    from qdk_chemistry.data import (
        AlgorithmRef,
        Circuit,
        Configuration,
        ModelOrbitals,
        StateVectorContainer,
        Wavefunction,
    )

except ImportError:
    raise ImportError("ERROR: qdk_chemistry is required. See README.md.")

try:
    from sparse_state_preparation import SparseStatePreparation
    from sparse_state_preparation.isometry import IsometryToSubspaceViaBatching
except ImportError:
    raise ImportError(
        "ERROR: sparse_state_preparation from https://zenodo.org/records/18234600 "
        "is required. See README.md.",
    )

try:
    from qualtran import QFxp
    from qualtran.bloqs.state_preparation.sparse_state_preparation_via_rotations import (
        SparseStatePreparationViaRotations,
    )
    from qualtran.resource_counting import QECGatesCost, QubitCount, get_cost_value
except ImportError:
    raise ImportError("ERROR: qualtran is required. See README.md.")

# number of phase bits and number of fractional bits required by the reference
# implementations but unused in the benchmarking
PHASE_BITSIZE = 6
NUM_FRAC = 6

# Gate counting helpers
_BASIS_GATES = [
    "x",
    "y",
    "z",
    "cx",
    "swap",
    "cz",
    "ccx",
    "ccz",
    "id",
    "h",
    "s",
    "sdg",
    "rz",
    "cswap",
]
_CLIFFORD_GATES = {"x", "y", "z", "cx", "cz", "h", "s", "sdg", "swap"}
_TOFFOLI_GATES = {"ccx", "ccz", "cswap"}


@dataclass
class ResourceEstimateData:
    """Resource estimate for a state preparation circuit or bloq.

    Args:
        logical_qubits: Total logical qubit count.
        toffoli_count: Number of Toffoli / And-bloq gates.
        rotation_count: Number of arbitrary-angle rotation gates.
        non_clifford_count: ``toffoli_count + rotation_count``.
        clifford_count: Number of Clifford gates.
    """

    logical_qubits: int
    toffoli_count: int
    rotation_count: int
    non_clifford_count: int
    clifford_count: int


def _to_qdk_wavefunction(bitstrings: list[str], coeffs: list[complex]) -> Wavefunction:
    """Convert MSB-first bitstrings and coefficients to a QDK ``Wavefunction``.

    QDK configuration strings list ``q[0]`` first, so reverse the benchmark's
    MSB-first strings at this boundary.
    """
    n_qubits = len(bitstrings[0])
    return Wavefunction(
        StateVectorContainer(
            np.array(coeffs),
            [Configuration.from_bitstring(bitstring[::-1]) for bitstring in bitstrings],
            ModelOrbitals(n_qubits),
        )
    )


def bitstring_from_qubit_occupations(occupations: np.ndarray) -> str:
    """Convert ``q[0]``-first occupations to an MSB-first bitstring."""
    return "".join(str(int(bit)) for bit in reversed(occupations))


def bitstring_from_basis_index(index: int, n_qubits: int) -> str:
    """Convert a computational-basis index to an MSB-first bitstring."""
    return f"{index:0{n_qubits}b}"


def estimate_bloq(bloq: Any) -> ResourceEstimateData:
    """Get resource estimates directly from a qualtran Bloq.

    Args:
        bloq (Any): A qualtran ``Bloq`` supporting ``QubitCount`` and
            ``QECGatesCost`` resource-counting protocols.

    Returns:
        ResourceEstimateData: Resource estimate for the bloq.
    """
    qubit_count = get_cost_value(bloq, QubitCount())
    gate_counts = get_cost_value(bloq, QECGatesCost())
    toffoli = int(gate_counts.toffoli + gate_counts.and_bloq)
    rotation = int(gate_counts.rotation)
    return ResourceEstimateData(
        logical_qubits=int(qubit_count),
        toffoli_count=toffoli,
        rotation_count=rotation,
        non_clifford_count=toffoli + rotation,
        clifford_count=int(gate_counts.clifford),
    )


def _count_qiskit_operations(circuit: QuantumCircuit) -> Counter[str]:
    """Count operations recursively using worst-case ``if_else`` branches."""
    counts: Counter[str] = Counter()
    for instruction in circuit.data:
        operation = instruction.operation
        if not isinstance(operation, ControlFlowOp):
            counts[operation.name] += 1
            continue
        if not isinstance(operation, IfElseOp):
            raise NotImplementedError(
                f"Unsupported Qiskit control-flow operation: {operation.name}"
            )

        branch_counts = [_count_qiskit_operations(block) for block in operation.blocks]
        operation_names = {
            name for branch_count in branch_counts for name in branch_count
        }
        for name in operation_names:
            counts[name] += max(branch_count[name] for branch_count in branch_counts)
    return counts


def estimate_qdk_circuit(circuit: Circuit) -> ResourceEstimateData:
    """Estimate resources for a qdk_chemistry Circuit.

    Args:
        circuit (Circuit): A ``qdk_chemistry`` ``Circuit`` object.

    Returns:
        ResourceEstimate: Resource estimate for the transpiled circuit.
    """
    qc = circuit.get_qiskit_circuit()
    qc = transpile(qc, basis_gates=_BASIS_GATES, optimization_level=0)
    ops = _count_qiskit_operations(qc)
    toffoli_count = sum(ops.get(g, 0) for g in _TOFFOLI_GATES)
    rotation_count = ops.get("rz", 0)
    clifford_count = sum(ops.get(g, 0) for g in _CLIFFORD_GATES)
    return ResourceEstimateData(
        logical_qubits=qc.num_qubits,
        toffoli_count=toffoli_count,
        rotation_count=rotation_count,
        non_clifford_count=toffoli_count + rotation_count,
        clifford_count=clifford_count,
    )


def dense_state_prep(n_qubits: int, sv: np.ndarray) -> ResourceEstimateData:
    """Estimate dense state prep resources from a pre-built statevector.

    Args:
        n_qubits (int): Number of qubits (log2 of the statevector length).
        sv (np.ndarray): Statevector of length ``2**n_qubits``.
            May be complex-valued; imaginary parts must be negligible.
            Will be L2-normalised before use.

    Returns:
        ResourceEstimate: Resource estimate; see ``estimate_qdk_circuit``.
    """
    if np.iscomplexobj(sv):
        if np.max(np.abs(sv.imag)) > 1e-10:
            raise ValueError(
                "Dense state preparation received non-negligible imaginary "
                f"amplitudes (max |imag| = {np.max(np.abs(sv.imag)):.2e})."
            )
        sv = sv.real.copy()
    norm = np.linalg.norm(sv)
    if norm > 0:
        sv = sv / norm

    indices = np.flatnonzero(sv)
    wavefunction = _to_qdk_wavefunction(
        [bitstring_from_basis_index(int(index), n_qubits) for index in indices],
        sv[indices].tolist(),
    )
    circuit = create("state_prep", "dense_pure_state").run(wavefunction)
    return estimate_qdk_circuit(circuit)


def _subtract_dense_estimate(
    combined: ResourceEstimateData, dense: ResourceEstimateData
) -> ResourceEstimateData:
    """Recover sparse-stage gate counts from a composed QDK estimate."""
    residuals = [
        combined.toffoli_count - dense.toffoli_count,
        combined.rotation_count - dense.rotation_count,
        combined.non_clifford_count - dense.non_clifford_count,
        combined.clifford_count - dense.clifford_count,
    ]
    if min(residuals) < 0:
        raise RuntimeError(
            "Dense circuit counts exceed the composed sparse-isometry counts."
        )
    return ResourceEstimateData(
        logical_qubits=combined.logical_qubits,
        toffoli_count=residuals[0],
        rotation_count=residuals[1],
        non_clifford_count=residuals[2],
        clifford_count=residuals[3],
    )


def _estimate_qdk_sparse_isometry(
    bitstrings: list[str],
    coeffs: list[complex],
    *,
    binary_encoding: bool,
) -> tuple[ResourceEstimateData, ResourceEstimateData]:
    """Estimate the public QDK sparse-isometry plugin and its dense substage."""
    wavefunction = _to_qdk_wavefunction(bitstrings, coeffs)

    state_prep = create(
        "state_prep",
        "sparse_isometry",
        binary_encoding=binary_encoding,
        dense_state_prep=AlgorithmRef("state_prep", "dense_pure_state"),
        include_negative_controls=True,
        measurement_based_uncompute=binary_encoding,
    )
    circuit = state_prep.run(wavefunction)
    dense_circuit = state_prep.create_dense(wavefunction)

    combined_est = estimate_qdk_circuit(circuit)
    dense_est = estimate_qdk_circuit(dense_circuit)
    sparse_est = _subtract_dense_estimate(combined_est, dense_est)
    return sparse_est, dense_est


def gf2x(
    bitstrings: list[str], coeffs: list[complex]
) -> tuple[ResourceEstimateData, ResourceEstimateData]:
    """Run GF2+X sparse isometry via qdk_chemistry.

    Args:
        bitstrings (list[str]): Computational-basis bitstrings representing
            the non-zero amplitudes of the target state.
        coeffs (list[complex]): Expansion coefficients corresponding to
            each bitstring.

    Returns:
        tuple[ResourceEstimate, ResourceEstimate]: A pair
            ``(sparse_est, dense_est)``.
    """
    return _estimate_qdk_sparse_isometry(bitstrings, coeffs, binary_encoding=False)


def gf2x_binary_encoding(
    bitstrings: list[str], coeffs: list[complex]
) -> tuple[ResourceEstimateData, ResourceEstimateData]:
    """Run GF2+X with binary encoding via qdk_chemistry.

    Uses the QDK sparse-isometry plugin with batched Toffoli-based binary
    encoding. QDK falls back to standard GF2+X when encoding offers no qubit
    advantage.

    Args:
        bitstrings (list[str]): Computational-basis bitstrings representing
            the non-zero amplitudes of the target state.
        coeffs (list[complex]): Expansion coefficients corresponding to
            each bitstring.

    Returns:
        tuple[ResourceEstimate, ResourceEstimate]: A pair
            ``(sparse_est, dense_est)``.
    """
    return _estimate_qdk_sparse_isometry(
        bitstrings,
        coeffs,
        binary_encoding=True,
    )


def Rupprecht2026(
    bitstrings: list[str],
    coeffs: list[complex],
    phase_bitsize: int = PHASE_BITSIZE,
    num_frac: int = NUM_FRAC,
) -> tuple[ResourceEstimateData, ResourceEstimateData]:
    """Rupprecht & Wolk 2026 batched isometry method.

    Estimates the isometry cost via qualtran ``IsometryToSubspaceViaBatching``
    and the dense state-preparation cost via ``dense_state_prep``.

    Args:
        bitstrings (list[str]): Computational-basis bitstrings representing
            the non-zero amplitudes of the target state.
        coeffs (list[complex]): Expansion coefficients corresponding to
            each bitstring.
        phase_bitsize (int): Total bit-size of the fixed-point phase register
            used inside ``SparseStatePreparation``.  Defaults to
            ``PHASE_BITSIZE`` (6).
        num_frac (int): Number of fractional bits in the fixed-point phase
            register.  Defaults to ``NUM_FRAC`` (6).

    Returns:
        tuple[ResourceEstimate, ResourceEstimate]: A pair
            ``(sparse_est, dense_est)`` where ``sparse_est`` comes from
            ``estimate_bloq`` and ``dense_est`` from ``dense_state_prep``.
    """
    states = np.array([[b == "1" for b in bs] for bs in bitstrings], dtype=bool)
    coeffs_arr = np.array(coeffs)

    sparse_prep = SparseStatePreparation(
        states,
        coeffs_arr,
        QFxp(bitsize=phase_bitsize, num_frac=num_frac),
        uncompute_in_isometry=True,
    )
    dense_bitsize = sparse_prep.isometry.subspace_bitsize
    dense_coeffs = sparse_prep._permuted_coefficients
    # Build the dense statevector directly — the Rupprecht dense register spans
    # all 2^dense_bitsize computational basis states (arbitrary electron counts),
    # which cannot be represented as a QDKWavefunction.  Bypass it entirely.
    sv = np.zeros(2**dense_bitsize, dtype=complex)
    for i, c in enumerate(dense_coeffs):
        sv[i] = c

    isometry_bloq = IsometryToSubspaceViaBatching(
        states=states,
        signs=sparse_prep.isometry.signs,
    )
    sparse_est = estimate_bloq(isometry_bloq)
    dense_est = dense_state_prep(dense_bitsize, sv)

    return sparse_est, dense_est


def _bitstrings_to_coefficient_map(
    bitstrings: list[str], coeffs: list[complex]
) -> dict[int, complex]:
    """Convert bitstrings to a format compatible with Ramacciotti2024.

    Args:
        bitstrings (list[str]): Computational-basis bitstrings.
        coeffs (list[complex]): Complex expansion coefficients corresponding
            to each bitstring.

    Returns:
        dict[int, complex]: Mapping from integer computational-basis index
            (``int(bitstring, 2)``) to its complex coefficient.
    """
    coef_map: dict[int, complex] = {}
    for bs, cf in zip(bitstrings, coeffs):
        idx = int(bs, 2)
        if idx in coef_map:
            raise ValueError(f"Duplicate bitstring encountered: {bs!r}")
        coef_map[idx] = cf
    return coef_map


def Ramacciotti2024(
    bitstrings: list[str], coeffs: list[complex], phase_bitsize: int = PHASE_BITSIZE
) -> tuple[ResourceEstimateData, ResourceEstimateData]:
    """Ramacciotti et al. 2024 permutation-based sparse state preparation.

    Estimates the basis-permutation isometry cost via qualtran
    ``SparseStatePreparationViaRotations`` and the dense state-preparation
    cost via ``dense_state_prep``.

    Args:
        bitstrings (list[str]): Computational-basis bitstrings representing
            the non-zero amplitudes of the target state.
        coeffs (list[complex]): Expansion coefficients corresponding to
            each bitstring.
        phase_bitsize (int): Bit-size of the phase register used inside
            ``SparseStatePreparationViaRotations``.  Defaults to
            ``PHASE_BITSIZE`` (6).

    Returns:
        tuple[ResourceEstimate, ResourceEstimate]: A pair
            ``(sparse_est, dense_est)`` where ``sparse_est`` comes from
            ``estimate_bloq`` and ``dense_est`` from ``dense_state_prep``.
    """
    num_qubits = len(bitstrings[0])

    coef_map = _bitstrings_to_coefficient_map(bitstrings, coeffs)

    sparse_prep = SparseStatePreparationViaRotations.from_coefficient_map(
        N=2**num_qubits, coeff_map=coef_map, phase_bitsize=phase_bitsize
    )
    isometry_bloq = sparse_prep._basis_permutation_bloq
    sparse_est = estimate_bloq(isometry_bloq)

    dense_bitsize = sparse_prep.dense_bitsize
    dense_coeffs = sparse_prep.nonzero_coeffs
    sv = np.zeros(2**dense_bitsize, dtype=complex)
    for i, c in enumerate(dense_coeffs):
        sv[i] = c
    dense_est = dense_state_prep(dense_bitsize, sv)

    return sparse_est, dense_est


def combine_estimates(
    sparse: ResourceEstimateData, dense: ResourceEstimateData
) -> ResourceEstimateData:
    """Combine sparse and dense resource estimates into a single estimate.

    Args:
        sparse: Resource estimate for the sparse isometry circuit.
        dense: Resource estimate for the dense state preparation circuit.

    Returns:
        Merged estimate where qubit count is the max and gate counts are summed.
    """
    return ResourceEstimateData(
        logical_qubits=max(sparse.logical_qubits, dense.logical_qubits),
        toffoli_count=sparse.toffoli_count + dense.toffoli_count,
        rotation_count=sparse.rotation_count + dense.rotation_count,
        non_clifford_count=sparse.non_clifford_count + dense.non_clifford_count,
        clifford_count=sparse.clifford_count + dense.clifford_count,
    )


@dataclass
class BenchmarkResult:
    """Resource estimates for one method on one wavefunction/system.

    Shared result type used by both the F2 and random-matrix benchmarks.

    Args:
        method: Method name (``"gf2x"``, ``"gf2x_binary_encoding"``,
            ``"Rupprecht2026"``, ``"Ramacciotti2024"``).
        num_qubits: Number of qubits (bitstring length).
        num_dets: Number of configurations (non-zero amplitudes).
        sparse: Resource estimate for the sparse isometry circuit.
        dense: Resource estimate for the dense state preparation circuit.
        source: Data source, e.g. ``"random"`` or ``"chemical"``.
        molecule: Molecule name when ``source == "chemical"``, else ``None``.
    """

    method: str
    num_qubits: int
    num_dets: int
    sparse: ResourceEstimateData
    dense: ResourceEstimateData
    source: str = "random"
    molecule: str | None = None

    @property
    def combined(self) -> ResourceEstimateData:
        """Combined sparse + dense resource estimate."""
        return combine_estimates(self.sparse, self.dense)
