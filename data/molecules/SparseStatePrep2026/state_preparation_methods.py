"""State preparation methods for sparse quantum wavefunctions.

All four methods accept MSB-first bitstrings and coefficients:

  - ``gf2x`` — GF2+X elimination-based sparse isometry via qdk_chemistry.
  - ``gf2x_binary_encoding`` — GF2+X with binary encoding via qdk_chemistry.
  - ``Rupprecht2026`` — Batched isometry from Rupprecht & Wolk (2026) via
    Qualtran bloqs.
  - ``Ramacciotti2024`` — Permutation-based sparse state preparation from
    Ramacciotti et al. (2024) via Qualtran.

Also provides helpers shared across methods: ``estimate_bloq`` and
``estimate_qdk_circuit``, and the shared ``ResourceEstimateData`` result type.
"""

# --------------------------------------------------------------------------------------------
# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT License. See LICENSE in the project root for license information.
# --------------------------------------------------------------------------------------------

import platform
from collections import Counter
from collections.abc import Callable
from dataclasses import dataclass
from functools import cache
from importlib.metadata import version as distribution_version
from typing import Any

import numpy as np
from qdk import TargetProfile
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
    from qdk_chemistry.utils.qsharp import (
        create_qsharp_context,
        use_qsharp_context,
    )

except ImportError as exc:
    raise ImportError("ERROR: qdk_chemistry is required. See README.md.") from exc
try:
    from sparse_state_preparation import SparseStatePreparation
    from sparse_state_preparation.isometry import IsometryToSubspaceViaBatching
except ImportError as exc:
    raise ImportError(
        "ERROR: sparse_state_preparation from https://zenodo.org/records/18234600 "
        "is required. See README.md.",
    ) from exc

try:
    from qualtran import Bloq, QFxp
    from qualtran.bloqs.mcmt import MultiTargetCNOT
    from qualtran.bloqs.state_preparation.sparse_state_preparation_via_rotations import (
        SparseStatePreparationViaRotations,
    )
    from qualtran.resource_counting import (
        GateCounts,
        QECGatesCost,
        QubitCount,
        get_bloq_callee_counts,
        get_cost_value,
    )
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


class _DecomposedMultiTargetCNOTGatesCost(QECGatesCost):
    """``QECGatesCost`` that resolves ``MultiTargetCNOT`` into its CNOT ladder.

    Qualtran's default ``QECGatesCost`` scores ``MultiTargetCNOT(k)`` as a single
    Clifford regardless of ``k``, whereas the QDK methods count every primitive 
    ``cx`` left after transpilation.
    """

    def compute(
        self, bloq: Bloq, get_callee_cost: Callable[[Bloq], GateCounts]
    ) -> GateCounts:
        if not isinstance(bloq, MultiTargetCNOT):
            return super().compute(bloq, get_callee_cost)
        totals = GateCounts()
        for callee, times in get_bloq_callee_counts(bloq, ignore_decomp_failure=False):
            totals += times * get_callee_cost(callee)
        return totals


@cache
def _adaptive_qsharp_context() -> Any:
    """Q# context allowing the mid-circuit measurement that `AND` uncompute uses."""
    return create_qsharp_context(target_profile=TargetProfile.Adaptive_RIF)


def benchmark_environment() -> dict[str, Any]:
    """Return the resolved package versions used by the benchmark."""
    return {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "qdk-chemistry": {
            "version": distribution_version("qdk-chemistry"),
        },
        "qdk": distribution_version("qdk"),
        "qsharp": distribution_version("qsharp"),
        "qualtran": distribution_version("qualtran"),
        "qiskit": distribution_version("qiskit"),
        "numpy": distribution_version("numpy"),
        "matplotlib": distribution_version("matplotlib"),
        "sparse-state-preparation": {
            "version": distribution_version("sparse-state-preparation"),
        },
    }


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
        return ResourceEstimateData(
            logical_qubits=max(self.sparse.logical_qubits, self.dense.logical_qubits),
            toffoli_count=self.sparse.toffoli_count + self.dense.toffoli_count,
            rotation_count=self.sparse.rotation_count + self.dense.rotation_count,
            non_clifford_count=self.sparse.non_clifford_count
            + self.dense.non_clifford_count,
            clifford_count=self.sparse.clifford_count + self.dense.clifford_count,
        )


def _to_qdk_wavefunction(bitstrings: list[str], coeffs: list[complex]) -> Wavefunction:
    """Convert MSB-first bitstrings and coefficients to a QDK ``Wavefunction``."""
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


def estimate_bloq(bloq: Any) -> ResourceEstimateData:
    """Get resource estimates directly from a qualtran Bloq.

    ``MultiTargetCNOT`` is expanded into its CNOT decomposition so that Clifford
    counts are comparable with the transpiled QDK circuits; see
    :class:`_DecomposedMultiTargetCNOTGatesCost`.

    Args:
        bloq (Any): A qualtran ``Bloq`` supporting ``QubitCount`` and
            ``QECGatesCost`` resource-counting protocols.

    Returns:
        ResourceEstimateData: Resource estimate for the bloq.
    """
    qubit_count = get_cost_value(bloq, QubitCount())
    gate_counts = get_cost_value(bloq, _DecomposedMultiTargetCNOTGatesCost())
    toffoli = int(gate_counts.toffoli + gate_counts.and_bloq + gate_counts.cswap)
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
        ResourceEstimateData: Resource estimate for the transpiled circuit.
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
        ResourceEstimateData: Resource estimate; see ``estimate_qdk_circuit``.
    """
    if np.iscomplexobj(sv):
        if np.max(np.abs(sv.imag)) > 1e-10:
            raise ValueError(
                "Dense state preparation received non-negligible imaginary "
                f"amplitudes (max |imag| = {np.max(np.abs(sv.imag)):.2e})."
            )
        sv = sv.real.copy()
    norm = np.linalg.norm(sv)
    if norm == 0:
        raise ValueError("Dense state preparation received a zero-norm statevector.")
    sv = sv / norm
    indices = np.flatnonzero(sv)
    wavefunction = _to_qdk_wavefunction(
        [f"{int(index):0{n_qubits}b}" for index in indices],
        sv[indices].tolist(),
    )
    circuit = create("state_prep", "dense_pure_state").run(wavefunction)
    return estimate_qdk_circuit(circuit)


def _estimate_qdk_sparse_isometry(
    bitstrings: list[str],
    coeffs: list[complex],
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
        measurement_based_uncompute=True,
    )
    with use_qsharp_context(_adaptive_qsharp_context()):
        circuit = state_prep.run(wavefunction)
        combined_est = estimate_qdk_circuit(circuit)
        dense_circuit = state_prep.create_dense(wavefunction)
        dense_est = estimate_qdk_circuit(dense_circuit)
    residuals = [
        combined_est.toffoli_count - dense_est.toffoli_count,
        combined_est.rotation_count - dense_est.rotation_count,
        combined_est.non_clifford_count - dense_est.non_clifford_count,
        combined_est.clifford_count - dense_est.clifford_count,
    ]
    if min(residuals) < 0:
        raise RuntimeError(
            "Dense circuit counts exceed the composed sparse-isometry counts."
        )
    sparse_est = ResourceEstimateData(
        logical_qubits=combined_est.logical_qubits,
        toffoli_count=residuals[0],
        rotation_count=residuals[1],
        non_clifford_count=residuals[2],
        clifford_count=residuals[3],
    )
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
        tuple[ResourceEstimateData, ResourceEstimateData]: A pair
            ``(sparse_est, dense_est)``.
    """
    return _estimate_qdk_sparse_isometry(
        bitstrings,
        coeffs,
        binary_encoding=False,
    )


def gf2x_binary_encoding(
    bitstrings: list[str],
    coeffs: list[complex],
) -> tuple[ResourceEstimateData, ResourceEstimateData]:
    """Run GF2+X with binary encoding via qdk_chemistry.

    Uses the QDK sparse-isometry plugin with batched Toffoli-based binary
    encoding. Measured uncomputation explicitly selects the Adaptive Q# profile
    for optimal resource estimates.

    Args:
        bitstrings (list[str]): Computational-basis bitstrings representing
            the non-zero amplitudes of the target state.
        coeffs (list[complex]): Expansion coefficients corresponding to
            each bitstring.

    Returns:
        tuple[ResourceEstimateData, ResourceEstimateData]: A pair
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
        tuple[ResourceEstimateData, ResourceEstimateData]: A pair
            ``(sparse_est, dense_est)`` where ``sparse_est`` comes from
            ``estimate_bloq`` and ``dense_est`` from ``dense_state_prep``.
    """
    states = np.array([[b == "1" for b in bs] for bs in bitstrings], dtype=bool)
    coeffs_arr = np.array(coeffs)

    sparse_prep = SparseStatePreparation(
        states,
        coeffs_arr,
        QFxp(bitsize=phase_bitsize, num_frac=num_frac),
        uncompute_in_isometry=False,
    )
    dense_bitsize = sparse_prep.isometry.subspace_bitsize
    dense_coeffs = getattr(sparse_prep, "_permuted_coefficients", None)
    if dense_coeffs is None:
        raise AttributeError(
            "sparse_state_preparation API changed: missing '_permuted_coefficients'"
        )
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
        tuple[ResourceEstimateData, ResourceEstimateData]: A pair
            ``(sparse_est, dense_est)`` where ``sparse_est`` comes from
            ``estimate_bloq`` and ``dense_est`` from ``dense_state_prep``.
    """
    num_qubits = len(bitstrings[0])

    coef_map: dict[int, complex] = {}
    for bs, cf in zip(bitstrings, coeffs):
        idx = int(bs, 2)
        if idx in coef_map:
            raise ValueError(f"Duplicate bitstring encountered: {bs!r}")
        coef_map[idx] = cf

    sparse_prep = SparseStatePreparationViaRotations.from_coefficient_map(
        N=2**num_qubits, coeff_map=coef_map, phase_bitsize=phase_bitsize
    )
    isometry_bloq = getattr(sparse_prep, "_basis_permutation_bloq", None)
    if isometry_bloq is None:
        raise AttributeError(
            "Qualtran API changed: missing '_basis_permutation_bloq' on SparseStatePreparationViaRotations"
        )
    sparse_est = estimate_bloq(isometry_bloq)

    dense_bitsize = sparse_prep.dense_bitsize
    dense_coeffs = sparse_prep.nonzero_coeffs
    sv = np.zeros(2**dense_bitsize, dtype=complex)
    for i, c in enumerate(dense_coeffs):
        sv[i] = c
    dense_est = dense_state_prep(dense_bitsize, sv)

    return sparse_est, dense_est
