"""
Domain models, typed data structures, and enumeration types for NLE Predictor.
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional
import numpy as np


class ModelType(str, Enum):
    """Supported asymmetric catalysis kinetic/thermodynamic models."""
    KAGAN_RESERVOIR = "kagan_ml2"
    KAGAN_ACTIVE_DIMER = "kagan_active"
    NOYORI_AGGREGATION = "noyori_daib"


@dataclass
class XYZStructure:
    """Represents a validated 3D molecular structure."""
    atoms: List[str]
    coordinates: np.ndarray  # Shape (N, 3), float64
    comment: str = ""
    atom_count: int = 0
    formula: str = ""
    molecular_weight: float = 0.0

    def __post_init__(self):
        if self.atom_count == 0:
            self.atom_count = len(self.atoms)
        if self.coordinates.ndim != 2 or self.coordinates.shape[1] != 3:
            raise ValueError(f"Coordinates must have shape (N, 3), got {self.coordinates.shape}")
        if len(self.atoms) != len(self.coordinates):
            raise ValueError(
                f"Mismatch: {len(self.atoms)} atom symbols vs {len(self.coordinates)} coordinate rows."
            )

    @property
    def centroid(self) -> np.ndarray:
        """Geometric centroid (mean coordinate vector)."""
        return np.mean(self.coordinates, axis=0)


@dataclass
class Conformer:
    """Individual conformer in an ensemble with energetic and statistical metrics."""
    index: int
    path: str
    energy_hartree: float
    energy_kcal: float = 0.0
    relative_energy_kcal: float = 0.0
    boltzmann_weight: float = 0.0

    def to_dict(self) -> Dict[str, Any]:
        return {
            "index": self.index,
            "path": self.path,
            "energy_hartree": self.energy_hartree,
            "energy_kcal": self.energy_kcal,
            "relative_energy_kcal": self.relative_energy_kcal,
            "boltzmann_weight": self.boltzmann_weight,
        }


@dataclass
class EnsembleResult:
    """Aggregated Boltzmann-weighted ensemble thermodynamics for a catalyst state."""
    label: str  # e.g., "RR" or "RS"
    conformers: List[Conformer]
    ensemble_gibbs_hartree: float
    ensemble_gibbs_kcal: float
    lowest_gibbs_hartree: float
    conformer_count: int
    conformational_entropy_cal_mol_k: float = 0.0

    def to_dict(self) -> Dict[str, Any]:
        return {
            "label": self.label,
            "conformer_count": self.conformer_count,
            "ensemble_gibbs_hartree": self.ensemble_gibbs_hartree,
            "ensemble_gibbs_kcal": self.ensemble_gibbs_kcal,
            "lowest_gibbs_hartree": self.lowest_gibbs_hartree,
            "conformational_entropy_cal_mol_k": self.conformational_entropy_cal_mol_k,
            "conformers": [c.to_dict() for c in self.conformers],
        }


@dataclass
class ThermodynamicParameters:
    """Thermodynamic equilibrium parameters governing enantiomer dimerization."""
    g_rr_hartree: float
    g_rs_hartree: float
    delta_g_hartree: float
    delta_g_kcal: float
    delta_g_kj: float
    k_homo: float
    k_hetero: float
    k_ratio: float  # K_hetero / K_homo
    temperature_k: float

    def to_dict(self) -> Dict[str, Any]:
        return {
            "g_rr_hartree": self.g_rr_hartree,
            "g_rs_hartree": self.g_rs_hartree,
            "delta_g_hartree": self.delta_g_hartree,
            "delta_g_kcal_per_mol": self.delta_g_kcal,
            "delta_g_kj_per_mol": self.delta_g_kj,
            "k_homo": self.k_homo,
            "k_hetero": self.k_hetero,
            "k_hetero_over_k_homo": self.k_ratio,
            "temperature_k": self.temperature_k,
        }


@dataclass
class NLESimulationResult:
    """Output of an NLE trajectory simulation across catalyst ee values."""
    model_name: str
    ee_catalyst: np.ndarray  # In range [0, 1]
    ee_product: np.ndarray   # In range [0, 1]
    monomers_r: np.ndarray
    monomers_s: np.ndarray
    classification: str       # "positive NLE", "negative NLE", "near-linear"
    product_at_50: float     # In range [0, 1]
    max_amplification: float
    max_deviation_from_linearity: float
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "model_name": self.model_name,
            "classification": self.classification,
            "product_ee_at_50_percent_catalyst": float(self.product_at_50 * 100),
            "max_amplification": float(self.max_amplification),
            "max_deviation_percent": float(self.max_deviation_from_linearity * 100),
            "metadata": self.metadata,
        }


@dataclass
class FittingResult:
    """Result of non-linear parameter estimation from experimental data."""
    delta_g_kcal: float
    delta_g_stderr: float
    k_ratio: float
    k_ratio_stderr: float
    r_squared: float
    rmse: float
    residuals: np.ndarray
    converged: bool = True
    metadata: Dict[str, Any] = field(default_factory=dict)
