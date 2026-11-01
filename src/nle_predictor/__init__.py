"""
NLE-Predictor: Automated Non-Linear Effect (NLE) Predictors
and Thermodynamic Intelligence for Asymmetric Catalysis.
"""

__version__ = "1.0.0"
__author__ = "Asif Raza"

from nle_predictor.core.constants import (
    DEFAULT_TEMP,
    GAS_CONSTANT_R,
    GAS_CONSTANT_R_SI,
    HARTREE_TO_KCAL,
    HARTREE_TO_KJ,
    KCAL_TO_KJ,
    standard_state_correction_gas_to_solution,
)
from nle_predictor.core.types import (
    Conformer,
    EnsembleResult,
    ModelType,
    NLESimulationResult,
    ThermodynamicParameters,
    XYZStructure,
)
from nle_predictor.core.xyz import (
    calculate_rmsd,
    read_xyz,
    validate_xyz_file,
    write_xyz,
)
from nle_predictor.models.base import NLEModel
from nle_predictor.models.kagan_ml2 import (
    NLECalculator,
    write_sensitivity_report,
)
from nle_predictor.qc.crest import ConformerSampler
from nle_predictor.qc.orca import (
    Config,
    DFTEvaluator,
    collect_software_versions,
    validate_orca_output,
)
from nle_predictor.thermo.boltzmann import (
    boltzmann_weighted_free_energy,
    boltzmann_weights,
    calculate_ensemble,
    conformational_entropy,
)

__all__ = [
    # Constants
    "DEFAULT_TEMP",
    "GAS_CONSTANT_R",
    "GAS_CONSTANT_R_SI",
    "HARTREE_TO_KCAL",
    "HARTREE_TO_KJ",
    "KCAL_TO_KJ",
    "standard_state_correction_gas_to_solution",
    # Types
    "XYZStructure",
    "Conformer",
    "EnsembleResult",
    "ThermodynamicParameters",
    "NLESimulationResult",
    "ModelType",
    # XYZ
    "read_xyz",
    "write_xyz",
    "validate_xyz_file",
    "calculate_rmsd",
    # Thermo
    "boltzmann_weighted_free_energy",
    "boltzmann_weights",
    "conformational_entropy",
    "calculate_ensemble",
    # Models
    "NLEModel",
    "NLECalculator",
    "write_sensitivity_report",
    # QC
    "Config",
    "DFTEvaluator",
    "ConformerSampler",
    "validate_orca_output",
    "collect_software_versions",
]
