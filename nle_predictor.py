"""
NLE-Predictor: Automated Non-Linear Effect Predictors and Thermodynamic Intelligence
for Asymmetric Catalysis.

Root runner and backwards-compatible facade delegating to `src/nle_predictor`.
"""

import os
import sys

# Ensure src package directory is in sys.path
_src_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "src"))
if _src_dir not in sys.path:
    sys.path.insert(0, _src_dir)

# Enable package-like submodule imports from src/nle_predictor
__path__ = [os.path.join(_src_dir, "nle_predictor")]
if "nle_predictor" in sys.modules:
    sys.modules["nle_predictor"].__path__ = __path__

# Re-export core constants and functions for 100% backward compatibility
from nle_predictor.cli import main
from nle_predictor.core.constants import (
    DEFAULT_TEMP,
    GAS_CONSTANT_R,
    HARTREE_TO_KCAL,
    HARTREE_TO_KJ,
    KCAL_TO_KJ,
    standard_state_correction_gas_to_solution,
)
from nle_predictor.core.types import (
    Conformer,
    EnsembleResult,
    FittingResult,
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
from nle_predictor.fitting.bootstrap import (
    BootstrapResult,
    run_bootstrap_monte_carlo,
)
from nle_predictor.fitting.experimental_fit import (
    fit_kagan_reservoir,
    load_experimental_csv,
)
from nle_predictor.fitting.temperature_fit import (
    VanTHoffFittingResult,
    fit_van_t_hoff_global,
    load_multi_temperature_csv,
)
from nle_predictor.models.base import NLEModel
from nle_predictor.models.diagnostics import (
    NLEDiagnostics,
    analyze_nle_diagnostics,
    calculate_area_under_asymmetry,
    calculate_curvature,
    calculate_initial_slope,
    calculate_max_amplification,
    find_inflection_points,
)
from nle_predictor.models.kagan_active import KaganActiveDimerModel
from nle_predictor.models.kagan_ml2 import (
    NLECalculator,
    write_sensitivity_report,
)
from nle_predictor.models.kagan_oligomers import KaganOligomerModel
from nle_predictor.models.noyori_daib import NoyoriDAIBModel
from nle_predictor.models.soai import SoaiAutocatalysisModel, SoaiTrajectory
from nle_predictor.qc.crest import ConformerSampler
from nle_predictor.qc.gallery import ChiralLigandGallery, LigandMetadata
from nle_predictor.qc.gaussian import (
    extract_gaussian_free_energy,
    validate_gaussian_output,
)
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
from nle_predictor.thermo.solvation import (
    SOLVENT_DATABASE,
    SolventProperties,
    compute_standard_state_dimerization_constant,
    concentration_correction_1atm_to_1m,
    get_solvent_properties,
)
from nle_predictor.viz.dashboard import generate_interactive_html
from nle_predictor.viz.plotting import (
    plot_amplification_contour,
    plot_monomer_speciation,
    plot_nle_with_experimental,
)

__all__ = [
    # Top-level CLI
    "main",
    # Constants
    "DEFAULT_TEMP",
    "GAS_CONSTANT_R",
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
    "FittingResult",
    "ModelType",
    # Geometry
    "read_xyz",
    "write_xyz",
    "validate_xyz_file",
    "calculate_rmsd",
    # Thermo & Solvation
    "boltzmann_weighted_free_energy",
    "boltzmann_weights",
    "conformational_entropy",
    "calculate_ensemble",
    "SOLVENT_DATABASE",
    "SolventProperties",
    "get_solvent_properties",
    "concentration_correction_1atm_to_1m",
    "compute_standard_state_dimerization_constant",
    # Models & Diagnostics
    "NLEModel",
    "NLECalculator",
    "KaganActiveDimerModel",
    "KaganOligomerModel",
    "NoyoriDAIBModel",
    "SoaiAutocatalysisModel",
    "SoaiTrajectory",
    "write_sensitivity_report",
    "NLEDiagnostics",
    "analyze_nle_diagnostics",
    "calculate_initial_slope",
    "calculate_max_amplification",
    "calculate_curvature",
    "find_inflection_points",
    "calculate_area_under_asymmetry",
    # Fitting
    "fit_kagan_reservoir",
    "load_experimental_csv",
    "VanTHoffFittingResult",
    "fit_van_t_hoff_global",
    "load_multi_temperature_csv",
    "BootstrapResult",
    "run_bootstrap_monte_carlo",
    # QC
    "Config",
    "DFTEvaluator",
    "ConformerSampler",
    "ChiralLigandGallery",
    "LigandMetadata",
    "validate_orca_output",
    "collect_software_versions",
    "validate_gaussian_output",
    "extract_gaussian_free_energy",
    # Viz
    "plot_nle_with_experimental",
    "plot_monomer_speciation",
    "plot_amplification_contour",
]


if __name__ == "__main__":
    main(sys.argv[1:])
