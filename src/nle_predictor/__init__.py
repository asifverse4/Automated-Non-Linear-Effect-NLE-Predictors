"""
NLE-Predictor: Automated Non-Linear Effect (NLE) Predictors
and Thermodynamic Intelligence for Asymmetric Catalysis.
"""

__version__ = "1.0.0"
__author__ = "Asif Raza"

from nle_predictor.cli import main as cli_main
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
    "cli_main",
    "generate_interactive_html",
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
    # Fitting & Regression
    "fit_kagan_reservoir",
    "load_experimental_csv",
    "FittingResult",
    "VanTHoffFittingResult",
    "fit_van_t_hoff_global",
    "load_multi_temperature_csv",
    "BootstrapResult",
    "run_bootstrap_monte_carlo",
    # Visualization
    "plot_nle_with_experimental",
    "plot_monomer_speciation",
    "plot_amplification_contour",
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
]
