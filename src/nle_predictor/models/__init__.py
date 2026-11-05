"""Mathematical models and diagnostic tools for non-linear effects in asymmetric catalysis."""

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
from nle_predictor.models.kagan_ml2 import NLECalculator, write_sensitivity_report
from nle_predictor.models.kagan_oligomers import KaganOligomerModel
from nle_predictor.models.noyori_daib import NoyoriDAIBModel

__all__ = [
    "NLEModel",
    "NLECalculator",
    "KaganActiveDimerModel",
    "KaganOligomerModel",
    "NoyoriDAIBModel",
    "write_sensitivity_report",
    "NLEDiagnostics",
    "analyze_nle_diagnostics",
    "calculate_initial_slope",
    "calculate_max_amplification",
    "calculate_curvature",
    "find_inflection_points",
    "calculate_area_under_asymmetry",
]
