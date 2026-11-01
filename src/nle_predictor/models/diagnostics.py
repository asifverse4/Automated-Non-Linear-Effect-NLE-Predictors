"""
Quantitative mathematical diagnostics for Non-Linear Effects in Asymmetric Catalysis.
Calculates initial slope (amplification at racemic state), curvature, inflection points,
and integral asymmetry metrics.
"""

from dataclasses import dataclass
from typing import List, Tuple
import numpy as np


@dataclass
class NLEDiagnostics:
    """Summary of quantitative curve diagnostics for an NLE trajectory."""
    initial_slope: float
    max_amplification: float
    catalyst_ee_at_max_amplification: float
    area_under_asymmetry: float
    inflection_points: List[float]
    classification: str
    is_hyperbolic: bool
    is_sigmoidal: bool

    def to_dict(self):
        return {
            "initial_slope_at_racemic": self.initial_slope,
            "max_amplification_factor": self.max_amplification,
            "catalyst_ee_at_max_amplification": self.catalyst_ee_at_max_amplification,
            "area_under_asymmetry": self.area_under_asymmetry,
            "inflection_points": self.inflection_points,
            "classification": self.classification,
            "is_hyperbolic": self.is_hyperbolic,
            "is_sigmoidal": self.is_sigmoidal,
        }


def calculate_initial_slope(ee_cat: np.ndarray, ee_prod: np.ndarray) -> float:
    """
    Calculate the initial slope d(ee_prod)/d(ee_cat) as ee_cat -> 0.
    Values > 1 indicate asymmetric amplification at low enantiopurity.
    """
    if len(ee_cat) < 3:
        raise ValueError("At least 3 points required for initial slope estimation.")
    # Fit cubic or quadratic polynomial near origin (first 5 points or 20% of data)
    n_pts = min(max(5, len(ee_cat) // 5), len(ee_cat))
    x_sub = ee_cat[:n_pts]
    y_sub = ee_prod[:n_pts]

    # Polynomial fit y = c1 * x + c2 * x^2 ... passing near origin
    poly_coeffs = np.polyfit(x_sub, y_sub, deg=min(2, n_pts - 1))
    # First derivative at x=0
    # For degree 2: p(x) = c0*x^2 + c1*x + c2 -> p'(0) = c1
    return float(poly_coeffs[-2]) if len(poly_coeffs) >= 2 else 1.0


def calculate_max_amplification(ee_cat: np.ndarray, ee_prod: np.ndarray) -> Tuple[float, float]:
    """
    Compute maximum amplification alpha = max(ee_prod / ee_cat) for ee_cat > 0
    and the corresponding catalyst ee.
    """
    mask = ee_cat > 1e-4
    if not np.any(mask):
        return 1.0, 0.0
    ratios = ee_prod[mask] / ee_cat[mask]
    max_idx = np.argmax(ratios)
    return float(ratios[max_idx]), float(ee_cat[mask][max_idx])


def calculate_curvature(ee_cat: np.ndarray, ee_prod: np.ndarray) -> np.ndarray:
    """
    Calculate geometric curvature kappa(x) = y'' / (1 + (y')^2)^(3/2).
    """
    dy_dx = np.gradient(ee_prod, ee_cat)
    d2y_dx2 = np.gradient(dy_dx, ee_cat)
    curvature = d2y_dx2 / np.power(1.0 + dy_dx**2, 1.5)
    return curvature


def find_inflection_points(
    ee_cat: np.ndarray,
    ee_prod: np.ndarray,
    threshold: float = 1e-3,
) -> List[float]:
    """
    Identify catalyst ee values where second derivative changes sign (d^2(ee_prod)/d(ee_cat)^2 = 0).
    Ignores numerical floating-point noise below threshold.
    """
    if len(ee_cat) < 5:
        return []
    dy_dx = np.gradient(ee_prod, ee_cat)
    d2y_dx2 = np.gradient(dy_dx, ee_cat)

    # If the entire second derivative is essentially zero, it's linear (no inflections)
    if np.max(np.abs(d2y_dx2)) < threshold:
        return []

    inflection_pts: List[float] = []
    for i in range(1, len(d2y_dx2) - 1):
        if (
            d2y_dx2[i - 1] * d2y_dx2[i + 1] < 0
            and max(abs(d2y_dx2[i - 1]), abs(d2y_dx2[i + 1])) > threshold
        ):
            # Linear zero-crossing interpolation
            x0, x1 = ee_cat[i - 1], ee_cat[i + 1]
            y0, y1 = d2y_dx2[i - 1], d2y_dx2[i + 1]
            if (y1 - y0) != 0:
                zero_x = x0 - y0 * (x1 - x0) / (y1 - y0)
                if 0.0 < zero_x < 1.0:
                    inflection_pts.append(round(float(zero_x), 4))
    return sorted(list(set(inflection_pts)))


def calculate_area_under_asymmetry(ee_cat: np.ndarray, ee_prod: np.ndarray) -> float:
    """
    Compute Area Under Asymmetry: integral from 0 to 1 of (ee_prod - ee_cat) dx.
    Positive indicates overall (+)NLE amplification; negative indicates (-)NLE depletion.
    """
    deviation = ee_prod - ee_cat
    return float(np.trapezoid(deviation, ee_cat))


def analyze_nle_diagnostics(ee_cat: np.ndarray, ee_prod: np.ndarray) -> NLEDiagnostics:
    """
    Run full diagnostic suite on an NLE curve.
    """
    s0 = calculate_initial_slope(ee_cat, ee_prod)
    max_amp, cat_ee_amp = calculate_max_amplification(ee_cat, ee_prod)
    aua = calculate_area_under_asymmetry(ee_cat, ee_prod)
    inflections = find_inflection_points(ee_cat, ee_prod)

    tol = 0.005
    if aua > tol:
        classification = "positive NLE (amplification)"
    elif aua < -tol:
        classification = "negative NLE (depletion)"
    else:
        classification = "linear (no effect)"

    is_sigmoidal = len(inflections) > 0
    is_hyperbolic = (not is_sigmoidal) and (max_amp > 1.05 or aua > tol)

    return NLEDiagnostics(
        initial_slope=s0,
        max_amplification=max_amp,
        catalyst_ee_at_max_amplification=cat_ee_amp,
        area_under_asymmetry=aua,
        inflection_points=inflections,
        classification=classification,
        is_hyperbolic=is_hyperbolic,
        is_sigmoidal=is_sigmoidal,
    )
