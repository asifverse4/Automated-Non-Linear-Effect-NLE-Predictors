"""
Bootstrap Monte Carlo Uncertainty Quantification for Non-Linear Effect Parameters.

Theoretical foundations:
B. Efron, R. J. Tibshirani, "An Introduction to the Bootstrap", Chapman & Hall/CRC, 1993.
J. R. Alston et al., Physical Organic Chemistry / Chemometrics parameter uncertainty propagation.

Replaces restrictive asymptotic Fisher information linear approximations with empirical
non-parametric resampling:
1. Residual resampling bootstrap (preserving catalyst ee grid design).
2. Case resampling bootstrap (non-parametric observational pair sampling).
3. Exact empirical percentile 95% confidence intervals, parameter bias, skewness,
   and 2D joint parameter correlation distributions.
"""

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, Optional, Tuple

import numpy as np
from scipy.optimize import least_squares

from nle_predictor.core.constants import DEFAULT_TEMP, GAS_CONSTANT_R
from nle_predictor.models.kagan_ml2 import NLECalculator

logger = logging.getLogger("NLE_Predictor")


@dataclass
class BootstrapResult:
    """Non-parametric bootstrap uncertainty quantification metrics."""
    n_bootstraps: int
    delta_g_samples: np.ndarray
    delta_g_mean: float
    delta_g_median: float
    delta_g_std: float
    delta_g_ci_95: Tuple[float, float]
    k_ratio_samples: np.ndarray
    k_ratio_median: float
    k_ratio_ci_95: Tuple[float, float]
    bias_delta_g: float
    fit_ee_0: bool = False
    ee_0_samples: Optional[np.ndarray] = None
    ee_0_ci_95: Optional[Tuple[float, float]] = None
    ee_0_median: Optional[float] = None
    parameter_correlation: float = 0.0
    converged_fraction: float = 1.0
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "n_bootstraps": self.n_bootstraps,
            "converged_fraction": self.converged_fraction,
            "delta_g_kcal_mol": {
                "mean": float(self.delta_g_mean),
                "median": float(self.delta_g_median),
                "std": float(self.delta_g_std),
                "ci_95_lower": float(self.delta_g_ci_95[0]),
                "ci_95_upper": float(self.delta_g_ci_95[1]),
                "bias": float(self.bias_delta_g),
            },
            "k_ratio": {
                "median": float(self.k_ratio_median),
                "ci_95_lower": float(self.k_ratio_ci_95[0]),
                "ci_95_upper": float(self.k_ratio_ci_95[1]),
            },
            "fit_ee_0": self.fit_ee_0,
            "metadata": self.metadata,
        }

    def percentile_interval(self, alpha: float = 0.05) -> Tuple[float, float]:
        """Compute arbitrary (1 - alpha) percentile confidence interval for Delta G."""
        lower_pct = 100.0 * (alpha / 2.0)
        upper_pct = 100.0 * (1.0 - alpha / 2.0)
        return (
            float(np.percentile(self.delta_g_samples, lower_pct)),
            float(np.percentile(self.delta_g_samples, upper_pct)),
        )


def run_bootstrap_monte_carlo(
    ee_cat_exp: np.ndarray,
    ee_prod_exp: np.ndarray,
    n_bootstraps: int = 500,
    method: str = "residual",
    temp: float = DEFAULT_TEMP,
    k_homo: float = 1e5,
    fit_ee_0: bool = False,
    seed: Optional[int] = 42,
) -> BootstrapResult:
    """
    Perform non-parametric bootstrap Monte Carlo resampling on experimental NLE data.

    Parameters
    ----------
    ee_cat_exp : np.ndarray
        Experimental catalyst enantiomeric excess values in [0, 1].
    ee_prod_exp : np.ndarray
        Experimental product enantiomeric excess values in [0, 1].
    n_bootstraps : int
        Number of bootstrap resampling replications (default 500).
    method : str
        Resampling strategy: 'residual' (recommended) or 'case'.
    temp : float
        Reaction temperature in Kelvin.
    k_homo : float
        Homochiral association constant.
    fit_ee_0 : bool
        Whether to simultaneously fit intrinsic monomer enantioselectivity ee_0.
    seed : Optional[int]
        Random seed for reproducible Monte Carlo draws.

    Returns
    -------
    BootstrapResult
        Empirical bootstrap sampling distributions, percentiles, and bias metrics.
    """
    if len(ee_cat_exp) != len(ee_prod_exp):
        raise ValueError("ee_cat_exp and ee_prod_exp must have identical lengths.")
    if len(ee_cat_exp) < 3:
        raise ValueError("Need at least 3 data points for bootstrap resampling.")
    if n_bootstraps < 10:
        raise ValueError("n_bootstraps must be >= 10.")
    if method not in ("residual", "case"):
        raise ValueError(f"Unknown bootstrap method '{method}', must be 'residual' or 'case'.")

    rng = np.random.default_rng(seed)
    calculator = NLECalculator(temp=temp)

    # 1. Base fit on original experimental dataset
    def fit_single(x_data: np.ndarray, y_data: np.ndarray, p0: np.ndarray) -> Optional[np.ndarray]:
        def obj(p):
            dg = p[0]
            ee0 = p[1] if fit_ee_0 else 1.0
            y_calc = calculator.simulate_kagan_model(x_data, delta_g_kcal=dg, k_homo_abs=k_homo, ee_0=ee0)
            return y_calc - y_data

        bnds = ([-20.0, 0.0], [20.0, 1.0]) if fit_ee_0 else ([-20.0], [20.0])
        res = least_squares(obj, x0=p0, bounds=bnds, method="trf", ftol=1e-7, xtol=1e-7)
        return res.x if res.success else None

    init_p = np.array([0.0, 0.99]) if fit_ee_0 else np.array([0.0])
    base_opt = fit_single(ee_cat_exp, ee_prod_exp, init_p)
    if base_opt is None:
        raise RuntimeError("Initial model fit failed to converge on the base experimental dataset.")

    base_dg = float(base_opt[0])
    base_ee0 = float(base_opt[1]) if fit_ee_0 else 1.0

    # Calculate model predictions and centered residuals
    base_pred = calculator.simulate_kagan_model(ee_cat_exp, delta_g_kcal=base_dg, k_homo_abs=k_homo, ee_0=base_ee0)
    residuals = ee_prod_exp - base_pred
    centered_residuals = residuals - np.mean(residuals)

    n_points = len(ee_cat_exp)
    dg_samples = []
    ee0_samples = []

    for _ in range(n_bootstraps):
        if method == "residual":
            # Resample residuals with replacement
            sampled_res = rng.choice(centered_residuals, size=n_points, replace=True)
            y_boot = np.clip(base_pred + sampled_res, -1.0, 1.0)
            x_boot = ee_cat_exp
        else:
            # Case resampling: resample (x, y) pairs with replacement
            indices = rng.choice(n_points, size=n_points, replace=True)
            x_boot = ee_cat_exp[indices]
            y_boot = ee_prod_exp[indices]

        sol_boot = fit_single(x_boot, y_boot, p0=base_opt)
        if sol_boot is not None and np.all(np.isfinite(sol_boot)):
            dg_samples.append(float(sol_boot[0]))
            if fit_ee_0:
                ee0_samples.append(float(sol_boot[1]))

    if len(dg_samples) < 5:
        raise RuntimeError("Bootstrap failure: fewer than 5 resampled fits converged.")

    dg_arr = np.array(dg_samples, dtype=float)
    dg_mean = float(np.mean(dg_arr))
    dg_median = float(np.median(dg_arr))
    dg_std = float(np.std(dg_arr, ddof=1))
    dg_ci = (float(np.percentile(dg_arr, 2.5)), float(np.percentile(dg_arr, 97.5)))
    bias_dg = float(dg_mean - base_dg)

    # K_hetero / K_homo ratio samples
    k_ratios = 2.0 * np.exp(-dg_arr / (GAS_CONSTANT_R * temp))
    k_ratio_med = float(np.median(k_ratios))
    k_ratio_ci = (float(np.percentile(k_ratios, 2.5)), float(np.percentile(k_ratios, 97.5)))

    corr = 0.0
    ee0_arr = None
    ee0_ci = None
    ee0_med = None
    if fit_ee_0 and len(ee0_samples) > 0:
        ee0_arr = np.array(ee0_samples, dtype=float)
        ee0_med = float(np.median(ee0_arr))
        ee0_ci = (float(np.percentile(ee0_arr, 2.5)), float(np.percentile(ee0_arr, 97.5)))
        if np.std(ee0_arr) > 1e-12 and dg_std > 1e-12:
            corr = float(np.corrcoef(dg_arr, ee0_arr)[0, 1])

    return BootstrapResult(
        n_bootstraps=n_bootstraps,
        delta_g_samples=dg_arr,
        delta_g_mean=dg_mean,
        delta_g_median=dg_median,
        delta_g_std=dg_std,
        delta_g_ci_95=dg_ci,
        k_ratio_samples=k_ratios,
        k_ratio_median=k_ratio_med,
        k_ratio_ci_95=k_ratio_ci,
        bias_delta_g=bias_dg,
        fit_ee_0=fit_ee_0,
        ee_0_samples=ee0_arr,
        ee_0_ci_95=ee0_ci,
        ee_0_median=ee0_med,
        parameter_correlation=corr,
        converged_fraction=float(len(dg_samples) / n_bootstraps),
        metadata={
            "method": method,
            "base_delta_g_kcal": base_dg,
            "base_ee_0": base_ee0,
            "temperature_k": temp,
            "k_homo": k_homo,
        },
    )
