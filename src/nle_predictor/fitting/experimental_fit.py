"""
Experimental NLE data regression and thermodynamic parameter estimation.
Extracts Delta G, K_hetero/K_homo, and intrinsic ee_0 with statistical confidence intervals.
"""

import csv
import logging
import os
from typing import Optional, Tuple
import numpy as np
from scipy.optimize import least_squares

from nle_predictor.core.constants import DEFAULT_TEMP, GAS_CONSTANT_R
from nle_predictor.core.types import FittingResult
from nle_predictor.models.kagan_ml2 import NLECalculator

logger = logging.getLogger("NLE_Predictor")


def load_experimental_csv(csv_path: str) -> Tuple[np.ndarray, np.ndarray]:
    """
    Load experimental catalyst ee vs product ee points from a CSV file.
    Accepts both percentage [0, 100] and fraction [0, 1.0] scales.
    Returns normalized arrays in range [0, 1.0].
    """
    if not os.path.isfile(csv_path):
        raise FileNotFoundError(f"Experimental CSV file not found: {csv_path}")

    cat_list = []
    prod_list = []
    with open(csv_path, "r", encoding="utf-8", errors="replace") as f:
        reader = csv.reader(f)
        header = next(reader, None)
        for row in reader:
            if not row or len(row) < 2:
                continue
            try:
                c = float(row[0].strip())
                p = float(row[1].strip())
                cat_list.append(c)
                prod_list.append(p)
            except ValueError:
                continue

    if len(cat_list) < 3:
        raise ValueError(f"Experimental CSV must contain at least 3 valid data pairs: {csv_path}")

    cat_arr = np.array(cat_list, dtype=np.float64)
    prod_arr = np.array(prod_list, dtype=np.float64)

    # Detect percentage scale
    if np.max(np.abs(cat_arr)) > 1.05 or np.max(np.abs(prod_arr)) > 1.05:
        cat_arr = cat_arr / 100.0
        prod_arr = prod_arr / 100.0

    return cat_arr, prod_arr


def fit_kagan_reservoir(
    ee_cat_exp: np.ndarray,
    ee_prod_exp: np.ndarray,
    temp: float = DEFAULT_TEMP,
    k_homo: float = 1e5,
    fit_ee_0: bool = False,
    initial_delta_g: float = 0.0,
) -> FittingResult:
    """
    Perform non-linear least squares regression to fit experimental data
    to the Kagan ML2 Reservoir model.

    Parameters
    ----------
    ee_cat_exp : np.ndarray
        Experimental catalyst ee (in range [0, 1.0]).
    ee_prod_exp : np.ndarray
        Experimental product ee (in range [0, 1.0]).
    temp : float
        Reaction temperature in Kelvin.
    k_homo : float
        Assumed absolute dimerization constant for homochiral species.
    fit_ee_0 : bool
        If True, simultaneously fit the intrinsic monomer enantioselectivity ee_0.
    initial_delta_g : float
        Initial guess for Delta G (RS - RR) in kcal/mol.

    Returns
    -------
    FittingResult
        Fitted parameters, standard errors, R^2, and residuals.
    """
    calculator = NLECalculator(temp=temp)

    def residuals(params):
        if fit_ee_0:
            dg, ee0 = params
        else:
            dg = params[0]
            ee0 = 1.0
        sim_prod = calculator.simulate_kagan_model(ee_cat_exp, delta_g_kcal=dg, k_homo_abs=k_homo, ee_0=ee0)
        return sim_prod - ee_prod_exp

    p0 = [initial_delta_g, 1.0] if fit_ee_0 else [initial_delta_g]
    bounds = ([-20.0, 0.0], [20.0, 1.0]) if fit_ee_0 else ([-20.0], [20.0])

    res = least_squares(residuals, p0, bounds=bounds, ftol=1e-8, xtol=1e-8)

    fitted_dg = float(res.x[0])
    fitted_ee0 = float(res.x[1]) if fit_ee_0 else 1.0

    # Calculate statistics
    res_vec = res.fun
    n_points = len(ee_cat_exp)
    n_params = len(res.x)
    ss_res = np.sum(res_vec**2)
    ss_tot = np.sum((ee_prod_exp - np.mean(ee_prod_exp))**2)
    r_squared = float(1.0 - (ss_res / ss_tot)) if ss_tot > 1e-12 else 1.0
    rmse = float(np.sqrt(ss_res / n_points))

    # Variance-covariance matrix
    # J^T J approximation
    jac = res.jac
    try:
        inv_jtj = np.linalg.pinv(np.dot(jac.T, jac))
        s_sq = ss_res / max(1, (n_points - n_params))
        cov = inv_jtj * s_sq
        dg_stderr = float(np.sqrt(cov[0, 0]))
    except Exception:
        dg_stderr = float("nan")

    # Ratio and its propagated error
    ratio = 2.0 * np.exp(-fitted_dg / (GAS_CONSTANT_R * temp))
    # d(ratio)/d(dg) = -ratio / (R*T)
    ratio_stderr = abs(ratio / (GAS_CONSTANT_R * temp)) * dg_stderr if np.isfinite(dg_stderr) else float("nan")

    return FittingResult(
        delta_g_kcal=fitted_dg,
        delta_g_stderr=dg_stderr,
        k_ratio=ratio,
        k_ratio_stderr=ratio_stderr,
        r_squared=r_squared,
        rmse=rmse,
        residuals=res_vec,
        converged=bool(res.success),
        metadata={
            "fitted_ee_0": fitted_ee0,
            "fit_ee_0_enabled": fit_ee_0,
            "temperature_k": temp,
            "k_homo": k_homo,
            "iterations": int(res.nfev),
        },
    )

