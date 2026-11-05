"""
Multi-Temperature van 't Hoff Deconvolution and Thermodynamic Parameter Estimation.

Physical Organic Chemistry Foundations:
Extracts Delta Delta H° (enthalpy difference) and Delta Delta S° (entropy difference)
governing homochiral vs heterochiral aggregation equilibria across multiple temperatures:
Delta Delta G°(T) = Delta Delta H° - T * Delta Delta S°
g(T) = exp(-Delta Delta G°(T) / (R * T))

Also evaluates:
1. Krug-Hunter-Grieger enthalpy-entropy compensation and correlation.
2. The isoinversion temperature T_inv = Delta Delta H° / Delta Delta S°, at which
   chiral discrimination vanishes (g = 1, perfectly linear NLE) or inverts.
3. Full asymptotic parameter covariance matrices and joint 95% confidence ellipses.
"""

import csv
import logging
import os
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import least_squares

from nle_predictor.core.constants import GAS_CONSTANT_R
from nle_predictor.models.kagan_ml2 import NLECalculator

logger = logging.getLogger("NLE_Predictor")


@dataclass
class VanTHoffFittingResult:
    """Thermodynamic deconvolution results across multi-temperature NLE curves."""
    delta_delta_h_kcal: float
    delta_delta_h_stderr: float
    delta_delta_s_cal_mol_k: float
    delta_delta_s_stderr: float
    covariance_matrix: np.ndarray
    correlation_coeff: float
    isoinversion_temp_k: Optional[float]
    isoinversion_stderr: Optional[float]
    fitted_temperatures: List[float]
    r_squared: float
    rmse: float
    residuals: np.ndarray
    g_at_temperatures: Dict[float, float]
    delta_g_at_temperatures: Dict[float, float]
    fit_ee_0: bool = False
    fitted_ee_0: float = 1.0
    ee_0_stderr: float = 0.0
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "delta_delta_h_kcal_mol": self.delta_delta_h_kcal,
            "delta_delta_h_stderr": self.delta_delta_h_stderr,
            "delta_delta_s_cal_mol_k": self.delta_delta_s_cal_mol_k,
            "delta_delta_s_stderr": self.delta_delta_s_stderr,
            "correlation_coefficient": self.correlation_coeff,
            "isoinversion_temperature_k": self.isoinversion_temp_k,
            "isoinversion_stderr_k": self.isoinversion_stderr,
            "r_squared": self.r_squared,
            "rmse": self.rmse,
            "fitted_temperatures": self.fitted_temperatures,
            "g_at_temperatures": self.g_at_temperatures,
            "delta_g_at_temperatures": self.delta_g_at_temperatures,
            "fitted_ee_0": self.fitted_ee_0 if self.fit_ee_0 else 1.0,
            "metadata": self.metadata,
        }

    def confidence_ellipse_2d(
        self,
        confidence: float = 0.95,
        n_points: int = 100,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute coordinates (H_coords, S_coords) of the joint confidence ellipse.

        Parameters
        ----------
        confidence : float
            Confidence level (default 0.95 for 95% ellipse).
        n_points : int
            Number of perimeter points.

        Returns
        -------
        Tuple[np.ndarray, np.ndarray]
            (H_points_in_kcal, S_points_in_cal_mol_k)
        """
        # Critical chi-squared value for 2 degrees of freedom
        # chi2(2, 0.95) = 5.991
        from scipy.stats import chi2
        critical_val = chi2.ppf(confidence, df=2)

        cov_sub = self.covariance_matrix[:2, :2]
        eigvals, eigvecs = np.linalg.eigh(cov_sub)
        eigvals = np.maximum(eigvals, 1e-12)

        theta = np.linspace(0, 2 * np.pi, n_points)
        circle = np.stack([np.cos(theta), np.sin(theta)])  # (2, N)

        # Scale by eigenvalues and critical value
        radii = np.sqrt(eigvals * critical_val)  # (2,)
        scaled = eigvecs @ (radii[:, None] * circle)

        h_center = self.delta_delta_h_kcal
        s_center = self.delta_delta_s_cal_mol_k

        return scaled[0] + h_center, scaled[1] + s_center


def load_multi_temperature_csv(csv_path: str) -> Dict[float, Tuple[np.ndarray, np.ndarray]]:
    """
    Load multi-temperature NLE datasets from CSV.
    Expected header: temperature,ee_cat,ee_prod (or temp,ee_catalyst,ee_product).
    Values can be in Kelvin or Celsius; numbers < 100 are treated as Celsius and converted to K.
    """
    if not os.path.isfile(csv_path):
        raise FileNotFoundError(f"Multi-temperature CSV not found: {csv_path}")

    temp_data: Dict[float, List[Tuple[float, float]]] = {}
    with open(csv_path, "r", encoding="utf-8", errors="replace") as f:
        reader = csv.reader(f)
        header = next(reader, None)
        if not header or len(header) < 3:
            raise ValueError("CSV header must have at least 3 columns: temperature, ee_cat, ee_prod.")

        for row in reader:
            if not row or len(row) < 3:
                continue
            try:
                t_val = float(row[0].strip())
                c_val = float(row[1].strip())
                p_val = float(row[2].strip())

                # Convert Celsius to Kelvin if needed
                if t_val < 150.0:
                    t_val = t_val + 273.15

                # Detect percentage scale
                if abs(c_val) > 1.05 or abs(p_val) > 1.05:
                    c_val = c_val / 100.0
                    p_val = p_val / 100.0

                if t_val not in temp_data:
                    temp_data[t_val] = []
                temp_data[t_val].append((c_val, p_val))
            except ValueError:
                continue

    result: Dict[float, Tuple[np.ndarray, np.ndarray]] = {}
    for t_val, pairs in temp_data.items():
        if len(pairs) < 2:
            continue
        arr = np.array(pairs, dtype=float)
        result[t_val] = (arr[:, 0], arr[:, 1])

    if len(result) < 2:
        raise ValueError(f"Need at least 2 distinct temperatures with data in {csv_path}")

    return result


def fit_van_t_hoff_global(
    temperature_datasets: Dict[float, Tuple[np.ndarray, np.ndarray]],
    k_homo: float = 1e5,
    fit_ee_0: bool = False,
    initial_dh: float = 0.0,
    initial_ds: float = 0.0,
) -> VanTHoffFittingResult:
    """
    Perform simultaneous global non-linear least squares fit across multi-temperature NLE curves.

    Parameters
    ----------
    temperature_datasets : Dict[float, Tuple[np.ndarray, np.ndarray]]
        Mapping of temperature in Kelvin -> (ee_cat_array, ee_prod_array).
    k_homo : float
        Homochiral association constant (default 1e5).
    fit_ee_0 : bool
        Whether to simultaneously fit intrinsic enantioselectivity ee_0.
    initial_dh : float
        Initial guess for Delta Delta H° in kcal/mol.
    initial_ds : float
        Initial guess for Delta Delta S° in cal/(mol·K).

    Returns
    -------
    VanTHoffFittingResult
        Complete thermodynamic deconvolution with confidence intervals and ellipse.
    """
    if len(temperature_datasets) < 2:
        raise ValueError("Must provide at least 2 temperature datasets for van 't Hoff deconvolution.")

    sorted_temps = sorted(temperature_datasets.keys())
    calculators = {t: NLECalculator(temp=t) for t in sorted_temps}

    # Flatten experimental arrays
    all_exp_prod: List[float] = []
    for t in sorted_temps:
        _, p_exp = temperature_datasets[t]
        all_exp_prod.extend(p_exp.tolist())
    y_exp = np.array(all_exp_prod, dtype=float)

    # Initial parameter vector
    # params = [Delta Delta H (kcal/mol), Delta Delta S (cal/mol·K), (optional ee_0)]
    if fit_ee_0:
        p0 = [initial_dh, initial_ds, 0.99]
        bounds = ([-50.0, -100.0, 0.0], [50.0, 100.0, 1.0])
    else:
        p0 = [initial_dh, initial_ds]
        bounds = ([-50.0, -100.0], [50.0, 100.0])

    def residual_fun(params: Sequence[float]) -> np.ndarray:
        dh = params[0]
        ds = params[1]
        ee0 = params[2] if fit_ee_0 else 1.0

        res_list: List[float] = []
        for t in sorted_temps:
            c_exp, p_exp = temperature_datasets[t]
            # Delta Delta G(T) = Delta Delta H - T * (Delta Delta S / 1000) in kcal/mol
            dg_t = dh - t * (ds / 1000.0)
            calc = calculators[t]
            p_sim = calc.simulate_kagan_model(c_exp, delta_g_kcal=dg_t, k_homo_abs=k_homo, ee_0=ee0)
            res_list.extend((p_sim - p_exp).tolist())

        return np.array(res_list, dtype=float)

    opt_res = least_squares(
        residual_fun,
        x0=p0,
        bounds=bounds,
        method="trf",
        ftol=1e-9,
        xtol=1e-9,
        gtol=1e-9,
    )

    fitted_dh = float(opt_res.x[0])
    fitted_ds = float(opt_res.x[1])
    fitted_ee0 = float(opt_res.x[2]) if fit_ee_0 else 1.0

    residuals = opt_res.fun
    n_obs = len(residuals)
    n_params = len(p0)
    dof = max(1, n_obs - n_params)
    mse = float(np.sum(residuals ** 2) / dof)
    rmse = float(np.sqrt(np.mean(residuals ** 2)))

    # Compute covariance matrix from Jacobian J: Cov = s^2 * (J^T J)^-1
    jac = opt_res.jac
    try:
        jt_j = jac.T @ jac
        inv_jt_j = np.linalg.pinv(jt_j)
        cov_matrix = mse * inv_jt_j
        stderrs = np.sqrt(np.diag(cov_matrix))
    except Exception:
        cov_matrix = np.zeros((n_params, n_params))
        stderrs = np.zeros(n_params)

    dh_err = float(stderrs[0])
    ds_err = float(stderrs[1])
    ee0_err = float(stderrs[2]) if fit_ee_0 else 0.0

    # Correlation coefficient r = Cov(H, S) / (sigma_H * sigma_S)
    denom = dh_err * ds_err
    corr = float(cov_matrix[0, 1] / denom) if denom > 1e-15 else 0.0

    # Isoinversion temperature: T_inv = (Delta Delta H * 1000) / Delta Delta S
    if abs(fitted_ds) > 1e-4:
        t_inv = (fitted_dh * 1000.0) / fitted_ds
        # Error propagation for ratio A/B:
        # (sigma_T / T)^2 = (sigma_H / H)^2 + (sigma_S / S)^2 - 2 Cov(H, S)/(H S)
        var_rel = (
            (dh_err / fitted_dh) ** 2
            + (ds_err / fitted_ds) ** 2
            - 2.0 * (cov_matrix[0, 1] / (fitted_dh * fitted_ds))
        )
        t_inv_err = float(abs(t_inv) * np.sqrt(max(0.0, var_rel)))
    else:
        t_inv = None
        t_inv_err = None

    # Calculate R^2
    ss_tot = float(np.sum((y_exp - np.mean(y_exp)) ** 2))
    ss_res = float(np.sum(residuals ** 2))
    r_squared = 1.0 - (ss_res / ss_tot) if ss_tot > 1e-15 else 1.0

    # Speciation at each temperature
    g_at_t: Dict[float, float] = {}
    dg_at_t: Dict[float, float] = {}
    for t in sorted_temps:
        dg_val = fitted_dh - t * (fitted_ds / 1000.0)
        dg_at_t[t] = float(dg_val)
        g_at_t[t] = float(np.exp(-dg_val / (GAS_CONSTANT_R * t)))

    return VanTHoffFittingResult(
        delta_delta_h_kcal=fitted_dh,
        delta_delta_h_stderr=dh_err,
        delta_delta_s_cal_mol_k=fitted_ds,
        delta_delta_s_stderr=ds_err,
        covariance_matrix=cov_matrix,
        correlation_coeff=corr,
        isoinversion_temp_k=t_inv,
        isoinversion_stderr=t_inv_err,
        fitted_temperatures=sorted_temps,
        r_squared=float(r_squared),
        rmse=rmse,
        residuals=residuals,
        g_at_temperatures=g_at_t,
        delta_g_at_temperatures=dg_at_t,
        fit_ee_0=fit_ee_0,
        fitted_ee_0=fitted_ee0,
        ee_0_stderr=ee0_err,
        metadata={
            "k_homo": k_homo,
            "n_temperatures": len(sorted_temps),
            "n_observations": n_obs,
            "degrees_of_freedom": dof,
        },
    )
