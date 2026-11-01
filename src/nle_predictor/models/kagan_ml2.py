"""
Kagan ML2 Reservoir Model for Non-Linear Effects in Asymmetric Catalysis.
Henri B. Kagan's foundational monomer-dimer reservoir equilibrium model.
"""

import csv
import logging
import os
from typing import Sequence, Tuple
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import root

from nle_predictor.core.constants import DEFAULT_TEMP, GAS_CONSTANT_R, HARTREE_TO_KCAL
from nle_predictor.core.types import NLESimulationResult
from nle_predictor.models.base import NLEModel

logger = logging.getLogger("NLE_Predictor")


class NLECalculator(NLEModel):
    """
    Calculates equilibrium constants and simulates NLE curves
    using a rigorous numerical solver for the ML2 (Reservoir) Kagan Model.
    """
    def __init__(self, temp: float = DEFAULT_TEMP):
        if temp <= 0:
            raise ValueError("Temperature must be greater than zero Kelvin.")
        self.temp = temp
        self.last_k_hetero = 0.0
        self.last_k_ratio = 1.0

    def calc_equilibrium_constant(self, g_rr: float, g_rs: float) -> float:
        """
        Calculates the thermodynamic parameters.
        Delta G = G_RS - G_RR in kcal/mol.
        Statistical symmetry factor of 2 is applied during association constant derivation.
        """
        delta_g_hartree = g_rs - g_rr
        delta_g_kcal = delta_g_hartree * HARTREE_TO_KCAL

        logger.info("--- Thermodynamic Results ---")
        logger.info(f"G(RR) = {g_rr:.6f} Eh")
        logger.info(f"G(RS) = {g_rs:.6f} Eh")
        logger.info(f"Delta G (RS - RR): {delta_g_kcal:.4f} kcal/mol")

        return delta_g_kcal

    def solve_mass_balance(
        self,
        ee_cat: float,
        k_homo: float,
        k_hetero: float,
    ) -> Tuple[float, float]:
        """
        Numerically solves the coupled non-linear mass balance equations for Monomers (M_R, M_S).
        Total Catalyst = C_R + C_S = 1.0
        Analytical Jacobian is used for quadratic convergence and numerical stability.
        """
        if not -1.0 <= ee_cat <= 1.0:
            raise ValueError("Catalyst ee must be between -1.0 and 1.0")
        if k_homo < 0 or k_hetero < 0:
            raise ValueError("Association constants must be non-negative")

        C_R = 0.5 * (1.0 + ee_cat)
        C_S = 0.5 * (1.0 - ee_cat)

        def equations(variables):
            M_R, M_S = variables
            eq1 = M_R + 2.0 * k_homo * (M_R**2) + k_hetero * M_R * M_S - C_R
            eq2 = M_S + 2.0 * k_homo * (M_S**2) + k_hetero * M_R * M_S - C_S
            return [eq1, eq2]

        def jacobian(variables):
            M_R, M_S = variables
            df1_dmr = 1.0 + 4.0 * k_homo * M_R + k_hetero * M_S
            df1_dms = k_hetero * M_R
            df2_dmr = k_hetero * M_S
            df2_dms = 1.0 + 4.0 * k_homo * M_S + k_hetero * M_R
            return [[df1_dmr, df1_dms], [df2_dmr, df2_dms]]

        # Initial guess based on pure homochiral limit
        mr_guess = (-1.0 + np.sqrt(1.0 + 8.0 * k_homo * C_R)) / (4.0 * k_homo) if (C_R > 0 and k_homo > 0) else C_R
        ms_guess = (-1.0 + np.sqrt(1.0 + 8.0 * k_homo * C_S)) / (4.0 * k_homo) if (C_S > 0 and k_homo > 0) else C_S

        sol = root(equations, [mr_guess, ms_guess], jac=jacobian, method="hybr")

        residual = np.max(np.abs(equations(sol.x))) if sol.success else np.inf
        if (
            not sol.success
            or not np.all(np.isfinite(sol.x))
            or np.any(sol.x < -1e-10)
            or residual > 1e-8
        ):
            raise RuntimeError(
                f"Mass-balance solve failed at ee_cat={ee_cat}: "
                f"{sol.message}; residual={residual:.3e}"
            )

        return max(0.0, float(sol.x[0])), max(0.0, float(sol.x[1]))

    def simulate_kagan_model(
        self,
        ee_cat_range: Sequence[float],
        delta_g_kcal: float,
        k_homo_abs: float = 1e5,
        ee_0: float = 1.0,
    ) -> np.ndarray:
        """
        Simulates the NLE curve.
        k_homo_abs: Absolute association constant for R + R -> RR (Default 1e5 for reservoir limit)
        ee_0: Intrinsic enantioselectivity of the active monomer (default 1.0 = 100%)
        """
        if self.temp <= 0:
            raise ValueError("Temperature must be greater than zero")
        if k_homo_abs < 0:
            raise ValueError("k_homo_abs must be non-negative")

        # Statistical factor of 2 accounts for hetero vs homo symmetry
        ratio = 2.0 * np.exp(-delta_g_kcal / (GAS_CONSTANT_R * self.temp))
        k_hetero_abs = k_homo_abs * ratio
        self.last_k_hetero = k_hetero_abs
        self.last_k_ratio = ratio

        logger.info(f"Assumed K_homo (RR binding): {k_homo_abs:.2e}")
        logger.info(f"Calculated K_hetero (RS binding): {k_hetero_abs:.2e}")
        logger.info(f"K_hetero / K_homo ratio: {ratio:.4f}")

        ee_prods = []
        for ee_cat in ee_cat_range:
            M_R, M_S = self.solve_mass_balance(ee_cat, k_homo_abs, k_hetero_abs)
            total_active = M_R + M_S
            if total_active <= 1e-12:
                ee_p = 0.0
            else:
                ee_p = ee_0 * ((M_R - M_S) / total_active)
            ee_prods.append(ee_p)

        return np.array(ee_prods)

    def simulate(
        self,
        ee_cat_range: Sequence[float],
        delta_g_kcal: float,
        k_homo: float = 1e5,
        ee_0: float = 1.0,
        **kwargs,
    ) -> NLESimulationResult:
        """Simulate and return a comprehensive typed NLESimulationResult."""
        ee_cat = np.asarray(ee_cat_range, dtype=float)
        ee_prod = self.simulate_kagan_model(ee_cat, delta_g_kcal, k_homo_abs=k_homo, ee_0=ee_0)

        monomers_r = []
        monomers_s = []
        for c in ee_cat:
            mr, ms = self.solve_mass_balance(c, k_homo, self.last_k_hetero)
            monomers_r.append(mr)
            monomers_s.append(ms)

        diff = ee_prod - ee_cat
        tol = 1e-6
        if np.max(diff) > tol:
            classification = "positive NLE"
        elif np.min(diff) < -tol:
            classification = "negative NLE"
        else:
            classification = "near-linear"

        product_at_50 = float(np.interp(0.5, ee_cat, ee_prod))
        amplification = np.divide(
            ee_prod[1:], ee_cat[1:], out=np.zeros_like(ee_prod[1:]), where=ee_cat[1:] > 0
        )
        max_amp = float(np.max(amplification)) if len(amplification) > 0 else 1.0

        return NLESimulationResult(
            model_name="Kagan ML2 Reservoir Model",
            ee_catalyst=ee_cat,
            ee_product=ee_prod,
            monomers_r=np.array(monomers_r),
            monomers_s=np.array(monomers_s),
            classification=classification,
            product_at_50=product_at_50,
            max_amplification=max_amp,
            max_deviation_from_linearity=float(np.max(np.abs(diff))),
            metadata={
                "k_homo": k_homo,
                "k_hetero": self.last_k_hetero,
                "k_ratio": self.last_k_ratio,
                "delta_g_kcal": delta_g_kcal,
                "temperature_k": self.temp,
            },
        )

    def export_data(self, ee_cat: np.ndarray, ee_prod: np.ndarray, filepath: str = "nle_results.csv") -> None:
        """Exports raw NLE data to a CSV for external plotting."""
        with open(filepath, mode="w", newline="", encoding="utf-8") as file:
            writer = csv.writer(file)
            writer.writerow(["ee_catalyst(%)", "ee_product(%)"])
            for cat, prod in zip(ee_cat * 100, ee_prod * 100):
                writer.writerow([f"{cat:.2f}", f"{prod:.2f}"])
        logger.info(f"Data successfully exported to {filepath}")

    def plot_nle_curve(self, ee_cat: np.ndarray, ee_prod: np.ndarray, save_path: str = "nle_curve.png") -> None:
        """Generates a publication-quality plot of the NLE curve."""
        plt.style.use("default")
        fig, ax = plt.subplots(figsize=(7, 6))

        # Ideal linear reference line
        ax.plot(
            ee_cat * 100, ee_cat * 100,
            color="#7f8c8d", linestyle="--", linewidth=1.5,
            label="Ideal Linear Effect (No Aggregation)",
        )

        # Predicted NLE curve
        ax.plot(
            ee_cat * 100, ee_prod * 100,
            color="#004488", linestyle="-", linewidth=2.5,
            label="Predicted NLE (ML$_2$ Reservoir)",
        )

        ax.set_title("Predicted Non-Linear Effect", fontsize=15, fontweight="bold", pad=12)
        ax.set_xlabel("Catalyst $ee$ (%)", fontsize=13)
        ax.set_ylabel("Product $ee$ (%)", fontsize=13)

        ax.tick_params(axis="both", which="major", labelsize=11)
        ax.grid(True, linestyle=":", alpha=0.6, color="gray")

        ax.set_xlim(0, 100)
        ax.set_ylim(0, 100)
        ax.set_aspect("equal", adjustable="box")

        ax.legend(fontsize=11, loc="lower right", framealpha=0.95, edgecolor="gray")

        plt.tight_layout()
        plt.savefig(save_path, dpi=600, bbox_inches="tight")
        plt.close(fig)
        logger.info(f"Publication-ready plot saved to {save_path}")


def write_sensitivity_report(
    output_dir: str,
    ee_cat: np.ndarray,
    delta_g_kcal: float,
    temp: float,
    k_homo: float,
) -> None:
    """Write a surrogate sensitivity table around the baseline prediction."""
    scenarios = [
        ("baseline", delta_g_kcal, temp, k_homo),
        ("delta_g_minus_0.5_kcal", delta_g_kcal - 0.5, temp, k_homo),
        ("delta_g_plus_0.5_kcal", delta_g_kcal + 0.5, temp, k_homo),
        ("temperature_minus_20K", delta_g_kcal, temp - 20, k_homo),
        ("temperature_plus_20K", delta_g_kcal, temp + 20, k_homo),
        ("k_homo_x0.1", delta_g_kcal, temp, k_homo * 0.1),
        ("k_homo_x10", delta_g_kcal, temp, k_homo * 10),
    ]
    report_path = os.path.join(output_dir, "nle_sensitivity.csv")
    with open(report_path, "w", newline="", encoding="utf-8") as sensitivity_file:
        writer = csv.writer(sensitivity_file)
        writer.writerow([
            "scenario",
            "delta_g_kcal_per_mol",
            "temperature_K",
            "K_homo",
            "product_ee_at_50_percent_catalyst_ee",
            "maximum_product_ee",
        ])
        for name, scenario_delta_g, scenario_temp, scenario_k_homo in scenarios:
            scenario_calculator = NLECalculator(temp=scenario_temp)
            product_ee = scenario_calculator.simulate_kagan_model(
                ee_cat,
                scenario_delta_g,
                k_homo_abs=scenario_k_homo,
            )
            writer.writerow([
                name,
                f"{scenario_delta_g:.6f}",
                f"{scenario_temp:.2f}",
                f"{scenario_k_homo:.6g}",
                f"{np.interp(0.5, ee_cat, product_ee) * 100:.4f}",
                f"{np.max(product_ee) * 100:.4f}",
            ])
    logger.info(f"Surrogate sensitivity report saved to {report_path}")
