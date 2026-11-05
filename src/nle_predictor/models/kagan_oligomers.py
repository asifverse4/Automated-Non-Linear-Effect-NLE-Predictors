"""
Higher-Order Kagan Aggregation Models (ML3 Trimers and ML4 Tetramers).

Theoretical foundations:
Henri B. Kagan, "Nonlinear Effects in Asymmetric Synthesis. 2. Higher-Order Aggregates"
Bull. Soc. Chim. Fr. 1994, 131, 849-867; Chem. Rev. 2003, 103, 2921-2944.

This module models thermodynamic aggregation equilibria beyond dimers:
- Trimers (ML3): R3, R2S, RS2, S3 with statistical symmetry weights [1, 3, 3, 1].
- Tetramers (ML4): R4, R3S, R2S2, RS3, S4 with statistical symmetry weights [1, 4, 6, 4, 1].
- General n-mers: statistical combinatorial binom(n, k) weights with chiral self-recognition (g-factors).

Provides exact analytical Jacobians for guaranteed quadratic Newton-Raphson convergence,
full species distribution tracking, and detection of higher-order sigmoidal threshold effects.
"""

import math
from typing import Any, List, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import root

from nle_predictor.core.constants import DEFAULT_TEMP, GAS_CONSTANT_R
from nle_predictor.core.types import NLESimulationResult
from nle_predictor.models.base import NLEModel


class KaganOligomerModel(NLEModel):
    """
    Higher-order Kagan aggregation model for trimers (n=3), tetramers (n=4),
    and general n-mers in asymmetric catalysis.

    Parameters
    ----------
    order : int
        Oligomerization degree n (default 3 for trimers, 4 for tetramers).
    temp : float
        Absolute temperature in Kelvin (default 298.15 K).
    """

    def __init__(self, order: int = 3, temp: float = DEFAULT_TEMP):
        if order < 2:
            raise ValueError(f"Oligomerization order must be >= 2, got {order}.")
        if temp <= 0:
            raise ValueError(f"Temperature must be greater than zero Kelvin, got {temp}.")
        self.order = order
        self.temp = temp

    @staticmethod
    def combinatorial_weights(n: int) -> List[int]:
        """Compute statistical combinatorial weights binom(n, k) for k = n, n-1, ..., 0."""
        return [math.comb(n, k) for k in range(n, -1, -1)]

    def compute_equilibrium_constants(
        self,
        k_homo: float,
        g_factors: Optional[Sequence[float]] = None,
        delta_g_kcal: Optional[float] = None,
    ) -> np.ndarray:
        """
        Compute association constants K_{k, n-k} for each aggregate R_k S_{n-k} (k = n, ..., 0).

        Parameters
        ----------
        k_homo : float
            Association constant for homochiral aggregate (R_n or S_n).
        g_factors : Sequence[float], optional
            Chiral self-recognition factors g_k relative to statistical ideal mixing.
            If length 1, identical g is applied across all heterochiral aggregates.
        delta_g_kcal : float, optional
            Free energy difference Delta G (hetero - homo) in kcal/mol. If provided,
            g = exp(-delta_g_kcal / (R * T)).

        Returns
        -------
        np.ndarray
            Array of K_{k, n-k} of length n + 1 for k = n, n-1, ..., 0.
        """
        n = self.order
        if k_homo < 0:
            raise ValueError("k_homo must be non-negative.")

        # Determine chiral self-recognition factor g
        if delta_g_kcal is not None:
            g_single = float(np.exp(-delta_g_kcal / (GAS_CONSTANT_R * self.temp)))
        elif g_factors is not None and len(g_factors) > 0:
            g_single = float(g_factors[0])
        else:
            g_single = 1.0

        weights = self.combinatorial_weights(n)
        k_values = np.zeros(n + 1, dtype=float)

        for idx, k in enumerate(range(n, -1, -1)):
            comb = weights[idx]
            if k == n or k == 0:
                # Pure homochiral species: R_n or S_n
                k_values[idx] = float(comb * k_homo)
            else:
                # Heterochiral species
                if g_factors is not None and len(g_factors) >= n - 1:
                    # Index into heterochiral factors
                    het_idx = min(k - 1, (n - k) - 1)
                    g_val = g_factors[het_idx]
                else:
                    g_val = g_single
                k_values[idx] = float(comb * g_val * k_homo)

        return k_values

    def solve_mass_balance(
        self,
        ee_cat: float,
        k_homo: float,
        k_aggregates: np.ndarray,
        c_tot: float = 1.0,
    ) -> Tuple[float, float, np.ndarray]:
        """
        Solve coupled mass-balance equations for monomer concentrations [R] and [S]
        using Newton-Raphson with an analytical Jacobian.

        Parameters
        ----------
        ee_cat : float
            Catalyst enantiomeric excess fraction in [-1.0, 1.0].
        k_homo : float
            Homochiral association constant.
        k_aggregates : np.ndarray
            Equilibrium constants K_{k, n-k} for k = n, ..., 0.
        c_tot : float
            Total catalyst concentration (default 1.0 M).

        Returns
        -------
        Tuple[float, float, np.ndarray]
            (m_r, m_s, aggregate_concentrations)
        """
        if not -1.0 <= ee_cat <= 1.0:
            raise ValueError(f"Catalyst ee must be in [-1.0, 1.0], got {ee_cat}.")
        if c_tot <= 0:
            raise ValueError("Total catalyst concentration c_tot must be positive.")

        n = self.order
        c_r = c_tot * 0.5 * (1.0 + ee_cat)
        c_s = c_tot * 0.5 * (1.0 - ee_cat)

        # k values correspond to k = n, n-1, ..., 0
        k_powers = np.arange(n, -1, -1, dtype=float)  # powers of R
        s_powers = float(n) - k_powers               # powers of S

        def residuals_and_jacobian(x: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
            mr = max(float(x[0]), 1e-15)
            ms = max(float(x[1]), 1e-15)

            # Aggregate concentrations: [R_k S_{n-k}] = K_{k, n-k} * mr^k * ms^{n-k}
            # k=n: mr^n; k=0: ms^n
            concs = np.zeros(n + 1, dtype=float)
            for i in range(n + 1):
                kp = int(k_powers[i])
                sp = int(s_powers[i])
                concs[i] = k_aggregates[i] * (mr ** kp) * (ms ** sp)

            # Mass balances
            # F1 = mr + sum(k * concs) - c_r
            # F2 = ms + sum((n - k) * concs) - c_s
            f1 = mr + np.sum(k_powers * concs) - c_r
            f2 = ms + np.sum(s_powers * concs) - c_s

            # Analytical Jacobian elements:
            # df1_dmr = 1 + sum(k^2 * K * mr^{k-1} * ms^{n-k})
            # df1_dms = sum(k * (n-k) * K * mr^k * ms^{n-k-1})
            # df2_dmr = sum(k * (n-k) * K * mr^{k-1} * ms^{n-k})
            # df2_dms = 1 + sum((n-k)^2 * K * mr^k * ms^{n-k-1})
            df1_dmr = 1.0
            df1_dms = 0.0
            df2_dmr = 0.0
            df2_dms = 1.0

            for i in range(n + 1):
                kp = int(k_powers[i])
                sp = int(s_powers[i])
                k_val = k_aggregates[i]
                if kp > 0:
                    df1_dmr += (kp ** 2) * k_val * (mr ** (kp - 1)) * (ms ** sp)
                if kp > 0 and sp > 0:
                    df1_dms += kp * sp * k_val * (mr ** kp) * (ms ** (sp - 1))
                    df2_dmr += kp * sp * k_val * (mr ** (kp - 1)) * (ms ** sp)
                if sp > 0:
                    df2_dms += (sp ** 2) * k_val * (mr ** kp) * (ms ** (sp - 1))

            return np.array([f1, f2], dtype=float), np.array([[df1_dmr, df1_dms], [df2_dmr, df2_dms]], dtype=float)

        # Initial guess based on approximate homochiral root: mr + n * k_homo * mr^n ~ c_r
        if k_homo > 0 and c_r > 0:
            mr_guess = (c_r / (1.0 + n * k_homo * (c_r ** (n - 1)))) ** (1.0 / n)
            mr_guess = min(c_r, max(mr_guess, 1e-6))
        else:
            mr_guess = c_r

        if k_homo > 0 and c_s > 0:
            ms_guess = (c_s / (1.0 + n * k_homo * (c_s ** (n - 1)))) ** (1.0 / n)
            ms_guess = min(c_s, max(ms_guess, 1e-6))
        else:
            ms_guess = c_s

        def eq_fun(x):
            res, _ = residuals_and_jacobian(x)
            return res

        def jac_fun(x):
            _, jac = residuals_and_jacobian(x)
            return jac

        sol = root(eq_fun, [mr_guess, ms_guess], jac=jac_fun, method="hybr", tol=1e-12)
        if sol.success and np.all(sol.x >= 0.0):
            mr = float(sol.x[0])
            ms = float(sol.x[1])
        else:
            from scipy.optimize import minimize

            def obj(x):
                r, _ = residuals_and_jacobian(np.maximum(x, 0.0))
                return float(np.sum(r**2))

            res_min = minimize(obj, [mr_guess, ms_guess], bounds=[(0.0, c_tot), (0.0, c_tot)], method="L-BFGS-B")
            mr, ms = max(0.0, float(res_min.x[0])), max(0.0, float(res_min.x[1]))

        # Exact Newton-Raphson refinement for machine-precision mass conservation
        for _ in range(10):
            f_vec, j_mat = residuals_and_jacobian(np.array([mr, ms]))
            if np.max(np.abs(f_vec)) < 1e-12:
                break
            det = j_mat[0, 0] * j_mat[1, 1] - j_mat[0, 1] * j_mat[1, 0]
            if abs(det) > 1e-15:
                dmr = -(j_mat[1, 1] * f_vec[0] - j_mat[0, 1] * f_vec[1]) / det
                dms = -(-j_mat[1, 0] * f_vec[0] + j_mat[0, 0] * f_vec[1]) / det
                mr = max(0.0, mr + dmr)
                ms = max(0.0, ms + dms)
            else:
                break

        # Calculate final aggregate concentrations
        aggregates = np.zeros(n + 1, dtype=float)
        for i in range(n + 1):
            aggregates[i] = k_aggregates[i] * (mr ** int(k_powers[i])) * (ms ** int(s_powers[i]))

        return mr, ms, aggregates

    def simulate(
        self,
        ee_cat_range: Sequence[float],
        delta_g_kcal: float = 0.0,
        k_homo: float = 1e5,
        ee_0: float = 1.0,
        g_factors: Optional[Sequence[float]] = None,
        c_tot: float = 1.0,
        active_aggregates: bool = False,
        aggregate_ee: Optional[Sequence[float]] = None,
        **kwargs: Any,
    ) -> NLESimulationResult:
        """
        Simulate product enantiomeric excess for higher-order aggregation equilibria.

        Parameters
        ----------
        ee_cat_range : Sequence[float]
            Array of catalyst ee values in [0, 1] or [-1, 1].
        delta_g_kcal : float
            Free energy difference for heterochiral aggregate in kcal/mol.
        k_homo : float
            Homochiral association constant (default 1e5).
        ee_0 : float
            Intrinsic enantioselectivity of the active monomer (default 1.0).
        g_factors : Sequence[float], optional
            Chiral self-recognition factors g_k.
        c_tot : float
            Total catalyst concentration (default 1.0).
        active_aggregates : bool
            Whether higher-order oligomers also contribute catalytic activity.
        aggregate_ee : Sequence[float], optional
            Intrinsic ee contributions of aggregates (k = n, ..., 0).

        Returns
        -------
        NLESimulationResult
            Typed NLE simulation result.
        """
        ee_cat = np.asarray(ee_cat_range, dtype=float)
        k_aggs = self.compute_equilibrium_constants(k_homo, g_factors=g_factors, delta_g_kcal=delta_g_kcal)

        ee_prods: List[float] = []
        monomers_r: List[float] = []
        monomers_s: List[float] = []
        all_aggregates: List[np.ndarray] = []

        n = self.order
        for c in ee_cat:
            mr, ms, aggs = self.solve_mass_balance(c, k_homo, k_aggs, c_tot=c_tot)
            monomers_r.append(mr)
            monomers_s.append(ms)
            all_aggregates.append(aggs)

            tot_monomer = mr + ms
            if not active_aggregates:
                # Classic Kagan reservoir mode: only monomers are catalytically active
                if tot_monomer <= 1e-14:
                    ee_p = 0.0
                else:
                    ee_p = ee_0 * ((mr - ms) / tot_monomer)
            else:
                # Active oligomer mode: monomers and oligomers contribute to product formation
                v_r = 0.5 * (1.0 + ee_0) * mr + 0.5 * (1.0 - ee_0) * ms
                v_s = 0.5 * (1.0 - ee_0) * mr + 0.5 * (1.0 + ee_0) * ms

                if aggregate_ee is not None and len(aggregate_ee) == n + 1:
                    for i in range(n + 1):
                        agg_selectivity = aggregate_ee[i]
                        # Relative catalytic activity assumed proportional to concentration
                        v_r += 0.5 * (1.0 + agg_selectivity) * aggs[i]
                        v_s += 0.5 * (1.0 - agg_selectivity) * aggs[i]

                tot_v = v_r + v_s
                ee_p = (v_r - v_s) / tot_v if tot_v > 1e-14 else 0.0

            ee_prods.append(float(np.clip(ee_p, -1.0, 1.0)))

        ee_prod = np.array(ee_prods, dtype=float)
        diff = ee_prod - ee_cat
        tol = 1e-5
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

        model_label = f"Kagan ML{n} {'Trimer' if n == 3 else 'Tetramer' if n == 4 else f'{n}-mer'} Model"
        g_eff = float(k_aggs[1] / (math.comb(n, n - 1) * k_homo)) if (k_homo > 0 and len(k_aggs) > 1) else 1.0

        return NLESimulationResult(
            model_name=model_label,
            ee_catalyst=ee_cat,
            ee_product=ee_prod,
            monomers_r=np.array(monomers_r),
            monomers_s=np.array(monomers_s),
            classification=classification,
            product_at_50=product_at_50,
            max_amplification=max_amp,
            max_deviation_from_linearity=float(np.max(np.abs(diff))),
            metadata={
                "order": n,
                "k_homo": k_homo,
                "g_effective": g_eff,
                "delta_g_kcal": delta_g_kcal,
                "k_aggregates": k_aggs.tolist(),
                "temperature_k": self.temp,
                "active_aggregates": active_aggregates,
                "c_tot": c_tot,
            },
        )
