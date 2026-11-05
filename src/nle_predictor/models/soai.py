"""
Asymmetric Autocatalysis and Soai Reaction Kinetic Simulation Engine.

Theoretical and experimental foundations:
- K. Soai, T. Shibata, H. Morioka, K. Choji, Nature 1995, 378, 767-768.
- D. G. Blackmond, C. R. McMillan, S. Ramdeehul, A. Schorm, J. M. Brown,
  J. Am. Chem. Soc. 2001, 123, 545-555.
- D. G. Blackmond, Proc. Natl. Acad. Sci. USA 2004, 101, 5732-5736.

This module models the kinetic amplification of enantiomeric excess in asymmetric autocatalysis,
featuring:
1. Reversible homochiral and heterochiral product oligomerization/dimerization (mutual inhibition).
2. Square-order and general nth-order autocatalytic propagation.
3. Full numerical ODE integration of temporal trajectories via scipy.integrate.solve_ivp.
4. Hyper-amplification from trace seed enantiomeric excess (ee_seed ~ 0.0001) to near enantiopurity.
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Sequence, Tuple

import numpy as np
from scipy.integrate import solve_ivp

from nle_predictor.core.types import NLESimulationResult
from nle_predictor.models.base import NLEModel


@dataclass
class SoaiTrajectory:
    """Time-resolved concentration and enantioselectivity trajectory of a Soai reaction."""
    time: np.ndarray
    substrate_a: np.ndarray
    substrate_b: np.ndarray
    product_r: np.ndarray
    product_s: np.ndarray
    ee_cumulative: np.ndarray
    ee_instantaneous: np.ndarray
    conversion: np.ndarray
    final_ee: float
    amplification_factor: float
    half_life_time: float
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "final_ee_percent": float(self.final_ee * 100),
            "amplification_factor": float(self.amplification_factor),
            "half_life_time_seconds": float(self.half_life_time),
            "time_points": len(self.time),
            "metadata": self.metadata,
        }


class SoaiAutocatalysisModel(NLEModel):
    """
    Kinetic simulation engine for asymmetric autocatalysis with mutual inhibition.

    Parameters
    ----------
    k_uncat : float
        Second-order rate constant for uncatalyzed background reaction A + B -> P (M^-1 s^-1).
    k_auto : float
        Rate constant for autocatalytic reaction (M^(-1-order) s^-1).
    auto_order : float
        Reaction order in catalyst product (typically 2.0 for dimer autocatalysis).
    k_hetero : float
        Association constant for inactive heterochiral dimer PR + PS <-> (PR PS).
    k_homo : float
        Association constant for homochiral dimers (PR)2 and (PS)2.
    ee_0 : float
        Intrinsic enantioselectivity of the active catalyst (default 1.0 = 100%).
    """

    def __init__(
        self,
        k_uncat: float = 1e-6,
        k_auto: float = 1e4,
        auto_order: float = 2.0,
        k_hetero: float = 1e4,
        k_homo: float = 10.0,
        ee_0: float = 1.0,
    ):
        if k_uncat < 0 or k_auto < 0:
            raise ValueError("Rate constants must be non-negative.")
        if auto_order <= 0:
            raise ValueError("Autocatalytic order must be positive.")
        if k_hetero < 0 or k_homo < 0:
            raise ValueError("Dimerization constants must be non-negative.")
        if not 0.0 <= ee_0 <= 1.0:
            raise ValueError("Intrinsic ee_0 must be in [0.0, 1.0].")

        self.k_uncat = k_uncat
        self.k_auto = k_auto
        self.auto_order = auto_order
        self.k_hetero = k_hetero
        self.k_homo = k_homo
        self.ee_0 = ee_0

    def solve_dimer_equilibrium(
        self,
        tot_r: float,
        tot_s: float,
    ) -> Tuple[float, float, float, float, float]:
        """
        Calculates free monomer concentrations [PR]_free and [PS]_free
        in rapid equilibrium with homochiral and heterochiral dimers:
        tot_r = [PR] + 2 [(PR)2] + [(PR PS)]
        tot_s = [PS] + 2 [(PS)2] + [(PR PS)]
        where [(PR)2] = K_homo [PR]^2, [(PS)2] = K_homo [PS]^2, [(PR PS)] = K_hetero [PR] [PS].
        """
        if self.k_hetero == 0.0 and self.k_homo == 0.0:
            return tot_r, tot_s, 0.0, 0.0, 0.0

        if tot_r <= 1e-15 and tot_s <= 1e-15:
            return 0.0, 0.0, 0.0, 0.0, 0.0

        kh = self.k_homo
        khet = self.k_hetero

        # Iterative Newton-Raphson for rapid solution inside ODE integrator
        mr = max(tot_r / (1.0 + 2.0 * kh * tot_r + khet * tot_s), 1e-15)
        ms = max(tot_s / (1.0 + 2.0 * kh * tot_s + khet * tot_r), 1e-15)

        for _ in range(8):
            f1 = mr + 2.0 * kh * (mr ** 2) + khet * mr * ms - tot_r
            f2 = ms + 2.0 * kh * (ms ** 2) + khet * mr * ms - tot_s

            if abs(f1) < 1e-12 and abs(f2) < 1e-12:
                break

            j11 = 1.0 + 4.0 * kh * mr + khet * ms
            j12 = khet * mr
            j21 = khet * ms
            j22 = 1.0 + 4.0 * kh * ms + khet * mr

            det = j11 * j22 - j12 * j21
            if abs(det) > 1e-16:
                mr = max(0.0, mr - (j22 * f1 - j12 * f2) / det)
                ms = max(0.0, ms - (-j21 * f1 + j11 * f2) / det)
            else:
                break

        dimer_rr = kh * (mr ** 2)
        dimer_ss = kh * (ms ** 2)
        dimer_rs = khet * mr * ms

        return mr, ms, dimer_rr, dimer_ss, dimer_rs

    def simulate_trajectory(
        self,
        init_a: float = 0.1,
        init_b: float = 0.1,
        seed_catalyst: float = 0.001,
        seed_ee: float = 0.01,
        t_span: Tuple[float, float] = (0.0, 3600.0),
        n_points: int = 500,
    ) -> SoaiTrajectory:
        """
        Simulate the full temporal kinetic trajectory of a single asymmetric autocatalysis reaction.

        Parameters
        ----------
        init_a : float
            Initial concentration of pyrimidine-5-carbaldehyde substrate [A]0 in M (default 0.1 M).
        init_b : float
            Initial concentration of dialkylzinc reagent [B]0 in M (default 0.1 M).
        seed_catalyst : float
            Total initial concentration of chiral seed product [P]0 in M (default 1e-3 M).
        seed_ee : float
            Enantiomeric excess of chiral seed product in [-1.0, 1.0] (default 0.01 = 1% ee).
        t_span : Tuple[float, float]
            Time interval for simulation in seconds (default 0 to 3600 s).
        n_points : int
            Number of evaluation points.

        Returns
        -------
        SoaiTrajectory
            Full time-resolved trajectory with concentrations, ee, conversion, and amplification.
        """
        if init_a <= 0 or init_b <= 0:
            raise ValueError("Substrate concentrations must be positive.")
        if seed_catalyst < 0:
            raise ValueError("Seed catalyst concentration cannot be negative.")
        if not -1.0 <= seed_ee <= 1.0:
            raise ValueError("Seed ee must be in [-1.0, 1.0].")

        # Initial concentrations of R and S product enantiomers
        p_r_0 = seed_catalyst * 0.5 * (1.0 + seed_ee)
        p_s_0 = seed_catalyst * 0.5 * (1.0 - seed_ee)

        # State vector: y = [A, B, P_R, P_S]
        y0 = [init_a, init_b, p_r_0, p_s_0]

        p = self.auto_order
        k0 = self.k_uncat
        ka = self.k_auto
        ee0 = self.ee_0

        def ode_system(t: float, y: np.ndarray) -> List[float]:
            a = max(0.0, float(y[0]))
            b = max(0.0, float(y[1]))
            tot_r = max(0.0, float(y[2]))
            tot_s = max(0.0, float(y[3]))

            # Rapid monomer-dimer speciation
            mr, ms, _, _, _ = self.solve_dimer_equilibrium(tot_r, tot_s)

            # Uncatalyzed background reaction (produces racemic mixture)
            v_uncat = k0 * a * b

            # Autocatalytic reactions driven by active monomers (or dimers)
            # v_R = ka * [A] * [B] * [PR]^p
            # v_S = ka * [A] * [B] * [PS]^p
            v_auto_r = ka * a * b * (mr ** p)
            v_auto_s = ka * a * b * (ms ** p)

            v_tot = v_uncat + v_auto_r + v_auto_s

            da_dt = -v_tot
            db_dt = -v_tot

            # Enantiomeric product formation
            # Uncatalyzed creates 50/50 R/S
            # Catalyzed by R produces R with selectivity ee0
            # Catalyzed by S produces S with selectivity ee0
            dp_r_dt = (
                0.5 * v_uncat
                + 0.5 * (1.0 + ee0) * v_auto_r
                + 0.5 * (1.0 - ee0) * v_auto_s
            )
            dp_s_dt = (
                0.5 * v_uncat
                + 0.5 * (1.0 - ee0) * v_auto_r
                + 0.5 * (1.0 + ee0) * v_auto_s
            )

            return [da_dt, db_dt, dp_r_dt, dp_s_dt]

        t_eval = np.linspace(t_span[0], t_span[1], n_points)
        sol = solve_ivp(
            ode_system,
            t_span,
            y0,
            t_eval=t_eval,
            method="Radau",
            rtol=1e-7,
            atol=1e-9,
        )

        if not sol.success:
            # Fallback to BDF or RK45 if Radau encountered a stiff boundary
            sol = solve_ivp(ode_system, t_span, y0, t_eval=t_eval, method="BDF", rtol=1e-6, atol=1e-8)

        t_arr = sol.t
        a_arr = np.maximum(0.0, sol.y[0])
        b_arr = np.maximum(0.0, sol.y[1])
        r_arr = np.maximum(0.0, sol.y[2])
        s_arr = np.maximum(0.0, sol.y[3])

        tot_p = r_arr + s_arr
        ee_cum = np.where(tot_p > 1e-12, (r_arr - s_arr) / tot_p, seed_ee)

        # Numerical derivative for instantaneous ee
        dr_dt = np.gradient(r_arr, t_arr)
        ds_dt = np.gradient(s_arr, t_arr)
        tot_rate = dr_dt + ds_dt
        ee_inst = np.where(tot_rate > 1e-12, (dr_dt - ds_dt) / tot_rate, seed_ee)

        conversion = 1.0 - (a_arr / init_a)
        final_ee = float(ee_cum[-1])
        amplification = (final_ee / seed_ee) if abs(seed_ee) > 1e-6 else 1.0

        # Estimate half-life (time to 50% conversion)
        idx_half = np.where(conversion >= 0.5)[0]
        half_life = float(t_arr[idx_half[0]]) if len(idx_half) > 0 else float(t_span[1])

        return SoaiTrajectory(
            time=t_arr,
            substrate_a=a_arr,
            substrate_b=b_arr,
            product_r=r_arr,
            product_s=s_arr,
            ee_cumulative=ee_cum,
            ee_instantaneous=ee_inst,
            conversion=conversion,
            final_ee=final_ee,
            amplification_factor=float(amplification),
            half_life_time=half_life,
            metadata={
                "k_uncat": self.k_uncat,
                "k_auto": self.k_auto,
                "auto_order": self.auto_order,
                "k_hetero": self.k_hetero,
                "k_homo": self.k_homo,
                "seed_catalyst_m": seed_catalyst,
                "seed_ee": seed_ee,
            },
        )

    def simulate(
        self,
        ee_cat_range: Sequence[float],
        delta_g_kcal: float = 0.0,
        k_homo: float = 1e5,
        ee_0: float = 1.0,
        seed_catalyst: float = 0.001,
        t_final: float = 3600.0,
        **kwargs: Any,
    ) -> NLESimulationResult:
        """
        Simulate asymmetric autocatalysis NLE curve across a range of seed catalyst enantiomeric excesses.

        Parameters
        ----------
        ee_cat_range : Sequence[float]
            Array of seed enantiomeric excess values.
        delta_g_kcal : float
            Free energy difference of dimerization (kcal/mol), overrides k_hetero if provided.
        k_homo : float
            Homochiral dimerization constant.
        ee_0 : float
            Intrinsic enantioselectivity.
        seed_catalyst : float
            Initial seed concentration (default 1e-3 M).
        t_final : float
            Final reaction time for each kinetic batch run (default 3600 s).

        Returns
        -------
        NLESimulationResult
            Typed simulation result containing NLE curve and diagnostics.
        """
        ee_cat = np.asarray(ee_cat_range, dtype=float)

        if delta_g_kcal != 0.0:
            # Overwrite k_hetero based on Delta G
            from nle_predictor.core.constants import GAS_CONSTANT_R
            self.k_homo = k_homo
            self.k_hetero = 2.0 * k_homo * float(np.exp(-delta_g_kcal / (GAS_CONSTANT_R * 298.15)))

        if ee_0 is not None:
            self.ee_0 = ee_0

        ee_prods: List[float] = []
        monomers_r: List[float] = []
        monomers_s: List[float] = []

        for c in ee_cat:
            traj = self.simulate_trajectory(
                seed_catalyst=seed_catalyst,
                seed_ee=c,
                t_span=(0.0, t_final),
                n_points=100,
            )
            ee_prods.append(traj.final_ee)
            monomers_r.append(float(traj.product_r[-1]))
            monomers_s.append(float(traj.product_s[-1]))

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

        return NLESimulationResult(
            model_name="Soai Asymmetric Autocatalysis Model",
            ee_catalyst=ee_cat,
            ee_product=ee_prod,
            monomers_r=np.array(monomers_r),
            monomers_s=np.array(monomers_s),
            classification=classification,
            product_at_50=product_at_50,
            max_amplification=max_amp,
            max_deviation_from_linearity=float(np.max(np.abs(diff))),
            metadata={
                "k_uncat": self.k_uncat,
                "k_auto": self.k_auto,
                "auto_order": self.auto_order,
                "k_hetero": self.k_hetero,
                "k_homo": self.k_homo,
                "seed_catalyst": seed_catalyst,
                "t_final_seconds": t_final,
            },
        )
