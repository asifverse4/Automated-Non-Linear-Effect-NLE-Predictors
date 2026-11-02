"""
Noyori DAIB Aggregation Model for Asymmetric Dialkylzinc Additions.
References:
1. M. Kitamura, S. Okada, S. Suga, R. Noyori, J. Am. Chem. Soc. 1989, 111, 4028-4036.
2. R. Noyori, Asymmetric Catalysis in Organic Synthesis, Wiley-Interscience, 1994.
"""

import logging
from typing import Optional, Sequence, Tuple
import numpy as np

from nle_predictor.core.constants import DEFAULT_TEMP, GAS_CONSTANT_R
from nle_predictor.core.types import NLESimulationResult
from nle_predictor.models.base import NLEModel
from nle_predictor.models.kagan_ml2 import NLECalculator

logger = logging.getLogger("NLE_Predictor")


class NoyoriDAIBModel(NLEModel):
    """
    Implements the Noyori DAIB catalyst aggregation model.

    In the Noyori mechanism, homochiral dimers (S,S) and (R,R) exist in equilibrium
    with monomeric alkylzinc complexes, while the heterochiral dimer (S,R)
    functions as an exceptionally stable, catalytically inert reservoir trap.

    Parameters
    ----------
    temp : float
        Reaction temperature in Kelvin.
    ee_0 : float
        Intrinsic enantioselectivity of the active monomeric catalyst (default 0.98 for DAIB).
    kd_homo : float
        Dimer dissociation constant for homochiral species: K_d(homo) = [M]^2 / [D_homo].
    kd_hetero : float
        Dimer dissociation constant for heterochiral species: K_d(hetero) = [M_S][M_R] / [D_hetero].
    """
    def __init__(
        self,
        temp: float = DEFAULT_TEMP,
        ee_0: float = 0.98,
        kd_homo: Optional[float] = 1e-3,
        kd_hetero: Optional[float] = 1e-6,
    ):
        if temp <= 0:
            raise ValueError("Temperature must be greater than zero Kelvin.")
        if not (-1.0 <= ee_0 <= 1.0):
            raise ValueError("Intrinsic enantioselectivity ee_0 must be in [-1.0, 1.0].")
        if kd_homo <= 0 or kd_hetero <= 0:
            raise ValueError("Dissociation constants must be strictly positive.")

        self.temp = temp
        self.ee_0 = ee_0
        self.kd_homo = kd_homo
        self.kd_hetero = kd_hetero
        self._calculator = NLECalculator(temp=temp)

    @property
    def k_homo_assoc(self) -> float:
        """Homochiral association constant K_assoc = 1 / K_d(homo)."""
        return 1.0 / self.kd_homo

    @property
    def k_hetero_assoc(self) -> float:
        """Heterochiral association constant K_assoc = 1 / K_d(hetero)."""
        return 1.0 / self.kd_hetero

    @property
    def delta_g_kcal(self) -> float:
        """
        Thermodynamic stability difference Delta G (RS - RR) implied by dissociation constants.
        K_hetero / K_homo = 2 * exp(-Delta G / RT)
        Delta G = -RT * ln((K_hetero / K_homo) / 2)
        """
        ratio = self.k_hetero_assoc / self.k_homo_assoc
        return float(-GAS_CONSTANT_R * self.temp * np.log(ratio / 2.0))

    def simulate(
        self,
        ee_cat_range: Sequence[float],
        delta_g_kcal: Optional[float] = None,
        k_homo: Optional[float] = None,
        ee_0: Optional[float] = None,
        **kwargs,
    ) -> NLESimulationResult:
        """
        Simulate the Noyori NLE curve across catalyst optical purities.

        Parameters
        ----------
        ee_cat_range : Sequence[float]
            Array of catalyst ee values (fraction in [0, 1]).
        delta_g_kcal : Optional[float]
            If provided, overrides kd_hetero using the thermodynamic relation.
        k_homo : Optional[float]
            If provided, overrides kd_homo (k_homo = 1 / kd_homo).
        ee_0 : Optional[float]
            Override intrinsic monomer selectivity.
        """
        ee_cat = np.asarray(ee_cat_range, dtype=float)
        ee0_val = self.ee_0 if ee_0 is None else ee_0

        if k_homo is not None:
            kh = k_homo
        else:
            kh = self.k_homo_assoc

        if delta_g_kcal is not None:
            dg = delta_g_kcal
            ratio = 2.0 * np.exp(-dg / (GAS_CONSTANT_R * self.temp))
            k_hetero = kh * ratio
        else:
            k_hetero = self.k_hetero_assoc
            dg = self.delta_g_kcal
            ratio = k_hetero / kh

        ee_prods = []
        monomers_r = []
        monomers_s = []

        for cat in ee_cat:
            mr, ms = self._calculator.solve_mass_balance(cat, kh, k_hetero)
            monomers_r.append(mr)
            monomers_s.append(ms)

            active_tot = mr + ms
            if active_tot <= 1e-14:
                ee_p = 0.0
            else:
                ee_p = ee0_val * ((mr - ms) / active_tot)
            ee_prods.append(ee_p)

        ee_prod_arr = np.array(ee_prods)
        diff = ee_prod_arr - ee_cat
        product_at_50 = float(np.interp(0.5, ee_cat, ee_prod_arr))

        amplification = np.divide(
            ee_prod_arr[1:], ee_cat[1:], out=np.zeros_like(ee_prod_arr[1:]), where=ee_cat[1:] > 0
        )
        max_amp = float(np.max(amplification)) if len(amplification) > 0 else 1.0

        tol = 1e-5
        if np.max(diff) > tol:
            classification = "positive NLE (Noyori asymmetric amplification)"
        elif np.min(diff) < -tol:
            classification = "negative NLE"
        else:
            classification = "near-linear"

        return NLESimulationResult(
            model_name="Noyori DAIB Aggregation Model",
            ee_catalyst=ee_cat,
            ee_product=ee_prod_arr,
            monomers_r=np.array(monomers_r),
            monomers_s=np.array(monomers_s),
            classification=classification,
            product_at_50=product_at_50,
            max_amplification=max_amp,
            max_deviation_from_linearity=float(np.max(np.abs(diff))),
            metadata={
                "k_homo_assoc": kh,
                "k_hetero_assoc": k_hetero,
                "kd_homo": 1.0 / kh,
                "kd_hetero": 1.0 / k_hetero,
                "k_ratio": ratio,
                "delta_g_kcal": dg,
                "ee_0": ee0_val,
                "temperature_k": self.temp,
            },
        )
