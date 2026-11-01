"""
Abstract base class and interface for Non-Linear Effect (NLE) mathematical models.
"""

from abc import ABC, abstractmethod
from typing import Optional, Sequence
import numpy as np

from nle_predictor.core.types import NLESimulationResult


class NLEModel(ABC):
    """Abstract base class for asymmetric catalysis non-linear effect models."""

    @abstractmethod
    def simulate(
        self,
        ee_cat_range: Sequence[float],
        delta_g_kcal: float,
        k_homo: float = 1e5,
        **kwargs,
    ) -> NLESimulationResult:
        """
        Simulate product enantiomeric excess as a function of catalyst enantiomeric excess.

        Parameters
        ----------
        ee_cat_range : Sequence[float]
            Array of catalyst ee values (in fraction, e.g. from 0.0 to 1.0).
        delta_g_kcal : float
            Free energy difference Delta G (RS - RR) in kcal/mol.
        k_homo : float
            Absolute homochiral dimerization association constant.
        **kwargs : dict
            Model-specific parameters.

        Returns
        -------
        NLESimulationResult
            Detailed simulation result with monomers, product ee, and diagnostics.
        """
        pass

