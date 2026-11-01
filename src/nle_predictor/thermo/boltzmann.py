"""
Boltzmann ensemble thermodynamic weighting and conformational entropy calculations.
"""

from typing import List, Sequence, Union
import numpy as np

from nle_predictor.core.constants import GAS_CONSTANT_R, HARTREE_TO_KCAL
from nle_predictor.core.types import Conformer, EnsembleResult


def boltzmann_weights(
    energies_hartree: Sequence[float],
    temperature: float,
) -> np.ndarray:
    """
    Calculate the fractional Boltzmann population distribution for an ensemble of conformers.

    p_i = exp(-Delta G_i / (R * T)) / sum_j exp(-Delta G_j / (R * T))

    Parameters
    ----------
    energies_hartree : Sequence[float]
        Conformer Gibbs free energies in Hartrees.
    temperature : float
        Temperature in Kelvin.

    Returns
    -------
    np.ndarray
        Normalized probabilities (summing to 1.0) for each conformer.
    """
    if temperature <= 0:
        raise ValueError("Temperature must be greater than zero")
    energies = np.asarray(energies_hartree, dtype=float)
    if energies.size == 0 or not np.all(np.isfinite(energies)):
        raise ValueError("Conformer energies must be a non-empty finite sequence")

    energies_kcal = energies * HARTREE_TO_KCAL
    rel_kcal = energies_kcal - np.min(energies_kcal)
    exp_terms = np.exp(-rel_kcal / (GAS_CONSTANT_R * temperature))
    partition_q = np.sum(exp_terms)
    return exp_terms / partition_q


def conformational_entropy(
    weights: np.ndarray,
) -> float:
    """
    Compute conformational mixing entropy in cal / (mol * K) from Boltzmann populations.

    S_conf = -R * sum_i (p_i * ln(p_i))

    Parameters
    ----------
    weights : np.ndarray
        Normalized Boltzmann population distribution.

    Returns
    -------
    float
        Conformational entropy in cal / (mol * K).
    """
    w = np.asarray(weights, dtype=float)
    w_nonzero = w[w > 1e-15]
    if len(w_nonzero) <= 1:
        return 0.0
    r_cal = GAS_CONSTANT_R * 1000.0  # kcal to cal
    return float(-r_cal * np.sum(w_nonzero * np.log(w_nonzero)))


def boltzmann_weighted_free_energy(
    energies_hartree: Sequence[float],
    temperature: float,
) -> float:
    """
    Return the ensemble Gibbs free energy in Hartrees from conformer Gibbs energies.

    G_ensemble = G_min - R*T*ln(sum_i exp(-(G_i - G_min) / (R*T)))

    Parameters
    ----------
    energies_hartree : Sequence[float]
        Conformer Gibbs free energies in Hartrees.
    temperature : float
        Temperature in Kelvin.

    Returns
    -------
    float
        Ensemble Gibbs free energy in Hartrees.
    """
    if temperature <= 0:
        raise ValueError("Temperature must be greater than zero")
    energies = np.asarray(energies_hartree, dtype=float)
    if energies.size == 0 or not np.all(np.isfinite(energies)):
        raise ValueError("Conformer energies must be a non-empty finite sequence")

    energies_kcal = energies * HARTREE_TO_KCAL
    minimum = np.min(energies_kcal)
    scaled = np.exp(-(energies_kcal - minimum) / (GAS_CONSTANT_R * temperature))
    ensemble_kcal = minimum - GAS_CONSTANT_R * temperature * np.log(np.sum(scaled))
    return float(ensemble_kcal / HARTREE_TO_KCAL)


def calculate_ensemble(
    label: str,
    conformer_paths: Sequence[str],
    energies_hartree: Sequence[float],
    temperature: float,
) -> EnsembleResult:
    """
    Construct a complete EnsembleResult with conformer breakdown, weights, and entropy.

    Parameters
    ----------
    label : str
        Catalyst dimer label (e.g., 'RR' or 'RS').
    conformer_paths : Sequence[str]
        File paths for the evaluated conformers.
    energies_hartree : Sequence[float]
        DFT Gibbs free energies in Hartrees.
    temperature : float
        Temperature in Kelvin.

    Returns
    -------
    EnsembleResult
        Rich ensemble thermodynamic summary.
    """
    energies = np.asarray(energies_hartree, dtype=float)
    weights = boltzmann_weights(energies, temperature)
    s_conf = conformational_entropy(weights)
    g_ensemble_hartree = boltzmann_weighted_free_energy(energies, temperature)
    g_ensemble_kcal = g_ensemble_hartree * HARTREE_TO_KCAL

    min_kcal = np.min(energies * HARTREE_TO_KCAL)
    conformers: List[Conformer] = []
    for idx, (path, e_hartree, w) in enumerate(zip(conformer_paths, energies, weights), start=1):
        e_kcal = e_hartree * HARTREE_TO_KCAL
        rel_kcal = e_kcal - min_kcal
        conformers.append(
            Conformer(
                index=idx,
                path=str(path),
                energy_hartree=float(e_hartree),
                energy_kcal=float(e_kcal),
                relative_energy_kcal=float(rel_kcal),
                boltzmann_weight=float(w),
            )
        )

    return EnsembleResult(
        label=label,
        conformers=conformers,
        ensemble_gibbs_hartree=float(g_ensemble_hartree),
        ensemble_gibbs_kcal=float(g_ensemble_kcal),
        lowest_gibbs_hartree=float(np.min(energies)),
        conformer_count=len(conformers),
        conformational_entropy_cal_mol_k=float(s_conf),
    )

