"""
Solvation thermochemistry, solvent dielectric mappings, and standard-state concentration corrections.
Converts gas-phase quantum chemistry free energies (1 atm) to solution-phase standard states (1 M).
"""

import math
from typing import Dict, NamedTuple, Optional

from nle_predictor.core.constants import DEFAULT_TEMP, GAS_CONSTANT_R, IDEAL_GAS_CONSTANT_L_ATM


class SolventProperties(NamedTuple):
    """Physical properties and dielectric constants for common reaction solvents."""
    name: str
    dielectric_constant: float  # Relative permittivity (epsilon)
    refractive_index: float
    crest_name: str
    orca_cpcm_name: str


# Comprehensive database of solvents used in asymmetric catalysis literature
SOLVENT_DATABASE: Dict[str, SolventProperties] = {
    "toluene": SolventProperties("toluene", 2.38, 1.496, "toluene", "Toluene"),
    "benzene": SolventProperties("benzene", 2.28, 1.501, "benzene", "Benzene"),
    "dichloromethane": SolventProperties("dichloromethane", 8.93, 1.424, "ch2cl2", "CH2Cl2"),
    "dcm": SolventProperties("dcm", 8.93, 1.424, "ch2cl2", "CH2Cl2"),
    "chloroform": SolventProperties("chloroform", 4.81, 1.446, "chcl3", "CHCl3"),
    "thf": SolventProperties("thf", 7.58, 1.407, "thf", "THF"),
    "tetrahydrofuran": SolventProperties("tetrahydrofuran", 7.58, 1.407, "thf", "THF"),
    "ether": SolventProperties("ether", 4.33, 1.353, "ether", "Et2O"),
    "diethylether": SolventProperties("diethylether", 4.33, 1.353, "ether", "Et2O"),
    "acetonitrile": SolventProperties("acetonitrile", 36.64, 1.344, "ch3cn", "Acetonitrile"),
    "meoh": SolventProperties("meoh", 32.70, 1.328, "methanol", "Methanol"),
    "methanol": SolventProperties("methanol", 32.70, 1.328, "methanol", "Methanol"),
    "hexane": SolventProperties("hexane", 1.88, 1.375, "hexane", "Hexane"),
    "water": SolventProperties("water", 80.1, 1.333, "water", "Water"),
}


def get_solvent_properties(solvent_name: str) -> Optional[SolventProperties]:
    """Retrieve dielectric and software alias properties for a solvent."""
    return SOLVENT_DATABASE.get(solvent_name.strip().lower())


def concentration_correction_1atm_to_1m(
    temp_k: float = DEFAULT_TEMP,
    delta_n: int = -1,
) -> float:
    """
    Calculate the free energy shift (in kcal/mol) for transitioning from the ideal
    gas standard state (1 atm) to the solution standard state (1 M = 1 mol/L).

    Delta G_corr = - delta_n * R * T * ln(V_m)
    where V_m = (R_atm * T / 1.0 atm) is the molar volume of ideal gas in liters.

    For dimerization (2 M -> 1 D, delta_n = -1):
    Delta G_corr = + R * T * ln(V_m) ~ +1.894 kcal/mol at 298.15 K.

    Parameters
    ----------
    temp_k : float
        Temperature in Kelvin.
    delta_n : int
        Change in moles of particles (moles product - moles reactant).
    """
    if temp_k <= 0:
        raise ValueError("Temperature must be greater than zero Kelvin.")
    v_molar = IDEAL_GAS_CONSTANT_L_ATM * temp_k
    return -float(delta_n) * GAS_CONSTANT_R * temp_k * math.log(v_molar)


def compute_standard_state_dimerization_constant(
    delta_g_dimer_gas_kcal: float,
    temp_k: float = DEFAULT_TEMP,
) -> float:
    """
    Compute the standard-state solution dimerization constant K_assoc (in M^-1)
    from computed gas-phase Delta G, applying the +1.89 kcal/mol standard state shift.
    """
    corr = concentration_correction_1atm_to_1m(temp_k, delta_n=-1)
    delta_g_soln = delta_g_dimer_gas_kcal + corr
    return float(math.exp(-delta_g_soln / (GAS_CONSTANT_R * temp_k)))

