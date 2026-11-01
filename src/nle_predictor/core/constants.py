"""
Physical and thermodynamic constants for quantum chemistry and kinetic modeling.
All values adhere to CODATA recommended values.
"""

import math

# Fundamental Thermodynamic and Conversion Constants
HARTREE_TO_KCAL = 627.5094740631  # kcal / mol per Hartree (CODATA 2018)
HARTREE_TO_KJ = 2625.4996394799   # kJ / mol per Hartree
KCAL_TO_KJ = 4.184                # Standard thermochemical conversion
GAS_CONSTANT_R = 0.00198720425864083  # kcal / (mol * K)
GAS_CONSTANT_R_SI = 8.31446261815324  # J / (mol * K)
DEFAULT_TEMP = 298.15                 # Standard temperature in Kelvin (25 °C)
STANDARD_PRESSURE_BAR = 1.0           # Standard pressure in bar
IDEAL_GAS_CONSTANT_L_ATM = 0.082057366081  # L * atm / (mol * K)

# Standard Atomic Weights (IUPAC CIAAW) for element mass calculation
ATOMIC_WEIGHTS = {
    "H": 1.008, "HE": 4.0026, "LI": 6.94, "BE": 9.0122, "B": 10.81,
    "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998, "NE": 20.180,
    "NA": 22.990, "MG": 24.305, "AL": 26.982, "SI": 28.085, "P": 30.974,
    "S": 32.06, "CL": 35.45, "AR": 39.948, "K": 39.098, "CA": 40.078,
    "SC": 44.956, "TI": 47.867, "V": 50.942, "CR": 51.996, "MN": 54.938,
    "FE": 55.845, "CO": 58.933, "NI": 58.693, "CU": 63.546, "ZN": 65.38,
    "GA": 69.723, "GE": 72.630, "AS": 74.922, "SE": 78.971, "BR": 79.904,
    "KR": 83.798, "RB": 85.468, "SR": 87.62, "Y": 88.906, "ZR": 91.224,
    "NB": 92.906, "MO": 95.95, "RU": 101.07, "RH": 102.91, "PD": 106.42,
    "AG": 107.87, "CD": 112.41, "IN": 114.82, "SN": 118.71, "SB": 121.76,
    "TE": 127.60, "I": 126.90, "XE": 131.29, "CS": 132.91, "BA": 137.33,
    "LA": 138.91, "CE": 140.12, "PR": 140.91, "ND": 144.24, "SM": 150.36,
    "EU": 151.96, "GD": 157.25, "TB": 158.93, "DY": 162.50, "HO": 164.93,
    "ER": 167.26, "TM": 168.93, "YB": 173.05, "LU": 174.97, "HF": 178.49,
    "TA": 180.95, "W": 183.84, "RE": 186.21, "OS": 190.23, "IR": 192.22,
    "PT": 195.08, "AU": 196.97, "HG": 200.59, "TL": 204.38, "PB": 207.2,
    "BI": 208.98, "TH": 232.04, "PA": 231.04, "U": 238.03,
}


def standard_state_correction_gas_to_solution(
    temp_k: float = DEFAULT_TEMP,
    delta_n: int = -1,
) -> float:
    """
    Calculate the Gibbs free energy adjustment (in kcal/mol) for converting from
    the ideal gas standard state (1 atm) to the solution standard state (1 M = 1 mol/L).

    For dimerization: 2 Monomer -> 1 Dimer (delta_n = 1 - 2 = -1),
    entropy is lost upon standard-state concentration shift:
    Delta G_correction = - delta_n * R * T * ln(V_m)
    where V_m = R_atm * T (approx 24.46 L/mol at 298.15 K).

    Parameters
    ----------
    temp_k : float
        Temperature in Kelvin.
    delta_n : int
        Change in moles of solute (product moles - reactant moles).
        Default is -1 for 2 M -> 1 D dimerizations.

    Returns
    -------
    float
        Free energy correction in kcal/mol to be added to Delta G_gas.
    """
    if temp_k <= 0:
        raise ValueError("Temperature must be greater than zero Kelvin.")
    # Molar volume of ideal gas at 1 atm and T
    molar_volume_liters = IDEAL_GAS_CONSTANT_L_ATM * temp_k
    # Delta G = - delta_n * R * T * ln(V_m)
    return -float(delta_n) * GAS_CONSTANT_R * temp_k * math.log(molar_volume_liters)
