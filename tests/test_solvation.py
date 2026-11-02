"""Tests for solvation thermochemistry and standard-state concentration corrections."""

import pytest

from nle_predictor.thermo.solvation import (
    SOLVENT_DATABASE,
    compute_standard_state_dimerization_constant,
    concentration_correction_1atm_to_1m,
    get_solvent_properties,
)


def test_solvent_database_lookup():
    tol = get_solvent_properties("Toluene")
    assert tol is not None
    assert tol.dielectric_constant == pytest.approx(2.38)
    assert tol.crest_name == "toluene"
    assert tol.orca_cpcm_name == "Toluene"

    dcm = get_solvent_properties("DCM")
    assert dcm is not None
    assert dcm.dielectric_constant == pytest.approx(8.93)
    assert dcm.orca_cpcm_name == "CH2Cl2"

    unknown = get_solvent_properties("nonexistent_solvent")
    assert unknown is None


def test_concentration_correction_1atm_to_1m():
    # At 298.15 K, V_m = 24.4654 L/mol
    # For dimerization (delta_n = -1):
    # Delta G_corr = -(-1) * R * T * ln(24.4654) ~ +1.894 kcal/mol
    corr = concentration_correction_1atm_to_1m(298.15, delta_n=-1)
    assert corr == pytest.approx(1.894, abs=0.01)

    # For dissociation (delta_n = +1):
    corr_diss = concentration_correction_1atm_to_1m(298.15, delta_n=1)
    assert corr_diss == pytest.approx(-1.894, abs=0.01)

    with pytest.raises(ValueError, match="Temperature"):
        concentration_correction_1atm_to_1m(-50)


def test_compute_standard_state_dimerization_constant():
    # If gas-phase dimerization Delta G is -5.0 kcal/mol,
    # with +1.894 kcal/mol solution shift, Delta G_soln is -3.106 kcal/mol
    # K_assoc = exp(-(-3.106) / (R*T)) ~ 188.8 M^-1
    k_soln = compute_standard_state_dimerization_constant(-5.0, temp_k=298.15)
    assert k_soln > 100.0
