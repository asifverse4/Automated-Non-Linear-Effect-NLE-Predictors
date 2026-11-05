"""
Unit tests for multi-temperature van 't Hoff deconvolution engine.
"""

import tempfile

import numpy as np
import pytest

from nle_predictor.fitting.temperature_fit import (
    fit_van_t_hoff_global,
    load_multi_temperature_csv,
)
from nle_predictor.models.kagan_ml2 import NLECalculator


def test_van_t_hoff_recovers_synthetic_thermodynamic_parameters():
    """Verify global multi-temperature regression recovers known Delta Delta H° and Delta Delta S°."""
    true_dh = -2.5   # kcal/mol
    true_ds = -6.0   # cal/(mol*K)
    temps = [260.0, 290.0, 320.0]

    ee_cat = np.linspace(0.0, 1.0, 9)
    datasets = {}

    for t in temps:
        calc = NLECalculator(temp=t)
        # Delta Delta G(T) = Delta Delta H - T * (Delta Delta S / 1000)
        dg_t = true_dh - t * (true_ds / 1000.0)
        ee_prod = calc.simulate_kagan_model(ee_cat, delta_g_kcal=dg_t, k_homo_abs=1e4, ee_0=1.0)
        datasets[t] = (ee_cat, ee_prod)

    res = fit_van_t_hoff_global(datasets, k_homo=1e4, initial_dh=0.0, initial_ds=0.0)

    assert res.r_squared > 0.9999
    assert res.rmse < 1e-4
    np.testing.assert_allclose(res.delta_delta_h_kcal, true_dh, atol=0.02)
    np.testing.assert_allclose(res.delta_delta_s_cal_mol_k, true_ds, atol=0.05)

    # Expected Isoinversion temperature T_inv = (-2.5 * 1000) / -6.0 = 416.67 K
    expected_t_inv = (true_dh * 1000.0) / true_ds
    assert res.isoinversion_temp_k == pytest.approx(expected_t_inv, abs=2.0)
    assert res.isoinversion_stderr is not None
    assert res.isoinversion_stderr > 0.0


def test_confidence_ellipse_geometry():
    """Confidence ellipse must return 2D coordinates enclosing the fitted center."""
    true_dh = -1.8
    true_ds = -4.0
    temps = [273.15, 298.15, 323.15]
    ee_cat = np.linspace(0.1, 0.9, 7)
    datasets = {}

    for t in temps:
        calc = NLECalculator(temp=t)
        dg_t = true_dh - t * (true_ds / 1000.0)
        # Add slight noise
        ee_prod = calc.simulate_kagan_model(ee_cat, delta_g_kcal=dg_t, k_homo_abs=1e4)
        datasets[t] = (ee_cat, ee_prod)

    res = fit_van_t_hoff_global(datasets, k_homo=1e4)
    h_pts, s_pts = res.confidence_ellipse_2d(confidence=0.95, n_points=50)

    assert len(h_pts) == 50
    assert len(s_pts) == 50
    # Center should be inside bounds of the ellipse
    assert np.min(h_pts) <= res.delta_delta_h_kcal <= np.max(h_pts)
    assert np.min(s_pts) <= res.delta_delta_s_cal_mol_k <= np.max(s_pts)


def test_load_multi_temperature_csv_roundtrip():
    """Verify loading and parsing of multi-temperature CSV file with Celsius detection."""
    csv_content = (
        "temperature,ee_cat,ee_prod\n"
        "-20,0.1,0.25\n"
        "-20,0.5,0.75\n"
        "-20,0.8,0.92\n"
        "25,0.1,0.18\n"
        "25,0.5,0.65\n"
        "25,0.8,0.88\n"
    )
    with tempfile.NamedTemporaryFile("w", delete=False, suffix=".csv") as tmp:
        tmp.write(csv_content)
        tmp_path = tmp.name

    data = load_multi_temperature_csv(tmp_path)
    # -20 C = 253.15 K, 25 C = 298.15 K
    assert len(data) == 2
    keys = list(data.keys())
    assert any(abs(k - 253.15) < 0.1 for k in keys)
    assert any(abs(k - 298.15) < 0.1 for k in keys)

    k_253 = [k for k in keys if abs(k - 253.15) < 0.1][0]
    t1_cat, t1_prod = data[k_253]
    assert len(t1_cat) == 3
    assert t1_cat[0] == pytest.approx(0.1)
    assert t1_prod[0] == pytest.approx(0.25)


def test_invalid_multi_temperature_inputs_raise():
    """Verify error checking for invalid multi-temperature configurations."""
    with pytest.raises(ValueError, match="Must provide at least 2 temperature datasets"):
        fit_van_t_hoff_global({298.15: (np.array([0.1, 0.5]), np.array([0.2, 0.6]))})
