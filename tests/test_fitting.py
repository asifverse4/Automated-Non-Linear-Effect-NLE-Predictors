"""Tests for experimental fitting and parameter estimation in nle_predictor.fitting."""

import numpy as np
import pytest

from nle_predictor.fitting.experimental_fit import (
    fit_kagan_reservoir,
    load_experimental_csv,
)
from nle_predictor.models.kagan_ml2 import NLECalculator


def test_fit_recovers_known_delta_g():
    true_dg = -1.5  # kcal/mol
    temp = 298.15
    k_homo = 1e4

    calculator = NLECalculator(temp=temp)
    ee_cat = np.array([0.05, 0.15, 0.30, 0.50, 0.70, 0.85, 0.95])
    ee_prod_true = calculator.simulate_kagan_model(ee_cat, delta_g_kcal=true_dg, k_homo_abs=k_homo)

    # Add small synthetic experimental measurement noise
    np.random.seed(42)
    ee_prod_noisy = np.clip(ee_prod_true + np.random.normal(0, 0.005, size=len(ee_cat)), 0.0, 1.0)

    fit_res = fit_kagan_reservoir(ee_cat, ee_prod_noisy, temp=temp, k_homo=k_homo, initial_delta_g=0.0)

    assert fit_res.converged
    assert fit_res.r_squared > 0.98
    assert fit_res.delta_g_kcal == pytest.approx(true_dg, abs=0.15)
    assert fit_res.rmse < 0.015
    assert np.isfinite(fit_res.delta_g_stderr)


def test_load_experimental_csv_both_scales(tmp_path):
    # Test percentage scale
    csv_pct = tmp_path / "exp_pct.csv"
    csv_pct.write_text("ee_cat,ee_prod\n10.0,25.5\n50.0,85.2\n90.0,98.1\n")

    cat_p, prod_p = load_experimental_csv(str(csv_pct))
    np.testing.assert_allclose(cat_p, [0.10, 0.50, 0.90])
    np.testing.assert_allclose(prod_p, [0.255, 0.852, 0.981])

    # Test decimal fraction scale
    csv_frac = tmp_path / "exp_frac.csv"
    csv_frac.write_text("ee_cat,ee_prod\n0.10,0.255\n0.50,0.852\n0.90,0.981\n")

    cat_f, prod_f = load_experimental_csv(str(csv_frac))
    np.testing.assert_allclose(cat_f, [0.10, 0.50, 0.90])
    np.testing.assert_allclose(prod_f, [0.255, 0.852, 0.981])


def test_fit_with_intrinsic_ee_0():
    calculator = NLECalculator()
    ee_cat = np.array([0.1, 0.3, 0.5, 0.7, 0.9, 1.0])
    # Dimer is +NLE with monomer reaching max 88%
    ee_prod = calculator.simulate_kagan_model(ee_cat, delta_g_kcal=-1.2, k_homo_abs=1e4, ee_0=0.88)

    fit_res = fit_kagan_reservoir(ee_cat, ee_prod, k_homo=1e4, fit_ee_0=True)
    assert fit_res.converged
    assert fit_res.metadata["fitted_ee_0"] == pytest.approx(0.88, abs=0.03)
