"""
Unit tests for non-parametric bootstrap Monte Carlo uncertainty quantification.
"""

import numpy as np
import pytest

from nle_predictor.fitting.bootstrap import run_bootstrap_monte_carlo
from nle_predictor.models.kagan_ml2 import NLECalculator


def test_residual_bootstrap_confidence_intervals():
    """Residual resampling bootstrap accurately encloses true Delta G."""
    calc = NLECalculator(temp=298.15)
    true_dg = -1.5  # kcal/mol
    ee_cats = np.linspace(0.05, 0.95, 10)
    ee_prods_clean = calc.simulate_kagan_model(ee_cats, delta_g_kcal=true_dg, k_homo_abs=1e4)

    # Add small deterministic noise
    rng = np.random.default_rng(123)
    noise = rng.normal(0.0, 0.01, size=len(ee_cats))
    ee_prods_noisy = np.clip(ee_prods_clean + noise, 0.0, 1.0)

    res = run_bootstrap_monte_carlo(
        ee_cats,
        ee_prods_noisy,
        n_bootstraps=100,
        method="residual",
        temp=298.15,
        k_homo=1e4,
        seed=42,
    )

    assert res.n_bootstraps == 100
    assert res.converged_fraction > 0.95
    assert len(res.delta_g_samples) >= 95

    # 95% CI must enclose the true parameter
    ci_low, ci_high = res.delta_g_ci_95
    assert ci_low < res.delta_g_median < ci_high
    assert ci_low <= true_dg <= ci_high

    # Custom percentile interval check
    p_low, p_high = res.percentile_interval(alpha=0.10)  # 90% CI
    assert ci_low <= p_low < p_high <= ci_high


def test_case_bootstrap_execution():
    """Case resampling bootstrap executes and computes empirical statistics."""
    calc = NLECalculator(temp=298.15)
    ee_cats = np.linspace(0.1, 0.9, 8)
    ee_prods = calc.simulate_kagan_model(ee_cats, delta_g_kcal=-0.8, k_homo_abs=1e4)

    res = run_bootstrap_monte_carlo(
        ee_cats,
        ee_prods,
        n_bootstraps=60,
        method="case",
        temp=298.15,
        k_homo=1e4,
        seed=101,
    )
    assert res.converged_fraction > 0.85
    assert abs(res.delta_g_median - (-0.8)) < 0.25


def test_bootstrap_with_ee_0_fitting():
    """Bootstrap with simultaneous ee_0 estimation yields both parameter distributions."""
    calc = NLECalculator(temp=298.15)
    ee_cats = np.linspace(0.1, 0.9, 9)
    ee_prods = calc.simulate_kagan_model(ee_cats, delta_g_kcal=-1.0, k_homo_abs=1e4, ee_0=0.92)

    res = run_bootstrap_monte_carlo(
        ee_cats,
        ee_prods,
        n_bootstraps=50,
        method="residual",
        fit_ee_0=True,
        seed=999,
    )

    assert res.fit_ee_0 is True
    assert res.ee_0_samples is not None
    assert len(res.ee_0_samples) == len(res.delta_g_samples)
    assert res.ee_0_ci_95 is not None
    assert 0.85 <= res.ee_0_median <= 0.98


def test_bootstrap_invalid_inputs_raise():
    """Verify input validation guards against invalid parameters."""
    with pytest.raises(ValueError, match="must have identical lengths"):
        run_bootstrap_monte_carlo(np.array([0.1, 0.2]), np.array([0.1]))

    with pytest.raises(ValueError, match="Need at least 3 data points"):
        run_bootstrap_monte_carlo(np.array([0.1, 0.2]), np.array([0.1, 0.2]))

    with pytest.raises(ValueError, match="Unknown bootstrap method"):
        run_bootstrap_monte_carlo(
            np.array([0.1, 0.2, 0.3]),
            np.array([0.1, 0.2, 0.3]),
            method="invalid_method",
        )
