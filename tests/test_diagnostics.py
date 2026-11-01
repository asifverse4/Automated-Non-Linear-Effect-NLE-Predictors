"""Tests for NLE curve diagnostics in nle_predictor.models.diagnostics."""

import numpy as np
import pytest

from nle_predictor.models.diagnostics import (
    analyze_nle_diagnostics,
    calculate_area_under_asymmetry,
    calculate_initial_slope,
    calculate_max_amplification,
    find_inflection_points,
)
from nle_predictor.models.kagan_ml2 import NLECalculator


def test_diagnostics_linear_system():
    x = np.linspace(0.0, 1.0, 51)
    y = x.copy()

    s0 = calculate_initial_slope(x, y)
    assert s0 == pytest.approx(1.0, abs=1e-3)

    amp, cat_ee = calculate_max_amplification(x, y)
    assert amp == pytest.approx(1.0, abs=1e-3)

    aua = calculate_area_under_asymmetry(x, y)
    assert aua == pytest.approx(0.0, abs=1e-6)

    diag = analyze_nle_diagnostics(x, y)
    assert "linear" in diag.classification
    assert not diag.is_sigmoidal


def test_diagnostics_positive_nle():
    calc = NLECalculator()
    x = np.linspace(0.0, 1.0, 101)
    # Delta G = -2.0 kcal/mol -> K_hetero >> K_homo -> strong (+)NLE
    y = calc.simulate_kagan_model(x, delta_g_kcal=-2.0, k_homo_abs=1e4)

    diag = analyze_nle_diagnostics(x, y)
    assert diag.initial_slope > 1.2
    assert diag.max_amplification > 1.2
    assert diag.area_under_asymmetry > 0.05
    assert "positive NLE" in diag.classification
    assert diag.is_hyperbolic


def test_diagnostics_negative_nle():
    calc = NLECalculator()
    x = np.linspace(0.0, 1.0, 101)
    # Delta G = +2.0 kcal/mol -> K_homo >> K_hetero -> (-)NLE
    y = calc.simulate_kagan_model(x, delta_g_kcal=2.0, k_homo_abs=1e4)

    diag = analyze_nle_diagnostics(x, y)
    assert diag.initial_slope < 0.8
    assert diag.area_under_asymmetry < -0.05
    assert "negative NLE" in diag.classification


def test_inflection_point_detection_on_sigmoidal():
    # Sigmoid function: y = 1 / (1 + exp(-10*(x - 0.5)))
    x = np.linspace(0.0, 1.0, 101)
    y = 1.0 / (1.0 + np.exp(-12.0 * (x - 0.5)))
    # normalize between 0 and 1
    y = (y - y[0]) / (y[-1] - y[0])

    inflections = find_inflection_points(x, y)
    assert len(inflections) == 1
    assert inflections[0] == pytest.approx(0.5, abs=0.05)

