"""Tests for publication visualization in nle_predictor.viz."""

import os

import numpy as np

from nle_predictor.models.kagan_ml2 import NLECalculator
from nle_predictor.viz.dashboard import generate_interactive_html
from nle_predictor.viz.plotting import (
    plot_amplification_contour,
    plot_monomer_speciation,
    plot_nle_with_experimental,
)


def test_plot_nle_standalone(tmp_path):
    calc = NLECalculator()
    res = calc.simulate(np.linspace(0.0, 1.0, 21), delta_g_kcal=-1.5)

    out_file = tmp_path / "nle_standalone.png"
    plot_nle_with_experimental(res, save_path=str(out_file))

    assert os.path.isfile(out_file)
    assert os.path.getsize(out_file) > 1000


def test_plot_nle_with_experimental_and_residuals(tmp_path):
    calc = NLECalculator()
    res = calc.simulate(np.linspace(0.0, 1.0, 21), delta_g_kcal=-1.5)

    exp_cat = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
    exp_prod = np.array([0.3, 0.65, 0.85, 0.95, 0.99])

    out_file = tmp_path / "nle_with_exp.png"
    plot_nle_with_experimental(res, exp_cat=exp_cat, exp_prod=exp_prod, save_path=str(out_file))

    assert os.path.isfile(out_file)
    assert os.path.getsize(out_file) > 1000


def test_plot_monomer_speciation(tmp_path):
    calc = NLECalculator()
    res = calc.simulate(np.linspace(0.0, 1.0, 21), delta_g_kcal=-1.5)

    out_file = tmp_path / "speciation.png"
    plot_monomer_speciation(res, save_path=str(out_file))

    assert os.path.isfile(out_file)
    assert os.path.getsize(out_file) > 1000


def test_plot_amplification_contour(tmp_path):
    calc = NLECalculator()
    out_file = tmp_path / "contour.png"
    plot_amplification_contour(calc, delta_g_range=(-2.0, 2.0), n_points=15, save_path=str(out_file))

    assert os.path.isfile(out_file)
    assert os.path.getsize(out_file) > 1000


def test_generate_interactive_html_dashboard(tmp_path):
    out_html = tmp_path / "dashboard.html"
    result_path = generate_interactive_html(str(out_html))

    assert os.path.isfile(result_path)
    assert os.path.getsize(result_path) > 5000
    with open(result_path, "r", encoding="utf-8") as f:
        html = f.read()
    assert "NLE-Predictor: Interactive Chemistry Dashboard" in html
    assert "solveMassBalance" in html
    assert "Kagan 1986" in html
    assert "Noyori DAIB" in html
