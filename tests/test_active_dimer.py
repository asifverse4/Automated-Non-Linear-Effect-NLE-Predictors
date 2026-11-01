"""Tests for the Kagan Active Dimer Model in nle_predictor.models.kagan_active."""

import numpy as np
import pytest

from nle_predictor.models.kagan_active import KaganActiveDimerModel
from nle_predictor.models.kagan_ml2 import NLECalculator


def test_active_dimer_reduces_to_reservoir_when_g_is_zero():
    cat_range = np.linspace(0.0, 1.0, 11)
    delta_g = -1.5

    reservoir_calc = NLECalculator()
    ee_reservoir = reservoir_calc.simulate_kagan_model(cat_range, delta_g, k_homo_abs=1e4)

    active_model = KaganActiveDimerModel(g=0.0, g_prime=0.0, ee_0=1.0, ee_d=1.0)
    result_active = active_model.simulate(cat_range, delta_g, k_homo=1e4)

    np.testing.assert_allclose(result_active.ee_product, ee_reservoir, atol=1e-7)


def test_active_dimer_endpoints():
    active_model = KaganActiveDimerModel(g=1.5, g_prime=0.2, ee_0=1.0, ee_d=1.0)
    res = active_model.simulate([0.0, 1.0], delta_g=0.0, k_homo=1e4)

    assert res.ee_product[0] == pytest.approx(0.0, abs=1e-8)
    assert res.ee_product[1] == pytest.approx(1.0, abs=1e-8)


def test_intrinsic_monomer_selectivity():
    # If monomer ee_0 is 80%, product ee at enantiopure catalyst should reflect monomer/dimer mixture
    active_model = KaganActiveDimerModel(g=0.0, ee_0=0.85)
    res = active_model.simulate([1.0], delta_g=-1.0, k_homo=1e4)
    assert res.ee_product[0] == pytest.approx(0.85, abs=1e-5)


def test_active_dimer_parameter_validation():
    with pytest.raises(ValueError, match="Relative rate"):
        KaganActiveDimerModel(g=-0.1)

    with pytest.raises(ValueError, match="Intrinsic ee"):
        KaganActiveDimerModel(ee_0=1.5)

