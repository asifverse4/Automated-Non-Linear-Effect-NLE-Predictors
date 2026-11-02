"""Tests for the Noyori DAIB Aggregation Model in nle_predictor.models.noyori_daib."""

import numpy as np
import pytest

from nle_predictor.models.noyori_daib import NoyoriDAIBModel


def test_noyori_daib_amplification():
    # In Noyori's DAIB system, Kd(hetero) is ~1000x smaller than Kd(homo)
    # causing dramatic asymmetric amplification
    model = NoyoriDAIBModel(ee_0=0.98, kd_homo=1e-3, kd_hetero=1e-6)
    assert model.k_homo_assoc == pytest.approx(1e3)
    assert model.k_hetero_assoc == pytest.approx(1e6)
    assert model.delta_g_kcal < -3.0  # Heterochiral dimer is strongly favored

    cat_ee = np.array([0.0, 0.10, 0.20, 0.50, 1.0])
    res = model.simulate(cat_ee)

    assert res.ee_product[0] == pytest.approx(0.0, abs=1e-8)
    assert res.ee_product[-1] == pytest.approx(0.98, abs=1e-4)

    # Even at 10% catalyst ee, DAIB gives >80% product ee
    assert res.ee_product[1] > 0.80
    assert "Noyori asymmetric amplification" in res.classification


def test_noyori_daib_endpoints_and_overrides():
    model = NoyoriDAIBModel(ee_0=1.0)
    # Override delta_g directly
    res = model.simulate([0.0, 0.5, 1.0], delta_g_kcal=-2.0, k_homo=1e4)

    assert res.ee_product[0] == pytest.approx(0.0, abs=1e-8)
    assert res.ee_product[2] == pytest.approx(1.0, abs=1e-8)
    assert res.product_at_50 > 0.85


def test_noyori_parameter_validation():
    with pytest.raises(ValueError, match="Temperature"):
        NoyoriDAIBModel(temp=-10.0)

    with pytest.raises(ValueError, match="Intrinsic enantioselectivity"):
        NoyoriDAIBModel(ee_0=1.5)

    with pytest.raises(ValueError, match="Dissociation constants"):
        NoyoriDAIBModel(kd_homo=-1.0)

