"""
Unit tests for higher-order Kagan aggregation models (ML3 trimer and ML4 tetramer).
"""

import numpy as np
import pytest

from nle_predictor.models.kagan_oligomers import KaganOligomerModel


def test_combinatorial_weights():
    """Verify statistical symmetry coefficients binom(n, k)."""
    assert KaganOligomerModel.combinatorial_weights(2) == [1, 2, 1]
    assert KaganOligomerModel.combinatorial_weights(3) == [1, 3, 3, 1]
    assert KaganOligomerModel.combinatorial_weights(4) == [1, 4, 6, 4, 1]
    assert KaganOligomerModel.combinatorial_weights(5) == [1, 5, 10, 10, 5, 1]


def test_trimer_statistical_limit_is_linear():
    """When g=1 (Delta G = 0), Kagan reservoir trimer model must yield linear NLE."""
    model = KaganOligomerModel(order=3)
    ee_cats = np.linspace(0.0, 1.0, 11)
    res = model.simulate(ee_cats, delta_g_kcal=0.0, k_homo=1e4, ee_0=1.0)

    assert res.classification == "near-linear"
    np.testing.assert_allclose(res.ee_product, ee_cats, atol=1e-3)


def test_tetramer_statistical_limit_is_linear():
    """When g=1 (Delta G = 0), Kagan reservoir tetramer model must yield linear NLE."""
    model = KaganOligomerModel(order=4)
    ee_cats = np.linspace(0.0, 1.0, 11)
    res = model.simulate(ee_cats, delta_g_kcal=0.0, k_homo=1e4, ee_0=1.0)

    assert res.classification == "near-linear"
    np.testing.assert_allclose(res.ee_product, ee_cats, atol=1e-3)


def test_trimer_positive_nle_amplification():
    """Heterochiral trimer stabilization (Delta G < 0) yields positive NLE."""
    model = KaganOligomerModel(order=3)
    ee_cats = np.linspace(0.0, 1.0, 11)
    # delta_g = -1.5 kcal/mol favors heterochiral trimer (R2S and RS2)
    res = model.simulate(ee_cats, delta_g_kcal=-1.5, k_homo=1e4, ee_0=1.0)

    assert res.classification == "positive NLE"
    # At 20% catalyst ee, product ee should be significantly amplified
    idx_20 = 2  # 0.2
    assert res.ee_product[idx_20] > ee_cats[idx_20]
    assert res.product_at_50 > 0.50


def test_trimer_negative_nle_depletion():
    """Homochiral trimer stabilization (Delta G > 0) yields negative NLE."""
    model = KaganOligomerModel(order=3)
    ee_cats = np.linspace(0.0, 1.0, 11)
    # delta_g = +1.5 kcal/mol favors homochiral trimer (R3, S3)
    res = model.simulate(ee_cats, delta_g_kcal=1.5, k_homo=1e4, ee_0=1.0)

    assert res.classification == "negative NLE"
    idx_50 = 5  # 0.5
    assert res.ee_product[idx_50] < ee_cats[idx_50]


def test_mass_conservation():
    """Monomer and aggregate concentrations must strictly sum to total C_R and C_S."""
    for n in [3, 4]:
        model = KaganOligomerModel(order=n)
        k_aggs = model.compute_equilibrium_constants(k_homo=5e3, delta_g_kcal=-1.0)
        c_tot = 1.5

        for ee in [0.0, 0.25, 0.5, 0.8, 1.0]:
            mr, ms, aggs = model.solve_mass_balance(ee, k_homo=5e3, k_aggregates=k_aggs, c_tot=c_tot)
            k_powers = np.arange(n, -1, -1)
            s_powers = n - k_powers

            total_r = mr + np.sum(k_powers * aggs)
            total_s = ms + np.sum(s_powers * aggs)

            expected_r = c_tot * 0.5 * (1.0 + ee)
            expected_s = c_tot * 0.5 * (1.0 - ee)

            np.testing.assert_allclose(total_r, expected_r, rtol=1e-5, atol=1e-6)
            np.testing.assert_allclose(total_s, expected_s, rtol=1e-5, atol=1e-6)


def test_endpoints_and_intrinsic_ee():
    """Check ee_cat=0 gives ee_prod=0, and ee_cat=1 gives ee_prod=ee_0."""
    model = KaganOligomerModel(order=3)
    ee_cats = np.array([0.0, 1.0])
    res = model.simulate(ee_cats, delta_g_kcal=-2.0, k_homo=1e4, ee_0=0.88)

    assert abs(res.ee_product[0]) < 1e-6
    np.testing.assert_allclose(res.ee_product[1], 0.88, atol=1e-4)


def test_active_aggregates_mode():
    """Active aggregate mode allows catalytic turnover from higher-order clusters."""
    model = KaganOligomerModel(order=3)
    ee_cats = np.linspace(0.0, 1.0, 5)
    # Aggregates [R3, R2S, RS2, S3] have chiral selectivities [0.9, 0.3, -0.3, -0.9]
    res = model.simulate(
        ee_cats,
        delta_g_kcal=-0.5,
        k_homo=1e3,
        ee_0=0.95,
        active_aggregates=True,
        aggregate_ee=[0.9, 0.3, -0.3, -0.9],
    )
    assert len(res.ee_product) == len(ee_cats)
    assert res.ee_product[0] == pytest.approx(0.0, abs=1e-5)
    assert res.ee_product[-1] > 0.8


def test_invalid_parameters_raise():
    """Verify robust parameter validation."""
    with pytest.raises(ValueError, match="order must be >= 2"):
        KaganOligomerModel(order=1)

    with pytest.raises(ValueError, match="Temperature must be greater than zero"):
        KaganOligomerModel(order=3, temp=-10)

    model = KaganOligomerModel(order=3)
    with pytest.raises(ValueError, match="k_homo must be non-negative"):
        model.compute_equilibrium_constants(k_homo=-5)

    with pytest.raises(ValueError, match="Catalyst ee must be in"):
        model.solve_mass_balance(ee_cat=1.5, k_homo=1e3, k_aggregates=np.ones(4))
