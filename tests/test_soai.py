"""
Unit tests for Soai asymmetric autocatalysis and hyper-amplification kinetic engine.
"""

import numpy as np
import pytest

from nle_predictor.models.soai import SoaiAutocatalysisModel


def test_racemic_seed_preserves_symmetry():
    """A strictly racemic seed (ee=0) must yield strictly racemic product (ee=0)."""
    model = SoaiAutocatalysisModel(
        k_uncat=1e-4,
        k_auto=5.0,
        auto_order=2.0,
        k_hetero=500.0,
        k_homo=10.0,
    )
    traj = model.simulate_trajectory(
        init_a=0.1,
        init_b=0.1,
        seed_catalyst=0.002,
        seed_ee=0.0,
        t_span=(0.0, 500.0),
        n_points=50,
    )
    assert abs(traj.final_ee) < 1e-6
    np.testing.assert_allclose(traj.product_r, traj.product_s, atol=1e-6)


def test_hyper_amplification_from_low_ee():
    """A low seed ee (2%) must hyper-amplify to high enantiomeric excess (>70%)."""
    model = SoaiAutocatalysisModel(
        k_uncat=1e-7,
        k_auto=5e4,
        auto_order=2.0,
        k_hetero=5e4,
        k_homo=10.0,
        ee_0=1.0,
    )
    traj = model.simulate_trajectory(
        init_a=0.1,
        init_b=0.1,
        seed_catalyst=0.002,
        seed_ee=0.02,  # 2% initial ee
        t_span=(0.0, 3600.0),
        n_points=100,
    )
    assert traj.final_ee > 0.70
    assert traj.amplification_factor > 30.0
    assert traj.conversion[-1] > 0.80


def test_mass_conservation_in_trajectory():
    """Substrate consumed must equal product formed throughout the trajectory."""
    model = SoaiAutocatalysisModel(k_uncat=1e-4, k_auto=2.0)
    init_a = 0.08
    seed_cat = 0.002
    traj = model.simulate_trajectory(
        init_a=init_a,
        init_b=0.08,
        seed_catalyst=seed_cat,
        seed_ee=0.05,
        t_span=(0.0, 1000.0),
        n_points=100,
    )
    # Total product formed at any time t is P_R(t) + P_S(t) - seed_cat
    # This must equal init_a - substrate_a(t)
    delta_substrate = init_a - traj.substrate_a
    delta_product = (traj.product_r + traj.product_s) - seed_cat
    np.testing.assert_allclose(delta_substrate, delta_product, rtol=1e-4, atol=1e-5)


def test_soai_nle_curve_simulation():
    """Simulate batch NLE curve and ensure positive NLE classification."""
    model = SoaiAutocatalysisModel(k_uncat=1e-7, k_auto=5e4, auto_order=2.0, k_hetero=5e4)
    ee_cats = np.array([0.0, 0.05, 0.1, 0.3, 0.5, 0.8, 1.0])
    res = model.simulate(ee_cats, seed_catalyst=0.002, t_final=3600.0)

    assert res.classification == "positive NLE"
    assert res.product_at_50 > 0.80
    assert res.max_amplification > 2.0


def test_soai_parameter_validation():
    """Ensure invalid kinetic parameters raise clear ValueError."""
    with pytest.raises(ValueError, match="Rate constants must be non-negative"):
        SoaiAutocatalysisModel(k_uncat=-1.0)

    with pytest.raises(ValueError, match="Autocatalytic order must be positive"):
        SoaiAutocatalysisModel(auto_order=0.0)

    with pytest.raises(ValueError, match="Intrinsic ee_0 must be in"):
        SoaiAutocatalysisModel(ee_0=1.5)

    model = SoaiAutocatalysisModel()
    with pytest.raises(ValueError, match="Substrate concentrations must be positive"):
        model.simulate_trajectory(init_a=-0.1)

    with pytest.raises(ValueError, match="Seed ee must be in"):
        model.simulate_trajectory(seed_ee=2.0)
