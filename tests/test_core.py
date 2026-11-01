"""Tests for nle_predictor.core: constants, types, and xyz operations."""

import numpy as np
import pytest

from nle_predictor.core.constants import (
    DEFAULT_TEMP,
    GAS_CONSTANT_R,
    HARTREE_TO_KCAL,
    standard_state_correction_gas_to_solution,
)
from nle_predictor.core.types import (
    Conformer,
    EnsembleResult,
    ModelType,
    XYZStructure,
)
from nle_predictor.core.xyz import (
    calculate_rmsd,
    read_xyz,
    validate_xyz_file,
    write_xyz,
)


def test_standard_state_correction():
    # At 298.15 K, delta_n = -1 (dimerization: 2M -> D)
    # Delta G_corr = -(-1) * R * T * ln(24.4654) ~ +1.89 kcal/mol
    corr = standard_state_correction_gas_to_solution(DEFAULT_TEMP, delta_n=-1)
    assert 1.85 < corr < 1.95

    with pytest.raises(ValueError, match="Temperature"):
        standard_state_correction_gas_to_solution(0.0)


def test_xyz_structure_centroid():
    atoms = ["C", "O"]
    coords = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.2]])
    structure = XYZStructure(atoms=atoms, coordinates=coords)
    assert structure.atom_count == 2
    np.testing.assert_allclose(structure.centroid, [0.0, 0.0, 0.6])


def test_read_and_write_xyz_roundtrip(tmp_path):
    file_path = tmp_path / "water.xyz"
    atoms = ["O", "H", "H"]
    coords = np.array([
        [0.000, 0.000, 0.117],
        [0.000, 0.757, -0.469],
        [0.000, -0.757, -0.469],
    ])
    write_xyz(str(file_path), atoms, coords, comment="Water Molecule")

    struct = read_xyz(str(file_path))
    assert struct.atom_count == 3
    assert struct.atoms == ["O", "H", "H"]
    assert struct.formula == "H2O"
    assert pytest.approx(struct.molecular_weight, rel=1e-3) == 18.015
    np.testing.assert_allclose(struct.coordinates, coords, atol=1e-5)


def test_kabsch_rmsd_invariance():
    # Identical coordinates should have 0 RMSD
    coords = np.array([
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [-1.0, -1.0, -1.0],
    ])
    assert calculate_rmsd(coords, coords) == pytest.approx(0.0, abs=1e-8)

    # Pure translation should yield 0 RMSD
    translated = coords + np.array([5.2, -3.1, 7.8])
    assert calculate_rmsd(coords, translated) == pytest.approx(0.0, abs=1e-7)

    # 90 degree rotation around z-axis
    rotation_matrix = np.array([
        [0.0, -1.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
    ])
    rotated = np.dot(coords, rotation_matrix)
    assert calculate_rmsd(coords, rotated) == pytest.approx(0.0, abs=1e-7)


def test_ensemble_result_dict_export():
    conf = Conformer(
        index=1,
        path="conf_001.xyz",
        energy_hartree=-500.0,
        energy_kcal=-313754.7,
        relative_energy_kcal=0.0,
        boltzmann_weight=1.0,
    )
    ensemble = EnsembleResult(
        label="RR",
        conformers=[conf],
        ensemble_gibbs_hartree=-500.0,
        ensemble_gibbs_kcal=-313754.7,
        lowest_gibbs_hartree=-500.0,
        conformer_count=1,
    )
    d = ensemble.to_dict()
    assert d["label"] == "RR"
    assert d["conformer_count"] == 1
    assert len(d["conformers"]) == 1

