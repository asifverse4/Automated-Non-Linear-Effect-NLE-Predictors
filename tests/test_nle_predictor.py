import numpy as np
import pytest

from nle_predictor import (
    DFTEvaluator,
    NLECalculator,
    boltzmann_weighted_free_energy,
    validate_orca_output,
    validate_xyz_file,
)


def test_validate_xyz_file_accepts_valid_structure(tmp_path):
    xyz_path = tmp_path / "dimer.xyz"
    xyz_path.write_text("2\nwater-like test\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74\n")

    assert validate_xyz_file(str(xyz_path)) == 2


def test_validate_xyz_file_rejects_bad_coordinates(tmp_path):
    xyz_path = tmp_path / "bad.xyz"
    xyz_path.write_text("2\ninvalid\nH 0.0 0.0\nH 0.0 0.0 0.74\n")

    with pytest.raises(ValueError, match="Invalid XYZ coordinate"):
        validate_xyz_file(str(xyz_path))


def test_extract_free_energy_uses_final_energy(tmp_path):
    output_path = tmp_path / "orca.out"
    output_path.write_text(
        "ORCA TERMINATED NORMALLY\n"
        "Number of imaginary frequencies        0\n"
        "SCF energy -100.0\n"
        "Final Gibbs free energy             -99.87654321 Eh\n"
    )

    energy = DFTEvaluator().extract_free_energy(str(output_path))

    assert energy == pytest.approx(-99.87654321)


def test_orca_validation_rejects_incomplete_output(tmp_path):
    output_path = tmp_path / "incomplete.out"
    output_path.write_text("Final Gibbs free energy -99.0 Eh\n")

    with pytest.raises(RuntimeError, match="terminate normally"):
        validate_orca_output(str(output_path))


def test_orca_validation_rejects_imaginary_frequency(tmp_path):
    output_path = tmp_path / "imaginary.out"
    output_path.write_text(
        "ORCA TERMINATED NORMALLY\n"
        "Number of imaginary frequencies        1\n"
    )

    with pytest.raises(RuntimeError, match="imaginary frequencies"):
        validate_orca_output(str(output_path))


def test_generate_input_preserves_electronic_state(tmp_path, monkeypatch):
    xyz_path = tmp_path / "charged_dimer.xyz"
    xyz_path.write_text("2\ncharged test\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74\n")
    monkeypatch.chdir(tmp_path)

    input_path = DFTEvaluator().generate_input(
        str(xyz_path), "charged", charge=-1, multiplicity=2
    )

    content = (tmp_path / input_path).read_text()
    assert "Charge -1" in content
    assert "Mult 2" in content


def test_boltzmann_ensemble_is_lower_than_its_lowest_member():
    ensemble = boltzmann_weighted_free_energy([-100.0, -99.99], 298.15)

    assert ensemble < -100.0


def test_boltzmann_ensemble_rejects_invalid_temperature():
    with pytest.raises(ValueError, match="Temperature"):
        boltzmann_weighted_free_energy([-100.0], 0)


def test_mass_balance_preserves_catalyst_material():
    calculator = NLECalculator()
    monomer_r, monomer_s = calculator.solve_mass_balance(0.25, 1e5, 1e2)
    residual_r = monomer_r + 2e5 * monomer_r**2 + 1e2 * monomer_r * monomer_s - 0.625
    residual_s = monomer_s + 2e5 * monomer_s**2 + 1e2 * monomer_r * monomer_s - 0.375

    assert monomer_r >= 0
    assert monomer_s >= 0
    assert max(abs(residual_r), abs(residual_s)) < 1e-8


def test_mass_balance_rejects_invalid_catalyst_ee():
    with pytest.raises(ValueError, match="Catalyst ee"):
        NLECalculator().solve_mass_balance(1.1, 1e5, 1e2)


def test_kagan_model_returns_expected_endpoints():
    calculator = NLECalculator()
    catalyst_ee = np.array([0.0, 1.0])

    product_ee = calculator.simulate_kagan_model(catalyst_ee, 0.0)

    assert product_ee[0] == pytest.approx(0.0, abs=1e-10)
    assert product_ee[1] == pytest.approx(1.0, abs=1e-10)