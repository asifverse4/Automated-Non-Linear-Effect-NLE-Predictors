"""Tests for Gaussian output validation and energy parsing in nle_predictor.qc.gaussian."""

import pytest

from nle_predictor.qc.gaussian import (
    extract_gaussian_free_energy,
    validate_gaussian_output,
)


def test_extract_gaussian_free_energy(tmp_path):
    log_file = tmp_path / "g16.log"
    log_file.write_text(
        " Harmonic frequencies (cm**-1), intensities (KM/Mole)\n"
        " NImag=0\n"
        " Zero-point correction=                           0.245100\n"
        " Thermal correction to Energy=                    0.261200\n"
        " Thermal correction to Enthalpy=                  0.262144\n"
        " Thermal correction to Gibbs Free Energy=         0.201500\n"
        " Sum of electronic and thermal Free Energies=    -875.123456\n"
        " Normal termination of Gaussian 16 at Sun Nov 01 12:00:00 2026.\n"
    )

    energy = extract_gaussian_free_energy(str(log_file))
    assert energy == pytest.approx(-875.123456)


def test_validate_gaussian_rejects_incomplete(tmp_path):
    log_file = tmp_path / "incomplete.log"
    log_file.write_text("Sum of electronic and thermal Free Energies= -875.123456\n")

    with pytest.raises(RuntimeError, match="terminate normally"):
        validate_gaussian_output(str(log_file))


def test_validate_gaussian_rejects_imaginary_modes(tmp_path):
    log_file = tmp_path / "ts.log"
    log_file.write_text(
        " NImag=1\n"
        " Sum of electronic and thermal Free Energies= -875.123456\n"
        " Normal termination of Gaussian 16\n"
    )

    with pytest.raises(RuntimeError, match="imaginary mode"):
        validate_gaussian_output(str(log_file))

