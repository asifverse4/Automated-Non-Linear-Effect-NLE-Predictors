"""
Gaussian 09/16 output parser and thermochemistry extractor.
Allows computational chemists using Gaussian to interface with NLE Predictor.
"""

import logging
import os
import re

logger = logging.getLogger("NLE_Predictor")


def validate_gaussian_output(log_file: str) -> None:
    """Validate that a Gaussian calculation terminated normally and has no imaginary modes."""
    if not os.path.isfile(log_file):
        raise FileNotFoundError(f"Gaussian output file not found: {log_file}")

    with open(log_file, "r", encoding="utf-8", errors="replace") as f:
        content = f.read()

    if "Normal termination of Gaussian" not in content:
        raise RuntimeError(f"Gaussian did not terminate normally: {log_file}")

    # Check NImag count in summary lines
    nimag_match = re.search(r"NImag\s*=\s*(\d+)", content, flags=re.IGNORECASE)
    if nimag_match and int(nimag_match.group(1)) > 0:
        raise RuntimeError(f"Gaussian output contains {nimag_match.group(1)} imaginary mode(s): {log_file}")


def extract_gaussian_free_energy(log_file: str) -> float:
    """
    Extract the 'Sum of electronic and thermal Free Energies' from a Gaussian log file in Hartrees.
    """
    validate_gaussian_output(log_file)
    gibbs_pattern = re.compile(
        r"Sum of electronic and thermal Free Energies\s*=\s*([-+]?\d+\.\d+)"
    )

    gibbs_energy = None
    with open(log_file, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            match = gibbs_pattern.search(line)
            if match:
                gibbs_energy = float(match.group(1))

    if gibbs_energy is None:
        raise ValueError(f"Could not find Gibbs free energy in Gaussian file: {log_file}")

    return gibbs_energy
