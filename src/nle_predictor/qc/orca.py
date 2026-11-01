"""
ORCA electronic structure input generation, execution management, and thermochemistry extraction.
"""

import logging
import os
import re
import shutil
import subprocess
from typing import Dict, Optional

from nle_predictor.core.xyz import validate_xyz_file

logger = logging.getLogger("NLE_Predictor")


class Config:
    """Configuration and dependency management for quantum chemistry binaries."""
    CREST_CMD = "crest"
    ORCA_CMD = "orca"

    @classmethod
    def check_dependencies(cls):
        """Checks if required computational chemistry software is in the system PATH."""
        missing = []
        if shutil.which(cls.CREST_CMD) is None:
            missing.append("CREST")
        if shutil.which(cls.ORCA_CMD) is None:
            missing.append("ORCA")
        return missing


def validate_orca_output(out_file: str) -> None:
    """Validate that an ORCA output completed and has no imaginary modes."""
    if not os.path.isfile(out_file):
        raise FileNotFoundError(f"ORCA output file not found: {out_file}")

    with open(out_file, "r", encoding="utf-8", errors="replace") as output:
        content = output.read()

    if "ORCA TERMINATED NORMALLY" not in content:
        raise RuntimeError(f"ORCA did not terminate normally: {out_file}")

    imaginary_match = re.search(
        r"(?:Number of imaginary frequencies|imaginary frequencies)\s*[:=]?\s*(\d+)",
        content,
        flags=re.IGNORECASE,
    )
    if imaginary_match and int(imaginary_match.group(1)) > 0:
        raise RuntimeError(f"ORCA output contains imaginary frequencies: {out_file}")


def collect_software_versions() -> Dict[str, str]:
    """Return version strings for the external chemistry executables."""
    versions = {}
    for name, command in (("crest", Config.CREST_CMD), ("orca", Config.ORCA_CMD)):
        if shutil.which(command) is None:
            versions[name] = "unavailable"
            continue
        try:
            result = subprocess.run(
                [command, "--version"],
                capture_output=True,
                text=True,
                timeout=15,
                check=False,
            )
            version_output = (result.stdout or result.stderr).strip().splitlines()
            versions[name] = version_output[0] if version_output else "unknown"
        except (OSError, subprocess.SubprocessError) as error:
            versions[name] = f"unavailable: {error.__class__.__name__}"
    return versions


class DFTEvaluator:
    """
    Handles DFT evaluations using ORCA to compute final Gibbs Free Energies.
    """
    def __init__(self, orca_path: str = Config.ORCA_CMD, cores: int = 4, mock_mode: bool = False):
        self.orca_path = orca_path
        self.cores = cores
        self.mock_mode = mock_mode

    def generate_input(
        self,
        xyz_file: str,
        run_name: str,
        functional: str = "B3LYP",
        basis: str = "def2-SVP",
        dispersion: str = "D4",
        solvent: Optional[str] = None,
        charge: int = 0,
        multiplicity: int = 1,
    ) -> str:
        """Generates an ORCA input file for geometry optimization and frequency calculations."""
        validate_xyz_file(xyz_file)
        inp_path = f"{run_name}.inp"

        coords = []
        with open(xyz_file, "r", encoding="utf-8") as f:
            lines = f.readlines()
            if len(lines) > 2:
                coords = lines[2:]

        header = f"! {functional} {basis} {dispersion} Opt Freq\n"
        if solvent:
            header += f"! CPCM({solvent})\n"

        if self.cores > 1:
            header += f"%pal nprocs {self.cores} end\n"

        with open(inp_path, "w", encoding="utf-8") as f:
            f.write(header)
            f.write(
                f"%coords\n  CTyp xyz\n  Charge {charge}\n  Mult {multiplicity}\n  coords\n"
            )
            for line in coords:
                if line.strip():
                    f.write(f"    {line.strip()}\n")
            f.write("  end\nend\n")

        return inp_path

    def extract_free_energy(self, out_file: str) -> float:
        """Parses ORCA output to extract the Final Gibbs Free Energy in Hartrees."""
        validate_orca_output(out_file)
        gibbs_energy = None
        energy_pattern = re.compile(
            r"Final Gibbs free energy\s+([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)"
        )
        with open(out_file, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                match = energy_pattern.search(line)
                if match:
                    gibbs_energy = float(match.group(1))

        if gibbs_energy is None:
            logger.error(f"Could not find Gibbs free energy in {out_file}")
            raise ValueError(f"Free energy extraction failed for {out_file}")
        return gibbs_energy

    def run_dft(self, inp_file: str, output_dir: str) -> str:
        """Executes the ORCA calculation."""
        out_file = inp_file.replace(".inp", ".out")
        out_path = os.path.join(output_dir, out_file)
        inp_path = os.path.join(output_dir, inp_file)

        logger.info(f"Running DFT evaluation: {inp_file}...")

        if self.mock_mode:
            logger.info(f"[MOCK] Simulating ORCA execution for {inp_file}")
            os.makedirs(output_dir, exist_ok=True)
            with open(out_path, "w", encoding="utf-8") as out_f:
                energy = -1000.00366 if "RS" in inp_file else -1000.00000
                out_f.write("ORCA TERMINATED NORMALLY\n")
                out_f.write("Number of imaginary frequencies        0\n")
                out_f.write(f"Final Gibbs free energy             {energy:.6f} Eh\n")
            return out_path

        try:
            with open(out_path, "w", encoding="utf-8") as out_f:
                subprocess.run(
                    [self.orca_path, inp_path],
                    cwd=output_dir,
                    stdout=out_f,
                    stderr=subprocess.STDOUT,
                    check=True,
                )
            logger.info(f"DFT evaluation complete for {inp_file}.")
            return out_path
        except subprocess.CalledProcessError:
            logger.error(f"ORCA calculation failed. Check {out_path}.")
            raise
