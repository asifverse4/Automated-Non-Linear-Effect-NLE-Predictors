"""
CREST (Conformer-Rotamer Ensemble Sampling Tool) execution and conformer extraction.
"""

import logging
import os
import shutil
import subprocess
from typing import List, Union

from nle_predictor.core.xyz import validate_xyz_file
from nle_predictor.qc.orca import Config

logger = logging.getLogger("NLE_Predictor")


class ConformerSampler:
    """
    Handles automated conformational sampling using CREST (xtb).
    Identifies the lowest energy conformers for a given chiral complex.
    """
    def __init__(self, crest_path: str = Config.CREST_CMD, cores: int = 4, mock_mode: bool = False):
        self.crest_path = crest_path
        self.cores = cores
        self.mock_mode = mock_mode

    def run_sampling(
        self,
        input_xyz: str,
        output_dir: str,
        solvent: str = None,
        max_conformers: int = 1,
    ) -> Union[str, List[str]]:
        """Runs CREST conformational search on a given XYZ file."""
        validate_xyz_file(input_xyz)
        if max_conformers <= 0:
            raise ValueError("max_conformers must be greater than zero")
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        logger.info(f"Starting conformer sampling for {input_xyz}...")

        base_name = os.path.splitext(os.path.basename(input_xyz))[0]
        run_dir = os.path.join(output_dir, f"crest_{base_name}")
        os.makedirs(run_dir, exist_ok=True)

        best_conformer_path = os.path.join(run_dir, "crest_best.xyz")

        if self.mock_mode:
            logger.info(f"[MOCK] Simulating CREST execution for {input_xyz}")
            conformer_paths = []
            for index in range(1, max_conformers + 1):
                conformer_path = os.path.join(run_dir, f"crest_conf_{index:03d}.xyz")
                with open(conformer_path, "w", encoding="utf-8") as conformer_file:
                    conformer_file.write(
                        f"2\nMock Conformer {base_name} {index}\n"
                        "H 0.0 0.0 0.0\nH 0.0 0.0 0.74\n"
                    )
                conformer_paths.append(conformer_path)
            shutil.copyfile(conformer_paths[0], best_conformer_path)
            return conformer_paths if max_conformers > 1 else best_conformer_path

        abs_input = os.path.abspath(input_xyz)
        cmd = [self.crest_path, abs_input, "-T", str(self.cores)]

        if solvent:
            cmd.extend(["-alpb", solvent])

        try:
            logger.info(f"Executing: {' '.join(cmd)} in {run_dir}")
            with open(os.path.join(run_dir, "crest.out"), "w", encoding="utf-8") as f_out:
                subprocess.run(cmd, cwd=run_dir, stdout=f_out, stderr=subprocess.STDOUT, check=True)
            logger.info(f"CREST sampling complete for {base_name}.")
            conformers_path = os.path.join(run_dir, "crest_conformers.xyz")
            if not os.path.exists(conformers_path):
                raise FileNotFoundError(
                    f"CREST completed but did not produce {conformers_path}"
                )

            with open(conformers_path, "r", encoding="utf-8") as conformers:
                lines = conformers.readlines()
            if len(lines) < 2:
                raise ValueError(f"Invalid CREST conformer file: {conformers_path}")
            conformer_paths = []
            cursor = 0
            while cursor < len(lines) and len(conformer_paths) < max_conformers:
                try:
                    atom_count = int(lines[cursor].strip())
                except ValueError as error:
                    raise ValueError(f"Invalid conformer atom count: {conformers_path}") from error
                block_size = atom_count + 2
                if atom_count <= 0 or cursor + block_size > len(lines):
                    raise ValueError(f"Incomplete CREST conformer file: {conformers_path}")
                conformer_path = os.path.join(
                    run_dir, f"crest_conf_{len(conformer_paths) + 1:03d}.xyz"
                )
                with open(conformer_path, "w", encoding="utf-8") as conformer_file:
                    conformer_file.writelines(lines[cursor:cursor + block_size])
                conformer_paths.append(conformer_path)
                cursor += block_size

            if not conformer_paths:
                raise ValueError(f"No conformers found in {conformers_path}")
            shutil.copyfile(conformer_paths[0], best_conformer_path)
            return conformer_paths if max_conformers > 1 else best_conformer_path
        except subprocess.CalledProcessError:
            logger.error(f"CREST failed for {base_name}. Check {run_dir}/crest.out")
            raise

