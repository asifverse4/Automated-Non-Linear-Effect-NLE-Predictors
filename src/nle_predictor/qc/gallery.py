"""
Chiral Ligand 3D Conformer Gallery and Automated Quantum Chemistry Input Generator.

Curated repository of benchmark homochiral and heterochiral catalyst aggregate structures:
1. DAIB (3-exo-(dimethylamino)isoborneol di-zinc dimer) - Noyori asymmetric alkylation
2. BINOL (1,1'-bi-2-naphthol titanium dimer) - Walsh/Kagan asymmetric Lewis acid catalysis
3. PyBOX (bis(oxazoline)pyridine copper dimer) - Evans/Pfaltz asymmetric cyclopropanation

Provides cross-platform input generation for:
- ORCA (5.0 / 6.0)
- Gaussian (09 / 16)
- Q-Chem (5.0 / 6.0)
with modern DFT functionals (B3LYP, wB97X-D4, M06-2X), def2 basis sets, and implicit solvation.
"""

import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from nle_predictor.core.types import XYZStructure
from nle_predictor.core.xyz import read_xyz, validate_xyz_file


@dataclass
class LigandMetadata:
    """Metadata for a curated benchmark chiral catalyst."""
    name: str
    full_name: str
    reaction_class: str
    default_charge: int
    default_multiplicity: int
    rr_file: str
    rs_file: str


class ChiralLigandGallery:
    """
    Gallery of benchmark chiral catalyst structures and multi-code DFT input generator.
    """

    _REGISTRY: Dict[str, LigandMetadata] = {
        "daib": LigandMetadata(
            name="daib",
            full_name="3-exo-(Dimethylamino)isoborneol di-zinc dimer",
            reaction_class="Asymmetric organozinc addition to aldehydes",
            default_charge=0,
            default_multiplicity=1,
            rr_file="dimer_RR.xyz",
            rs_file="dimer_RS.xyz",
        ),
        "binol": LigandMetadata(
            name="binol",
            full_name="1,1'-Bi-2-naphthol titanium di-isopropoxide dimer",
            reaction_class="Asymmetric titanium Lewis acid catalysis",
            default_charge=0,
            default_multiplicity=1,
            rr_file="binol_RR.xyz",
            rs_file="binol_RS.xyz",
        ),
        "pybox": LigandMetadata(
            name="pybox",
            full_name="Bis(oxazolinyl)pyridine copper(II) chloride dimer",
            reaction_class="Asymmetric Lewis acid and radical catalysis",
            default_charge=0,
            default_multiplicity=1,
            rr_file="pybox_RR.xyz",
            rs_file="pybox_RS.xyz",
        ),
    }

    def __init__(self, structures_dir: Optional[str] = None):
        if structures_dir is None:
            # Default to bundled benchmark directory relative to module
            # gallery.py -> qc -> nle_predictor -> src -> repo_root
            repo_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
            self.structures_dir = os.path.join(repo_root, "examples", "benchmarks", "structures")
        else:
            self.structures_dir = structures_dir

    def list_ligands(self) -> List[str]:
        """Return list of available benchmark ligand identifiers."""
        return sorted(list(self._REGISTRY.keys()))

    def get_metadata(self, ligand: str) -> LigandMetadata:
        """Retrieve structural metadata for a ligand."""
        lig_key = ligand.lower().strip()
        if lig_key not in self._REGISTRY:
            raise KeyError(f"Ligand '{ligand}' not recognized. Available: {self.list_ligands()}")
        return self._REGISTRY[lig_key]

    def get_structure_path(self, ligand: str, chirality: str) -> str:
        """Get absolute path to the XYZ file for a ligand and chirality ('RR' or 'RS')."""
        meta = self.get_metadata(ligand)
        chir = chirality.upper().strip()
        if chir in ("RR", "HOMOCHIRAL", "HOMO"):
            filename = meta.rr_file
        elif chir in ("RS", "HETEROCHIRAL", "HETERO", "MESO"):
            filename = meta.rs_file
        else:
            raise ValueError(f"Chirality must be 'RR' or 'RS', got '{chirality}'.")

        path = os.path.join(self.structures_dir, filename)
        if not os.path.isfile(path):
            raise FileNotFoundError(f"Structure coordinate file not found: {path}")
        return path

    def load_structure(self, ligand: str, chirality: str) -> XYZStructure:
        """Load and validate the 3D XYZ structure of a gallery member."""
        path = self.get_structure_path(ligand, chirality)
        validate_xyz_file(path)
        return read_xyz(path)

    def generate_input_text(
        self,
        ligand: str,
        chirality: str,
        software: str = "orca",
        functional: str = "B3LYP",
        basis: str = "def2-TZVP",
        dispersion: str = "D4",
        solvent: Optional[str] = "toluene",
        n_cores: int = 8,
        memory_mb: int = 4000,
    ) -> str:
        """
        Generate electronic structure calculation input text.

        Parameters
        ----------
        ligand : str
            Ligand key ('daib', 'binol', 'pybox').
        chirality : str
            Chirality state ('RR' or 'RS').
        software : str
            Target quantum chemistry package ('orca', 'gaussian', 'qchem').
        functional : str
            DFT exchange-correlation functional.
        basis : str
            Gaussian orbital basis set.
        dispersion : str
            Empirical dispersion correction ('D4', 'D3BJ', 'D3').
        solvent : Optional[str]
            Implicit continuum dielectric solvent (e.g. 'toluene', 'dichloromethane', 'THF').
        n_cores : int
            Number of CPU cores for parallelization.
        memory_mb : int
            Memory allocation per core in MB.

        Returns
        -------
        str
            Complete, formatted quantum chemistry input file content.
        """
        meta = self.get_metadata(ligand)
        struct = self.load_structure(ligand, chirality)
        sw = software.lower().strip()

        # Format coordinate block
        coord_lines = [
            f"  {atom:<4} {coord[0]:>14.8f} {coord[1]:>14.8f} {coord[2]:>14.8f}"
            for atom, coord in zip(struct.atoms, struct.coordinates, strict=True)
        ]
        coords_str = "\n".join(coord_lines)

        chg = meta.default_charge
        mult = meta.default_multiplicity

        if sw == "orca":
            solv_str = f" ! CPCM({solvent})" if solvent else ""
            header = (
                f"! {functional} {basis} {dispersion} Opt Freq{solv_str}\n"
                f"%pal nprocs {n_cores} end\n"
                f"%maxcore {memory_mb}\n"
                f"* xyz {chg} {mult}\n"
            )
            footer = "*\n"
            return header + coords_str + "\n" + footer

        elif sw in ("gaussian", "g16", "g09"):
            disp_str = f" empiricaldispersion=g{dispersion.lower()}" if dispersion else ""
            solv_str = f" scrf=(smd,solvent={solvent})" if solvent else ""
            route = f"#p {functional}/{basis}{disp_str}{solv_str} opt freq"
            header = (
                f"%nprocshared={n_cores}\n"
                f"%mem={memory_mb * n_cores}MB\n"
                f"{route}\n\n"
                f"{meta.full_name} [{chirality.upper()}] geometry optimization\n\n"
                f"{chg} {mult}\n"
            )
            return header + coords_str + "\n\n"

        elif sw == "qchem":
            solv_rem = "  solvent_method pcm\n" if solvent else ""
            pcm_block = f"\n$pcm\n  theory cpcm\n  solvent {solvent}\n$end\n" if solvent else ""
            text = (
                f"$rem\n"
                f"  jobtype opt\n"
                f"  method {functional.lower()}\n"
                f"  basis {basis.lower()}\n"
                f"  dft_d {dispersion.lower()}\n"
                f"{solv_rem}"
                f"$end\n\n"
                f"$molecule\n"
                f"  {chg} {mult}\n"
                f"{coords_str}\n"
                f"$end\n"
                f"{pcm_block}"
            )
            return text
        else:
            raise ValueError(f"Unsupported software package '{software}'. Choose 'orca', 'gaussian', or 'qchem'.")

    def export_pair(
        self,
        ligand: str,
        output_dir: str,
        software: str = "orca",
        **kwargs,
    ) -> Tuple[str, str]:
        """
        Export both homochiral (RR) and heterochiral (RS) input files to output_dir.

        Returns
        -------
        Tuple[str, str]
            Paths to generated (rr_input_path, rs_input_path).
        """
        os.makedirs(output_dir, exist_ok=True)
        ext = "inp" if software.lower() == "orca" else "gjf" if "gauss" in software.lower() else "in"

        rr_path = os.path.join(output_dir, f"{ligand}_RR.{ext}")
        rs_path = os.path.join(output_dir, f"{ligand}_RS.{ext}")

        rr_text = self.generate_input_text(ligand, "RR", software=software, **kwargs)
        rs_text = self.generate_input_text(ligand, "RS", software=software, **kwargs)

        with open(rr_path, "w", encoding="utf-8") as f:
            f.write(rr_text)
        with open(rs_path, "w", encoding="utf-8") as f:
            f.write(rs_text)

        return rr_path, rs_path
