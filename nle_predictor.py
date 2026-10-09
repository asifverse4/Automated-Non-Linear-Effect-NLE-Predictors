import os
import subprocess
import re
import csv
import hashlib
import json
import logging
import argparse
import shutil
import platform
from datetime import datetime, timezone
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root

# Configuration and Thermodynamic Constants
HARTREE_TO_KCAL = 627.509
GAS_CONSTANT_R = 0.00198720425864083  # kcal/(K*mol)
DEFAULT_TEMP = 298.15 # Kelvin

def setup_logger():
    """Sets up a professional logger for the module."""
    logger = logging.getLogger("NLE_Predictor")
    logger.setLevel(logging.INFO)
    if not logger.handlers:
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - [%(levelname)s] - %(message)s', datefmt='%H:%M:%S')
        ch.setFormatter(formatter)
        logger.addHandler(ch)
    return logger

logger = setup_logger()


def validate_xyz_file(xyz_path):
    """Validate an XYZ file and return its atom count."""
    if not os.path.isfile(xyz_path):
        raise FileNotFoundError(f"XYZ input file not found: {xyz_path}")

    with open(xyz_path, 'r') as xyz_file:
        lines = xyz_file.readlines()

    if len(lines) < 2:
        raise ValueError(f"XYZ file is incomplete: {xyz_path}")

    try:
        atom_count = int(lines[0].strip())
    except ValueError as error:
        raise ValueError(f"First XYZ line must be an atom count: {xyz_path}") from error

    if atom_count <= 0 or len(lines) < atom_count + 2:
        raise ValueError(f"XYZ file does not contain {atom_count} atoms: {xyz_path}")

    for line_number, line in enumerate(lines[2:atom_count + 2], start=3):
        fields = line.split()
        if len(fields) < 4:
            raise ValueError(f"Invalid XYZ coordinate at {xyz_path}:{line_number}")
        try:
            [float(value) for value in fields[1:4]]
        except ValueError as error:
            raise ValueError(f"Invalid XYZ coordinate at {xyz_path}:{line_number}") from error

    return atom_count


def sha256_file(file_path):
    """Return the SHA-256 digest for a file."""
    digest = hashlib.sha256()
    with open(file_path, 'rb') as input_file:
        for chunk in iter(lambda: input_file.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def write_run_manifest(output_dir, args, input_files, missing_dependencies):
    """Write machine-readable provenance for a workflow run."""
    manifest = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "arguments": vars(args),
        "mock_mode": bool(args.mock),
        "missing_dependencies": missing_dependencies,
        "input_files": {
            path: {"sha256": sha256_file(path), "size_bytes": os.path.getsize(path)}
            for path in input_files
        },
    }
    manifest_path = os.path.join(output_dir, "run_manifest.json")
    with open(manifest_path, 'w') as manifest_file:
        json.dump(manifest, manifest_file, indent=2)
        manifest_file.write('\n')
    logger.info(f"Run manifest saved to {manifest_path}")


def write_prediction_report(
    output_dir,
    args,
    g_rr,
    g_rs,
    delta_g_kcal,
    k_homo,
    k_hetero,
    ee_cat,
    ee_prod,
):
    """Write a structured summary of the thermodynamic prediction."""
    product_at_50 = float(np.interp(0.5, ee_cat, ee_prod))
    amplification = np.divide(
        ee_prod[1:], ee_cat[1:], out=np.zeros_like(ee_prod[1:]), where=ee_cat[1:] > 0
    )
    difference = ee_prod - ee_cat
    tolerance = 1e-6
    if np.max(difference) > tolerance:
        classification = "positive NLE"
    elif np.min(difference) < -tolerance:
        classification = "negative NLE"
    else:
        classification = "near-linear"

    report = {
        "model": "Kagan ML2 reservoir",
        "temperature_K": args.temp,
        "electronic_state": {
            "charge": args.charge,
            "multiplicity": args.multiplicity,
        },
        "electronic_structure": {
            "functional": args.functional,
            "basis": args.basis,
            "dispersion": args.dispersion,
            "solvent": args.solvent,
        },
        "free_energy": {
            "RR_hartree": g_rr,
            "RS_hartree": g_rs,
            "RS_minus_RR_kcal_per_mol": delta_g_kcal,
        },
        "association_constants": {
            "K_homo": k_homo,
            "K_hetero": k_hetero,
            "K_hetero_over_K_homo": k_hetero / k_homo if k_homo else None,
        },
        "nle_assessment": {
            "classification": classification,
            "product_ee_at_50_percent_catalyst_ee_percent": product_at_50 * 100,
            "maximum_product_ee_percent": float(np.max(ee_prod) * 100),
            "maximum_ee_amplification": float(np.max(amplification)),
        },
        "data_file": "nle_results.csv",
        "plot_file": "nle_curve.png",
    }
    report_path = os.path.join(output_dir, "nle_report.json")
    with open(report_path, 'w') as report_file:
        json.dump(report, report_file, indent=2)
        report_file.write('\n')

    logger.info("NLE assessment: %s", classification)
    logger.info("Predicted product ee at 50%% catalyst ee: %.2f%%", product_at_50 * 100)
    logger.info(f"Prediction report saved to {report_path}")

class Config:
    """Configuration and dependency management for quantum chemistry binaries."""
    CREST_CMD = "crest"
    ORCA_CMD = "orca"
    
    @classmethod
    def check_dependencies(cls):
        """Checks if required computational chemistry software is in the system PATH."""
        missing = []
        if shutil.which(cls.CREST_CMD) is None: missing.append('CREST')
        if shutil.which(cls.ORCA_CMD) is None: missing.append('ORCA')
        return missing

class ConformerSampler:
    """
    Handles automated conformational sampling using CREST (xtb).
    Identifies the lowest energy conformer for a given chiral complex.
    """
    def __init__(self, crest_path=Config.CREST_CMD, cores=4, mock_mode=False):
        self.crest_path = crest_path
        self.cores = cores
        self.mock_mode = mock_mode

    def run_sampling(self, input_xyz, output_dir, solvent=None):
        """Runs CREST conformational search on a given XYZ file."""
        validate_xyz_file(input_xyz)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            
        logger.info(f"Starting conformer sampling for {input_xyz}...")
        
        base_name = os.path.splitext(os.path.basename(input_xyz))[0]
        run_dir = os.path.join(output_dir, f"crest_{base_name}")
        os.makedirs(run_dir, exist_ok=True)
        
        best_conformer_path = os.path.join(run_dir, "crest_best.xyz")
        
        if self.mock_mode:
            logger.info(f"[MOCK] Simulating CREST execution for {input_xyz}")
            with open(best_conformer_path, 'w') as f:
                f.write(f"2\nMock Conformer {base_name}\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74\n")
            return best_conformer_path

        # Construct CREST command
        abs_input = os.path.abspath(input_xyz)
        cmd = [self.crest_path, abs_input, f"-T", str(self.cores)]
        
        if solvent:
            # CREST ALPB model
            cmd.extend(["-alpb", solvent])
            
        try:
            logger.info(f"Executing: {' '.join(cmd)} in {run_dir}")
            with open(os.path.join(run_dir, "crest.out"), "w") as f_out:
                subprocess.run(cmd, cwd=run_dir, stdout=f_out, stderr=subprocess.STDOUT, check=True)
            logger.info(f"CREST sampling complete for {base_name}.")
            conformers_path = os.path.join(run_dir, "crest_conformers.xyz")
            if not os.path.exists(conformers_path):
                raise FileNotFoundError(
                    f"CREST completed but did not produce {conformers_path}"
                )

            # CREST writes conformers in increasing energy order. Keep the
            # first XYZ block as the reproducible input for the DFT step.
            with open(conformers_path, 'r') as conformers:
                lines = conformers.readlines()
            if len(lines) < 2:
                raise ValueError(f"Invalid CREST conformer file: {conformers_path}")
            atom_count = int(lines[0].strip())
            block_size = atom_count + 2
            if len(lines) < block_size:
                raise ValueError(f"Incomplete CREST conformer file: {conformers_path}")
            with open(best_conformer_path, 'w') as best_conformer:
                best_conformer.writelines(lines[:block_size])
            return best_conformer_path
        except subprocess.CalledProcessError as e:
            logger.error(f"CREST failed for {base_name}. Check {run_dir}/crest.out")
            raise

class DFTEvaluator:
    """
    Handles DFT evaluations using ORCA to compute final Gibbs Free Energies.
    """
    def __init__(self, orca_path=Config.ORCA_CMD, cores=4, mock_mode=False):
        self.orca_path = orca_path
        self.cores = cores
        self.mock_mode = mock_mode

    def generate_input(
        self,
        xyz_file,
        run_name,
        functional="B3LYP",
        basis="def2-SVP",
        dispersion="D4",
        solvent=None,
        charge=0,
        multiplicity=1,
    ):
        """Generates an ORCA input file for geometry optimization and frequency calculations."""
        validate_xyz_file(xyz_file)
        inp_path = f"{run_name}.inp"
        
        # Read coordinates
        coords = []
        with open(xyz_file, 'r') as f:
            lines = f.readlines()
            if len(lines) > 2:
                coords = lines[2:] # Skip atom count and comment line
        
        # Construct header
        header = f"! {functional} {basis} {dispersion} Opt Freq\n"
        if solvent:
            header += f"! CPCM({solvent})\n"
            
        if self.cores > 1:
            header += f"%pal nprocs {self.cores} end\n"
            
        with open(inp_path, 'w') as f:
            f.write(header)
            f.write(
                f"%coords\n  CTyp xyz\n  Charge {charge}\n  Mult {multiplicity}\n  coords\n"
            )
            for line in coords:
                if line.strip():
                    f.write(f"    {line.strip()}\n")
            f.write("  end\nend\n")
            
        return inp_path

    def extract_free_energy(self, out_file):
        """Parses ORCA output to extract the Final Gibbs Free Energy in Hartrees."""
        gibbs_energy = None
        energy_pattern = re.compile(
            r"Final Gibbs free energy\s+([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)"
        )
        with open(out_file, 'r') as f:
            for line in f:
                match = energy_pattern.search(line)
                if match:
                    gibbs_energy = float(match.group(1))
        
        if gibbs_energy is None:
            logger.error(f"Could not find Gibbs free energy in {out_file}")
            raise ValueError(f"Free energy extraction failed for {out_file}")
        return gibbs_energy

    def run_dft(self, inp_file, output_dir):
        """Executes the ORCA calculation."""
        out_file = inp_file.replace('.inp', '.out')
        out_path = os.path.join(output_dir, out_file)
        inp_path = os.path.join(output_dir, inp_file)
        
        logger.info(f"Running DFT evaluation: {inp_file}...")
        
        if self.mock_mode:
            logger.info(f"[MOCK] Simulating ORCA execution for {inp_file}")
            os.makedirs(output_dir, exist_ok=True)
            with open(out_path, 'w') as out_f:
                energy = -1000.00366 if "RS" in inp_file else -1000.00000
                out_f.write(f"Final Gibbs free energy             {energy:.6f} Eh\n")
            return out_path

        try:
            # Note: ORCA requires full path to the executable to run smoothly in some environments
            with open(out_path, 'w') as out_f:
                subprocess.run([self.orca_path, inp_path], cwd=output_dir, stdout=out_f, stderr=subprocess.STDOUT, check=True)
            logger.info(f"DFT evaluation complete for {inp_file}.")
            return out_path
        except subprocess.CalledProcessError:
            logger.error(f"ORCA calculation failed. Check {out_path}.")
            raise

class NLECalculator:
    """
    Calculates equilibrium constants and simulates NLE curves 
    using a rigorous numerical solver for the ML2 (Reservoir) Kagan Model.
    """
    def __init__(self, temp=DEFAULT_TEMP):
        self.temp = temp
        
    def calc_equilibrium_constant(self, g_rr, g_rs):
        """
        Calculates the thermodynamic parameters.
        Delta G = G_RS - G_RR
        The statistical factor of 2 is handled in the association constants later.
        """
        delta_g_hartree = g_rs - g_rr
        delta_g_kcal = delta_g_hartree * HARTREE_TO_KCAL
        
        logger.info(f"--- Thermodynamic Results ---")
        logger.info(f"G(RR) = {g_rr:.6f} Eh")
        logger.info(f"G(RS) = {g_rs:.6f} Eh")
        logger.info(f"Delta G (RS - RR): {delta_g_kcal:.4f} kcal/mol")
        
        return delta_g_kcal

    def solve_mass_balance(self, ee_cat, k_homo, k_hetero):
        """
        Numerically solves the mass balance equations for Monomers (M_R, M_S).
        Total Catalyst = C_R + C_S = 1.0
        """
        if not -1.0 <= ee_cat <= 1.0:
            raise ValueError("Catalyst ee must be between -1.0 and 1.0")
        if k_homo < 0 or k_hetero < 0:
            raise ValueError("Association constants must be non-negative")

        C_R = 0.5 * (1.0 + ee_cat)
        C_S = 0.5 * (1.0 - ee_cat)
        
        def equations(vars):
            M_R, M_S = vars
            # Mass balance for R monomer
            eq1 = M_R + 2 * k_homo * (M_R**2) + k_hetero * M_R * M_S - C_R
            # Mass balance for S monomer
            eq2 = M_S + 2 * k_homo * (M_S**2) + k_hetero * M_R * M_S - C_S
            return [eq1, eq2]
        
        # Initial guesses based on pure homochiral limit (ignoring heterochiral)
        # M + 2*K*M^2 - C = 0 => 2*K*M^2 + M - C = 0
        mr_guess = (-1 + np.sqrt(1 + 8 * k_homo * C_R)) / (4 * k_homo) if C_R > 0 else 0.0
        ms_guess = (-1 + np.sqrt(1 + 8 * k_homo * C_S)) / (4 * k_homo) if C_S > 0 else 0.0
        
        sol = root(equations, [mr_guess, ms_guess], method='hybr')

        residual = np.max(np.abs(equations(sol.x))) if sol.success else np.inf
        if (
            not sol.success
            or not np.all(np.isfinite(sol.x))
            or np.any(sol.x < -1e-10)
            or residual > 1e-8
        ):
            raise RuntimeError(
                f"Mass-balance solve failed at ee_cat={ee_cat}: "
                f"{sol.message}; residual={residual:.3e}"
            )

        return max(0.0, sol.x[0]), max(0.0, sol.x[1])

    def simulate_kagan_model(self, ee_cat_range, delta_g_kcal, k_homo_abs=1e5):
        """
        Simulates the NLE curve.
        k_homo_abs: Absolute association constant for R + R -> RR (Default 1e5 for reservoir limit)
        """
        if self.temp <= 0:
            raise ValueError("Temperature must be greater than zero")
        if k_homo_abs < 0:
            raise ValueError("k_homo_abs must be non-negative")

        # Calculate k_hetero relative to k_homo
        # K_RS / K_RR = 2 * exp(-DeltaG / RT)  --> factor of 2 is the symmetry number
        ratio = 2.0 * np.exp(-delta_g_kcal / (GAS_CONSTANT_R * self.temp))
        k_hetero_abs = k_homo_abs * ratio
        self.last_k_hetero = k_hetero_abs
        self.last_k_ratio = ratio
        
        logger.info(f"Assumed K_homo (RR binding): {k_homo_abs:.2e}")
        logger.info(f"Calculated K_hetero (RS binding): {k_hetero_abs:.2e}")
        logger.info(f"K_hetero / K_homo ratio: {ratio:.4f}")

        ee_prods = []
        for ee_cat in ee_cat_range:
            M_R, M_S = self.solve_mass_balance(ee_cat, k_homo_abs, k_hetero_abs)
            
            # ee of the product is dictated by the active monomer concentrations
            total_active = M_R + M_S
            if total_active <= 1e-12:
                ee_p = 0.0
            else:
                ee_p = (M_R - M_S) / total_active
                
            ee_prods.append(ee_p)
            
        return np.array(ee_prods)

    def export_data(self, ee_cat, ee_prod, filepath="nle_results.csv"):
        """Exports the raw NLE data to a CSV for external plotting."""
        with open(filepath, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(["ee_catalyst(%)", "ee_product(%)"])
            for cat, prod in zip(ee_cat * 100, ee_prod * 100):
                writer.writerow([f"{cat:.2f}", f"{prod:.2f}"])
        logger.info(f"Data successfully exported to {filepath}")

    def plot_nle_curve(self, ee_cat, ee_prod, save_path="nle_curve.png"):
        """Generates a publication-quality plot of the NLE curve."""
        plt.style.use('default')
        fig, ax = plt.subplots(figsize=(7, 6))
        
        # Ideal linear relationship
        ax.plot(ee_cat * 100, ee_cat * 100, color='gray', linestyle='--', linewidth=1.5, label='Ideal Linear Effect')
        
        # Predicted NLE
        ax.plot(ee_cat * 100, ee_prod * 100, color='#004488', linestyle='-', linewidth=2.5, label='Predicted NLE (ML$_2$ Model)')
        
        # Aesthetics
        ax.set_title('Predicted Non-Linear Effect', fontsize=16, fontweight='bold', pad=15)
        ax.set_xlabel('Catalyst $ee$ (%)', fontsize=14)
        ax.set_ylabel('Product $ee$ (%)', fontsize=14)
        
        ax.tick_params(axis='both', which='major', labelsize=12)
        ax.grid(True, linestyle=':', alpha=0.6, color='gray')
        
        ax.set_xlim(0, 100)
        ax.set_ylim(0, 100)
        
        # Format axes purely square
        ax.set_aspect('equal', adjustable='box')
        
        ax.legend(fontsize=12, loc='lower right', framealpha=1.0, edgecolor='black')
        
        plt.tight_layout()
        plt.savefig(save_path, dpi=600, bbox_inches='tight')
        logger.info(f"Publication-ready plot saved to {save_path}")

def setup_dummy_files(rr_path, rs_path):
    """Creates dummy XYZ files if the user didn't provide any, purely for testing."""
    for path in [rr_path, rs_path]:
        parent = os.path.dirname(path)
        if parent:
            os.makedirs(parent, exist_ok=True)
    dummy_xyz = "2\nDummy\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74\n"
    
    for path in [rr_path, rs_path]:
        if not os.path.exists(path):
            with open(path, 'w') as f:
                f.write(dummy_xyz)
            logger.info(f"Created dummy input file: {path}")

def main():
    parser = argparse.ArgumentParser(
        description="Automated Non-Linear Effect (NLE) Predictor for Asymmetric Catalysis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required inputs (made optional via default mock files for safety)
    parser.add_argument('--rr', type=str, default='inputs/dimer_RR.xyz', help="Path to homochiral (R,R) dimer XYZ file")
    parser.add_argument('--rs', type=str, default='inputs/dimer_RS.xyz', help="Path to heterochiral (R,S) dimer XYZ file")
    
    # Chemical Parameters
    parser.add_argument('--temp', type=float, default=298.15, help="Temperature in Kelvin")
    parser.add_argument('--solvent', type=str, default=None, help="Solvent name (e.g., toluene, water) for CREST and ORCA CPCM")
    parser.add_argument('--khomo', type=float, default=1e5, help="Assumed absolute dimerization constant for RR (ML2 reservoir)")
    parser.add_argument('--charge', type=int, default=0, help="Molecular charge for ORCA (default: 0)")
    parser.add_argument('--multiplicity', type=int, default=1, help="Spin multiplicity for ORCA (default: 1)")
    
    # Computational Details
    parser.add_argument('--functional', type=str, default="B3LYP", help="DFT Functional")
    parser.add_argument('--basis', type=str, default="def2-SVP", help="DFT Basis Set")
    parser.add_argument('--dispersion', type=str, default="D4", help="DFT Dispersion Correction")
    parser.add_argument('--cores', type=int, default=4, help="Number of CPU cores for CREST and ORCA")
    
    # Workflow directives
    parser.add_argument('--outdir', type=str, default="calc_outputs", help="Directory to store intermediate and final calculations")
    parser.add_argument('--mock', action='store_true', help="Force mock execution (skip actual CREST/ORCA runs)")
    parser.add_argument('--dry-run', action='store_true', help="Validate inputs and print the planned commands without running chemistry")
    
    args = parser.parse_args()

    logger.info("="*55)
    logger.info(" Automated Non-Linear Effect (NLE) Predictor")
    logger.info("="*55)
    
    # 0. Check dependencies
    missing_deps = Config.check_dependencies()
    if missing_deps and not args.mock and not args.dry_run:
        parser.error(
            "Missing computational binaries: " + ", ".join(missing_deps) +
            ". Install them or rerun with --mock for a synthetic demonstration."
        )
    mock_mode = args.mock

    if args.temp <= 0:
        parser.error("--temp must be greater than zero")
    if args.cores <= 0:
        parser.error("--cores must be greater than zero")
    if args.khomo < 0:
        parser.error("--khomo must be non-negative")
    if args.multiplicity <= 0:
        parser.error("--multiplicity must be greater than zero")
    
    # 1. Setup Inputs
    if not os.path.exists(args.rr) or not os.path.exists(args.rs):
        if not args.mock:
            parser.error("Both --rr and --rs must point to existing XYZ files")
        logger.info("Mock input files not found. Generating dummy XYZ files for testing.")
        setup_dummy_files(args.rr, args.rs)

    validate_xyz_file(args.rr)
    validate_xyz_file(args.rs)

    if args.dry_run:
        solvent_flag = f" -alpb {args.solvent}" if args.solvent else ""
        logger.info("Dry run: no calculations or output files will be created.")
        logger.info(
            "CREST RR: %s %s -T %s%s",
            Config.CREST_CMD,
            os.path.abspath(args.rr),
            args.cores,
            solvent_flag,
        )
        logger.info(
            "CREST RS: %s %s -T %s%s",
            Config.CREST_CMD,
            os.path.abspath(args.rs),
            args.cores,
            solvent_flag,
        )
        logger.info(
            "ORCA: %s [generated opt_RR.inp and opt_RS.inp] (%s/%s, charge=%s, multiplicity=%s)",
            Config.ORCA_CMD,
            args.functional,
            args.basis,
            args.charge,
            args.multiplicity,
        )
        return

    os.makedirs(args.outdir, exist_ok=True)
    
    # 2. Conformer Sampling (CREST)
    logger.info("\n--- [Step 1] Conformational Sampling ---")
    sampler = ConformerSampler(cores=args.cores, mock_mode=mock_mode)
    best_rr = sampler.run_sampling(args.rr, output_dir=args.outdir, solvent=args.solvent)
    best_rs = sampler.run_sampling(args.rs, output_dir=args.outdir, solvent=args.solvent)
    
    # 3. DFT Free Energy Evaluation (ORCA)
    logger.info("\n--- [Step 2] DFT Free Energy Evaluation ---")
    dft = DFTEvaluator(cores=args.cores, mock_mode=mock_mode)
    
    # Generate inputs
    inp_rr_name = dft.generate_input(
        best_rr,
        "opt_RR",
        args.functional,
        args.basis,
        args.dispersion,
        args.solvent,
        args.charge,
        args.multiplicity,
    )
    inp_rs_name = dft.generate_input(
        best_rs,
        "opt_RS",
        args.functional,
        args.basis,
        args.dispersion,
        args.solvent,
        args.charge,
        args.multiplicity,
    )
    
    # Move them to working directory
    inp_rr_path = os.path.join(args.outdir, inp_rr_name)
    inp_rs_path = os.path.join(args.outdir, inp_rs_name)
    shutil.move(inp_rr_name, inp_rr_path)
    shutil.move(inp_rs_name, inp_rs_path)
    
    # Execute ORCA
    out_rr = dft.run_dft(inp_rr_name, args.outdir)
    out_rs = dft.run_dft(inp_rs_name, args.outdir)
    
    # Extract Data
    g_rr = dft.extract_free_energy(out_rr)
    g_rs = dft.extract_free_energy(out_rs)
    
    # 4. Thermodynamic & NLE Calculation
    logger.info("\n--- [Step 3] NLE Prediction & Visualization ---")
    nle_calc = NLECalculator(temp=args.temp)
    delta_g = nle_calc.calc_equilibrium_constant(g_rr, g_rs)
    
    # Simulate data points (0% to 100% ee)
    ee_cat_range = np.linspace(0, 1.0, 101)
    ee_prod_range = nle_calc.simulate_kagan_model(ee_cat_range, delta_g, k_homo_abs=args.khomo)
    
    # 5. Export and Plot
    nle_calc.export_data(
        ee_cat_range,
        ee_prod_range,
        filepath=os.path.join(args.outdir, "nle_results.csv")
    )
    nle_calc.plot_nle_curve(
        ee_cat_range,
        ee_prod_range,
        save_path=os.path.join(args.outdir, "nle_curve.png")
    )
    write_prediction_report(
        args.outdir,
        args,
        g_rr,
        g_rs,
        delta_g,
        args.khomo,
        nle_calc.last_k_hetero,
        ee_cat_range,
        ee_prod_range,
    )
    write_run_manifest(args.outdir, args, [args.rr, args.rs], missing_deps)
    
    logger.info("="*55)
    logger.info("Workflow completed successfully!")
    logger.info("="*55)

if __name__ == "__main__":
    main()