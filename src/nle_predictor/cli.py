"""
Unified Command Line Interface for NLE-Predictor.
Supports quantum chemistry workflows, rapid analytical simulations,
experimental curve fitting, and mathematical diagnostics.
"""

import argparse
import json
import logging
import os
import sys
from typing import List, Optional
import numpy as np

from nle_predictor.core.constants import DEFAULT_TEMP
from nle_predictor.core.xyz import validate_xyz_file
from nle_predictor.fitting.experimental_fit import fit_kagan_reservoir, load_experimental_csv
from nle_predictor.models.diagnostics import analyze_nle_diagnostics
from nle_predictor.models.kagan_active import KaganActiveDimerModel
from nle_predictor.models.kagan_ml2 import (
    NLECalculator,
    write_sensitivity_report,
)
from nle_predictor.qc.crest import ConformerSampler
from nle_predictor.qc.orca import (
    Config,
    DFTEvaluator,
    collect_software_versions,
)
from nle_predictor.thermo.boltzmann import boltzmann_weighted_free_energy
from nle_predictor.viz.plotting import (
    plot_monomer_speciation,
    plot_nle_with_experimental,
)

logger = logging.getLogger("NLE_Predictor")


def build_parser() -> argparse.ArgumentParser:
    """Build the top-level argument parser with subcommands."""
    parser = argparse.ArgumentParser(
        prog="nle-predictor",
        description="Automated Non-Linear Effect (NLE) Predictors and Thermodynamic Intelligence for Asymmetric Catalysis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    subparsers = parser.add_subparsers(dest="subcommand", help="Available subcommands")

    # 1. RUN Subcommand (Quantum chemistry workflow)
    run_parser = subparsers.add_parser(
        "run",
        help="Run full conformer sampling and DFT workflow",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    _add_run_arguments(run_parser)

    # 2. SIMULATE Subcommand (Rapid analytical/numerical simulation)
    sim_parser = subparsers.add_parser(
        "simulate",
        help="Simulate NLE curve directly from thermodynamic Delta G or K constants",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    sim_parser.add_argument("--delta-g", type=float, default=-1.5, help="Delta G (RS - RR) in kcal/mol")
    sim_parser.add_argument("--khomo", type=float, default=1e5, help="Homochiral association constant")
    sim_parser.add_argument("--temp", type=float, default=DEFAULT_TEMP, help="Temperature in Kelvin")
    sim_parser.add_argument("--ee0", type=float, default=1.0, help="Intrinsic monomer enantioselectivity (0.0 to 1.0)")
    sim_parser.add_argument("--g", type=float, default=0.0, help="Relative dimer reactivity (k_dimer / k_monomer)")
    sim_parser.add_argument("--outdir", type=str, default="sim_outputs", help="Output directory")

    # 3. FIT Subcommand (Experimental data fitting)
    fit_parser = subparsers.add_parser(
        "fit",
        help="Fit experimental (ee_cat, ee_prod) CSV to extract thermodynamic parameters",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    fit_parser.add_argument("--data", type=str, required=True, help="Path to experimental CSV file")
    fit_parser.add_argument("--temp", type=float, default=DEFAULT_TEMP, help="Reaction temperature in Kelvin")
    fit_parser.add_argument("--khomo", type=float, default=1e5, help="Assumed K_homo constant")
    fit_parser.add_argument("--fit-ee0", action="store_true", help="Simultaneously fit intrinsic monomer selectivity ee_0")
    fit_parser.add_argument("--outdir", type=str, default="fit_outputs", help="Output directory")

    # 4. DIAGNOSTICS Subcommand
    diag_parser = subparsers.add_parser(
        "diagnostics",
        help="Compute mathematical curvature, initial slope, and amplification diagnostics",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    diag_parser.add_argument("--data", type=str, required=True, help="Path to NLE CSV data file")

    return parser


def _add_run_arguments(parser: argparse.ArgumentParser) -> None:
    """Helper to attach standard quantum workflow arguments."""
    parser.add_argument("--rr", type=str, default="inputs/dimer_RR.xyz", help="Path to homochiral (R,R) dimer XYZ file")
    parser.add_argument("--rs", type=str, default="inputs/dimer_RS.xyz", help="Path to heterochiral (R,S) dimer XYZ file")
    parser.add_argument("--temp", type=float, default=DEFAULT_TEMP, help="Temperature in Kelvin")
    parser.add_argument("--solvent", type=str, default=None, help="Solvent name for CREST and ORCA CPCM")
    parser.add_argument("--khomo", type=float, default=1e5, help="Assumed absolute dimerization constant for RR")
    parser.add_argument("--charge", type=int, default=0, help="Molecular charge for ORCA")
    parser.add_argument("--multiplicity", type=int, default=1, help="Spin multiplicity for ORCA")
    parser.add_argument("--functional", type=str, default="B3LYP", help="DFT Functional")
    parser.add_argument("--basis", type=str, default="def2-SVP", help="DFT Basis Set")
    parser.add_argument("--dispersion", type=str, default="D4", help="DFT Dispersion Correction")
    parser.add_argument("--cores", type=int, default=4, help="Number of CPU cores for CREST and ORCA")
    parser.add_argument("--conformers", type=int, default=1, help="Number of low-energy conformers per dimer")
    parser.add_argument("--outdir", type=str, default="calc_outputs", help="Directory to store calculations and results")
    parser.add_argument("--mock", action="store_true", help="Force mock execution (skip external binaries)")
    parser.add_argument("--dry-run", action="store_true", help="Preview planned chemistry commands without running")
    parser.add_argument("--sensitivity", action="store_true", help="Generate surrogate parameter sensitivity table")


def execute_simulate(args: argparse.Namespace) -> None:
    """Execute rapid analytical NLE simulation."""
    os.makedirs(args.outdir, exist_ok=True)
    ee_cat = np.linspace(0.0, 1.0, 101)

    if args.g > 0:
        model = KaganActiveDimerModel(temp=args.temp, g=args.g, ee_0=args.ee0)
        res = model.simulate(ee_cat, delta_g_kcal=args.delta_g, k_homo=args.khomo)
    else:
        calc = NLECalculator(temp=args.temp)
        res = calc.simulate(ee_cat, delta_g_kcal=args.delta_g, k_homo=args.khomo, ee_0=args.ee0)

    # Save CSV
    csv_path = os.path.join(args.outdir, "simulated_nle.csv")
    calc = NLECalculator(temp=args.temp)
    calc.export_data(res.ee_catalyst, res.ee_product, filepath=csv_path)

    # Plot
    plot_path = os.path.join(args.outdir, "simulated_curve.png")
    plot_nle_with_experimental(res, save_path=plot_path)

    # Speciation
    spec_path = os.path.join(args.outdir, "monomer_speciation.png")
    plot_monomer_speciation(res, save_path=spec_path)

    # Diagnostics
    diag = analyze_nle_diagnostics(res.ee_catalyst, res.ee_product)

    print(f"\n--- NLE Simulation Summary ---")
    print(f"Model: {res.model_name}")
    print(f"Classification: {res.classification}")
    print(f"Initial Slope at Racemic: {diag.initial_slope:.3f}")
    print(f"Max Amplification: {diag.max_amplification:.3f} (at ee_cat = {diag.catalyst_ee_at_max_amplification*100:.1f}%)")
    print(f"Area Under Asymmetry: {diag.area_under_asymmetry:.4f}")
    print(f"Results saved to: {args.outdir}\n")


def execute_fit(args: argparse.Namespace) -> None:
    """Execute experimental data fitting."""
    os.makedirs(args.outdir, exist_ok=True)
    cat_exp, prod_exp = load_experimental_csv(args.data)

    fit_res = fit_kagan_reservoir(
        cat_exp,
        prod_exp,
        temp=args.temp,
        k_homo=args.khomo,
        fit_ee_0=args.fit_ee0,
    )

    print("\n" + "=" * 50)
    print(" EXPERIMENTAL NLE REGRESSION RESULTS")
    print("=" * 50)
    print(f"Data Source:           {args.data} ({len(cat_exp)} points)")
    print(f"Fitted Delta G (RS-RR): {fit_res.delta_g_kcal:.4f} ± {fit_res.delta_g_stderr:.4f} kcal/mol")
    print(f"Fitted K_hetero/K_homo: {fit_res.k_ratio:.4f} ± {fit_res.k_ratio_stderr:.4f}")
    if args.fit_ee0:
        print(f"Fitted Intrinsic ee_0: {fit_res.metadata['fitted_ee_0']*100:.2f}%")
    print(f"Coefficient of Det (R²): {fit_res.r_squared:.4f}")
    print(f"Root Mean Square Error: {fit_res.rmse*100:.2f}%")
    print("=" * 50)

    # Simulate dense fitted curve
    cat_dense = np.linspace(0.0, 1.0, 101)
    calc = NLECalculator(temp=args.temp)
    ee_0_val = fit_res.metadata["fitted_ee_0"] if args.fit_ee0 else 1.0
    sim_result = calc.simulate(cat_dense, delta_g_kcal=fit_res.delta_g_kcal, k_homo=args.khomo, ee_0=ee_0_val)

    # Plot fitted curve with experimental points and residuals
    fig_path = os.path.join(args.outdir, "fitted_nle_curve.png")
    plot_nle_with_experimental(
        sim_result,
        exp_cat=cat_exp,
        exp_prod=prod_exp,
        save_path=fig_path,
        title=rf"Experimental Fit: $\Delta G = {fit_res.delta_g_kcal:.2f}$ kcal/mol ($R^2 = {fit_res.r_squared:.3f}$)",
    )

    # Save summary JSON
    summary_path = os.path.join(args.outdir, "fit_summary.json")
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump({
            "delta_g_kcal": fit_res.delta_g_kcal,
            "delta_g_stderr": fit_res.delta_g_stderr,
            "k_ratio": fit_res.k_ratio,
            "k_ratio_stderr": fit_res.k_ratio_stderr,
            "r_squared": fit_res.r_squared,
            "rmse": fit_res.rmse,
            "fitted_ee_0": fit_res.metadata.get("fitted_ee_0", 1.0),
        }, f, indent=2)

    print(f"Fitted plot and regression summary saved to: {args.outdir}\n")


def execute_diagnostics(args: argparse.Namespace) -> None:
    """Compute and print curve diagnostics from a CSV file."""
    cat, prod = load_experimental_csv(args.data)
    diag = analyze_nle_diagnostics(cat, prod)

    print("\n" + "=" * 50)
    print(" NLE CURVE MATHEMATICAL DIAGNOSTICS")
    print("=" * 50)
    print(f"Classification:            {diag.classification}")
    print(f"Initial Slope at Racemic:  {diag.initial_slope:.4f}")
    print(f"Max Amplification Factor:  {diag.max_amplification:.4f}")
    print(f"Catalyst ee at Max Amp:    {diag.catalyst_ee_at_max_amplification*100:.2f}%")
    print(f"Area Under Asymmetry (AUA):{diag.area_under_asymmetry:.4f}")
    print(f"Is Hyperbolic (+NLE):      {diag.is_hyperbolic}")
    print(f"Is Sigmoidal:              {diag.is_sigmoidal}")
    if diag.inflection_points:
        print(f"Inflection Point(s):       {diag.inflection_points}")
    print("=" * 50 + "\n")


def main(argv: Optional[List[str]] = None) -> None:
    """Main CLI entry point."""
    if argv is None:
        argv = sys.argv[1:]

    parser = build_parser()

    # If first argument is not a subcommand and arguments look like standard flags
    # (e.g. `nle_predictor.py --mock` or `nle_predictor.py --rr ...`), default to 'run'
    subcommands = {"run", "simulate", "fit", "diagnostics"}
    if argv and argv[0] not in subcommands and not argv[0].startswith("-h") and not argv[0] == "--help":
        argv = ["run"] + argv
    elif not argv:
        # Default to showing help
        parser.print_help()
        return

    args = parser.parse_args(argv)

    if args.subcommand == "simulate":
        execute_simulate(args)
    elif args.subcommand == "fit":
        execute_fit(args)
    elif args.subcommand == "diagnostics":
        execute_diagnostics(args)
    elif args.subcommand == "run":
        # Run quantum workflow via root nle_predictor
        import nle_predictor
        # If nle_predictor is package or script
        if hasattr(nle_predictor, "main"):
            # run legacy main
            pass
        from nle_predictor.models.kagan_ml2 import write_sensitivity_report
        # Execute workflow
        _execute_qc_workflow(args)


def _execute_qc_workflow(args: argparse.Namespace) -> None:
    """Core QC workflow orchestrator."""
    import shutil
    from datetime import datetime, timezone

    logger.info("=" * 55)
    logger.info(" Automated Non-Linear Effect (NLE) Predictor")
    logger.info("=" * 55)

    missing_deps = Config.check_dependencies()
    if missing_deps and not args.mock and not args.dry_run:
        raise SystemExit(
            "Missing computational binaries: " + ", ".join(missing_deps) +
            ". Install them or rerun with --mock for a synthetic demonstration."
        )

    if args.temp <= 0:
        raise SystemExit("--temp must be greater than zero")
    if args.sensitivity and args.temp <= 20:
        raise SystemExit("--temp must be greater than 20 K when using --sensitivity")
    if args.cores <= 0:
        raise SystemExit("--cores must be greater than zero")
    if args.conformers <= 0:
        raise SystemExit("--conformers must be greater than zero")
    if args.khomo < 0:
        raise SystemExit("--khomo must be non-negative")
    if args.multiplicity <= 0:
        raise SystemExit("--multiplicity must be greater than zero")

    if not os.path.exists(args.rr) or not os.path.exists(args.rs):
        if not args.mock:
            raise SystemExit("Both --rr and --rs must point to existing XYZ files")
        logger.info("Mock input files not found. Generating dummy XYZ files for testing.")
        for p in [args.rr, args.rs]:
            parent = os.path.dirname(p)
            if parent:
                os.makedirs(parent, exist_ok=True)
            with open(p, "w", encoding="utf-8") as f:
                f.write("2\nDummy\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74\n")

    validate_xyz_file(args.rr)
    validate_xyz_file(args.rs)

    if args.dry_run:
        logger.info("Dry run: no calculations or output files will be created.")
        return

    os.makedirs(args.outdir, exist_ok=True)

    # 1. Conformational Sampling
    sampler = ConformerSampler(cores=args.cores, mock_mode=args.mock)
    sampled_rr = sampler.run_sampling(args.rr, output_dir=args.outdir, solvent=args.solvent, max_conformers=args.conformers)
    sampled_rs = sampler.run_sampling(args.rs, output_dir=args.outdir, solvent=args.solvent, max_conformers=args.conformers)
    rr_confs = sampled_rr if isinstance(sampled_rr, list) else [sampled_rr]
    rs_confs = sampled_rs if isinstance(sampled_rs, list) else [sampled_rs]

    # 2. DFT Free Energies
    dft = DFTEvaluator(cores=args.cores, mock_mode=args.mock)
    rr_energies = []
    rs_energies = []
    for label, confs, energies in (("RR", rr_confs, rr_energies), ("RS", rs_confs, rs_energies)):
        for index, conf in enumerate(confs, start=1):
            run_name = f"opt_{label}_conf{index:03d}"
            input_name = dft.generate_input(
                conf, run_name, args.functional, args.basis, args.dispersion,
                args.solvent, args.charge, args.multiplicity,
            )
            input_path = os.path.join(args.outdir, input_name)
            shutil.move(input_name, input_path)
            output_path = dft.run_dft(input_name, args.outdir)
            energies.append(dft.extract_free_energy(output_path))

    g_rr = boltzmann_weighted_free_energy(rr_energies, args.temp)
    g_rs = boltzmann_weighted_free_energy(rs_energies, args.temp)

    # 3. NLE Simulation
    calc = NLECalculator(temp=args.temp)
    delta_g = calc.calc_equilibrium_constant(g_rr, g_rs)
    ee_cat = np.linspace(0.0, 1.0, 101)
    ee_prod = calc.simulate_kagan_model(ee_cat, delta_g, k_homo_abs=args.khomo)

    # Export & Plot
    calc.export_data(ee_cat, ee_prod, filepath=os.path.join(args.outdir, "nle_results.csv"))
    calc.plot_nle_curve(ee_cat, ee_prod, save_path=os.path.join(args.outdir, "nle_curve.png"))

    # Speciation Plot
    sim_res = calc.simulate(ee_cat, delta_g_kcal=delta_g, k_homo=args.khomo)
    plot_monomer_speciation(sim_res, save_path=os.path.join(args.outdir, "speciation_curve.png"))

    if args.sensitivity:
        write_sensitivity_report(args.outdir, ee_cat, delta_g, args.temp, args.khomo)

    logger.info("Workflow completed successfully!")


if __name__ == "__main__":
    main()
