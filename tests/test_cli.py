"""Tests for the unified NLE-Predictor CLI interface in nle_predictor.cli."""

import os

from nle_predictor.cli import main


def test_cli_simulate(tmp_path):
    outdir = tmp_path / "sim"
    main(["simulate", "--delta-g", "-1.2", "--outdir", str(outdir)])

    assert os.path.isfile(outdir / "simulated_nle.csv")
    assert os.path.isfile(outdir / "simulated_curve.png")
    assert os.path.isfile(outdir / "monomer_speciation.png")


def test_cli_fit(tmp_path):
    csv_file = tmp_path / "data.csv"
    csv_file.write_text("ee_cat,ee_prod\n10,35\n30,70\n50,88\n70,95\n90,99\n")

    outdir = tmp_path / "fit"
    main(["fit", "--data", str(csv_file), "--outdir", str(outdir)])

    assert os.path.isfile(outdir / "fitted_nle_curve.png")
    assert os.path.isfile(outdir / "fit_summary.json")


def test_cli_diagnostics(tmp_path, capsys):
    csv_file = tmp_path / "data.csv"
    csv_file.write_text("ee_cat,ee_prod\n0,0\n20,50\n50,85\n80,96\n100,100\n")

    main(["diagnostics", "--data", str(csv_file)])
    captured = capsys.readouterr()
    assert "NLE CURVE MATHEMATICAL DIAGNOSTICS" in captured.out
    assert "Initial Slope at Racemic" in captured.out


def test_cli_run_mock(tmp_path):
    outdir = tmp_path / "calc"
    rr_file = tmp_path / "rr.xyz"
    rs_file = tmp_path / "rs.xyz"
    content = "2\nTest\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74\n"
    rr_file.write_text(content)
    rs_file.write_text(content)

    main([
        "run",
        "--rr", str(rr_file),
        "--rs", str(rs_file),
        "--outdir", str(outdir),
        "--mock",
    ])

    assert os.path.isfile(outdir / "nle_results.csv")
    assert os.path.isfile(outdir / "nle_curve.png")
    assert os.path.isfile(outdir / "speciation_curve.png")


def test_cli_empty_args_or_interactive_in_non_tty(capsys):
    main([])
    captured = capsys.readouterr()
    assert "nle-predictor" in captured.out
    assert "usage:" in captured.out


def test_pyinstaller_spec_configuration():
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    spec_path = os.path.join(repo_root, "nle_predictor.spec")
    assert os.path.isfile(spec_path)
    with open(spec_path, "r", encoding="utf-8") as f:
        content = f.read()
    assert "src/nle_predictor/__main__.py" in content
    assert "nle-predictor" in content


def test_cli_dashboard(tmp_path):
    dash_out = tmp_path / "custom_dash.html"
    main(["dashboard", "--output", str(dash_out)])

    assert os.path.isfile(dash_out)
    assert os.path.getsize(dash_out) > 5000


def test_cli_oligomer(tmp_path):
    outdir = tmp_path / "oligomer"
    main(["oligomer", "--order", "3", "--delta-g", "-1.5", "--outdir", str(outdir)])

    assert os.path.isfile(outdir / "oligomer_order3_nle.csv")
    assert os.path.isfile(outdir / "oligomer_order3_curve.png")


def test_cli_soai_trajectory(tmp_path):
    outdir = tmp_path / "soai_traj"
    main(["soai", "--seed-ee", "0.05", "--t-final", "500", "--outdir", str(outdir)])

    assert os.path.isfile(outdir / "soai_kinetic_trajectory.csv")


def test_cli_vthoff(tmp_path):
    csv_file = tmp_path / "multi_t.csv"
    csv_file.write_text(
        "temperature,ee_cat,ee_prod\n"
        "273.15,0.1,0.25\n"
        "273.15,0.5,0.75\n"
        "298.15,0.1,0.20\n"
        "298.15,0.5,0.70\n"
    )
    outdir = tmp_path / "vthoff"
    main(["vthoff", "--data", str(csv_file), "--outdir", str(outdir)])

    assert os.path.isfile(outdir / "vthoff_summary.json")


def test_cli_bootstrap(tmp_path):
    csv_file = tmp_path / "boot_data.csv"
    csv_file.write_text("ee_cat,ee_prod\n0.1,0.2\n0.3,0.5\n0.5,0.7\n0.7,0.85\n0.9,0.98\n")
    outdir = tmp_path / "boot"
    main(["bootstrap", "--data", str(csv_file), "--n-bootstraps", "20", "--outdir", str(outdir)])

    assert os.path.isfile(outdir / "bootstrap_summary.json")


def test_cli_gallery(tmp_path):
    outdir = tmp_path / "inputs"
    main(["gallery", "--ligand", "binol", "--software", "orca", "--outdir", str(outdir)])

    assert os.path.isfile(outdir / "binol_RR.inp")
    assert os.path.isfile(outdir / "binol_RS.inp")




