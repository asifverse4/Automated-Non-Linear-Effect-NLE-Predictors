"""Tests for the unified NLE-Predictor CLI interface in nle_predictor.cli."""

import os
import pytest

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
