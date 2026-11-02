"""Tests for cross-platform portability, module execution, and environment specifications."""

import os
import subprocess
import sys


def test_python_m_nle_predictor_module_invocation(tmp_path):
    outdir = str(tmp_path / "mod_sim")
    cmd = [
        sys.executable,
        "-m",
        "nle_predictor",
        "simulate",
        "--delta-g",
        "-1.5",
        "--outdir",
        outdir,
    ]
    env = os.environ.copy()
    env["PYTHONPATH"] = "src"
    res = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        env=env,
        stdin=subprocess.DEVNULL,
    )

    assert res.returncode == 0
    assert os.path.isfile(os.path.join(outdir, "simulated_nle.csv"))
    assert os.path.isfile(os.path.join(outdir, "simulated_curve.png"))


def test_root_script_execution_with_legacy_and_subcommand(tmp_path):
    # Test subcommand via root script
    outdir_sub = str(tmp_path / "root_sim")
    res_sub = subprocess.run(
        [sys.executable, "nle_predictor.py", "simulate", "--delta-g", "-1.8", "--outdir", outdir_sub],
        capture_output=True,
        text=True,
        stdin=subprocess.DEVNULL,
    )
    assert res_sub.returncode == 0
    assert os.path.isfile(os.path.join(outdir_sub, "simulated_nle.csv"))

    # Test legacy flags without subcommand via root script
    outdir_mock = str(tmp_path / "root_mock")
    res_mock = subprocess.run(
        [sys.executable, "nle_predictor.py", "--mock", "--outdir", outdir_mock],
        capture_output=True,
        text=True,
        stdin=subprocess.DEVNULL,
    )
    assert res_mock.returncode == 0
    assert os.path.isfile(os.path.join(outdir_mock, "nle_results.csv"))


def test_environment_specifications_exist():
    assert os.path.isfile("environment.yml")
    assert os.path.isfile("environment-dev.yml")
    assert os.path.isfile("Dockerfile")
    assert os.path.isfile("docker-compose.yml")
    assert os.path.isfile("Makefile")
    assert os.path.isfile("make.ps1")

    # Verify Dockerfile syntax contains entrypoint
    with open("Dockerfile", "r", encoding="utf-8") as f:
        docker_content = f.read()
    assert "ENTRYPOINT [\"nle-predictor\"]" in docker_content
    assert "MPLBACKEND=Agg" in docker_content

    # Verify environment.yml contains dependencies
    with open("environment.yml", "r", encoding="utf-8") as f:
        env_content = f.read()
    assert "name: nle-predictor" in env_content
    assert "conda-forge" in env_content
