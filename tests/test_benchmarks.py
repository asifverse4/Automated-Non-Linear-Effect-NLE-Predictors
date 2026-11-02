"""Tests for benchmark chemical structures and literature datasets."""

import os

from nle_predictor.core.xyz import read_xyz, validate_xyz_file
from nle_predictor.fitting.experimental_fit import fit_kagan_reservoir, load_experimental_csv


def test_benchmark_structures():
    bench_dir = os.path.join("examples", "benchmarks", "structures")
    rr_file = os.path.join(bench_dir, "dimer_RR.xyz")
    rs_file = os.path.join(bench_dir, "dimer_RS.xyz")

    assert os.path.isfile(rr_file)
    assert os.path.isfile(rs_file)

    assert validate_xyz_file(rr_file) == 26
    assert validate_xyz_file(rs_file) == 26

    struct_rr = read_xyz(rr_file)
    struct_rs = read_xyz(rs_file)

    assert struct_rr.atom_count == 26
    assert struct_rs.atom_count == 26
    assert "Zn" in struct_rr.atoms
    assert struct_rr.formula == "C8H12N2O2Zn2"
    assert struct_rs.formula == "C8H12N2O2Zn2"


def test_kagan_literature_dataset_fit():
    csv_file = os.path.join("examples", "benchmarks", "datasets", "kagan_amplification.csv")
    assert os.path.isfile(csv_file)

    cat, prod = load_experimental_csv(csv_file)
    assert len(cat) == 13

    fit_res = fit_kagan_reservoir(cat, prod)
    assert fit_res.converged
    assert fit_res.r_squared > 0.95
    # Strong positive NLE indicates Delta G (RS - RR) < 0
    assert fit_res.delta_g_kcal < -0.5


def test_negative_nle_literature_dataset_fit():
    csv_file = os.path.join("examples", "benchmarks", "datasets", "negative_depletion.csv")
    assert os.path.isfile(csv_file)

    cat, prod = load_experimental_csv(csv_file)
    fit_res = fit_kagan_reservoir(cat, prod)
    assert fit_res.converged
    # Asymmetric depletion indicates Delta G (RS - RR) > 0
    assert fit_res.delta_g_kcal > 0.2


def test_walsh_binol_titanium_dataset_fit():
    csv_file = os.path.join("examples", "benchmarks", "datasets", "walsh_binol_titanium.csv")
    assert os.path.isfile(csv_file)

    cat, prod = load_experimental_csv(csv_file)
    assert len(cat) == 8

    fit_res = fit_kagan_reservoir(cat, prod)
    assert fit_res.converged
    assert fit_res.r_squared > 0.95
    assert fit_res.delta_g_kcal < -0.5


def test_noyori_daib_experimental_dataset():
    csv_file = os.path.join("examples", "benchmarks", "datasets", "noyori_daib_experimental.csv")
    assert os.path.isfile(csv_file)

    cat, prod = load_experimental_csv(csv_file)
    assert len(cat) == 8
    # 8% catalyst ee produces >= 90% product ee
    assert prod[1] >= 0.90

