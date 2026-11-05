"""
Unit tests for Chiral Ligand 3D Conformer Gallery and QC input generation.
"""

import os
import tempfile

import pytest

from nle_predictor.qc.gallery import ChiralLigandGallery


def test_gallery_available_ligands():
    """Verify registry of benchmark chiral ligands."""
    gallery = ChiralLigandGallery()
    ligands = gallery.list_ligands()
    assert "daib" in ligands
    assert "binol" in ligands
    assert "pybox" in ligands


def test_gallery_structure_loading():
    """Verify loading and geometric consistency of homochiral and heterochiral conformers."""
    gallery = ChiralLigandGallery()
    for lig in ["daib", "binol", "pybox"]:
        struct_rr = gallery.load_structure(lig, "RR")
        struct_rs = gallery.load_structure(lig, "RS")

        assert struct_rr.atom_count > 0
        assert struct_rs.atom_count > 0
        # Diastereomeric pairs should have identical stoichiometry/atom counts
        assert struct_rr.atom_count == struct_rs.atom_count
        assert struct_rr.coordinates.shape == (struct_rr.atom_count, 3)


def test_orca_input_generation():
    """Verify ORCA input format matches quantum chemical specification."""
    gallery = ChiralLigandGallery()
    text = gallery.generate_input_text(
        "binol",
        "RR",
        software="orca",
        functional="wB97X-D4",
        basis="def2-TZVP",
        dispersion="",
        solvent="dichloromethane",
        n_cores=16,
        memory_mb=6000,
    )
    assert "! wB97X-D4 def2-TZVP" in text
    assert "! CPCM(dichloromethane)" in text
    assert "%pal nprocs 16 end" in text
    assert "%maxcore 6000" in text
    assert "* xyz 0 1" in text
    assert "Ti" in text


def test_gaussian_input_generation():
    """Verify Gaussian 16 input syntax including route, title, and coordinates."""
    gallery = ChiralLigandGallery()
    text = gallery.generate_input_text(
        "pybox",
        "RS",
        software="gaussian",
        functional="B3LYP",
        basis="def2-SVP",
        dispersion="D3BJ",
        solvent="THF",
        n_cores=8,
    )
    assert "%nprocshared=8" in text
    assert "#p B3LYP/def2-SVP empiricaldispersion=gd3bj scrf=(smd,solvent=THF) opt freq" in text
    assert "0 1" in text
    assert "Cu" in text


def test_qchem_input_generation():
    """Verify Q-Chem input blocks."""
    gallery = ChiralLigandGallery()
    text = gallery.generate_input_text(
        "daib",
        "RR",
        software="qchem",
        functional="M06-2X",
        basis="def2-SV(P)",
        solvent="toluene",
    )
    assert "$rem" in text
    assert "jobtype opt" in text
    assert "method m06-2x" in text
    assert "$molecule" in text
    assert "$pcm" in text


def test_export_pair_writes_files():
    """Verify export_pair writes both files to disk."""
    gallery = ChiralLigandGallery()
    with tempfile.TemporaryDirectory() as tmpdir:
        rr_path, rs_path = gallery.export_pair("daib", output_dir=tmpdir, software="orca")
        assert os.path.isfile(rr_path)
        assert os.path.isfile(rs_path)
        assert os.path.getsize(rr_path) > 100
        assert os.path.getsize(rs_path) > 100


def test_invalid_parameters_raise():
    """Verify error checking on invalid ligand, chirality, or software."""
    gallery = ChiralLigandGallery()
    with pytest.raises(KeyError, match="not recognized"):
        gallery.get_metadata("unknown_ligand")

    with pytest.raises(ValueError, match="Chirality must be 'RR' or 'RS'"):
        gallery.get_structure_path("binol", "INVALID")

    with pytest.raises(ValueError, match="Unsupported software package"):
        gallery.generate_input_text("daib", "RR", software="turbomole")
