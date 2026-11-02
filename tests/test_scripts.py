"""Tests for manuscript figure generation pipeline."""

import os

from scripts.generate_paper_figures import (
    ensure_output_dir,
    generate_figure_1_kagan_models,
    generate_figure_2_speciation,
    generate_figure_3_experimental_benchmarks,
    generate_figure_4_contour_landscape,
    generate_figure_5_noyori_daib,
)


def test_paper_figures_generation(tmp_path):
    fig_dir = str(tmp_path / "figures")
    ensure_output_dir(fig_dir)

    generate_figure_1_kagan_models(fig_dir)
    generate_figure_2_speciation(fig_dir)
    generate_figure_3_experimental_benchmarks(fig_dir)
    generate_figure_4_contour_landscape(fig_dir)
    generate_figure_5_noyori_daib(fig_dir)

    for fig_name in [
        "figure_1_kagan_models",
        "figure_2_equilibrium_speciation",
        "figure_3_experimental_benchmarks",
        "figure_4_amplification_contour",
        "figure_5_noyori_daib",
    ]:
        for ext in [".png", ".pdf", ".svg"]:
            path = os.path.join(fig_dir, f"{fig_name}{ext}")
            assert os.path.isfile(path)
            assert os.path.getsize(path) > 1000


def test_manuscript_sources_and_submission_bundle(tmp_path):
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    manuscript_dir = os.path.join(repo_root, "manuscript")

    # Verify key manuscript files exist
    for filename in [
        "MANUSCRIPT.md",
        "manuscript.tex",
        "manuscript_generic.tex",
        "references.bib",
        "SUPPLEMENTARY_INFORMATION.md",
        "supplementary_information.tex",
        "interactive_dashboard.html",
    ]:
        path = os.path.join(manuscript_dir, filename)
        assert os.path.isfile(path), f"Missing {filename}"
        assert os.path.getsize(path) > 500

    # Test submission bundle packaging into temporary zip
    import zipfile
    bundle_zip = str(tmp_path / "test_submission_bundle.zip")
    with zipfile.ZipFile(bundle_zip, "w", zipfile.ZIP_DEFLATED) as z:
        for root, _, files in os.walk(manuscript_dir):
            for file in files:
                if not file.endswith(".zip"):
                    full_p = os.path.join(root, file)
                    rel_p = os.path.relpath(full_p, manuscript_dir)
                    z.write(full_p, rel_p)

    assert os.path.isfile(bundle_zip)
    with zipfile.ZipFile(bundle_zip, "r") as z:
        names = z.namelist()
        assert "manuscript.tex" in names
        assert "manuscript_generic.tex" in names
        assert "interactive_dashboard.html" in names
        assert "references.bib" in names
        assert any("figure_1_kagan_models" in n for n in names)


