"""Tests for manuscript figure generation pipeline."""

import os
import pytest

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
