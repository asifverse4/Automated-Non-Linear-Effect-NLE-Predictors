"""
Automated Publication Figure Generation Pipeline for NLE-Predictor Manuscript.
Generates all manuscript figures in 600 DPI PNG, vector PDF, and SVG formats.
Complies with JACS / Nature Computational Science publication standards.
"""

import os
import sys

# Ensure src is on path when running as standalone script
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

import matplotlib

matplotlib.use("Agg")  # Headless backend for CI/scripts
import matplotlib.pyplot as plt
import numpy as np

from nle_predictor.core.constants import DEFAULT_TEMP
from nle_predictor.fitting.experimental_fit import fit_kagan_reservoir, load_experimental_csv
from nle_predictor.models.kagan_active import KaganActiveDimerModel
from nle_predictor.models.kagan_ml2 import NLECalculator
from nle_predictor.models.noyori_daib import NoyoriDAIBModel
from nle_predictor.viz.plotting import (
    plot_amplification_contour,
    setup_publication_style,
)


def ensure_output_dir(dir_path: str = "manuscript/figures") -> str:
    """Ensure output figure directory exists."""
    os.makedirs(dir_path, exist_ok=True)
    return dir_path


def generate_figure_1_kagan_models(fig_dir: str) -> None:
    """Figure 1: Kagan ML2 Reservoir and Active Dimer NLE Trajectories."""
    setup_publication_style()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5.2))

    cat_ee = np.linspace(0.0, 1.0, 101)
    calc = NLECalculator(temp=DEFAULT_TEMP)

    # Subplot A: Reservoir model with varying Delta G
    ax1.plot(cat_ee * 100, cat_ee * 100, color="gray", linestyle="--", linewidth=1.5, label="Linear ($ee_p = ee_c$)")
    dg_values = [
        (-3.0, "#084594", r"$\Delta G = -3.0$ kcal/mol (Strong (+)NLE)"),
        (-1.5, "#2171b5", r"$\Delta G = -1.5$ kcal/mol ((+)NLE)"),
        (-0.41, "#4292c6", r"$\Delta G = -0.41$ kcal/mol (Statistical Linear)"),
        (+1.5, "#cb181d", r"$\Delta G = +1.5$ kcal/mol ((-)NLE Depletion)"),
    ]
    for dg, color, label in dg_values:
        prod = calc.simulate_kagan_model(cat_ee, delta_g_kcal=dg, k_homo_abs=1e4)
        ax1.plot(cat_ee * 100, prod * 100, color=color, linewidth=2.2, label=label)

    ax1.set_title("(a) Kagan ML$_2$ Reservoir Model", fontsize=13, fontweight="bold", pad=10)
    ax1.set_xlabel("Catalyst $ee$ (%)", fontsize=12, fontweight="bold")
    ax1.set_ylabel("Product $ee$ (%)", fontsize=12, fontweight="bold")
    ax1.set_xlim(0, 100)
    ax1.set_ylim(0, 100)
    ax1.grid(True, linestyle=":", alpha=0.5)
    ax1.legend(loc="lower right", fontsize=9.5, framealpha=0.95)

    # Subplot B: Active dimer model with varying kinetic activity ratio g
    ax2.plot(cat_ee * 100, cat_ee * 100, color="gray", linestyle="--", linewidth=1.5, label="Linear ($ee_p = ee_c$)")
    g_values = [
        (0.0, "#084594", "$g = 0.0$ (Pure Reservoir)"),
        (0.2, "#238b45", "$g = 0.2$ ($k_D/k_M = 0.2$)"),
        (0.8, "#d95f0e", "$g = 0.8$ ($k_D/k_M = 0.8$)"),
        (2.0, "#7a0177", "$g = 2.0$ (Dimer-Dominant)"),
    ]
    for g_val, color, label in g_values:
        active_model = KaganActiveDimerModel(g=g_val, g_prime=0.0, ee_0=1.0, ee_d=1.0)
        res = active_model.simulate(cat_ee, delta_g_kcal=-1.5, k_homo=1e4)
        ax2.plot(cat_ee * 100, res.ee_product * 100, color=color, linewidth=2.2, label=label)

    ax2.set_title(r"(b) Kagan Active Dimer Model ($\Delta G = -1.5$ kcal/mol)", fontsize=13, fontweight="bold", pad=10)
    ax2.set_xlabel("Catalyst $ee$ (%)", fontsize=12, fontweight="bold")
    ax2.set_ylabel("Product $ee$ (%)", fontsize=12, fontweight="bold")
    ax2.set_xlim(0, 100)
    ax2.set_ylim(0, 100)
    ax2.grid(True, linestyle=":", alpha=0.5)
    ax2.legend(loc="lower right", fontsize=9.5, framealpha=0.95)

    plt.tight_layout()
    _save_all_formats(fig, os.path.join(fig_dir, "figure_1_kagan_models"))
    plt.close(fig)


def generate_figure_2_speciation(fig_dir: str) -> None:
    """Figure 2: Catalyst Equilibrium Speciation and Trapping."""
    setup_publication_style()
    calc = NLECalculator()
    cat_ee = np.linspace(0.0, 1.0, 101)
    k_homo = 1e4
    delta_g = -1.8  # Strong (+)NLE
    ratio = 2.0 * np.exp(-delta_g / (0.001987204 * DEFAULT_TEMP))
    k_hetero = k_homo * ratio

    m_r = []
    m_s = []
    d_rr = []
    d_ss = []
    d_rs = []

    for c in cat_ee:
        mr, ms = calc.solve_mass_balance(c, k_homo, k_hetero)
        m_r.append(mr)
        m_s.append(ms)
        d_rr.append(k_homo * (mr ** 2))
        d_ss.append(k_homo * (ms ** 2))
        d_rs.append(k_hetero * mr * ms)

    fig, ax = plt.subplots(figsize=(7, 5.5))
    x_pct = cat_ee * 100

    ax.plot(x_pct, m_r, color="#1b9e77", linewidth=2.5, label="Major Active Monomer $[M_R]$")
    ax.plot(x_pct, m_s, color="#d95f02", linestyle="-.", linewidth=2.0, label="Minor Active Monomer $[M_S]$")
    ax.plot(x_pct, d_rs, color="#7570b3", linewidth=2.5, label="Heterochiral Dimer $[(R,S)]$ (Meso Trap)")
    ax.plot(x_pct, d_rr, color="#e7298a", linestyle=":", linewidth=2.0, label="Homochiral Dimer $[(R,R)]$")

    ax.set_title(r"Equilibrium Speciation Profile ($\Delta G = -1.8$ kcal/mol)", fontsize=13, fontweight="bold", pad=12)
    ax.set_xlabel("Catalyst $ee$ (%)", fontsize=12, fontweight="bold")
    ax.set_ylabel("Species Concentration Fraction", fontsize=12, fontweight="bold")
    ax.set_xlim(0, 100)
    ax.set_ylim(bottom=0.0)
    ax.grid(True, linestyle=":", alpha=0.5)
    ax.legend(loc="center right", fontsize=10.5, framealpha=0.95)

    plt.tight_layout()
    _save_all_formats(fig, os.path.join(fig_dir, "figure_2_equilibrium_speciation"))
    plt.close(fig)


def generate_figure_3_experimental_benchmarks(fig_dir: str) -> None:
    """Figure 3: Experimental Data Fitting on Kagan and Blackmond Datasets."""
    setup_publication_style()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.5))

    cat_dense = np.linspace(0.0, 1.0, 101)
    calc = NLECalculator()

    # Kagan Diethylzinc Fit
    kagan_csv = os.path.join("examples", "benchmarks", "datasets", "kagan_amplification.csv")
    c_kagan, p_kagan = load_experimental_csv(kagan_csv)
    res_kagan = fit_kagan_reservoir(c_kagan, p_kagan)
    sim_kagan = calc.simulate_kagan_model(cat_dense, delta_g_kcal=res_kagan.delta_g_kcal)

    ax1.plot(cat_dense * 100, cat_dense * 100, color="gray", linestyle="--", label="Linear Line")
    ax1.plot(cat_dense * 100, sim_kagan * 100, color="#08519c", linewidth=2.5,
             label=rf"Fitted Kagan ML$_2$ Model\n($\Delta G = {res_kagan.delta_g_kcal:.2f}$ kcal/mol, $R^2 = {res_kagan.r_squared:.3f}$)")
    ax1.scatter(c_kagan * 100, p_kagan * 100, color="#e6550d", edgecolors="black", s=60, zorder=5, label="Kagan (1986) Data")
    ax1.set_title("(a) Kagan Diethylzinc Addition", fontsize=12, fontweight="bold", pad=10)
    ax1.set_xlabel("Catalyst $ee$ (%)", fontsize=12, fontweight="bold")
    ax1.set_ylabel("Product $ee$ (%)", fontsize=12, fontweight="bold")
    ax1.set_xlim(0, 100)
    ax1.set_ylim(0, 100)
    ax1.grid(True, linestyle=":", alpha=0.5)
    ax1.legend(loc="lower right", fontsize=9.5)

    # Blackmond Aldol Fit
    blackmond_csv = os.path.join("examples", "benchmarks", "datasets", "blackmond_aldol.csv")
    c_bm, p_bm = load_experimental_csv(blackmond_csv)
    res_bm = fit_kagan_reservoir(c_bm, p_bm)
    sim_bm = calc.simulate_kagan_model(cat_dense, delta_g_kcal=res_bm.delta_g_kcal)

    ax2.plot(cat_dense * 100, cat_dense * 100, color="gray", linestyle="--", label="Linear Line")
    ax2.plot(cat_dense * 100, sim_bm * 100, color="#08519c", linewidth=2.5,
             label=rf"Fitted Kagan ML$_2$ Model\n($\Delta G = {res_bm.delta_g_kcal:.2f}$ kcal/mol, $R^2 = {res_bm.r_squared:.3f}$)")
    ax2.scatter(c_bm * 100, p_bm * 100, color="#31a354", edgecolors="black", s=60, zorder=5, label="Blackmond (2004) Data")
    ax2.set_title("(b) Blackmond Proline-Catalyzed Aldol", fontsize=12, fontweight="bold", pad=10)
    ax2.set_xlabel("Catalyst $ee$ (%)", fontsize=12, fontweight="bold")
    ax2.set_ylabel("Product $ee$ (%)", fontsize=12, fontweight="bold")
    ax2.set_xlim(0, 100)
    ax2.set_ylim(0, 100)
    ax2.grid(True, linestyle=":", alpha=0.5)
    ax2.legend(loc="lower right", fontsize=9.5)

    plt.tight_layout()
    _save_all_formats(fig, os.path.join(fig_dir, "figure_3_experimental_benchmarks"))
    plt.close(fig)


def generate_figure_4_contour_landscape(fig_dir: str) -> None:
    """Figure 4: 2D Amplification Contour Map."""
    calc = NLECalculator()
    fig = plot_amplification_contour(calc, delta_g_range=(-3.0, 3.0), n_points=50)
    _save_all_formats(fig, os.path.join(fig_dir, "figure_4_amplification_contour"))
    plt.close(fig)


def generate_figure_5_noyori_daib(fig_dir: str) -> None:
    """Figure 5: Noyori DAIB Aggregation Model."""
    setup_publication_style()
    fig, ax = plt.subplots(figsize=(6.5, 5.5))

    cat_ee = np.linspace(0.0, 1.0, 101)
    ax.plot(cat_ee * 100, cat_ee * 100, color="gray", linestyle="--", label="Ideal Linear Line")

    # Classic Noyori DAIB parameters
    noyori = NoyoriDAIBModel(ee_0=0.98, kd_homo=1e-3, kd_hetero=1e-6)
    res_noyori = noyori.simulate(cat_ee)

    ax.plot(cat_ee * 100, res_noyori.ee_product * 100, color="#6a51a3", linewidth=2.8,
            label="Noyori DAIB Model\n($K_d(hetero) = 10^{-6}$, $K_d(homo) = 10^{-3}$)")

    # Experimental landmark points from Noyori (1989)
    noyori_exp_cat = np.array([0.08, 0.15, 0.30, 0.50, 0.70, 0.90])
    noyori_exp_prod = np.array([0.92, 0.95, 0.97, 0.98, 0.98, 0.98])
    ax.scatter(noyori_exp_cat * 100, noyori_exp_prod * 100, color="#ff7f00", edgecolors="black", s=65, zorder=5, label="Noyori (1989) Data Points")

    ax.set_title("Noyori DAIB Extreme Asymmetric Amplification", fontsize=13, fontweight="bold", pad=12)
    ax.set_xlabel("Catalyst $ee$ (%)", fontsize=12, fontweight="bold")
    ax.set_ylabel("Product $ee$ (%)", fontsize=12, fontweight="bold")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.grid(True, linestyle=":", alpha=0.5)
    ax.legend(loc="lower right", fontsize=10, framealpha=0.95)

    plt.tight_layout()
    _save_all_formats(fig, os.path.join(fig_dir, "figure_5_noyori_daib"))
    plt.close(fig)


def _save_all_formats(fig: plt.Figure, base_path: str) -> None:
    """Save figure in PNG (600 DPI), vector PDF, and SVG formats."""
    fig.savefig(f"{base_path}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{base_path}.pdf", bbox_inches="tight")
    fig.savefig(f"{base_path}.svg", bbox_inches="tight")
    print(f"Generated: {base_path} (.png, .pdf, .svg)")


def main():
    print("=" * 60)
    print(" NLE-Predictor: Publication Figure Generation Pipeline")
    print("=" * 60)
    fig_dir = ensure_output_dir("manuscript/figures")

    print("\nGenerating Figure 1: Kagan ML2 Reservoir & Active Dimer Models...")
    generate_figure_1_kagan_models(fig_dir)

    print("Generating Figure 2: Catalyst Speciation & Enantiomeric Trapping...")
    generate_figure_2_speciation(fig_dir)

    print("Generating Figure 3: Experimental Benchmark Regressions...")
    generate_figure_3_experimental_benchmarks(fig_dir)

    print("Generating Figure 4: 2D Amplification Contour Landscape...")
    generate_figure_4_contour_landscape(fig_dir)

    print("Generating Figure 5: Noyori DAIB Aggregation Model...")
    generate_figure_5_noyori_daib(fig_dir)

    print("\n" + "=" * 60)
    print(f"All manuscript figures successfully rendered to: {fig_dir}")
    print("=" * 60)


if __name__ == "__main__":
    main()
