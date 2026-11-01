"""
Publication-ready scientific plotting suite tailored for asymmetric catalysis literature.
Complies with JACS / Angewandte Chemie graphic standards.
"""

from typing import Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np

from nle_predictor.core.types import NLESimulationResult


def setup_publication_style():
    """Apply clean, modern typography and styling suitable for chemical journals."""
    plt.style.use("default")
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 11,
        "axes.labelsize": 12,
        "axes.labelweight": "bold",
        "axes.titlesize": 13,
        "axes.titleweight": "bold",
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
        "figure.titlesize": 14,
        "lines.linewidth": 2.0,
        "lines.markersize": 6,
        "axes.edgecolor": "#333333",
        "axes.linewidth": 1.2,
    })


def plot_nle_with_experimental(
    sim_result: NLESimulationResult,
    exp_cat: Optional[np.ndarray] = None,
    exp_prod: Optional[np.ndarray] = None,
    save_path: Optional[str] = None,
    title: str = "Asymmetric Amplification & Non-Linear Effect",
) -> plt.Figure:
    """
    Plot simulated NLE curve with optional experimental data scatter points
    and an optional residual panel if experimental data is supplied.
    """
    setup_publication_style()

    has_exp = exp_cat is not None and exp_prod is not None and len(exp_cat) > 0

    if has_exp:
        fig, (ax_main, ax_res) = plt.subplots(
            2, 1, figsize=(6.5, 7.5), gridspec_kw={"height_ratios": [3.5, 1]}, sharex=True
        )
    else:
        fig, ax_main = plt.subplots(figsize=(6, 5.5))
        ax_res = None

    cat_pct = sim_result.ee_catalyst * 100.0
    prod_pct = sim_result.ee_product * 100.0

    # Ideal linear reference line
    ax_main.plot(
        cat_pct, cat_pct,
        color="#888888", linestyle="--", linewidth=1.5,
        label="Ideal Linear Line (No Aggregation)",
    )

    # Theoretical simulation curve
    ax_main.plot(
        cat_pct, prod_pct,
        color="#004488", linewidth=2.5,
        label=f"{sim_result.model_name} ({sim_result.classification})",
    )

    # Experimental points
    if has_exp:
        exp_c_pct = exp_cat * 100.0
        exp_p_pct = exp_prod * 100.0
        ax_main.scatter(
            exp_c_pct, exp_p_pct,
            color="#D95F02", edgecolors="black", s=50, zorder=5,
            label="Experimental Data",
        )

        # Residuals in bottom panel
        pred_interp = np.interp(exp_cat, sim_result.ee_catalyst, sim_result.ee_product)
        residuals_pct = (exp_prod - pred_interp) * 100.0

        ax_res.axhline(0, color="gray", linestyle="--", linewidth=1.0)
        ax_res.scatter(exp_c_pct, residuals_pct, color="#D95F02", edgecolors="black", s=40)
        ax_res.set_ylabel("Residual (%)", fontsize=10, fontweight="bold")
        ax_res.set_xlabel("Catalyst $ee$ (%)", fontsize=12, fontweight="bold")
        ax_res.set_xlim(0, 100)
        ax_res.set_ylim(-15, 15)
        ax_res.grid(True, linestyle=":", alpha=0.5)

    ax_main.set_xlim(0, 100)
    ax_main.set_ylim(0, 100)
    if not has_exp:
        ax_main.set_xlabel("Catalyst $ee$ (%)", fontsize=12, fontweight="bold")
        ax_main.set_aspect("equal", adjustable="box")
    ax_main.set_ylabel("Product $ee$ (%)", fontsize=12, fontweight="bold")
    ax_main.set_title(title, pad=12)
    ax_main.grid(True, linestyle=":", alpha=0.5)
    ax_main.legend(loc="lower right", framealpha=0.95, edgecolor="gray")

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=600, bbox_inches="tight")
        plt.close(fig)
    return fig


def plot_monomer_speciation(
    sim_result: NLESimulationResult,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """
    Plot concentration speciation of monomer pools [M_R] and [M_S]
    as catalyst enantiomeric excess increases from 0% to 100%.
    """
    setup_publication_style()
    fig, ax = plt.subplots(figsize=(6.5, 5))

    cat_pct = sim_result.ee_catalyst * 100.0
    ax.plot(cat_pct, sim_result.monomers_r, color="#1B9E77", label="Active Monomer $[M_R]$")
    ax.plot(cat_pct, sim_result.monomers_s, color="#D95F02", linestyle="-.", label="Active Monomer $[M_S]$")
    ax.plot(
        cat_pct, sim_result.monomers_r + sim_result.monomers_s,
        color="#7570B3", linestyle=":", label="Total Active Pool $[M_R + M_S]$",
    )

    ax.set_xlim(0, 100)
    ax.set_ylim(bottom=0.0)
    ax.set_xlabel("Catalyst $ee$ (%)", fontsize=12, fontweight="bold")
    ax.set_ylabel("Monomer Fraction", fontsize=12, fontweight="bold")
    ax.set_title("Catalyst Equilibrium Speciation", pad=12)
    ax.grid(True, linestyle=":", alpha=0.5)
    ax.legend(loc="upper left", framealpha=0.95, edgecolor="gray")

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=600, bbox_inches="tight")
        plt.close(fig)
    return fig


def plot_amplification_contour(
    calculator,
    delta_g_range: Sequence[float] = (-3.0, 3.0),
    n_points: int = 40,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """
    Generate a 2D contour map of Product ee across a grid of (Catalyst ee, Delta G).
    """
    setup_publication_style()
    fig, ax = plt.subplots(figsize=(7, 5.5))

    cat_ee = np.linspace(0.0, 1.0, n_points)
    delta_gs = np.linspace(delta_g_range[0], delta_g_range[1], n_points)
    X, Y = np.meshgrid(cat_ee * 100.0, delta_gs)
    Z = np.zeros_like(X)

    for i, dg in enumerate(delta_gs):
        prod = calculator.simulate_kagan_model(cat_ee, delta_g_kcal=dg, k_homo_abs=1e4)
        Z[i, :] = prod * 100.0

    contour = ax.contourf(X, Y, Z, levels=20, cmap="viridis")
    cbar = fig.colorbar(contour, ax=ax)
    cbar.set_label("Product $ee$ (%)", fontweight="bold")

    # Overlay Delta G = 0 neutral threshold line
    ax.axhline(0, color="white", linestyle="--", linewidth=1.5, label=r"$\Delta G = 0$ (Linear limit)")

    ax.set_xlabel("Catalyst $ee$ (%)", fontweight="bold")
    ax.set_ylabel(r"$\Delta G (RS - RR)$ (kcal/mol)", fontweight="bold")
    ax.set_title(r"NLE Landscape: Product $ee$ vs $\Delta G$", pad=12)
    ax.legend(loc="upper right", framealpha=0.8)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=600, bbox_inches="tight")
        plt.close(fig)
    return fig
