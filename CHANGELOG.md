# Changelog

All notable changes to the **NLE-Predictor** project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [1.0.0] - 2026-11-02

### Added
- **Unified Scientific Package Architecture:**
  - Modern PEP 517/621 modular package structure under `src/nle_predictor/`.
  - Subpackages: `core`, `thermo`, `models`, `fitting`, `qc`, and `viz`.
  - Console script entry point `nle-predictor = "nle_predictor.cli:main"` and universal module invocation `python -m nle_predictor`.
  - 100% backwards-compatible root `nle_predictor.py` facade re-exporting all public APIs and handling legacy flags.

- **Kinetic & Thermodynamic Modeling Engines:**
  - **Kagan $\text{ML}_2$ Reservoir Model:** Exact analytical Jacobian solver for non-linear mass-balance equilibria, guaranteeing quadratic convergence even at extreme association constants ($K \ge 10^7$).
  - **Kagan Active Dimer Model:** Flexible kinetic simulation supporting differential turnover rates ($g = k_D / k_M$, $g' = k_H / k_M$) and intrinsic stereoselectivities ($ee_0, ee_D$).
  - **Noyori DAIB Aggregation Model:** Non-linear model reproducing asymmetric amplification in dialkylzinc additions with active monomer dissociation equilibria.
  - **Mathematical Diagnostics Profiling:** Automated calculation of initial slope ($S_0$), maximum amplification ($\alpha_{\text{max}}$), Area Under Asymmetry (AUA), mean curvature ($\kappa$), and numerical inflection point detection.
  - **Non-Linear Least Squares Regression:** Robust experimental fitting engine using `scipy.optimize.curve_fit` extracting $\Delta G$, $K_{\text{hetero}}/K_{\text{homo}}$, and intrinsic $ee_0$ with standard errors, 95% confidence intervals, and $R^2$ statistics.

- **Quantum Chemical & Solvation Tooling:**
  - Standard-state concentration correction ($\Delta G^\circ = +1.894\text{ kcal/mol}$ for $2\text{ M} \to 1\text{ D}$ in 1 M solution).
  - Solvent database covering 14 common reaction solvents with dielectric constants and standard-state shifts.
  - Conformer ensemble Boltzmann weighting with conformational entropy ($S_{\text{conf}}$).
  - Robust parser and validator for ORCA and Gaussian 09/16 output logs.
  - Atomic structure reader/writer (`XYZStructure`) with Hill chemical formula calculation and Kabsch RMSD alignment algorithm.

- **Academic Publication & Manuscript Artifacts:**
  - Full scientific preprint in Markdown (`manuscript/MANUSCRIPT.md`) and Supplementary Information (`manuscript/SUPPLEMENTARY_INFORMATION.md`).
  - JACS/ACS template (`manuscript/manuscript.tex`) and generic multi-journal LaTeX template (`manuscript/manuscript_generic.tex`) with BibTeX references (`manuscript/references.bib`).
  - Automated high-resolution publication figure generation script (`scripts/generate_paper_figures.py`) outputting 600 DPI PNG, vector PDF, and SVG for Figures 1–5.
  - Interactive research tutorial Jupyter notebook (`notebooks/nle_research_tutorial.ipynb`).

- **Multi-Platform Portability & DevOps:**
  - Standalone Windows desktop executable (`nle-predictor.exe`) bundled into `nle-predictor-windows-x64.zip` with zero runtime Python dependencies.
  - Interactive desktop console mode with formatted guided menus for simulation, fitting, Noyori analysis, and figure generation.
  - Containerized reproducible environment (`Dockerfile`, `docker-compose.yml`, `.dockerignore`).
  - Multi-platform Conda environments (`environment.yml` and `environment-dev.yml`).
  - Cross-platform task automation runners (`Makefile` for Linux/macOS and `make.ps1` for Windows PowerShell).
  - Automated journal submission packaging target (`make paper-bundle` / `.\make.ps1 paper-bundle`).
  - GitHub Actions CI matrix testing across Ubuntu, macOS, and Windows with Python 3.10–3.13 (`.github/workflows/ci.yml`).
  - Automated GitHub release pipeline packaging wheels, sdist, and manuscript bundles (`.github/workflows/release.yml`).
  - Comprehensive 53-test verification suite covering 100% of submodules.

