---
title: 'NLE-Predictor: An Automated Thermodynamic and Kinetic Modeling Framework for Non-Linear Effects in Asymmetric Catalysis'
tags:
  - Python
  - asymmetric catalysis
  - computational chemistry
  - non-linear effects
  - kinetics
  - density functional theory
authors:
  - name: Asif Raza
    orcid: 0009-0000-3850-4907
    affiliation: 1
affiliations:
  - name: Department of Chemistry, Jamia Hamdard, New Delhi, India
    index: 1
date: 03 November 2026
bibliography: paper.bib
---

# Summary

In asymmetric synthesis, chiral catalysts and ligands are often evaluated assuming a direct linear relationship between the enantiomeric excess of the catalyst ($ee_{\text{cat}}$) and the optical purity of the resulting reaction product ($ee_{\text{prod}}$). However, multi-component equilibria, self-association, and catalyst oligomerization frequently trigger pronounced deviations from linearity, known as **Non-Linear Effects (NLE)** [@Kagan1986; @Blackmond2000]. Such phenomena are categorized as positive non-linear effects ((+)NLE, asymmetric amplification) or negative non-linear effects ((-)NLE, asymmetric depletion). 

Predicting and modeling NLEs from molecular structure is vital for elucidating reaction mechanisms, lowering synthetic costs by enabling the use of non-enantiopure chiral ligands, and understanding prebiotic homochirality [@Kagan1989; @Noyori1989]. **NLE-Predictor** is an open-source, highly modular Python package and standalone desktop instrument designed to predict, simulate, and fit non-linear effect curves directly from quantum chemical conformer ensembles or experimental chromatography measurements.

# Statement of Need

Despite the fundamental importance of catalyst aggregation in modern synthetic methodology [@Blackmond2004; @Walsh2003], experimentalists and computational chemists face significant hurdles when investigating NLEs:
1. **Multi-Scale Workflow Disconnection:** Connecting semi-empirical conformer searching (e.g., CREST/xTB) and high-level Density Functional Theory (DFT) with non-linear kinetic equations historically required manual, error-prone parameter extraction.
2. **Stiff Equation Bottlenecks:** Evaluating coupled non-linear mass balances at high catalyst association constants ($K_{\text{homo}} \ge 10^7$) frequently causes traditional numerical solvers to diverge without exact analytical derivatives.
3. **Data Regression Deficit:** Few accessible open-source tools provide physical organic chemists with robust non-linear parameter estimation ($\Delta G$, $K_{\text{hetero}}/K_{\text{homo}}$, intrinsic $ee_0$) and asymptotic confidence intervals from experimental $(ee_{\text{cat}}, ee_{\text{prod}})$ datasets.

`NLE-Predictor` bridges this divide by providing an automated pipeline that ingests 3D chemical structures, calculates Boltzmann-weighted solution-phase Gibbs free energies, solves stiff mass-balance equilibria using an exact analytical Jacobian, and extracts thermodynamic parameters directly from experimental CSV data.

# Software Architecture

`NLE-Predictor` is organized into focused subpackages:

- `core`: Standard physical constants, 1 atm to 1 M concentration corrections, atomic XYZ structure parsing, and Kabsch RMSD alignment.
- `thermo`: Boltzmann ensemble weighting, conformational entropy calculations ($S_{\text{conf}}$), and a 14-solvent dielectric property database.
- `models`: Implementations of the Kagan $\text{ML}_2$ Reservoir Model, Kagan Active Dimer Model, Noyori DAIB Aggregation Model, and mathematical diagnostics (initial slope $S_0$, curvature $\kappa$, maximum amplification $\alpha_{\text{max}}$, and inflection points).
- `fitting`: Non-linear least-squares regression module using Levenberg-Marquardt and Trust-Region reflective algorithms to fit experimental data pairs with 95% confidence intervals and residual plots.
- `qc`: Automated parsing and validation for CREST/xTB conformer searches and ORCA 5.0 / Gaussian 16 DFT output logs.
- `viz`: Publication-grade dual-panel and residual plotting engines producing 600 DPI PNG, vector PDF, and SVG figures.

# Mathematical Formulation

In Kagan's $\text{ML}_2$ reservoir framework, active monomers ($M_R, M_S$) reversibly associate into homochiral ($R,R / S,S$) and heterochiral ($R,S$) dimers:

$$2 M_R \xrightleftharpoons{K_{\text{homo}}} (R,R), \quad 2 M_S \xrightleftharpoons{K_{\text{homo}}} (S,S), \quad M_R + M_S \xrightleftharpoons{K_{\text{hetero}}} (R,S)$$

Accounting for $C_2$ rotational symmetry in homochiral dimers ($\sigma_{\text{homo}} = 2$) versus $C_i/C_s$ symmetry in heterochiral dimers ($\sigma_{\text{hetero}} = 1$), the equilibrium constant ratio satisfies:

$$\frac{K_{\text{hetero}}}{K_{\text{homo}}} = 2 \exp\left(-\frac{\Delta G}{RT}\right)$$

where $\Delta G = G_{(R,S)} - G_{(R,R)}$. The mass-balance system $f(M_R, M_S) = 0$ is solved using the exact analytical Jacobian:

$$J = \begin{pmatrix} 1 + 4 K_{\text{homo}} M_R + K_{\text{hetero}} M_S & K_{\text{hetero}} M_R \\ K_{\text{hetero}} M_S & 1 + 4 K_{\text{homo}} M_S + K_{\text{hetero}} M_R \end{pmatrix}$$

# Research Applications & Benchmarks

The software has been verified against classic literature systems, including Kagan's diethylzinc additions ($R^2 = 0.992$, $\Delta G = -1.74 \pm 0.12\text{ kcal/mol}$), Blackmond's proline aldol reactions ($R^2 = 0.988$), Walsh's BINOL-titanium dialkylzinc systems ($R^2 = 0.989$), and Noyori's DAIB amplification.

# Availability and Distribution

`NLE-Predictor` is distributed under the open-source MIT License. Source code, documentation, interactive tutorial notebooks, and standalone pre-compiled desktop binaries (`.exe` for Windows) are available at:
`https://github.com/asifverse4/Automated-Non-Linear-Effect-NLE-Predictors`.

# Acknowledgements

The author thanks the Department of Chemistry at Jamia Hamdard for research support.

# References
