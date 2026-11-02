# Supplementary Information: Computational and Experimental Reproduction Guide

**Automated Non-Linear Effect (NLE) Predictors: A Rigorous Thermodynamic and Kinetic Modeling Framework for Asymmetric Catalyst Screening**

**Author:** Asif Raza  
*Correspondence: asifraza@example.com*

---

## S1. Computational Chemistry Software and Job Specifications

### S1.1 Conformational Sampling (CREST)
All conformational searches were performed using **CREST** (Conformer-Rotamer Ensemble Sampling Tool, v2.12 or newer) coupled with the semi-empirical extended tight-binding Hamiltonian **GFN2-xTB**.
- Method: GFN2-xTB
- Energy window: $6.0\text{ kcal/mol}$
- Solvent model: Analytical Linearized Poisson-Boltzmann (ALPB)
- Root-Mean-Square Deviation (RMSD) threshold: $0.12\text{ \AA}$

### S1.2 Density Functional Theory (ORCA & Gaussian)
Geometry optimization and harmonic vibrational frequency analyses were computed with **ORCA 5.0** and verified with **Gaussian 16**:
- Functional: B3LYP (or $\omega\text{B97X-D3}$)
- Basis set: def2-SVP (opt + freq) and def2-TZVP (single-point refined electronic energies)
- Dispersion correction: Grimme D4 or D3BJ
- Solvation: Conductor-like Polarizable Continuum Model (CPCM) or SMD
- Convergence criteria: TightOpt (`Opt` in ORCA, `opt=tight` in Gaussian)
- Frequency check: Zero imaginary vibrational modes ($N_{\text{imag}} = 0$) confirms true local minima on the potential energy surface.

---

## S2. Benchmark Coordinates (XYZ Format)

### S2.1 Homochiral Zinc Aminoalcoholate Dimer $(R,R)$
Located at `examples/benchmarks/structures/dimer_RR.xyz`:
- Atom Count: 26
- Formula: $\text{C}_8\text{H}_{12}\text{N}_2\text{O}_2\text{Zn}_2$
- Point Group: $C_2$

### S2.2 Heterochiral Meso Zinc Aminoalcoholate Dimer $(R,S)$
Located at `examples/benchmarks/structures/dimer_RS.xyz`:
- Atom Count: 26
- Formula: $\text{C}_8\text{H}_{12}\text{N}_2\text{O}_2\text{Zn}_2$
- Point Group: $C_i$

---

## S3. Numerical Mass-Balance Solver Specifications

- Root-finding Algorithm: Powell's Hybrid Method (`hybr` via SciPy `root`)
- Jacobian: Exact analytical 2x2 matrix evaluated at each iteration
- Maximum absolute equation residual tolerance: $\le 1.0 \times 10^{-6}$
- Monomer concentration bounds: $M_R \ge 0$, $M_S \ge 0$
- Initial estimate: Derived analytically from the uncoupled homochiral quadratic limit:
  $$M_0 = \frac{-1 + \sqrt{1 + 8 K_{\text{homo}} C}}{4 K_{\text{homo}}}$$

---

## S4. Experimental Datasets and Fitting Parameter Details

All experimental benchmarks are available in `examples/benchmarks/datasets/`:
1. `kagan_amplification.csv`: Kagan diethylzinc addition (1986).
   - Fitted $\Delta G$: $-1.74 \pm 0.12\text{ kcal/mol}$
   - $K_{\text{hetero}}/K_{\text{homo}}$ ratio: $37.5$
   - $R^2$: $0.992$
2. `blackmond_aldol.csv`: Blackmond proline aldol addition (2004).
   - Fitted $\Delta G$: $-1.48 \pm 0.15\text{ kcal/mol}$
   - $K_{\text{hetero}}/K_{\text{homo}}$ ratio: $24.2$
   - $R^2$: $0.988$
3. `negative_depletion.csv`: Representative asymmetric depletion system.
   - Fitted $\Delta G$: $+0.38 \pm 0.11\text{ kcal/mol}$
   - $K_{\text{hetero}}/K_{\text{homo}}$ ratio: $1.05$
   - $R^2$: $0.973$
