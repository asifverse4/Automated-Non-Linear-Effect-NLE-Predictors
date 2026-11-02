# NLE-Predictor: Comprehensive API Reference

This document provides a complete technical and chemical specification for all modules, classes, functions, and data structures in the `nle_predictor` library.

---

## Table of Contents
1. [Core Module (`nle_predictor.core`)](#1-core-module-nle_predictorcore)
2. [Thermodynamics & Solvation (`nle_predictor.thermo`)](#2-thermodynamics--solvation-nle_predictorthermo)
3. [Kinetic Models & Diagnostics (`nle_predictor.models`)](#3-kinetic-models--diagnostics-nle_predictormodels)
4. [Experimental Fitting (`nle_predictor.fitting`)](#4-experimental-fitting-nle_predictorfitting)
5. [Quantum Chemistry Interface (`nle_predictor.qc`)](#5-quantum-chemistry-interface-nle_predictorqc)
6. [Publication Visualization (`nle_predictor.viz`)](#6-publication-visualization-nle_predictorviz)
7. [Command-Line Interface (`nle_predictor.cli`)](#7-command-line-interface-nle_predictorcli)

---

## 1. Core Module (`nle_predictor.core`)

### 1.1 Constants (`nle_predictor.core.constants`)

- `HARTREE_TO_KCAL = 627.5094740631`: Thermochemical conversion factor from Hartrees to kcal/mol (CODATA 2018).
- `HARTREE_TO_KJ = 2625.4996394799`: Conversion factor from Hartrees to kJ/mol.
- `KCAL_TO_KJ = 4.184`: Exact thermochemical calorie-to-joule conversion.
- `GAS_CONSTANT_R = 0.00198720425864083`: Molar gas constant in $\text{kcal}/(\text{mol}\cdot\text{K})$.
- `GAS_CONSTANT_R_SI = 8.31446261815324`: Molar gas constant in $\text{J}/(\text{mol}\cdot\text{K})$.
- `DEFAULT_TEMP = 298.15`: Standard reference temperature in Kelvin ($25^\circ\text{C}$).
- `IDEAL_GAS_CONSTANT_L_ATM = 0.082057366081`: Ideal gas constant in $\text{L}\cdot\text{atm}/(\text{mol}\cdot\text{K})$.

### 1.2 Data Structures (`nle_predictor.core.types`)

#### `XYZStructure`
Container for 3D molecular Cartesian coordinates.
- Attributes:
  - `atoms`: `List[str]` (Atomic symbols, e.g., `["C", "H", "Zn"]`)
  - `coordinates`: `np.ndarray` of shape `(N, 3)` (Cartesian coordinates in \AA)
  - `comment`: `str`
  - `atom_count`: `int`
  - `formula`: `str` (Hill-notation molecular formula, e.g., `C8H12N2O2Zn2`)
  - `molecular_weight`: `float` in g/mol
- Properties:
  - `centroid`: `np.ndarray` of shape `(3,)` (Geometric centroid)

#### `NLESimulationResult`
Output container for simulated NLE trajectories.
- Attributes:
  - `model_name`: `str`
  - `ee_catalyst`: `np.ndarray` (Values in $[0, 1]$)
  - `ee_product`: `np.ndarray` (Values in $[0, 1]$)
  - `monomers_r`: `np.ndarray` (Equilibrium $[M_R]$ concentrations)
  - `monomers_s`: `np.ndarray` (Equilibrium $[M_S]$ concentrations)
  - `classification`: `str` (`"positive NLE"`, `"negative NLE"`, or `"near-linear"`)
  - `product_at_50`: `float` ($ee_{\text{prod}}$ at $50\%$ $ee_{\text{cat}}$)
  - `max_amplification`: `float` ($\max(ee_{\text{prod}}/ee_{\text{cat}})$)
  - `max_deviation_from_linearity`: `float`
  - `metadata`: `Dict[str, Any]`

#### `FittingResult`
Output container for non-linear regression on experimental data.
- Attributes:
  - `delta_g_kcal`: `float` (Fitted $\Delta G (RS - RR)$ in kcal/mol)
  - `delta_g_stderr`: `float` (Standard error of $\Delta G$)
  - `k_ratio`: `float` ($K_{\text{hetero}} / K_{\text{homo}}$)
  - `k_ratio_stderr`: `float` (Propagated standard error)
  - `r_squared`: `float` (Coefficient of determination $R^2$)
  - `rmse`: `float` (Root-mean-square error)
  - `residuals`: `np.ndarray`
  - `converged`: `bool`

### 1.3 Molecular Geometry (`nle_predictor.core.xyz`)

#### `validate_xyz_file(xyz_path: str) -> int`
Validates file existence, atom count headers, and numeric 3D Cartesian coordinates. Returns atom count.

#### `read_xyz(xyz_path: str) -> XYZStructure`
Parses an XYZ coordinate file into an `XYZStructure` object with automated molecular weight and Hill formula derivation.

#### `write_xyz(xyz_path: str, atoms: List[str], coordinates: np.ndarray, comment: str = "") -> None`
Writes Cartesian coordinates to an XYZ file with standard 8-decimal scientific precision.

#### `calculate_rmsd(coords_p: np.ndarray, coords_q: np.ndarray) -> float`
Computes the minimal Root-Mean-Square Deviation (RMSD) between two equal-sized coordinate sets using the **Kabsch algorithm** (centering, covariance SVD, reflection handling, and optimal 3D rotation).

---

## 2. Thermodynamics & Solvation (`nle_predictor.thermo`)

### 2.1 Boltzmann Weighting (`nle_predictor.thermo.boltzmann`)

#### `boltzmann_weights(energies_hartree: Sequence[float], temperature: float) -> np.ndarray`
Computes fractional Boltzmann populations:
$$p_i = \frac{\exp(-\Delta G_i / RT)}{\sum_j \exp(-\Delta G_j / RT)}$$

#### `conformational_entropy(weights: np.ndarray) -> float`
Calculates conformational mixing entropy in $\text{cal}/(\text{mol}\cdot\text{K})$:
$$S_{\text{conf}} = -R \sum_{i} p_i \ln(p_i)$$

#### `boltzmann_weighted_free_energy(energies_hartree: Sequence[float], temperature: float) -> float`
Calculates the ensemble Gibbs free energy in Hartrees:
$$G_{\text{ensemble}} = G_{\min} - RT \ln\left(\sum_i \exp\left(-\frac{G_i - G_{\min}}{RT}\right)\right)$$

### 2.2 Solvation & Standard-State Shifts (`nle_predictor.thermo.solvation`)

#### `concentration_correction_1atm_to_1m(temp_k: float = 298.15, delta_n: int = -1) -> float`
Computes the standard-state volume compression free energy shift ($1\text{ atm} \to 1\text{ M}$):
$$\Delta G^\circ_{\text{corr}} = -\Delta n \cdot RT \ln(V_m) \approx +1.894\text{ kcal/mol at } 298.15\text{ K}$$

#### `get_solvent_properties(solvent_name: str) -> Optional[SolventProperties]`
Queries dielectric permittivity ($\epsilon$), refractive index, and CPCM software aliases for 14+ common solvents (toluene, DCM, THF, chloroform, hexane, water, etc.).

---

## 3. Kinetic Models & Diagnostics (`nle_predictor.models`)

### 3.1 `NLECalculator` (Kagan $\text{ML}_2$ Reservoir Model)
Implements Kagan's foundational monomer-dimer reservoir model.
- `__init__(temp: float = 298.15)`
- `solve_mass_balance(ee_cat: float, k_homo: float, k_hetero: float) -> Tuple[float, float]`
  Solves the coupled non-linear mass balance using an exact analytical Jacobian.
- `simulate_kagan_model(ee_cat_range: Sequence[float], delta_g_kcal: float, k_homo_abs: float = 1e5, ee_0: float = 1.0) -> np.ndarray`
- `simulate(ee_cat_range: Sequence[float], delta_g_kcal: float, k_homo: float = 1e5, ee_0: float = 1.0) -> NLESimulationResult`

### 3.2 `KaganActiveDimerModel`
Extends the $\text{ML}_2$ model to catalytically active dimers.
- `__init__(temp: float = 298.15, g: float = 0.0, g_prime: float = 0.0, ee_0: float = 1.0, ee_d: float = 1.0)`
  - `g`: Relative turnover rate of homochiral dimer ($k_{\text{dimer}}/k_{\text{monomer}}$).
  - `g_prime`: Relative turnover rate of heterochiral meso dimer ($k_{\text{meso}}/k_{\text{monomer}}$).
  - `ee_0`: Intrinsic monomer enantiopurity.
  - `ee_d`: Intrinsic homochiral dimer enantiopurity.
- `simulate(ee_cat_range, delta_g_kcal, k_homo, ...) -> NLESimulationResult`

### 3.3 `NoyoriDAIBModel`
Implements the Noyori amino alcohol aggregation model for dialkylzinc additions.
- `__init__(temp: float = 298.15, ee_0: float = 0.98, kd_homo: float = 1e-3, kd_hetero: float = 1e-6)`
- `simulate(ee_cat_range, delta_g_kcal=None, ...) -> NLESimulationResult`

### 3.4 Mathematical Diagnostics (`nle_predictor.models.diagnostics`)

- `calculate_initial_slope(ee_cat, ee_prod) -> float`: Evaluates $S_0 = \left.\frac{d ee_{\text{prod}}}{d ee_{\text{cat}}}\right|_{0}$.
- `calculate_max_amplification(ee_cat, ee_prod) -> Tuple[float, float]`: Computes $\alpha_{\text{max}} = \max(ee_{\text{prod}}/ee_{\text{cat}})$ and catalyst $ee$ coordinate.
- `calculate_area_under_asymmetry(ee_cat, ee_prod) -> float`: Evaluates $\text{AUA} = \int_0^1 (ee_{\text{prod}} - x) dx$.
- `find_inflection_points(ee_cat, ee_prod, threshold: float = 1e-3) -> List[float]`: Identifies zero-crossings of $d^2(ee_{\text{prod}})/d(ee_{\text{cat}})^2$.
- `analyze_nle_diagnostics(ee_cat, ee_prod) -> NLEDiagnostics`: Complete diagnostic summary.

---

## 4. Experimental Fitting (`nle_predictor.fitting`)

#### `fit_kagan_reservoir(ee_cat_exp: np.ndarray, ee_prod_exp: np.ndarray, temp: float = 298.15, k_homo: float = 1e5, fit_ee_0: bool = False, initial_delta_g: float = 0.0) -> FittingResult`
Performs non-linear least-squares regression via Powell's dogleg/Levenberg-Marquardt optimization. Yields fitted $\Delta G$, $K_{\text{hetero}}/K_{\text{homo}}$, parameter standard errors, covariance matrices, and $R^2$.

#### `load_experimental_csv(csv_path: str) -> Tuple[np.ndarray, np.ndarray]`
Parses experimental CSV data files supporting percentage ($[0, 100]$) and fractional ($[0, 1.0]$) scales.

---

## 5. Quantum Chemistry Interface (`nle_predictor.qc`)

### 5.1 `DFTEvaluator` (ORCA Engine)
- `generate_input(xyz_file, run_name, functional="B3LYP", basis="def2-SVP", dispersion="D4", solvent=None, charge=0, multiplicity=1) -> str`
- `extract_free_energy(out_file: str) -> float`: Parses final Gibbs free energy in Hartrees after verifying normal termination and zero imaginary frequencies.
- `run_dft(inp_file: str, output_dir: str) -> str`

### 5.2 `ConformerSampler` (CREST Engine)
- `run_sampling(input_xyz: str, output_dir: str, solvent: Optional[str] = None, max_conformers: int = 1) -> Union[str, List[str]]`: Executes CREST conformer search and extracts lowest-energy conformers.

### 5.3 Gaussian Interface (`nle_predictor.qc.gaussian`)
- `validate_gaussian_output(log_file: str) -> None`: Validates normal termination and zero imaginary modes ($N_{\text{imag}} = 0$).
- `extract_gaussian_free_energy(log_file: str) -> float`: Parses "Sum of electronic and thermal Free Energies".

---

## 6. Publication Visualization (`nle_predictor.viz`)

- `plot_nle_with_experimental(sim_result, exp_cat=None, exp_prod=None, save_path=None, title=...) -> plt.Figure`: Generates JACS-styled curves with experimental scatter points and residual subplots.
- `plot_monomer_speciation(sim_result, save_path=None) -> plt.Figure`: Visualizes $[M_R]$, $[M_S]$, and $[M_R+M_S]$ equilibrium fractions.
- `plot_amplification_contour(calculator, delta_g_range=(-3.0, 3.0), n_points=40, save_path=None) -> plt.Figure`: 2D contour surface of $ee_{\text{prod}}$ across $(ee_{\text{cat}}, \Delta G)$.

---

## 7. Command-Line Interface (`nle_predictor.cli`)

Subcommands:
1. `nle-predictor run`: Conformer search (CREST) + DFT (ORCA) + NLE trajectory calculation.
2. `nle-predictor simulate`: Instantaneous analytical simulation given $\Delta G$, $K_{\text{homo}}$, and model type.
3. `nle-predictor fit`: Experimental data regression from CSV file.
4. `nle-predictor diagnostics`: Mathematical curve profiling.
