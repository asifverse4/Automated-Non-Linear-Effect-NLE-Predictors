# NLE // ATLAS

### Thermodynamic intelligence for asymmetric catalyst design

<p align="center">
  <img src="assets/nle-workflow.svg" alt="Animated three-dimensional NLE Atlas computational workflow" width="100%">
</p>

<p align="center">
  <strong>A research-grade computational instrument for non-linear effects.</strong><br>
  Conformer ensembles. Free-energy surfaces. Enantioselectivity forecasts.
</p>

<p align="center">
  <img src="https://img.shields.io/badge/Python-3.10%2B-3776AB?style=flat-square&logo=python&logoColor=white" alt="Python 3.10 or newer">
  <img src="https://img.shields.io/badge/CREST-xTB-176B87?style=flat-square" alt="CREST and xTB">
  <img src="https://img.shields.io/badge/ORCA-DFT-5B4B8A?style=flat-square" alt="ORCA DFT">
  <img src="https://img.shields.io/badge/Model-Kagan%20ML2-C45A2A?style=flat-square" alt="Kagan ML2 model">
</p>

## Overview

NLE Predictor is a reproducible computational screening instrument for
asymmetric catalyst design. It compares homochiral (`R,R`) and heterochiral
(`R,S`) catalyst dimers, evaluates low-energy conformer ensembles, estimates
their thermodynamic stability, and propagates that difference through a
Kagan-style mass-balance model.

Each run preserves the source structures, conformer archive, calculation
inputs and outputs, ensemble free energies, structured JSON reports,
catalyst-ee/product-ee data, and a publication-ready plot.

<table>
  <tr>
    <td><strong>ENSEMBLE ENGINE</strong><br>CREST-ranked conformers with Boltzmann-weighted Gibbs energies</td>
    <td><strong>QUANTUM CORE</strong><br>ORCA optimization and frequency calculations</td>
    <td><strong>SELECTIVITY MODEL</strong><br>Kagan ML2 mass balance with solver safeguards</td>
  </tr>
</table>

### Built for serious screening

- **Explicit execution modes:** real runs never silently become mock runs;
  `--dry-run` previews the calculation plan without executing chemistry.
- **Electronic-state control:** configure molecular charge and spin
  multiplicity for ORCA input generation.
- **Numerical safeguards:** validate XYZ structures, reject invalid
  concentrations, and enforce mass-balance residual tolerances.
- **Machine-readable results:** `nle_report.json` records free energies,
  association constants, NLE classification, 50% catalyst performance, and
  amplification metrics.
- **Audit trail:** `run_manifest.json` records arguments, input hashes, Python
  version, platform, and timestamp.

## Workflow

```mermaid
flowchart LR
    A[RR and RS XYZ structures] --> B[CREST / xTB<br/>conformer sampling]
    B --> C[Lowest-energy<br/>conformers]
    C --> D[ORCA optimization<br/>and frequencies]
    D --> E[Gibbs free energies]
    E --> F[K_homo and K_hetero]
    F --> G[Kagan ML2<br/>mass balance]
    G --> H[CSV data and<br/>NLE curve]
```

### What is calculated

For each catalyst ee value, the model solves the coupled monomer balances

$$
C_R = M_R + 2K_{homo}M_R^2 + K_{hetero}M_RM_S
$$

$$
C_S = M_S + 2K_{homo}M_S^2 + K_{hetero}M_RM_S
$$

and reports the product ee implied by the active monomer pool. The relative
free energy is converted using the configured temperature and
$R = 0.00198720425864083$ kcal mol$^{-1}$ K$^{-1}$.

## Quick start

### 1. Install Python dependencies

```bash
python -m pip install -r requirements.txt
```

For development and testing:

```bash
python -m pip install -r requirements-dev.txt
```

### 2. Install computational chemistry software

Install CREST with xTB and ORCA according to their respective documentation.
Both executables must be available on `PATH`:

```bash
crest --version
orca
```

### 3. Run a real calculation

```bash
python nle_predictor.py \
  --rr inputs/dimer_RR.xyz \
  --rs inputs/dimer_RS.xyz \
  --outdir calc_outputs \
  --solvent toluene \
  --cores 8 \
  --conformers 5 \
  --charge 0 \
  --multiplicity 1
```

Preview the plan without running CREST or ORCA:

```bash
python nle_predictor.py \
  --rr inputs/dimer_RR.xyz \
  --rs inputs/dimer_RS.xyz \
  --dry-run
```

### 4. Run the pipeline demonstration

The mock mode exercises file generation, parsing, numerical solving, CSV
export, and plotting without running CREST or ORCA:

```bash
python nle_predictor.py --mock --outdir calc_outputs
```

Mock energies are synthetic and must not be interpreted as chemical results.
Mock mode must be requested explicitly; a real run stops if CREST or ORCA is
not available.

## Inputs

The `--rr` and `--rs` arguments accept XYZ files containing the homochiral and
heterochiral catalyst dimers. The generated ORCA input currently assumes a
neutral, singlet system (`Charge 0`, `Mult 1`); edit the input-generation
parameters before using charged or open-shell complexes.

| Option | Default | Purpose |
| --- | --- | --- |
| `--rr` | `inputs/dimer_RR.xyz` | Homochiral dimer structure |
| `--rs` | `inputs/dimer_RS.xyz` | Heterochiral dimer structure |
| `--temp` | `298.15` | Temperature in kelvin |
| `--solvent` | none | CREST solvent and ORCA CPCM solvent |
| `--khomo` | `1e5` | Absolute homochiral association constant |
| `--charge` | `0` | Molecular charge passed to ORCA |
| `--multiplicity` | `1` | Spin multiplicity passed to ORCA |
| `--functional` | `B3LYP` | ORCA density functional |
| `--basis` | `def2-SVP` | ORCA basis set |
| `--dispersion` | `D4` | Dispersion correction |
| `--cores` | `4` | CREST and ORCA parallel workers |
| `--conformers` | `1` | Low-energy conformers evaluated per dimer |
| `--outdir` | `calc_outputs` | Calculation and result directory |
| `--mock` | off | Explicitly use deterministic synthetic chemistry outputs |
| `--dry-run` | off | Validate and preview commands without running chemistry |

## Outputs

Every generated artifact is placed below `--outdir`:

```text
calc_outputs/
|-- crest_dimer_RR.xyz/
|   |-- crest.out
|   |-- crest_best.xyz
|   `|-- crest_conf_001.xyz ...
|-- crest_dimer_RS.xyz/
|   |-- crest.out
|   |-- crest_best.xyz
|   `|-- crest_conf_001.xyz ...
|-- opt_RR_conf001.inp / .out
|-- opt_RS_conf001.inp / .out
|-- nle_results.csv
|-- nle_report.json
|-- run_manifest.json
`-- nle_curve.png
```

- `nle_results.csv` contains catalyst ee (%) and product ee (%), ready for
  OriginLab, GraphPad Prism, Excel, or custom analysis.
- `nle_report.json` summarizes the thermodynamics and classifies the result as
  positive NLE, negative NLE, or near-linear.
- `--conformers N` preserves and evaluates up to `N` CREST conformers per
  dimer; the reported free energies are Boltzmann-weighted ensemble values.
- `nle_curve.png` compares the predicted NLE against the ideal linear baseline
  at 600 dpi.
- `opt_*.inp` and `opt_*.out` preserve the ORCA calculation inputs and logs.
- `crest_*/crest_best.xyz` preserves the conformer passed to ORCA.
- `run_manifest.json` records arguments, input hashes, platform, and Python
  version for reproducibility.
- Console logs report $G_{RR}$, $G_{RS}$, relative stability, `K_homo`, and
  `K_hetero`.

## Reproducibility notes

- Record the exact CREST/xTB and ORCA versions used for publication work.
- Preserve the input XYZ files, command line, solvent, temperature, functional,
  basis, dispersion model, and core count with the generated `calc_outputs/`.
- Review ORCA convergence and frequency results before accepting a free energy.
- Real runs fail when a mass-balance point does not converge or has a large
  residual; no approximate chemistry result is returned.
- Dummy hydrogen structures are created only in explicit `--mock` mode when
  input files are missing; real runs require both valid XYZ files.

## Testing

Run the focused regression suite with:

```bash
python -m pytest -q
```

The tests cover XYZ validation, Gibbs-energy parsing, mass-balance residuals,
invalid input rejection, and model endpoint behavior.

## Project layout

```text
.
|-- nle_predictor.py    # Workflow, thermodynamics, solver, export, plotting
|-- requirements.txt    # Python runtime dependencies
|-- requirements-dev.txt # Test and development dependencies
|-- README.md           # Usage and reproducibility documentation
|-- tests/
|   `-- test_nle_predictor.py
`-- assets/
    `-- nle-workflow.svg # Animated project visual
```

## Scope

This project is a workflow scaffold for computational screening. It does not
replace chemical validation, conformer quality assessment, ORCA convergence
review, or experimental measurement. Treat predicted NLE curves as hypotheses
to prioritize, not as standalone evidence of catalytic performance.