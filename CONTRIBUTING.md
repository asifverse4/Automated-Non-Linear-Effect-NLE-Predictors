# Contributing to NLE // ATLAS

Thank you for helping improve NLE // ATLAS. Contributions should strengthen scientific reproducibility, computational safety, or usability for asymmetric-catalysis research.

## Before You Start

- Search existing issues before opening a new one.
- For substantial changes, open an issue first to describe the proposed design.
- Never include private molecular data, credentials, proprietary calculation files, or unpublished results without permission.
- Keep scientific claims proportional to the evidence available.

## Development Setup

```bash
python -m pip install -r requirements-dev.txt
python -m pytest -q
```

The test suite must pass before submitting a change. Real CREST/xTB and ORCA calculations are not required for unit tests; use explicit `--mock` mode for workflow demonstrations.

## Making Changes

- Keep changes focused and consistent with the existing Python style.
- Add or update tests for behavior changes.
- Preserve deterministic mock behavior.
- Update `README.md` when CLI options, outputs, or scientific assumptions change.
- Document new approximations and clearly distinguish them from experimental or quantum-chemistry validation.
- Do not silently convert a failed real calculation into mock output.

## Pull Requests

A useful pull request includes:

- A concise description of the scientific or engineering motivation.
- The exact tests and commands run.
- Any changes to model assumptions or thermodynamic conventions.
- Limitations, benchmark gaps, or known follow-up work.
- Representative output only when it is synthetic, public, or explicitly approved.

Keep commits small and descriptive. Review generated files for accidental input structures, calculation outputs, caches, or secrets before pushing.

## Reporting Bugs

Include the command, Python version, operating system, dependency versions, relevant log excerpts, and a minimal reproducible input where possible. Remove confidential molecular or computational data before sharing.

## Scientific Contributions

Changes to the Kagan model, association-constant convention, conformer weighting, or free-energy treatment require tests and a short explanation of the underlying reference or derivation. A passing test confirms software behavior; it does not by itself establish chemical validity.
