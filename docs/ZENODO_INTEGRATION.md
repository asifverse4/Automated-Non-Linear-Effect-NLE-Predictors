# Zenodo DOI Integration Guide

This guide outlines how to link **NLE-Predictor** with **Zenodo** (hosted by CERN) to mint permanent, citable Digital Object Identifiers (DOIs) for each software release.

---

## Why Zenodo?

1. **Permanent Academic Archive:** Zenodo preserves snapshots of the codebase in the CERN data center.
2. **Citable DOI:** Assigns a persistent DOI (e.g., `10.5281/zenodo.XXXXXXX`) indexed across DataCite, Crossref, Google Scholar, and Web of Science.
3. **Automated Release Archiving:** Whenever a GitHub release or tag is published, Zenodo automatically archives the release assets (wheels, source tarballs, manuscripts, standalone binaries).

---

## Setup Instructions (Takes 60 Seconds)

### Step 1: Log in to Zenodo via GitHub
1. Navigate to [https://zenodo.org](https://zenodo.org).
2. Click **Log in** at the top right and select **Log in with GitHub**.
3. Authorize Zenodo to access your public GitHub repositories.

### Step 2: Enable NLE-Predictor Repository
1. Navigate to your Zenodo GitHub settings:
   👉 **[https://zenodo.org/account/settings/github/](https://zenodo.org/account/settings/github/)**
2. Locate `asifverse4/Automated-Non-Linear-Effect-NLE-Predictors` in the repository list.
3. Toggle the switch next to the repository to **ON** (Enabled).

### Step 3: Trigger DOI Minting
Once enabled, Zenodo will automatically archive the repository and mint a DOI on the next GitHub release or when a release tag (e.g., `v1.0.0` or `v1.0.1`) is published.

---

## Metadata Configuration

The repository contains pre-configured metadata files that Zenodo parses automatically:
- **[`.zenodo.json`](../.zenodo.json):** Provides the official title, author details (Asif Raza, Jamia Hamdard), abstract, license (MIT), and keywords.
- **[`CITATION.cff`](../CITATION.cff):** Citation File Format metadata compliant with GitHub's native "Cite this repository" widget.

---

## Adding the Badge to README

Once your repository is enabled on Zenodo and the first release is archived:
1. Zenodo will generate a badge markdown snippet:
   ```markdown
   [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
   ```
2. Place this badge in the top badge section of [`README.md`](../README.md).
