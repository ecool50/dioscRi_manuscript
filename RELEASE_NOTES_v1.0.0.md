# dioscRi manuscript v1.0.0

**Submission release.** Snapshot of the analysis code and numerical results for the manuscript:

> *dioscRi enables transferable prediction of clinical outcomes in multi-parameter cytometry data*

This release is the immutable reference cited in the manuscript's Code Availability section. For the `dioscRi` R package itself, see [ecool50/dioscRi](https://github.com/ecool50/dioscRi) (also tagged `v1.0.0`).

## What's in this release

This release contains **code and data only**. The manuscript LaTeX source and response letter are not part of this repository.

### Per-dataset analyses (`analysis_files/`)
- `BioHEART/` — primary discovery + validation analysis (12 notebooks).
- `Breast_Cancer_Wagner_2019/`, `CMV_study_SDY519/`, `COVID_19_PBMC_Mathew_2020/` — external-dataset analyses.
- All hardcoded paths replaced with the `DEEPLEARNING_CYTOF_ROOT` environment variable for portability.

### Revision experiments (`revision/scripts/`)
Scripts added during the revision cycle, organized by topic:
- `uncertainty/` — stratified bootstrap CIs (BioHEART, Breast, COVID, CMV) and the combined uncertainty supplementary figure.
- `ablation/` — 24-variant ablation (normalization × hierarchy × model × features).
- `downsampling/` — cells-per-sample sensitivity (1K–50K).
- `batch_composition/` — Frobenius / batch-structure analyses with panel labels (Supp Figs S6, S7).
- `nref/` — N reference-samples sensitivity (Supp Fig S8).
- `cytogpnet/` — cytoGPNet baseline (data conversion + end-to-end runner; response-only).
- `helpers/` — environment setup and SCE extraction.

All revision scripts use `BIOHEART_ROOT` (env var) for upstream inputs and `here::here()` for outputs, with no machine-specific paths.

### Figures
- Not tracked in this repository — figures are regenerable outputs of the revision scripts (`plot_*.R` under `revision/scripts/`).

### Reproducibility infrastructure
- `revision/results/sessionInfo.txt` — full dependency snapshot (R 4.5.0 with the 18+ packages used in the analyses).
- Environment specification table and step-by-step reproduction workflow in `README.md`.
- `CITATION.cff` — machine-readable citation file.

## Pinned dependencies

| Component | Version |
|---|---|
| R | 4.5.0 |
| `keras3` (R) | 1.2.0 |
| `tensorflow` (Python) | 2.16.2 |
| Python | 3.10.15 (virtualenv `r-reticulate`) |
| cytoGPNet env | PyTorch 2.0, gpytorch 1.10 (separate conda env, see `revision/scripts/cytogpnet/`) |

## Data

Raw and processed inputs are archived on [Zenodo](https://zenodo.org/records/15694581).

## License

GPL-3.0
