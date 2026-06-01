# dioscRi: manuscript analysis code

Code and reproduction materials for the manuscript

> *dioscRi enables transferable prediction of clinical outcomes in multi-parameter cytometry data*

This repository contains everything needed to reproduce the analyses, figures, and supplementary materials, including the experiments added for the revision (uncertainty quantification, ablation, downsampling sensitivity, batch composition, reference-sample selection, and the cytoGPNet baseline).

For the **dioscRi R package itself**, see [github.com/ecool50/dioscRi](https://github.com/ecool50/dioscRi). The package and this repository are versioned together; the release tag `v1.0.0` (created at acceptance) is the immutable reference for this manuscript.

---

## Repository layout

This repository contains **code and data only**. The manuscript LaTeX source, response letter, and other write-up artefacts are not tracked here.

```
dioscRi_manuscript/
├── analysis_files/              Per-dataset analysis notebooks
│   ├── BioHEART/                  Bioheart_Prediction.Rmd is the primary analysis;
│   │                              other Rmds (Lambda_experiment, bottleneck_experiment,
│   │                              Number_of_ref_samples_experiment, Number_of_clusters_experiment,
│   │                              cyCombine / iMUBAC / MG variants) are supporting experiments
│   │                              referenced in supplementary figures and are not required
│   │                              for the main results.
│   ├── Breast_Cancer_Wagner_2019/ breast_cancer_prediction.Rmd + cluster-number experiment
│   ├── CMV_study_SDY519/          cmv_prediction.Rmd + cluster-number experiment
│   └── COVID_19_PBMC_Mathew_2020/ bcell_prediction.Rmd + cluster-number experiment
└── revision/                    Reproducible scripts and outputs added during the revision cycle
    ├── scripts/                 Scripts organised by topic
    │   ├── uncertainty/           Bootstrap CIs + Supp Fig 9 plotter
    │   ├── ablation/              24-variant ablation (Supp Fig 11)
    │   ├── downsampling/          Cells-per-sample sensitivity (Supp Fig 10)
    │   ├── batch_composition/     Frobenius / batch-structure analyses (Supp Figs 6, 7;
    │   │                          batch_composition_check.R is a precursor diagnostic)
    │   ├── nref/                  Reference-sample sensitivity (Supp Fig 8 + plot_nref_sensitivity.R
    │   │                          for the per-tier AUC distribution)
    │   ├── cytogpnet/             cytoGPNet baseline runner + plot_cytogpnet_comparison.R
    │   │                          (Response Fig R1, not in the manuscript itself)
    │   └── helpers/               setup_environment.R, fetch_data.R,
    │                              extract_discovery_from_full_sce.R (one-off data prep)
    └── results/                 Output CSVs/RDS, supplementary numerical results
        └── sessionInfo.txt        Dependency snapshot (R 4.5.0)
```

Figures are regenerated locally by the relevant scripts (`revision/scripts/.../plot_*.R`); they are not tracked in this repository.

## Data

Raw and processed inputs are archived on **[Zenodo](https://zenodo.org/records/15694581)** as a single ~7.2 GB archive. The helper

```bash
export BIOHEART_ROOT=/path/to/bioheart_analysis
Rscript revision/scripts/helpers/fetch_data.R
```

downloads it into `$BIOHEART_ROOT/data/`, verifies the md5, and unzips it (resumable via `--force`). After unzipping, the layout is preserved as `dioscRi_analysis_data/` under `$BIOHEART_ROOT/data/`.

External cytometry inputs used by some analyses live outside the manuscript repo (e.g. CMV data from the DeepLearningCyTOF distribution); paths are documented at the top of the relevant scripts.

## Environment

Pinned versions used for the manuscript:

| Component | Version |
|---|---|
| R | 4.5.0 |
| `keras3` (R) | 1.2.0 |
| `tensorflow` (Python) | 2.16.2 |
| Python | 3.10.15 (virtualenv `r-reticulate`) |

One-time setup (R packages + matched TensorFlow/Keras):

```r
Rscript revision/scripts/helpers/setup_environment.R
```

A full `sessionInfo()` snapshot is at `revision/results/sessionInfo.txt`.

For the cytoGPNet baseline only, a separate conda environment is needed (PyTorch 2.0 + gpytorch 1.10, named `cytoGPNet`). See `revision/scripts/cytogpnet/run_cytogpnet.sh`.

## Path convention

All scripts in this repository resolve their data locations via environment variables, so nothing is hardcoded to one machine:

- `BIOHEART_ROOT`, used by `revision/scripts/`. Points to the upstream analysis directory containing raw data, SCE objects, and model weights. Defaults to `~/Documents/Academic/PhD/bioheart_analysis`.
- `DEEPLEARNING_CYTOF_ROOT`, used by `analysis_files/{CMV_study_SDY519,Breast_Cancer_Wagner_2019,COVID_19_PBMC_Mathew_2020}/*.Rmd`. Points to the parent directory containing the DeepLearningCyTOF data distribution. Defaults to `~/Documents/PhD/DeepLearning_CyTOF`.
- `here::here()`, resolves to this repository root (`dioscRi_manuscript/`). All script *outputs* are written under `revision/results/` or `revision/figures/` relative to this root.

```bash
export BIOHEART_ROOT=/path/to/bioheart_analysis
export DEEPLEARNING_CYTOF_ROOT=/path/to/DeepLearning_CyTOF
Rscript revision/scripts/uncertainty/bioheart_bootstrap_ci.R
Rscript -e 'rmarkdown::render("analysis_files/CMV_study_SDY519/cmv_prediction.Rmd")'
```

## Step-by-step reproduction

```bash
# 1. Clone the package and this repository side-by-side
git clone https://github.com/ecool50/dioscRi.git
git clone https://github.com/ecool50/dioscRi_manuscript.git

# 2. Install the R environment + matched TF/Keras
Rscript dioscRi_manuscript/revision/scripts/helpers/setup_environment.R

# 3. Install the dioscRi package
Rscript -e 'devtools::install_github("ecool50/dioscRi")'

# 4. Download data (~7.2 GB) from Zenodo into $BIOHEART_ROOT/data/
export BIOHEART_ROOT=/path/to/bioheart_analysis
Rscript dioscRi_manuscript/revision/scripts/helpers/fetch_data.R

# 5. Reproduce the main analyses (per dataset)
#    BioHEART-CT, the primary discovery + validation analysis:
Rscript -e 'rmarkdown::render("dioscRi_manuscript/analysis_files/BioHEART/Bioheart_Prediction.Rmd")'
#    The other three datasets (Breast Cancer, CMV, COVID-19) each have a
#    *_prediction.Rmd that uses DEEPLEARNING_CYTOF_ROOT for raw inputs:
export DEEPLEARNING_CYTOF_ROOT=/path/to/DeepLearning_CyTOF
Rscript -e 'rmarkdown::render("dioscRi_manuscript/analysis_files/Breast_Cancer_Wagner_2019/breast_cancer_prediction.Rmd")'
Rscript -e 'rmarkdown::render("dioscRi_manuscript/analysis_files/CMV_study_SDY519/cmv_prediction.Rmd")'
Rscript -e 'rmarkdown::render("dioscRi_manuscript/analysis_files/COVID_19_PBMC_Mathew_2020/bcell_prediction.Rmd")'

# 6. Reproduce the revision experiments (each script is standalone)
Rscript dioscRi_manuscript/revision/scripts/uncertainty/bioheart_bootstrap_ci.R
Rscript dioscRi_manuscript/revision/scripts/uncertainty/breast_cancer_bootstrap_ci.R
Rscript dioscRi_manuscript/revision/scripts/uncertainty/covid_bootstrap_ci.R
Rscript dioscRi_manuscript/revision/scripts/uncertainty/cmv_bootstrap_ci.R
Rscript dioscRi_manuscript/revision/scripts/uncertainty/plot_uncertainty_supp_figure.R   # Supp Fig 9

Rscript dioscRi_manuscript/revision/scripts/ablation/bioheart_ablation.R                  # Supp Fig 11
Rscript dioscRi_manuscript/revision/scripts/downsampling/bioheart_downsampling_sensitivity.R  # Supp Fig 10
Rscript dioscRi_manuscript/revision/scripts/batch_composition/plot_s6_bioheart_frobenius.R    # Supp Fig 6
Rscript dioscRi_manuscript/revision/scripts/batch_composition/plot_s7_cmv_frobenius.R         # Supp Fig 7
Rscript dioscRi_manuscript/revision/scripts/nref/plot_s8_tier_combined.R                      # Supp Fig 8
Rscript dioscRi_manuscript/revision/scripts/nref/plot_nref_sensitivity.R                      # per-tier AUC detail

# 7. Reproduce the cytoGPNet baseline (separate conda env)
bash dioscRi_manuscript/revision/scripts/cytogpnet/run_cytogpnet.sh bioheart
Rscript dioscRi_manuscript/revision/scripts/cytogpnet/plot_cytogpnet_comparison.R         # Response Fig R1
```

## License

Code in this repository is released under GPL-3 (matching the `dioscRi` package).

## Contact

- Issues with scripts in this repository: [open an issue](https://github.com/ecool50/dioscRi_manuscript/issues)
- Issues with the `dioscRi` package itself: [package issues](https://github.com/ecool50/dioscRi/issues)
- General correspondence: ewil3501@uni.sydney.edu.au
