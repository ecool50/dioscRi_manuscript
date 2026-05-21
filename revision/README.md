# dioscRi Revision Materials

Everything needed for the reviewer response: code, results, figures, and
response documents. Organized by reviewer concern so each item maps cleanly to
a point in the response letter.

## Layout

```
revision/
├── docs/                       # Plan, tracker, response letter
│   ├── reviewer_response_plan.md       Plan + reviewer-by-reviewer status
│   ├── revision_tracker.md             Tracker mapping points to scripts/results
│   ├── response_to_reviewers.md        Full response letter (markdown)
│   ├── response_to_reviewers.docx
│   ├── reviewer1_summary.{md,pdf}
│   ├── revision_summary.{md,pdf}
│   ├── bootstrap_ci_justification.md   Methodology defense for uncertainty CIs
│   ├── dioscRI_Reviewer_Comments.docx  Original reviewer comments
│   ├── reviewer_response_plan.docx
│   └── latex/                          revisions.tex + response_to_reviewers.tex
├── scripts/                    # Code (R + Python + shell)
│   ├── uncertainty/    R2 Major 10  — bootstrap CIs + repeated splits
│   ├── ablation/       R2 Major 8   — ablation study (10 variants)
│   ├── downsampling/   R2 Major 4   — cells-per-sample sensitivity (1K–50K)
│   ├── batch_composition/  R1 Major 1 / R2 Major 2 — reference sample batches
│   ├── nref/           R2 Major 3   — K=2 vs K=6–7 sensitivity
│   ├── cytogpnet/      R2 Major 6   — cytoGPNet baseline comparison
│   └── helpers/        setup_environment.R, extract_discovery_from_full_sce.R
├── results/                    # Output CSVs / pickled objects
│   ├── uncertainty/    bootstrap_ci.csv + bootstrap_distribution.csv per dataset
│   ├── ablation/       bioheart_ablation_results.csv + no_norm / basic_norm
│   ├── downsampling/   per-cell-count results + per-sample ARI + F1
│   ├── batch_composition/  reference_sample_batch_composition.csv
│   ├── nref/           K-sensitivity sweeps across datasets
    └── cytogpnet/bioheart/  Pickled inputs + .pth checkpoints + predictions
```

Figures are regenerated locally by the relevant `plot_*.R` scripts and are not tracked in the repository.

## Path conventions

Scripts resolve two roots:

- **`bioheart_root`** — upstream analysis directory containing raw data,
  pre-normalized SCE objects, and trained model weights. Defaults to
  `~/Documents/Academic/PhD/bioheart_analysis`. Override with
  `BIOHEART_ROOT` env var when running on a different machine.
- **`revision_root`** — this folder (`dioscRi_manuscript/revision/`),
  resolved via `here::here("revision")` (R) or relative to the script
  location (Python/shell).

**Inputs** (raw data, SCEs, model weights) load from `bioheart_root`.
**Outputs** (CSVs, figures) write to `revision_root`. Each script also keeps
`base_dir <- bioheart_root` as a legacy alias.

## Reproducing

```bash
# One-time
Rscript scripts/helpers/setup_environment.R    # installs here, keras3 1.2.0, TF 2.16.2

# Set BIOHEART_ROOT if your upstream data lives elsewhere
export BIOHEART_ROOT=/path/to/bioheart_analysis

# Run individual experiments
Rscript scripts/uncertainty/bioheart_bootstrap_ci.R
Rscript scripts/uncertainty/breast_cancer_bootstrap_ci.R
Rscript scripts/uncertainty/covid_bootstrap_ci.R
Rscript scripts/uncertainty/cmv_bootstrap_ci.R
Rscript scripts/uncertainty/plot_uncertainty_supp_figure.R   # combined fig

Rscript scripts/ablation/bioheart_ablation.R
Rscript scripts/downsampling/bioheart_downsampling_sensitivity.R
Rscript scripts/batch_composition/batch_composition_check.R

# cytoGPNet — end-to-end (clone, convert, pretrain, train, test, AUC)
#   One-time: conda env create -f .cytoGPNet/environment.yml
bash scripts/cytogpnet/run_cytogpnet.sh bioheart
#   EPOCHS=50 bash scripts/cytogpnet/run_cytogpnet.sh bioheart  # override
```

## Cross-reference

`docs/revision_tracker.md` is the live status table. `docs/reviewer_response_plan.md`
has the per-point response strategy. Both reference scripts and results by their
paths within this folder.
