#!/usr/bin/env bash
#
# Run the cytoGPNet baseline on a converted dioscRi dataset.
#
# Pipeline:
#   1. (one-time) clone the upstream cytoGPNet repo at https://github.com/llin-lab/cytoGPNet
#   2. (one-time) convert dioscRi data to cytoGPNet pickle format
#   3. stage fold<N>/ subdir (cytoGPNet's train_simplified.py / test.py expect it)
#   4. pretrain the autoencoder           -> ae_output/simpleAE_epoch<E>.pth
#   5. fine-tune GP + attention + clf     -> model_output/<component>_epoch<E>.pth
#   6. test                               -> model_output/test_predictions<F>.csv
#   7. compute AUC and append to results/cytogpnet/auc_summary.csv
#
# Prereqs:
#   conda env create -f $CYTOGPNET_DIR/environment.yml   # one-time, named "cytoGPNet"
#
# Usage:
#   bash run_cytogpnet.sh [bioheart]                 # dataset (default: bioheart)
#   EPOCHS=50 bash run_cytogpnet.sh bioheart         # override epochs (default: 100)
#   CYTOGPNET_DIR=/path/to/cytoGPNet ./run...        # use existing clone

set -euo pipefail

DATASET="${1:-bioheart}"
EPOCHS="${EPOCHS:-100}"
FOLD="${FOLD:-1}"
CONDA_ENV="${CONDA_ENV:-cytoGPNet}"

REVISION_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
DATA_ROOT="${REVISION_ROOT}/results/cytogpnet/${DATASET}"
CYTOGPNET_DIR="${CYTOGPNET_DIR:-${REVISION_ROOT}/.cytoGPNet}"
CYTOGPNET_CODE="${CYTOGPNET_DIR}/cytoGPNet model"

echo "=========================================="
echo "  cytoGPNet runner"
echo "  dataset:        $DATASET"
echo "  epochs:         $EPOCHS"
echo "  fold:           $FOLD"
echo "  data root:      $DATA_ROOT"
echo "  cytoGPNet repo: $CYTOGPNET_DIR"
echo "  conda env:      $CONDA_ENV"
echo "=========================================="

# Run a python command inside the cytoGPNet conda env. `conda run` keeps the
# parent shell clean (no `conda activate` needed in this script).
crun() { conda run --no-capture-output -n "$CONDA_ENV" "$@"; }

# ---- 1. Clone upstream cytoGPNet if absent ----
if [ ! -d "$CYTOGPNET_DIR" ]; then
  echo "[setup] cloning cytoGPNet to $CYTOGPNET_DIR"
  git clone https://github.com/llin-lab/cytoGPNet.git "$CYTOGPNET_DIR"
fi
if [ ! -f "$CYTOGPNET_CODE/pretrain.py" ]; then
  echo "ERROR: expected $CYTOGPNET_CODE/pretrain.py -- is CYTOGPNET_DIR pointing at the repo root?"
  exit 1
fi

# ---- 2. Convert dioscRi data to cytoGPNet pickle format ----
CONVERT_SCRIPT="$REVISION_ROOT/scripts/cytogpnet/convert_${DATASET}_for_cytogpnet.py"
if [ ! -f "$CONVERT_SCRIPT" ]; then
  echo "ERROR: no conversion script for dataset '$DATASET' at $CONVERT_SCRIPT"
  echo "       Only 'bioheart' is wired up today. Add a converter and rerun."
  exit 1
fi
if [ ! -f "$DATA_ROOT/train_Data.obj" ] || [ ! -f "$DATA_ROOT/test_Data.obj" ]; then
  echo "[1/4] converting $DATASET to cytoGPNet pickle format"
  crun python "$CONVERT_SCRIPT"
else
  echo "[1/4] $DATASET pickles already present, skipping conversion"
fi

# ---- 3. Stage fold<N>/ subdir (cytoGPNet expects this for train/test) ----
FOLD_DIR="$DATA_ROOT/fold${FOLD}"
mkdir -p "$FOLD_DIR"
ln -sf "$DATA_ROOT/train_Data.obj" "$FOLD_DIR/train_Data.obj"
ln -sf "$DATA_ROOT/test_Data.obj"  "$FOLD_DIR/test_Data.obj"

# ---- 4. Pretrain the autoencoder ----
AE_OUT="$DATA_ROOT/ae_output"
AE_CKPT="$AE_OUT/simpleAE_epoch${EPOCHS}.pth"
mkdir -p "$AE_OUT"
if [ -f "$AE_CKPT" ]; then
  echo "[2/4] AE checkpoint already exists at $AE_CKPT, skipping pretrain (delete to force rerun)"
else
  echo "[2/4] pretraining AE ($EPOCHS epochs)"
  (
    cd "$CYTOGPNET_CODE"
    crun python pretrain.py \
      --data-dir "$DATA_ROOT" \
      --save-dir "$AE_OUT" \
      --max-epochs "$EPOCHS"
  )
fi

# ---- 5. Fine-tune GP + attention + classifier ----
MODEL_OUT="$DATA_ROOT/model_output"
mkdir -p "$MODEL_OUT"
echo "[3/4] fine-tuning GP + attention + classifier ($EPOCHS epochs)"
(
  cd "$CYTOGPNET_CODE"
  crun python train_simplified.py \
    --data-dir "$DATA_ROOT" \
    -fold "$FOLD" \
    --save-dir "$MODEL_OUT" \
    --pretrained-file "$AE_OUT/simpleAE_epoch${EPOCHS}.pth" \
    --max-epochs "$EPOCHS"
)

# ---- 6. Test ----
echo "[4/4] testing"
(
  cd "$CYTOGPNET_CODE"
  crun python test.py \
    --data-dir "$DATA_ROOT" \
    -fold "$FOLD" \
    --save-dir "$MODEL_OUT" \
    --max-epochs "$EPOCHS"
)

# ---- 7. Compute AUC and log it ----
PRED_CSV="$MODEL_OUT/test_predictions${FOLD}.csv"
SUMMARY_CSV="$REVISION_ROOT/results/cytogpnet/auc_summary.csv"
crun python - <<PY
import csv, os, pandas as pd
from datetime import datetime
from sklearn.metrics import roc_auc_score

df = pd.read_csv("$PRED_CSV")
auc = roc_auc_score(df["label"], df["prediction_score"])

print()
print(f"=== cytoGPNet result ===")
print(f"  dataset:      $DATASET")
print(f"  fold:         $FOLD")
print(f"  epochs:       $EPOCHS")
print(f"  n test:       {len(df)}")
print(f"  prevalence:   {df['label'].mean():.3f}")
print(f"  AUC:          {auc:.4f}")
print(f"  predictions:  $PRED_CSV")

summary = "$SUMMARY_CSV"
write_header = not os.path.exists(summary)
with open(summary, "a", newline="") as f:
    w = csv.writer(f)
    if write_header:
        w.writerow(["timestamp", "dataset", "fold", "epochs", "n_test", "prevalence", "auc"])
    w.writerow([datetime.now().isoformat(timespec="seconds"),
                "$DATASET", "$FOLD", "$EPOCHS", len(df), f"{df['label'].mean():.4f}", f"{auc:.4f}"])
print(f"  appended to:  {summary}")
PY
