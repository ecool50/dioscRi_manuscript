#!/usr/bin/env python3
"""
Convert BioHEART-CT data to cytoGPNet format.

cytoGPNet expects a pickle .obj file with:
  - expr_list: np.array [n_samples, n_channels, n_cells, n_markers]
    (n_channels = 1 for cross-sectional data)
  - cytof_files: pd.DataFrame with patient metadata including outcome
  - marker_names: list of marker names

We create separate train (discovery) and test (validation) .obj files.
"""

import pandas as pd
import numpy as np
import pickle
import os

bioheart_root = os.environ.get(
    'BIOHEART_ROOT',
    os.path.expanduser('~/Documents/Academic/PhD/bioheart_analysis')
)
# Outputs land under revision/results/cytogpnet/bioheart/
revision_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
out_dir = os.path.join(revision_root, 'results/cytogpnet/bioheart')
os.makedirs(out_dir, exist_ok=True)
base_dir = bioheart_root  # alias: upstream inputs (raw_data, misc/Study_3_2019)

# Markers (underscore format matching CSV)
markers = [
    'HLA_DR', 'CD3', 'CD4', 'CD8a', 'CD25', 'CD127', 'FoxP3', 'CD27',
    'KLRG1', 'CD56', 'CD45RO', 'CD45RA', 'CD192_CCR2', 'CD194_CCR4',
    'CD196_CCR6', 'CD39', 'CD38', 'Ki67', 'CD183_CXCR3', 'CCR7',
    'CD19', 'CD20', 'IgD', 'CD14', 'CD304', 'CD141', 'CD1c_PE'
]
num_markers = len(markers)

# --- Load discovery cohort (training) ---
print("Loading discovery cohort...")
# Try the clean copy first, fall back to /tmp
csv_path = os.path.join(base_dir, 'data/raw_data/bioheart_b4_clean.csv')
if not os.path.exists(csv_path):
    csv_path = '/tmp/bioheart_b4.csv'
if not os.path.exists(csv_path):
    # Try original
    csv_path = os.path.join(base_dir, 'data/raw_data/bioheart_ct_cytof_data_b4_mg.csv')

df_disc = pd.read_csv(csv_path)
print(f"  Discovery: {len(df_disc)} cells, {df_disc['sample_id'].nunique()} samples")

# --- Load validation cohort ---
print("Loading validation cohort...")
val_csv = os.path.join(base_dir, 'misc/Study_3_2019/all_batches_processed_mg_10K.csv')
df_val = pd.read_csv(val_csv)
print(f"  Validation: {len(df_val)} cells, {df_val['sample_id'].nunique()} samples")

# --- arcsinh transform (cofactor 5 for CyTOF, same as our pipeline) ---
print("Applying arcsinh transform...")
for m in markers:
    df_disc[m] = np.arcsinh(df_disc[m] / 5.0)
    df_val[m] = np.arcsinh(df_val[m] / 5.0)


def convert_to_obj(df, markers, outcome_col, outcome_map, obj_path, max_cells=10000):
    """
    Convert a DataFrame to cytoGPNet .obj format.

    Args:
        df: DataFrame with sample_id column and marker columns
        markers: list of marker column names
        outcome_col: column name for the outcome variable
        outcome_map: dict mapping outcome values to 0/1
        obj_path: output path for the .obj file
        max_cells: pad/truncate to this many cells per sample
    """
    samples = df['sample_id'].unique()
    n_samples = len(samples)

    expr_list = np.zeros((n_samples, 1, max_cells, len(markers)), dtype=np.float32)
    metadata_rows = []

    for i, sid in enumerate(samples):
        sample_df = df[df['sample_id'] == sid]
        expr = sample_df[markers].values
        n_cells = expr.shape[0]

        if n_cells >= max_cells:
            # Subsample without replacement
            idx = np.random.choice(n_cells, max_cells, replace=False)
            expr_list[i, 0, :, :] = expr[idx]
        else:
            # Pad with zeros
            expr_list[i, 0, :n_cells, :] = expr

        # Get outcome for this sample
        outcome_val = sample_df[outcome_col].iloc[0]
        if isinstance(outcome_map, dict):
            y = outcome_map.get(outcome_val, outcome_val)
        else:
            y = outcome_val

        metadata_rows.append({
            'patient_id': sid,
            'label': int(y),
            'y': int(y)
        })

    cytof_files = pd.DataFrame(metadata_rows)

    all_data = {
        'expr_list': expr_list,
        'cytof_files': cytof_files,
        'marker_names': markers
    }

    with open(obj_path, 'wb') as f:
        pickle.dump(all_data, f)

    print(f"  Saved to {obj_path}")
    print(f"  Shape: {expr_list.shape}")
    print(f"  Outcome distribution: {cytof_files['y'].value_counts().to_dict()}")


# --- Convert discovery (training) ---
print("\nConverting discovery cohort...")
# Gensini > 0 = CAD+ (1), Gensini == 0 = CAD- (0)
df_disc['y'] = (df_disc['Gensini'] > 0).astype(int)
convert_to_obj(
    df=df_disc,
    markers=markers,
    outcome_col='y',
    outcome_map=None,
    obj_path=os.path.join(out_dir, 'train_Data.obj'),
    max_cells=10000
)

# --- Convert validation (testing) ---
print("\nConverting validation cohort...")
# Gensini_bin is already 0/1
convert_to_obj(
    df=df_val,
    markers=markers,
    outcome_col='Gensini_bin',
    outcome_map=None,
    obj_path=os.path.join(out_dir, 'test_Data.obj'),
    max_cells=10000
)

# --- Verify ---
print("\n=== Verification ===")
for name in ['train_Data.obj', 'test_Data.obj']:
    path = os.path.join(out_dir, name)
    with open(path, 'rb') as f:
        data = pickle.load(f)
    print(f"\n{name}:")
    print(f"  expr_list shape: {data['expr_list'].shape}")
    print(f"  cytof_files shape: {data['cytof_files'].shape}")
    print(f"  markers: {data['marker_names'][:5]}...")
    print(f"  outcome: {data['cytof_files']['y'].value_counts().to_dict()}")
    print(f"  any NaN: {np.any(np.isnan(data['expr_list']))}")

print("\nDone. Files saved to:", out_dir)
