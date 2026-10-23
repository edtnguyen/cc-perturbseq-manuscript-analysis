# Notebook mirror: `oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/PAH.PerturbSeq/present_figures/IP_GFP_assigned_guides_per_cell_histogram.ipynb`

- Source archive path: `/oak/stanford/groups/engreitz/Users/tri/CodeBase/cc-perturbseq-figure-reproduction/recovery/sherlock_notebooks/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/PAH.PerturbSeq/present_figures/IP_GFP_assigned_guides_per_cell_histogram.ipynb`
- SHA-256: `f8950f39b2009b18eab3762bd3806d79168587e449db6d7c66f869af968f75ee`
- Kernel language: `python`
- Full-fidelity notebook storage: Git LFS

## Markdown cell 0

# IP_GFP: number of SCEPTRE-assigned guides per cell

This notebook calculates the number of **SCEPTRE-assigned guides** per IP_GFP cell and exports:

1. one row per cell;
2. a histogram-ready count table;
3. a compact summary table.

The default input is the final K30/density-threshold-2.0 H5MU. If that file is not found, the notebook falls back to the prepared PerturbNMF H5AD.

Because these objects contain the guide-positive analysis subset, the distribution normally starts at one guide per cell. To include zero-guide cells, use the full RNA-QC object with assignments instead.

## Code cell 1

```python
from __future__ import annotations

from pathlib import Path
import json
import sys

import anndata as ad
from IPython.display import display
import muon as mu
import numpy as np
import pandas as pd
import scipy.sparse as sp

print("Python:", sys.version.split()[0])

```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/muon/_core/preproc.py:31: FutureWarning: `__version__` is deprecated, use `importlib.metadata.version('scanpy')` instead
  if Version(scanpy.__version__) < Version("1.10"):

```

### Text output (stream)

```text
Python: 3.11.14

```

## Code cell 2

```python
# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------
BASE = Path(
    "/oak/stanford/groups/engreitz/Users/tri/10xPerturbSeq/"
    "analysis/250507_ePAH1_NovaSeq/CRISPR_count/"
    "IP_GV_unified_sceptre_cnmf_no_celltype"
)

H5MU_PATH = (
    BASE
    / "PerturbNMF"
    / "results"
    / "IP_GV_skcnmf_K30_dt2_0_n100"
    / "Inference"
    / "adata"
    / "cNMF_30_2_0.h5mu"
)

H5AD_FALLBACK = (
    BASE
    / "PerturbNMF"
    / "input"
    / "IP_GFP_GV_GFP_guide_positive_PerturbNMF_input.h5ad"
)

SAMPLE_NAME = "IP_GFP"

OUTPUT_DIR = (
    BASE
    / "PerturbNMF"
    / "results"
    / "IP_GV_skcnmf_K30_dt2_0_n100"
    / "Evaluation"
    / "30_2_0"
    / "guide_count_histogram"
)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("Preferred input:", H5MU_PATH)
print("Fallback input:", H5AD_FALLBACK)
print("Output directory:", OUTPUT_DIR)

```

### Text output (stream)

```text
Preferred input: /oak/stanford/groups/engreitz/Users/tri/10xPerturbSeq/analysis/250507_ePAH1_NovaSeq/CRISPR_count/IP_GV_unified_sceptre_cnmf_no_celltype/PerturbNMF/results/IP_GV_skcnmf_K30_dt2_0_n100/Inference/adata/cNMF_30_2_0.h5mu
Fallback input: /oak/stanford/groups/engreitz/Users/tri/10xPerturbSeq/analysis/250507_ePAH1_NovaSeq/CRISPR_count/IP_GV_unified_sceptre_cnmf_no_celltype/PerturbNMF/input/IP_GFP_GV_GFP_guide_positive_PerturbNMF_input.h5ad
Output directory: /oak/stanford/groups/engreitz/Users/tri/10xPerturbSeq/analysis/250507_ePAH1_NovaSeq/CRISPR_count/IP_GV_unified_sceptre_cnmf_no_celltype/PerturbNMF/results/IP_GV_skcnmf_K30_dt2_0_n100/Evaluation/30_2_0/guide_count_histogram

```

## Code cell 3

```python
def choose_first(container, keys):
    for key in keys:
        if key in container:
            return key
    return None


def row_sums(matrix) -> np.ndarray:
    if sp.issparse(matrix):
        return np.asarray(matrix.sum(axis=1)).ravel()
    return np.asarray(matrix).sum(axis=1)


# ---------------------------------------------------------------------
# Load the object
# ---------------------------------------------------------------------
if H5MU_PATH.exists():
    source_path = H5MU_PATH
    source_kind = "h5mu"
    mdata = mu.read_h5mu(source_path)

    if "cNMF" in mdata.mod:
        obj = mdata["cNMF"]
        modality_name = "cNMF"
    elif "rna" in mdata.mod:
        obj = mdata["rna"]
        modality_name = "rna"
    else:
        raise KeyError(
            f"No cNMF or RNA modality found. Modalities: {list(mdata.mod)}"
        )

elif H5AD_FALLBACK.exists():
    source_path = H5AD_FALLBACK
    source_kind = "h5ad"
    obj = ad.read_h5ad(source_path)
    modality_name = "AnnData"

else:
    raise FileNotFoundError(
        f"Neither input exists:\n{H5MU_PATH}\n{H5AD_FALLBACK}"
    )

print("Loaded:", source_path)
print("Object:", modality_name, obj.shape)

# ---------------------------------------------------------------------
# Resolve sample labels
# ---------------------------------------------------------------------
sample_key = choose_first(
    obj.obs.columns,
    ["sample", "sample_id", "mouse"],
)
if sample_key is None:
    raise KeyError(
        "Could not find a sample column among: sample, sample_id, mouse"
    )

sample_values = obj.obs[sample_key].astype(str)
ip_mask = (sample_values == SAMPLE_NAME).to_numpy()

if ip_mask.sum() == 0:
    raise ValueError(
        f"No cells have {sample_key} == {SAMPLE_NAME!r}. "
        f"Observed values: {sorted(sample_values.unique().tolist())}"
    )

# ---------------------------------------------------------------------
# Resolve the binary guide-assignment matrix
# ---------------------------------------------------------------------
guide_key = choose_first(
    obj.obsm.keys(),
    ["guide_assignment", "guide_assignments"],
)
if guide_key is None:
    raise KeyError(
        "Could not find obsm['guide_assignment'] or "
        "obsm['guide_assignments']."
    )

guide_assignment = obj.obsm[guide_key]

if sp.issparse(guide_assignment):
    guide_assignment = guide_assignment.tocsr().copy()
    # Count assigned guides, regardless of the stored nonzero value.
    binary_assignment = guide_assignment.copy()
    binary_assignment.data = np.ones(
        binary_assignment.nnz,
        dtype=np.int8,
    )
else:
    binary_assignment = (
        np.asarray(guide_assignment) != 0
    ).astype(np.int8)

n_assigned_guides = row_sums(binary_assignment).astype(int)

if len(n_assigned_guides) != obj.n_obs:
    raise ValueError(
        "Guide-assignment row count does not match the number of cells."
    )

# Optional consistency check against obs.
if "n_assigned_guides" in obj.obs:
    stored = pd.to_numeric(
        obj.obs["n_assigned_guides"],
        errors="coerce",
    ).to_numpy()
    comparable = np.isfinite(stored)
    n_disagree = int(
        np.sum(
            stored[comparable].astype(int)
            != n_assigned_guides[comparable]
        )
    )
    print(
        "Guide-count consistency check against obs['n_assigned_guides']:",
        f"{n_disagree} disagreements",
    )

print(f"{SAMPLE_NAME} cells:", int(ip_mask.sum()))
print("Guide-assignment key:", guide_key)

```

### Text output (stream)

```text
Loaded: /oak/stanford/groups/engreitz/Users/tri/10xPerturbSeq/analysis/250507_ePAH1_NovaSeq/CRISPR_count/IP_GV_unified_sceptre_cnmf_no_celltype/PerturbNMF/results/IP_GV_skcnmf_K30_dt2_0_n100/Inference/adata/cNMF_30_2_0.h5mu
Object: cNMF (4095, 30)
Guide-count consistency check against obs['n_assigned_guides']: 0 disagreements
IP_GFP cells: 2849
Guide-assignment key: guide_assignment

```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)

```

## Code cell 4

```python
# ---------------------------------------------------------------------
# Cell-level table
# ---------------------------------------------------------------------
cell_table = pd.DataFrame(
    {
        "cell_id": obj.obs_names.astype(str),
        "sample": sample_values.to_numpy(),
        "n_assigned_guides": n_assigned_guides,
    }
)

# Preserve useful existing columns when available.
for column in [
    "barcode",
    "n_assigned_targets",
    "assigned_guide_umis",
    "assigned_guide_ids",
    "assigned_targets",
    "total_counts",
    "n_genes_by_counts",
    "pct_counts_mt",
]:
    if column in obj.obs and column not in cell_table:
        cell_table[column] = obj.obs[column].to_numpy()

ip_cells = cell_table.loc[
    cell_table["sample"] == SAMPLE_NAME
].copy()

# ---------------------------------------------------------------------
# Histogram-ready table
# ---------------------------------------------------------------------
histogram = (
    ip_cells["n_assigned_guides"]
    .value_counts()
    .sort_index()
    .rename_axis("n_assigned_guides")
    .rename("n_cells")
    .reset_index()
)

# Include missing integer bins, if any.
all_bins = pd.DataFrame(
    {
        "n_assigned_guides": np.arange(
            int(histogram["n_assigned_guides"].min()),
            int(histogram["n_assigned_guides"].max()) + 1,
        )
    }
)
histogram = all_bins.merge(
    histogram,
    on="n_assigned_guides",
    how="left",
)
histogram["n_cells"] = (
    histogram["n_cells"].fillna(0).astype(int)
)
histogram["fraction_cells"] = (
    histogram["n_cells"] / len(ip_cells)
)
histogram["percent_cells"] = (
    100.0 * histogram["fraction_cells"]
)
histogram["cumulative_fraction_leq"] = (
    histogram["fraction_cells"].cumsum()
)
histogram["cumulative_percent_leq"] = (
    100.0 * histogram["cumulative_fraction_leq"]
)

# ---------------------------------------------------------------------
# Compact summary
# ---------------------------------------------------------------------
guide_counts = ip_cells["n_assigned_guides"].to_numpy()

summary = pd.DataFrame(
    [
        {
            "sample": SAMPLE_NAME,
            "n_cells": int(len(ip_cells)),
            "mean_guides_per_cell": float(np.mean(guide_counts)),
            "median_guides_per_cell": float(np.median(guide_counts)),
            "standard_deviation": float(
                np.std(guide_counts, ddof=1)
            ),
            "minimum": int(np.min(guide_counts)),
            "maximum": int(np.max(guide_counts)),
            "fraction_single_guide": float(
                np.mean(guide_counts == 1)
            ),
            "fraction_multiguide": float(
                np.mean(guide_counts > 1)
            ),
            "fraction_three_or_more": float(
                np.mean(guide_counts >= 3)
            ),
        }
    ]
)

display(histogram)
display(summary)

```

### Text output (text/plain)

```text
    n_assigned_guides  n_cells  fraction_cells  percent_cells  \
0                   1     1300        0.456300      45.630046   
1                   2      581        0.203931      20.393120   
2                   3      369        0.129519      12.951913   
3                   4      235        0.082485       8.248508   
4                   5      145        0.050895       5.089505   
5                   6      102        0.035802       3.580204   
6                   7       57        0.020007       2.000702   
7                   8       33        0.011583       1.158301   
8                   9       18        0.006318       0.631801   
9                  10        6        0.002106       0.210600   
10                 11        1        0.000351       0.035100   
11                 12        2        0.000702       0.070200   

    cumulative_fraction_leq  cumulative_percent_leq  
0                  0.456300               45.630046  
1                  0.660232               66.023166  
2                  0.789751               78.975079  
3                  0.872236               87.223587  
4                  0.923131               92.313092  
5                  0.958933               95.893296  
6                  0.978940               97.893998  
7                  0.990523               99.052299  
8                  0.996841               99.684100  
9                  0.998947               99.894700  
10                 0.999298               99.929800  
11                 1.000000              100.000000  
```

### Text output (text/plain)

```text
   sample  n_cells  mean_guides_per_cell  median_guides_per_cell  \
0  IP_GFP     2849              2.374868                     2.0   

   standard_deviation  minimum  maximum  fraction_single_guide  \
0            1.807603        1       12                 0.4563   

   fraction_multiguide  fraction_three_or_more  
0               0.5437                0.339768  
```

## Code cell 5

```python
# ---------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------
cell_path = (
    OUTPUT_DIR / "IP_GFP_assigned_guides_per_cell.csv.gz"
)
histogram_path = (
    OUTPUT_DIR / "IP_GFP_assigned_guides_per_cell_histogram.csv"
)
summary_path = (
    OUTPUT_DIR / "IP_GFP_assigned_guides_per_cell_summary.csv"
)
metadata_path = (
    OUTPUT_DIR / "IP_GFP_assigned_guides_per_cell_metadata.json"
)

ip_cells.to_csv(
    cell_path,
    index=False,
    compression="gzip",
)
histogram.to_csv(
    histogram_path,
    index=False,
)
summary.to_csv(
    summary_path,
    index=False,
)

metadata = {
    "source_path": str(source_path),
    "source_kind": source_kind,
    "modality": modality_name,
    "sample_column": sample_key,
    "sample_value": SAMPLE_NAME,
    "guide_assignment_key": guide_key,
    "guide_count_definition": (
        "number of nonzero SCEPTRE guide assignments per cell"
    ),
    "n_cells": int(len(ip_cells)),
    "includes_zero_guide_cells": bool(
        (ip_cells["n_assigned_guides"] == 0).any()
    ),
}
metadata_path.write_text(
    json.dumps(metadata, indent=2) + "\n"
)

print("Cell-level CSV:", cell_path)
print("Histogram CSV:", histogram_path)
print("Summary CSV:", summary_path)
print("Metadata:", metadata_path)

```

### Text output (stream)

```text
Cell-level CSV: /oak/stanford/groups/engreitz/Users/tri/10xPerturbSeq/analysis/250507_ePAH1_NovaSeq/CRISPR_count/IP_GV_unified_sceptre_cnmf_no_celltype/PerturbNMF/results/IP_GV_skcnmf_K30_dt2_0_n100/Evaluation/30_2_0/guide_count_histogram/IP_GFP_assigned_guides_per_cell.csv.gz
Histogram CSV: /oak/stanford/groups/engreitz/Users/tri/10xPerturbSeq/analysis/250507_ePAH1_NovaSeq/CRISPR_count/IP_GV_unified_sceptre_cnmf_no_celltype/PerturbNMF/results/IP_GV_skcnmf_K30_dt2_0_n100/Evaluation/30_2_0/guide_count_histogram/IP_GFP_assigned_guides_per_cell_histogram.csv
Summary CSV: /oak/stanford/groups/engreitz/Users/tri/10xPerturbSeq/analysis/250507_ePAH1_NovaSeq/CRISPR_count/IP_GV_unified_sceptre_cnmf_no_celltype/PerturbNMF/results/IP_GV_skcnmf_K30_dt2_0_n100/Evaluation/30_2_0/guide_count_histogram/IP_GFP_assigned_guides_per_cell_summary.csv
Metadata: /oak/stanford/groups/engreitz/Users/tri/10xPerturbSeq/analysis/250507_ePAH1_NovaSeq/CRISPR_count/IP_GV_unified_sceptre_cnmf_no_celltype/PerturbNMF/results/IP_GV_skcnmf_K30_dt2_0_n100/Evaluation/30_2_0/guide_count_histogram/IP_GFP_assigned_guides_per_cell_metadata.json

```

