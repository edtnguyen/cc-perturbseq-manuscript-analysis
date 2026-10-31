# Notebook mirror: `oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/PAH.PerturbSeq/IP_GV_K30_dt2_CRT_notebook_bundle/IP_GV_K30_dt2_CRT_calibration_and_discovery.ipynb`

- Source archive path: `/oak/stanford/groups/engreitz/Users/tri/CodeBase/cc-perturbseq-figure-reproduction/recovery/sherlock_notebooks/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/PAH.PerturbSeq/IP_GV_K30_dt2_CRT_notebook_bundle/IP_GV_K30_dt2_CRT_calibration_and_discovery.ipynb`
- SHA-256: `c09e054a3397568c210cf9c10504e574207f8f7bed7c30c39ea14ef98ad1b171`
- Kernel language: `python`
- Full-fidelity notebook storage: Git LFS

## Markdown cell 0

# IP/GV target–program CRT

**Selected cNMF model:** `K=30`, density threshold `2.0`, 100 cNMF replicates.

The pilot metrics do not provide a strong reason to reject density threshold 2.0 at K=30:

- total explained variance is effectively unchanged;
- the fraction of programs with significant GO terms is the same;
- Reactome coverage is slightly higher at threshold 2.0;
- the mouse-associated fraction is the same;
- cross-mouse perturbation concordance is modestly lower but remains in the same range.

This notebook first compares CRT calibration across covariate sets using pseudo-targets made from non-targeting guides. It then optionally runs the full target–program analysis.

Used the **Python (programDE IP-GV)** kernel created by `setup_programDE_ipgv_mamba.sh`.


## Code cell 1

```python
import os

for variable in [
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
]:
    os.environ.setdefault(variable, "1")

os.environ.setdefault("PYTHONNOUSERSITE", "1")

```

### Text output (text/plain)

```text
'1'
```

## Code cell 2

```python
from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
import inspect
import json
import math
import re
import subprocess
import sys
import time
from typing import Any

import anndata as ad
from IPython.display import display
import matplotlib.pyplot as plt
import muon as mu
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import kstest, norm
from statsmodels.stats.multitest import multipletests

from src.sceptre import (
    build_ntc_group_inputs,
    compute_guide_set_null_pvals,
    crt_pvals_for_ntc_groups_ensemble,
    crt_pvals_for_ntc_groups_ensemble_skew,
    make_ntc_groups_ensemble,
    prepare_crt_inputs,
    run_all_genes_union_crt,
)
from src.visualization import qq_plot_ntc_pvals

try:
    from src.sceptre import limit_threading
    limit_threading()
except (ImportError, AttributeError):
    pass

print("Python:", sys.version.split()[0])
print("prepare_crt_inputs:", inspect.signature(prepare_crt_inputs))
print("run_all_genes_union_crt:", inspect.signature(run_all_genes_union_crt))

```

### Text output (stream)

```text
Python: 3.11.14
prepare_crt_inputs: (adata: Any, usage_key: str = 'cnmf_usage', covar_key: str = 'covar', guide_assignment_key: str = 'guide_assignment', guide_names_key: str = 'guide_names', guide2gene_key: str = 'guide2gene', eps_quantile: float = 0.0001, add_intercept: bool = True, standardize: bool = True, numeric_as_category_threshold: Optional[int] = 20, batch_key: str = 'batch', clamp_threads: bool = True) -> src.sceptre.pipeline.CRTInputs
run_all_genes_union_crt: (inputs: src.sceptre.pipeline.CRTInputs, genes: Optional[Iterable[str]] = None, B: int = 1023, n_jobs: int = 8, base_seed: int = 123, propensity_model: Callable = <function fit_propensity_logistic at 0x7f94b12c6840>, backend: str = 'loky', resampling_method: str = 'bernoulli_index', resampling_kwargs: Optional[Dict[str, Any]] = None, null_method: str = 'crt', null_kwargs: Optional[Dict[str, Any]] = None, qq_crossfit: bool = False, calibrate_skew_normal: bool = False, skew_normal_side_code: int = 0, return_skew_normal: bool = False, return_raw_pvals: bool = False, return_format: str = 'dict', test_stat: str = 'ols', test_stat_kwargs: Optional[Dict[str, Any]] = None) -> Union[Tuple, Dict[str, Any]]

```

## Markdown cell 3

## Configuration

## Code cell 4

```python
# ---------------------------------------------------------------------
# Final selected cNMF model
# ---------------------------------------------------------------------
K = 30
DENSITY_THRESHOLD = 2.0

BASE = Path(
    "/oak/stanford/groups/engreitz/Users/tri/10xPerturbSeq/"
    "analysis/250507_ePAH1_NovaSeq/CRISPR_count/"
    "IP_GV_unified_sceptre_cnmf_no_celltype"
)

RESULTS_ROOT = BASE / "PerturbNMF" / "results"
EXPECTED_RUN_NAME = "IP_GV_skcnmf_K30_dt2_0_n100"
EXPECTED_H5MU = (
    RESULTS_ROOT
    / EXPECTED_RUN_NAME
    / "Inference"
    / "adata"
    / "cNMF_30_2_0.h5mu"
)

# Auto-discover only when the expected path is absent.
if EXPECTED_H5MU.exists():
    H5MU_PATH = EXPECTED_H5MU
else:
    candidates = sorted(
        RESULTS_ROOT.glob(
            "*K30*dt2_0*n100/Inference/adata/cNMF_30_2_0.h5mu"
        )
    )
    if len(candidates) != 1:
        raise FileNotFoundError(
            f"Expected {EXPECTED_H5MU}, and auto-discovery found "
            f"{len(candidates)} candidates: {candidates}"
        )
    H5MU_PATH = candidates[0]

RUN_ROOT = H5MU_PATH.parents[2]
CRT_ROOT = (
    RUN_ROOT
    / "Evaluation"
    / "30_2_0"
    / "CRT_notebook"
)
CALIBRATION_ROOT = CRT_ROOT / "calibration"
DISCOVERY_ROOT = CRT_ROOT / "discovery"

for directory in [CRT_ROOT, CALIBRATION_ROOT, DISCOVERY_ROOT]:
    directory.mkdir(parents=True, exist_ok=True)

CONTROL_TARGET_LABELS = ["non-targeting"]
CATEGORICAL_KEY = "sample"

# Calibration defaults
CALIBRATION_B = 1024
GUIDES_PER_PSEUDO_TARGET = 4
NTC_ENSEMBLES = 20
NTC_MAX_GROUPS = 7
NTC_FREQUENCY_BINS = 5
NTC_GROUP_SEED = 7
NTC_CRT_SEED = 23
CALIBRATION_COMPUTE_SKEW = True

# Full production run defaults
DISCOVERY_B = 5000
N_JOBS = min(16, os.cpu_count() or 1)

print("H5MU:", H5MU_PATH)
print("Run root:", RUN_ROOT)
print("CRT output:", CRT_ROOT)
print("n_jobs:", N_JOBS)

```

### Text output (stream)

```text
H5MU: /oak/stanford/groups/engreitz/Users/tri/10xPerturbSeq/analysis/250507_ePAH1_NovaSeq/CRISPR_count/IP_GV_unified_sceptre_cnmf_no_celltype/PerturbNMF/results/IP_GV_skcnmf_K30_dt2_0_n100/Inference/adata/cNMF_30_2_0.h5mu
Run root: /oak/stanford/groups/engreitz/Users/tri/10xPerturbSeq/analysis/250507_ePAH1_NovaSeq/CRISPR_count/IP_GV_unified_sceptre_cnmf_no_celltype/PerturbNMF/results/IP_GV_skcnmf_K30_dt2_0_n100
CRT output: /oak/stanford/groups/engreitz/Users/tri/10xPerturbSeq/analysis/250507_ePAH1_NovaSeq/CRISPR_count/IP_GV_unified_sceptre_cnmf_no_celltype/PerturbNMF/results/IP_GV_skcnmf_K30_dt2_0_n100/Evaluation/30_2_0/CRT_notebook
n_jobs: 16

```

## Markdown cell 5

## Load and validate the final H5MU

## Code cell 6

```python
mdata = mu.read(H5MU_PATH)

print(mdata)
print("Modalities:", list(mdata.mod.keys()))

if "rna" not in mdata.mod or "cNMF" not in mdata.mod:
    raise KeyError("Expected RNA and cNMF modalities.")

rna = mdata["rna"]
programs = mdata["cNMF"]

if programs.n_vars != K:
    raise ValueError(
        f"Expected K={K} programs, found {programs.n_vars}."
    )
if rna.n_obs != programs.n_obs:
    raise ValueError("RNA and cNMF modalities have different cell counts.")
if not rna.obs_names.equals(programs.obs_names):
    raise ValueError("RNA and cNMF cell orders differ.")

print("RNA shape:", rna.shape)
print("Program shape:", programs.shape)

```

### Text output (stream)

```text
MuData object with n_obs × n_vars = 4095 × 16896
  2 modalities
    rna:	4095 x 16866
      obs:	'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'total_counts_mt', 'log1p_total_counts_mt', 'pct_counts_mt', 'total_counts_ribo', 'log1p_total_counts_ribo', 'pct_counts_ribo', 'total_counts_hb', 'log1p_total_counts_hb', 'pct_counts_hb', 'sample_id', 'barcode', 'passes_rna_qc', 'raw_total_guide_umis', 'raw_n_detected_guides', 'total_guide_umis', 'n_detected_guides', 'n_assigned_guides', 'n_assigned_targets', 'assigned_guide_umis', 'assigned_guide_ids', 'assigned_targets', 'guide_assignment_status', 'has_sceptre_guide', 'sample', 'mouse'
      var:	'gene_ids', 'feature_types', 'genome', 'pattern', 'read', 'sequence', 'target_gene_id', 'target_gene_name', 'feature_name', 'symbol', 'mt', 'ribo', 'hb', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts', 'n_cells_combined', 'var_names'
      uns:	'guide_names', 'guide_targets', 'perturbnmf_preparation', 'target_names'
      obsm:	'X_pca', 'X_umap', 'assigned_guide_counts', 'guide_assignment', 'guide_assignments', 'guide_counts', 'target_assigned_guide_counts', 'target_assignments'
      layers:	'counts'
    cNMF:	4095 x 30
      obs:	'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'total_counts_mt', 'log1p_total_counts_mt', 'pct_counts_mt', 'total_counts_ribo', 'log1p_total_counts_ribo', 'pct_counts_ribo', 'total_counts_hb', 'log1p_total_counts_hb', 'pct_counts_hb', 'sample_id', 'barcode', 'passes_rna_qc', 'raw_total_guide_umis', 'raw_n_detected_guides', 'total_guide_umis', 'n_detected_guides', 'n_assigned_guides', 'n_assigned_targets', 'assigned_guide_umis', 'assigned_guide_ids', 'assigned_targets', 'guide_assignment_status', 'has_sceptre_guide', 'sample', 'mouse'
      uns:	'guide_names', 'guide_targets', 'var_names'
      obsm:	'X_pca', 'X_umap', 'guide_assignment'
      varm:	'loadings'
Modalities: ['rna', 'cNMF']
RNA shape: (4095, 16866)
Program shape: (4095, 30)

```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)

```

## Code cell 7

```python
def resolve_obs_series(key: str) -> pd.Series:
    """Resolve an obs column from cNMF first, then RNA."""
    if key in programs.obs:
        values = programs.obs[key]
    elif key in rna.obs:
        values = rna.obs.loc[programs.obs_names, key]
    else:
        raise KeyError(
            f"{key!r} is absent from both cNMF.obs and RNA.obs."
        )
    return pd.Series(
        values.to_numpy(),
        index=programs.obs_names,
        name=key,
    )


# Resolve guide metadata from cNMF first, with RNA fallback.
if "guide_assignment" in programs.obsm:
    guide_assignment = programs.obsm["guide_assignment"]
elif "guide_assignments" in programs.obsm:
    guide_assignment = programs.obsm["guide_assignments"]
elif "guide_assignment" in rna.obsm:
    guide_assignment = rna.obsm["guide_assignment"]
elif "guide_assignments" in rna.obsm:
    guide_assignment = rna.obsm["guide_assignments"]
else:
    raise KeyError("No guide-assignment matrix found.")

if sp.issparse(guide_assignment):
    guide_assignment = guide_assignment.tocsr().copy()
else:
    guide_assignment = sp.csr_matrix(np.asarray(guide_assignment))

guide_assignment.eliminate_zeros()
guide_assignment.data = np.ones(
    guide_assignment.nnz,
    dtype=np.int8,
)

guide_names = None
guide_targets = None
for modality in [programs, rna]:
    if guide_names is None and "guide_names" in modality.uns:
        guide_names = np.asarray(modality.uns["guide_names"]).astype(str)
    if guide_targets is None and "guide_targets" in modality.uns:
        guide_targets = np.asarray(
            modality.uns["guide_targets"]
        ).astype(str)

if guide_names is None or guide_targets is None:
    raise KeyError("guide_names or guide_targets are missing.")
if guide_assignment.shape != (programs.n_obs, len(guide_names)):
    raise ValueError(
        f"Guide assignment has shape {guide_assignment.shape}; "
        f"expected {(programs.n_obs, len(guide_names))}."
    )
if len(guide_targets) != len(guide_names):
    raise ValueError("guide_names and guide_targets lengths differ.")

usage = np.asarray(programs.X, dtype=float)
if usage.shape != (programs.n_obs, K):
    raise ValueError("Unexpected program-usage shape.")
if not np.isfinite(usage).all() or np.any(usage < 0):
    raise ValueError("Program usage contains invalid values.")

usage = np.maximum(usage, 1e-8)
usage /= usage.sum(axis=1, keepdims=True)

sample = resolve_obs_series(CATEGORICAL_KEY).astype(str)
sample_counts = sample.value_counts()

target_counts = (
    pd.Series(guide_targets, name="target")
    .value_counts()
    .rename("n_guides")
)

assigned_per_cell = np.asarray(
    guide_assignment.sum(axis=1)
).ravel()

print("Sample counts:")
display(sample_counts.to_frame())

print("Guide counts by target:")
display(target_counts.to_frame())

print(
    "Assigned guides/cell:",
    pd.Series(assigned_per_cell).describe().to_dict(),
)

if not set(CONTROL_TARGET_LABELS).intersection(set(guide_targets)):
    raise ValueError(
        f"Control labels {CONTROL_TARGET_LABELS} do not occur in guide_targets."
    )

```

### Text output (stream)

```text
Sample counts:

```

### Text output (text/plain)

```text
        count
sample       
IP_GFP   2849
GV_GFP   1246
```

### Text output (stream)

```text
Guide counts by target:

```

### Text output (text/plain)

```text
               n_guides
target                 
non-targeting        30
Cdkn1a                5
Notch1                5
Acvrl1                4
Bmpr2                 4
Cav1                  4
Eif2ak4               4
Eng                   4
Atp13a3               4
Gdf2                  4
Kcnk3                 4
Smad9                 4
Kdr                   4
Sox17                 4
Tbx4                  4
Smad4                 4
Trp53                 4
```

### Text output (stream)

```text
Assigned guides/cell: {'count': 4095.0, 'mean': 2.2515262515262515, 'std': 1.7212235880871438, 'min': 1.0, '25%': 1.0, '50%': 2.0, '75%': 3.0, 'max': 12.0}

```

## Code cell 8

```python
available_obs = sorted(
    set(programs.obs.columns).union(rna.obs.columns)
)
print("Available obs columns:")
print(available_obs)

candidate_keys = [
    "total_counts",
    "n_genes_by_counts",
    "pct_counts_mt",
    "pct_counts_ribo",
    "pct_counts_hb",
    "n_assigned_guides",
    "total_guide_umis",
    "assigned_guide_umis",
]

summary_rows = []
for key in candidate_keys:
    if key in available_obs:
        values = pd.to_numeric(
            resolve_obs_series(key),
            errors="coerce",
        )
        summary_rows.append(
            {
                "key": key,
                "missing": int(values.isna().sum()),
                "n_unique": int(values.nunique(dropna=True)),
                "min": float(values.min()),
                "median": float(values.median()),
                "max": float(values.max()),
            }
        )

display(pd.DataFrame(summary_rows))

```

### Text output (stream)

```text
Available obs columns:
['assigned_guide_ids', 'assigned_guide_umis', 'assigned_targets', 'barcode', 'guide_assignment_status', 'has_sceptre_guide', 'log1p_n_genes_by_counts', 'log1p_total_counts', 'log1p_total_counts_hb', 'log1p_total_counts_mt', 'log1p_total_counts_ribo', 'mouse', 'n_assigned_guides', 'n_assigned_targets', 'n_detected_guides', 'n_genes_by_counts', 'passes_rna_qc', 'pct_counts_hb', 'pct_counts_mt', 'pct_counts_ribo', 'raw_n_detected_guides', 'raw_total_guide_umis', 'sample', 'sample_id', 'total_counts', 'total_counts_hb', 'total_counts_mt', 'total_counts_ribo', 'total_guide_umis']

```

### Text output (text/plain)

```text
                   key  missing  n_unique         min       median  \
0         total_counts        0      3634  533.000000  8339.000000   
1    n_genes_by_counts        0      2667  420.000000  3226.000000   
2        pct_counts_mt        0      4085    0.000000     3.077370   
3      pct_counts_ribo        0      4091    0.986193     6.745182   
4        pct_counts_hb        0      2266    0.000000     0.007552   
5    n_assigned_guides        0        12    1.000000     2.000000   
6     total_guide_umis        0      1251    1.000000    63.000000   
7  assigned_guide_umis        0      1231    1.000000    60.000000   

            max  
0  47697.000000  
1   7976.000000  
2     14.712064  
3     29.503580  
4     25.227612  
5     12.000000  
6   7962.000000  
7   7811.000000  
```

## Markdown cell 9

## Candidate covariate sets

The CRT is run separately within each mouse, so `sample` is not included as a within-mouse covariate.

The initial variants are intentionally limited:

- **none**: intercept only;
- **rna_qc**: RNA depth, detected genes, mitochondrial and ribosomal percentages;
- **rna_qc_guides**: RNA QC plus assigned-guide number and guide UMI burden;
- **rna_qc_guides_hb**: adds hemoglobin percentage when available.

Avoid including redundant columns such as both total and assigned guide UMIs in the same model. Prefer the simplest set that gives acceptable calibration.


## Code cell 10

```python
def existing(keys: list[str]) -> list[str]:
    return [key for key in keys if key in available_obs]


COVARIATE_VARIANTS = OrderedDict(
    {
        "none": {
            "covariates": [],
            "log_covariates": [],
        },
        "rna_qc": {
            "covariates": existing(
                ["pct_counts_mt", "pct_counts_ribo"]
            ),
            "log_covariates": existing(
                ["total_counts", "n_genes_by_counts"]
            ),
        },
        "rna_qc_guides": {
            "covariates": existing(
                ["pct_counts_mt", "pct_counts_ribo"]
            ),
            "log_covariates": existing(
                [
                    "total_counts",
                    "n_genes_by_counts",
                    "n_assigned_guides",
                    "total_guide_umis",
                ]
            ),
        },
        "rna_qc_guides_hb": {
            "covariates": existing(
                [
                    "pct_counts_mt",
                    "pct_counts_ribo",
                    "pct_counts_hb",
                ]
            ),
            "log_covariates": existing(
                [
                    "total_counts",
                    "n_genes_by_counts",
                    "n_assigned_guides",
                    "total_guide_umis",
                ]
            ),
        },
    }
)

display(
    pd.DataFrame(
        [
            {
                "variant": name,
                "covariates": ", ".join(spec["covariates"]),
                "log_covariates": ", ".join(
                    spec["log_covariates"]
                ),
            }
            for name, spec in COVARIATE_VARIANTS.items()
        ]
    )
)

```

### Text output (text/plain)

```text
            variant                                     covariates  \
0              none                                                  
1            rna_qc                 pct_counts_mt, pct_counts_ribo   
2     rna_qc_guides                 pct_counts_mt, pct_counts_ribo   
3  rna_qc_guides_hb  pct_counts_mt, pct_counts_ribo, pct_counts_hb   

                                      log_covariates  
0                                                     
1                    total_counts, n_genes_by_counts  
2  total_counts, n_genes_by_counts, n_assigned_gu...  
3  total_counts, n_genes_by_counts, n_assigned_gu...  
```

## Markdown cell 11

## CRT helper functions

## Code cell 12

```python
def safe_name(value: object) -> str:
    return (
        re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_")
        or "condition"
    )


def build_crt_adata(
    covariates: list[str],
    log_covariates: list[str],
) -> ad.AnnData:
    adata = programs.copy()
    adata.obsm["cnmf_usage"] = usage.copy()
    adata.obsm["guide_assignment"] = guide_assignment.copy()
    adata.uns["program_names"] = list(adata.var_names)
    adata.uns["guide_names"] = guide_names.tolist()
    adata.uns["guide2gene"] = dict(
        zip(guide_names.tolist(), guide_targets.tolist())
    )

    if CATEGORICAL_KEY not in adata.obs:
        adata.obs[CATEGORICAL_KEY] = resolve_obs_series(
            CATEGORICAL_KEY
        ).to_numpy()

    adata.obs["guides_per_cell"] = assigned_per_cell

    covar = pd.DataFrame(index=adata.obs_names)

    for key in covariates:
        values = pd.to_numeric(
            resolve_obs_series(key),
            errors="coerce",
        )
        if values.isna().any():
            raise ValueError(
                f"Covariate {key!r} has missing/non-numeric values."
            )
        covar[key] = values

    for key in log_covariates:
        values = pd.to_numeric(
            resolve_obs_series(key),
            errors="coerce",
        )
        if values.isna().any():
            raise ValueError(
                f"Log covariate {key!r} has missing/non-numeric values."
            )
        if (values < 0).any():
            raise ValueError(
                f"Log covariate {key!r} contains negative values."
            )
        covar[f"log_{key}"] = np.log1p(values)

    # Remove globally constant covariates. A second condition-specific
    # check is performed below.
    keep = [
        column
        for column in covar.columns
        if covar[column].nunique(dropna=False) > 1
    ]
    adata.obsm["covar"] = covar[keep]
    return adata


def condition_subset(
    adata: ad.AnnData,
    condition: str,
) -> ad.AnnData:
    mask = (
        adata.obs[CATEGORICAL_KEY].astype(str)
        == str(condition)
    )
    subset = adata[mask].copy()

    covar = subset.obsm["covar"]
    if not isinstance(covar, pd.DataFrame):
        covar = pd.DataFrame(
            np.asarray(covar),
            index=subset.obs_names,
        )

    keep = [
        column
        for column in covar.columns
        if covar[column].nunique(dropna=False) > 1
    ]
    subset.obsm["covar"] = covar[keep]
    return subset


def prepare_inputs(adata: ad.AnnData):
    kwargs = dict(
        adata=adata,
        usage_key="cnmf_usage",
        covar_key="covar",
        guide_assignment_key="guide_assignment",
        guide2gene_key="guide2gene",
    )
    signature = inspect.signature(prepare_crt_inputs)

    # Treat numeric covariates as continuous even when they have <=20
    # distinct values, such as guides per cell.
    if "numeric_as_category_threshold" in signature.parameters:
        kwargs["numeric_as_category_threshold"] = None

    return prepare_crt_inputs(**kwargs)


def flatten_ntc_frames(frames: dict) -> np.ndarray:
    arrays = []
    for frame in frames.values():
        values = frame.to_numpy(dtype=float, copy=False).ravel()
        arrays.append(values[np.isfinite(values)])
    if not arrays:
        return np.empty(0, dtype=float)
    return np.concatenate(arrays)


def qq_xy(pvalues: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    values = np.asarray(pvalues, dtype=float)
    values = values[np.isfinite(values)]
    values = np.clip(values, 1e-300, 1.0)
    values.sort()
    n = values.size
    if n == 0:
        return np.empty(0), np.empty(0)
    expected = (np.arange(1, n + 1) - 0.5) / n
    return -np.log10(expected), -np.log10(values)


def calibration_metrics(pvalues: np.ndarray) -> dict[str, float]:
    values = np.asarray(pvalues, dtype=float)
    values = values[np.isfinite(values)]

    if values.size == 0:
        return {"n_pvalues": 0}

    ks = kstest(values, "uniform")
    result = {
        "n_pvalues": int(values.size),
        "mean_p": float(values.mean()),
        "median_p": float(np.median(values)),
        "ks_stat": float(ks.statistic),
        "ks_pvalue": float(ks.pvalue),
    }

    for threshold in [0.05, 0.01, 0.001]:
        fraction = float(np.mean(values <= threshold))
        result[f"fraction_p_le_{threshold:g}"] = fraction
        result[f"inflation_ratio_{threshold:g}"] = (
            fraction / threshold
        )

    return result


def save_ntc_frames(
    frames: dict,
    path: Path,
) -> None:
    rows = []
    for ensemble, frame in frames.items():
        current = frame.copy()
        current.index.name = "group_id"
        current = current.reset_index()
        current.insert(0, "ensemble", ensemble)
        rows.append(current)

    if rows:
        pd.concat(rows, ignore_index=True).to_csv(
            path,
            sep="\t",
            index=False,
        )


@dataclass
class CalibrationResult:
    variant: str
    condition: str
    raw_pvalues: np.ndarray
    skew_pvalues: np.ndarray
    crt_null_pvalues: np.ndarray
    output_dir: Path


def plot_calibration_result(
    result: CalibrationResult,
    title: str | None = None,
):
    figure, axis = plt.subplots(figsize=(7, 6))

    plotted = []
    for label, values, style in [
        ("NTC pseudo-target empirical", result.raw_pvalues, "-"),
        ("NTC pseudo-target skew", result.skew_pvalues, "--"),
        ("CRT-null leave-one-out", result.crt_null_pvalues, ":"),
    ]:
        x, y = qq_xy(values)
        if x.size:
            axis.plot(x, y, style, label=label)
            plotted.extend([x.max(), y.max()])

    limit = max([1.0] + [float(value) for value in plotted])
    axis.plot([0, limit], [0, limit], label="y=x")
    axis.set_xlim(0, limit)
    axis.set_ylim(0, limit)
    axis.set_xlabel("Expected -log10(p)")
    axis.set_ylabel("Observed -log10(p)")
    axis.set_title(
        title
        or f"{result.variant}: {result.condition}"
    )
    axis.legend(frameon=False)
    figure.tight_layout()
    return figure, axis

```

## Code cell 13

```python
def run_ntc_calibration(
    variant_name: str,
    condition: str,
    *,
    B: int = CALIBRATION_B,
    compute_skew: bool = CALIBRATION_COMPUTE_SKEW,
    overwrite: bool = False,
) -> CalibrationResult:
    specification = COVARIATE_VARIANTS[variant_name]
    output_dir = (
        CALIBRATION_ROOT
        / variant_name
        / safe_name(condition)
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    raw_path = output_dir / "ntc_empirical_pvalues.npy"
    skew_path = output_dir / "ntc_skew_pvalues.npy"
    null_path = output_dir / "crt_null_pvalues.npy"

    if (
        not overwrite
        and raw_path.exists()
        and null_path.exists()
        and (not compute_skew or skew_path.exists())
    ):
        return CalibrationResult(
            variant=variant_name,
            condition=condition,
            raw_pvalues=np.load(raw_path),
            skew_pvalues=(
                np.load(skew_path)
                if skew_path.exists()
                else np.empty(0)
            ),
            crt_null_pvalues=np.load(null_path),
            output_dir=output_dir,
        )

    adata = build_crt_adata(
        specification["covariates"],
        specification["log_covariates"],
    )
    adata_condition = condition_subset(adata, condition)
    print(
        f"{variant_name} / {condition}: "
        f"{adata_condition.n_obs:,} cells, "
        f"{adata_condition.obsm['covar'].shape[1]} covariates"
    )

    inputs = prepare_inputs(adata_condition)

    (
        ntc_guides,
        guide_frequency,
        guide_to_bin,
        real_gene_bin_signatures,
    ) = build_ntc_group_inputs(
        inputs=inputs,
        ntc_label=CONTROL_TARGET_LABELS,
        group_size=GUIDES_PER_PSEUDO_TARGET,
        n_bins=NTC_FREQUENCY_BINS,
    )

    groups = make_ntc_groups_ensemble(
        ntc_guides=ntc_guides,
        ntc_freq=guide_frequency,
        real_gene_bin_sigs=real_gene_bin_signatures,
        guide_to_bin=guide_to_bin,
        n_ensemble=NTC_ENSEMBLES,
        seed0=NTC_GROUP_SEED,
        group_size=GUIDES_PER_PSEUDO_TARGET,
        max_groups=NTC_MAX_GROUPS,
    )

    empirical_frames = crt_pvals_for_ntc_groups_ensemble(
        inputs=inputs,
        ntc_groups_ens=groups,
        B=B,
        seed0=NTC_CRT_SEED,
    )
    empirical = flatten_ntc_frames(empirical_frames)

    if compute_skew:
        skew_frames = crt_pvals_for_ntc_groups_ensemble_skew(
            inputs=inputs,
            ntc_groups_ens=groups,
            B=B,
            seed0=NTC_CRT_SEED,
        )
        skew = flatten_ntc_frames(skew_frames)
    else:
        skew_frames = {}
        skew = np.empty(0, dtype=float)

    guide_to_column = {
        guide: index
        for index, guide in enumerate(inputs.guide_names)
    }
    null_chunks = []
    for ensemble_groups in groups:
        for guides in ensemble_groups.values():
            null_chunks.append(
                compute_guide_set_null_pvals(
                    guide_idx=[
                        guide_to_column[guide]
                        for guide in guides
                    ],
                    inputs=inputs,
                    B=B,
                ).ravel()
            )
    null_pvalues = (
        np.concatenate(null_chunks)
        if null_chunks
        else np.empty(0, dtype=float)
    )

    np.save(raw_path, empirical)
    np.save(skew_path, skew)
    np.save(null_path, null_pvalues)

    save_ntc_frames(
        empirical_frames,
        output_dir / "ntc_empirical_pvalues_wide.tsv",
    )
    if skew_frames:
        save_ntc_frames(
            skew_frames,
            output_dir / "ntc_skew_pvalues_wide.tsv",
        )

    config = {
        "variant": variant_name,
        "condition": condition,
        "B": B,
        "covariates": specification["covariates"],
        "log_covariates": specification["log_covariates"],
        "guides_per_pseudo_target": GUIDES_PER_PSEUDO_TARGET,
        "ntc_ensembles": NTC_ENSEMBLES,
        "ntc_max_groups": NTC_MAX_GROUPS,
        "ntc_frequency_bins": NTC_FREQUENCY_BINS,
        "ntc_group_seed": NTC_GROUP_SEED,
        "ntc_crt_seed": NTC_CRT_SEED,
    }
    (output_dir / "config.json").write_text(
        json.dumps(config, indent=2) + "\n"
    )

    result = CalibrationResult(
        variant=variant_name,
        condition=condition,
        raw_pvalues=empirical,
        skew_pvalues=skew,
        crt_null_pvalues=null_pvalues,
        output_dir=output_dir,
    )
    figure, _ = plot_calibration_result(result)
    figure.savefig(
        output_dir / "calibration_qq.png",
        dpi=180,
        bbox_inches="tight",
    )
    plt.close(figure)
    return result

```

## Markdown cell 14

## Run the calibration grid

Set `RUN_CALIBRATION_GRID = True` and run this cell. Results are cached on disk, so rerunning the cell does not recompute completed variants unless `OVERWRITE_CALIBRATION = True`.


## Code cell 15

```python
RUN_CALIBRATION_GRID = True
OVERWRITE_CALIBRATION = False

conditions = sample_counts.index.astype(str).tolist()
calibration_results: dict[tuple[str, str], CalibrationResult] = {}

if RUN_CALIBRATION_GRID:
    started = time.time()

    for variant_name in COVARIATE_VARIANTS:
        for condition in conditions:
            result = run_ntc_calibration(
                variant_name,
                condition,
                overwrite=OVERWRITE_CALIBRATION,
            )
            calibration_results[
                (variant_name, condition)
            ] = result

            figure, _ = plot_calibration_result(result)
            display(figure)
            plt.close(figure)

    print(
        "Elapsed minutes:",
        round((time.time() - started) / 60, 2),
    )
else:
    print(
        "Calibration is disabled. Set RUN_CALIBRATION_GRID=True "
        "when ready."
    )

```

### Text output (stream)

```text
none / IP_GFP: 2,849 cells, 0 covariates

```

### Text output (text/plain)

```text
<Figure size 700x600 with 1 Axes>
```

### Text output (stream)

```text
none / GV_GFP: 1,246 cells, 0 covariates

```

### Text output (text/plain)

```text
<Figure size 700x600 with 1 Axes>
```

### Text output (stream)

```text
rna_qc / IP_GFP: 2,849 cells, 4 covariates

```

### Text output (text/plain)

```text
<Figure size 700x600 with 1 Axes>
```

### Text output (stream)

```text
rna_qc / GV_GFP: 1,246 cells, 4 covariates

```

### Text output (text/plain)

```text
<Figure size 700x600 with 1 Axes>
```

### Text output (stream)

```text
rna_qc_guides / IP_GFP: 2,849 cells, 6 covariates

```

### Text output (text/plain)

```text
<Figure size 700x600 with 1 Axes>
```

### Text output (stream)

```text
rna_qc_guides / GV_GFP: 1,246 cells, 6 covariates

```

### Text output (text/plain)

```text
<Figure size 700x600 with 1 Axes>
```

### Text output (stream)

```text
rna_qc_guides_hb / IP_GFP: 2,849 cells, 7 covariates

```

### Text output (text/plain)

```text
<Figure size 700x600 with 1 Axes>
```

### Text output (stream)

```text
rna_qc_guides_hb / GV_GFP: 1,246 cells, 7 covariates

```

### Text output (text/plain)

```text
<Figure size 700x600 with 1 Axes>
```

### Text output (stream)

```text
Elapsed minutes: 6.31

```

## Markdown cell 16

## Summarize and compare calibration

## Code cell 17

```python
# Load all completed calibration results from disk.
calibration_results = {}

for variant_name in COVARIATE_VARIANTS:
    for condition in conditions:
        directory = (
            CALIBRATION_ROOT
            / variant_name
            / safe_name(condition)
        )
        raw_path = directory / "ntc_empirical_pvalues.npy"
        skew_path = directory / "ntc_skew_pvalues.npy"
        null_path = directory / "crt_null_pvalues.npy"

        if raw_path.exists() and null_path.exists():
            calibration_results[
                (variant_name, condition)
            ] = CalibrationResult(
                variant=variant_name,
                condition=condition,
                raw_pvalues=np.load(raw_path),
                skew_pvalues=(
                    np.load(skew_path)
                    if skew_path.exists()
                    else np.empty(0)
                ),
                crt_null_pvalues=np.load(null_path),
                output_dir=directory,
            )

rows = []
for (variant_name, condition), result in calibration_results.items():
    for pvalue_type, values in [
        ("empirical", result.raw_pvalues),
        ("skew", result.skew_pvalues),
        ("crt_null", result.crt_null_pvalues),
    ]:
        if values.size == 0:
            continue
        row = {
            "variant": variant_name,
            "condition": condition,
            "pvalue_type": pvalue_type,
        }
        row.update(calibration_metrics(values))
        rows.append(row)

calibration_summary = pd.DataFrame(rows)

if calibration_summary.empty:
    print(
        "No calibration results found. Run the calibration grid first."
    )
else:
    # Distance is only a sorting aid, not an automatic selection rule.
    calibration_summary["calibration_distance"] = (
        calibration_summary["ks_stat"]
        + (
            calibration_summary["inflation_ratio_0.05"] - 1
        ).abs()
        + (
            calibration_summary["inflation_ratio_0.01"] - 1
        ).abs()
    )

    calibration_summary.to_csv(
        CALIBRATION_ROOT / "calibration_summary.tsv",
        sep="\t",
        index=False,
    )

    display(
        calibration_summary.sort_values(
            [
                "pvalue_type",
                "condition",
                "calibration_distance",
            ]
        )
    )

```

### Text output (text/plain)

```text
             variant condition pvalue_type  n_pvalues    mean_p  median_p  \
5               none    GV_GFP    crt_null    3686400  0.500488  0.500488   
11            rna_qc    GV_GFP    crt_null    3686400  0.500488  0.500488   
17     rna_qc_guides    GV_GFP    crt_null    3686400  0.500488  0.500488   
23  rna_qc_guides_hb    GV_GFP    crt_null    3686400  0.500488  0.500488   
2               none    IP_GFP    crt_null    4147200  0.500488  0.500488   
8             rna_qc    IP_GFP    crt_null    4147200  0.500488  0.500488   
14     rna_qc_guides    IP_GFP    crt_null    4147200  0.500488  0.500488   
20  rna_qc_guides_hb    IP_GFP    crt_null    4147200  0.500488  0.500488   
15     rna_qc_guides    GV_GFP   empirical       3600  0.499764  0.497561   
21  rna_qc_guides_hb    GV_GFP   empirical       3600  0.499658  0.498537   
9             rna_qc    GV_GFP   empirical       3600  0.461455  0.447805   
3               none    GV_GFP   empirical       3600  0.461142  0.450732   
12     rna_qc_guides    IP_GFP   empirical       4050  0.513006  0.516098   
18  rna_qc_guides_hb    IP_GFP   empirical       4050  0.513004  0.514146   
0               none    IP_GFP   empirical       4050  0.426242  0.408780   
6             rna_qc    IP_GFP   empirical       4050  0.404490  0.371707   
22  rna_qc_guides_hb    GV_GFP        skew       3600  0.494178  0.488024   
16     rna_qc_guides    GV_GFP        skew       3600  0.494344  0.487297   
10            rna_qc    GV_GFP        skew       3600  0.457792  0.442107   
4               none    GV_GFP        skew       3600  0.459601  0.448848   
13     rna_qc_guides    IP_GFP        skew       4050  0.514084  0.517813   
19  rna_qc_guides_hb    IP_GFP        skew       4050  0.514188  0.514214   
1               none    IP_GFP        skew       4050  0.425145  0.403889   
7             rna_qc    IP_GFP        skew       4050  0.403693  0.371063   

     ks_stat     ks_pvalue  fraction_p_le_0.05  inflation_ratio_0.05  \
5   0.000977  1.766502e-03            0.049805              0.996094   
11  0.000977  1.766502e-03            0.049805              0.996094   
17  0.000977  1.766502e-03            0.049805              0.996094   
23  0.000977  1.766502e-03            0.049805              0.996094   
2   0.000977  7.335168e-04            0.049805              0.996094   
8   0.000977  7.335168e-04            0.049805              0.996094   
14  0.000977  7.335168e-04            0.049805              0.996094   
20  0.000977  7.335168e-04            0.049805              0.996094   
15  0.010467  8.210704e-01            0.046111              0.922222   
21  0.008598  9.507904e-01            0.044444              0.888889   
9   0.065332  8.388709e-14            0.099167              1.983333   
3   0.069573  1.348477e-15            0.109722              2.194444   
12  0.021722  4.312660e-02            0.045926              0.918519   
18  0.022228  3.600608e-02            0.046420              0.928395   
0   0.116007  6.110772e-48            0.150370              3.007407   
6   0.148720  1.173856e-78            0.169877              3.397531   
22  0.017328  2.272117e-01            0.046667              0.933333   
16  0.016125  3.032623e-01            0.047778              0.955556   
10  0.071910  1.239538e-16            0.100833              2.016667   
4   0.071108  2.836193e-16            0.110833              2.216667   
13  0.024137  1.755960e-02            0.044938              0.898765   
19  0.025681  9.406686e-03            0.045926              0.918519   
1   0.115963  6.638570e-48            0.152099              3.041975   
7   0.147269  3.945623e-77            0.172346              3.446914   

    fraction_p_le_0.01  inflation_ratio_0.01  fraction_p_le_0.001  \
5             0.009766              0.976562             0.000977   
11            0.009766              0.976562             0.000977   
17            0.009766              0.976562             0.000977   
23            0.009766              0.976562             0.000977   
2             0.009766              0.976562             0.000977   
8             0.009766              0.976562             0.000977   
14            0.009766              0.976562             0.000977   
20            0.009766              0.976562             0.000977   
15            0.007778              0.777778             0.000556   
21            0.007778              0.777778             0.000278   
9             0.038333              3.833333             0.010833   
3             0.051944              5.194444             0.017778   
12            0.010123              1.012346             0.000494   
18            0.009383              0.938272             0.000247   
0             0.086914              8.691358             0.048889   
6             0.087654              8.765432             0.038519   
22            0.010833              1.083333             0.000833   
16            0.011667              1.166667             0.000556   
10            0.045556              4.555556             0.016944   
4             0.054722              5.472222             0.027222   
13            0.010123              1.012346             0.000741   
19            0.010370              1.037037             0.000494   
1             0.084444              8.444444             0.050864   
7             0.087901              8.790123             0.041481   

    inflation_ratio_0.001  calibration_distance  
5                0.976562              0.028320  
11               0.976562              0.028320  
17               0.976562              0.028320  
23               0.976562              0.028320  
2                0.976562              0.028320  
8                0.976562              0.028320  
14               0.976562              0.028320  
20               0.976562              0.028320  
15               0.555556              0.310467  
21               0.277778              0.341931  
9               10.833333              3.881999  
3               17.777778              5.458462  
12               0.493827              0.115550  
18               0.246914              0.155562  
0               48.888889              9.814773  
6               38.518519             10.311683  
22               0.833333              0.167328  
16               0.555556              0.227236  
10              16.944444              4.644132  
4               27.222222              5.759997  
13               0.740741              0.137717  
19               0.493827              0.144199  
1               50.864198              9.602383  
7               41.481481             10.384306  
```

## Code cell 18

```python
# Overlay candidate variants for each mouse and p-value family.
if not calibration_results:
    print("No calibration results to plot.")
else:
    for condition in conditions:
        for pvalue_type in ["empirical", "skew"]:
            figure, axis = plt.subplots(figsize=(7, 6))
            maxima = [1.0]

            for variant_name in COVARIATE_VARIANTS:
                key = (variant_name, condition)
                if key not in calibration_results:
                    continue

                result = calibration_results[key]
                values = (
                    result.raw_pvalues
                    if pvalue_type == "empirical"
                    else result.skew_pvalues
                )
                x, y = qq_xy(values)
                if x.size == 0:
                    continue

                axis.plot(x, y, label=variant_name)
                maxima.extend([float(x.max()), float(y.max())])

            limit = max(maxima)
            axis.plot([0, limit], [0, limit], label="y=x")
            axis.set_xlim(0, limit)
            axis.set_ylim(0, limit)
            axis.set_xlabel("Expected -log10(p)")
            axis.set_ylabel("Observed -log10(p)")
            axis.set_title(
                f"{condition}: {pvalue_type} NTC calibration"
            )
            axis.legend(frameon=False)
            figure.tight_layout()

            figure.savefig(
                CALIBRATION_ROOT
                / f"compare_{safe_name(condition)}_{pvalue_type}.png",
                dpi=180,
                bbox_inches="tight",
            )
            display(figure)
            plt.close(figure)

```

### Text output (text/plain)

```text
<Figure size 700x600 with 1 Axes>
```

### Text output (text/plain)

```text
<Figure size 700x600 with 1 Axes>
```

### Text output (text/plain)

```text
<Figure size 700x600 with 1 Axes>
```

### Text output (text/plain)

```text
<Figure size 700x600 with 1 Axes>
```

## Markdown cell 19

## Choose the final calibration specification

Guidelines:

- `fraction_p_le_0.05` should be near 0.05;
- `fraction_p_le_0.01` should be near 0.01;
- inflation ratios substantially above 1 indicate anti-conservative p-values;
- values substantially below 1 indicate conservative/deflated p-values;
- prefer a simpler covariate set when calibration is similar;
- verify calibration in **both mice**.

The default below is a starting point, not an automatic conclusion.


## Code cell 20

```python
SELECTED_VARIANT = "rna_qc_guides"
PRIMARY_PVALUE = "skew"  # "skew" or "empirical"

if SELECTED_VARIANT not in COVARIATE_VARIANTS:
    raise KeyError(SELECTED_VARIANT)
if PRIMARY_PVALUE not in {"skew", "empirical"}:
    raise ValueError(PRIMARY_PVALUE)

print("Selected variant:", SELECTED_VARIANT)
print("Specification:", COVARIATE_VARIANTS[SELECTED_VARIANT])
print("Primary p-value:", PRIMARY_PVALUE)

if not calibration_summary.empty:
    display(
        calibration_summary.loc[
            calibration_summary["variant"]
            == SELECTED_VARIANT
        ].sort_values(["condition", "pvalue_type"])
    )

```

### Text output (stream)

```text
Selected variant: rna_qc_guides
Specification: {'covariates': ['pct_counts_mt', 'pct_counts_ribo'], 'log_covariates': ['total_counts', 'n_genes_by_counts', 'n_assigned_guides', 'total_guide_umis']}
Primary p-value: skew

```

### Text output (text/plain)

```text
          variant condition pvalue_type  n_pvalues    mean_p  median_p  \
17  rna_qc_guides    GV_GFP    crt_null    3686400  0.500488  0.500488   
15  rna_qc_guides    GV_GFP   empirical       3600  0.499764  0.497561   
16  rna_qc_guides    GV_GFP        skew       3600  0.494344  0.487297   
14  rna_qc_guides    IP_GFP    crt_null    4147200  0.500488  0.500488   
12  rna_qc_guides    IP_GFP   empirical       4050  0.513006  0.516098   
13  rna_qc_guides    IP_GFP        skew       4050  0.514084  0.517813   

     ks_stat  ks_pvalue  fraction_p_le_0.05  inflation_ratio_0.05  \
17  0.000977   0.001767            0.049805              0.996094   
15  0.010467   0.821070            0.046111              0.922222   
16  0.016125   0.303262            0.047778              0.955556   
14  0.000977   0.000734            0.049805              0.996094   
12  0.021722   0.043127            0.045926              0.918519   
13  0.024137   0.017560            0.044938              0.898765   

    fraction_p_le_0.01  inflation_ratio_0.01  fraction_p_le_0.001  \
17            0.009766              0.976562             0.000977   
15            0.007778              0.777778             0.000556   
16            0.011667              1.166667             0.000556   
14            0.009766              0.976562             0.000977   
12            0.010123              1.012346             0.000494   
13            0.010123              1.012346             0.000741   

    inflation_ratio_0.001  calibration_distance  
17               0.976562              0.028320  
15               0.555556              0.310467  
16               0.555556              0.227236  
14               0.976562              0.028320  
12               0.493827              0.115550  
13               0.740741              0.137717  
```

## Markdown cell 21

## Full target–program CRT

## Code cell 22

```python
def long_results_table(
    out: dict[str, Any],
    *,
    condition: str,
    primary_pvalue: str,
) -> pd.DataFrame:
    empirical = out["pvals_raw_df"]

    skew = out.get("pvals_df")
    if skew is None:
        skew = out.get("pvals_skew_df")

    beta = out["betas_df"]

    def melt(frame: pd.DataFrame, value_name: str):
        index_name = frame.index.name or "target_name"
        current = frame.copy()
        current.index.name = "target_name"
        return (
            current.reset_index()
            .melt(
                id_vars="target_name",
                var_name="program_name",
                value_name=value_name,
            )
        )

    table = melt(empirical, "p_empirical").merge(
        melt(beta, "beta_clr"),
        on=["target_name", "program_name"],
        how="inner",
    )

    if skew is not None:
        table = table.merge(
            melt(skew, "p_skew"),
            on=["target_name", "program_name"],
            how="left",
        )
    else:
        table["p_skew"] = np.nan

    table["q_empirical_BH"] = multipletests(
        table["p_empirical"].fillna(1.0),
        method="fdr_bh",
    )[1]

    if table["p_skew"].notna().any():
        table["q_skew_BH"] = multipletests(
            table["p_skew"].fillna(1.0),
            method="fdr_bh",
        )[1]
    else:
        table["q_skew_BH"] = np.nan

    # beta_clr is a CLR log-ratio effect, not an ordinary expression log2FC.
    table["approx_relative_log2fc"] = (
        (K / (K - 1))
        * table["beta_clr"]
        / np.log(2)
    )

    table["condition"] = condition
    table["primary_pvalue_type"] = primary_pvalue

    if primary_pvalue == "skew":
        table["primary_pvalue"] = table["p_skew"]
        table["primary_q_BH"] = table["q_skew_BH"]
    else:
        table["primary_pvalue"] = table["p_empirical"]
        table["primary_q_BH"] = table["q_empirical_BH"]

    return table


def signed_stouffer_two_mouse(
    table_a: pd.DataFrame,
    table_b: pd.DataFrame,
) -> pd.DataFrame:
    key = ["target_name", "program_name"]

    left = table_a[
        key + ["beta_clr", "primary_pvalue", "primary_q_BH"]
    ].rename(
        columns={
            "beta_clr": "beta_IP_GFP",
            "primary_pvalue": "p_IP_GFP",
            "primary_q_BH": "q_IP_GFP",
        }
    )
    right = table_b[
        key + ["beta_clr", "primary_pvalue", "primary_q_BH"]
    ].rename(
        columns={
            "beta_clr": "beta_GV_GFP",
            "primary_pvalue": "p_GV_GFP",
            "primary_q_BH": "q_GV_GFP",
        }
    )

    merged = left.merge(right, on=key, how="inner")
    merged["same_sign"] = (
        np.sign(merged["beta_IP_GFP"])
        == np.sign(merged["beta_GV_GFP"])
    )

    p_ip = merged["p_IP_GFP"].clip(1e-300, 1.0)
    p_gv = merged["p_GV_GFP"].clip(1e-300, 1.0)

    z_ip = np.sign(merged["beta_IP_GFP"]) * norm.isf(
        p_ip / 2
    )
    z_gv = np.sign(merged["beta_GV_GFP"]) * norm.isf(
        p_gv / 2
    )
    z_combined = (z_ip + z_gv) / math.sqrt(2)

    merged["signed_stouffer_z"] = z_combined
    merged["signed_stouffer_p"] = (
        2 * norm.sf(np.abs(z_combined))
    )
    merged["signed_stouffer_q_BH"] = multipletests(
        merged["signed_stouffer_p"],
        method="fdr_bh",
    )[1]
    return merged


def run_full_crt_for_condition(
    condition: str,
    *,
    variant_name: str,
    B: int,
    primary_pvalue: str,
    overwrite: bool = False,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    specification = COVARIATE_VARIANTS[variant_name]
    output_dir = (
        DISCOVERY_ROOT
        / variant_name
        / safe_name(condition)
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    table_path = output_dir / "target_program_CRT.tsv"
    if table_path.exists() and not overwrite:
        return pd.read_csv(table_path, sep="\t"), {}

    adata = build_crt_adata(
        specification["covariates"],
        specification["log_covariates"],
    )
    adata_condition = condition_subset(adata, condition)
    inputs = prepare_inputs(adata_condition)

    print(
        f"Running full CRT for {condition}: "
        f"{adata_condition.n_obs:,} cells, B={B:,}"
    )

    out = run_all_genes_union_crt(
        inputs=inputs,
        B=B,
        n_jobs=N_JOBS,
        calibrate_skew_normal=True,
        return_raw_pvals=True,
        return_skew_normal=True,
    )

    # Recreate matched NTC pseudo-targets for the final QQ plot.
    (
        ntc_guides,
        guide_frequency,
        guide_to_bin,
        real_gene_bin_signatures,
    ) = build_ntc_group_inputs(
        inputs=inputs,
        ntc_label=CONTROL_TARGET_LABELS,
        group_size=GUIDES_PER_PSEUDO_TARGET,
        n_bins=NTC_FREQUENCY_BINS,
    )
    groups = make_ntc_groups_ensemble(
        ntc_guides=ntc_guides,
        ntc_freq=guide_frequency,
        real_gene_bin_sigs=real_gene_bin_signatures,
        guide_to_bin=guide_to_bin,
        n_ensemble=NTC_ENSEMBLES,
        seed0=NTC_GROUP_SEED,
        group_size=GUIDES_PER_PSEUDO_TARGET,
        max_groups=NTC_MAX_GROUPS,
    )

    ntc_empirical = crt_pvals_for_ntc_groups_ensemble(
        inputs=inputs,
        ntc_groups_ens=groups,
        B=B,
        seed0=NTC_CRT_SEED,
    )
    ntc_skew = crt_pvals_for_ntc_groups_ensemble_skew(
        inputs=inputs,
        ntc_groups_ens=groups,
        B=B,
        seed0=NTC_CRT_SEED,
    )

    guide_to_column = {
        guide: index
        for index, guide in enumerate(inputs.guide_names)
    }
    null_pvalues = np.concatenate(
        [
            compute_guide_set_null_pvals(
                guide_idx=[
                    guide_to_column[guide]
                    for guide in guides
                ],
                inputs=inputs,
                B=B,
            ).ravel()
            for ensemble_groups in groups
            for guides in ensemble_groups.values()
        ]
    )

    skew_frame = out.get("pvals_df")
    if skew_frame is None:
        skew_frame = out.get("pvals_skew_df")

    axis = qq_plot_ntc_pvals(
        pvals_raw_df=out["pvals_raw_df"],
        guide2gene=adata_condition.uns["guide2gene"],
        ntc_genes=CONTROL_TARGET_LABELS,
        pvals_skew_df=skew_frame,
        null_pvals=null_pvalues,
        ntc_group_pvals_ens=ntc_empirical,
        ntc_group_pvals_skew_ens=ntc_skew,
        show_ntc_ensemble_band=True,
        show_all_pvals=True,
        title=f"K30, dt2.0 CRT calibration: {condition}",
    )
    figure = axis.figure
    figure.tight_layout()
    figure.savefig(
        output_dir / "CRT_calibration_QQ.png",
        dpi=180,
        bbox_inches="tight",
    )
    display(figure)
    plt.close(figure)

    table = long_results_table(
        out,
        condition=condition,
        primary_pvalue=primary_pvalue,
    )
    table.to_csv(table_path, sep="\t", index=False)

    metadata = {
        "K": K,
        "density_threshold": DENSITY_THRESHOLD,
        "condition": condition,
        "variant": variant_name,
        "covariates": specification["covariates"],
        "log_covariates": specification["log_covariates"],
        "B": B,
        "n_jobs": N_JOBS,
        "primary_pvalue": primary_pvalue,
        "effect_scale": "CLR log-ratio beta",
    }
    (output_dir / "config.json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )
    return table, out

```

## Markdown cell 23

### Execute the full run

Set `RUN_FULL_DISCOVERY = True` only after reviewing the calibration grid and setting `SELECTED_VARIANT` and `PRIMARY_PVALUE`.

The production default is 5,000 permutations. This can be computationally intensive in an interactive Jupyter session; use a sufficiently large Sherlock Jupyter allocation.


## Code cell 24

```python
RUN_FULL_DISCOVERY = True
OVERWRITE_DISCOVERY = False

discovery_tables: dict[str, pd.DataFrame] = {}

if RUN_FULL_DISCOVERY:
    started = time.time()

    for condition in conditions:
        table, _ = run_full_crt_for_condition(
            condition,
            variant_name=SELECTED_VARIANT,
            B=DISCOVERY_B,
            primary_pvalue=PRIMARY_PVALUE,
            overwrite=OVERWRITE_DISCOVERY,
        )
        discovery_tables[condition] = table
        display(
            table.sort_values("primary_q_BH").head(20)
        )

    print(
        "Elapsed minutes:",
        round((time.time() - started) / 60, 2),
    )
else:
    print(
        "Full discovery is disabled. Set RUN_FULL_DISCOVERY=True "
        "after selecting the calibration specification."
    )

```

### Text output (stream)

```text
Running full CRT for IP_GFP: 2,849 cells, B=5,000

```

### Text output (stream)

```text
/tmp/ipykernel_21869/3815331543.py:242: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  figure.tight_layout()
/tmp/ipykernel_21869/3815331543.py:243: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  figure.savefig(
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/IPython/core/pylabtools.py:170: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  fig.canvas.print_figure(bytes_io, **kw)

```

### Text output (text/plain)

```text
<Figure size 600x500 with 1 Axes>
```

### Text output (text/plain)

```text
    target_name program_name  p_empirical  beta_clr    p_skew  q_empirical_BH  \
28        Smad4            2     0.000200 -1.435775  0.000049        0.033993   
77          Kdr            5     0.000200  1.414975  0.000085        0.033993   
146      Notch1            9     0.000200 -1.518763  0.000349        0.033993   
348       Kcnk3           21     0.001600 -1.187703  0.000594        0.203959   
132       Sox17            8     0.002000 -1.456698  0.002798        0.203959   
120     Atp13a3            8     0.002999 -1.088562  0.002660        0.254949   
459      Acvrl1           28     0.004199  1.215362  0.004129        0.267696   
141     Eif2ak4            9     0.006599  1.235948  0.004149        0.305939   
27       Notch1            2     0.005199 -1.014715  0.004902        0.294608   
503      Notch1           30     0.004199  0.933932  0.005567        0.267696   
142         Eng            9     0.008398 -1.225768  0.007660        0.345162   
209     Eif2ak4           13     0.008798 -1.147002  0.010057        0.345162   
9           Kdr            1     0.005799 -0.856077  0.010600        0.295741   
167        Tbx4           10     0.014197 -1.011529  0.015587        0.441912   
336       Sox17           20     0.015397  1.246134  0.016012        0.441912   
116        Tbx4            7     0.017996  1.085852  0.015145        0.469106   
462        Cav1           28     0.020796 -0.763041  0.013979        0.482085   
261         Eng           16     0.022995  0.920033  0.015149        0.488652   
278         Eng           17     0.014197 -0.986488  0.014520        0.441912   
91          Eng            6     0.012398  1.114963  0.016333        0.441912   

     q_skew_BH  approx_relative_log2fc condition primary_pvalue_type  \
28    0.021787               -2.142812    IP_GFP                skew   
77    0.021787                2.111770    IP_GFP                skew   
146   0.059312               -2.266667    IP_GFP                skew   
348   0.075795               -1.772580    IP_GFP                skew   
132   0.237831               -2.174039    IP_GFP                skew   
120   0.237831               -1.624617    IP_GFP                skew   
459   0.264478                1.813858    IP_GFP                skew   
141   0.264478                1.844583    IP_GFP                skew   
27    0.277790               -1.514405    IP_GFP                skew   
503   0.283906                1.393840    IP_GFP                skew   
142   0.355136               -1.829389    IP_GFP                skew   
209   0.415865               -1.711835    IP_GFP                skew   
9     0.415865               -1.277647    IP_GFP                skew   
167   0.416500               -1.509649    IP_GFP                skew   
336   0.416500                1.859784    IP_GFP                skew   
116   0.416500                1.620572    IP_GFP                skew   
462   0.416500               -1.138795    IP_GFP                skew   
261   0.416500                1.373097    IP_GFP                skew   
278   0.416500               -1.472277    IP_GFP                skew   
91    0.416500                1.664020    IP_GFP                skew   

     primary_pvalue  primary_q_BH  
28         0.000049      0.021787  
77         0.000085      0.021787  
146        0.000349      0.059312  
348        0.000594      0.075795  
132        0.002798      0.237831  
120        0.002660      0.237831  
459        0.004129      0.264478  
141        0.004149      0.264478  
27         0.004902      0.277790  
503        0.005567      0.283906  
142        0.007660      0.355136  
209        0.010057      0.415865  
9          0.010600      0.415865  
167        0.015587      0.416500  
336        0.016012      0.416500  
116        0.015145      0.416500  
462        0.013979      0.416500  
261        0.015149      0.416500  
278        0.014520      0.416500  
91         0.016333      0.416500  
```

### Text output (stream)

```text
Running full CRT for GV_GFP: 1,246 cells, B=5,000

```

### Text output (stream)

```text
/tmp/ipykernel_21869/3815331543.py:242: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  figure.tight_layout()
/tmp/ipykernel_21869/3815331543.py:243: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  figure.savefig(
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/IPython/core/pylabtools.py:170: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  fig.canvas.print_figure(bytes_io, **kw)

```

### Text output (text/plain)

```text
<Figure size 600x500 with 1 Axes>
```

### Text output (text/plain)

```text
       target_name program_name  p_empirical  beta_clr    p_skew  \
348          Kcnk3           21     0.002000 -1.922923  0.001264   
3             Cav1            1     0.014997 -1.206155  0.016902   
48            Tbx4            3     0.013197 -1.412683  0.020890   
7             Gdf2            1     0.013597  1.446777  0.011345   
460        Atp13a3           28     0.029794  1.090428  0.027983   
272         Acvrl1           17     0.033193 -1.768372  0.031964   
268          Sox17           16     0.007998  2.051879  0.010103   
456           Tbx4           27     0.016997  1.247864  0.023373   
237  non-targeting           14     0.018596 -0.970403  0.019378   
85          Acvrl1            6     0.030994  2.165996  0.034446   
409        Atp13a3           25     0.026395  1.093609  0.030906   
338          Trp53           20     0.015397 -1.951276  0.006709   
119         Acvrl1            8     0.079184  1.545133  0.032168   
393          Bmpr2           24     0.032194  1.597395  0.034570   
92            Gdf2            6     0.022196 -2.073819  0.020178   
239        Atp13a3           15     0.014397 -1.360228  0.013093   
274          Bmpr2           17     0.017596  1.745996  0.016675   
12           Smad9            1     0.016597  1.476967  0.009529   
279           Gdf2           17     0.031194  1.722615  0.028197   
19           Bmpr2            2     0.030994 -1.886793  0.029609   

     q_empirical_BH  q_skew_BH  approx_relative_log2fc condition  \
348        0.790342   0.644685               -2.869853    GV_GFP   
3          0.790342   0.801391               -1.800118    GV_GFP   
48         0.790342   0.801391               -2.108349    GV_GFP   
7          0.790342   0.801391                2.159232    GV_GFP   
460        0.803599   0.801391                1.627402    GV_GFP   
272        0.803599   0.801391               -2.639195    GV_GFP   
268        0.790342   0.801391                3.062313    GV_GFP   
456        0.790342   0.801391                1.862367    GV_GFP   
237        0.790342   0.801391               -1.448271    GV_GFP   
85         0.803599   0.801391                3.232627    GV_GFP   
409        0.803599   0.801391                1.632149    GV_GFP   
338        0.790342   0.801391               -2.912169    GV_GFP   
119        0.834379   0.801391                2.306023    GV_GFP   
393        0.803599   0.801391                2.384022    GV_GFP   
92         0.803599   0.801391               -3.095057    GV_GFP   
239        0.790342   0.801391               -2.030063    GV_GFP   
274        0.790342   0.801391                2.605800    GV_GFP   
12         0.790342   0.801391                2.204290    GV_GFP   
279        0.803599   0.801391                2.570904    GV_GFP   
19         0.803599   0.801391               -2.815931    GV_GFP   

    primary_pvalue_type  primary_pvalue  primary_q_BH  
348                skew        0.001264      0.644685  
3                  skew        0.016902      0.801391  
48                 skew        0.020890      0.801391  
7                  skew        0.011345      0.801391  
460                skew        0.027983      0.801391  
272                skew        0.031964      0.801391  
268                skew        0.010103      0.801391  
456                skew        0.023373      0.801391  
237                skew        0.019378      0.801391  
85                 skew        0.034446      0.801391  
409                skew        0.030906      0.801391  
338                skew        0.006709      0.801391  
119                skew        0.032168      0.801391  
393                skew        0.034570      0.801391  
92                 skew        0.020178      0.801391  
239                skew        0.013093      0.801391  
274                skew        0.016675      0.801391  
12                 skew        0.009529      0.801391  
279                skew        0.028197      0.801391  
19                 skew        0.029609      0.801391  
```

### Text output (stream)

```text
Elapsed minutes: 14.04

```

## Markdown cell 25

## Cross-mouse summary

## Code cell 26

```python
# Load completed per-mouse tables from disk.
discovery_tables = {}

for condition in conditions:
    path = (
        DISCOVERY_ROOT
        / SELECTED_VARIANT
        / safe_name(condition)
        / "target_program_CRT.tsv"
    )
    if path.exists():
        discovery_tables[condition] = pd.read_csv(
            path,
            sep="\t",
        )

if {"IP_GFP", "GV_GFP"}.issubset(discovery_tables):
    cross_mouse = signed_stouffer_two_mouse(
        discovery_tables["IP_GFP"],
        discovery_tables["GV_GFP"],
    )
    cross_mouse.to_csv(
        DISCOVERY_ROOT
        / SELECTED_VARIANT
        / "cross_mouse_target_program_summary.tsv",
        sep="\t",
        index=False,
    )

    display(
        cross_mouse.sort_values(
            "signed_stouffer_q_BH"
        ).head(30)
    )

    print(
        "Same-sign fraction:",
        float(cross_mouse["same_sign"].mean()),
    )
else:
    print(
        "Both per-mouse discovery tables are required before "
        "creating the cross-mouse summary."
    )

```

### Text output (text/plain)

```text
    target_name  program_name  beta_IP_GFP  p_IP_GFP  q_IP_GFP  beta_GV_GFP  \
348       Kcnk3            21    -1.187703  0.000594  0.075795    -1.922923   
146      Notch1             9    -1.518763  0.000349  0.059312    -0.944759   
338       Trp53            20    -0.857691  0.055520  0.643523    -1.951276   
423       Trp53            25    -0.859412  0.044762  0.605575    -1.230048   
141     Eif2ak4             9     1.235948  0.004149  0.264478     0.753258   
28        Smad4             2    -1.435775  0.000049  0.021787     0.106286   
209     Eif2ak4            13    -1.147002  0.010057  0.415865    -0.956724   
274       Bmpr2            17     0.755403  0.081705  0.693913     1.745996   
73      Eif2ak4             5     0.866946  0.020545  0.423077     0.572211   
120     Atp13a3             8    -1.088562  0.002660  0.237831    -0.367209   
336       Sox17            20     1.246134  0.016012  0.416500     1.007969   
393       Bmpr2            24     0.595373  0.175750  0.794827     1.597395   
77          Kdr             5     1.414975  0.000085  0.021787    -0.147726   
503      Notch1            30     0.933932  0.005567  0.283906     0.430120   
49        Trp53             3     0.943357  0.020739  0.423077     1.033564   
48         Tbx4             3    -0.481140  0.234967  0.806136    -1.412683   
456        Tbx4            27     0.477920  0.229066  0.806136     1.247864   
23          Eng             2    -0.799938  0.050664  0.630215    -1.243817   
167        Tbx4            10    -1.011529  0.015587  0.416500    -0.740270   
7          Gdf2             1     0.220988  0.370680  0.875344     1.446777   
227         Eng            14     0.864510  0.037745  0.605575     0.833457   
145         Kdr             9     0.607393  0.142955  0.794827     0.974279   
142         Eng             9    -1.225768  0.007660  0.355136    -0.237975   
139        Cav1             9     0.512230  0.219236  0.806136     1.138123   
476      Acvrl1            29     0.369643  0.356497  0.875344     1.560007   
340      Acvrl1            21    -0.539523  0.177295  0.794827    -1.457551   
116        Tbx4             7     1.085852  0.015145  0.416500     0.389694   
114       Smad9             7    -1.101619  0.059904  0.647221    -1.199778   
147       Smad4             9    -0.898665  0.022940  0.449973    -0.482597   
409     Atp13a3            25     0.285732  0.437038  0.902912     1.093609   

     p_GV_GFP  q_GV_GFP  same_sign  signed_stouffer_z  signed_stouffer_p  \
348  0.001264  0.644685       True          -4.708013           0.000003   
146  0.154500  0.878319       True          -3.535402           0.000407   
338  0.006709  0.801391       True          -3.270878           0.001072   
423  0.054258  0.849658       True          -2.780103           0.005434   
141  0.270285  0.883623       True           2.806533           0.005008   
28   0.912639  0.975440      False          -2.795177           0.005187   
209  0.109869  0.878319       True          -2.950506           0.003173   
274  0.016675  0.801391       True           2.923662           0.003459   
73   0.217847  0.878319       True           2.509182           0.012101   
120  0.529408  0.975440       True          -2.569203           0.010193   
336  0.236491  0.878319       True           2.540247           0.011077   
393  0.034570  0.801391       True           2.451765           0.014216   
77   0.735712  0.975440      False           2.539269           0.011108   
503  0.349404  0.895457       True           2.621959           0.008743   
49   0.133191  0.878319       True           2.697144           0.006994   
48   0.020890  0.801391       True          -2.473197           0.013391   
456  0.023373  0.801391       True           2.453695           0.014140   
23   0.078993  0.853343       True          -2.623981           0.008691   
167  0.245077  0.878319       True          -2.532027           0.011341   
7    0.011345  0.801391       True           2.423322           0.015379   
227  0.191034  0.878319       True           2.393653           0.016682   
145  0.147407  0.878319       True           2.060246           0.039375   
142  0.743552  0.975440       True          -2.116980           0.034262   
139  0.085243  0.853343       True           2.085687           0.037007   
476  0.044922  0.849658       True           2.070016           0.038451   
340  0.065526  0.849658       True          -2.256227           0.024056   
116  0.572555  0.975440       True           2.116503           0.034302   
114  0.163560  0.878319       True          -2.315560           0.020582   
147  0.428154  0.918684       True          -2.168548           0.030117   
409  0.030906  0.801391       True           2.075704           0.037921   

     signed_stouffer_q_BH  
348              0.001276  
146              0.103825  
338              0.182264  
423              0.346428  
141              0.346428  
28               0.346428  
209              0.346428  
274              0.346428  
73               0.381581  
120              0.381581  
336              0.381581  
393              0.381581  
77               0.381581  
503              0.381581  
49               0.381581  
48               0.381581  
456              0.381581  
23               0.381581  
167              0.381581  
7                0.392172  
227              0.405123  
145              0.427261  
142              0.427261  
139              0.427261  
476              0.427261  
340              0.427261  
116              0.427261  
114              0.427261  
147              0.427261  
409              0.427261  
```

### Text output (stream)

```text
Same-sign fraction: 0.5235294117647059

```

## Markdown cell 27

## Interpretation notes

- The CRT coefficient is a **CLR log-ratio effect**, not a conventional gene-expression log2 fold change.
- Empirical and skew-calibrated p-values are both retained.
- BH correction is applied within each mouse across all target–program tests.
- The signed-Stouffer table is an exploratory two-mouse summary. With only two biological replicates, it should not be treated as a substitute for a larger replicate-level analysis.
- Because this is a high-MOI experiment, inspect NTC calibration, guide burden, and cross-mouse direction consistency before prioritizing associations.


