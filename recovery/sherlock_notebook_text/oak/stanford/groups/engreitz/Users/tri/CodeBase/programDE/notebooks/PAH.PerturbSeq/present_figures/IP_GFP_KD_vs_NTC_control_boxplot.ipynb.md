# Notebook mirror: `oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/PAH.PerturbSeq/present_figures/IP_GFP_KD_vs_NTC_control_boxplot.ipynb`

- Source archive path: `/oak/stanford/groups/engreitz/Users/tri/CodeBase/cc-perturbseq-figure-reproduction/recovery/sherlock_notebooks/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/PAH.PerturbSeq/present_figures/IP_GFP_KD_vs_NTC_control_boxplot.ipynb`
- SHA-256: `236af21e32180faff24a4bf0569de8c457420c4bd93c75d044b0b1852efe1d47`
- Kernel language: `python`
- Full-fidelity notebook storage: Git LFS

## Markdown cell 0

# IP_GFP target-gene knockdown versus NTC pseudo-target control

This notebook makes the simple two-column knockdown plot requested for the **IP_GFP mouse only**.

For each real target gene:

\[
\mathrm{KD}\% =
100\left(
1 -
\frac{\mathrm{target\ gene\ pseudobulk\ CPM\ in\ target-guide\ cells}}
{\mathrm{target\ gene\ pseudobulk\ CPM\ in\ NTC-positive,\ target-negative\ cells}}
\right).
\]

The high-MOI reference permits other perturbations in both groups.

For the matched NTC control, the notebook repeatedly:

1. selects the same number of NTC guides as the number of guides targeting that gene;
2. treats cells carrying any selected NTC guide as a pseudo-target group;
3. compares them with cells carrying an unselected NTC guide but none of the selected NTC guides;
4. measures the **same target gene**;
5. excludes cells carrying a true guide for that target gene.

The mean apparent KD across the random NTC groups is one control point for that gene. Thus, the final boxplot has one target point and one matched NTC-control point per evaluable target gene.

This is a descriptive screen-QC plot. It is not a biological-replicate significance test.

## Code cell 1

```python
# Run before importing NumPy/SciPy.
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

from pathlib import Path
import hashlib
import json
import re
import sys

from IPython.display import display
import matplotlib.pyplot as plt
import muon as mu
import numpy as np
import pandas as pd
import scipy.sparse as sp

# Keep text editable when the PDF is opened in Adobe Illustrator.
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["font.family"] = "DejaVu Sans"

print("Python:", sys.version.split()[0])

```

### Text output (stream)

```text
Python: 3.11.14

```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/muon/_core/preproc.py:31: FutureWarning: `__version__` is deprecated, use `importlib.metadata.version('scanpy')` instead
  if Version(scanpy.__version__) < Version("1.10"):

```

## Markdown cell 3

## Configuration

## Code cell 4

```python
# ---------------------------------------------------------------------
# Input and output
# ---------------------------------------------------------------------
K = 30
DENSITY_THRESHOLD = 2.0
SAMPLE_NAME = "IP_GFP"

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
OUTPUT_DIR = (
    RUN_ROOT
    / "Evaluation"
    / "30_2_0"
    / "CRT_notebook"
    / "knockdown_estimation"
    / "IP_GFP_target_vs_NTC_control"
)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------
# Analysis settings
# ---------------------------------------------------------------------
CONTROL_TARGET_LABELS = {"non-targeting"}

# Number of random NTC pseudo-target groups averaged for each control point.
# Set to 1 for exactly one random group per target gene.
N_NTC_DRAWS = 100
RANDOM_SEED = 20260621

# A target/control estimate is retained for plotting only if both groups
# have at least this many cells.
MIN_CELLS_PER_GROUP = 10
MIN_VALID_NTC_DRAWS = 10

print("H5MU:", H5MU_PATH)
print("Output:", OUTPUT_DIR)
print("NTC random draws per target:", N_NTC_DRAWS)

```

### Text output (stream)

```text
H5MU: /oak/stanford/groups/engreitz/Users/tri/10xPerturbSeq/analysis/250507_ePAH1_NovaSeq/CRISPR_count/IP_GV_unified_sceptre_cnmf_no_celltype/PerturbNMF/results/IP_GV_skcnmf_K30_dt2_0_n100/Inference/adata/cNMF_30_2_0.h5mu
Output: /oak/stanford/groups/engreitz/Users/tri/10xPerturbSeq/analysis/250507_ePAH1_NovaSeq/CRISPR_count/IP_GV_unified_sceptre_cnmf_no_celltype/PerturbNMF/results/IP_GV_skcnmf_K30_dt2_0_n100/Evaluation/30_2_0/CRT_notebook/knockdown_estimation/IP_GFP_target_vs_NTC_control
NTC random draws per target: 100

```

## Markdown cell 5

## Load RNA counts and guide assignments

## Code cell 6

```python
mdata = mu.read(H5MU_PATH)

if "rna" not in mdata.mod or "cNMF" not in mdata.mod:
    raise KeyError("Expected 'rna' and 'cNMF' modalities.")

rna = mdata["rna"]
programs = mdata["cNMF"]

if rna.n_obs != programs.n_obs:
    raise ValueError("RNA and cNMF modalities have different cell counts.")
if not rna.obs_names.equals(programs.obs_names):
    raise ValueError("RNA and cNMF cell order differs.")


def resolve_obs_series(key: str) -> pd.Series:
    if key in programs.obs:
        values = programs.obs[key]
    elif key in rna.obs:
        values = rna.obs.loc[programs.obs_names, key]
    else:
        raise KeyError(f"{key!r} is absent from RNA.obs and cNMF.obs.")
    return pd.Series(
        values.to_numpy(),
        index=programs.obs_names,
        name=key,
    )


# Guide assignment
if "guide_assignment" in programs.obsm:
    guide_assignment = programs.obsm["guide_assignment"]
elif "guide_assignments" in programs.obsm:
    guide_assignment = programs.obsm["guide_assignments"]
elif "guide_assignment" in rna.obsm:
    guide_assignment = rna.obsm["guide_assignment"]
elif "guide_assignments" in rna.obsm:
    guide_assignment = rna.obsm["guide_assignments"]
else:
    raise KeyError("No guide-assignment matrix was found.")

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
        guide_targets = np.asarray(modality.uns["guide_targets"]).astype(str)

if guide_names is None or guide_targets is None:
    raise KeyError("guide_names or guide_targets are missing.")
if guide_assignment.shape != (rna.n_obs, len(guide_names)):
    raise ValueError(
        f"Guide matrix shape {guide_assignment.shape} does not match "
        f"{rna.n_obs} cells x {len(guide_names)} guides."
    )

sample = resolve_obs_series("sample").astype(str)
ip_mask = (sample == SAMPLE_NAME).to_numpy()

if ip_mask.sum() == 0:
    raise ValueError(f"No cells have sample == {SAMPLE_NAME!r}.")

print(mdata)
print(f"{SAMPLE_NAME} cells:", int(ip_mask.sum()))
print("Guides:", len(guide_names))

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
IP_GFP cells: 2849
Guides: 96

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
def _matrix_values_for_check(matrix, max_values: int = 200_000):
    if sp.issparse(matrix):
        values = np.asarray(matrix.data)
    else:
        values = np.asarray(matrix).ravel()
    values = values[np.isfinite(values)]
    if values.size > max_values:
        positions = np.linspace(
            0,
            values.size - 1,
            max_values,
            dtype=int,
        )
        values = values[positions]
    return values


def _is_nonnegative_integer_like(
    matrix,
    tolerance: float = 1e-3,
) -> bool:
    values = _matrix_values_for_check(matrix)
    if values.size == 0:
        return True
    if values.min() < 0:
        return False
    return bool(
        np.max(np.abs(values - np.rint(values))) <= tolerance
    )


def resolve_rna_count_matrix(adata):
    candidates = []
    for layer_name in ["counts", "raw_counts", "count"]:
        if layer_name in adata.layers:
            candidates.append(
                (f"rna.layers[{layer_name!r}]", adata.layers[layer_name])
            )
    candidates.append(("rna.X", adata.X))

    for source_name, matrix in candidates:
        if _is_nonnegative_integer_like(matrix):
            if sp.issparse(matrix):
                return matrix.tocsr().astype(np.float64), source_name
            return np.asarray(matrix, dtype=np.float64), source_name

    raise ValueError(
        "No nonnegative integer-like RNA count matrix was found."
    )


def row_sum(matrix) -> np.ndarray:
    if sp.issparse(matrix):
        return np.asarray(matrix.sum(axis=1)).ravel()
    return np.asarray(matrix).sum(axis=1)


def gene_vector(matrix, columns: list[int]) -> np.ndarray:
    if not columns:
        return np.zeros(matrix.shape[0], dtype=float)
    if sp.issparse(matrix):
        return np.asarray(matrix[:, columns].sum(axis=1)).ravel()
    return np.asarray(matrix[:, columns]).sum(axis=1)


def pseudobulk_cpm(
    gene_counts: np.ndarray,
    library_counts: np.ndarray,
    mask: np.ndarray,
) -> float:
    indices = np.flatnonzero(mask)
    if indices.size == 0:
        return np.nan
    library_sum = float(library_counts[indices].sum())
    if library_sum <= 0:
        return np.nan
    gene_sum = float(gene_counts[indices].sum())
    return 1e6 * gene_sum / library_sum


def kd_percent(
    gene_counts: np.ndarray,
    library_counts: np.ndarray,
    exposed_mask: np.ndarray,
    reference_mask: np.ndarray,
) -> tuple[float, float, float]:
    exposed_cpm = pseudobulk_cpm(
        gene_counts,
        library_counts,
        exposed_mask,
    )
    reference_cpm = pseudobulk_cpm(
        gene_counts,
        library_counts,
        reference_mask,
    )
    if not np.isfinite(reference_cpm) or reference_cpm <= 0:
        return np.nan, exposed_cpm, reference_cpm
    value = 100.0 * (1.0 - exposed_cpm / reference_cpm)
    return float(value), exposed_cpm, reference_cpm


def stable_target_seed(target: str) -> int:
    digest = hashlib.sha256(target.encode("utf-8")).digest()
    offset = int.from_bytes(digest[:4], "little")
    return (RANDOM_SEED + offset) % (2**32)

```

## Code cell 8

```python
rna_counts, rna_count_source = resolve_rna_count_matrix(rna)
library_counts = row_sum(rna_counts)

if np.any(library_counts <= 0):
    raise ValueError("At least one RNA count row has zero library size.")

if "symbol" in rna.var:
    rna_symbols = rna.var["symbol"].astype(str)
else:
    rna_symbols = pd.Series(
        rna.var_names.astype(str),
        index=rna.var_names,
        name="symbol",
    )

symbol_to_columns: dict[str, list[int]] = {}
for column_index, symbol in enumerate(rna_symbols):
    key = str(symbol).strip().upper()
    if key and key not in {"NAN", "NONE"}:
        symbol_to_columns.setdefault(key, []).append(column_index)

control_keys = {
    str(label).strip().lower()
    for label in CONTROL_TARGET_LABELS
}
normalized_targets = np.asarray(
    [str(target).strip().upper() for target in guide_targets]
)

ntc_columns = np.flatnonzero(
    np.asarray(
        [
            str(target).strip().lower() in control_keys
            for target in guide_targets
        ]
    )
)
if len(ntc_columns) == 0:
    raise ValueError("No NTC guides were found.")

ntc_positive = np.asarray(
    guide_assignment[:, ntc_columns].sum(axis=1)
).ravel() > 0

real_targets = sorted(
    {
        str(target)
        for target in guide_targets
        if str(target).strip().lower() not in control_keys
    },
    key=str.upper,
)

print("RNA count source:", rna_count_source)
print("Real targets:", len(real_targets))
print("NTC guides:", len(ntc_columns))

```

### Text output (stream)

```text
RNA count source: rna.layers['counts']
Real targets: 16
NTC guides: 30

```

## Markdown cell 9

## Calculate real KD and matched NTC control

## Code cell 10

```python
per_gene_rows = []
null_rows = []

all_ntc_columns = np.asarray(ntc_columns, dtype=int)

for target in real_targets:
    target_key = target.strip().upper()
    target_guide_columns = np.flatnonzero(
        normalized_targets == target_key
    )
    n_target_guides = int(len(target_guide_columns))

    gene_columns = symbol_to_columns.get(target_key, [])
    if gene_columns:
        gene_counts = gene_vector(rna_counts, gene_columns)
        matched_rna_features = ";".join(
            rna.var_names[np.asarray(gene_columns)].astype(str)
        )
    else:
        gene_counts = np.zeros(rna.n_obs, dtype=float)
        matched_rna_features = ""

    target_positive = np.asarray(
        guide_assignment[:, target_guide_columns].sum(axis=1)
    ).ravel() > 0

    # Real target comparison:
    # target-guide cells versus NTC-positive cells lacking this target.
    real_target_mask = ip_mask & target_positive
    real_reference_mask = (
        ip_mask
        & ntc_positive
        & ~target_positive
    )

    real_kd, real_cpm, real_reference_cpm = kd_percent(
        gene_counts,
        library_counts,
        real_target_mask,
        real_reference_mask,
    )

    # Matched NTC pseudo-target controls.
    # One control point per target gene is the mean across valid random draws.
    rng = np.random.default_rng(stable_target_seed(target))
    current_null_values = []

    if 0 < n_target_guides < len(all_ntc_columns):
        for draw in range(N_NTC_DRAWS):
            selected_ntc = np.sort(
                rng.choice(
                    all_ntc_columns,
                    size=n_target_guides,
                    replace=False,
                )
            )
            unselected_ntc = np.setdiff1d(
                all_ntc_columns,
                selected_ntc,
                assume_unique=True,
            )

            selected_positive = np.asarray(
                guide_assignment[:, selected_ntc].sum(axis=1)
            ).ravel() > 0
            unselected_positive = np.asarray(
                guide_assignment[:, unselected_ntc].sum(axis=1)
            ).ravel() > 0

            # Exclude true target-guide cells from both pseudo groups so the
            # target gene readout is not affected by its real perturbation.
            pseudo_target_mask = (
                ip_mask
                & selected_positive
                & ~target_positive
            )
            pseudo_reference_mask = (
                ip_mask
                & unselected_positive
                & ~selected_positive
                & ~target_positive
            )

            control_kd, control_cpm, control_reference_cpm = kd_percent(
                gene_counts,
                library_counts,
                pseudo_target_mask,
                pseudo_reference_mask,
            )

            valid_draw = bool(
                gene_columns
                and np.isfinite(control_kd)
                and pseudo_target_mask.sum() >= MIN_CELLS_PER_GROUP
                and pseudo_reference_mask.sum() >= MIN_CELLS_PER_GROUP
            )

            null_rows.append(
                {
                    "target_name": target,
                    "draw": draw,
                    "n_guides_in_pseudo_target": n_target_guides,
                    "selected_NTC_guides": ";".join(
                        guide_names[selected_ntc].astype(str)
                    ),
                    "n_pseudo_target_cells": int(
                        pseudo_target_mask.sum()
                    ),
                    "n_pseudo_reference_cells": int(
                        pseudo_reference_mask.sum()
                    ),
                    "pseudo_target_CPM": control_cpm,
                    "pseudo_reference_CPM": control_reference_cpm,
                    "apparent_KD_percent": control_kd,
                    "valid_draw": valid_draw,
                }
            )

            if valid_draw:
                current_null_values.append(control_kd)

    current_null = np.asarray(current_null_values, dtype=float)
    n_valid = int(len(current_null))

    if n_valid:
        control_mean = float(np.mean(current_null))
        control_median = float(np.median(current_null))
        control_sd = float(np.std(current_null, ddof=1)) if n_valid > 1 else 0.0
        control_low, control_high = np.percentile(
            current_null,
            [2.5, 97.5],
        )
    else:
        control_mean = np.nan
        control_median = np.nan
        control_sd = np.nan
        control_low = np.nan
        control_high = np.nan

    plot_eligible = bool(
        gene_columns
        and np.isfinite(real_kd)
        and real_target_mask.sum() >= MIN_CELLS_PER_GROUP
        and real_reference_mask.sum() >= MIN_CELLS_PER_GROUP
        and n_valid >= MIN_VALID_NTC_DRAWS
    )

    status = []
    if not gene_columns:
        status.append("target_gene_not_found")
    if real_target_mask.sum() < MIN_CELLS_PER_GROUP:
        status.append("few_real_target_cells")
    if real_reference_mask.sum() < MIN_CELLS_PER_GROUP:
        status.append("few_real_reference_cells")
    if not np.isfinite(real_kd):
        status.append("real_KD_not_estimable")
    if n_valid < MIN_VALID_NTC_DRAWS:
        status.append("few_valid_NTC_draws")
    if not status:
        status.append("ok")

    per_gene_rows.append(
        {
            "target_name": target,
            "target_gene_found": bool(gene_columns),
            "matched_rna_features": matched_rna_features,
            "n_target_guides": n_target_guides,
            "n_real_target_cells": int(real_target_mask.sum()),
            "n_real_reference_cells": int(
                real_reference_mask.sum()
            ),
            "real_target_CPM": real_cpm,
            "real_reference_CPM": real_reference_cpm,
            "real_target_KD_percent": real_kd,
            "n_valid_NTC_draws": n_valid,
            "NTC_control_mean_apparent_KD_percent": control_mean,
            "NTC_control_median_apparent_KD_percent": control_median,
            "NTC_control_SD": control_sd,
            "NTC_control_2.5_percentile": float(control_low),
            "NTC_control_97.5_percentile": float(control_high),
            "plot_eligible": plot_eligible,
            "status": ";".join(status),
        }
    )

per_gene = pd.DataFrame(per_gene_rows).sort_values("target_name")
ntc_draws = pd.DataFrame(null_rows).sort_values(
    ["target_name", "draw"]
)

display(per_gene)

```

### Text output (text/plain)

```text
   target_name  target_gene_found matched_rna_features  n_target_guides  \
0       Acvrl1               True               Acvrl1                4   
1      Atp13a3               True              Atp13a3                4   
2        Bmpr2               True                Bmpr2                4   
3         Cav1               True                 Cav1                4   
4       Cdkn1a               True               Cdkn1a                5   
5      Eif2ak4               True              Eif2ak4                4   
6          Eng               True                  Eng                4   
7         Gdf2              False                                     4   
8        Kcnk3               True                Kcnk3                4   
9          Kdr               True                  Kdr                4   
10      Notch1               True               Notch1                5   
11       Smad4               True                Smad4                4   
12       Smad9               True                Smad9                4   
13       Sox17               True                Sox17                4   
14        Tbx4               True                 Tbx4                4   
15       Trp53               True                Trp53                4   

    n_real_target_cells  n_real_reference_cells  real_target_CPM  \
0                   193                    1309       922.991382   
1                   365                    1231        28.028365   
2                   215                    1292      1144.439676   
3                   316                    1264      1098.659351   
4                   365                    1206       506.912687   
5                   306                    1262        15.193297   
6                   287                    1268       510.379928   
7                   221                    1297         0.000000   
8                   291                    1255        23.862127   
9                   326                    1253       344.589938   
10                  415                    1217       118.084989   
11                  430                    1238        80.928185   
12                  198                    1307         2.437896   
13                  178                    1317       110.752212   
14                  276                    1278        12.733153   
15                  207                    1313        42.021838   

    real_reference_CPM  real_target_KD_percent  n_valid_NTC_draws  \
0           878.399751               -5.076462                100   
1            39.599989               29.221282                100   
2          1149.725285                0.459728                100   
3          1339.920660               18.005641                100   
4           796.643406               36.368935                100   
5            22.553061               32.633106                100   
6           719.615598               29.076033                100   
7             0.000000                     NaN                  0   
8            22.960140               -3.928487                100   
9           690.447396               50.091790                100   
10          175.971037               32.895213                100   
11           99.908269               18.997511                100   
12            3.620489               32.663905                100   
13          148.968577               25.653977                100   
14           19.073602               33.242011                100   
15           65.067522               35.418107                100   

    NTC_control_mean_apparent_KD_percent  \
0                              -2.884108   
1                               1.834374   
2                              -0.899253   
3                               1.822939   
4                              -5.694640   
5                              -3.265609   
6                               0.289653   
7                                    NaN   
8                              20.979788   
9                              -7.716546   
10                             -0.307983   
11                              2.906180   
12                            -12.870944   
13                              6.778285   
14                              6.754431   
15                              1.126001   

    NTC_control_median_apparent_KD_percent  NTC_control_SD  \
0                                -2.711233        3.485485   
1                                 3.157555       16.406973   
2                                -1.407852        4.525554   
3                                 2.090488        3.150057   
4                                -5.247780        7.760845   
5                                -2.355002       16.457441   
6                                 0.727622        4.855451   
7                                      NaN             NaN   
8                                19.459971       27.235677   
9                                -8.821202        8.051350   
10                                0.114329        5.928657   
11                                2.599318        7.357682   
12                              -18.178325       46.151594   
13                                7.560684       10.022211   
14                               11.573501       23.517075   
15                                0.513966        7.407364   

    NTC_control_2.5_percentile  NTC_control_97.5_percentile  plot_eligible  \
0                    -9.557554                     3.638213           True   
1                   -33.361376                    29.426252           True   
2                    -9.904905                     7.550394           True   
3                    -5.068955                     6.644198           True   
4                   -18.641629                     8.043998           True   
5                   -35.914366                    24.536582           True   
6                    -8.800231                     8.753027           True   
7                          NaN                          NaN          False   
8                   -40.356619                    72.897391           True   
9                   -23.505722                     7.932990           True   
10                  -11.598187                    10.204023           True   
11                  -11.024560                    19.738070           True   
12                  -91.397520                    73.072657           True   
13                  -13.334559                    24.784836           True   
14                  -41.665926                    38.318514           True   
15                  -12.102575                    16.738140           True   

                                               status  
0                                                  ok  
1                                                  ok  
2                                                  ok  
3                                                  ok  
4                                                  ok  
5                                                  ok  
6                                                  ok  
7   target_gene_not_found;real_KD_not_estimable;fe...  
8                                                  ok  
9                                                  ok  
10                                                 ok  
11                                                 ok  
12                                                 ok  
13                                                 ok  
14                                                 ok  
15                                                 ok  
```

## Markdown cell 11

## Export CSV files

## Code cell 12

```python
per_gene_path = (
    OUTPUT_DIR / "IP_GFP_target_vs_NTC_control_per_gene.csv"
)
null_draws_path = (
    OUTPUT_DIR / "IP_GFP_NTC_control_random_draws.csv"
)
plot_data_path = (
    OUTPUT_DIR / "IP_GFP_target_vs_NTC_control_plot_data.csv"
)
group_summary_path = (
    OUTPUT_DIR / "IP_GFP_target_vs_NTC_control_group_summary.csv"
)

per_gene.to_csv(per_gene_path, index=False)
ntc_draws.to_csv(null_draws_path, index=False)

eligible = per_gene.loc[per_gene["plot_eligible"]].copy()

target_plot = eligible[
    ["target_name", "real_target_KD_percent"]
].rename(columns={"real_target_KD_percent": "KD_percent"})
target_plot["group"] = "Targeting guides"

control_plot = eligible[
    ["target_name", "NTC_control_mean_apparent_KD_percent"]
].rename(
    columns={
        "NTC_control_mean_apparent_KD_percent": "KD_percent"
    }
)
control_plot["group"] = "NTC pseudo-target"

plot_data = pd.concat(
    [target_plot, control_plot],
    ignore_index=True,
)
plot_data["group"] = pd.Categorical(
    plot_data["group"],
    categories=["Targeting guides", "NTC pseudo-target"],
    ordered=True,
)
plot_data.to_csv(plot_data_path, index=False)

group_summary = (
    plot_data.groupby("group", observed=True)["KD_percent"]
    .agg(
        n="count",
        mean="mean",
        median="median",
        standard_deviation="std",
        minimum="min",
        maximum="max",
    )
    .reset_index()
)
group_summary.to_csv(group_summary_path, index=False)

config = {
    "sample": SAMPLE_NAME,
    "formula": (
        "100 * (1 - exposed pseudobulk CPM / reference pseudobulk CPM)"
    ),
    "real_reference": (
        "IP_GFP cells with at least one NTC guide and no guide "
        "targeting the tested gene; other perturbations allowed"
    ),
    "control_definition": (
        "same gene measured in cells with a random NTC guide group "
        "versus cells with remaining NTC guides; true target-guide "
        "cells excluded; other perturbations allowed"
    ),
    "NTC_draws_per_target": N_NTC_DRAWS,
    "control_point_summary": "mean apparent KD across valid NTC draws",
    "random_seed": RANDOM_SEED,
    "minimum_cells_per_group": MIN_CELLS_PER_GROUP,
    "minimum_valid_NTC_draws": MIN_VALID_NTC_DRAWS,
    "RNA_count_source": rna_count_source,
}
(OUTPUT_DIR / "IP_GFP_target_vs_NTC_control_config.json").write_text(
    json.dumps(config, indent=2) + "\n"
)

print("Per-gene CSV:", per_gene_path)
print("Full NTC draws:", null_draws_path)
print("Plot-data CSV:", plot_data_path)
print("Group summary:", group_summary_path)
display(group_summary)

```

### Text output (stream)

```text
Per-gene CSV: /oak/stanford/groups/engreitz/Users/tri/10xPerturbSeq/analysis/250507_ePAH1_NovaSeq/CRISPR_count/IP_GV_unified_sceptre_cnmf_no_celltype/PerturbNMF/results/IP_GV_skcnmf_K30_dt2_0_n100/Evaluation/30_2_0/CRT_notebook/knockdown_estimation/IP_GFP_target_vs_NTC_control/IP_GFP_target_vs_NTC_control_per_gene.csv
Full NTC draws: /oak/stanford/groups/engreitz/Users/tri/10xPerturbSeq/analysis/250507_ePAH1_NovaSeq/CRISPR_count/IP_GV_unified_sceptre_cnmf_no_celltype/PerturbNMF/results/IP_GV_skcnmf_K30_dt2_0_n100/Evaluation/30_2_0/CRT_notebook/knockdown_estimation/IP_GFP_target_vs_NTC_control/IP_GFP_NTC_control_random_draws.csv
Plot-data CSV: /oak/stanford/groups/engreitz/Users/tri/10xPerturbSeq/analysis/250507_ePAH1_NovaSeq/CRISPR_count/IP_GV_unified_sceptre_cnmf_no_celltype/PerturbNMF/results/IP_GV_skcnmf_K30_dt2_0_n100/Evaluation/30_2_0/CRT_notebook/knockdown_estimation/IP_GFP_target_vs_NTC_control/IP_GFP_target_vs_NTC_control_plot_data.csv
Group summary: /oak/stanford/groups/engreitz/Users/tri/10xPerturbSeq/analysis/250507_ePAH1_NovaSeq/CRISPR_count/IP_GV_unified_sceptre_cnmf_no_celltype/PerturbNMF/results/IP_GV_skcnmf_K30_dt2_0_n100/Evaluation/30_2_0/CRT_notebook/knockdown_estimation/IP_GFP_target_vs_NTC_control/IP_GFP_target_vs_NTC_control_group_summary.csv

```

### Text output (text/plain)

```text
               group   n       mean     median  standard_deviation    minimum  \
0   Targeting guides  15  24.381486  29.221282           15.974090  -5.076462   
1  NTC pseudo-target  15   0.590171   0.289653            7.655536 -12.870944   

     maximum  
0  50.091790  
1  20.979788  
```

## Markdown cell 13

## Two-column boxplot

## Code cell 14

```python
if plot_data.empty:
    raise ValueError(
        "No target genes passed the plotting filters. "
        "Inspect the per-gene CSV and lower MIN_CELLS_PER_GROUP only "
        "if scientifically justified."
    )

order = ["Targeting guides", "NTC pseudo-target"]
values = [
    plot_data.loc[
        plot_data["group"].astype(str) == group,
        "KD_percent",
    ].to_numpy(dtype=float)
    for group in order
]

target_color = "#17365D"   # dark blue
control_color = "#4D4D4D"  # dark gray
box_colors = ["#DCE6F1", "#E3E3E3"]

fig, ax = plt.subplots(figsize=(4.6, 5.2))

box = ax.boxplot(
    values,
    positions=[0, 1],
    widths=0.55,
    patch_artist=True,
    showfliers=False,
    medianprops={
        "color": "black",
        "linewidth": 1.5,
    },
    whiskerprops={
        "color": "black",
        "linewidth": 1.0,
    },
    capprops={
        "color": "black",
        "linewidth": 1.0,
    },
)

for patch, fill in zip(box["boxes"], box_colors):
    patch.set_facecolor(fill)
    patch.set_edgecolor("black")
    patch.set_linewidth(1.0)

plot_rng = np.random.default_rng(11)
for x, group, color in zip(
    [0, 1],
    order,
    [target_color, control_color],
):
    current = plot_data.loc[
        plot_data["group"].astype(str) == group,
        "KD_percent",
    ].to_numpy(dtype=float)
    jitter = plot_rng.uniform(-0.09, 0.09, size=len(current))
    ax.scatter(
        np.full(len(current), x) + jitter,
        current,
        s=34,
        color=color,
        edgecolor="white",
        linewidth=0.5,
        alpha=0.95,
        zorder=3,
    )

# White diamonds mark the group means.
means = [float(np.mean(current)) for current in values]
ax.scatter(
    [0, 1],
    means,
    marker="D",
    s=50,
    facecolor="white",
    edgecolor="black",
    linewidth=1.0,
    zorder=4,
)

ax.axhline(
    0,
    color="#7F7F7F",
    linestyle="--",
    linewidth=1.0,
)
ax.set_xticks([0, 1])
ax.set_xticklabels(order)
ax.set_ylabel("Estimated target-gene knockdown (%)")
ax.set_title("IP_GFP knockdown efficiency")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.tick_params(axis="both", labelsize=10)

summary_text = (
    f"Mean target KD: {means[0]:.1f}%\n"
    f"Mean NTC apparent KD: {means[1]:.1f}%\n"
    f"n = {len(values[0])} matched genes"
)
ax.text(
    0.03,
    0.97,
    summary_text,
    transform=ax.transAxes,
    ha="left",
    va="top",
    fontsize=9,
)

fig.tight_layout()

pdf_path = OUTPUT_DIR / "IP_GFP_KD_vs_NTC_control_boxplot.pdf"
png_path = OUTPUT_DIR / "IP_GFP_KD_vs_NTC_control_boxplot.png"

fig.savefig(pdf_path, bbox_inches="tight")
fig.savefig(png_path, dpi=300, bbox_inches="tight")
plt.show()

print("Editable PDF:", pdf_path)
print("PNG:", png_path)

```

### Text output (text/plain)

```text
<Figure size 460x520 with 1 Axes>
```

### Text output (stream)

```text
Editable PDF: /oak/stanford/groups/engreitz/Users/tri/10xPerturbSeq/analysis/250507_ePAH1_NovaSeq/CRISPR_count/IP_GV_unified_sceptre_cnmf_no_celltype/PerturbNMF/results/IP_GV_skcnmf_K30_dt2_0_n100/Evaluation/30_2_0/CRT_notebook/knockdown_estimation/IP_GFP_target_vs_NTC_control/IP_GFP_KD_vs_NTC_control_boxplot.pdf
PNG: /oak/stanford/groups/engreitz/Users/tri/10xPerturbSeq/analysis/250507_ePAH1_NovaSeq/CRISPR_count/IP_GV_unified_sceptre_cnmf_no_celltype/PerturbNMF/results/IP_GV_skcnmf_K30_dt2_0_n100/Evaluation/30_2_0/CRT_notebook/knockdown_estimation/IP_GFP_target_vs_NTC_control/IP_GFP_KD_vs_NTC_control_boxplot.png

```

## Markdown cell 15

## Interpretation

- Each **targeting-guides** point is one target gene.
- Each **NTC pseudo-target** point is the mean apparent KD for the same gene across random, guide-count-matched NTC groups.
- Positive values indicate lower expression in the exposed group.
- Values near zero in the NTC column indicate that random NTC grouping does not systematically reduce the matched gene.
- Negative values indicate apparent higher expression rather than knockdown.
- The white diamond is the arithmetic mean; the line inside each box is the median.

