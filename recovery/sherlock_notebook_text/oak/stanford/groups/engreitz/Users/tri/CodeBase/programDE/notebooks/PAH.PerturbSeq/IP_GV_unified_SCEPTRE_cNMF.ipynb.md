# Notebook mirror: `oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/PAH.PerturbSeq/IP_GV_unified_SCEPTRE_cNMF.ipynb`

- Source archive path: `/oak/stanford/groups/engreitz/Users/tri/CodeBase/cc-perturbseq-figure-reproduction/recovery/sherlock_notebooks/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/PAH.PerturbSeq/IP_GV_unified_SCEPTRE_cNMF.ipynb`
- SHA-256: `cb4ad4093bce8935fc3c9df54421932cae9f83b620fc0bf37fd5049c59e06ef5`
- Kernel language: `python`
- Full-fidelity notebook storage: Git LFS

## Markdown cell 0

# Unified IP_GFP + GV_GFP processing, SCEPTRE assignment, and cNMF input — no cell-type filter

This notebook:

1. reads the 10x Gene Expression + CRISPR Guide Capture matrices for `IP_GFP` and `GV_GFP`;
2. applies RNA QC independently to each mouse, without clustering or cell-type filtering;
3. runs SCEPTRE guide assignment separately for each mouse using every RNA-QC-passed cell;
4. preserves a full QC-passed object containing assigned and unassigned cells;
5. keeps cells with at least one SCEPTRE-assigned guide for the cNMF input;
6. concatenates IP_GFP and GV_GFP while retaining raw guide UMI counts, binary guide assignments,
   assigned-guide UMI counts, guide targets, and sample/barcode metadata;
7. writes a canonical `.h5mu`, a cNMF-ready RNA `.h5ad`, and a standalone guide `.h5ad`.

`rna.X` and `rna.layers["counts"]` remain raw RNA counts. In the guide modality,
`guide.X` and `guide.layers["counts"]` are raw guide UMI counts,
`guide.layers["sceptre_assignment"]` is the binary assignment matrix, and
`guide.layers["assigned_counts"]` contains raw guide counts masked by the assignments.
Original 10x barcodes are retained in `obs["barcode"]`; sample-prefixed cell IDs are introduced
only after SCEPTRE assignment.

## Code cell 1

```python
# Run this before importing numpy/scipy/R to avoid nested thread oversubscription.
import os

for var in [
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "RCPP_PARALLEL_NUM_THREADS",
]:
    os.environ[var] = "1"
os.environ["OMP_THREAD_LIMIT"] = "1"

print("Thread limits set.")
```

### Text output (stream)

```text
Thread limits set.

```

## Code cell 2

```python
from pathlib import Path
import json
import re
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import mudata
from mudata import MuData
from scipy import sparse
import scipy.io as sio
from IPython.display import display

warnings.filterwarnings("default")
sc.settings.verbosity = 2
sc.settings.set_figure_params(dpi=100, fontsize=10, dpi_save=300, figsize=(5, 4))

# ---------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------
BASE = Path("/oak/stanford/groups/engreitz/Users/tri/10xPerturbSeq/analysis/250507_ePAH1_NovaSeq/CRISPR_count")

SAMPLES = {
    "IP_GFP": {"data_dir": BASE / "IP_GFP"},
    "GV_GFP": {"data_dir": BASE / "GV_GFP"},
}

# The 10x H5 appears to contain target_gene_name metadata. This CSV is only a fallback.
# It is assumed to describe the same guide library used for both mice.
SHARED_GUIDE_METADATA = BASE / "IP_GFP" / "crispr" / "ePAH1_gRNA.csv"

OUTPUT_ROOT = BASE / "IP_GV_unified_sceptre_cnmf_no_celltype"
OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------
# Same expression-QC thresholds as the earlier IP notebooks
# ---------------------------------------------------------------------
MIN_GENES = 250
MAX_GENES = 8000
MAX_TOTAL_COUNTS = 50000
MAX_PCT_MT = 15


# Apply this once after combining the assigned-guide cells from both mice.
MIN_CELLS_PER_GENE_COMBINED = 10

# SCEPTRE settings. "high" + "mixture" matches the attached rpy2 example.
SCEPTRE_MOI = "high"
ASSIGN_GRNA_METHOD = "mixture"
ASSIGN_PARALLEL = False   # Set True if desired; logs are written per sample.

# Re-write Python-to-R transfer matrices when rerunning preprocessing.
OVERWRITE_TRANSFER_FILES = True

# Saves all RNA-QC-passed cells with SCEPTRE assignments before guide-positive filtering.
# Set False to reduce disk use after the workflow is working as expected.
SAVE_ALL_QC_WITH_ASSIGNMENTS = True

print("OUTPUT_ROOT:", OUTPUT_ROOT)
for sample_id, cfg in SAMPLES.items():
    print(sample_id, "->", cfg["data_dir"])
```

### Text output (stream)

```text
OUTPUT_ROOT: /oak/stanford/groups/engreitz/Users/tri/10xPerturbSeq/analysis/250507_ePAH1_NovaSeq/CRISPR_count/IP_GV_unified_sceptre_cnmf_no_celltype
IP_GFP -> /oak/stanford/groups/engreitz/Users/tri/10xPerturbSeq/analysis/250507_ePAH1_NovaSeq/CRISPR_count/IP_GFP
GV_GFP -> /oak/stanford/groups/engreitz/Users/tri/10xPerturbSeq/analysis/250507_ePAH1_NovaSeq/CRISPR_count/GV_GFP

```

## Markdown cell 3

## RNA QC and SCEPTRE transfer-file helpers

## Code cell 4

```python
def as_csr(x):
    """Return a CSR sparse matrix without densifying."""
    if sparse.issparse(x):
        return x.tocsr()
    return sparse.csr_matrix(x)




def make_anndata_index_names_write_safe(adata):
    """Avoid H5AD index-name/column-name collisions without changing IDs."""
    adata.obs.index.name = None
    adata.var.index.name = None
    return adata


def write_lines(values, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for value in values:
            handle.write(f"{value}\n")


def first_existing_column(df, candidates):
    for col in candidates:
        if col in df.columns:
            return col
    return None


def normalize_target_name(value):
    value = str(value).strip()
    if value.lower() in {"", "nan", "none", "na"}:
        return np.nan
    if re.search(r"(^|[-_ ])non[-_ ]?target(ing)?($|[-_ ])", value, flags=re.I):
        return "non-targeting"
    return value


def read_guide_metadata(path):
    path = Path(path)
    if not path.exists():
        return None
    return pd.read_csv(path, sep=None, engine="python")


def resolve_guide_targets(guide, fallback_metadata=None):
    """Resolve guide ID -> target using embedded 10x metadata first, then a CSV/TSV fallback."""
    targets = pd.Series(np.nan, index=guide.var_names, dtype="object")

    embedded_target_col = first_existing_column(
        guide.var,
        ["grna_target", "target_gene_name", "intended_target_name", "target"],
    )
    if embedded_target_col is not None:
        embedded = guide.var[embedded_target_col].map(normalize_target_name)
        embedded.index = guide.var_names
        targets.loc[embedded.notna()] = embedded.loc[embedded.notna()]

    meta = read_guide_metadata(fallback_metadata) if fallback_metadata is not None else None
    if meta is not None:
        id_col = first_existing_column(meta, ["id", "guide_id", "grna_id", "feature_id"])
        target_col = first_existing_column(
            meta,
            ["target_gene_name", "intended_target_name", "grna_target", "target"],
        )
        if id_col is None or target_col is None:
            raise KeyError(
                f"Guide metadata {fallback_metadata} must contain a guide-ID column and a target column. "
                f"Observed columns: {list(meta.columns)}"
            )
        fallback_map = (
            meta[[id_col, target_col]]
            .dropna(subset=[id_col])
            .assign(**{id_col: lambda x: x[id_col].astype(str)})
            .drop_duplicates(id_col)
            .set_index(id_col)[target_col]
            .map(normalize_target_name)
        )

        identifiers_to_try = {
            "var_names": pd.Series(guide.var_names, index=guide.var_names, dtype="object")
        }
        for col in ["grna_id", "gene_ids", "feature_name"]:
            if col in guide.var.columns:
                identifiers_to_try[col] = guide.var[col].astype(str)
                identifiers_to_try[col].index = guide.var_names

        for _, identifiers in identifiers_to_try.items():
            unresolved = targets.isna()
            if not unresolved.any():
                break
            mapped = identifiers.map(fallback_map)
            targets.loc[unresolved & mapped.notna()] = mapped.loc[unresolved & mapped.notna()]

    if targets.isna().any():
        missing = targets.index[targets.isna()].tolist()[:20]
        raise ValueError(
            "Could not resolve grna_target for every guide. "
            f"First unresolved guide IDs: {missing}. "
            "Set SHARED_GUIDE_METADATA to the correct guide metadata file or inspect guide.var."
        )

    return targets.astype(str)


def load_and_qc_sample(sample_id, cfg):
    """
    Read one 10x feature-barcode H5, split RNA and guide modalities, apply RNA QC,
    and retain all RNA-QC-passed cells for SCEPTRE assignment.
    """
    data_dir = Path(cfg["data_dir"])
    h5_path = data_dir / "filtered_feature_bc_matrix.h5"

    if not h5_path.exists():
        raise FileNotFoundError(h5_path)

    print(f"\n=== {sample_id}: reading {h5_path} ===")
    all_features = sc.read_10x_h5(h5_path, gex_only=False)
    all_features.obs_names = all_features.obs_names.astype(str)
    all_features.obs.index.name = None
    all_features.var["feature_name"] = all_features.var_names.astype(str)

    feature_type_col = first_existing_column(
        all_features.var,
        ["feature_types", "feature_type"],
    )
    if feature_type_col is None:
        raise KeyError(
            f"No feature-type column found in {h5_path}. "
            f"var columns: {list(all_features.var.columns)}"
        )

    feature_types = all_features.var[feature_type_col].astype(str)
    gex_mask = feature_types.eq("Gene Expression")
    guide_mask = feature_types.str.contains("Guide Capture", case=False, regex=False)

    if int(gex_mask.sum()) == 0 or int(guide_mask.sum()) == 0:
        raise ValueError(
            f"Expected Gene Expression and Guide Capture features in {h5_path}. "
            f"Observed feature types: {feature_types.value_counts().to_dict()}"
        )

    gex = all_features[:, gex_mask].copy()
    guide = all_features[:, guide_mask].copy()
    del all_features

    # Keep gene symbols as RNA var_names for cNMF interpretation.
    # to_numpy() avoids inheriting the Series name as the index name.
    gex.var["symbol"] = gex.var["feature_name"].astype(str)
    gex.var_names = pd.Index(gex.var["symbol"].astype(str).to_numpy(), name=None)
    gex.var_names_make_unique()
    gex.var.index.name = None

    # Prefer stable 10x feature IDs for guide var_names.
    guide.var["feature_name"] = guide.var["feature_name"].astype(str)
    guide_id_col = first_existing_column(
        guide.var,
        ["gene_ids", "id", "feature_id"],
    )
    if guide_id_col is not None:
        candidate_ids = guide.var[guide_id_col].astype(str)
    else:
        candidate_ids = guide.var["feature_name"].astype(str)

    guide.var_names = pd.Index(candidate_ids.astype(str).to_numpy(), name=None)
    guide.var_names_make_unique()
    guide.var.index.name = None
    guide.var["grna_id"] = guide.var_names.astype(str)

    gex.X = as_csr(gex.X)
    guide.X = as_csr(guide.X)

    n_raw_cells = int(gex.n_obs)

    # Compute RNA QC only from Gene Expression counts.
    symbols = gex.var["symbol"].astype(str)
    gex.var["mt"] = symbols.str.match(r"(?i)^mt-")
    gex.var["ribo"] = symbols.str.match(r"(?i)^(rps|rpl)")
    gex.var["hb"] = symbols.str.match(r"(?i)^hb(?!p)")

    sc.pp.calculate_qc_metrics(
        gex,
        qc_vars=["mt", "ribo", "hb"],
        percent_top=None,
        log1p=True,
        inplace=True,
    )

    qc_keep = (
        (gex.obs["n_genes_by_counts"] >= MIN_GENES)
        & (gex.obs["n_genes_by_counts"] < MAX_GENES)
        & (gex.obs["total_counts"] < MAX_TOTAL_COUNTS)
        & (gex.obs["pct_counts_mt"] < MAX_PCT_MT)
    )

    gex = gex[qc_keep].copy()
    guide = guide[gex.obs_names].copy()
    n_after_expression_qc = int(gex.n_obs)

    if n_after_expression_qc == 0:
        raise ValueError(f"{sample_id}: no cells remain after RNA QC.")

    # Keep cells with zero or low guide counts at this stage. SCEPTRE needs the
    # full guide-count distribution among QC-passed cells for assignment.
    # Remove only guide features with no UMIs anywhere in this sample.
    guide_total_by_feature = np.asarray(guide.X.sum(axis=0)).ravel()
    nonzero_guide_feature = guide_total_by_feature > 0
    if int(nonzero_guide_feature.sum()) == 0:
        raise ValueError(
            f"{sample_id}: no Guide Capture feature has any UMI after RNA QC."
        )
    guide = guide[:, nonzero_guide_feature].copy()

    # Keep identical cell metadata in both modalities, including RNA QC metrics.
    guide.obs = gex.obs.copy()

    guide.var["grna_target"] = resolve_guide_targets(
        guide,
        fallback_metadata=SHARED_GUIDE_METADATA,
    ).reindex(guide.var_names).to_numpy()

    # Preserve raw counts in X and explicit count layers.
    gex.X = as_csr(gex.X)
    guide.X = as_csr(guide.X)
    gex.layers["counts"] = gex.X.copy()
    guide.layers["counts"] = guide.X.copy()

    raw_guide_umis = np.asarray(guide.X.sum(axis=1)).ravel()
    raw_detected_guides = np.asarray((guide.X > 0).sum(axis=1)).ravel()

    for obj in (gex, guide):
        obj.obs["sample_id"] = sample_id
        obj.obs["barcode"] = obj.obs_names.astype(str)
        obj.obs["passes_rna_qc"] = True

    # Guide-count summaries are useful before assignment and are copied to both modalities.
    for obj in (gex, guide):
        obj.obs["raw_total_guide_umis"] = raw_guide_umis
        obj.obs["raw_n_detected_guides"] = raw_detected_guides

    summary = {
        "sample_id": sample_id,
        "raw_cells": n_raw_cells,
        "after_expression_qc": n_after_expression_qc,
        "rna_qc_retention_rate": (
            float(n_after_expression_qc / n_raw_cells) if n_raw_cells else np.nan
        ),
        "rna_features": int(gex.n_vars),
        "guide_features_nonzero": int(guide.n_vars),
        "cells_with_any_raw_guide_umi": int((raw_guide_umis > 0).sum()),
        "fraction_with_any_raw_guide_umi": float((raw_guide_umis > 0).mean()),
    }
    return gex, guide, summary


def export_sample_for_sceptre(sample_id, cfg):
    gex, guide, summary = load_and_qc_sample(sample_id, cfg)
    sample_dir = OUTPUT_ROOT / sample_id
    transfer_dir = sample_dir / "sceptre_transfer"
    sample_dir.mkdir(parents=True, exist_ok=True)
    transfer_dir.mkdir(parents=True, exist_ok=True)

    # Save all RNA-QC-passed cells before SCEPTRE so the workflow can resume after a kernel restart.
    # Clear axis index names defensively before H5AD serialization.
    make_anndata_index_names_write_safe(gex)
    make_anndata_index_names_write_safe(guide)
    gex.write_h5ad(sample_dir / "rna_after_expression_qc_all_cells.h5ad", compression="gzip")
    guide.write_h5ad(sample_dir / "guide_after_expression_qc_all_cells.h5ad", compression="gzip")

    required_transfer_files = [
        transfer_dir / "response_matrix.mtx",
        transfer_dir / "grna_matrix.mtx",
        transfer_dir / "gene_names.txt",
        transfer_dir / "response_names.txt",
        transfer_dir / "guide_names.txt",
        transfer_dir / "cell_names.txt",
        transfer_dir / "guide_mapping.csv",
    ]
    should_write = OVERWRITE_TRANSFER_FILES or not all(path.exists() for path in required_transfer_files)

    if should_write:
        print(f"{sample_id}: writing SCEPTRE transfer files to {transfer_dir}")
        sio.mmwrite(transfer_dir / "response_matrix.mtx", as_csr(gex.X).T.tocoo())
        sio.mmwrite(transfer_dir / "grna_matrix.mtx", as_csr(guide.X).T.tocoo())
        write_lines(gex.var_names.astype(str), transfer_dir / "gene_names.txt")
        write_lines(gex.var["symbol"].astype(str), transfer_dir / "response_names.txt")
        write_lines(guide.var_names.astype(str), transfer_dir / "guide_names.txt")
        write_lines(gex.obs_names.astype(str), transfer_dir / "cell_names.txt")
        (
            guide.var[["grna_target"]]
            .assign(grna_id=guide.var_names.astype(str))
            [["grna_id", "grna_target"]]
            .to_csv(transfer_dir / "guide_mapping.csv", index=False)
        )
    else:
        print(f"{sample_id}: transfer files already exist; skipping Matrix Market rewrite.")

    with (sample_dir / "python_preprocessing_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2)

    return summary
```

## Code cell 5

```python
preprocessing_summaries = []
for sample_id, cfg in SAMPLES.items():
    preprocessing_summaries.append(export_sample_for_sceptre(sample_id, cfg))

preprocessing_summary_df = pd.DataFrame(preprocessing_summaries)
preprocessing_summary_df.to_csv(OUTPUT_ROOT / "preprocessing_summary.tsv", sep="\t", index=False)
display(preprocessing_summary_df)
```

### Text output (stream)

```text

=== IP_GFP: reading /oak/stanford/groups/engreitz/Users/tri/10xPerturbSeq/analysis/250507_ePAH1_NovaSeq/CRISPR_count/IP_GFP/filtered_feature_bc_matrix.h5 ===
reading /oak/stanford/groups/engreitz/Users/tri/10xPerturbSeq/analysis/250507_ePAH1_NovaSeq/CRISPR_count/IP_GFP/filtered_feature_bc_matrix.h5
 (0:00:02)

```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/anndata/_core/anndata.py:1758: UserWarning: Variable names are not unique. To make them unique, call `.var_names_make_unique`.
  utils.warn_names_duplicates("var")
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/anndata/_core/anndata.py:1758: UserWarning: Variable names are not unique. To make them unique, call `.var_names_make_unique`.
  utils.warn_names_duplicates("var")

```

### Text output (stream)

```text
IP_GFP: writing SCEPTRE transfer files to /oak/stanford/groups/engreitz/Users/tri/10xPerturbSeq/analysis/250507_ePAH1_NovaSeq/CRISPR_count/IP_GV_unified_sceptre_cnmf_no_celltype/IP_GFP/sceptre_transfer

=== GV_GFP: reading /oak/stanford/groups/engreitz/Users/tri/10xPerturbSeq/analysis/250507_ePAH1_NovaSeq/CRISPR_count/GV_GFP/filtered_feature_bc_matrix.h5 ===
reading /oak/stanford/groups/engreitz/Users/tri/10xPerturbSeq/analysis/250507_ePAH1_NovaSeq/CRISPR_count/GV_GFP/filtered_feature_bc_matrix.h5
 (0:00:02)

```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/anndata/_core/anndata.py:1758: UserWarning: Variable names are not unique. To make them unique, call `.var_names_make_unique`.
  utils.warn_names_duplicates("var")
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/anndata/_core/anndata.py:1758: UserWarning: Variable names are not unique. To make them unique, call `.var_names_make_unique`.
  utils.warn_names_duplicates("var")

```

### Text output (stream)

```text
GV_GFP: writing SCEPTRE transfer files to /oak/stanford/groups/engreitz/Users/tri/10xPerturbSeq/analysis/250507_ePAH1_NovaSeq/CRISPR_count/IP_GV_unified_sceptre_cnmf_no_celltype/GV_GFP/sceptre_transfer

```

### Text output (text/plain)

```text
  sample_id  raw_cells  after_expression_qc  rna_qc_retention_rate  \
0    IP_GFP      19508                19273               0.987954   
1    GV_GFP      14650                14475               0.988055   

   rna_features  guide_features_nonzero  cells_with_any_raw_guide_umi  \
0         33696                      96                         10094   
1         33696                      96                          4168   

   fraction_with_any_raw_guide_umi  
0                         0.523738  
1                         0.287945  
```

## Markdown cell 6

## Run SCEPTRE guide assignment through rpy2

## Code cell 7

```python
%load_ext rpy2.ipython

cpus_per_task = int(
    os.environ.get("SLURM_CPUS_PER_TASK")
    or os.environ.get("SLURM_CPUS_ON_NODE")
    or 1
)

SAMPLE_IDS_R = list(SAMPLES.keys())
SCEPTRE_OUT_ROOT_STR = str(OUTPUT_ROOT)
R_N_PROCESSORS = max(1, cpus_per_task - 2)
SCEPTRE_MOI_R = SCEPTRE_MOI
ASSIGN_GRNA_METHOD_R = ASSIGN_GRNA_METHOD
ASSIGN_PARALLEL_R = ASSIGN_PARALLEL

print("Samples:", SAMPLE_IDS_R)
print("SCEPTRE MOI:", SCEPTRE_MOI_R)
print("Assignment method:", ASSIGN_GRNA_METHOD_R)
print("Parallel assignment:", ASSIGN_PARALLEL_R)
print("R processors:", R_N_PROCESSORS)
```

### Text output (stream)

```text
Samples: ['IP_GFP', 'GV_GFP']
SCEPTRE MOI: high
Assignment method: mixture
Parallel assignment: False
R processors: 6

```

## Code cell 8

```python
%%R -i SCEPTRE_OUT_ROOT_STR -i SAMPLE_IDS_R -i R_N_PROCESSORS -i SCEPTRE_MOI_R -i ASSIGN_GRNA_METHOD_R -i ASSIGN_PARALLEL_R

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  BLIS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1",
  RCPP_PARALLEL_NUM_THREADS = "1",
  OMP_THREAD_LIMIT = "1"
)

suppressWarnings({
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
  }
})

if (!requireNamespace("sceptre", quietly = TRUE)) {
  stop("The R package `sceptre` is not installed in this Jupyter/R environment.")
}
if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("The R package `Matrix` is not installed.")
}

library(sceptre)
library(Matrix)

out_root <- as.character(SCEPTRE_OUT_ROOT_STR)
sample_ids <- as.character(SAMPLE_IDS_R)
n_processors <- as.integer(R_N_PROCESSORS)
moi_value <- as.character(SCEPTRE_MOI_R)
assign_method <- as.character(ASSIGN_GRNA_METHOD_R)
parallel_flag <- as.logical(ASSIGN_PARALLEL_R)

run_assignment <- function(sample_id) {
  sample_dir <- file.path(out_root, sample_id)
  transfer_dir <- file.path(sample_dir, "sceptre_transfer")
  result_dir <- file.path(sample_dir, "sceptre_assignment")
  log_dir <- file.path(result_dir, "logs")
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  cat("\n========================================\n")
  cat("SCEPTRE guide assignment:", sample_id, "\n")
  cat("========================================\n")

  response_mat <- Matrix::readMM(file.path(transfer_dir, "response_matrix.mtx"))
  grna_mat <- Matrix::readMM(file.path(transfer_dir, "grna_matrix.mtx"))
  response_mat <- as(response_mat, "dgCMatrix")
  grna_mat <- as(grna_mat, "dgCMatrix")

  gene_names <- readLines(file.path(transfer_dir, "gene_names.txt"))
  response_names <- readLines(file.path(transfer_dir, "response_names.txt"))
  guide_names <- readLines(file.path(transfer_dir, "guide_names.txt"))
  cell_names <- readLines(file.path(transfer_dir, "cell_names.txt"))
  guide_mapping <- read.csv(file.path(transfer_dir, "guide_mapping.csv"), stringsAsFactors = FALSE)

  stopifnot(
    nrow(response_mat) == length(gene_names),
    nrow(grna_mat) == length(guide_names),
    ncol(response_mat) == length(cell_names),
    ncol(grna_mat) == length(cell_names)
  )

  rownames(response_mat) <- gene_names
  colnames(response_mat) <- cell_names
  rownames(grna_mat) <- guide_names
  colnames(grna_mat) <- cell_names

  sceptre_object <- sceptre::import_data(
    response_matrix = response_mat,
    grna_matrix = grna_mat,
    grna_target_data_frame = guide_mapping,
    moi = moi_value,
    response_names = response_names
  )
  sceptre_object <- sceptre::set_analysis_parameters(sceptre_object)

  if (parallel_flag) {
    sceptre_object <- sceptre::assign_grnas(
      sceptre_object = sceptre_object,
      method = assign_method,
      parallel = TRUE,
      n_processors = n_processors,
      log_dir = log_dir
    )
  } else {
    sceptre_object <- sceptre::assign_grnas(
      sceptre_object = sceptre_object,
      method = assign_method,
      parallel = FALSE
    )
  }

  assignment_mat <- sceptre::get_grna_assignments(
    sceptre_object = sceptre_object,
    apply_cellwise_qc = FALSE
  )

  # get_grna_assignments() is guides x cells. Some SCEPTRE versions omit one
  # or both dimname vectors, so restore missing names from the exact matrices
  # supplied to import_data(). Do not merely call as.character(NULL), because
  # that would silently write an empty names file.
  expected_dim <- c(length(guide_names), length(cell_names))
  observed_dim <- dim(assignment_mat)
  if (!identical(as.integer(observed_dim), as.integer(expected_dim))) {
    stop(
      sprintf(
        paste0(
          "Unexpected assignment matrix dimensions for %s: got %d x %d, ",
          "expected %d guides x %d cells."
        ),
        sample_id,
        observed_dim[1], observed_dim[2],
        expected_dim[1], expected_dim[2]
      )
    )
  }

  assignment_guide_names <- rownames(assignment_mat)
  if (is.null(assignment_guide_names)) {
    assignment_guide_names <- as.character(guide_names)
    cat("Restored missing assignment guide names from guide_names.txt.\n")
  } else {
    assignment_guide_names <- as.character(assignment_guide_names)
    if (!setequal(assignment_guide_names, as.character(guide_names))) {
      stop("SCEPTRE assignment row names do not match the imported guide names.")
    }
  }

  assignment_cell_names <- colnames(assignment_mat)
  if (is.null(assignment_cell_names)) {
    assignment_cell_names <- as.character(cell_names)
    cat("Restored missing assignment cell names from cell_names.txt.\n")
  } else {
    assignment_cell_names <- as.character(assignment_cell_names)
    if (!setequal(assignment_cell_names, as.character(cell_names))) {
      stop("SCEPTRE assignment column names do not match the imported cell names.")
    }
  }

  dimnames(assignment_mat) <- list(
    assignment_guide_names,
    assignment_cell_names
  )

  # Save the expensive assignment result before performing file exports.
  saveRDS(
    sceptre_object,
    file.path(result_dir, "sceptre_assignment_object.rds")
  )

  # Export a numeric sparse matrix plus explicit names. Matrix Market drops dimnames.
  assignment_numeric <- Matrix::Matrix(assignment_mat * 1, sparse = TRUE)
  Matrix::writeMM(
    assignment_numeric,
    file.path(result_dir, "grna_assignment_matrix.mtx")
  )
  writeLines(
    assignment_guide_names,
    file.path(result_dir, "assignment_guide_names.txt")
  )
  writeLines(
    assignment_cell_names,
    file.path(result_dir, "assignment_cell_names.txt")
  )

  assignment_summary <- data.frame(
    barcode = assignment_cell_names,
    n_assigned_guides = as.integer(Matrix::colSums(assignment_mat))
  )
  write.csv(
    assignment_summary,
    file.path(result_dir, "assignment_cell_summary.csv"),
    row.names = FALSE
  )

  triplets <- summary(assignment_mat)
  if (nrow(triplets) > 0) {
    assignment_long <- data.frame(
      barcode = assignment_cell_names[triplets$j],
      grna_id = assignment_guide_names[triplets$i]
    )
  } else {
    assignment_long <- data.frame(barcode = character(0), grna_id = character(0))
  }
  write.csv(
    assignment_long,
    file.path(result_dir, "grna_assignments_long.csv"),
    row.names = FALSE
  )

  capture.output(sessionInfo(), file = file.path(result_dir, "R_sessionInfo.txt"))

  tryCatch({
    png(file.path(result_dir, "sceptre_assignment_plot.png"), width = 2400, height = 2000, res = 300)
    print(sceptre::plot_assign_grnas(sceptre_object, return_indiv_plots = TRUE))
    dev.off()
  }, error = function(e) {
    cat("Assignment plot failed:", conditionMessage(e), "\n")
    if (dev.cur() > 1) dev.off()
  })

  cat("Assigned cells:", sum(Matrix::colSums(assignment_mat) > 0), "/", ncol(assignment_mat), "\n")
  invisible(file.path(result_dir, "grna_assignment_matrix.mtx"))
}

assignment_files <- vapply(sample_ids, run_assignment, FUN.VALUE = character(1))
print(assignment_files)
```

### Text output (stream)

```text

========================================
SCEPTRE guide assignment: IP_GFP 
========================================
Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`.

Performing gRNA-to-cell assignments for gRNA Cdkn1a_2 (1 of 96)
Performing gRNA-to-cell assignments for gRNA Cdkn1a_C3 (5 of 96)
Performing gRNA-to-cell assignments for gRNA Notch1_C1 (10 of 96)
Performing gRNA-to-cell assignments for gRNA EIF2AK4_1 (15 of 96)
Performing gRNA-to-cell assignments for gRNA SMAD9_1 (20 of 96)
Performing gRNA-to-cell assignments for gRNA BMPR2_2 (25 of 96)
Performing gRNA-to-cell assignments for gRNA ENG_2 (30 of 96)
Performing gRNA-to-cell assignments for gRNA SOX17_2 (35 of 96)
Performing gRNA-to-cell assignments for gRNA ACVRL1_3 (40 of 96)
Performing gRNA-to-cell assignments for gRNA GDF2_3 (45 of 96)
Performing gRNA-to-cell assignments for gRNA TBX4_3 (50 of 96)
Performing gRNA-to-cell assignments for gRNA ATP13A3_4 (55 of 96)
Performing gRNA-to-cell assignments for gRNA KCNK3_4 (60 of 96)
Performing gRNA-to-cell assignments for gRNA SMAD4_4 (65 of 96)
Performing gRNA-to-cell assignments for gRNA non-targeting_4 (70 of 96)
Performing gRNA-to-cell assignments for gRNA non-targeting_9 (75 of 96)
Performing gRNA-to-cell assignments for gRNA non-targeting_14 (80 of 96)
Performing gRNA-to-cell assignments for gRNA non-targeting_19 (85 of 96)
Performing gRNA-to-cell assignments for gRNA non-targeting_24 (90 of 96)
Performing gRNA-to-cell assignments for gRNA non-targeting_29 (95 of 96)
Restored missing assignment cell names from cell_names.txt.
[[1]]

[[2]]

[[3]]

Assigned cells: 2849 / 19273 

========================================
SCEPTRE guide assignment: GV_GFP 
========================================
Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`.

Performing gRNA-to-cell assignments for gRNA Cdkn1a_2 (1 of 96)
Performing gRNA-to-cell assignments for gRNA Cdkn1a_C3 (5 of 96)
Performing gRNA-to-cell assignments for gRNA Notch1_C1 (10 of 96)
Performing gRNA-to-cell assignments for gRNA EIF2AK4_1 (15 of 96)
Performing gRNA-to-cell assignments for gRNA SMAD9_1 (20 of 96)
Performing gRNA-to-cell assignments for gRNA BMPR2_2 (25 of 96)
Performing gRNA-to-cell assignments for gRNA ENG_2 (30 of 96)
Performing gRNA-to-cell assignments for gRNA SOX17_2 (35 of 96)
Performing gRNA-to-cell assignments for gRNA ACVRL1_3 (40 of 96)
Performing gRNA-to-cell assignments for gRNA GDF2_3 (45 of 96)
Performing gRNA-to-cell assignments for gRNA TBX4_3 (50 of 96)
Performing gRNA-to-cell assignments for gRNA ATP13A3_4 (55 of 96)
Performing gRNA-to-cell assignments for gRNA KCNK3_4 (60 of 96)
Performing gRNA-to-cell assignments for gRNA SMAD4_4 (65 of 96)
Performing gRNA-to-cell assignments for gRNA non-targeting_4 (70 of 96)
Performing gRNA-to-cell assignments for gRNA non-targeting_9 (75 of 96)
Performing gRNA-to-cell assignments for gRNA non-targeting_14 (80 of 96)
Performing gRNA-to-cell assignments for gRNA non-targeting_19 (85 of 96)
Performing gRNA-to-cell assignments for gRNA non-targeting_24 (90 of 96)
Performing gRNA-to-cell assignments for gRNA non-targeting_29 (95 of 96)
Restored missing assignment cell names from cell_names.txt.
[[1]]

[[2]]

[[3]]

Assigned cells: 1246 / 14475 
                                                                                                                                                                                        IP_GFP 
"/oak/stanford/groups/engreitz/Users/tri/10xPerturbSeq/analysis/250507_ePAH1_NovaSeq/CRISPR_count/IP_GV_unified_sceptre_cnmf_no_celltype/IP_GFP/sceptre_assignment/grna_assignment_matrix.mtx" 
                                                                                                                                                                                        GV_GFP 
"/oak/stanford/groups/engreitz/Users/tri/10xPerturbSeq/analysis/250507_ePAH1_NovaSeq/CRISPR_count/IP_GV_unified_sceptre_cnmf_no_celltype/GV_GFP/sceptre_assignment/grna_assignment_matrix.mtx" 

```

### Text output (text/plain)

```text
'as(<dgTMatrix>, "dgCMatrix")' is deprecated.
Use 'as(., "CsparseMatrix")' instead.
See help("Deprecated") and help("Matrix-deprecated").

```

## Markdown cell 9

## Import assignments, filter guide-positive cells, combine mice, and write outputs

## Code cell 10

```python
def read_text_index(path):
    with Path(path).open() as handle:
        values = [line.rstrip("\n") for line in handle]
    index = pd.Index(values, dtype="object")
    if not index.is_unique:
        raise ValueError(f"Names in {path} are not unique.")
    return index


def row_labels_from_sparse(binary_matrix, names, unique=False):
    matrix = as_csr(binary_matrix)
    names = np.asarray(names, dtype=object)
    labels = []
    for row in range(matrix.shape[0]):
        row_names = names[matrix.indices[matrix.indptr[row]:matrix.indptr[row + 1]]]
        if unique:
            row_names = list(dict.fromkeys(row_names.tolist()))
        labels.append(";".join(map(str, row_names)))
    return labels


def add_guide_matrices_to_rna(rna, guide_counts, guide_assignments, assigned_guide_counts, guide):
    """Mirror guide- and target-level matrices into an RNA AnnData object."""
    guide_names = guide.var_names.astype(str).to_numpy()
    guide_targets = guide.var["grna_target"].astype(str).to_numpy()

    rna.obsm["guide_counts"] = as_csr(guide_counts).copy()
    rna.obsm["guide_assignments"] = as_csr(guide_assignments).copy()
    rna.obsm["assigned_guide_counts"] = as_csr(assigned_guide_counts).copy()
    rna.uns["guide_names"] = guide_names
    rna.uns["guide_targets"] = guide_targets

    target_names = pd.Index(pd.unique(guide_targets), dtype="object")
    target_codes = target_names.get_indexer(guide_targets)
    guide_to_target = sparse.csr_matrix(
        (
            np.ones(guide.n_vars, dtype=np.int8),
            (np.arange(guide.n_vars), target_codes),
        ),
        shape=(guide.n_vars, len(target_names)),
    )

    target_assignments = (as_csr(guide_assignments) @ guide_to_target).tocsr()
    if target_assignments.nnz:
        target_assignments.data[:] = 1
    target_assigned_counts = (
        as_csr(assigned_guide_counts) @ guide_to_target
    ).tocsr()

    rna.obsm["target_assignments"] = target_assignments
    rna.obsm["target_assigned_guide_counts"] = target_assigned_counts
    rna.uns["target_names"] = target_names.astype(str).to_numpy()


def attach_assignment_and_filter(sample_id):
    """
    Attach SCEPTRE assignments to all RNA-QC-passed cells, save that complete object,
    then return and save the subset with at least one assigned guide.
    """
    sample_dir = OUTPUT_ROOT / sample_id
    assignment_dir = sample_dir / "sceptre_assignment"

    gex = sc.read_h5ad(sample_dir / "rna_after_expression_qc_all_cells.h5ad")
    guide = sc.read_h5ad(sample_dir / "guide_after_expression_qc_all_cells.h5ad")
    gex.X = as_csr(gex.X)
    guide.X = as_csr(guide.X)

    # Ensure explicit raw-count layers are present after reload.
    gex.layers["counts"] = as_csr(gex.layers["counts"] if "counts" in gex.layers else gex.X)
    guide.layers["counts"] = as_csr(
        guide.layers["counts"] if "counts" in guide.layers else guide.X
    )

    assignment_guides = read_text_index(
        assignment_dir / "assignment_guide_names.txt"
    )
    assignment_cells = read_text_index(
        assignment_dir / "assignment_cell_names.txt"
    )
    assignment_gc = as_csr(
        sio.mmread(assignment_dir / "grna_assignment_matrix.mtx")
    )

    expected_shape = (len(assignment_guides), len(assignment_cells))
    if assignment_gc.shape != expected_shape:
        raise ValueError(
            f"{sample_id}: assignment matrix shape {assignment_gc.shape} does not "
            f"match {expected_shape[0]} guides x {expected_shape[1]} cells."
        )

    guide_order = assignment_guides.get_indexer(guide.var_names.astype(str))
    cell_order = assignment_cells.get_indexer(guide.obs_names.astype(str))
    if (guide_order < 0).any():
        missing = guide.var_names[np.flatnonzero(guide_order < 0)[:20]].tolist()
        raise ValueError(
            f"{sample_id}: guides missing from SCEPTRE output: {missing}"
        )
    if (cell_order < 0).any():
        missing = guide.obs_names[np.flatnonzero(cell_order < 0)[:20]].tolist()
        raise ValueError(
            f"{sample_id}: cells missing from SCEPTRE output: {missing}"
        )

    # SCEPTRE returns guides x cells. Reorder and transpose to cells x guides.
    assignment = (
        assignment_gc[guide_order, :][:, cell_order]
        .T
        .tocsr()
        .astype(np.int8)
    )
    if assignment.nnz:
        assignment.data[:] = 1

    raw_guide_counts = as_csr(guide.X)
    assigned_guide_counts = raw_guide_counts.multiply(assignment).tocsr()

    total_guide_umis = np.asarray(raw_guide_counts.sum(axis=1)).ravel()
    detected_guides = np.asarray((raw_guide_counts > 0).sum(axis=1)).ravel()
    n_assigned_guides = np.asarray(assignment.sum(axis=1)).ravel().astype(int)
    assigned_guide_umis = np.asarray(
        assigned_guide_counts.sum(axis=1)
    ).ravel()

    guide_names = guide.var_names.astype(str).to_numpy()
    guide_targets = guide.var["grna_target"].astype(str).to_numpy()
    assigned_guide_ids = row_labels_from_sparse(assignment, guide_names)
    assigned_targets = row_labels_from_sparse(
        assignment,
        guide_targets,
        unique=True,
    )

    n_assigned_targets = np.array(
        [0 if value == "" else len(value.split(";")) for value in assigned_targets],
        dtype=int,
    )
    guide_status = np.select(
        [n_assigned_guides == 0, n_assigned_guides == 1],
        ["no_assignment", "single_guide"],
        default="multiple_guides",
    )

    metrics = pd.DataFrame(
        {
            "total_guide_umis": total_guide_umis,
            "n_detected_guides": detected_guides,
            "n_assigned_guides": n_assigned_guides,
            "n_assigned_targets": n_assigned_targets,
            "assigned_guide_umis": assigned_guide_umis,
            "assigned_guide_ids": assigned_guide_ids,
            "assigned_targets": assigned_targets,
            "guide_assignment_status": guide_status,
            "has_sceptre_guide": n_assigned_guides > 0,
        },
        index=gex.obs_names,
    )

    for column in metrics.columns:
        gex.obs[column] = metrics[column]
        guide.obs[column] = metrics[column]

    guide.layers["sceptre_assignment"] = assignment
    guide.layers["assigned_counts"] = assigned_guide_counts
    add_guide_matrices_to_rna(
        gex,
        raw_guide_counts,
        assignment,
        assigned_guide_counts,
        guide,
    )

    # Save an enriched long assignment table.
    assignment_coo = assignment.tocoo()
    if assignment_coo.nnz:
        assigned_umis_long = np.asarray(
            raw_guide_counts[assignment_coo.row, assignment_coo.col]
        ).reshape(-1)
        long_df = pd.DataFrame(
            {
                "sample_id": sample_id,
                "barcode": gex.obs_names.to_numpy()[assignment_coo.row],
                "grna_id": guide_names[assignment_coo.col],
                "grna_target": guide_targets[assignment_coo.col],
                "guide_umi_count": assigned_umis_long,
            }
        )
    else:
        long_df = pd.DataFrame(
            columns=[
                "sample_id",
                "barcode",
                "grna_id",
                "grna_target",
                "guide_umi_count",
            ]
        )
    long_df.to_csv(sample_dir / "grna_assignments.csv", index=False)

    # Optionally preserve all RNA-QC-passed cells, including no-assignment cells.
    make_anndata_index_names_write_safe(gex)
    make_anndata_index_names_write_safe(guide)
    if SAVE_ALL_QC_WITH_ASSIGNMENTS:
        all_qc_mdata = MuData({"rna": gex, "guide": guide})
        all_qc_mdata.uns["description"] = (
            "All RNA-QC-passed cells; includes cells with and without a "
            "SCEPTRE-assigned guide."
        )
        all_qc_mdata.update()
        all_qc_mdata.write(
            sample_dir / f"{sample_id}_all_qc_with_sceptre.h5mu"
        )
        del all_qc_mdata

    # The cNMF branch keeps only cells with one or more SCEPTRE-assigned guides.
    keep = gex.obs["has_sceptre_guide"].to_numpy(dtype=bool)
    n_before = int(gex.n_obs)
    gex = gex[keep].copy()
    guide = guide[keep].copy()

    if gex.n_obs == 0:
        raise ValueError(
            f"{sample_id}: SCEPTRE assigned no guides to any RNA-QC-passed cell."
        )

    # Make cell IDs globally unique only after assignment.
    raw_barcodes = gex.obs_names.astype(str).to_numpy()
    global_ids = pd.Index(
        [f"{sample_id}:{barcode}" for barcode in raw_barcodes],
        name=None,
    )
    gex.obs_names = global_ids
    guide.obs_names = global_ids
    gex.obs["barcode"] = raw_barcodes
    guide.obs["barcode"] = raw_barcodes
    gex.obs["sample_id"] = sample_id
    guide.obs["sample_id"] = sample_id

    # Reattach guide matrices after subsetting and renaming.
    guide_counts_positive = as_csr(guide.X)
    guide_assignments_positive = as_csr(
        guide.layers["sceptre_assignment"]
    ).astype(np.int8)
    if guide_assignments_positive.nnz:
        guide_assignments_positive.data[:] = 1
    assigned_counts_positive = guide_counts_positive.multiply(
        guide_assignments_positive
    ).tocsr()

    guide.layers["counts"] = guide_counts_positive.copy()
    guide.layers["sceptre_assignment"] = guide_assignments_positive
    guide.layers["assigned_counts"] = assigned_counts_positive
    gex.layers["counts"] = as_csr(gex.X).copy()

    # AnnData subsets obsm automatically, but rebuild these matrices explicitly
    # to guarantee that their guide/target column order matches the guide modality.
    for key in [
        "guide_counts",
        "guide_assignments",
        "assigned_guide_counts",
        "target_assignments",
        "target_assigned_guide_counts",
    ]:
        if key in gex.obsm:
            del gex.obsm[key]
    for key in ["guide_names", "guide_targets", "target_names"]:
        if key in gex.uns:
            del gex.uns[key]

    add_guide_matrices_to_rna(
        gex,
        guide_counts_positive,
        guide_assignments_positive,
        assigned_counts_positive,
        guide,
    )

    sample_mdata = MuData({"rna": gex, "guide": guide})
    sample_mdata.uns["description"] = (
        "RNA-QC-passed cells with at least one SCEPTRE-assigned guide."
    )
    sample_mdata.update()
    sample_mdata.write(sample_dir / f"{sample_id}_guide_positive.h5mu")
    del sample_mdata

    # Remove sample-specific guide matrices before RNA concatenation. Their
    # column dimensions can differ when a guide has zero counts in one sample.
    # The combined matrices are rebuilt after the guide modalities are aligned.
    for key in [
        "guide_counts",
        "guide_assignments",
        "assigned_guide_counts",
        "target_assignments",
        "target_assigned_guide_counts",
    ]:
        if key in gex.obsm:
            del gex.obsm[key]
    for key in ["guide_names", "guide_targets", "target_names"]:
        if key in gex.uns:
            del gex.uns[key]

    gex.obs.to_csv(
        sample_dir / f"{sample_id}_guide_positive_cell_metadata.tsv",
        sep="\t",
        index=True,
    )

    summary = {
        "sample_id": sample_id,
        "cells_entering_sceptre": n_before,
        "cells_with_any_raw_guide_umi": int((total_guide_umis > 0).sum()),
        "cells_with_assigned_guide": int(gex.n_obs),
        "assignment_rate_all_qc": (
            float(gex.n_obs / n_before) if n_before else np.nan
        ),
        "assignment_rate_among_raw_guide_positive": (
            float(gex.n_obs / (total_guide_umis > 0).sum())
            if (total_guide_umis > 0).sum()
            else np.nan
        ),
        "single_guide_cells": int((gex.obs["n_assigned_guides"] == 1).sum()),
        "multiple_guide_cells": int((gex.obs["n_assigned_guides"] > 1).sum()),
        "guides_in_matrix": int(guide.n_vars),
    }
    return gex, guide, summary


sample_objects = {}
assignment_summaries = []
for sample_id in SAMPLES:
    gex_i, guide_i, summary_i = attach_assignment_and_filter(sample_id)
    sample_objects[sample_id] = {"rna": gex_i, "guide": guide_i}
    assignment_summaries.append(summary_i)

assignment_summary_df = pd.DataFrame(assignment_summaries)
assignment_summary_df.to_csv(OUTPUT_ROOT / "assignment_summary.tsv", sep="\t", index=False)
display(assignment_summary_df)
```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)

```

### Text output (text/plain)

```text
  sample_id  cells_entering_sceptre  cells_with_any_raw_guide_umi  \
0    IP_GFP                   19273                         10094   
1    GV_GFP                   14475                          4168   

   cells_with_assigned_guide  assignment_rate_all_qc  \
0                       2849                0.147823   
1                       1246                0.086079   

   assignment_rate_among_raw_guide_positive  single_guide_cells  \
0                                  0.282247                1300   
1                                  0.298944                 672   

   multiple_guide_cells  guides_in_matrix  
0                  1549                96  
1                   574                96  
```

## Code cell 11

```python
# Concatenate RNA and guide modalities separately. Cell IDs are already unique.
rna_combined = ad.concat(
    [sample_objects[sid]["rna"] for sid in SAMPLES],
    axis=0,
    join="outer",
    merge="first",
    uns_merge="first",
    index_unique=None,
    fill_value=0,
)
guide_combined = ad.concat(
    [sample_objects[sid]["guide"] for sid in SAMPLES],
    axis=0,
    join="outer",
    merge="first",
    uns_merge="first",
    index_unique=None,
    fill_value=0,
)

rna_combined.X = as_csr(rna_combined.X)
guide_combined.X = as_csr(guide_combined.X)

if not rna_combined.obs_names.is_unique:
    raise ValueError("Combined RNA cell IDs are not unique.")
if not guide_combined.obs_names.is_unique:
    raise ValueError("Combined guide cell IDs are not unique.")

# Force identical cell order across modalities.
guide_combined = guide_combined[rna_combined.obs_names].copy()

# Rebuild guide target metadata and explicitly check for cross-sample conflicts.
target_map = {}
for sample_id in SAMPLES:
    g = sample_objects[sample_id]["guide"]
    for grna_id, target in zip(g.var_names.astype(str), g.var["grna_target"].astype(str)):
        if grna_id in target_map and target_map[grna_id] != target:
            raise ValueError(
                f"Conflicting targets for guide {grna_id}: {target_map[grna_id]} vs {target}"
            )
        target_map[grna_id] = target

guide_combined.var["grna_id"] = guide_combined.var_names.astype(str)
guide_combined.var["grna_target"] = [target_map[x] for x in guide_combined.var_names.astype(str)]

# Filter RNA genes once on the combined guide-positive matrix, rather than separately per mouse.
n_cells_by_gene = np.asarray((rna_combined.X > 0).sum(axis=0)).ravel()
keep_genes = n_cells_by_gene >= MIN_CELLS_PER_GENE_COMBINED
rna_combined = rna_combined[:, keep_genes].copy()
rna_combined.var["n_cells_combined"] = n_cells_by_gene[keep_genes]
rna_combined.layers["counts"] = as_csr(rna_combined.X).copy()

# Standardize guide layers after outer concatenation.
guide_counts = as_csr(guide_combined.X)
guide_assignments = as_csr(
    guide_combined.layers["sceptre_assignment"]
).astype(np.int8)
if guide_assignments.nnz:
    guide_assignments.data[:] = 1
assigned_guide_counts = guide_counts.multiply(guide_assignments).tocsr()
guide_combined.layers["counts"] = guide_counts.copy()
guide_combined.layers["sceptre_assignment"] = guide_assignments
guide_combined.layers["assigned_counts"] = assigned_guide_counts

# Mirror guide-level matrices into RNA.obsm for compatibility with the prior cNMF evaluation workflow.
rna_combined.obsm["guide_counts"] = guide_counts.copy()
rna_combined.obsm["guide_assignments"] = guide_assignments.copy()
rna_combined.obsm["assigned_guide_counts"] = assigned_guide_counts.copy()
rna_combined.uns["guide_names"] = guide_combined.var_names.astype(str).to_numpy()
rna_combined.uns["guide_targets"] = guide_combined.var["grna_target"].astype(str).to_numpy()

# Also preserve target-level assignment and assigned-UMI matrices.
guide_targets = guide_combined.var["grna_target"].astype(str).to_numpy()
target_names = pd.Index(pd.unique(guide_targets))
target_codes = target_names.get_indexer(guide_targets)
guide_to_target = sparse.csr_matrix(
    (
        np.ones(guide_combined.n_vars, dtype=np.int8),
        (np.arange(guide_combined.n_vars), target_codes),
    ),
    shape=(guide_combined.n_vars, len(target_names)),
)
target_assignments = (guide_assignments @ guide_to_target).tocsr()
target_assignments.data[:] = 1
target_assigned_counts = (assigned_guide_counts @ guide_to_target).tocsr()

rna_combined.obsm["target_assignments"] = target_assignments
rna_combined.obsm["target_assigned_guide_counts"] = target_assigned_counts
rna_combined.uns["target_names"] = target_names.astype(str).to_numpy()

# Canonical multimodal object and a directly usable RNA h5ad for cNMF.
combined_dir = OUTPUT_ROOT / "combined"
combined_dir.mkdir(parents=True, exist_ok=True)

combined_h5mu = combined_dir / "IP_GFP_GV_GFP_guide_positive_raw_counts.h5mu"
cnmf_h5ad = combined_dir / "IP_GFP_GV_GFP_guide_positive_cNMF_input.h5ad"
guide_h5ad = combined_dir / "IP_GFP_GV_GFP_guide_modality.h5ad"

mdata_combined = MuData({"rna": rna_combined, "guide": guide_combined})
mdata_combined.uns["pipeline_parameters"] = {
    "min_genes": int(MIN_GENES),
    "max_genes": int(MAX_GENES),
    "max_total_counts": int(MAX_TOTAL_COUNTS),
    "max_pct_mt": float(MAX_PCT_MT),
    "cell_type_filter": "none",
    "sceptre_input_cells": "all_rna_qc_passed_cells",
    "saved_all_qc_with_assignments": bool(SAVE_ALL_QC_WITH_ASSIGNMENTS),
    "min_cells_per_gene_combined": int(MIN_CELLS_PER_GENE_COMBINED),
    "sceptre_moi": str(SCEPTRE_MOI),
    "sceptre_assignment_method": str(ASSIGN_GRNA_METHOD),
}
make_anndata_index_names_write_safe(rna_combined)
make_anndata_index_names_write_safe(guide_combined)
mdata_combined.update()
mdata_combined.write(combined_h5mu)
rna_combined.write_h5ad(cnmf_h5ad, compression="gzip")
guide_combined.write_h5ad(guide_h5ad, compression="gzip")

rna_combined.obs.to_csv(
    combined_dir / "IP_GFP_GV_GFP_guide_positive_cell_metadata.tsv",
    sep="\t",
    index=True,
)
guide_combined.var.to_csv(
    combined_dir / "IP_GFP_GV_GFP_guide_metadata.tsv",
    sep="\t",
    index=True,
)

combined_assignment_tables = []
for sample_id in SAMPLES:
    table_path = OUTPUT_ROOT / sample_id / "grna_assignments.csv"
    if table_path.exists():
        combined_assignment_tables.append(pd.read_csv(table_path))
if combined_assignment_tables:
    pd.concat(combined_assignment_tables, ignore_index=True).to_csv(
        combined_dir / "IP_GFP_GV_GFP_grna_assignments_long.csv",
        index=False,
    )

print("Saved canonical h5mu:", combined_h5mu)
print("Saved cNMF h5ad:", cnmf_h5ad)
print("Saved guide modality h5ad:", guide_h5ad)
print("Combined RNA shape:", rna_combined.shape)
print("Combined guide shape:", guide_combined.shape)
```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)

```

### Text output (stream)

```text
Saved canonical h5mu: /oak/stanford/groups/engreitz/Users/tri/10xPerturbSeq/analysis/250507_ePAH1_NovaSeq/CRISPR_count/IP_GV_unified_sceptre_cnmf_no_celltype/combined/IP_GFP_GV_GFP_guide_positive_raw_counts.h5mu
Saved cNMF h5ad: /oak/stanford/groups/engreitz/Users/tri/10xPerturbSeq/analysis/250507_ePAH1_NovaSeq/CRISPR_count/IP_GV_unified_sceptre_cnmf_no_celltype/combined/IP_GFP_GV_GFP_guide_positive_cNMF_input.h5ad
Saved guide modality h5ad: /oak/stanford/groups/engreitz/Users/tri/10xPerturbSeq/analysis/250507_ePAH1_NovaSeq/CRISPR_count/IP_GV_unified_sceptre_cnmf_no_celltype/combined/IP_GFP_GV_GFP_guide_modality.h5ad
Combined RNA shape: (4095, 16866)
Combined guide shape: (4095, 96)

```

## Markdown cell 12

## Final validation

## Code cell 13

```python
# Basic consistency checks before starting cNMF.
assert rna_combined.obs_names.equals(guide_combined.obs_names)
assert rna_combined.n_obs == guide_assignments.shape[0]
assert guide_combined.n_vars == guide_assignments.shape[1]
assert np.all(np.asarray(guide_assignments.sum(axis=1)).ravel() > 0)
assert np.all(rna_combined.obs["n_assigned_guides"].to_numpy() > 0)

validation = pd.DataFrame(
    {
        "n_cells": [rna_combined.n_obs],
        "n_rna_genes": [rna_combined.n_vars],
        "n_guides": [guide_combined.n_vars],
        "n_targets": [len(target_names)],
        "median_total_guide_umis": [float(np.median(rna_combined.obs["total_guide_umis"]))],
        "median_assigned_guide_umis": [float(np.median(rna_combined.obs["assigned_guide_umis"]))],
        "fraction_multiple_guides": [float((rna_combined.obs["n_assigned_guides"] > 1).mean())],
    }
)
validation.to_csv(OUTPUT_ROOT / "combined" / "validation_summary.tsv", sep="\t", index=False)
display(validation)

display(
    rna_combined.obs.groupby("sample_id", observed=True)
    .agg(
        cells=("sample_id", "size"),
        median_total_guide_umis=("total_guide_umis", "median"),
        median_assigned_guide_umis=("assigned_guide_umis", "median"),
        mean_assigned_guides=("n_assigned_guides", "mean"),
    )
)

print("RNA X contains raw counts. The combined object includes every RNA-QC-passed cell with at least one SCEPTRE-assigned guide; no cell-type filter was applied.")
```

### Text output (text/plain)

```text
   n_cells  n_rna_genes  n_guides  n_targets  median_total_guide_umis  \
0     4095        16866        96         17                     63.0   

   median_assigned_guide_umis  fraction_multiple_guides  
0                        60.0                  0.518437  
```

### Text output (text/plain)

```text
           cells  median_total_guide_umis  median_assigned_guide_umis  \
sample_id                                                               
GV_GFP      1246                     43.5                        41.0   
IP_GFP      2849                     72.0                        68.0   

           mean_assigned_guides  
sample_id                        
GV_GFP                 1.969502  
IP_GFP                 2.374868  
```

### Text output (stream)

```text
RNA X contains raw counts. The combined object includes every RNA-QC-passed cell with at least one SCEPTRE-assigned guide; no cell-type filter was applied.

```

## Code cell 14

```python
from pathlib import Path
import anndata as ad
import pandas as pd
import scipy.sparse as sp

INPUT_H5AD = Path(
    "/oak/stanford/groups/engreitz/Users/tri/10xPerturbSeq/"
    "analysis/250507_ePAH1_NovaSeq/CRISPR_count/"
    "IP_GV_unified_sceptre_cnmf_no_celltype/combined/"
    "IP_GFP_GV_GFP_guide_positive_cNMF_input.h5ad"
)

print("Input exists:", INPUT_H5AD.exists())
print("Input path:", INPUT_H5AD)

adata = ad.read_h5ad(INPUT_H5AD, backed="r")

print("\nShape")
print("cells:", adata.n_obs)
print("genes:", adata.n_vars)

print("\nX")
print("type:", type(adata.X))
print("dtype:", getattr(adata.X, "dtype", None))

print("\nobs columns")
print(adata.obs.columns.tolist())

for key in ["sample", "sample_id", "n_assigned_guides"]:
    if key in adata.obs:
        print(f"\nobs[{key!r}]")
        print(adata.obs[key].value_counts(dropna=False).head(20))

print("\nvar columns")
print(adata.var.columns.tolist())
print("first var names:", adata.var_names[:10].tolist())
print("unique var names:", adata.var_names.is_unique)

print("\nobsm")
for key in adata.obsm.keys():
    value = adata.obsm[key]
    print(
        key,
        "shape =", getattr(value, "shape", None),
        "sparse =", sp.issparse(value),
    )

print("\nuns")
print(list(adata.uns.keys()))

if "guide_names" in adata.uns:
    print("number of guides:", len(adata.uns["guide_names"]))
    print("first guides:", list(adata.uns["guide_names"][:10]))

if "guide_targets" in adata.uns:
    targets = pd.Series(adata.uns["guide_targets"], dtype=str)
    print("\nGuide target counts")
    print(targets.value_counts(dropna=False).head(40))

for key in ["guide_assignment", "guide_assignments"]:
    if key in adata.obsm:
        x = adata.obsm[key]
        print(
            f"\n{key}:",
            "shape =", x.shape,
            "sparse =", sp.issparse(x),
        )

adata.file.close()

print("\nPerturbNMF locations")
for repo in [
    Path("/oak/stanford/groups/engreitz/Users/ymo/Tools/PerturbNMF"),
    Path("/oak/stanford/groups/engreitz/Users/tri/Tools/PerturbNMF"),
]:
    print(repo, "exists =", repo.exists())
```

## Code cell 15

```python
from pathlib import Path
import anndata as ad
import pandas as pd
import scipy.sparse as sp

INPUT_H5AD = Path(
    "/oak/stanford/groups/engreitz/Users/tri/10xPerturbSeq/"
    "analysis/250507_ePAH1_NovaSeq/CRISPR_count/"
    "IP_GV_unified_sceptre_cnmf_no_celltype/combined/"
    "IP_GFP_GV_GFP_guide_positive_cNMF_input.h5ad"
)

print("Input exists:", INPUT_H5AD.exists())
print("Input path:", INPUT_H5AD)

adata = ad.read_h5ad(INPUT_H5AD, backed="r")

print("\nShape")
print("cells:", adata.n_obs)
print("genes:", adata.n_vars)

print("\nX")
print("type:", type(adata.X))
print("dtype:", getattr(adata.X, "dtype", None))

print("\nobs columns")
print(adata.obs.columns.tolist())

for key in ["sample", "sample_id", "n_assigned_guides"]:
    if key in adata.obs:
        print(f"\nobs[{key!r}]")
        print(adata.obs[key].value_counts(dropna=False).head(20))

print("\nvar columns")
print(adata.var.columns.tolist())
print("first var names:", adata.var_names[:10].tolist())
print("unique var names:", adata.var_names.is_unique)

print("\nobsm")
for key in adata.obsm.keys():
    value = adata.obsm[key]
    print(
        key,
        "shape =", getattr(value, "shape", None),
        "sparse =", sp.issparse(value),
    )

print("\nuns")
print(list(adata.uns.keys()))

if "guide_names" in adata.uns:
    print("number of guides:", len(adata.uns["guide_names"]))
    print("first guides:", list(adata.uns["guide_names"][:10]))

if "guide_targets" in adata.uns:
    targets = pd.Series(adata.uns["guide_targets"], dtype=str)
    print("\nGuide target counts")
    print(targets.value_counts(dropna=False).head(40))

for key in ["guide_assignment", "guide_assignments"]:
    if key in adata.obsm:
        x = adata.obsm[key]
        print(
            f"\n{key}:",
            "shape =", x.shape,
            "sparse =", sp.issparse(x),
        )

adata.file.close()

print("\nPerturbNMF locations")
for repo in [
    Path("/oak/stanford/groups/engreitz/Users/ymo/Tools/PerturbNMF"),
    Path("/oak/stanford/groups/engreitz/Users/tri/Tools/PerturbNMF"),
]:
    print(repo, "exists =", repo.exists())
```

### Text output (stream)

```text
Input exists: True
Input path: /oak/stanford/groups/engreitz/Users/tri/10xPerturbSeq/analysis/250507_ePAH1_NovaSeq/CRISPR_count/IP_GV_unified_sceptre_cnmf_no_celltype/combined/IP_GFP_GV_GFP_guide_positive_cNMF_input.h5ad

Shape
cells: 4095
genes: 16866

X
type: <class 'anndata._core.sparse_dataset._CSRDataset'>
dtype: float32

obs columns
['n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'total_counts_mt', 'log1p_total_counts_mt', 'pct_counts_mt', 'total_counts_ribo', 'log1p_total_counts_ribo', 'pct_counts_ribo', 'total_counts_hb', 'log1p_total_counts_hb', 'pct_counts_hb', 'sample_id', 'barcode', 'passes_rna_qc', 'raw_total_guide_umis', 'raw_n_detected_guides', 'total_guide_umis', 'n_detected_guides', 'n_assigned_guides', 'n_assigned_targets', 'assigned_guide_umis', 'assigned_guide_ids', 'assigned_targets', 'guide_assignment_status', 'has_sceptre_guide']

obs['sample_id']
sample_id
IP_GFP    2849
GV_GFP    1246
Name: count, dtype: int64

obs['n_assigned_guides']
n_assigned_guides
1     1972
2      849
3      526
4      300
5      182
6      126
7       65
8       41
9       22
10       9
12       2
11       1
Name: count, dtype: int64

var columns
['gene_ids', 'feature_types', 'genome', 'pattern', 'read', 'sequence', 'target_gene_id', 'target_gene_name', 'feature_name', 'symbol', 'mt', 'ribo', 'hb', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts', 'n_cells_combined']
first var names: ['Xkr4', 'Gm19938', 'Rp1', 'Sox17', 'Gm37587', 'Mrpl15', 'Lypla1', 'Tcea1', 'Rgs20', 'Atp6v1h']
unique var names: True

obsm
assigned_guide_counts shape = (4095, 96) sparse = True
guide_assignments shape = (4095, 96) sparse = True
guide_counts shape = (4095, 96) sparse = True
target_assigned_guide_counts shape = (4095, 17) sparse = True
target_assignments shape = (4095, 17) sparse = True

uns
['guide_names', 'guide_targets', 'target_names']
number of guides: 96
first guides: ['Cdkn1a_2', 'Cdkn1a_3', 'Cdkn1a_C1', 'Cdkn1a_C2', 'Cdkn1a_C3', 'Notch1_1', 'Notch1_2', 'Notch1_3', 'Notch1_4', 'Notch1_C1']

Guide target counts
non-targeting    30
Cdkn1a            5
Notch1            5
Acvrl1            4
Bmpr2             4
Cav1              4
Eif2ak4           4
Eng               4
Atp13a3           4
Gdf2              4
Kcnk3             4
Smad9             4
Kdr               4
Sox17             4
Tbx4              4
Smad4             4
Trp53             4
Name: count, dtype: int64

guide_assignments: shape = (4095, 96) sparse = True

PerturbNMF locations
/oak/stanford/groups/engreitz/Users/ymo/Tools/PerturbNMF exists = True
/oak/stanford/groups/engreitz/Users/tri/Tools/PerturbNMF exists = False

```

