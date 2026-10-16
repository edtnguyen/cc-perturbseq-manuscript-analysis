# Notebook mirror: `oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/Hippo_small_molecule.ipynb`

- Source archive path: `/oak/stanford/groups/engreitz/Users/tri/CodeBase/cc-perturbseq-figure-reproduction/recovery/sherlock_notebooks/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/Hippo_small_molecule.ipynb`
- SHA-256: `ab714e8811bd5e4e8408f6f4d23fb815c87a1bf946d9b324006e9449b7f92af9`
- Kernel language: `python`
- Full-fidelity notebook storage: Git LFS

## Code cell 0

```python
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

# Avoid nested CPU overuse inside R parallel workers.
os.environ["OMP_THREAD_LIMIT"] = "1"

print("Thread env vars set:")
for var in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS", "NUMEXPR_NUM_THREADS"]:
    print(var, os.environ.get(var))
```

### Text output (stream)

```text
Thread env vars set:
OMP_NUM_THREADS 1
OPENBLAS_NUM_THREADS 1
MKL_NUM_THREADS 1
BLIS_NUM_THREADS 1
NUMEXPR_NUM_THREADS 1

```

## Code cell 1

```python
# %% [markdown]
# # SCEPTRE trans-gene effects for Hippo perturbations
# 
# This notebook is a cleaned, parameterized version of the uploaded `ipsc_ec_sceptre.ipynb` workflow. Instead of testing only target-gene/cis pairs, it builds **trans discovery pairs**:
# 
# ```text
# Hippo perturbation target × all expressed response genes
# ```
# 
# For each requested day, it exports matrices from the `h5mu`, runs SCEPTRE through `rpy2`, and saves one result table per day plus a combined table across days.
# 
# Default Hippo targets are editable at the top. The default response gene set is all genes detected in at least `MIN_DETECTED_CELLS` cells, excluding the target gene itself.

# %%
%reset -f
import os
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import muon as mu
from scipy import sparse
import scipy.io as sio

warnings.filterwarnings("ignore")
print("Python OK")

# %%
# ---------------------------------------------------------------------
# User settings
# ---------------------------------------------------------------------
DATA_DIR = Path("/oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/2.preprocessing_tonys_pipeline/2.normalized_data/normalized_data")
WORKING_DIR = Path("/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks")
GUIDE_META_PATH = WORKING_DIR / "300genes_guide_metadata_v43.tsv"

# Where all SCEPTRE transfer files and results will go.
SCEPTRE_OUT_ROOT = Path("/scratch/users/edtng/sceptre_hippo_trans_gene_effects")
SCEPTRE_OUT_ROOT.mkdir(parents=True, exist_ok=True)

# Day h5mu files.
DAY_TO_H5MU = {
    "d0": DATA_DIR / "d0_normalized_gex_and_guide.h5mu",
    "d1": DATA_DIR / "d1_normalized_gex_and_guide.h5mu",
    "d2": DATA_DIR / "d2_normalized_gex_and_guide.h5mu",
    "d3": DATA_DIR / "d3_normalized_gex_and_guide.h5mu",
}

# Run all days by default. For a first test, set DAYS_TO_RUN = ["d2"]
DAYS_TO_RUN = ["d0", "d1", "d2", "d3"]

# Edit this if Jesse/Kyle want a different Hippo target set.
# These were the Hippo-like targets detected in the earlier cNMF workbook analyses.
HIPPO_TARGETS = ["LATS2", "NF2", "SAV1", "TAOK1", "WWC3", "AMOTL2"]

# Expression matrix source in the GEX modality.
RAW_COUNT_LAYER = "raw_counts"
GENE_SYMBOL_COL = "symbol"

# Response genes to test. This is important for runtime.
# Use 0 for truly all genes, but >= 20 or >= 50 is usually safer.
MIN_DETECTED_CELLS = 50
MIN_TOTAL_COUNTS = 0
EXCLUDE_SELF_TARGET = True

# Keep all cells by default, as in your original notebook.
FILTER_CELLS_WITH_NO_GUIDES = False

# Matrix export is slow; set True if input files changed.
OVERWRITE_TRANSFER_FILES = False

print("SCEPTRE_OUT_ROOT:", SCEPTRE_OUT_ROOT)
print("DAYS_TO_RUN:", DAYS_TO_RUN)
print("HIPPO_TARGETS:", HIPPO_TARGETS)

# %%
def _get_obs_col(obs, candidates, default=None):
    for c in candidates:
        if c in obs.columns:
            return obs[c]
    if default is not None:
        return pd.Series(default, index=obs.index)
    raise KeyError(f"None of these columns were found in obs: {candidates}")


def export_day_for_sceptre(day: str):
    """Export matrices and discovery pairs for one day."""
    h5mu_path = DAY_TO_H5MU[day]
    day_dir = SCEPTRE_OUT_ROOT / day
    day_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n=== {day}: loading {h5mu_path} ===")
    mdata = mu.read(h5mu_path)

    if "GEX" not in mdata.mod or "guide" not in mdata.mod:
        raise KeyError(f"Expected modalities `GEX` and `guide`; found {list(mdata.mod.keys())}")

    gex = mdata["GEX"]
    guide = mdata["guide"]

    if RAW_COUNT_LAYER not in gex.layers:
        raise KeyError(f"{RAW_COUNT_LAYER!r} not in GEX layers: {list(gex.layers.keys())}")
    if GENE_SYMBOL_COL not in gex.var.columns:
        raise KeyError(f"{GENE_SYMBOL_COL!r} not in GEX var columns: {list(gex.var.columns)}")

    if FILTER_CELLS_WITH_NO_GUIDES:
        if "guides_per_cell" not in gex.obs.columns:
            raise KeyError("FILTER_CELLS_WITH_NO_GUIDES=True but `guides_per_cell` missing from GEX.obs")
        keep_cells = gex.obs["guides_per_cell"].to_numpy() > 0
        mdata = mdata[keep_cells].copy()
        gex = mdata["GEX"]
        guide = mdata["guide"]
        print(f"Filtered cells with no guides; n cells now {gex.n_obs:,}")

    X_raw = gex.layers[RAW_COUNT_LAYER]
    if not sparse.issparse(X_raw):
        X_raw = sparse.csr_matrix(X_raw)

    # Response gene metadata and filtering.
    detected_cells = np.asarray((X_raw > 0).sum(axis=0)).ravel()
    total_counts = np.asarray(X_raw.sum(axis=0)).ravel()
    response_meta = pd.DataFrame({
        "response_id": gex.var_names.astype(str),
        "response_symbol": gex.var[GENE_SYMBOL_COL].astype(str).to_numpy(),
        "detected_cells": detected_cells,
        "total_counts": total_counts,
    })
    response_meta["keep_response"] = (
        (response_meta["detected_cells"] >= MIN_DETECTED_CELLS) &
        (response_meta["total_counts"] >= MIN_TOTAL_COUNTS)
    )
    response_meta.to_csv(day_dir / "response_gene_metadata.csv", index=False)

    kept_responses = response_meta.loc[response_meta["keep_response"], ["response_id", "response_symbol"]].copy()
    print(f"Response genes kept: {kept_responses.shape[0]:,} / {response_meta.shape[0]:,}")

    # Guide mapping from guide IDs to intended target names.
    guide_meta_df = pd.read_csv(GUIDE_META_PATH, sep="\t")
    required_meta_cols = {"guide_id", "intended_target_name"}
    missing = required_meta_cols - set(guide_meta_df.columns)
    if missing:
        raise KeyError(f"Guide metadata missing columns: {missing}")

    lookup = guide_meta_df.set_index("guide_id")
    guide_mapping = pd.DataFrame(index=guide.var_names.astype(str))
    guide_mapping.index.name = "grna_id"
    guide_mapping["grna_target"] = lookup.reindex(guide_mapping.index)["intended_target_name"].to_numpy()
    guide_mapping = guide_mapping.reset_index()

    if guide_mapping["grna_target"].isna().any():
        bad = guide_mapping.loc[guide_mapping["grna_target"].isna(), "grna_id"].head(10).tolist()
        raise ValueError(f"Some guides have no intended_target_name in metadata. Examples: {bad}")

    guide_mapping.to_csv(day_dir / "guide_mapping.csv", index=False)

    available_targets = sorted(set(guide_mapping["grna_target"].astype(str)))
    hippo_available = [t for t in HIPPO_TARGETS if t in available_targets]
    hippo_missing = [t for t in HIPPO_TARGETS if t not in available_targets]

    print("Hippo targets available in guide matrix:", hippo_available)
    if hippo_missing:
        print("WARNING: Hippo targets missing from guide matrix:", hippo_missing)

    # Discovery pairs: each Hippo target vs every kept response gene.
    pair_records = []
    for target in hippo_available:
        sub = kept_responses.copy()
        if EXCLUDE_SELF_TARGET:
            sub = sub[sub["response_symbol"].astype(str) != target]
        sub["grna_target"] = target
        pair_records.append(sub[["grna_target", "response_id", "response_symbol"]])

    if not pair_records:
        raise ValueError(f"No discovery pairs for {day}; no Hippo targets available.")

    discovery_pairs = pd.concat(pair_records, ignore_index=True)
    discovery_pairs.to_csv(day_dir / "discovery_pairs_hippo_vs_genes.csv", index=False)
    print(f"Discovery pairs: {discovery_pairs.shape[0]:,}")
    print(discovery_pairs.head())

    # Gene symbol map, useful for merging results after SCEPTRE.
    gene_symbol_map = response_meta[["response_id", "response_symbol", "detected_cells", "total_counts"]].copy()
    gene_symbol_map.to_csv(day_dir / "gene_symbol_map.csv", index=False)

    # Covariates. Use same logic as your original notebook, with fallbacks.
    obs = gex.obs
    covariates_df = pd.DataFrame(index=gex.obs_names.astype(str))
    covariates_df["batch"] = _get_obs_col(obs, ["sample_id", "batch", "pool", "stype"], default="batch1").astype(str).to_numpy()
    covariates_df["grna_n_nonzero"] = _get_obs_col(obs, ["guides_per_cell"], default=0).to_numpy()
    covariates_df["grna_n_umis"] = _get_obs_col(obs, ["guide_umi_counts"], default=0).to_numpy()
    covariates_df["response_n_nonzero"] = _get_obs_col(obs, ["n_genes_by_counts"], default=0).to_numpy()
    covariates_df["response_n_umis"] = _get_obs_col(obs, ["total_counts"], default=0).to_numpy()
    covariates_df["response_p_mito"] = _get_obs_col(obs, ["pct_counts_mt"], default=0).to_numpy()
    covariates_df.to_csv(day_dir / "covariates.csv")

    # Matrices and names. Matrix Market loses row/col names, so write names separately.
    need_write = OVERWRITE_TRANSFER_FILES or not (day_dir / "response_matrix.mtx").exists()
    if need_write:
        print("Writing Matrix Market files...")
        sio.mmwrite(day_dir / "response_matrix.mtx", X_raw.T)        # genes x cells
        sio.mmwrite(day_dir / "grna_matrix.mtx", guide.X.T)          # guides x cells
        pd.Series(gex.var_names.astype(str)).to_csv(day_dir / "gene_names.csv", index=False, header=False)
        pd.Series(guide.var_names.astype(str)).to_csv(day_dir / "guide_names.csv", index=False, header=False)
        pd.Series(gex.obs_names.astype(str)).to_csv(day_dir / "cell_names.csv", index=False, header=False)
    else:
        print("Matrix Market files already exist; skipping rewrite.")

    return {
        "day": day,
        "n_cells": gex.n_obs,
        "n_genes": gex.n_vars,
        "n_kept_response_genes": kept_responses.shape[0],
        "n_discovery_pairs": discovery_pairs.shape[0],
        "hippo_available": hippo_available,
        "hippo_missing": hippo_missing,
    }


export_summaries = []
for day in DAYS_TO_RUN:
    export_summaries.append(export_day_for_sceptre(day))

export_summary_df = pd.DataFrame(export_summaries)
export_summary_df.to_csv(SCEPTRE_OUT_ROOT / "python_export_summary.tsv", sep="\t", index=False)
display(export_summary_df)
```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/tqdm/auto.py:21: TqdmWarning: IProgress not found. Please update jupyter and ipywidgets. See https://ipywidgets.readthedocs.io/en/stable/user_install.html
  from .autonotebook import tqdm as notebook_tqdm
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/muon/_core/preproc.py:31: FutureWarning: `__version__` is deprecated, use `importlib.metadata.version('scanpy')` instead
  if Version(scanpy.__version__) < Version("1.10"):

```

### Text output (stream)

```text
Python OK
SCEPTRE_OUT_ROOT: /scratch/users/edtng/sceptre_hippo_trans_gene_effects
DAYS_TO_RUN: ['d0', 'd1', 'd2', 'd3']
HIPPO_TARGETS: ['LATS2', 'NF2', 'SAV1', 'TAOK1', 'WWC3', 'AMOTL2']

=== d0: loading /oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/2.preprocessing_tonys_pipeline/2.normalized_data/normalized_data/d0_normalized_gex_and_guide.h5mu ===
Response genes kept: 16,392 / 16,392
Hippo targets available in guide matrix: ['LATS2', 'NF2', 'SAV1', 'TAOK1', 'WWC3', 'AMOTL2']
Discovery pairs: 98,346
  grna_target         response_id  response_symbol
0       LATS2   ENSG00000227232.5           WASH7P
1       LATS2   ENSG00000238009.6  ENSG00000238009
2       LATS2   ENSG00000241860.7  ENSG00000241860
3       LATS2   ENSG00000279457.4           WASH9P
4       LATS2  ENSG00000228463.10  ENSG00000228463
Matrix Market files already exist; skipping rewrite.

=== d1: loading /oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/2.preprocessing_tonys_pipeline/2.normalized_data/normalized_data/d1_normalized_gex_and_guide.h5mu ===
Response genes kept: 13,383 / 13,383
Hippo targets available in guide matrix: ['LATS2', 'NF2', 'SAV1', 'TAOK1', 'WWC3', 'AMOTL2']
Discovery pairs: 80,292
  grna_target        response_id  response_symbol
0       LATS2  ENSG00000227232.5           WASH7P
1       LATS2  ENSG00000238009.6  ENSG00000238009
2       LATS2  ENSG00000268903.1  ENSG00000268903
3       LATS2  ENSG00000241860.7  ENSG00000241860
4       LATS2  ENSG00000279457.4           WASH9P
Matrix Market files already exist; skipping rewrite.

=== d2: loading /oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/2.preprocessing_tonys_pipeline/2.normalized_data/normalized_data/d2_normalized_gex_and_guide.h5mu ===
Response genes kept: 12,732 / 12,732
Hippo targets available in guide matrix: ['LATS2', 'NF2', 'SAV1', 'TAOK1', 'WWC3', 'AMOTL2']
Discovery pairs: 76,386
  grna_target         response_id  response_symbol
0       LATS2   ENSG00000227232.5           WASH7P
1       LATS2   ENSG00000238009.6  ENSG00000238009
2       LATS2   ENSG00000241860.7  ENSG00000241860
3       LATS2   ENSG00000279457.4           WASH9P
4       LATS2  ENSG00000228463.10  ENSG00000228463
Matrix Market files already exist; skipping rewrite.

=== d3: loading /oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/2.preprocessing_tonys_pipeline/2.normalized_data/normalized_data/d3_normalized_gex_and_guide.h5mu ===
Response genes kept: 14,712 / 14,712
Hippo targets available in guide matrix: ['LATS2', 'NF2', 'SAV1', 'TAOK1', 'WWC3', 'AMOTL2']
Discovery pairs: 88,266
  grna_target         response_id  response_symbol
0       LATS2   ENSG00000227232.5           WASH7P
1       LATS2   ENSG00000238009.6  ENSG00000238009
2       LATS2   ENSG00000241860.7  ENSG00000241860
3       LATS2   ENSG00000279457.4           WASH9P
4       LATS2  ENSG00000228463.10  ENSG00000228463
Matrix Market files already exist; skipping rewrite.

```

### Text output (text/plain)

```text
  day  n_cells  n_genes  n_kept_response_genes  n_discovery_pairs  \
0  d0    38609    16392                  16392              98346   
1  d1    16925    13383                  13383              80292   
2  d2    21980    12732                  12732              76386   
3  d3    44102    14712                  14712              88266   

                           hippo_available hippo_missing  
0  [LATS2, NF2, SAV1, TAOK1, WWC3, AMOTL2]            []  
1  [LATS2, NF2, SAV1, TAOK1, WWC3, AMOTL2]            []  
2  [LATS2, NF2, SAV1, TAOK1, WWC3, AMOTL2]            []  
3  [LATS2, NF2, SAV1, TAOK1, WWC3, AMOTL2]            []  
```

## Code cell 2

```python
# %%
# Cell 2 — load rpy2 after thread env vars are set.
%load_ext rpy2.ipython

# %%
# Cell 3 — safe SCEPTRE run settings.
# Use the SLURM allocation but leave some overhead for code-server/Jupyter.

import os

# Thread caps. These should already be set after kernel restart,
# but setting them again here is harmless.
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

cpus_per_task = int(
    os.environ.get("SLURM_CPUS_PER_TASK")
    or os.environ.get("SLURM_CPUS_ON_NODE")
    or 1
)

SCEPTRE_OUT_ROOT_STR = str(SCEPTRE_OUT_ROOT)

# Full D2 run first
DAYS_TO_RUN_R = ["d2"]

# Conservative parallel settings
R_PARALLEL = True
R_N_PROCESSORS = cpus_per_task - 2

ASSIGN_GRNA_METHOD = "mixture"
P_MITO_THRESHOLD = 0.25

# Full discovery, no calibration for first full pass
RUN_CALIBRATION_CHECK = False
MAX_DISCOVERY_PAIRS_FOR_TEST = 0

print("SLURM_CPUS_PER_TASK:", os.environ.get("SLURM_CPUS_PER_TASK"))
print("SLURM_CPUS_ON_NODE:", os.environ.get("SLURM_CPUS_ON_NODE"))
print("R_N_PROCESSORS:", R_N_PROCESSORS)
print("DAYS_TO_RUN_R:", DAYS_TO_RUN_R)
print("RUN_CALIBRATION_CHECK:", RUN_CALIBRATION_CHECK)
print("MAX_DISCOVERY_PAIRS_FOR_TEST:", MAX_DISCOVERY_PAIRS_FOR_TEST)

```

### Text output (stream)

```text
SLURM_CPUS_PER_TASK: 16
SLURM_CPUS_ON_NODE: 16
R_N_PROCESSORS: 14
DAYS_TO_RUN_R: ['d2']
RUN_CALIBRATION_CHECK: False
MAX_DISCOVERY_PAIRS_FOR_TEST: 0

```

## Code cell 3

```python
%%R -i SCEPTRE_OUT_ROOT_STR -i DAYS_TO_RUN_R -i R_PARALLEL -i R_N_PROCESSORS -i ASSIGN_GRNA_METHOD -i P_MITO_THRESHOLD -i RUN_CALIBRATION_CHECK -i MAX_DISCOVERY_PAIRS_FOR_TEST

# Cell 4 — define a parallel-safe SCEPTRE runner.
# IMPORTANT: %%R must be the first line of this cell.

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

library(sceptre)
library(Matrix)

out_root <- as.character(SCEPTRE_OUT_ROOT_STR)
days <- as.character(DAYS_TO_RUN_R)
parallel_flag <- as.logical(R_PARALLEL)
n_processors <- as.integer(R_N_PROCESSORS)
assign_method <- as.character(ASSIGN_GRNA_METHOD)
p_mito_threshold <- as.numeric(P_MITO_THRESHOLD)
run_calibration <- as.logical(RUN_CALIBRATION_CHECK)
max_pairs <- as.integer(MAX_DISCOVERY_PAIRS_FOR_TEST)

cat("R parallel flag:", parallel_flag, "\n")
cat("R n_processors:", n_processors, "\n")
cat("Thread caps:\n")
print(Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS", "NUMEXPR_NUM_THREADS")))

run_sceptre_day <- function(day) {
  data_dir <- file.path(out_root, day)
  result_dir <- file.path(data_dir, "sceptre_results")
  log_dir <- file.path(data_dir, "sceptre_logs")
  dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

  cat("\n============================\n")
  cat("Running SCEPTRE for day:", day, "\n")
  cat("data_dir:", data_dir, "\n")
  cat("============================\n")

  response_mat <- readMM(file.path(data_dir, "response_matrix.mtx"))
  grna_mat <- readMM(file.path(data_dir, "grna_matrix.mtx"))

  gene_names <- read.csv(file.path(data_dir, "gene_names.csv"), header = FALSE)$V1
  guide_names <- read.csv(file.path(data_dir, "guide_names.csv"), header = FALSE)$V1
  cell_names <- read.csv(file.path(data_dir, "cell_names.csv"), header = FALSE)$V1

  rownames(response_mat) <- gene_names
  colnames(response_mat) <- cell_names
  rownames(grna_mat) <- guide_names
  colnames(grna_mat) <- cell_names

  covariates <- read.csv(file.path(data_dir, "covariates.csv"), row.names = 1)
  guide_mapping <- read.csv(file.path(data_dir, "guide_mapping.csv"))
  discovery_pairs_full <- read.csv(file.path(data_dir, "discovery_pairs_hippo_vs_genes.csv"))

  if (max_pairs > 0 && nrow(discovery_pairs_full) > max_pairs) {
    cat("SMOKE TEST: keeping first", max_pairs, "pairs out of", nrow(discovery_pairs_full), "\n")
    discovery_pairs_full <- discovery_pairs_full[seq_len(max_pairs), , drop = FALSE]
  }

  discovery_pairs <- unique(discovery_pairs_full[, c("grna_target", "response_id")])
  cat("Discovery pairs:", nrow(discovery_pairs), "\n")
  cat("Targets:", paste(sort(unique(discovery_pairs$grna_target)), collapse = ", "), "\n")

  if ("response_p_mito" %in% colnames(covariates)) {
    names(covariates)[names(covariates) == "response_p_mito"] <- "mt_percent"
  }
  reserved_cols <- c("response_n_nonzero", "response_n_umis", "grna_n_nonzero", "grna_n_umis")
  covariates <- covariates[, !(colnames(covariates) %in% reserved_cols), drop = FALSE]

  sceptre_object <- import_data(
    response_matrix = response_mat,
    grna_matrix = grna_mat,
    grna_target_data_frame = guide_mapping,
    extra_covariates = covariates,
    moi = "high"
  )

  sceptre_object <- set_analysis_parameters(
    sceptre_object = sceptre_object,
    discovery_pairs = discovery_pairs,
    side = "both",
    control_group = "complement"
  )

  # Keep this serial. It is usually not the bottleneck, and parallelizing it inside rpy2 can oversubscribe.
  sceptre_object <- assign_grnas(
    sceptre_object = sceptre_object,
    method = assign_method,
    parallel = FALSE
  )

  sceptre_object <- run_qc(
    sceptre_object = sceptre_object,
    p_mito_threshold = p_mito_threshold
  )

  tryCatch({
    png(file.path(result_dir, "sceptre_QC_plot.png"), width = 2400, height = 2000, res = 300)
    print(plot(sceptre_object))
    dev.off()
  }, error = function(e) {
    cat("QC plot failed:", conditionMessage(e), "\n")
    if (dev.cur() > 1) dev.off()
  })

  if (run_calibration) {
    cat("Running calibration check...\n")
    sceptre_object <- run_calibration_check(
      sceptre_object = sceptre_object,
      parallel = parallel_flag,
      n_processors = n_processors,
      log_dir = log_dir
    )
    calibration_result <- get_result(sceptre_object, analysis = "run_calibration_check")
    write.csv(calibration_result, file.path(result_dir, "calibration_check_result.csv"), row.names = FALSE)
  } else {
    cat("Skipping calibration check.\n")
  }

  cat("Running discovery analysis...\n")
  sceptre_object <- run_discovery_analysis(
    sceptre_object = sceptre_object,
    output_amount = 3,
    parallel = parallel_flag,
    n_processors = n_processors,
    log_dir = log_dir
  )

  discovery_result <- get_result(sceptre_object, analysis = "run_discovery_analysis")
  discovery_result$day <- day

  gene_map_path <- file.path(data_dir, "gene_symbol_map.csv")
  if (file.exists(gene_map_path)) {
    gene_map <- read.csv(gene_map_path)
    discovery_result <- merge(discovery_result, gene_map, by = "response_id", all.x = TRUE)
  }

  out_csv <- file.path(result_dir, paste0("sceptre_hippo_vs_genes_", day, ".csv.gz"))
  con <- gzfile(out_csv, "wt")
  write.csv(discovery_result, con, row.names = FALSE)
  close(con)

  cat("Saved discovery result to:", out_csv, "\n")
  invisible(out_csv)
}

cat("Defined run_sceptre_day().\n")
cat("Days to run:", paste(days, collapse = ", "), "\n")



```

### Text output (stream)

```text
R parallel flag: TRUE 
R n_processors: 14 
Thread caps:
     OMP_NUM_THREADS OPENBLAS_NUM_THREADS      MKL_NUM_THREADS 
                 "1"                  "1"                  "1" 
    BLIS_NUM_THREADS  NUMEXPR_NUM_THREADS 
                 "1"                  "1" 
Defined run_sceptre_day().
Days to run: d2 

```

## Code cell 4

```python
%%R
# Cell 5 — run the analysis after the function is defined.
# Start with DAYS_TO_RUN_R = c("d2") and MAX_DISCOVERY_PAIRS_FOR_TEST = 1000.
result_files <- c()
for (day in days) {
  result_files <- c(result_files, run_sceptre_day(day))
}
cat("\nAll result files:\n")
print(result_files)
```

### Text output (stream)

```text

============================
Running SCEPTRE for day: d2 
data_dir: /scratch/users/edtng/sceptre_hippo_trans_gene_effects/d2 
============================
Discovery pairs: 76386 
Targets: AMOTL2, LATS2, NF2, SAV1, TAOK1, WWC3 
Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`.

Performing gRNA-to-cell assignments for gRNA ACAA1_+_38178488.23-P1P2 (1 of 1926)
Performing gRNA-to-cell assignments for gRNA ACAA1_+_38178620.23-P1P2 (5 of 1926)
Performing gRNA-to-cell assignments for gRNA ACTL7B_-_111617995.23-P1P2 (10 of 1926)
Performing gRNA-to-cell assignments for gRNA ACVR1B_+_52345762.23-P1P2 (15 of 1926)
Performing gRNA-to-cell assignments for gRNA ACVR2A_+_148602343.23-P1P2 (20 of 1926)
Performing gRNA-to-cell assignments for gRNA ACVR2B_+_38495339.23-P1P2 (25 of 1926)
Performing gRNA-to-cell assignments for gRNA ACVR2B_-_38495819.23-P1P2 (30 of 1926)
Performing gRNA-to-cell assignments for gRNA ADCY2_-_7396429.23-P1P2 (35 of 1926)
Performing gRNA-to-cell assignments for gRNA AMOTL2_+_134093279.23-P1P2 (40 of 1926)
Performing gRNA-to-cell assignments for gRNA ANP32E_+_150208342.23-P1P2 (45 of 1926)
Performing gRNA-to-cell assignments for gRNA APAF1_-_99039071.23-P1P2 (50 of 1926)
Performing gRNA-to-cell assignments for gRNA ARIH1_+_72766725.23-P1P2 (55 of 1926)
Performing gRNA-to-cell assignments for gRNA ARIH1_-_72766755.23-P1P2 (60 of 1926)
Performing gRNA-to-cell assignments for gRNA ARTN_-_44399429.23-P1 (65 of 1926)
Performing gRNA-to-cell assignments for gRNA ASCC3_-_101329193.23-P1P2 (70 of 1926)
Performing gRNA-to-cell assignments for gRNA ASH2L_-_37963059.23-P1P2 (75 of 1926)
Performing gRNA-to-cell assignments for gRNA ATMIN_+_81069661.23-P1P2 (80 of 1926)
Performing gRNA-to-cell assignments for gRNA ATP5F1_+_111992112.23-P1P2 (85 of 1926)
Performing gRNA-to-cell assignments for gRNA ATP5F1_-_111992129.23-P1P2 (90 of 1926)
Performing gRNA-to-cell assignments for gRNA ATP6AP1_-_153657075.23-P1P2 (95 of 1926)
Performing gRNA-to-cell assignments for gRNA ATXN2L_-_28834668.23-P1P2 (100 of 1926)
Performing gRNA-to-cell assignments for gRNA BAD_+_64052102.23-P1P2 (105 of 1926)
Performing gRNA-to-cell assignments for gRNA BAG6_+_31620398.23-P1P2 (110 of 1926)
Performing gRNA-to-cell assignments for gRNA BAX_+_49458268.23-P1P2 (115 of 1926)
Performing gRNA-to-cell assignments for gRNA BAX_-_49458323.23-P1P2 (120 of 1926)
Performing gRNA-to-cell assignments for gRNA BCAS2_-_115124045.23-P1P2 (125 of 1926)
Performing gRNA-to-cell assignments for gRNA BMPR1A_-_88516386.23-P1P2 (130 of 1926)
Performing gRNA-to-cell assignments for gRNA BRD3_-_136932911.23-P1P2 (135 of 1926)
Performing gRNA-to-cell assignments for gRNA C17orf49_+_6918114.23-P2 (140 of 1926)
Performing gRNA-to-cell assignments for gRNA C5orf22_+_31532436.23-P1P2 (145 of 1926)
Performing gRNA-to-cell assignments for gRNA C5orf22_-_31532650.23-P1P2 (150 of 1926)
Performing gRNA-to-cell assignments for gRNA CACFD1_-_136325160.23-P1P2 (155 of 1926)
Performing gRNA-to-cell assignments for gRNA CACNG1_+_65040756.23-P1P2 (160 of 1926)
Performing gRNA-to-cell assignments for gRNA CAPN5_+_76778225.23-P1P2 (165 of 1926)
Performing gRNA-to-cell assignments for gRNA CASP3_+_185570560.23-P1P2 (170 of 1926)
Performing gRNA-to-cell assignments for gRNA CBLL1_+_107384793.23-P1P2 (175 of 1926)
Performing gRNA-to-cell assignments for gRNA CBLL1_-_107384635.23-P1P2 (180 of 1926)
Performing gRNA-to-cell assignments for gRNA CCM2_-_45039941.23-P1P2 (185 of 1926)
Performing gRNA-to-cell assignments for gRNA CCNC_-_100016267.23-P1P2 (190 of 1926)
Performing gRNA-to-cell assignments for gRNA CCNT2_-_135676416.23-P1P2 (195 of 1926)
Performing gRNA-to-cell assignments for gRNA CCT6B_+_33288381.23-P1P2 (200 of 1926)
Performing gRNA-to-cell assignments for gRNA CD2BP2_+_30366413.23-P1P2 (205 of 1926)
Performing gRNA-to-cell assignments for gRNA CD2BP2_-_30366630.23-P1P2 (210 of 1926)
Performing gRNA-to-cell assignments for gRNA CD99L2_-_150067120.23-P1P2 (215 of 1926)
Performing gRNA-to-cell assignments for gRNA CDH5_-_66400630.23-P1P2 (220 of 1926)
Performing gRNA-to-cell assignments for gRNA CDKN1A_-_36646538.23-P1P2 (225 of 1926)
Performing gRNA-to-cell assignments for gRNA CENPC_+_68411190.23-P1P2 (230 of 1926)
Performing gRNA-to-cell assignments for gRNA C16orf93_+_30773828.23-P1P2 (235 of 1926)
Performing gRNA-to-cell assignments for gRNA C16orf93_-_30773391.23-P1P2 (240 of 1926)
Performing gRNA-to-cell assignments for gRNA CHMP3_-_86790451.23-P1P2 (245 of 1926)
Performing gRNA-to-cell assignments for gRNA CHORDC1_-_89955894.23-P1P2 (250 of 1926)
Performing gRNA-to-cell assignments for gRNA CNEP1R1_+_50059284.23-P1P2 (255 of 1926)
Performing gRNA-to-cell assignments for gRNA CNOT2_+_70637651.23-P1P2 (260 of 1926)
Performing gRNA-to-cell assignments for gRNA COPS4_+_83956339.23-P1P2 (265 of 1926)
Performing gRNA-to-cell assignments for gRNA COPS4_-_83956572.23-P1P2 (270 of 1926)
Performing gRNA-to-cell assignments for gRNA COPS6_-_99686676.23-P1P2 (275 of 1926)
Performing gRNA-to-cell assignments for gRNA CSNK2B_-_31634088.23-P1P2 (280 of 1926)
Performing gRNA-to-cell assignments for gRNA CTR9_+_10772881.23-P1P2 (285 of 1926)
Performing gRNA-to-cell assignments for gRNA CUL3_-_225449667.23-P1P2 (290 of 1926)
Performing gRNA-to-cell assignments for gRNA CYC1_+_145150267.23-P1P2 (295 of 1926)
Performing gRNA-to-cell assignments for gRNA CYC1_-_145149996.23-P1P2 (300 of 1926)
Performing gRNA-to-cell assignments for gRNA CYCS_+_25164869.23-P1P2 (305 of 1926)
Performing gRNA-to-cell assignments for gRNA DACH1_-_72441101.23-P1P2 (310 of 1926)
Performing gRNA-to-cell assignments for gRNA DBR1_+_137893673.23-P1P2 (315 of 1926)
Performing gRNA-to-cell assignments for gRNA DCTN2_-_57940910.23-P1P2 (320 of 1926)
Performing gRNA-to-cell assignments for gRNA DDX42_+_61851532.23-P1P2 (325 of 1926)
Performing gRNA-to-cell assignments for gRNA DDX42_-_61851812.23-P1P2 (330 of 1926)
Performing gRNA-to-cell assignments for gRNA DGCR8_-_20068158.23-P1P2 (335 of 1926)
Performing gRNA-to-cell assignments for gRNA DNAJC17_-_41099215.23-P1P2 (340 of 1926)
Performing gRNA-to-cell assignments for gRNA DNM2_+_10829157.23-P1P2 (345 of 1926)
Performing gRNA-to-cell assignments for gRNA DPM1_+_49575018.23-P1P2 (350 of 1926)
Performing gRNA-to-cell assignments for gRNA DPM2_+_130700056.23-P1P2 (355 of 1926)
Performing gRNA-to-cell assignments for gRNA DPM2_-_130700076.23-P1P2 (360 of 1926)
Performing gRNA-to-cell assignments for gRNA DPM3_+_155113020.23-P1P2 (365 of 1926)
Performing gRNA-to-cell assignments for gRNA DROSHA_-_31531982.23-P1P2 (370 of 1926)
Performing gRNA-to-cell assignments for gRNA EEF2_+_3985081.23-P1 (375 of 1926)
Performing gRNA-to-cell assignments for gRNA EEF2_-_3985443.23-P1 (380 of 1926)
Performing gRNA-to-cell assignments for gRNA EHMT2_-_31865199.23-P1P2 (385 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF3F_+_8009347.23-P1P2 (390 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF3M_-_32605405.23-P1P2 (395 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF4E_+_99849831.23-P1P2 (400 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF4G1_-_184032364.23-P1P2 (405 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF4G1_-_184033042.23-P1P2 (410 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF4G2_+_10830447.23-P1 (415 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF5A_-_7211217.23-P1P2 (420 of 1926)
Performing gRNA-to-cell assignments for gRNA TCEB3_-_24069975.23-P1P2 (425 of 1926)
Performing gRNA-to-cell assignments for gRNA EMC2_+_109456080.23-P1P2 (430 of 1926)
Performing gRNA-to-cell assignments for gRNA EME2_-_1823220.23-P1P2 (435 of 1926)
Performing gRNA-to-cell assignments for gRNA EME2_-_1823252.23-P1P2 (440 of 1926)
Performing gRNA-to-cell assignments for gRNA EOMES_-_27763342.23-P1P2 (445 of 1926)
Performing gRNA-to-cell assignments for gRNA ERCC2_+_45873802.23-P1P2 (450 of 1926)
Performing gRNA-to-cell assignments for gRNA ERO1L_+_53162425.23-P1P2 (455 of 1926)
Performing gRNA-to-cell assignments for gRNA ETS1_+_128392044.23-P1P2 (460 of 1926)
Performing gRNA-to-cell assignments for gRNA ETV2_+_36134592.23-P1P2 (465 of 1926)
Performing gRNA-to-cell assignments for gRNA ETV2_-_36134789.23-P1P2 (470 of 1926)
Performing gRNA-to-cell assignments for gRNA EWSR1_-_29664347.23-P1P2 (475 of 1926)
Performing gRNA-to-cell assignments for gRNA EXOSC8_-_37574979.23-P1P2 (480 of 1926)
Performing gRNA-to-cell assignments for gRNA FBXW11_-_171433300.23-P1P2 (485 of 1926)
Performing gRNA-to-cell assignments for gRNA FEN1_+_61560782.23-P1P2 (490 of 1926)
Performing gRNA-to-cell assignments for gRNA FUT11_+_75532080.23-P1P2 (495 of 1926)
Performing gRNA-to-cell assignments for gRNA FUT11_-_75532461.23-P1P2 (500 of 1926)
Performing gRNA-to-cell assignments for gRNA FUT7_-_139926495.23-P1P2 (505 of 1926)
Performing gRNA-to-cell assignments for gRNA G3BP1_-_151151499.23-P1P2 (510 of 1926)
Performing gRNA-to-cell assignments for gRNA GCK_+_44229021.23-P1P2 (515 of 1926)
Performing gRNA-to-cell assignments for gRNA GOLT1B_+_21654777.23-P1P2 (520 of 1926)
Performing gRNA-to-cell assignments for gRNA GOSR2_+_45000529.23-P1P2 (525 of 1926)
Performing gRNA-to-cell assignments for gRNA GOSR2_-_45000948.23-P1P2 (530 of 1926)
Performing gRNA-to-cell assignments for gRNA GPI_+_34856284.23-P1P2 (535 of 1926)
Performing gRNA-to-cell assignments for gRNA GPR107_-_132816044.23-P1P2 (540 of 1926)
Performing gRNA-to-cell assignments for gRNA GSPT1_-_12009714.23-P1P2 (545 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF2A1_+_81687215.23-P1P2 (550 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF2A2_+_59949355.23-P1P2 (555 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF2A2_-_59949657.23-P1P2 (560 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF2H4_-_30876185.23-P1P2 (565 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF3C2_+_27579604.23-P1P2 (570 of 1926)
Performing gRNA-to-cell assignments for gRNA HAND1_+_153857732.23-P1P2 (575 of 1926)
Performing gRNA-to-cell assignments for gRNA HDAC2_-_114292253.23-P1P2 (580 of 1926)
Performing gRNA-to-cell assignments for gRNA HGF_+_81399312.23-P1P2 (585 of 1926)
Performing gRNA-to-cell assignments for gRNA HGF_-_81399351.23-P1P2 (590 of 1926)
Performing gRNA-to-cell assignments for gRNA HINT2_+_35814995.23-P1P2 (595 of 1926)
Performing gRNA-to-cell assignments for gRNA HNF4A_+_43029971.23-P1P2 (600 of 1926)
Performing gRNA-to-cell assignments for gRNA HNRNPM_-_8509885.23-P1 (605 of 1926)
Performing gRNA-to-cell assignments for gRNA HNRNPU_+_245027753.23-P1P2 (610 of 1926)
Performing gRNA-to-cell assignments for gRNA HSD17B12_+_43702292.23-P1P2 (615 of 1926)
Performing gRNA-to-cell assignments for gRNA HSD17B12_-_43702305.23-P1P2 (620 of 1926)
Performing gRNA-to-cell assignments for gRNA HSPA5_-_128003534.23-P1P2 (625 of 1926)
Performing gRNA-to-cell assignments for gRNA HUWE1_-_53710540.23-P1P2 (630 of 1926)
Performing gRNA-to-cell assignments for gRNA INO80E_-_30007588.23-P1P2 (635 of 1926)
Performing gRNA-to-cell assignments for gRNA INTS10_+_19675000.23-P1P2 (640 of 1926)
Performing gRNA-to-cell assignments for gRNA VWA9_+_65903046.23-P1P2 (645 of 1926)
Performing gRNA-to-cell assignments for gRNA VWA9_-_65903523.23-P1P2 (650 of 1926)
Performing gRNA-to-cell assignments for gRNA C7orf26_+_6630358.23-P1P2 (655 of 1926)
Performing gRNA-to-cell assignments for gRNA IPO8_+_30848771.23-P1P2 (660 of 1926)
Performing gRNA-to-cell assignments for gRNA ITGB1BP1_+_9563543.23-P1P2 (665 of 1926)
Performing gRNA-to-cell assignments for gRNA SSFA2_+_182756743.23-P1P2 (670 of 1926)
Performing gRNA-to-cell assignments for gRNA KCTD10_+_109914946.23-P1P2 (675 of 1926)
Performing gRNA-to-cell assignments for gRNA KCTD10_-_109915112.23-P1P2 (680 of 1926)
Performing gRNA-to-cell assignments for gRNA KDM1A_-_23346409.23-P1P2 (685 of 1926)
Performing gRNA-to-cell assignments for gRNA KDM3B_-_137688431.23-P1P2 (690 of 1926)
Performing gRNA-to-cell assignments for gRNA KDR_-_55991659.23-P1P2 (695 of 1926)
Performing gRNA-to-cell assignments for gRNA KIAA0232_+_6784651.23-P1P2 (700 of 1926)
Performing gRNA-to-cell assignments for gRNA LAGE3_+_153706791.23-P1P2 (705 of 1926)
Performing gRNA-to-cell assignments for gRNA LAGE3_+_153707248.23-P1P2 (710 of 1926)
Performing gRNA-to-cell assignments for gRNA LAMTOR4_-_99746921.23-P1P2 (715 of 1926)
Performing gRNA-to-cell assignments for gRNA LARP4_-_50794864.23-P1P2 (720 of 1926)
Performing gRNA-to-cell assignments for gRNA LATS2_-_21635336.23-P1P2 (725 of 1926)
Performing gRNA-to-cell assignments for gRNA LEMD3_+_65563441.23-P1P2 (730 of 1926)
Performing gRNA-to-cell assignments for gRNA LENG8_+_54960428.23-P1P2 (735 of 1926)
Performing gRNA-to-cell assignments for gRNA LENG8_-_54960616.23-P1P2 (740 of 1926)
Performing gRNA-to-cell assignments for gRNA LHFPL2_-_77944519.23-P1P2 (745 of 1926)
Performing gRNA-to-cell assignments for gRNA LSM12_-_42143603.23-P1P2 (750 of 1926)
Performing gRNA-to-cell assignments for gRNA MANBAL_+_35918136.23-P1P2 (755 of 1926)
Performing gRNA-to-cell assignments for gRNA MAU2_+_19431732.23-P1P2 (760 of 1926)
Performing gRNA-to-cell assignments for gRNA MBD1_+_47807678.23-P1P2 (765 of 1926)
Performing gRNA-to-cell assignments for gRNA MBD1_-_47807877.23-P1P2 (770 of 1926)
Performing gRNA-to-cell assignments for gRNA MEAF6_-_37980319.23-P1P2 (775 of 1926)
Performing gRNA-to-cell assignments for gRNA MED1_-_37607140.23-P1P2 (780 of 1926)
Performing gRNA-to-cell assignments for gRNA MED12_+_70338541.23-P1P2 (785 of 1926)
Performing gRNA-to-cell assignments for gRNA MED15_+_20862189.23-P1 (790 of 1926)
Performing gRNA-to-cell assignments for gRNA MED19_+_57479335.23-P1P2 (795 of 1926)
Performing gRNA-to-cell assignments for gRNA MED19_-_57479437.23-P1P2 (800 of 1926)
Performing gRNA-to-cell assignments for gRNA MEMO1_-_32264882.23-P1 (805 of 1926)
Performing gRNA-to-cell assignments for gRNA MIXL1_+_226411425.23-P1P2 (810 of 1926)
Performing gRNA-to-cell assignments for gRNA MRPL28_+_420487.23-P1P2 (815 of 1926)
Performing gRNA-to-cell assignments for gRNA MRPL32_+_42972018.23-P1P2 (820 of 1926)
Performing gRNA-to-cell assignments for gRNA MRPL47_+_179322181.23-P1P2 (825 of 1926)
Performing gRNA-to-cell assignments for gRNA MRPL47_-_179322213.23-P1P2 (830 of 1926)
Performing gRNA-to-cell assignments for gRNA LUZP6_-_135612750.23-P1P2 (835 of 1926)
Performing gRNA-to-cell assignments for gRNA NAA20_-_19998023.23-P1P2 (840 of 1926)
Performing gRNA-to-cell assignments for gRNA NAA25_+_112546538.23-P1P2 (845 of 1926)
Performing gRNA-to-cell assignments for gRNA NCBP1_+_100396082.23-P1P2 (850 of 1926)
Performing gRNA-to-cell assignments for gRNA NCEH1_+_172428630.23-P1P2 (855 of 1926)
Performing gRNA-to-cell assignments for gRNA NCEH1_-_172428841.23-P1P2 (860 of 1926)
Performing gRNA-to-cell assignments for gRNA NEDD8_-_24701363.23-P1P2 (865 of 1926)
Performing gRNA-to-cell assignments for gRNA NF2_+_29999635.23-P1P2 (870 of 1926)
Performing gRNA-to-cell assignments for gRNA NFIB_+_14314559.23-P1P2 (875 of 1926)
Performing gRNA-to-cell assignments for gRNA NOTCH1_+_139418329.23-P2 (880 of 1926)
Performing gRNA-to-cell assignments for gRNA NOTCH1_-_139418349.23-P2 (885 of 1926)
Performing gRNA-to-cell assignments for gRNA NOTCH1_-_139440503.23-P1 (890 of 1926)
Performing gRNA-to-cell assignments for gRNA NUFIP2_+_27621073.23-P1P2 (895 of 1926)
Performing gRNA-to-cell assignments for gRNA NUP93_-_56764090.23-P1P2 (900 of 1926)
Performing gRNA-to-cell assignments for gRNA OGT_+_70752948.23-P1P2 (905 of 1926)
Performing gRNA-to-cell assignments for gRNA OSBP_+_59383170.23-P1P2 (910 of 1926)
Performing gRNA-to-cell assignments for gRNA OSBP_-_59383348.23-P1P2 (915 of 1926)
Performing gRNA-to-cell assignments for gRNA PAM16_-_4401236.23-P1P2 (920 of 1926)
Performing gRNA-to-cell assignments for gRNA PAXBP1_-_34143614.23-P1P2 (925 of 1926)
Performing gRNA-to-cell assignments for gRNA PBLD_-_70092436.23-P1 (930 of 1926)
Performing gRNA-to-cell assignments for gRNA PCBP2_+_53835614.23-P2 (935 of 1926)
Performing gRNA-to-cell assignments for gRNA PEG10_+_94292752.23-P2 (940 of 1926)
Performing gRNA-to-cell assignments for gRNA PEG10_-_94293027.23-P2 (945 of 1926)
Performing gRNA-to-cell assignments for gRNA PGRMC2_+_129208913.23-P1P2 (950 of 1926)
Performing gRNA-to-cell assignments for gRNA PIK3CA_-_178866211.23-P1P2 (955 of 1926)
Performing gRNA-to-cell assignments for gRNA PITPNC1_-_65373433.23-P1P2 (960 of 1926)
Performing gRNA-to-cell assignments for gRNA PLCG1_+_39766176.23-P1P2 (965 of 1926)
Performing gRNA-to-cell assignments for gRNA PLEKHG1_+_151042273.23-P2 (970 of 1926)
Performing gRNA-to-cell assignments for gRNA PLEKHG1_-_151042422.23-P2 (975 of 1926)
Performing gRNA-to-cell assignments for gRNA PLRG1_-_155471154.23-P1P2 (980 of 1926)
Performing gRNA-to-cell assignments for gRNA PLXNA2_-_208417275.23-P1P2 (985 of 1926)
Performing gRNA-to-cell assignments for gRNA POLE4_-_75185725.23-P1P2 (990 of 1926)
Performing gRNA-to-cell assignments for gRNA POLR2M_+_57998962.23-P1P2 (995 of 1926)
Performing gRNA-to-cell assignments for gRNA POLR3H_+_41940045.23-P1P2 (1000 of 1926)
Performing gRNA-to-cell assignments for gRNA POLR3H_+_41940414.23-P1P2 (1005 of 1926)
Performing gRNA-to-cell assignments for gRNA POLR3K_-_103144.23-P1P2 (1010 of 1926)
Performing gRNA-to-cell assignments for gRNA PPHLN1_-_42720007.23-P1P2 (1015 of 1926)
Performing gRNA-to-cell assignments for gRNA PPIL4_+_149867146.23-P1P2 (1020 of 1926)
Performing gRNA-to-cell assignments for gRNA PPP2R1A_+_52693317.23-P1 (1025 of 1926)
Performing gRNA-to-cell assignments for gRNA PRKCD_+_53195217.23-P1P2 (1030 of 1926)
Performing gRNA-to-cell assignments for gRNA PRKCD_-_53195332.23-P1P2 (1035 of 1926)
Performing gRNA-to-cell assignments for gRNA PRKD1_-_30396996.23-P1P2 (1040 of 1926)
Performing gRNA-to-cell assignments for gRNA PRMT5_+_23398554.23-P1 (1045 of 1926)
Performing gRNA-to-cell assignments for gRNA C2orf53_-_27362230.23-P1P2 (1050 of 1926)
Performing gRNA-to-cell assignments for gRNA PRSS8_+_31147010.23-P1P2 (1055 of 1926)
Performing gRNA-to-cell assignments for gRNA PSMB3_+_36909062.23-P1P2 (1060 of 1926)
Performing gRNA-to-cell assignments for gRNA PSMB3_-_36909230.23-P1P2 (1065 of 1926)
Performing gRNA-to-cell assignments for gRNA PSMD4_-_151227289.23-P1P2 (1070 of 1926)
Performing gRNA-to-cell assignments for gRNA PSMG4_-_3259227.23-P1P2 (1075 of 1926)
Performing gRNA-to-cell assignments for gRNA PTAR1_+_72374835.23-P1P2 (1080 of 1926)
Performing gRNA-to-cell assignments for gRNA PTBP1_+_797463.23-P1P2 (1085 of 1926)
Performing gRNA-to-cell assignments for gRNA PTEN_+_89623388.23-P1P2 (1090 of 1926)
Performing gRNA-to-cell assignments for gRNA PTEN_-_89623526.23-P1P2 (1095 of 1926)
Performing gRNA-to-cell assignments for gRNA PTMS_-_6875591.23-P1P2 (1100 of 1926)
Performing gRNA-to-cell assignments for gRNA PYGO2_-_154934135.23-P1P2 (1105 of 1926)
Performing gRNA-to-cell assignments for gRNA RAB18_-_27793292.23-P1P2 (1110 of 1926)
Performing gRNA-to-cell assignments for gRNA RAB6A_+_73471983.23-P1P2 (1115 of 1926)
Performing gRNA-to-cell assignments for gRNA RABGGTB_+_76251918.23-P1P2 (1120 of 1926)
Performing gRNA-to-cell assignments for gRNA RABGGTB_-_76252008.23-P1P2 (1125 of 1926)
Performing gRNA-to-cell assignments for gRNA GNB2L1_-_180670840.23-P1P2 (1130 of 1926)
Performing gRNA-to-cell assignments for gRNA RAD21_-_117886635.23-P1P2 (1135 of 1926)
Performing gRNA-to-cell assignments for gRNA RANGRF_+_8192257.23-P1P2 (1140 of 1926)
Performing gRNA-to-cell assignments for gRNA RBBP5_+_205091050.23-P1P2 (1145 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM14-RBM4_+_66406280.23-P2 (1150 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM33_+_155437802.23-P1P2 (1155 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM4_+_66406319.23-P2 (1160 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM7_+_114271352.23-P1P2 (1165 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM7_-_114271700.23-P1P2 (1170 of 1926)
Performing gRNA-to-cell assignments for gRNA RBX1_-_41347852.23-P1P2 (1175 of 1926)
Performing gRNA-to-cell assignments for gRNA KIAA1551_+_32112399.23-P1P2 (1180 of 1926)
Performing gRNA-to-cell assignments for gRNA RFX3_-_3525750.23-P1P2 (1185 of 1926)
Performing gRNA-to-cell assignments for gRNA RGP1_+_35749328.23-P1P2 (1190 of 1926)
Performing gRNA-to-cell assignments for gRNA RHOBTB3_+_95067102.23-P1P2 (1195 of 1926)
Performing gRNA-to-cell assignments for gRNA RHOBTB3_-_95067164.23-P1P2 (1200 of 1926)
Performing gRNA-to-cell assignments for gRNA RIBC1_-_53449606.23-P1P2 (1205 of 1926)
Performing gRNA-to-cell assignments for gRNA KIAA1432_+_5629173.23-P1P2 (1210 of 1926)
Performing gRNA-to-cell assignments for gRNA RLF_+_40627522.23-P1P2 (1215 of 1926)
Performing gRNA-to-cell assignments for gRNA RNF115_-_145611148.23-P1P2 (1220 of 1926)
Performing gRNA-to-cell assignments for gRNA RNF20_+_104296142.23-P1P2 (1225 of 1926)
Performing gRNA-to-cell assignments for gRNA RNF20_-_104296250.23-P1P2 (1230 of 1926)
Performing gRNA-to-cell assignments for gRNA RNF40_-_30773681.23-P1P2 (1235 of 1926)
Performing gRNA-to-cell assignments for gRNA RPL37A_-_217363600.23-P1P2 (1240 of 1926)
Performing gRNA-to-cell assignments for gRNA RPS15A_-_18801590.23-P1P2 (1245 of 1926)
Performing gRNA-to-cell assignments for gRNA SAV1_+_51134992.23-P1P2 (1250 of 1926)
Performing gRNA-to-cell assignments for gRNA SCAF4_+_33104038.23-P1P2 (1255 of 1926)
Performing gRNA-to-cell assignments for gRNA SCAF4_-_33104343.23-P1P2 (1260 of 1926)
Performing gRNA-to-cell assignments for gRNA SCAF8_+_155054976.23-P1P2 (1265 of 1926)
Performing gRNA-to-cell assignments for gRNA EPT1_+_26569050.23-P1P2 (1270 of 1926)
Performing gRNA-to-cell assignments for gRNA SERF2_+_44069659.23-P2 (1275 of 1926)
Performing gRNA-to-cell assignments for gRNA SERPING1_+_57365137.23-P1P2 (1280 of 1926)
Performing gRNA-to-cell assignments for gRNA SETD1B_+_122242041.23-P1P2 (1285 of 1926)
Performing gRNA-to-cell assignments for gRNA SETD1B_-_122242294.23-P1P2 (1290 of 1926)
Performing gRNA-to-cell assignments for gRNA SETDB1_-_150898815.23-P1P2 (1295 of 1926)
Performing gRNA-to-cell assignments for gRNA SF1_-_64545937.23-P1P2 (1300 of 1926)
Performing gRNA-to-cell assignments for gRNA CCDC101_-_28565268.23-P1P2 (1305 of 1926)
Performing gRNA-to-cell assignments for gRNA SKA3_+_21750519.23-P1P2 (1310 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR61_+_78591497.23-P1P2 (1315 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR61_-_78591737.23-P1P2 (1320 of 1926)
Performing gRNA-to-cell assignments for gRNA SLC29A1_-_44187438.23-P2 (1325 of 1926)
Performing gRNA-to-cell assignments for gRNA SLC35C1_-_45826711.23-P1P2 (1330 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD1_+_146403428.23-P1P2 (1335 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD3_+_67358229.23-P1P2 (1340 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD4_+_48556623.23-P1P2 (1345 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD4_-_48556657.23-P1P2 (1350 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD5_-_135468694.23-P1P2 (1355 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCA4_-_11071741.23-P1P2 (1360 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCA5_+_144435243.23-P1P2 (1365 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCB1_+_24129228.23-P1 (1370 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCC1_+_47822975.23-P1P2 (1375 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCC1_-_47822893.23-P1P2 (1380 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCD3_-_150974212.23-P1 (1385 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCE1_-_38803669.23-P1P2 (1390 of 1926)
Performing gRNA-to-cell assignments for gRNA SMC2_-_106856631.23-P1P2 (1395 of 1926)
Performing gRNA-to-cell assignments for gRNA SMG5_+_156252537.23-P1P2 (1400 of 1926)
Performing gRNA-to-cell assignments for gRNA SNRNP25_+_103852.23-P1P2 (1405 of 1926)
Performing gRNA-to-cell assignments for gRNA SNRNP25_-_104242.23-P1P2 (1410 of 1926)
Performing gRNA-to-cell assignments for gRNA SNRNP27_-_70121181.23-P1P2 (1415 of 1926)
Performing gRNA-to-cell assignments for gRNA SNUPN_-_75917960.23-P1P2 (1420 of 1926)
Performing gRNA-to-cell assignments for gRNA SOX4_+_21594213.23-P1 (1425 of 1926)
Performing gRNA-to-cell assignments for gRNA SPIN1_+_91003737.23-P1P2 (1430 of 1926)
Performing gRNA-to-cell assignments for gRNA SPPL2B_+_2328680.23-P1P2 (1435 of 1926)
Performing gRNA-to-cell assignments for gRNA SPPL2B_-_2328912.23-P1P2 (1440 of 1926)
Performing gRNA-to-cell assignments for gRNA SPTY2D1_-_18655846.23-P1P2 (1445 of 1926)
Performing gRNA-to-cell assignments for gRNA SRP19_+_112197188.23-P1P2 (1450 of 1926)
Performing gRNA-to-cell assignments for gRNA SRP54_-_35452244.23-P1P2 (1455 of 1926)
Performing gRNA-to-cell assignments for gRNA SRP68_+_74068540.23-P1P2 (1460 of 1926)
Performing gRNA-to-cell assignments for gRNA SRPR_+_126138730.23-P1P2 (1465 of 1926)
Performing gRNA-to-cell assignments for gRNA SRPR_-_126138742.23-P1P2 (1470 of 1926)
Performing gRNA-to-cell assignments for gRNA SRPRB_-_133524719.23-P1P2 (1475 of 1926)
Performing gRNA-to-cell assignments for gRNA SRRM2_-_2802714.23-P1P2 (1480 of 1926)
Performing gRNA-to-cell assignments for gRNA SRSF6_+_42086993.23-P1 (1485 of 1926)
Performing gRNA-to-cell assignments for gRNA SSB_+_170655450.23-P1P2 (1490 of 1926)
Performing gRNA-to-cell assignments for gRNA SSRP1_+_57102902.23-P1P2 (1495 of 1926)
Performing gRNA-to-cell assignments for gRNA SSRP1_-_57103355.23-P1P2 (1500 of 1926)
Performing gRNA-to-cell assignments for gRNA STAG2_-_123095658.23-P1 (1505 of 1926)
Performing gRNA-to-cell assignments for gRNA STK40_+_36851507.23-P1P2 (1510 of 1926)
Performing gRNA-to-cell assignments for gRNA SUPT20H_+_37633798.23-P1P2 (1515 of 1926)
Performing gRNA-to-cell assignments for gRNA SUZ12_+_30264157.23-P1P2 (1520 of 1926)
Performing gRNA-to-cell assignments for gRNA SYVN1_+_64901973.23-P1P2 (1525 of 1926)
Performing gRNA-to-cell assignments for gRNA SYVN1_-_64901981.23-P1P2 (1530 of 1926)
Performing gRNA-to-cell assignments for gRNA TADA2B_-_7045121.23-P1P2 (1535 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF10_+_6633406.23-P1P2 (1540 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF2_+_120845032.23-P1P2 (1545 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF3_-_7860540.23-P1P2 (1550 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF4_+_60640814.23-P1P2 (1555 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF4_-_60640498.23-P1P2 (1560 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF6L_-_62538954.23-P1P2 (1565 of 1926)
Performing gRNA-to-cell assignments for gRNA TAOK1_-_27717907.23-P1P2 (1570 of 1926)
Performing gRNA-to-cell assignments for gRNA FAM208A_+_56717273.23-P1P2 (1575 of 1926)
Performing gRNA-to-cell assignments for gRNA TBC1D9_+_141677481.23-P1P2 (1580 of 1926)
Performing gRNA-to-cell assignments for gRNA TBP_+_170863524.23-P1P2 (1585 of 1926)
Performing gRNA-to-cell assignments for gRNA TBP_-_170863481.23-P1P2 (1590 of 1926)
Performing gRNA-to-cell assignments for gRNA TBX3_-_115121998.23-P1P2 (1595 of 1926)
Performing gRNA-to-cell assignments for gRNA T_-_166581333.23-P1P2 (1600 of 1926)
Performing gRNA-to-cell assignments for gRNA TGFBR2_+_30648526.23-P1P2 (1605 of 1926)
Performing gRNA-to-cell assignments for gRNA THBD_+_23029196.23-P2 (1610 of 1926)
Performing gRNA-to-cell assignments for gRNA TLE6_+_2977592.23-P1P2 (1615 of 1926)
Performing gRNA-to-cell assignments for gRNA TLE6_-_2977643.23-P1P2 (1620 of 1926)
Performing gRNA-to-cell assignments for gRNA MESDC1_-_81293340.23-P1P2 (1625 of 1926)
Performing gRNA-to-cell assignments for gRNA TMEM41B_+_9336156.23-P1P2 (1630 of 1926)
Performing gRNA-to-cell assignments for gRNA TMTC3_-_88536127.23-P1P2 (1635 of 1926)
Performing gRNA-to-cell assignments for gRNA TOP1_+_39657538.23-P1P2 (1640 of 1926)
Performing gRNA-to-cell assignments for gRNA TP53_+_7590531.23-P1P2 (1645 of 1926)
Performing gRNA-to-cell assignments for gRNA TP53_-_7590704.23-P1P2 (1650 of 1926)
Performing gRNA-to-cell assignments for gRNA TRIM24_-_138145308.23-P1P2 (1655 of 1926)
Performing gRNA-to-cell assignments for gRNA TRRAP_-_98476465.23-P1P2 (1660 of 1926)
Performing gRNA-to-cell assignments for gRNA TSTD2_+_100395594.23-P1P2 (1665 of 1926)
Performing gRNA-to-cell assignments for gRNA UBE2E1_+_23847541.23-P1P2 (1670 of 1926)
Performing gRNA-to-cell assignments for gRNA UFM1_+_38924011.23-P1P2 (1675 of 1926)
Performing gRNA-to-cell assignments for gRNA UFM1_-_38924421.23-P1P2 (1680 of 1926)
Performing gRNA-to-cell assignments for gRNA UPF2_+_12084799.23-P1P2 (1685 of 1926)
Performing gRNA-to-cell assignments for gRNA UXT_-_47518505.23-P1P2 (1690 of 1926)
Performing gRNA-to-cell assignments for gRNA VHL_-_10183525.23-P1P2 (1695 of 1926)
Performing gRNA-to-cell assignments for gRNA VPS54_+_64246104.23-P1P2 (1700 of 1926)
Performing gRNA-to-cell assignments for gRNA WAC_+_28821698.23-P1P2 (1705 of 1926)
Performing gRNA-to-cell assignments for gRNA WAC_-_28822074.23-P1P2 (1710 of 1926)
Performing gRNA-to-cell assignments for gRNA WAPAL_+_88281533.23-P1P2 (1715 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR12_-_203776383.23-P1P2 (1720 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR77_+_111991860.23-P1P2 (1725 of 1926)
Performing gRNA-to-cell assignments for gRNA WIZ_-_15543016.23-P1P2 (1730 of 1926)
Performing gRNA-to-cell assignments for gRNA WWC3_+_9983328.23-P1P2 (1735 of 1926)
Performing gRNA-to-cell assignments for gRNA WWC3_-_9983510.23-P1P2 (1740 of 1926)
Performing gRNA-to-cell assignments for gRNA XPO5_-_43543670.23-P1P2 (1745 of 1926)
Performing gRNA-to-cell assignments for gRNA XRN1_-_142166648.23-P1P2 (1750 of 1926)
Performing gRNA-to-cell assignments for gRNA YARS2_+_32908787.23-P1P2 (1755 of 1926)
Performing gRNA-to-cell assignments for gRNA ZBTB22_+_33285454.23-P1P2 (1760 of 1926)
Performing gRNA-to-cell assignments for gRNA ZBTB48_+_6640093.23-P1P2 (1765 of 1926)
Performing gRNA-to-cell assignments for gRNA ZBTB48_-_6640138.23-P1P2 (1770 of 1926)
Performing gRNA-to-cell assignments for gRNA ZC3H8_-_113012435.23-P1P2 (1775 of 1926)
Performing gRNA-to-cell assignments for gRNA ZNF281_+_200379151.23-P1P2 (1780 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_00117 (1785 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_00233 (1790 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_00518 (1795 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_00656 (1800 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_01226 (1805 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_01730 (1810 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_02009 (1815 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_02288 (1820 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_02777 (1825 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_03026 (1830 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_03252 (1835 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_03620 (1840 of 1926)
Performing gRNA-to-cell assignments for gRNA CSNK2B_SGTN_04 (1845 of 1926)
Performing gRNA-to-cell assignments for gRNA INO80E_SGTN_25 (1850 of 1926)
Performing gRNA-to-cell assignments for gRNA PIK3CA_SGTN_26 (1855 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_1044 (1860 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_1414 (1865 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_1554 (1870 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_1785 (1875 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_2009 (1880 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_2271 (1885 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_2425 (1890 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_2614 (1895 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_293 (1900 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_541 (1905 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_725 (1910 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_860 (1915 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR61_SGTN_14 (1920 of 1926)
Performing gRNA-to-cell assignments for gRNA TP53_SGTN_01 (1925 of 1926)
Skipping calibration check.
Running discovery analysis...
Warning: The calibration check (`run_calibration_check()`) should be run before the discovery analysis.

Running discovery_analysis in parallel. Change directories to /scratch/users/edtng/sceptre_hippo_trans_gene_effects/d2/sceptre_logs/sceptre_logs/ and view the files discovery_analysis_*.out for progress updates. 
```

## Code cell 5

```python
# %% [markdown]
# ## Combine and annotate SCEPTRE outputs
# 
# This cell merges SCEPTRE results with response gene symbols and saves one combined table across days.

# %%
combined = []
for day in DAYS_TO_RUN:
    day_dir = SCEPTRE_OUT_ROOT / day
    result_csv = day_dir / "sceptre_results" / f"sceptre_hippo_vs_genes_{day}.csv"
    if not result_csv.exists():
        print(f"Missing result for {day}: {result_csv}")
        continue

    res = pd.read_csv(result_csv)
    res["day"] = day

    meta = pd.read_csv(day_dir / "gene_symbol_map.csv")
    res = res.merge(meta, on="response_id", how="left")
    res["is_self_target"] = res["grna_target"].astype(str) == res["response_symbol"].astype(str)
    combined.append(res)

if not combined:
    raise ValueError("No SCEPTRE result files found to combine.")

combined_df = pd.concat(combined, ignore_index=True)

# Sort by a likely p-value column if present.
pval_cols = [c for c in ["p_value", "p.value", "pval", "p_val"] if c in combined_df.columns]
if pval_cols:
    combined_df = combined_df.sort_values(pval_cols[0])

combined_path = SCEPTRE_OUT_ROOT / "sceptre_hippo_vs_genes_all_days.annotated.csv.gz"
combined_df.to_csv(combined_path, index=False, compression="gzip")

print("Combined shape:", combined_df.shape)
print("Saved:", combined_path)
print("Columns:", list(combined_df.columns))
display(combined_df.head())

# %%
# Quick summaries.
summary_cols = ["day", "grna_target"]
if "p_value" in combined_df.columns:
    sig = combined_df[combined_df["p_value"] < 0.05].copy()
    print("Nominal p < 0.05 counts:")
    display(sig.groupby(summary_cols).size().rename("n_nominal_p_lt_0_05").reset_index())

# If SCEPTRE returns an adjusted p/q-value column, summarize that too.
q_cols = [c for c in combined_df.columns if c.lower() in {"q_value", "qval", "fdr", "padj", "adjusted_p_value", "adjusted_p.value"}]
print("Candidate adjusted p/q columns:", q_cols)
for q in q_cols[:3]:
    sig_q = combined_df[combined_df[q] < 0.1]
    print(f"\n{q} < 0.1 counts:")
    display(sig_q.groupby(summary_cols).size().rename(f"n_{q}_lt_0_1").reset_index())

# Strongest hits per target/day by p-value.
if "p_value" in combined_df.columns:
    top_hits = (
        combined_df
        .sort_values("p_value")
        .groupby(["day", "grna_target"], observed=True)
        .head(25)
    )
    top_hits.to_csv(SCEPTRE_OUT_ROOT / "sceptre_hippo_vs_genes_top25_per_day_target.csv", index=False)
    display(top_hits[["day", "grna_target", "response_id", "response_symbol", "p_value"] + [c for c in ["log_2_fold_change", "z_value", "statistic"] if c in top_hits.columns]].head(50))

# %% [markdown]
# ## Notes
# 
# - This notebook tests **trans effects**, i.e. each Hippo target against many response genes.
# - The discovery-pair table has required SCEPTRE columns `grna_target` and `response_id`.
# - `control_group="complement"` matches your original notebook.
# - For first runs, use `DAYS_TO_RUN=["d2"]` and `MAX_DISCOVERY_PAIRS_FOR_TEST=1000` to check that the pipeline works.
# - For full runs, set `MAX_DISCOVERY_PAIRS_FOR_TEST=0`.
# - If runtime is too long, increase `MIN_DETECTED_CELLS`, reduce `HIPPO_TARGETS`, or split days/targets.
```

