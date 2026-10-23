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

# %% [markdown]
# ## Run SCEPTRE in R via rpy2
# 
# This cell runs the SCEPTRE high-MOI workflow for each exported day. It uses one discovery set per day containing all Hippo target × response gene pairs.
# 
# For a quick smoke test, set `MAX_DISCOVERY_PAIRS_FOR_TEST` to e.g. `1000`. For the full run, set it to `0`.

# %%
%load_ext rpy2.ipython

# %%
# R/SCEPTRE settings passed into the %%R cell below.
SCEPTRE_OUT_ROOT_STR = str(SCEPTRE_OUT_ROOT)
R_PARALLEL = True
R_N_PROCESSORS = 8
ASSIGN_GRNA_METHOD = "mixture"
P_MITO_THRESHOLD = 0.25
RUN_CALIBRATION_CHECK = False
MAX_DISCOVERY_PAIRS_FOR_TEST = 0   # set to 1000 for a smoke test, 0 for all pairs


# %%
%%R -i SCEPTRE_OUT_ROOT_STR -i DAYS_TO_RUN -i R_PARALLEL -i R_N_PROCESSORS -i ASSIGN_GRNA_METHOD -i P_MITO_THRESHOLD -i RUN_CALIBRATION_CHECK -i MAX_DISCOVERY_PAIRS_FOR_TEST
library(sceptre)
library(Matrix)

out_root <- SCEPTRE_OUT_ROOT_STR
days <- as.character(DAYS_TO_RUN)
parallel_flag <- as.logical(R_PARALLEL)
n_processors <- as.integer(R_N_PROCESSORS)
assign_method <- as.character(ASSIGN_GRNA_METHOD)
p_mito_threshold <- as.numeric(P_MITO_THRESHOLD)
run_calibration <- as.logical(RUN_CALIBRATION_CHECK)
max_pairs <- as.integer(MAX_DISCOVERY_PAIRS_FOR_TEST)

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

  cat("Loading matrices...\n")
  response_mat <- readMM(file.path(data_dir, "response_matrix.mtx"))
  grna_mat <- readMM(file.path(data_dir, "grna_matrix.mtx"))

  gene_names <- read.csv(file.path(data_dir, "gene_names.csv"), header = FALSE)$V1
  guide_names <- read.csv(file.path(data_dir, "guide_names.csv"), header = FALSE)$V1
  cell_names <- read.csv(file.path(data_dir, "cell_names.csv"), header = FALSE)$V1

  rownames(response_mat) <- gene_names
  colnames(response_mat) <- cell_names
  rownames(grna_mat) <- guide_names
  colnames(grna_mat) <- cell_names

  cat("response_mat dim:", paste(dim(response_mat), collapse = " x "), "\n")
  cat("grna_mat dim:", paste(dim(grna_mat), collapse = " x "), "\n")

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

  # Rename mitochondrial percentage covariate to avoid SCEPTRE reserved name behavior.
  if ("response_p_mito" %in% colnames(covariates)) {
    names(covariates)[names(covariates) == "response_p_mito"] <- "mt_percent"
  }
  reserved_cols <- c("response_n_nonzero", "response_n_umis", "grna_n_nonzero", "grna_n_umis")
  covariates <- covariates[, !(colnames(covariates) %in% reserved_cols), drop = FALSE]

  cat("Importing data into SCEPTRE...\n")
  sceptre_object <- import_data(
    response_matrix = response_mat,
    grna_matrix = grna_mat,
    grna_target_data_frame = guide_mapping,
    extra_covariates = covariates,
    moi = "high"
  )

  cat("Setting analysis parameters...\n")
  sceptre_object <- set_analysis_parameters(
    sceptre_object = sceptre_object,
    discovery_pairs = discovery_pairs,
    side = "both",
    control_group = "complement"
  )

  cat("Assigning gRNAs...\n")
  sceptre_object <- assign_grnas(
    sceptre_object = sceptre_object,
    method = assign_method,
    parallel = parallel_flag,
    n_processors = n_processors
  )

  cat("Running QC...\n")
  sceptre_object <- run_qc(
    sceptre_object = sceptre_object,
    p_mito_threshold = p_mito_threshold
  )

  # Save basic plots. If plotting fails, continue.
  try({
    png(file.path(result_dir, "sceptre_QC_plot.png"), width = 2400, height = 2000, res = 300)
    print(plot(sceptre_object))
    dev.off()
  }, silent = TRUE)

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
    cat("Skipping calibration check. Set RUN_CALIBRATION_CHECK=True if needed.\n")
  }

  cat("Running discovery analysis...\n")
  sceptre_object <- run_discovery_analysis(
    sceptre_object = sceptre_object,
    output_amount = 3,
    print_progress = TRUE,
    parallel = parallel_flag,
    n_processors = n_processors,
    log_dir = log_dir
  )

  discovery_result <- get_result(sceptre_object, analysis = "run_discovery_analysis")
  discovery_result$day <- day

  out_csv <- file.path(result_dir, paste0("sceptre_hippo_vs_genes_", day, ".csv"))
  write.csv(discovery_result, out_csv, row.names = FALSE)
  cat("Saved discovery result to:", out_csv, "\n")

  return(out_csv)
}

result_files <- c()
for (day in days) {
  result_files <- c(result_files, run_sceptre_day(day))
}

cat("\nAll result files:\n")
print(result_files)


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