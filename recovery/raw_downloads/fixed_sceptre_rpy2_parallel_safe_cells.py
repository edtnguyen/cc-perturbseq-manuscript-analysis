# %% [markdown]
# # Parallel-safe SCEPTRE rpy2 run cells
#
# These cells limit BLAS/OpenMP threads to avoid oversubscribing Sherlock nodes.
# Restart the kernel before running Cell 1.

# %%
# Cell 1 — run this before importing numpy/scipy/rpy2/scanpy, then restart the kernel if those were already imported.
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

# %%
# Cell 2 — load rpy2 after thread env vars are set.
%load_ext rpy2.ipython

# %%
# Cell 3 — safe SCEPTRE run settings.
# Use the SLURM allocation but leave some overhead for code-server/Jupyter.
import os

requested_cpus = int(os.environ.get("SLURM_CPUS_PER_TASK") or os.environ.get("SLURM_CPUS_ON_NODE") or 1)

# In code-server/Jupyter, do not use all allocated CPUs; leave 1-2 cores for the server/kernel overhead.
if requested_cpus >= 4:
    safe_workers = max(1, requested_cpus - 2)
else:
    safe_workers = 1

SCEPTRE_OUT_ROOT_STR = str(SCEPTRE_OUT_ROOT)
DAYS_TO_RUN_R = ["d2"]  # start with one day
R_PARALLEL = safe_workers > 1
R_N_PROCESSORS = safe_workers
ASSIGN_GRNA_METHOD = "mixture"
P_MITO_THRESHOLD = 0.25
RUN_CALIBRATION_CHECK = False  # turn on only after smoke test succeeds
MAX_DISCOVERY_PAIRS_FOR_TEST = 1000  # smoke test first; 0 for full run

print("SLURM_JOB_ID:", os.environ.get("SLURM_JOB_ID"))
print("SLURM_CPUS_PER_TASK:", os.environ.get("SLURM_CPUS_PER_TASK"))
print("SLURM_CPUS_ON_NODE:", os.environ.get("SLURM_CPUS_ON_NODE"))
print("requested_cpus inferred:", requested_cpus)
print("R_PARALLEL:", R_PARALLEL)
print("R_N_PROCESSORS:", R_N_PROCESSORS)
print("DAYS_TO_RUN_R:", DAYS_TO_RUN_R)
print("MAX_DISCOVERY_PAIRS_FOR_TEST:", MAX_DISCOVERY_PAIRS_FOR_TEST)

# %%
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

# %%
%%R
# Cell 5 — run the analysis after the function is defined.
# Start with DAYS_TO_RUN_R = c("d2") and MAX_DISCOVERY_PAIRS_FOR_TEST = 1000.
result_files <- c()
for (day in days) {
  result_files <- c(result_files, run_sceptre_day(day))
}
cat("\nAll result files:\n")
print(result_files)
