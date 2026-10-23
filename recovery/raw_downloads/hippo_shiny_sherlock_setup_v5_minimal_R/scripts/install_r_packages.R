#!/usr/bin/env Rscript
# Install R packages needed for the Shiny UMAP viewer.
# Run on Sherlock after loading R, e.g.:
#   module load R
#   module load hdf5  # or whichever HDF5 module is available on Sherlock
#   Rscript scripts/install_r_packages.R

packages <- c(
  "shiny",
  "bslib",
  "ggplot2",
  "hdf5r",
  "viridisLite"
)

repos <- Sys.getenv("CRAN_REPO", unset = "https://cloud.r-project.org")
lib <- Sys.getenv("R_LIBS_USER")
if (!nzchar(lib)) {
  lib <- file.path(Sys.getenv("HOME"), "R", paste0("x86_64-pc-linux-gnu-library/", paste(R.version$major, strsplit(R.version$minor, "\\.")[[1]][1], sep = ".")))
  Sys.setenv(R_LIBS_USER = lib)
}
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))

missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  message("Installing into: ", lib)
  message("Packages: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = repos, Ncpus = as.integer(Sys.getenv("NCPUS", unset = "1")))
}

still_missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(still_missing) > 0) {
  stop("Still missing packages after install: ", paste(still_missing, collapse = ", "))
}
message("R package install/check complete.")
