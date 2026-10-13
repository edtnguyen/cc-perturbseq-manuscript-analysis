#!/usr/bin/env Rscript
repos <- Sys.getenv("CRAN_REPO", unset = "https://cloud.r-project.org")
lib <- Sys.getenv("R_LIBS_USER")
if (!nzchar(lib)) {
  lib <- file.path(Sys.getenv("HOME"), "R", paste0("x86_64-pc-linux-gnu-library/", paste(R.version$major, strsplit(R.version$minor, "\\.")[[1]][1], sep = ".")))
  Sys.setenv(R_LIBS_USER = lib)
}
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
cat("R_LIBS_USER=", lib, "\n", sep = "")
cat("PATH=", Sys.getenv("PATH"), "\n", sep = "")
cat("LD_LIBRARY_PATH=", Sys.getenv("LD_LIBRARY_PATH"), "\n", sep = "")
cat("PKG_CONFIG_PATH=", Sys.getenv("PKG_CONFIG_PATH"), "\n", sep = "")
cat("h5cc=", Sys.which("h5cc"), "\n", sep = "")
install.packages("hdf5r", repos = repos, Ncpus = 1, verbose = TRUE)
cat("hdf5r installed: ", requireNamespace("hdf5r", quietly = TRUE), "\n", sep = "")
