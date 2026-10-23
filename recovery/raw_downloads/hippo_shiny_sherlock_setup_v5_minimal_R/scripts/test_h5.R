#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(hdf5r))
h5_path <- Sys.getenv("SHINY_UMAP_H5", unset = "/scratch/users/edtng/Hippo_Shiny_adata/data/umap_gene_program_data.h5")
stopifnot(file.exists(h5_path))
h5 <- H5File$new(h5_path, mode = "r")
on.exit(h5$close_all(), add = TRUE)
read_chr <- function(path) as.character(h5[[path]][])
cat("HDF5:", h5_path, "\n")
cat("n cells:", length(h5[["obs/cell_id"]][]), "\n")
cat("n genes:", length(h5[["genes/names"]][]), "\n")
cat("n programs:", length(h5[["programs/names"]][]), "\n")
cat("First genes:", paste(head(read_chr("genes/names"), 10), collapse = ", "), "\n")
cat("First programs:", paste(head(read_chr("programs/names"), 10), collapse = ", "), "\n")
if ("obs_annotation_names" %in% names(h5)) {
  obs_names <- read_chr("obs_annotation_names")
  obs_names <- obs_names[obs_names != ""]
  cat("obs annotations:", paste(obs_names, collapse = ", "), "\n")
  if ("day" %in% obs_names) {
    day <- read_chr("obs_annotations/day")
    cat("day counts:\n")
    print(sort(table(day, useNA = "ifany"), decreasing = FALSE))
  } else {
    cat("WARNING: obs annotation `day` is not present in this HDF5 bundle.\n")
  }
} else {
  cat("WARNING: no obs annotations are present in this HDF5 bundle.\n")
}
