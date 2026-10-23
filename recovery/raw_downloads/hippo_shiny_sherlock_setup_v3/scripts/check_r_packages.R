#!/usr/bin/env Rscript
packages <- c("shiny", "bslib", "ggplot2", "hdf5r", "viridisLite", "svglite")
status <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
print(data.frame(package = packages, installed = as.logical(status)))
if (!all(status)) quit(save = "no", status = 1)
