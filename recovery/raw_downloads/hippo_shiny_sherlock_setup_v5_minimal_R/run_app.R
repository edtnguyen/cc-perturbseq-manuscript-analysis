#!/usr/bin/env Rscript
# Run the Shiny UMAP viewer on Sherlock.

root <- Sys.getenv("HIPPO_SHINY_ROOT", unset = "/scratch/users/edtng/Hippo_Shiny_adata")
app_dir <- Sys.getenv("SHINY_UMAP_APP_DIR", unset = file.path(root, "app"))
h5_path <- Sys.getenv("SHINY_UMAP_H5", unset = file.path(root, "data", "umap_gene_program_data.h5"))
port <- as.integer(Sys.getenv("PORT", unset = "3838"))
host <- Sys.getenv("HOST", unset = "0.0.0.0")

Sys.setenv(SHINY_UMAP_H5 = h5_path)
message("App dir: ", app_dir)
message("HDF5:    ", h5_path)
message("Host:    ", host)
message("Port:    ", port)

if (!dir.exists(app_dir)) stop("App directory does not exist: ", app_dir)
if (!file.exists(h5_path)) stop("HDF5 bundle does not exist: ", h5_path)

shiny::runApp(app_dir, host = host, port = port, launch.browser = FALSE)
