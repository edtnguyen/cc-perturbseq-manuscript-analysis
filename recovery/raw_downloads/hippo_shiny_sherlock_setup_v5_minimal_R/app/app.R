# Shiny UMAP viewer for gene expression, cNMF program scores, and cell metadata.
# Default data bundle:
#   /scratch/users/edtng/Hippo_Shiny_adata/data/umap_gene_program_data.h5
# Override with:
#   Sys.setenv(SHINY_UMAP_H5 = "/path/to/umap_gene_program_data.h5")

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(ggplot2)
  library(hdf5r)
  library(viridisLite)
})

HIPPO_SHINY_ROOT <- Sys.getenv(
  "HIPPO_SHINY_ROOT",
  unset = "/scratch/users/edtng/Hippo_Shiny_adata"
)
DATA_H5 <- Sys.getenv(
  "SHINY_UMAP_H5",
  unset = file.path(HIPPO_SHINY_ROOT, "data", "umap_gene_program_data.h5")
)

if (!file.exists(DATA_H5)) {
  stop(
    "Cannot find HDF5 data bundle: ", DATA_H5, "\n",
    "Set Sys.setenv(SHINY_UMAP_H5 = '/path/to/umap_gene_program_data.h5') ",
    "or create the default scratch file under ", HIPPO_SHINY_ROOT, "/data/."
  )
}

message("Using HDF5 data bundle: ", normalizePath(DATA_H5, mustWork = TRUE))

h5 <- H5File$new(DATA_H5, mode = "r")
onStop(function() {
  try(h5$close_all(), silent = TRUE)
})

read_chr <- function(path) {
  as.character(h5[[path]][])
}

cell_id <- read_chr("obs/cell_id")
umap_df <- data.frame(
  cell_id = cell_id,
  UMAP_1 = as.numeric(h5[["obs/umap_1"]][]),
  UMAP_2 = as.numeric(h5[["obs/umap_2"]][]),
  stringsAsFactors = FALSE
)

obs_annotation_names <- character(0)
if ("obs_annotation_names" %in% names(h5)) {
  obs_annotation_names <- read_chr("obs_annotation_names")
}
obs_annotation_names <- obs_annotation_names[obs_annotation_names != ""]

for (col in obs_annotation_names) {
  umap_df[[col]] <- read_chr(paste0("obs_annotations/", col))
}

program_ids <- as.integer(h5[["programs/ids"]][])
program_names <- read_chr("programs/names")
program_labels <- read_chr("programs/labels")
program_colors <- read_chr("programs/colors")
program_choices <- setNames(program_names, paste0(program_names, " — ", program_labels))

gene_names <- read_chr("genes/names")
gene_indptr <- as.numeric(h5[["genes/indptr"]][])
gene_indices_ds <- h5[["genes/indices"]]
gene_data_ds <- h5[["genes/data"]]
gene_shape <- as.numeric(h5[["genes/shape"]][])
n_genes <- gene_shape[1]
n_cells <- gene_shape[2]

if (n_cells != nrow(umap_df)) {
  stop("Gene expression cell count does not match UMAP cell count.")
}

metadata_default <- if ("day" %in% obs_annotation_names) "day" else if (length(obs_annotation_names) > 0) obs_annotation_names[1] else NULL
feature_type_choices <- c("Program score" = "program", "Gene expression" = "gene")
if (length(obs_annotation_names) > 0) {
  feature_type_choices <- c(feature_type_choices, "Cell metadata / day" = "metadata")
}

read_h5_1d_range <- function(dataset, start0, end0_exclusive) {
  n <- end0_exclusive - start0
  if (n <= 0) {
    return(numeric(0))
  }

  # hdf5r uses 1-based R indexing. The exported sparse CSR pointers are 0-based.
  dataset[seq.int(start0 + 1, end0_exclusive)]
}

get_gene_vector <- function(gene) {
  idx <- match(gene, gene_names)
  if (is.na(idx)) {
    stop("Gene not found: ", gene)
  }

  # HDF5 sparse arrays are exported from Python as 0-based CSR arrays.
  start0 <- gene_indptr[idx]
  end0 <- gene_indptr[idx + 1]

  value <- numeric(n_cells)
  if (end0 > start0) {
    cell_pos <- as.integer(read_h5_1d_range(gene_indices_ds, start0, end0)) + 1L
    vals <- as.numeric(read_h5_1d_range(gene_data_ds, start0, end0))
    value[cell_pos] <- vals
  }
  value
}

get_program_vector <- function(program_name) {
  idx <- match(program_name, program_names)
  if (is.na(idx)) {
    stop("Program not found: ", program_name)
  }
  as.numeric(h5[["programs/scores"]][idx, ])
}

clip_values <- function(x, lower_q = 0.00, upper_q = 0.99) {
  finite_x <- x[is.finite(x)]
  if (length(finite_x) == 0) {
    return(x)
  }
  lo <- as.numeric(stats::quantile(finite_x, probs = lower_q, na.rm = TRUE, names = FALSE))
  hi <- as.numeric(stats::quantile(finite_x, probs = upper_q, na.rm = TRUE, names = FALSE))
  pmin(pmax(x, lo), hi)
}

safe_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  gsub("_+", "_", x)
}

annotation_filter_ui <- function() {
  if (length(obs_annotation_names) == 0) {
    return(NULL)
  }
  tagList(
    selectInput(
      "filter_column",
      "Optional cell annotation filter",
      choices = c("None", obs_annotation_names),
      selected = "None"
    ),
    uiOutput("filter_values_ui")
  )
}

ui <- page_sidebar(
  title = "UMAP gene/program viewer",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  sidebar = sidebar(
    radioButtons(
      "feature_type",
      "Display",
      choices = feature_type_choices,
      selected = "program"
    ),
    conditionalPanel(
      condition = "input.feature_type == 'program'",
      selectInput(
        "program",
        "Program",
        choices = program_choices,
        selected = program_names[1]
      )
    ),
    conditionalPanel(
      condition = "input.feature_type == 'gene'",
      selectizeInput(
        "gene",
        "Gene",
        choices = NULL,
        selected = NULL,
        options = list(placeholder = "Type a gene symbol")
      ),
      checkboxInput("log1p_gene", "Show log1p(expression)", value = FALSE)
    ),
    conditionalPanel(
      condition = "input.feature_type == 'metadata'",
      selectInput(
        "metadata_col",
        "Cell metadata column",
        choices = obs_annotation_names,
        selected = metadata_default
      )
    ),
    conditionalPanel(
      condition = "input.feature_type != 'metadata'",
      sliderInput(
        "upper_clip",
        "Upper color clip quantile",
        min = 0.90,
        max = 1.00,
        value = 0.99,
        step = 0.005
      )
    ),
    numericInput("point_size", "Point size", value = 0.12, min = 0.01, max = 2, step = 0.01),
    sliderInput("point_alpha", "Point alpha", min = 0.05, max = 1, value = 0.85, step = 0.05),
    annotation_filter_ui(),
    hr(),
    downloadButton("download_svg", "Download SVG"),
    downloadButton("download_png", "Download PNG")
  ),
  card(
    card_header(textOutput("plot_title")),
    plotOutput("umap_plot", height = "760px"),
    verbatimTextOutput("feature_summary")
  )
)

server <- function(input, output, session) {
  updateSelectizeInput(
    session,
    "gene",
    choices = gene_names,
    selected = if ("KDR" %in% gene_names) "KDR" else gene_names[1],
    server = TRUE
  )

  output$filter_values_ui <- renderUI({
    req(input$filter_column)
    if (identical(input$filter_column, "None")) {
      return(NULL)
    }
    vals <- sort(unique(umap_df[[input$filter_column]]))
    selectizeInput(
      "filter_values",
      "Values to include",
      choices = vals,
      selected = vals,
      multiple = TRUE,
      options = list(plugins = list("remove_button"))
    )
  })

  feature_data <- reactive({
    if (identical(input$feature_type, "program")) {
      req(input$program)
      value <- get_program_vector(input$program)
      label <- program_labels[match(input$program, program_names)]
      list(
        value = value,
        is_numeric = TRUE,
        feature_name = input$program,
        feature_label = paste0(input$program, " — ", label),
        legend_title = "Program score"
      )
    } else if (identical(input$feature_type, "gene")) {
      req(input$gene)
      value <- get_gene_vector(input$gene)
      if (isTRUE(input$log1p_gene)) {
        value <- log1p(value)
      }
      list(
        value = value,
        is_numeric = TRUE,
        feature_name = input$gene,
        feature_label = input$gene,
        legend_title = if (isTRUE(input$log1p_gene)) "log1p expression" else "Expression"
      )
    } else {
      req(input$metadata_col)
      if (!input$metadata_col %in% obs_annotation_names) {
        stop("Metadata column not found in exported data: ", input$metadata_col)
      }
      value <- as.character(umap_df[[input$metadata_col]])
      value[is.na(value) | value == ""] <- "NA"
      list(
        value = value,
        is_numeric = FALSE,
        feature_name = paste0("metadata_", input$metadata_col),
        feature_label = paste0("Cell metadata: ", input$metadata_col),
        legend_title = input$metadata_col
      )
    }
  })

  plot_df <- reactive({
    fd <- feature_data()
    df <- umap_df
    df$value <- fd$value
    if (isTRUE(fd$is_numeric)) {
      df$value_plot <- clip_values(df$value, lower_q = 0.00, upper_q = input$upper_clip)
    } else {
      df$value_plot <- factor(df$value)
    }

    if (!is.null(input$filter_column) && !identical(input$filter_column, "None")) {
      req(input$filter_values)
      df <- df[df[[input$filter_column]] %in% input$filter_values, , drop = FALSE]
    }
    df
  })

  output$plot_title <- renderText({
    feature_data()$feature_label
  })

  make_umap_plot <- reactive({
    fd <- feature_data()
    df <- plot_df()

    p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = value_plot)) +
      geom_point(size = input$point_size, alpha = input$point_alpha, stroke = 0) +
      coord_equal() +
      labs(title = fd$feature_label) +
      theme_void(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9)
      )

    if (isTRUE(fd$is_numeric)) {
      p + scale_color_viridis_c(option = "magma", name = fd$legend_title)
    } else {
      p + scale_color_discrete(name = fd$legend_title, drop = FALSE, na.translate = FALSE)
    }
  })

  output$umap_plot <- renderPlot({
    make_umap_plot()
  }, res = 140)

  output$feature_summary <- renderPrint({
    fd <- feature_data()
    df <- plot_df()
    x <- df$value
    cat("Feature:", fd$feature_label, "\n")
    cat("Cells shown:", nrow(df), "of", nrow(umap_df), "\n")

    if (isTRUE(fd$is_numeric)) {
      cat("Nonzero cells:", sum(x > 0, na.rm = TRUE), "\n")
      cat("Min:", signif(min(x, na.rm = TRUE), 4), "\n")
      cat("Mean:", signif(mean(x, na.rm = TRUE), 4), "\n")
      cat("Max:", signif(max(x, na.rm = TRUE), 4), "\n")
      cat("Upper clip quantile:", input$upper_clip, "\n")
    } else {
      tab <- sort(table(x, useNA = "ifany"), decreasing = TRUE)
      cat("Category counts:\n")
      print(tab)
    }
  })

  output$download_svg <- downloadHandler(
    filename = function() {
      paste0("umap_", safe_filename(feature_data()$feature_name), ".svg")
    },
    content = function(file) {
      grDevices::svg(file, width = 7.2, height = 6.2, onefile = FALSE)
      on.exit(grDevices::dev.off(), add = TRUE)
      print(make_umap_plot())
    }
  )

  output$download_png <- downloadHandler(
    filename = function() {
      paste0("umap_", safe_filename(feature_data()$feature_name), ".png")
    },
    content = function(file) {
      ggsave(file, plot = make_umap_plot(), width = 7.2, height = 6.2, dpi = 300)
    }
  )
}

shinyApp(ui, server)
