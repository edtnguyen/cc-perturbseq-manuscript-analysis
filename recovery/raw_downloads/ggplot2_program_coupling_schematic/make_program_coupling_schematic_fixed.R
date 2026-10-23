suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

# ============================================================
# Robust path handling
# ============================================================
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  hit <- grep(file_arg, args, value = TRUE)
  if (length(hit) > 0) {
    return(dirname(normalizePath(sub(file_arg, "", hit[1]))))
  }
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  }
  return(getwd())
}

indir  <- get_script_dir()
outdir <- indir

nodes_path <- file.path(indir, "program_coupling_nodes.csv")
edges_path <- file.path(indir, "program_coupling_edges.csv")

nodes <- read.csv(nodes_path, stringsAsFactors = FALSE)
edges <- read.csv(edges_path, stringsAsFactors = FALSE)

# Convert pipe separators to line breaks.
nodes$long_label <- gsub("\\|", "\n", nodes$long_label)

# Join node coordinates to edges.
node_from <- nodes[, c("node_id", "x", "y")]
colnames(node_from) <- c("from_node", "x_from", "y_from")
node_to <- nodes[, c("node_id", "x", "y")]
colnames(node_to) <- c("to_node", "x_to", "y_to")

edges <- merge(edges, node_from, by = "from_node", all.x = TRUE)
edges <- merge(edges, node_to,   by = "to_node",   all.x = TRUE)

edges$label_x <- (edges$x_from + edges$x_to) / 2 + edges$label_nudge_x
edges$label_y <- (edges$y_from + edges$y_to) / 2 + edges$label_nudge_y
edges$linetype <- ifelse(is.na(edges$linetype), "solid", edges$linetype)

# Manual day-band positions.
day_bands <- data.frame(
  day  = c("D0", "D1", "D2", "D3"),
  xmin = c(-0.55, 0.78, 3.00, 5.32),
  xmax = c( 0.55, 2.32, 4.50, 6.78),
  ymin = c(-1.95, -1.95, -1.95, -1.95),
  ymax = c( 1.95,  1.95,  1.95,  1.95),
  fill = c("#F8F8F8", "#FFF8F2", "#F6FAFE", "#FBFBFB")
)

day_text <- data.frame(
  day = c("D0", "D1", "D2", "D3"),
  x   = c(0.0, 1.55, 3.75, 6.05),
  y   = c(2.32, 2.32, 2.32, 2.32)
)

nodes_start <- subset(nodes, node_id == "iPSC")
nodes_prog  <- subset(nodes, node_id != "iPSC")

# ============================================================
# Plot
# ============================================================
p <- ggplot() +
  geom_rect(
    data = day_bands,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
    color = NA,
    alpha = 0.65
  ) +
  geom_vline(xintercept = c(0.68, 2.66, 5.05), color = "#E5E5E5", linewidth = 0.45) +

  # Edges. This intentionally uses geom_segment instead of geom_curve.
  # Previous error came from passing a vector to geom_curve(curvature = ...).
  geom_segment(
    data = edges,
    aes(
      x = x_from, y = y_from,
      xend = x_to, yend = y_to,
      color = color,
      linewidth = linewidth,
      alpha = alpha,
      linetype = linetype
    ),
    lineend = "round",
    arrow = arrow(length = unit(0.15, "inches"), type = "closed")
  ) +

  # Correlation labels.
  geom_label(
    data = subset(edges, show_label == 1),
    aes(x = label_x, y = label_y, label = label, color = color),
    fill = "white",
    linewidth = 0,
    label.padding = unit(0.10, "lines"),
    size = 3.7,
    fontface = "bold"
  ) +

  # Nodes: rectangles.
  geom_rect(
    data = nodes,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_color, color = border_color),
    linewidth = 0.75
  ) +

  # Program short labels.
  geom_text(
    data = nodes_prog,
    aes(x = x, y = y + 0.15, label = short_label, color = text_color),
    size = 6.8,
    fontface = "bold",
    lineheight = 0.9
  ) +

  # Program long labels.
  geom_text(
    data = nodes_prog,
    aes(x = x, y = y - 0.13, label = long_label, color = text_color),
    size = 3.6,
    lineheight = 0.92
  ) +

  # iPSC label.
  geom_text(
    data = nodes_start,
    aes(x = x, y = y, label = short_label),
    size = 7.0,
    fontface = "bold",
    color = "#333333"
  ) +

  # Day labels.
  geom_text(
    data = day_text,
    aes(x = x, y = y, label = day),
    size = 7.5,
    fontface = "bold",
    color = "#111111"
  ) +

  annotate(
    "text",
    x = 3.25, y = 2.80,
    label = "Hypothesized coupling of program effects across D1–D3 endothelial differentiation",
    size = 5.2,
    fontface = "bold",
    color = "#111111"
  ) +
  annotate(
    "text",
    x = 3.25, y = 2.56,
    label = "Links show Pearson correlation of target-level perturbation effects between selected programs",
    size = 3.5,
    color = "#4B4B4B"
  ) +

  # Legend-like guide.
  annotate("segment", x = 5.45, xend = 5.82, y = -2.42, yend = -2.42,
           linewidth = 1.4, color = "#EB8C63", alpha = 0.95,
           arrow = arrow(length = unit(0.11, "inches"), type = "closed")) +
  annotate("text", x = 5.95, y = -2.42, label = "positive coupling", hjust = 0,
           size = 3.4, color = "#555555") +
  annotate("segment", x = 5.45, xend = 5.82, y = -2.67, yend = -2.67,
           linewidth = 1.4, color = "#4A9AD4", alpha = 0.95,
           arrow = arrow(length = unit(0.11, "inches"), type = "closed")) +
  annotate("text", x = 5.95, y = -2.67, label = "negative coupling", hjust = 0,
           size = 3.4, color = "#555555") +
  annotate("text", x = -0.30, y = -2.67, label = "line thickness ~ |r|; dashed = weaker link", hjust = 0,
           size = 3.4, color = "#666666") +

  scale_color_identity() +
  scale_fill_identity() +
  scale_linewidth_identity() +
  scale_alpha_identity() +
  scale_linetype_identity() +
  coord_cartesian(xlim = c(-0.72, 7.05), ylim = c(-2.86, 2.96), clip = "off") +
  theme_void(base_size = 12) +
  theme(
    plot.margin = margin(t = 12, r = 18, b = 8, l = 16),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  )

# ============================================================
# Export
# ============================================================
out_png <- file.path(outdir, "program_coupling_schematic_ggplot2_fixed.png")
out_pdf <- file.path(outdir, "program_coupling_schematic_ggplot2_fixed.pdf")
out_svg <- file.path(outdir, "program_coupling_schematic_ggplot2_fixed.svg")

ggsave(out_png, p, width = 12.4, height = 5.5, dpi = 350, bg = "white")
ggsave(out_pdf, p, width = 12.4, height = 5.5, bg = "white")

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(out_svg, p, width = 12.4, height = 5.5, bg = "white", device = svglite::svglite)
}

message("Saved:")
message(out_png)
message(out_pdf)
if (file.exists(out_svg)) message(out_svg)
