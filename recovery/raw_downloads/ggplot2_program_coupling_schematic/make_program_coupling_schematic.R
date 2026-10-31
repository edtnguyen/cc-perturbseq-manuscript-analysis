suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

# =========================
# Inputs
# =========================
indir  <- "/Users/edtng/Downloads/ggplot2_program_coupling_schematic"
outdir <- indir

nodes_path <- file.path(indir, "program_coupling_nodes.csv")
edges_path <- file.path(indir, "program_coupling_edges.csv")

nodes <- read.csv(nodes_path, stringsAsFactors = FALSE)
edges <- read.csv(edges_path, stringsAsFactors = FALSE)

# =========================
# Helper data
# =========================
node_from <- nodes[, c("node_id", "x", "y")]
colnames(node_from) <- c("from_node", "x_from", "y_from")
node_to <- nodes[, c("node_id", "x", "y")]
colnames(node_to) <- c("to_node", "x_to", "y_to")

edges <- merge(edges, node_from, by = "from_node", all.x = TRUE)
edges <- merge(edges, node_to,   by = "to_node",   all.x = TRUE)

edges$label_x <- (edges$x_from + edges$x_to) / 2 + edges$label_nudge_x
edges$label_y <- (edges$y_from + edges$y_to) / 2 + edges$label_nudge_y

# Manual day-band positions. Adjust if you move node x coordinates.
day_bands <- data.frame(
  day  = c("D0", "D1", "D2", "D3"),
  xmin = c(-0.55, 0.85, 3.25, 5.65),
  xmax = c( 0.55, 2.55, 4.95, 7.35),
  ymin = c(-2.05, -2.05, -2.05, -2.05),
  ymax = c( 2.05,  2.05,  2.05,  2.05),
  fill = c("#F7F7F7", "#FFF8F2", "#F6FAFE", "#FBFBFB")
)

day_text <- data.frame(
  day = c("D0", "D1", "D2", "D3"),
  x   = c(0.0, 1.7, 4.1, 6.5),
  y   = c(2.38, 2.38, 2.38, 2.38)
)

nodes_box <- subset(nodes, node_id != "iPSC")
node_w <- 1.42
node_h <- 0.96

# =========================
# Plot
# =========================
p <- ggplot() +
  # light day bands
  geom_rect(
    data = day_bands,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = day_bands$fill,
    color = NA
  ) +
  # subtle separators between days
  geom_vline(xintercept = c(0.70, 2.90, 5.30), color = "#E3E3E3", linewidth = 0.45) +
  # developmental/correlation arrows
  geom_curve(
    data = edges,
    aes(
      x = x_from, y = y_from,
      xend = x_to, yend = y_to,
      color = color,
      linewidth = size,
      alpha = alpha
    ),
    curvature = edges$curvature,
    lineend = "round",
    arrow = arrow(length = unit(0.16, "inches"), type = "closed")
  ) +
  # r labels on colored links
  geom_label(
    data = subset(edges, show_label == 1),
    aes(x = label_x, y = label_y, label = label, fill = color),
    label.size = 0,
    label.padding = unit(0.10, "lines"),
    size = 4.0,
    fontface = "bold",
    color = "#3A3A3A"
  ) +
  # node rectangles
  geom_tile(
    data = nodes_box,
    aes(x = x, y = y),
    width = node_w,
    height = node_h,
    fill = nodes_box$fill_color,
    color = nodes_box$border_color,
    linewidth = 0.75
  ) +
  # short labels
  geom_text(
    data = nodes_box,
    aes(x = x, y = y + 0.16, label = short_label, color = text_color),
    size = 10,
    fontface = "bold",
    lineheight = 0.9
  ) +
  # long labels
  geom_text(
    data = nodes_box,
    aes(x = x, y = y - 0.13, label = long_label, color = text_color),
    size = 4.9,
    lineheight = 0.92
  ) +
  # starting state label
  geom_text(
    data = subset(nodes, node_id == "iPSC"),
    aes(x = x, y = y, label = short_label),
    size = 11,
    fontface = "bold",
    color = "#333333"
  ) +
  # day labels
  geom_text(
    data = day_text,
    aes(x = x, y = y, label = day),
    size = 12,
    fontface = "bold",
    color = "#111111"
  ) +
  # figure annotations
  annotate(
    "text",
    x = 3.25, y = 2.90,
    label = "Hypothesized coupling of program effects across D1–D3 endothelial differentiation",
    size = 7.0, fontface = "bold", color = "#111111"
  ) +
  annotate(
    "text",
    x = 3.25, y = 2.63,
    label = "Colored links show Pearson correlation of target-level perturbation effects between selected programs",
    size = 4.3, color = "#4B4B4B"
  ) +
  # small legend-like guide
  annotate("segment", x = 5.95, xend = 6.32, y = -2.50, yend = -2.50,
           linewidth = 1.5, color = "#EB8C63", alpha = 0.95,
           arrow = arrow(length = unit(0.12, "inches"), type = "closed")) +
  annotate("text", x = 6.80, y = -2.50, label = "positive coupling", hjust = 0,
           size = 4.0, color = "#555555") +
  annotate("segment", x = 5.95, xend = 6.32, y = -2.77, yend = -2.77,
           linewidth = 1.5, color = "#4A9AD4", alpha = 0.95,
           arrow = arrow(length = unit(0.12, "inches"), type = "closed")) +
  annotate("text", x = 6.80, y = -2.77, label = "negative coupling", hjust = 0,
           size = 4.0, color = "#555555") +
  annotate("text", x = 0.0, y = -2.77, label = "line thickness ~ |r|", hjust = 0,
           size = 4.0, color = "#666666") +
  scale_color_identity() +
  scale_fill_identity() +
  scale_linewidth_identity() +
  scale_alpha_identity() +
  coord_cartesian(xlim = c(-0.75, 7.55), ylim = c(-2.95, 3.05), clip = "off") +
  theme_void(base_size = 12) +
  theme(
    plot.margin = margin(t = 12, r = 18, b = 8, l = 16),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  )

# =========================
# Export
# =========================
out_png <- file.path(outdir, "program_coupling_schematic_ggplot2.png")
out_pdf <- file.path(outdir, "program_coupling_schematic_ggplot2.pdf")
out_svg <- file.path(outdir, "program_coupling_schematic_ggplot2.svg")

ggsave(out_png, p, width = 12.4, height = 5.5, dpi = 350, bg = "white")
ggsave(out_pdf, p, width = 12.4, height = 5.5, bg = "white", device = cairo_pdf)

# SVG export if svglite is available
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(out_svg, p, width = 12.4, height = 5.5, bg = "white", device = svglite::svglite)
}

message("Saved: ")
message(out_png)
message(out_pdf)
if (file.exists(out_svg)) message(out_svg)
