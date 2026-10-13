# make_hippo_d1_p13_p30_publication_ggplot2.R
# Real ggplot2 version of the Hippo D1 P13/P30 publication-style plots.
#
# Input:
#   Hippo_D1_D2_D3_scores_source.csv
#
# Outputs:
#   PNG and PDF files in ./ggplot2_outputs/
#
# To run:
#   install.packages(c("ggplot2", "dplyr", "tidyr", "readr", "ggrepel", "scales"))
#   source("make_hippo_d1_p13_p30_publication_ggplot2.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggrepel)
  library(scales)
})

input_csv <- "Hippo_D1_D2_D3_scores_source.csv"
outdir <- "ggplot2_outputs"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

df <- read_csv(input_csv, show_col_types = FALSE)

# Expected columns:
# target, D1_P13, D1_P30, D1_split, D2_P15, D2_P17, D2_split, D3_endpoint, set

hippo_core <- df %>%
  filter(set == "Hippo core") %>%
  arrange(D1_split)

hippo_labels <- hippo_core$target

theme_pub <- function(base_size = 15) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
      axis.title = element_text(size = base_size, face = "bold"),
      axis.text = element_text(size = base_size - 2, color = "black"),
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 1),
      legend.key.size = unit(0.85, "cm"),
      legend.background = element_rect(fill = "white", color = "grey40", linewidth = 0.5),
      legend.box.background = element_rect(fill = "white", color = "grey40", linewidth = 0.5),
      legend.margin = margin(8, 8, 8, 8),
      plot.margin = margin(12, 16, 12, 12)
    )
}

point_colors <- c(
  "Other targets" = "grey72",
  "Hippo core" = "#D95F02",
  "YAP/TEAD" = "#1B9E77"
)

point_shapes <- c(
  "Other targets" = 16,
  "Hippo core" = 8,
  "YAP/TEAD" = 18
)

point_sizes <- c(
  "Other targets" = 2.2,
  "Hippo core" = 5.5,
  "YAP/TEAD" = 4.2
)

save_plot <- function(p, name, width = 7.8, height = 6.3) {
  ggsave(file.path(outdir, paste0(name, ".png")), p, width = width, height = height, dpi = 450)
  ggsave(file.path(outdir, paste0(name, ".pdf")), p, width = width, height = height)
}

format_p <- function(p) {
  ifelse(p < 1e-3, scales::scientific(p, digits = 2), sprintf("%.3f", p))
}

scatter_hippo <- function(data, x, y, xlab, ylab, title, name) {
  xvec <- data[[x]]
  yvec <- data[[y]]
  ok <- is.finite(xvec) & is.finite(yvec)

  r <- cor(xvec[ok], yvec[ok], method = "pearson")
  pval <- cor.test(xvec[ok], yvec[ok], method = "pearson")$p.value
  n <- sum(ok)

  stat_label <- paste0(
    "Pearson r = ", sprintf("%.2f", r), "\n",
    "p = ", format_p(pval), "\n",
    "n = ", n
  )

  p <- ggplot(data, aes(x = .data[[x]], y = .data[[y]])) +
    geom_hline(yintercept = 0, linewidth = 0.45, linetype = "dashed", color = "grey45") +
    geom_vline(xintercept = 0, linewidth = 0.45, linetype = "dashed", color = "grey45") +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, color = "black") +
    geom_point(
      aes(color = set, shape = set, size = set),
      alpha = 0.85,
      stroke = 1.15
    ) +
    geom_text_repel(
      data = data %>% filter(target %in% hippo_labels),
      aes(label = target),
      size = 4.9,
      fontface = "bold",
      min.segment.length = 0,
      segment.size = 0.35,
      box.padding = 0.38,
      point.padding = 0.24,
      max.overlaps = Inf
    ) +
    annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = stat_label,
      hjust = -0.08,
      vjust = 1.12,
      size = 4.35,
      lineheight = 1.05
    ) +
    scale_color_manual(values = point_colors, breaks = c("Other targets", "Hippo core", "YAP/TEAD")) +
    scale_shape_manual(values = point_shapes, breaks = c("Other targets", "Hippo core", "YAP/TEAD")) +
    scale_size_manual(values = point_sizes, breaks = c("Other targets", "Hippo core", "YAP/TEAD")) +
    labs(x = xlab, y = ylab, title = title) +
    theme_pub(base_size = 15) +
    guides(
      color = guide_legend(override.aes = list(size = 5.0, alpha = 1)),
      shape = guide_legend(override.aes = list(size = 5.0, alpha = 1)),
      size = "none"
    )

  save_plot(p, name, width = 8.2, height = 6.6)
  return(p)
}

p1 <- scatter_hippo(
  df,
  x = "D1_P13",
  y = "D1_P30",
  xlab = "D1 P13 perturbation effect",
  ylab = "D1 P30 perturbation effect",
  title = "Hippo targets on the D1 P13/P30 plane",
  name = "01_D1_P13_vs_P30_ggplot2"
)

p2 <- scatter_hippo(
  df,
  x = "D1_split",
  y = "D2_split",
  xlab = "D1 split effect: P30 - P13",
  ylab = "D2 split effect: P15 - P17",
  title = "Hippo-associated D1 shift aligns with the D2 fate axis",
  name = "02_D1_split_vs_D2_split_ggplot2"
)

p3 <- scatter_hippo(
  df,
  x = "D1_split",
  y = "D3_endpoint",
  xlab = "D1 split effect: P30 - P13",
  ylab = "D3 EC endpoint score",
  title = "Hippo-associated D1 shift foreshadows D3 endpoint loss",
  name = "03_D1_split_vs_D3_endpoint_ggplot2"
)

hippo_long <- hippo_core %>%
  select(target, D1_P13, D1_P30, D1_split) %>%
  pivot_longer(
    cols = c(D1_P13, D1_P30, D1_split),
    names_to = "score",
    values_to = "effect"
  ) %>%
  mutate(
    target = factor(target, levels = hippo_core$target),
    score = factor(
      score,
      levels = c("D1_P13", "D1_P30", "D1_split"),
      labels = c("D1 P13", "D1 P30", "D1 split: P30 - P13")
    )
  )

p4 <- ggplot(hippo_long, aes(x = effect, y = target, fill = score)) +
  geom_vline(xintercept = 0, linewidth = 0.5, linetype = "dashed", color = "grey45") +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  labs(
    x = "Perturbation effect",
    y = NULL,
    title = "D1 P13/P30 effects for Hippo-core targets"
  ) +
  theme_pub(base_size = 15) +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(face = "bold", size = 14)
  )

save_plot(p4, "04_Hippo_D1_P13_P30_barplot_ggplot2", width = 8.4, height = 5.5)

hist_counts <- hist(df$D1_split, breaks = 35, plot = FALSE)$counts
ymax <- max(hist_counts, na.rm = TRUE)

p5 <- ggplot(df, aes(x = D1_split)) +
  geom_histogram(bins = 35, fill = "grey72", color = "white", linewidth = 0.3) +
  geom_vline(xintercept = 0, linewidth = 0.45, linetype = "dashed", color = "grey45") +
  geom_vline(
    data = hippo_core,
    aes(xintercept = D1_split),
    linewidth = 0.8,
    color = "#D95F02"
  ) +
  geom_text_repel(
    data = hippo_core,
    aes(x = D1_split, y = ymax * 0.96, label = target),
    angle = 90,
    size = 4.7,
    fontface = "bold",
    direction = "x",
    segment.size = 0.3,
    max.overlaps = Inf
  ) +
  labs(
    x = "D1 split effect: P30 - P13",
    y = "Number of targets",
    title = "Hippo-core targets on the D1 fate-split distribution"
  ) +
  theme_pub(base_size = 15) +
  theme(legend.position = "none")

save_plot(p5, "05_D1_split_distribution_Hippo_overlay_ggplot2", width = 8.2, height = 6.0)

heat_df <- hippo_core %>%
  select(target, D1_P13, D1_P30, D1_split, D2_P15, D2_P17, D2_split, D3_endpoint) %>%
  mutate(target = factor(target, levels = hippo_core$target)) %>%
  pivot_longer(
    cols = -target,
    names_to = "score",
    values_to = "effect"
  ) %>%
  mutate(
    score = factor(
      score,
      levels = c("D1_P13", "D1_P30", "D1_split", "D2_P15", "D2_P17", "D2_split", "D3_endpoint"),
      labels = c("D1 P13", "D1 P30", "D1 P30-P13", "D2 P15", "D2 P17", "D2 P15-P17", "D3 EC endpoint")
    )
  )

max_abs <- max(abs(heat_df$effect), na.rm = TRUE)

p6 <- ggplot(heat_df, aes(x = score, y = target, fill = effect)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.2f", effect)), size = 3.9, fontface = "bold") +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-max_abs, max_abs),
    name = "Effect"
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Hippo-core target effects across D1, D2, and D3"
  ) +
  theme_pub(base_size = 15) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(face = "bold", size = 14),
    legend.position = "right"
  )

save_plot(p6, "06_Hippo_D1_D2_D3_score_matrix_ggplot2", width = 9.3, height = 5.5)

message("Wrote ggplot2 figures to: ", outdir)
