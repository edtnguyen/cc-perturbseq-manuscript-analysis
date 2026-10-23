
library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)

out_dir <- "/mnt/data/unbiased_hippo_d3_live"

target_rank <- read_csv(file.path(out_dir, "source_figure1_target_rank.csv"), show_col_types = FALSE)
pca_df <- read_csv(file.path(out_dir, "source_figure3_target_effectspace_pca.csv"), show_col_types = FALSE)
heat_long <- read_csv(file.path(out_dir, "source_figure4_top40_targets_all_programs_heatmap.csv"), show_col_types = FALSE)
null_df <- read_csv(file.path(out_dir, "source_figure5_permutation_null.csv"), show_col_types = FALSE)
sum_df <- read_csv(file.path(out_dir, "source_figure5_permutation_summary.csv"), show_col_types = FALSE)

hippo_color <- "#D55E00"
other_color <- "#BDBDBD"

# Fig1
target_rank$rank <- seq_len(nrow(target_rank))
p1 <- ggplot(target_rank, aes(x = global_d3_burden_score, y = rank)) +
  geom_point(aes(color = is_hippo), size = 2.2, alpha = 0.9) +
  scale_color_manual(values = c("TRUE" = hippo_color, "FALSE" = other_color), guide = "none") +
  scale_y_reverse() +
  geom_text_repel(data = subset(target_rank, target_name %in% c(head(target_rank$target_name,15), target_rank$target_name[target_rank$is_hippo])),
                  aes(label = target_name, color = is_hippo), size = 2.8, show.legend = FALSE) +
  labs(title = "Unbiased D3 target burden rank",
       x = "Global D3 burden score (sum of |significant program effects|)",
       y = "Target rank") +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "figure1_unbiased_d3_target_rank_fromR.png"), p1, width = 7, height = 8, dpi = 300)

# Fig3
p3 <- ggplot(pca_df, aes(PC1, PC2)) +
  geom_point(data = subset(pca_df, !is_hippo), color = other_color, size = 1.8, alpha = 0.65) +
  geom_point(data = subset(pca_df, is_hippo), color = hippo_color, size = 2.8) +
  geom_text_repel(data = subset(pca_df, is_hippo | target_name %in% head(target_rank$target_name,10)),
                  aes(label = target_name, color = is_hippo), size = 2.8, show.legend = FALSE) +
  scale_color_manual(values = c("TRUE" = hippo_color, "FALSE" = "#4D4D4D"), guide = "none") +
  labs(title = "Target effect-space PCA (all targets, D3)",
       x = paste0("PC1 (14.2% var)"),
       y = paste0("PC2 (10.1% var)")) +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "figure3_target_effect_space_pca_fromR.png"), p3, width = 8.2, height = 6.2, dpi = 300)

# Fig4
heat_long$target_name <- factor(heat_long$target_name, levels = rev(unique(heat_long$target_name)))
heat_long$program <- factor(heat_long$program, levels = unique(heat_long$program))
p4 <- ggplot(heat_long, aes(program, target_name, fill = row_scaled_effect)) +
  geom_tile() +
  scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0, limits = c(-2.5, 2.5)) +
  labs(title = "Top targets × all programs heatmap (unsupervised)",
       x = NULL, y = NULL, fill = "Row-scaled\nprogram effect") +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(file.path(out_dir, "figure4_top_targets_all_programs_heatmap_fromR.png"), p4, width = 12, height = 8.5, dpi = 300)

# Fig5
obs <- sum_df$observed_mean_pairwise_distance[1]
p5 <- ggplot(null_df, aes(x = null_mean_pairwise_distance)) +
  geom_histogram(bins = 50, fill = "lightgray", color = "white") +
  geom_vline(xintercept = obs, color = hippo_color, linewidth = 1.2) +
  annotate("text", x = obs, y = Inf, label = paste0("Observed Hippo\nmean distance = ", round(obs,2), "\nempirical p = ", signif(sum_df$empirical_p[1],3)),
           hjust = -0.05, vjust = 1.2, color = hippo_color, size = 3.3) +
  labs(title = "Hippo compactness permutation test",
       x = "Mean pairwise Euclidean distance in D3 effect space",
       y = "Random gene-set count") +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "figure5_hippo_compactness_permutation_test_fromR.png"), p5, width = 8, height = 5.6, dpi = 300)
