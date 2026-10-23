
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(tidyr)
})

out_dir <- "."
pal <- c(
  "Hippo/YAP" = "#6A3D9A",
  "TGF/BMP-SMAD" = "#D55E00",
  "EC/vascular" = "#0072B2",
  "Mesoderm/dev TF" = "#009E73",
  "Chromatin/TF" = "#8C564B",
  "RNA/translation" = "#4DBBD5",
  "Apoptosis/stress" = "#E64B35",
  "Trafficking/metabolism" = "#7E6148",
  "DNA repair/cell cycle" = "#BCBD22",
  "Signaling/adaptors" = "#F0E442",
  "Other" = "#BDBDBD"
)

theme_pub <- theme_bw(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
    axis.text = element_text(color = "black"),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )

target <- read_csv("source_fig1_target_global_burden_rank.csv", show_col_types = FALSE)
enrich <- read_csv("source_fig3_tophit_enrichment.csv", show_col_types = FALSE)
coh <- read_csv("source_fig4_pathway_strength_vs_coherence.csv", show_col_types = FALSE)
heat <- read_csv("source_fig5_top_targets_all_program_heatmap.csv", show_col_types = FALSE)

# Figure 1
fig1 <- target %>%
  arrange(desc(global_burden_score)) %>%
  slice_head(n = 45) %>%
  mutate(target_name = factor(target_name, levels = rev(target_name))) %>%
  ggplot(aes(x = global_burden_score, y = target_name, fill = pathway_class)) +
  geom_col(color = "white", linewidth = 0.2) +
  scale_fill_manual(values = pal, drop = FALSE) +
  labs(
    title = "Global D3 target ranking identifies Hippo/YAP among strongest perturbation classes",
    x = "Global D3 perturbation burden\nΣ |significant program effects| across all cNMF programs",
    y = NULL,
    fill = "Regulator class"
  ) +
  theme_pub
ggsave("figure1_global_d3_target_burden_rank_R.png", fig1, width = 7.0, height = 9.0, dpi = 400)
ggsave("figure1_global_d3_target_burden_rank_R.pdf", fig1, width = 7.0, height = 9.0)
ggsave("figure1_global_d3_target_burden_rank_R.svg", fig1, width = 7.0, height = 9.0)

# Figure 2
class_order <- target %>% group_by(pathway_class) %>% summarize(med = median(global_burden_score), .groups="drop") %>% arrange(desc(med)) %>% pull(pathway_class)
fig2 <- target %>%
  mutate(pathway_class = factor(pathway_class, levels = class_order)) %>%
  ggplot(aes(x = pathway_class, y = global_burden_score, fill = pathway_class)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.45, linewidth = 0.25) +
  geom_jitter(aes(color = pathway_class), width = 0.16, height = 0, size = 1.7, alpha = 0.75, show.legend = FALSE) +
  scale_fill_manual(values = pal, drop = FALSE) +
  scale_color_manual(values = pal, drop = FALSE) +
  labs(title = "Pathway-level distribution of global D3 perturbation burden",
       x = NULL, y = "Global D3 perturbation burden", fill = "Regulator class") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "none")
ggsave("figure2_pathway_global_burden_distribution_R.png", fig2, width = 8.0, height = 4.8, dpi = 400)
ggsave("figure2_pathway_global_burden_distribution_R.pdf", fig2, width = 8.0, height = 4.8)
ggsave("figure2_pathway_global_burden_distribution_R.svg", fig2, width = 8.0, height = 4.8)

# Figure 3
fig3 <- enrich %>% filter(n_top20 > 0) %>%
  arrange(fold_enrichment_top20) %>%
  mutate(label = paste0(pathway_class, " (n=", n_top20, ")"),
         label = factor(label, levels = label)) %>%
  ggplot(aes(x = fold_enrichment_top20, y = label, fill = pathway_class)) +
  geom_col(color = "white") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_fill_manual(values = pal, drop = FALSE) +
  labs(title = "Hippo/YAP is enriched among the strongest global D3 hits",
       x = "Fold enrichment among top 20 global D3 hits", y = NULL) +
  theme_pub +
  theme(legend.position = "none")
ggsave("figure3_top20_pathway_enrichment_R.png", fig3, width = 6.7, height = 4.2, dpi = 400)
ggsave("figure3_top20_pathway_enrichment_R.pdf", fig3, width = 6.7, height = 4.2)
ggsave("figure3_top20_pathway_enrichment_R.svg", fig3, width = 6.7, height = 4.2)

# Figure 4
fig4 <- coh %>% filter(!is.na(within_pathway_coherence), pathway_class != "Other") %>%
  ggplot(aes(x = mean_global_burden, y = within_pathway_coherence, color = pathway_class, size = n_targets, label = pathway_class)) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_point(alpha = 0.9) +
  geom_text_repel(size = 2.4, show.legend = FALSE, max.overlaps = Inf) +
  scale_color_manual(values = pal, drop = FALSE) +
  scale_size_continuous(range = c(2,6)) +
  labs(title = "Hippo/YAP combines strong global burden with high pathway coherence",
       x = "Mean global D3 perturbation burden",
       y = "Within-pathway coherence\nmean pairwise correlation across all 50 program effects",
       color = "Regulator class", size = "# targets") +
  theme_pub
ggsave("figure4_pathway_strength_vs_coherence_R.png", fig4, width = 6.5, height = 4.9, dpi = 400)
ggsave("figure4_pathway_strength_vs_coherence_R.pdf", fig4, width = 6.5, height = 4.9)
ggsave("figure4_pathway_strength_vs_coherence_R.svg", fig4, width = 6.5, height = 4.9)

# Figure 5
top_targets <- target %>% arrange(desc(global_burden_score)) %>% slice_head(n = 35) %>% pull(target_name)
program_order <- heat %>% distinct(program_name, program_label) %>% pull(program_label)
heat2 <- heat %>%
  filter(target_name %in% top_targets) %>%
  mutate(target_label = paste0(rank, ". ", target_name),
         target_label = factor(target_label, levels = rev(unique(paste0(target %>% filter(target_name %in% top_targets) %>% arrange(desc(global_burden_score)) %>% pull(rank), ". ", top_targets)))),
         program_label = factor(program_label, levels = unique(program_label)),
         effect_clipped = pmax(pmin(sig_effect, 4), -4))
fig5 <- ggplot(heat2, aes(x = program_label, y = target_label, fill = effect_clipped)) +
  geom_tile() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-4,4)) +
  labs(title = "Unsupervised all-program D3 perturbation signatures of top global targets",
       x = "cNMF programs", y = NULL, fill = "Significant\neffect") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
        legend.position = "right")
ggsave("figure5_unsupervised_top_targets_all_program_heatmap_R.png", fig5, width = 10.0, height = 7.0, dpi = 400)
ggsave("figure5_unsupervised_top_targets_all_program_heatmap_R.pdf", fig5, width = 10.0, height = 7.0)
ggsave("figure5_unsupervised_top_targets_all_program_heatmap_R.svg", fig5, width = 10.0, height = 7.0)
