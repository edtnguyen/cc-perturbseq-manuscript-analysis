
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

out_dir <- "/mnt/data/tgfb_bmp_d2_revision"

group_cols <- c(
  "TGFβ/Activin arm" = "#7B3294",
  "BMP arm" = "#008837",
  "Shared TGF/BMP-SMAD" = "#D95F02",
  "Hippo/YAP" = "#1F78B4",
  "EC specification" = "#33A02C",
  "Chromatin/cofactor" = "#6A3D9A",
  "RNA processing/translation" = "#A6CEE3",
  "Apoptosis/stress" = "#E31A1C",
  "Trafficking/metabolism" = "#FDBF6F",
  "Other" = "#BDBDBD"
)

theme_pub <- theme_bw(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
    strip.background = element_rect(fill = "grey95", color = "grey80"),
    strip.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )

# Figure 2: D2 state map
state <- read_csv(file.path(out_dir, "figure2_d2_state_map_source.csv"), show_col_types = FALSE)
label_genes <- c("ACVR1B","ACVR2A","ACVR2B","BMPR1A","SMAD1","SMAD3","SMAD4","SMAD5","TGFBR2",
                 "NF2","LATS2","AMOTL2","ETV2","KDR","MED12","SETDB1","HAND1")
fig2 <- ggplot(state, aes(x = D2_EC_induction, y = D2_P19_minus_P14P17)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_point(aes(fill = P19, color = regulator_group), shape = 21, size = 3, stroke = 0.9, alpha = 0.9) +
  geom_text_repel(data = filter(state, target_name %in% label_genes),
                  aes(label = target_name), size = 2.4, max.overlaps = Inf, seed = 1) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       name = "D2 P19 effect") +
  scale_color_manual(values = group_cols, name = "Regulator class") +
  labs(title = "D2 TGF/BMP perturbation-state map",
       x = "D2 EC-induction effect\nmean(P15, P23)",
       y = "D2 posterior-skew effect\nP19 − mean(P14, P17)") +
  theme_pub
ggsave(file.path(out_dir, "figure2_d2_tgfb_bmp_state_map_R.png"), fig2, width = 6.2, height = 4.8, dpi = 400)
ggsave(file.path(out_dir, "figure2_d2_tgfb_bmp_state_map_R.pdf"), fig2, width = 6.2, height = 4.8)
ggsave(file.path(out_dir, "figure2_d2_tgfb_bmp_state_map_R.svg"), fig2, width = 6.2, height = 4.8)

# Figure 3: D2 heatmap
heat <- read_csv(file.path(out_dir, "figure3_d2_selected_program_heatmap_source.csv"), show_col_types = FALSE)
heat_long <- heat %>%
  pivot_longer(cols = -target_name, names_to = "program", values_to = "effect") %>%
  mutate(program = factor(program, levels = c("P19","P14","P17","P15","P23","P34_D2"),
                          labels = c("P19 posterior\nmesoderm","P14 outflow/\nBMP-like","P17 Lat-Meso/\nepithelial","P15 endothelial","P23 EC induction","P34 YAP/\ncontractile")),
         target_name = factor(target_name, levels = unique(heat$target_name)))
fig3 <- ggplot(heat_long, aes(x = program, y = target_name, fill = effect)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(abs(effect) >= 0.75, sprintf("%.1f", effect), "")), size = 2.1) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       name = "Significant effect\n(0 if NS)") +
  labs(title = "D2 program effects distinguish posterior, lateral/outflow, and EC-induction axes",
       x = NULL, y = NULL) +
  theme_pub +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(out_dir, "figure3_d2_selected_program_heatmap_R.png"), fig3, width = 7.2, height = 5.2, dpi = 400)
ggsave(file.path(out_dir, "figure3_d2_selected_program_heatmap_R.pdf"), fig3, width = 7.2, height = 5.2)
ggsave(file.path(out_dir, "figure3_d2_selected_program_heatmap_R.svg"), fig3, width = 7.2, height = 5.2)

# Figure 4: rank
bar <- read_csv(file.path(out_dir, "figure4_d2_failure_rank_displayed_source.csv"), show_col_types = FALSE) %>%
  mutate(target_name = factor(target_name, levels = target_name))
fig4 <- ggplot(bar, aes(x = D2_P19_minus_EC_induction, y = target_name, fill = regulator_group)) +
  geom_col(color = "black", linewidth = 0.25) +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
  scale_fill_manual(values = group_cols) +
  labs(title = "D2 P19-vs-EC failure score",
       x = "P19 − mean(P15, P23)", y = NULL, fill = "Regulator class") +
  theme_pub
ggsave(file.path(out_dir, "figure4_d2_p19_vs_ec_failure_rank_R.png"), fig4, width = 6.2, height = 6.0, dpi = 400)
ggsave(file.path(out_dir, "figure4_d2_p19_vs_ec_failure_rank_R.pdf"), fig4, width = 6.2, height = 6.0)
ggsave(file.path(out_dir, "figure4_d2_p19_vs_ec_failure_rank_R.svg"), fig4, width = 6.2, height = 6.0)

# Figure 5: lead-lag
lead <- read_csv(file.path(out_dir, "figure5_d2_failure_vs_d3_endpoint_source.csv"), show_col_types = FALSE)
lead_long <- lead %>%
  pivot_longer(cols = c(D2_P19_minus_EC_induction, D2_P19_minus_P14P17),
               names_to = "score_type", values_to = "D2_score") %>%
  mutate(score_type = factor(score_type,
                             levels = c("D2_P19_minus_EC_induction","D2_P19_minus_P14P17"),
                             labels = c("P19 − mean(P15,P23)", "P19 − mean(P14,P17)")))
fig5 <- ggplot(lead_long, aes(x = D2_score, y = D3_EC_endpoint_score, color = regulator_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_point(size = 2.5, alpha = 0.85) +
  geom_text_repel(data = filter(lead_long, target_name %in% label_genes),
                  aes(label = target_name), size = 2.2, max.overlaps = Inf, seed = 1) +
  facet_wrap(~score_type, nrow = 1, scales = "free_x") +
  scale_color_manual(values = group_cols) +
  labs(title = "D2 posterior-skew signatures are associated with later EC endpoint loss",
       x = "D2 score", y = "D3 EC endpoint score\nmean(P2,P3,P8,P11) − mean(P27,P34,P35)",
       color = "Regulator class") +
  theme_pub
ggsave(file.path(out_dir, "figure5_d2_failure_vs_d3_endpoint_R.png"), fig5, width = 10.5, height = 4.5, dpi = 400)
ggsave(file.path(out_dir, "figure5_d2_failure_vs_d3_endpoint_R.pdf"), fig5, width = 10.5, height = 4.5)
ggsave(file.path(out_dir, "figure5_d2_failure_vs_d3_endpoint_R.svg"), fig5, width = 10.5, height = 4.5)

# Figure 6: state heatmap
stateh <- read_csv(file.path(out_dir, "figure6_stage_specific_state_signature_heatmap_source.csv"), show_col_types = FALSE)
stateh_long <- stateh %>% pivot_longer(cols = -target_name, names_to = "state", values_to = "value")
state_levels <- colnames(stateh)[-1]
stateh_long$state <- factor(stateh_long$state, levels = state_levels)
stateh_long$target_name <- factor(stateh_long$target_name, levels = unique(stateh$target_name))
fig6 <- ggplot(stateh_long, aes(x = state, y = target_name, fill = value)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = ifelse(abs(value) >= 1, sprintf("%.1f", value), "")), size = 1.8) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       name = "Effect/score") +
  labs(title = "Stage-specific perturbation signatures separate D2 transition and D3 endpoint failure modes",
       x = NULL, y = NULL) +
  theme_pub +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(out_dir, "figure6_stage_specific_state_signature_heatmap_R.png"), fig6, width = 8.5, height = 5.4, dpi = 400)
ggsave(file.path(out_dir, "figure6_stage_specific_state_signature_heatmap_R.pdf"), fig6, width = 8.5, height = 5.4)
ggsave(file.path(out_dir, "figure6_stage_specific_state_signature_heatmap_R.svg"), fig6, width = 8.5, height = 5.4)

# Figure 7: D3 map
d3 <- read_csv(file.path(out_dir, "figure7_d3_endpoint_signature_map_source.csv"), show_col_types = FALSE)
fig7 <- ggplot(d3, aes(x = D3_EC_endpoint, y = D3_YAP_contractile)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_point(aes(fill = D3_ECM_mesenchymal, color = regulator_group), shape = 21, size = 3, stroke = 0.9, alpha = 0.9) +
  geom_text_repel(data = filter(d3, target_name %in% label_genes),
                  aes(label = target_name), size = 2.4, max.overlaps = Inf, seed = 1) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       name = "D3 ECM/mesenchymal") +
  scale_color_manual(values = group_cols, name = "Regulator class") +
  labs(title = "D3 endpoint map keeps YAP/ECM failure separate from D2 P19",
       x = "D3 EC endpoint effect\nmean(P2,P3,P8,P11)",
       y = "D3 P34 YAP/contractile effect") +
  theme_pub
ggsave(file.path(out_dir, "figure7_d3_endpoint_signature_map_R.png"), fig7, width = 6.2, height = 4.8, dpi = 400)
ggsave(file.path(out_dir, "figure7_d3_endpoint_signature_map_R.pdf"), fig7, width = 6.2, height = 4.8)
ggsave(file.path(out_dir, "figure7_d3_endpoint_signature_map_R.svg"), fig7, width = 6.2, height = 4.8)
