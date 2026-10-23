
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(tidyr)
})

data_dir <- "."
out_dir <- data_dir

reg_cols <- c(
  "TGFβ/Activin arm" = "#d55e00",
  "BMP/SMAD1/5 arm" = "#e69f00",
  "Shared TGF/BMP" = "#cc79a7",
  "Hippo/YAP" = "#6a3d9a",
  "EC specification" = "#0072b2",
  "Chromatin/cofactor" = "#009e73",
  "RNA processing" = "#56b4e9",
  "Other" = "#bdbdbd"
)

theme_pub <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

prog <- read_csv(file.path(data_dir, "global_d2_program_discovery.csv"), show_col_types = FALSE)
scores <- read_csv(file.path(data_dir, "target_state_scores_d2_d3.csv"), show_col_types = FALSE)
corrd3 <- read_csv(file.path(data_dir, "d3_program_correlations_with_d2_p19_vs_ec.csv"), show_col_types = FALSE)
selected <- read_csv(file.path(data_dir, "selected_regulator_state_scores.csv"), show_col_types = FALSE)

key_programs <- c(14,15,17,19,23,34,33,30,13,28)

fig1 <- ggplot(prog, aes(x = mean_score_D2, y = abs_burden_D2)) +
  geom_point(aes(size = n_sig_D2, color = corr_D2effect_D3_EC), alpha = 0.9) +
  geom_text_repel(data = filter(prog, program %in% key_programs),
                  aes(label = paste0("P", program)), size = 3, max.overlaps = Inf) +
  scale_color_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
                        name = "corr(D2 program,\nD3 EC)") +
  scale_size_continuous(name = "# significant\nregulators", range = c(2, 8)) +
  labs(title = "Global D2 program screen identifies P19", x = "Mean activity at D2", y = "D2 perturbation burden\nsum |significant effects|") +
  theme_pub
ggsave(file.path(out_dir, "R_figure1_global_d2_program_discovery.png"), fig1, width = 7.2, height = 5.0, dpi = 400)

fig3 <- ggplot(scores, aes(x = D2_EC_induction, y = D2_P19_vs_latmeso_outflow)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
  geom_point(aes(fill = D2_P19, color = regulator_class), shape = 21, size = 3, stroke = 1) +
  geom_text_repel(data = filter(scores, target_name %in% c("ACVR2A","ACVR2B","BMPR1A","SMAD1","SMAD4","SMAD5","ACVR1B","TGFBR2","SMAD3","ETV2","KDR","NF2","LATS2","AMOTL2")),
                  aes(label = target_name), size = 2.7, max.overlaps = Inf) +
  scale_fill_viridis_c(name = "D2 P19 effect") +
  scale_color_manual(values = reg_cols, name = "Regulator class") +
  labs(title = "D2 state map separates posterior P19 skew from EC-induction",
       x = "D2 EC-induction effect\nmean(P15, P23)",
       y = "D2 posterior-skew effect\nP19 - mean(P14, P17)") +
  theme_pub
ggsave(file.path(out_dir, "R_figure3_d2_state_map_p19_p14_ec.png"), fig3, width = 6.8, height = 5.2, dpi = 400)

topn <- scores %>% arrange(desc(D2_P19_vs_EC)) %>% slice_head(n = 24) %>% mutate(target_name = factor(target_name, levels = rev(target_name)))
fig4 <- ggplot(topn, aes(x = D2_P19_vs_EC, y = target_name, fill = regulator_class)) +
  geom_col() +
  geom_vline(xintercept = 0, color = "grey40") +
  scale_fill_manual(values = reg_cols) +
  labs(title = "Regulators ranked by D2 P19-high / EC-induction-low signature",
       x = "D2 P19-vs-EC score\nP19 - mean(P15, P23)",
       y = NULL, fill = "Class") +
  theme_pub
ggsave(file.path(out_dir, "R_figure4_d2_p19_vs_ec_failure_rank.png"), fig4, width = 6.0, height = 6.6, dpi = 400)

corrplot <- bind_rows(slice_head(corrd3, n = 12), slice_tail(corrd3, n = 12)) %>%
  arrange(corr_D2P19vsEC_to_D3_program) %>%
  mutate(program_label = factor(program_label, levels = program_label),
         direction = ifelse(corr_D2P19vsEC_to_D3_program < 0, "anti-correlated", "correlated"))
fig5 <- ggplot(corrplot, aes(x = corr_D2P19vsEC_to_D3_program, y = program_label, fill = direction)) +
  geom_col() +
  geom_vline(xintercept = 0, color = "grey40") +
  scale_fill_manual(values = c("anti-correlated" = "#0072B2", "correlated" = "#D55E00")) +
  labs(title = "What follows D2 P19 skew?",
       x = "Correlation with D2 P19-vs-EC score\nacross regulators", y = NULL, fill = NULL) +
  theme_pub
ggsave(file.path(out_dir, "R_figure5_d3_outcomes_correlated_with_d2_p19_skew.png"), fig5, width = 7.0, height = 6.2, dpi = 400)

heat_states <- c("D2_P19","D2_P14","D2_P17","D2_EC_induction","D2_P19_vs_EC","D3_EC_endpoint","D3_P33","D3_P24","D3_P34","D3_endpoint_score")
state_labels <- c("D2 P19","D2 P14","D2 P17","D2 EC induction","D2 P19-vs-EC","D3 EC endpoint","D3 P33 mesenchymal","D3 P24 HOX","D3 P34 YAP","D3 endpoint score")
heat <- selected %>%
  arrange(desc(D2_P19_vs_EC)) %>%
  select(target_name, all_of(heat_states)) %>%
  pivot_longer(-target_name, names_to = "state", values_to = "effect") %>%
  mutate(state = factor(state, levels = heat_states, labels = state_labels),
         target_name = factor(target_name, levels = unique(selected[order(-selected$D2_P19_vs_EC), "target_name"]$target_name)))
fig7 <- ggplot(heat, aes(x = state, y = target_name, fill = effect)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, name = "Effect/score") +
  labs(title = "Top D2 P19-skew regulators and downstream state effects", x = NULL, y = NULL) +
  theme_pub + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "R_figure7_selected_regulators_p19_downstream_heatmap.png"), fig7, width = 8.5, height = 6.8, dpi = 400)
