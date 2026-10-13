
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
})

data_dir <- "./tgfb_bmp_d2_celllevel_exports"
out_dir <- file.path(data_dir, "ggplot2_figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

target_order <- c("non-targeting", "ACVR1B", "TGFBR2", "SMAD3", "BMPR1A", "SMAD1", "SMAD5", "ACVR2A", "ACVR2B", "SMAD4")

target_cols <- c(
  "non-targeting" = "grey45",
  "ACVR1B" = "#7B3294",
  "TGFBR2" = "#7B3294",
  "SMAD3" = "#7B3294",
  "BMPR1A" = "#008837",
  "SMAD1" = "#008837",
  "SMAD5" = "#008837",
  "ACVR2A" = "#D95F02",
  "ACVR2B" = "#D95F02",
  "SMAD4" = "#D95F02"
)

branch_cols <- c(
  "EC-induction-biased" = "#009E73",
  "P19-high / EC-induction-low" = "#CC79A7"
)

program_cols <- c(
  "Program 14" = "#E69F00",
  "Program 17" = "#D55E00",
  "Program 19" = "#CC79A7",
  "Program 15" = "#009E73",
  "Program 23" = "#0072B2",
  "Program 34" = "#6A3D9A"
)

theme_pub <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )

# D2 scores by target
d2 <- read_csv(file.path(data_dir, "d2_cell_scores_by_target.csv"), show_col_types = FALSE) %>%
  mutate(target = factor(target, levels = target_order))

fig1 <- ggplot(d2, aes(x = target, y = D2_P19_vs_EC_score, fill = target)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.55, color = "black", linewidth = 0.25) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.85, linewidth = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  scale_fill_manual(values = target_cols, guide = "none") +
  labs(
    title = "D2 P19-vs-EC induction score by TGF/BMP target-carrier group",
    x = "Guide target carried by cell",
    y = "D2 P19-vs-EC score\nzP19 - mean(zP15,zP23)"
  ) +
  theme_pub
ggsave(file.path(out_dir, "celllevel_fig1_d2_p19_vs_ec_score.png"), fig1, width = 7.0, height = 4.4, dpi = 400)
ggsave(file.path(out_dir, "celllevel_fig1_d2_p19_vs_ec_score.pdf"), fig1, width = 7.0, height = 4.4)
ggsave(file.path(out_dir, "celllevel_fig1_d2_p19_vs_ec_score.svg"), fig1, width = 7.0, height = 4.4)

# Program score distributions
prog <- read_csv(file.path(data_dir, "d2_program_scores_long_by_target.csv"), show_col_types = FALSE) %>%
  mutate(target = factor(target, levels = target_order),
         program = factor(program, levels = c("Program 19","Program 14","Program 17","Program 15","Program 23","Program 34")))

fig2 <- ggplot(prog, aes(x = target, y = z_score, fill = target)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.55, color = "black", linewidth = 0.25) +
  geom_boxplot(width = 0.14, outlier.shape = NA, alpha = 0.85, linewidth = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  facet_wrap(~program, nrow = 2, scales = "free_y") +
  scale_fill_manual(values = target_cols, guide = "none") +
  labs(
    title = "D2 program-score shifts by TGF/BMP target-carrier group",
    x = "Guide target carried by cell",
    y = "Program score, z-scored within D2"
  ) +
  theme_pub +
  theme(strip.background = element_rect(fill = "grey95", color = "grey80"),
        strip.text = element_text(face = "bold", size = 10))
ggsave(file.path(out_dir, "celllevel_fig2_d2_program_scores_by_target.png"), fig2, width = 12.0, height = 6.0, dpi = 400)
ggsave(file.path(out_dir, "celllevel_fig2_d2_program_scores_by_target.pdf"), fig2, width = 12.0, height = 6.0)
ggsave(file.path(out_dir, "celllevel_fig2_d2_program_scores_by_target.svg"), fig2, width = 12.0, height = 6.0)

# Branch fractions
branch <- read_csv(file.path(data_dir, "d2_p19_vs_ec_branch_fraction_by_target.csv"), show_col_types = FALSE) %>%
  mutate(target = factor(target, levels = target_order),
         D2_binary_branch = factor(D2_binary_branch, levels = c("EC-induction-biased","P19-high / EC-induction-low")))
fig3 <- ggplot(branch, aes(x = target, y = fraction, fill = D2_binary_branch)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.25) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
  scale_fill_manual(values = branch_cols) +
  labs(
    title = "D2 P19-vs-EC branch fractions by target-carrier group",
    x = "Guide target carried by cell",
    y = "Fraction of D2 carrier cells",
    fill = "D2 branch bias"
  ) +
  theme_pub
ggsave(file.path(out_dir, "celllevel_fig3_d2_p19_vs_ec_branch_fraction.png"), fig3, width = 7.0, height = 4.6, dpi = 400)
ggsave(file.path(out_dir, "celllevel_fig3_d2_p19_vs_ec_branch_fraction.pdf"), fig3, width = 7.0, height = 4.6)
ggsave(file.path(out_dir, "celllevel_fig3_d2_p19_vs_ec_branch_fraction.svg"), fig3, width = 7.0, height = 4.6)

# Guide-level concordance
guide <- read_csv(file.path(data_dir, "d2_guide_level_concordance.csv"), show_col_types = FALSE) %>%
  mutate(target = factor(target, levels = target_order))
fig4 <- ggplot(guide, aes(x = target, y = D2_P19_vs_EC_score_mean, color = target)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_jitter(aes(size = n_cells_d2), width = 0.15, height = 0, alpha = 0.85) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.45, color = "black", linewidth = 0.35) +
  scale_color_manual(values = target_cols, guide = "none") +
  scale_size_continuous(range = c(1.5, 5.0)) +
  labs(
    title = "Guide-level concordance of D2 P19-vs-EC score",
    x = "Guide target",
    y = "Mean D2 P19-vs-EC score per guide",
    size = "# D2 carrier cells"
  ) +
  theme_pub
ggsave(file.path(out_dir, "celllevel_fig4_guide_level_concordance_p19_vs_ec.png"), fig4, width = 7.0, height = 4.6, dpi = 400)
ggsave(file.path(out_dir, "celllevel_fig4_guide_level_concordance_p19_vs_ec.pdf"), fig4, width = 7.0, height = 4.6)
ggsave(file.path(out_dir, "celllevel_fig4_guide_level_concordance_p19_vs_ec.svg"), fig4, width = 7.0, height = 4.6)

# Lead-lag target-level
lead <- read_csv(file.path(data_dir, "d2_tgfb_scores_vs_d3_endpoint_leadlag.csv"), show_col_types = FALSE) %>%
  mutate(target = factor(target, levels = target_order))
fig5 <- ggplot(lead, aes(x = D2_P19_vs_EC_score_mean, y = D3_endpoint_score_mean, label = target, color = target)) +
  geom_hline(yintercept = 0, color = "grey75", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "grey75", linewidth = 0.4) +
  geom_point(size = 3.0) +
  geom_text_repel(size = 3.1, max.overlaps = Inf, show.legend = FALSE) +
  scale_color_manual(values = target_cols, guide = "none") +
  labs(
    title = "D2 P19-vs-EC bias versus D3 EC endpoint score",
    x = "Mean D2 P19-vs-EC score",
    y = "Mean D3 EC endpoint score"
  ) +
  theme_pub
ggsave(file.path(out_dir, "celllevel_fig5_d2_p19_vs_ec_to_d3_endpoint.png"), fig5, width = 5.8, height = 4.6, dpi = 400)
ggsave(file.path(out_dir, "celllevel_fig5_d2_p19_vs_ec_to_d3_endpoint.pdf"), fig5, width = 5.8, height = 4.6)
ggsave(file.path(out_dir, "celllevel_fig5_d2_p19_vs_ec_to_d3_endpoint.svg"), fig5, width = 5.8, height = 4.6)

message("Done. Figures saved to: ", out_dir)
