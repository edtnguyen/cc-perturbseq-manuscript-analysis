
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2); library(ggrepel); library(tidyr)
})

data_dir <- "/mnt/data/tgfb_d2_signature_updates"
out_dir <- data_dir
class_cols <- c(
  "TGF/BMP/SMAD" = "#D55E00", "Hippo/YAP" = "#6A3D9A", "EC specification" = "#0072B2",
  "Chromatin/cofactor" = "#1B9E77", "RNA processing/translation" = "#E7298A",
  "Mesoderm/developmental" = "#A6761D", "Apoptosis/stress" = "#E6AB02",
  "Trafficking/metabolism" = "#66A61E", "Other" = "#BDBDBD")
theme_pub <- theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "grey90", linewidth = 0.3), axis.title = element_text(face = "bold"), axis.text = element_text(color = "black"), plot.title = element_text(face = "bold", hjust = 0.5, size = 11), legend.title = element_text(face = "bold"), legend.key = element_blank())

d2 <- read_csv(file.path(data_dir, "D2_state_scores.csv"), show_col_types = FALSE)
lead <- read_csv(file.path(data_dir, "D2_to_D3_scores.csv"), show_col_types = FALSE)
heat <- read_csv(file.path(data_dir, "D2_selected_program_effects.csv"), show_col_types = FALSE)
state <- read_csv(file.path(data_dir, "state_signature_heatmap.csv"), show_col_types = FALSE)

# This R script recreates the same figures rendered in Python. Run with: Rscript make_tgfb_d2_program_signature_figures_ggplot2.R
# Figure 1
label_genes <- c("ACVR2A","BMPR1A","ACVR2B","SMAD1","SMAD4","SMAD5","TGFBR2","SETDB1","KDR")
fig1 <- ggplot(d2, aes(x=D2_EC_induction, y=D2_P19_minus_P14, fill=P19_posterior_mesoderm, color=regulator_class)) + geom_hline(yintercept=0, color="grey60", linewidth=0.4) + geom_vline(xintercept=0, color="grey60", linewidth=0.4) + geom_point(shape=21, size=2.7, stroke=0.9, alpha=0.9) + geom_text_repel(data=d2 %>% filter(target_name %in% label_genes), aes(label=target_name), size=2.7, show.legend=FALSE, max.overlaps=Inf, seed=1) + scale_color_manual(values=class_cols, drop=FALSE) + scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B", midpoint=0, limits=c(-4,4), oob=scales::squish) + labs(title="D2 TGF/BMP-SMAD perturbations shift toward P19 and away from P14/EC induction", x="D2 EC-induction effect\nmean(P15, P23)", y="D2 P19/P14 imbalance\nP19 - P14", fill="P19 effect\n(D2)", color="Regulator class") + theme_pub

ggsave(file.path(out_dir, "figure_tgfb_d2_state_map_clean_R.png"), fig1, width=7.0, height=4.7, dpi=400)
ggsave(file.path(out_dir, "figure_tgfb_d2_state_map_clean_R.pdf"), fig1, width=7.0, height=4.7)
ggsave(file.path(out_dir, "figure_tgfb_d2_state_map_clean_R.svg"), fig1, width=7.0, height=4.7)

# Figure 2
fig2 <- ggplot(heat, aes(x=program_label, y=target_name, fill=effect, size=pmin(abs(effect),4))) + geom_point(shape=21, color="grey30", stroke=0.25) + scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B", midpoint=0, limits=c(-5,5), oob=scales::squish) + scale_size_continuous(range=c(1.3,5.8), breaks=c(1,2,4), name="|effect|") + labs(title="D2 program effects for TGF/BMP-SMAD and comparison regulators", x=NULL, y=NULL, fill="Significant effect\n(non-significant = 0)") + theme_pub + theme(axis.text.x=element_text(angle=35, hjust=1), axis.text.y=element_text(face="bold"))
ggsave(file.path(out_dir, "figure_tgfb_d2_selected_program_heatmap_R.png"), fig2, width=9.2, height=5.8, dpi=400)
ggsave(file.path(out_dir, "figure_tgfb_d2_selected_program_heatmap_R.pdf"), fig2, width=9.2, height=5.8)
ggsave(file.path(out_dir, "figure_tgfb_d2_selected_program_heatmap_R.svg"), fig2, width=9.2, height=5.8)

# Figure 3
rank_d2 <- d2 %>% arrange(desc(D2_P19_failure_score)) %>% slice_head(n=22) %>% mutate(target_name=factor(target_name, levels=rev(target_name)))
fig3 <- ggplot(rank_d2, aes(x=D2_P19_failure_score, y=target_name, fill=regulator_class)) + geom_col(color="black", linewidth=0.25, width=0.75) + scale_fill_manual(values=class_cols, drop=FALSE) + labs(title="TGF/BMP-SMAD perturbations are the strongest D2 P19-shift regulators", x="D2 P19 failure score\nP19 - mean(P14, P15, P23)", y=NULL, fill="Class") + theme_pub
ggsave(file.path(out_dir, "figure_tgfb_d2_p19_failure_rank_R.png"), fig3, width=7.2, height=5.9, dpi=400)
ggsave(file.path(out_dir, "figure_tgfb_d2_p19_failure_rank_R.pdf"), fig3, width=7.2, height=5.9)
ggsave(file.path(out_dir, "figure_tgfb_d2_p19_failure_rank_R.svg"), fig3, width=7.2, height=5.9)

# Figure 4
corr <- cor(lead$D2_P19_failure_score, lead$D3_endpoint_score, use="pairwise.complete.obs")
lead_labels <- c("ACVR2A","ACVR2B","BMPR1A","SMAD1","SMAD4","SMAD5","TGFBR2","NF2","LATS2","ETV2","KDR","SETDB1")
fig4 <- ggplot(lead, aes(x=D2_P19_failure_score, y=D3_endpoint_score, color=regulator_class)) + geom_hline(yintercept=0, color="grey60", linewidth=0.4) + geom_vline(xintercept=0, color="grey60", linewidth=0.4) + geom_point(alpha=0.85, size=2.0) + geom_smooth(method="lm", se=FALSE, color="black", linetype="dashed", linewidth=0.6) + geom_text_repel(data=lead %>% filter(target_name %in% lead_labels), aes(label=target_name), size=2.7, show.legend=FALSE, max.overlaps=Inf, seed=2) + annotate("text", x=-Inf, y=Inf, hjust=-0.05, vjust=1.2, label=paste0("Pearson r = ", sprintf("%.2f", corr)), size=3.0) + scale_color_manual(values=class_cols, drop=FALSE) + labs(title="D2 P19-shift signature is associated with later EC endpoint loss", x="D2 P19 failure score\nP19 - mean(P14, P15, P23)", y="D3 EC endpoint score\nmean(P2,P3,P8,P11) - mean(P27,P34,P35)", color="Class") + theme_pub
ggsave(file.path(out_dir, "figure_tgfb_d2_failure_vs_d3_endpoint_R.png"), fig4, width=7.0, height=5.2, dpi=400)
ggsave(file.path(out_dir, "figure_tgfb_d2_failure_vs_d3_endpoint_R.pdf"), fig4, width=7.0, height=5.2)
ggsave(file.path(out_dir, "figure_tgfb_d2_failure_vs_d3_endpoint_R.svg"), fig4, width=7.0, height=5.2)

# State heatmap and D3 map versions can be recreated analogously from state_signature_heatmap.csv and D2_to_D3_scores.csv.
