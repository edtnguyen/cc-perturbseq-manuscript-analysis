
suppressPackageStartupMessages({library(readr); library(dplyr); library(ggplot2); library(ggrepel); library(scales)})
data_dir <- "./tgfb_d2_celllevel_csv_exports"
out_dir <- file.path(data_dir, "ggplot2_figures"); dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
target_order <- c("non-targeting","ACVR1B","ACVR2A","ACVR2B","BMPR1A","SMAD1","SMAD3","SMAD4","SMAD5","TGFBR2")
target_cols <- c("non-targeting"="grey45", "ACVR1B"="#D55E00", "ACVR2A"="#D55E00", "ACVR2B"="#D55E00", "BMPR1A"="#D55E00", "SMAD1"="#D55E00", "SMAD3"="#D55E00", "SMAD4"="#D55E00", "SMAD5"="#D55E00", "TGFBR2"="#D55E00")
theme_pub <- theme_bw(base_size=11) + theme(panel.grid.minor=element_blank(), panel.grid.major=element_line(color="grey90", linewidth=0.3), axis.title=element_text(face="bold"), plot.title=element_text(face="bold", hjust=0.5), axis.text.x=element_text(angle=45,hjust=1), legend.title=element_text(face="bold"))

# Fig 1: D2 P19 failure score distribution
cell <- read_csv(file.path(data_dir,"d2_tgfb_cell_scores_by_target.csv"), show_col_types=FALSE) %>% mutate(target=factor(target, levels=target_order))
fig1 <- ggplot(cell, aes(x=target, y=d2_p19_failure_score, fill=target)) + geom_violin(trim=TRUE, scale="width", alpha=0.55, color="black", linewidth=0.25) + geom_boxplot(width=0.15, outlier.shape=NA, alpha=0.85, linewidth=0.25) + geom_hline(yintercept=0, linetype="dashed", linewidth=0.5) + scale_fill_manual(values=target_cols, guide="none") + labs(title="D2 P19 failure score by TGF/BMP-SMAD target carrier", x="Guide target carried by cell", y="D2 P19 failure score\nzP19 - mean(zP14,zP15,zP23)") + theme_pub
ggsave(file.path(out_dir,"fig1_d2_tgfb_p19_failure_score_by_target.png"), fig1, width=7.2, height=4.4, dpi=400)
ggsave(file.path(out_dir,"fig1_d2_tgfb_p19_failure_score_by_target.pdf"), fig1, width=7.2, height=4.4)
ggsave(file.path(out_dir,"fig1_d2_tgfb_p19_failure_score_by_target.svg"), fig1, width=7.2, height=4.4)

# Fig 2: program distributions
prog <- read_csv(file.path(data_dir,"d2_tgfb_program_scores_long_by_target.csv"), show_col_types=FALSE) %>% mutate(target=factor(target, levels=target_order), program=factor(program, levels=c("Program 14","Program 19","Program 17","Program 15","Program 23","Program 34")))
fig2 <- ggplot(prog, aes(x=target, y=z_score, fill=target)) + geom_violin(trim=TRUE, scale="width", alpha=0.55, color="black", linewidth=0.25) + geom_boxplot(width=0.14, outlier.shape=NA, alpha=0.85, linewidth=0.25) + geom_hline(yintercept=0, linetype="dashed", linewidth=0.4) + facet_wrap(~program, nrow=1, scales="free_y") + scale_fill_manual(values=target_cols, guide="none") + labs(title="D2 program-score shifts by TGF/BMP-SMAD target carrier", x="Guide target carried by cell", y="Program score, z-scored within D2") + theme_pub + theme(strip.text=element_text(face="bold"))
ggsave(file.path(out_dir,"fig2_d2_tgfb_program_scores_by_target.png"), fig2, width=14.0, height=4.0, dpi=400)
ggsave(file.path(out_dir,"fig2_d2_tgfb_program_scores_by_target.pdf"), fig2, width=14.0, height=4.0)
ggsave(file.path(out_dir,"fig2_d2_tgfb_program_scores_by_target.svg"), fig2, width=14.0, height=4.0)

# Fig 3: branch fractions
branch <- read_csv(file.path(data_dir,"d2_tgfb_binary_branch_fraction_by_target.csv"), show_col_types=FALSE) %>% mutate(target=factor(target, levels=target_order), binary_branch=factor(binary_branch, levels=c("EC/lateral-mesoderm biased","P19-biased / failed-transition")))
fig3 <- ggplot(branch, aes(x=target, y=fraction, fill=binary_branch)) + geom_col(width=0.7, color="black", linewidth=0.25) + scale_y_continuous(labels=percent_format(accuracy=1), expand=c(0,0)) + scale_fill_manual(values=c("EC/lateral-mesoderm biased"="#0072B2", "P19-biased / failed-transition"="#D55E00")) + labs(title="D2 P19-biased fraction by target carrier", x="Guide target carried by cell", y="Fraction of D2 carrier cells", fill="D2 branch bias") + theme_pub
ggsave(file.path(out_dir,"fig3_d2_tgfb_p19_branch_fraction.png"), fig3, width=7.4, height=4.5, dpi=400)
ggsave(file.path(out_dir,"fig3_d2_tgfb_p19_branch_fraction.pdf"), fig3, width=7.4, height=4.5)
ggsave(file.path(out_dir,"fig3_d2_tgfb_p19_branch_fraction.svg"), fig3, width=7.4, height=4.5)

# Fig 4: guide-level concordance
guide <- read_csv(file.path(data_dir,"d2_tgfb_guide_level_concordance.csv"), show_col_types=FALSE) %>% mutate(target=factor(target, levels=target_order))
fig4 <- ggplot(guide, aes(x=target, y=d2_p19_failure_score_mean, color=target)) + geom_hline(yintercept=0, linetype="dashed", linewidth=0.5) + geom_jitter(aes(size=n_cells), width=0.15, height=0, alpha=0.85) + stat_summary(fun=mean, geom="crossbar", width=0.45, color="black", linewidth=0.35) + scale_color_manual(values=target_cols, guide="none") + scale_size_continuous(range=c(1.5,5.0)) + labs(title="Guide-level concordance of D2 P19 failure score", x="Guide target", y="Mean D2 P19 failure score per guide", size="# D2 carrier cells") + theme_pub
ggsave(file.path(out_dir,"fig4_d2_tgfb_guide_level_concordance.png"), fig4, width=7.2, height=4.5, dpi=400)
ggsave(file.path(out_dir,"fig4_d2_tgfb_guide_level_concordance.pdf"), fig4, width=7.2, height=4.5)
ggsave(file.path(out_dir,"fig4_d2_tgfb_guide_level_concordance.svg"), fig4, width=7.2, height=4.5)

# Fig 5: lead-lag target-level
lead <- read_csv(file.path(data_dir,"d2_tgfb_p19_failure_vs_d3_endpoint_leadlag.csv"), show_col_types=FALSE) %>% mutate(target=factor(target, levels=target_order))
fig5 <- ggplot(lead, aes(x=d2_p19_failure_score_mean, y=d3_endpoint_score_mean, label=target, color=target)) + geom_hline(yintercept=0, color="grey75", linewidth=0.4) + geom_vline(xintercept=0, color="grey75", linewidth=0.4) + geom_point(size=3.0) + geom_text_repel(size=3.1, show.legend=FALSE, max.overlaps=Inf) + scale_color_manual(values=target_cols, guide="none") + labs(title="D2 P19 failure bias versus D3 EC endpoint score", x="Mean D2 P19 failure score", y="Mean D3 EC endpoint score") + theme_pub
ggsave(file.path(out_dir,"fig5_tgfb_d2_p19_failure_vs_d3_endpoint.png"), fig5, width=5.6, height=4.5, dpi=400)
ggsave(file.path(out_dir,"fig5_tgfb_d2_p19_failure_vs_d3_endpoint.pdf"), fig5, width=5.6, height=4.5)
ggsave(file.path(out_dir,"fig5_tgfb_d2_p19_failure_vs_d3_endpoint.svg"), fig5, width=5.6, height=4.5)
