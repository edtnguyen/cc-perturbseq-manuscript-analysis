## 2025/1/16 Plotting Genome-wide screen analysis
## Dulguun
## Replace "/Users/dulguun/oak" with "/oak/stanford/groups/engreitz"

if(!"MAGeCKFlute" %in% installed.packages()) BiocManager::install("MAGeCKFlute")
if(!"clusterProfiler" %in% installed.packages()) BiocManager::install("clusterProfiler")
if(!"ggplot2" %in% installed.packages()) BiocManager::install("ggplot2")

library(MAGeCKFlute)
library(clusterProfiler)
library(ggplot2)

#install.packages("readxl")
#install.packages("plotly")
library(readxl)

setwd("/Users/dulguun/oak/Users/dulguun/230916_Genome-wide_Screen_Analysis")
OUTDIR = "/Users/dulguun/oak/Users/dulguun/230916_Genome-wide_Screen_Analysis/figures_250115"

####################
### VOLCANO PLOT ###
####################

## Plotting overall Mageck results
data_m = "/Users/dulguun/oak/Users/dulguun/230916_Genome-wide_Screen_Analysis/mageck/Pos_vs_Neg/Pos_vs_Neg.gene_summary_GeneNameCorrected.txt"
gdata = ReadRRA(data_m)
head(gdata)

p1 = VolcanoView(gdata, x = "Score", y = "FDR", Label = "id", x_cutoff = 0, ylab = "-Log10(FDR)")
p1
ggsave(paste0(OUTDIR,"/D3pos_vs_D3neg_mageck_overall_volcano_plot.pdf"), plot = p1, width = 8, height = 6, dpi = 300)

## Plotting U-test results
data_u <- read.csv('full_result_files_from_Chad/Pos vs Neg/combined_utest_P1P2_adj_test.txt', sep = '')
data_u$hits <- data_u$p_value_corrected < 0.05 # Add a column to indicate if FDR < 0.05

plot_u <- ggplot(data_u, aes(x = log2FC, y = -log10(p_value_corrected))) +
  geom_point(alpha = 0.4, size = 1.75, aes(color = hits)) +
  geom_hline(yintercept = fdr_threshold, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("black", "red")) +
  labs(x = "Log2FC", y = "-Log10(FDR)") +
  theme_bw(base_size = 16) +
  ggtitle("D3 pos vs D3 neg (U-test)") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(paste0(OUTDIR, "/D3pos_vs_D3neg_u_test_volcano_plot.pdf"), plot = plot_u, width = 8, height = 6, dpi = 300)

## Probably not super necessary plots
## Plotting negative selection Mageck only
data <- read.csv('full_result_files_from_Chad/Pos vs Neg/Pos_vs_Neg.gene_summary.txt', sep = '')
data$neg.hits <- data$neg.fdr < 0.05 # Add a column to indicate if FDR < 0.05
data$pos.hits <- data$pos.fdr < 0.05 # Add a column to indicate if FDR < 0.05
fdr_threshold <- -log10(0.05)

plot_neg <- ggplot(data, aes(x = neg.lfc, y = -log10(neg.fdr))) +
  geom_point(alpha = 0.4, size = 1.75, aes(color = neg.hits)) +
  geom_hline(yintercept = fdr_threshold, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("black", "red")) +
  labs(x = "Log2FC", y = "-Log10(FDR)") +
  theme_bw(base_size = 16) +
  ggtitle("D3 pos vs D3 neg (Mageck negative selection)") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("figures_250115/D3pos_vs_D3neg_mageck_negative_selection_volcano_plot.pdf", plot = plot_neg, width = 8, height = 6, dpi = 300)

## Plotting positive selection Mageck only
plot_pos <- ggplot(data, aes(x = pos.lfc, y = -log10(pos.fdr))) +
  geom_point(alpha = 0.4, size = 1.75, aes(color = pos.hits)) +
  geom_hline(yintercept = fdr_threshold, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("black", "red")) +
  labs(x = "Log2FC", y = "-Log10(FDR)") +
  theme_bw(base_size = 16) +
  ggtitle("D3 pos vs D3 neg (Mageck positive selection)") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("figures_250115/D3pos_vs_D3neg_mageck_positive_selection_volcano_plot.pdf", plot = plot_pos, width = 8, height = 6, dpi = 300)
