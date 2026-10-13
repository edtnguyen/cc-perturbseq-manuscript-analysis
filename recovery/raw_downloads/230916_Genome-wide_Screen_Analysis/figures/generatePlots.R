library(ggplot2)

setwd("/Users/dulguun/oak/Users/dulguun/230916_Genome-wide_Screen_Analysis/figures")
fileName = file.path("Pos vs Neg", "Pos_vs_Neg.gene_summary.txt")
data <- read.csv(fileName, sep = '\t')

## VOLCANO PLOT ##
data$log10_fdr <- -log10(data$neg.fdr)
data$log10_p <- -log10(data$neg.p.value)
data$log10_pos_fdr <- -log10(data$pos.fdr)
fdr_threshold <- -log10(0.05)

plot <- ggplot(data, aes(x = neg.lfc, y = log10_fdr)) +
  geom_point(alpha = 0.4, size = 1.75) +
  labs(x = "Log2FC", y = "-Log10(FDR)") +
  ggtitle("D3 vs D0") +
  theme(plot.title = element_text(hjust = 0.5))

plot2 <- ggplot(data, aes(x = neg.lfc, y = log10_p)) +
  geom_point(alpha = 0.4, size = 1.75) +
  labs(x = "Log2FC", y = "-Log10(P)") +
  ggtitle("D3 vs D0") +
  theme(plot.title = element_text(hjust = 0.5))

plot3 <- ggplot(data, aes(x = pos.lfc, y = log10_pos_fdr)) +
  geom_point(alpha = 0.4, size = 1.75) +
  geom_hline(yintercept = fdr_threshold, linetype = "dashed", color = "red") +
  labs(x = "Log2FC", y = "-Log10(FDR)") +
  ggtitle("D3 vs D0") +
  theme(plot.title = element_text(hjust = 0.5))
