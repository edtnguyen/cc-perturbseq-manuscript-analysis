## Replace "/Users/dulguun/oak" with "/oak/stanford/groups/engreitz"

library(ggplot2)
library(readxl)
#install.packages("readxl")
#install.packages("plotly")


setwd("/Users/dulguun/oak/Users/dulguun/231023_Validation_Screen_Files")
fileName = "d3pos_d3neg_mageckLFC"
data <- read.csv(paste(fileName, '.csv',sep = ''))

####################
### VOLCANO PLOT ###
####################

data$hits <- data$FDR < 0.05 # Add a column to indicate if FDR < 0.05
fdr_threshold <- -log10(0.05)
highlight_gene <- "CDH5"
highlight_coords <- subset(data, Gene == highlight_gene)


plot_fdr <- ggplot(data, aes(x = Average.U.Test.FC, y = -log10(FDR.corrected))) +
  geom_point(alpha = 0.4, size = 1.75, aes(color = hits)) +
  geom_point(data = highlight_coords, color = "purple", alpha = 0.4, size = 1.75) +
  geom_hline(yintercept = fdr_threshold, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_text(aes(label = ifelse(Gene == highlight_gene, Gene, "")),
            vjust = -1, hjust = 0.5, color = "red") +
  scale_color_manual(values = c("black", "red")) +
  labs(x = "Log2FC", y = "-Log10(FDR)") +
  theme_bw(base_size = 16) +
  ggtitle("D3 pos vs D3 neg") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(paste(fileName,"_volcanoplot.pdf",sep = ""), plot = plot_fdr, width = 8, height = 6, dpi = 300)

#####################
### BR1 vs BR2 FC ###
#####################

data$U.Test.FC.BR1 <- as.numeric(data$U.Test.FC.BR1)
data$U.Test.FC.BR2 <- as.numeric(data$U.Test.FC.BR2)

temp <- na.omit(data[, c("U.Test.FC.BR1", "U.Test.FC.BR2")])

lm_fit <- lm(U.Test.FC.BR2 ~ U.Test.FC.BR1, data = temp)
r_squared <- summary(lm_fit)$r.squared

plot <- ggplot(data, aes(x = U.Test.FC.BR1, y = U.Test.FC.BR2)) +
  geom_point(alpha = 0.4, size = 1.75) +
  theme_bw(base_size = 16) +
  labs(x = "BR1 log2FC", y = "BR2 Log2FC") +
  ggtitle("Bio-replicate log2FC Comparison") +
  annotate("text", x = -Inf, y = Inf, label = paste("R²:", round(r_squared, 2)),
           hjust = -0.2, vjust = 2, size = 5, color = "red") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(paste(fileName,"_BioReplicate_Comparison.pdf",sep = ""), plot = plot, width = 6, height = 6, dpi = 300)

#################################
### Original vs Validation FC ###
#################################

data$Original.FC.U.Test <- as.numeric(data$Original.FC.U.Test)
data$Original.FC.Mageck <- as.numeric(data$Original.FC.Mageck)
data$Average.U.Test.FC <- as.numeric(data$Average.U.Test.FC)

# U-test
temp <- na.omit(data[, c("Original.FC.U.Test", "Average.U.Test.FC")])
#correlation <- cor(temp$Original.FC, temp$Average.U.Test.FC)

lm_fit <- lm(Average.U.Test.FC ~ Original.FC.U.Test, data = temp)
r_squared <- summary(lm_fit)$r.squared # Extract R-squared value

plot <- ggplot(data, aes(x = Original.FC.U.Test, y = Average.U.Test.FC)) +
  geom_point(alpha = 0.4, size = 1.75, aes(color = hits)) +
  scale_color_manual(values = c("black", "red")) +
  theme_bw(base_size = 16) +
  labs(x = "Original Screen log2FC", y = "Validation Screen Log2FC") +
  ggtitle("D3 pos vs D3 neg Comparison") +
  annotate("text", x = -Inf, y = Inf, label = paste("R²:", round(r_squared, 2)),
           hjust = -0.2, vjust = 2, size = 5, color = "red") +
  theme(plot.title = element_text(hjust = 0.5))

# Mageck
temp2 <- na.omit(data[, c("Original.FC.Mageck", "Average.U.Test.FC")])
lm_fit2 <- lm(Average.U.Test.FC ~ Original.FC.Mageck, data = temp2)
r_squared2 <- summary(lm_fit2)$r.squared # Extract R-squared value

plot2 <- ggplot(data, aes(x = Original.FC.Mageck, y = Average.U.Test.FC)) +
  geom_point(alpha = 0.4, size = 1.75, aes(color = hits)) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"),
                     labels = c("TRUE" = "Hit in Validation Screen", "FALSE" = "Non-Hit in Validation Screen")) +
  theme_bw(base_size = 16) +
  labs(x = "Original Screen log2FC", y = "Validation Screen Log2FC") +
  ggtitle("D3 pos vs D3 neg Comparison") +
  annotate("text", x = -Inf, y = Inf, label = paste("R²:", round(r_squared2, 2)),
           hjust = -0.2, vjust = 2, size = 5, color = "red") +
  scale_x_continuous(breaks = c(-4, -3, -2, -1, 0, 1, 2)) +  # Set specific breaks
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") + # Adding a diagonal line
  theme(plot.title = element_text(hjust = 0.5))

ggsave(paste(fileName,"_Original_Comparison_2.pdf",sep = ""), plot = plot2, width = 10, height = 6, dpi = 300)

# Highlight original hits
gene_type <- read_excel("validation_genes_type.xlsx")
colnames(gene_type)[colnames(gene_type) == "Target"] <- "Gene" # Rename the column 'Target' to 'Gene'
data_merged <- merge(data, gene_type, by = "Gene", all.x = TRUE)

plot <- ggplot(data_merged, aes(x = Original.FC.Mageck, y = Average.U.Test.FC)) +
  geom_point(alpha = 0.4, size = 1.75, aes(color = Type)) +
  #geom_point(aes(color = ifelse(type == "typeA", "typeA", "other"))) + # Highlight typeA
  #scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"),
                     #labels = c("TRUE" = "Hit in Validation Screen", "FALSE" = "Non-Hit in Validation Screen")) +
  scale_color_manual(values = c("genome_wide_hit" = "blue", "other" = "black")) +
  theme_bw(base_size = 16) +
  labs(x = "Original Screen log2FC", y = "Validation Screen Log2FC") +
  ggtitle("D3 pos vs D3 neg Comparison") +
  annotate("text", x = -Inf, y = Inf, label = paste("R²:", round(r_squared2, 2)),
           hjust = -0.2, vjust = 2, size = 5, color = "red") +
  scale_x_continuous(breaks = c(-4, -3, -2, -1, 0, 1, 2)) +  # Set specific breaks
  #geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") + # Adding a diagonal line
  theme(plot.title = element_text(hjust = 0.5))
plot

# Highlight certain genes
highlight_gene <- "TAOK1"
highlight_coords <- subset(data_merged, Gene == highlight_gene)

plot <- ggplot(data_merged, aes(x = Original.FC.Mageck, y = Average.U.Test.FC)) +
  geom_point(alpha = 0.4, size = 1.75, aes(color = Type)) +
  #geom_point(data = highlight_coords, color = "red", alpha = 0.4, size = 1.75) +
  #geom_text(aes(label = ifelse(Gene == highlight_gene, Gene, "")),
            #vjust = -1, hjust = 0.5, color = "red") +
  scale_color_manual(values = c("genome_wide_hit" = "blue", "other" = "black")) +
  theme_bw(base_size = 16) +
  labs(x = "Original Screen log2FC", y = "Validation Screen Log2FC") +
  ggtitle("D3 pos vs D3 neg Comparison") +
  annotate("text", x = -Inf, y = Inf, label = paste("R²:", round(r_squared2, 2)), hjust = -0.2, vjust = 2, size = 5, color = "red") +
  scale_x_continuous(breaks = c(-4, -3, -2, -1, 0, 1, 2)) +  # Set specific breaks
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") + # Adding a diagonal line
  theme(plot.title = element_text(hjust = 0.5))
plot
ggsave(paste(fileName,"_Original_Comparison_GW_hits.pdf",sep = ""), plot = plot, width = 9, height = 6, dpi = 300)

highlight_genes <- c("TAOK1", "TCEB3")
highlight_coords <- subset(data_merged, Gene %in% highlight_genes)
plot <- ggplot(data_merged, aes(x = Original.FC.Mageck, y = Average.U.Test.FC)) +
  geom_point(alpha = 0.4, size = 1.75, aes(color = Type)) +
  geom_point(data = highlight_coords, color = "red", alpha = 0.4, size = 1.75) +
  geom_text(aes(label = ifelse(Gene %in% highlight_genes, Gene, "")), vjust = -1, hjust = 0.5, color = "red") + # Label highlighted genes
  scale_color_manual(values = c("genome_wide_hit" = "blue", "other" = "black")) +
  theme_bw(base_size = 16) +
  labs(x = "Original Screen log2FC", y = "Validation Screen Log2FC") +
  ggtitle("D3 pos vs D3 neg Comparison") +
  annotate("text", x = -Inf, y = Inf, label = paste("R²:", round(r_squared2, 2)), hjust = -0.2, vjust = 2, size = 5, color = "red") +
  scale_x_continuous(breaks = c(-4, -3, -2, -1, 0, 1, 2)) +  # Set specific breaks
  #geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") + # Adding a diagonal line
  theme(plot.title = element_text(hjust = 0.5))
plot
ggsave(paste(fileName,"_Original_Comparison_TAOK1_TCEB3.pdf",sep = ""), plot = plot, width = 10, height = 6, dpi = 300)

highlight_gene <- "SKA3"
highlight_coords <- subset(data_merged, Gene == highlight_gene)
plot <- ggplot(data_merged, aes(x = Original.FC.Mageck, y = Average.U.Test.FC)) +
  geom_point(alpha = 0.4, size = 1.75, aes(color = Type)) +
  #geom_point(data = highlight_coords, color = "red", alpha = 0.4, size = 1.75) +
  geom_text(aes(label = ifelse(Gene == highlight_gene, Gene, "")), vjust = -1, hjust = 0.5, color = "red") + # Label highlighted genes
  scale_color_manual(values = c("genome_wide_hit" = "blue", "other" = "black")) +
  theme_bw(base_size = 16) +
  labs(x = "Original Screen log2FC", y = "Validation Screen Log2FC") +
  ggtitle("D3 pos vs D3 neg Comparison") +
  annotate("text", x = -Inf, y = Inf, label = paste("R²:", round(r_squared2, 2)), hjust = -0.2, vjust = 2, size = 5, color = "red") +
  scale_x_continuous(breaks = c(-4, -3, -2, -1, 0, 1, 2)) +  # Set specific breaks
  #geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") + # Adding a diagonal line
  theme(plot.title = element_text(hjust = 0.5))
plot
ggsave(paste(fileName,"_Original_Comparison_SKA3.pdf",sep = ""), plot = plot, width = 9, height = 6, dpi = 300)

highlight_genes <- c("WDR12", "CENPC", "NUP93", "SMC2", "RPS15A", "SKA3")
highlight_coords <- subset(data_merged, Gene %in% highlight_genes)
plot <- ggplot(data_merged, aes(x = Original.FC.Mageck, y = Average.U.Test.FC)) +
  geom_point(alpha = 0.4, size = 1.75, aes(color = Type)) +
  #geom_point(data = highlight_coords, color = "red", alpha = 0.4, size = 1.75) +
  geom_text(aes(label = ifelse(Gene %in% highlight_genes, Gene, "")), vjust = -1, hjust = 0.5, color = "red") + # Label highlighted genes
  scale_color_manual(values = c("genome_wide_hit" = "blue", "other" = "black")) +
  theme_bw(base_size = 16) +
  labs(x = "Original Screen log2FC", y = "Validation Screen Log2FC") +
  ggtitle("D3 pos vs D3 neg Comparison") +
  annotate("text", x = -Inf, y = Inf, label = paste("R²:", round(r_squared2, 2)), hjust = -0.2, vjust = 2, size = 5, color = "red") +
  scale_x_continuous(breaks = c(-4, -3, -2, -1, 0, 1, 2)) +  # Set specific breaks
  #geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") + # Adding a diagonal line
  theme(plot.title = element_text(hjust = 0.5))
plot
ggsave(paste(fileName,"_Original_Comparison_pos-hits_unvalidated.pdf",sep = ""), plot = plot, width = 11, height = 6, dpi = 300)

highlight_genes <- c("LAMTOR4", "GSPT1", "PGRMC2")
highlight_coords <- subset(data_merged, Gene %in% highlight_genes)
plot <- ggplot(data_merged, aes(x = Original.FC.Mageck, y = Average.U.Test.FC)) +
  geom_point(alpha = 0.4, size = 1.75, aes(color = Type)) +
  geom_point(data = highlight_coords, color = "red", alpha = 0.6, size = 1.75) +
  geom_text(aes(label = ifelse(Gene %in% highlight_genes, Gene, "")), vjust = -1, hjust = 0.5, color = "red") + # Label highlighted genes
  scale_color_manual(values = c("genome_wide_hit" = "blue", "other" = "black")) +
  theme_bw(base_size = 16) +
  labs(x = "Original Screen log2FC", y = "Validation Screen Log2FC") +
  ggtitle("D3 pos vs D3 neg Comparison") +
  annotate("text", x = -Inf, y = Inf, label = paste("R²:", round(r_squared2, 2)), hjust = -0.2, vjust = 2, size = 5, color = "red") +
  scale_x_continuous(breaks = c(-4, -3, -2, -1, 0, 1, 2)) +  # Set specific breaks
  #geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") + # Adding a diagonal line
  theme(plot.title = element_text(hjust = 0.5))
plot
ggsave(paste(fileName,"_Original_Comparison_non-hits.pdf",sep = ""), plot = plot, width = 11, height = 6, dpi = 300)

highlight_genes <- c("EIF3M", "C7ORF26", "SYVN1", "UPF2", "TP53", "CCM2", "MESDC1", "PIK3CA", "INO80E", "SUZ12")
highlight_coords <- subset(data_merged, Gene %in% highlight_genes)
plot <- ggplot(data_merged, aes(x = Original.FC.Mageck, y = Average.U.Test.FC)) +
  geom_point(alpha = 0.4, size = 1.75, aes(color = Type)) +
  geom_point(data = highlight_coords, color = "red", alpha = 0.6, size = 1.75) +
  geom_text(aes(label = ifelse(Gene %in% highlight_genes, Gene, "")), vjust = -1, hjust = 0.5, color = "red") + # Label highlighted genes
  scale_color_manual(values = c("genome_wide_hit" = "blue", "other" = "black")) +
  theme_bw(base_size = 16) +
  labs(x = "Original Screen log2FC", y = "Validation Screen Log2FC") +
  ggtitle("D3 pos vs D3 neg Comparison") +
  annotate("text", x = -Inf, y = Inf, label = paste("R²:", round(r_squared2, 2)), hjust = -0.2, vjust = 2, size = 5, color = "red") +
  scale_x_continuous(breaks = c(-4, -3, -2, -1, 0, 1, 2)) +  # Set specific breaks
  #geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") + # Adding a diagonal line
  theme(plot.title = element_text(hjust = 0.5))
plot
ggsave(paste(fileName,"_Original_Comparison_tri_genes.pdf",sep = ""), plot = plot, width = 9, height = 6, dpi = 300)



# Re-plot by omitting outliers
# Omit certain rows using filter() from dplyr
filtered_data <- data_merged %>% filter(!Gene %in% c("WDR12", "CENPC", "NUP93", "SMC2", "RPS15A", "SKA3"))
temp <- na.omit(filtered_data[, c("Original.FC.Mageck", "Average.U.Test.FC")])
lm_fit <- lm(Average.U.Test.FC ~ Original.FC.Mageck, data = temp)
r_squared <- summary(lm_fit)$r.squared # Extract R-squared value
plot <- ggplot(filtered_data, aes(x = Original.FC.Mageck, y = Average.U.Test.FC)) +
  geom_point(alpha = 0.4, size = 1.75, aes(color = Type)) +
  scale_color_manual(values = c("genome_wide_hit" = "blue", "other" = "black")) +
  theme_bw(base_size = 16) +
  labs(x = "Original Screen log2FC", y = "Validation Screen Log2FC") +
  ggtitle("D3 pos vs D3 neg Comparison") +
  annotate("text", x = -Inf, y = Inf, label = paste("R²:", round(r_squared, 2)), hjust = -0.2, vjust = 2, size = 5, color = "red") +
  scale_x_continuous(breaks = c(-4, -3, -2, -1, 0, 1, 2)) +  # Set specific breaks
  #geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") + # Adding a diagonal line
  theme(plot.title = element_text(hjust = 0.5))
plot

#################################
### Original U-test vs Mageck ###
#################################

temp <- na.omit(data[, c("Original.FC.U.Test", "Original.FC.Mageck")])
#correlation <- cor(temp$Original.FC, temp$Average.U.Test.FC)
lm_fit <- lm(Original.FC.U.Test ~ Original.FC.Mageck, data = temp)
r_squared <- summary(lm_fit)$r.squared # Extract R-squared value

plot <- ggplot(data, aes(x = Original.FC.U.Test, y = Original.FC.Mageck)) +
  geom_point(alpha = 0.4, size = 1.75, aes(color = hits)) +
  scale_color_manual(values = c("black", "red")) +
  theme_bw(base_size = 16) +
  labs(x = "Original U-test log2FC", y = "Origina Mageck Log2FC") +
  ggtitle("D3 pos vs D3 neg Comparison") +
  annotate("text", x = -Inf, y = Inf, label = paste("R²:", round(r_squared, 2)),
           hjust = -0.2, vjust = 2, size = 5, color = "red") +
  theme(plot.title = element_text(hjust = 0.5))
plot

ggsave(paste(fileName,"_Original_U-test_vs_Mageck.pdf",sep = ""), plot = plot, width = 8, height = 6, dpi = 300)




plot <- ggplot(data, aes(x = Average.U.Test.FC, y =  X.log..U.Test.P.)) +
  geom_point(alpha = 0.4, size = 1.75, aes(color = NonHit)) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"),
                     labels = c("TRUE" = "Non-Hit", "FALSE" = "Hit")) +
  labs(x = "Log2FC", y = "-Log10(P)", color = "Genes Tested") +
  theme_bw(base_size = 16) +
  labs(x = "Log2FC", y = "-Log10(P)") +
  ggtitle("D3 vs D0") +
  theme(plot.title = element_text(hjust = 0.5))









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
