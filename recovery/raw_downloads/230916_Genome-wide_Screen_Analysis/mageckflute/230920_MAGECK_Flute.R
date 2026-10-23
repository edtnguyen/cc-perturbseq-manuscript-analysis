if(!"MAGeCKFlute" %in% installed.packages()) BiocManager::install("MAGeCKFlute")
if(!"clusterProfiler" %in% installed.packages()) BiocManager::install("clusterProfiler")
if(!"ggplot2" %in% installed.packages()) BiocManager::install("ggplot2")

library(MAGeCKFlute)
library(clusterProfiler)
library(ggplot2)

OUTDIR = "/Users/dulguun/oak/Users/dulguun/230916_Genome-wide_Screen_Analysis/mageckflute"
#################################################
## path to the gene summary file (required)
file1 = "/Users/dulguun/oak/Users/dulguun/230916_Genome-wide_Screen_Analysis/mageck/Pos_vs_Neg/Pos_vs_Neg.gene_summary_GeneNameCorrected.txt"

## path to the sgRNA summary file (optional)
file2 = "/Users/dulguun/oak/Users/dulguun/230916_Genome-wide_Screen_Analysis/mageck/Pos_vs_Neg/Pos_vs_Neg.sgrna_summary_GeneNameCorrected.txt"

# Run FluteRRA with only gene summary file
FluteRRA(file1, proj="Pos_vs_Neg", organism="hsa", outdir = "/Users/dulguun/oak/Users/dulguun/230916_Genome-wide_Screen_Analysis/Flute")

# Run FluteRRA with both gene summary file and sgRNA summary file
FluteRRA(file1, file2, proj="Pos_vs_Neg_sgRNA", organism="hsa", outdir = "/Users/dulguun/oak/Users/dulguun/230916_Genome-wide_Screen_Analysis/Flute")

gdata = ReadRRA(file1)
head(gdata)
sdata = ReadsgRRA(file2)
head(sdata)
#################################################
## Volcano plot

# LogFDR vs RRA score
# Score = Log2FC
gdata$LogFDR = -log10(gdata$FDR)
p1 = ScatterView(gdata, x = "Score", y = "LogFDR", label = "id", 
                 model = "volcano", top = 5)
p1
ggsave(paste0(OUTDIR,"/Pos_vs_Neg_volcano.png"), width = 8, height = 6)

# Log2FC vs p-value
p2 = VolcanoView(gdata, x = "Score", y = "FDR", Label = "id", x_cutoff = 0)
p2
ggsave(paste0(OUTDIR,"/Pos_vs_Neg_volcano_2.png"), width = 8, height = 6)

sdata$LogFDR = -log10(sdata$FDR)
p1 = ScatterView(sdata, x = "LFC", y = "LogFDR", label = "Gene", 
                 model = "volcano", top = 5)
p1
ggsave(paste0(OUTDIR,"/Pos_vs_Neg_sgrna_volcano.png"), width = 8, height = 6)

#################################################
## Rank plot
# Rank all the genes based on their scores and label genes in the rank plot.

genes_to_label = c("CDH5","KDR","SMAD4","ACVR2A","SMAD1","DBR1","SKA3","CENPC","RPS15A","SMC2","WDR12","NUP93","CNOT2","CBLL1")
top_genes <- c("DBR1", "SMAD1", "ACVR2A", "CDH5", "SMAD5", "LEMD3", "SMARCE1", "PLCG1", "BMPR1A", "NF2", "DROSHA", "KDR", "LENG8", "SMARCC1", "TOP1", "SMARCA4", "ATXN2L", "NAA25", "AMOTL2", "TMTC3", "STAG2", "SMAD4", "SRPRB", "ETV2", "CHORDC1", "PRKCD", "ASCC3", "HAND1", "CUL3", "SRPR", "ARIH1", "SRP54", "TBX3", "PAXBP1", "NEDD8-P1P2", "MED1", "KIAA1432", "MIXL1", "TAOK1", "ERCC2", "OGT", "KIAA1551", "LATS2", "C5orf22", "IPO8", "DGCR8", "CSNK2B", "FBXW11", "CNEP1R1", "RNF20", "RAD21", "SCAF8", "G3BP1", "GTF2A2", "CASP3", "BAX", "TAF4", "MAU2", "RBM33", "MEMO1", "RABGGTB", "APAF1", "VHL", "CCNC", "DPY30", "TRRAP", "SMARCB1", "SMARCD3", "KCTD10", "CD2BP2", "NEDD8", "MED12", "PTBP1", "GNB2L1", "ACVR2B", "KDM1A", "PPIL4", "HUWE1", "EIF4G1", "RAB6A", "EIF5A", "SETD1B", "SRP68", "GTF2H4", "TADA2B", "SUPT20H", "PSMG4", "RBM4", "ASH2L", "SRP19", "NUFIP2", "PRMT5", "SMARCA5", "VWA9", "SNRNP25", "T", "ZBTB48", "SETDB1", "RNF40", "DPM2", "WDR77", "RBX1", "SNRNP27", "TCEB3", "SSFA2", "RBBP5", "SF1", "PPP2R1A", "HSD17B12", "GTF3C2", "POLR3H", "PLXNA2", "DPM1", "HDAC2", "NAA20", "WIZ", "CYCS", "ATP6AP1", "TAF2", "HNRNPM", "ATP5F1", "SAV1", "LARP4", "PPHLN1", "RBM14-P2", "TAF3", "SPTY2D1", "KDM3B", "EIF4G2", "XRN1", "EOMES", "RBM7", "ZC3H8", "ACVR1B", "EXOSC8", "WAPAL", "TBP", "LHFPL2", "SMAD3", "CCDC101", "BAD", "DDX42", "LAGE3", "RAB18", "FAM208A", "SRSF6", "SCAF4", "WAC", "POLR3K", "VPS54", "CAPN5", "BRD3", "EEF2", "NCBP1", "DPM3", "TBC1D9", "TAF10", "MRPL32", "RNF115", "GPI", "MANBAL", "GTF2A1", "INTS10", "PTEN", "EWSR1", "POLR2M", "EMC2", "TRIM24", "XPO5", "FUT7", "COPS6", "CCNT2", "PSMB3", "SLC29A1", "PBLD", "TAF6L", "SSB", "SPPL2B", "BCAS2", "COPS4", "LSM12", "C2orf53", "ZBTB22", "MBD1", "EPT1", "EIF3F", "PYGO2", "SLC35C1", "WWC3", "HINT2", "GCK", "RLF", "OSBP", "CACNG1", "SOX4", "RIBC1", "C16orf93", "RGP1", "EHMT2", "WDR61", "GPR107", "ATMIN", "NCEH1", "SRRM2", "MED19", "ACAA1", "SSRP1", "PEG10", "SPIN1", "RFX3", "CACFD1", "HNRNPU", "KIAA0232", "PTAR1", "EIF4E", "MED15", "UFM1", "ACTL7B", "TLE6", "SNUPN", "DNAJC17", "CCT6B", "BAG6", "CD99L2", "SMG5", "PSMD4", "PITPNC1", "TSTD2", "PRKD1", "ARTN", "PAM16", "ERO1L", "POLE4", "STK40", "UXT", "CTR9", "CHMP1")
bottom_genes <- c("NUP93", "CNOT2", "SMC2", "SKA3", "RPS15A", "CENPC", "WDR12", "CBLL")

gdata$Rank = rank(gdata$Score)
gdata$group <- "other"
gdata$group[gdata$id %in% top_genes] <- "top"
gdata$group[gdata$id %in% bottom_genes] <- "bottom"

p1 = ScatterView(gdata, x = "Rank", y = "Score", label = "id", 
                 auto_cut_y = TRUE, ylab = "Log2FC", 
                 groups = c("top", "bottom"), toplabels = genes_to_label, group_col = c("darkgreen", "red"), top = 8)
p1
ggsave(paste0(OUTDIR,"/Pos_vs_Neg_rank_4.png"), width = 5, height = 8)

#################################################
### EC total vs D0 ### 
#################################################
## path to the gene summary file (required)
file3 = "/Users/dulguun/oak/Users/dulguun/230916_Genome-wide_Screen_Analysis/mageck/EC_total_vs_D0/EC_total_vs_D0.gene_summary_GeneNameCorrected.txt"

## Run FluteRRA with only gene summary file
FluteRRA(file3, proj="EC_total_vs_D0", organism="hsa", outdir = "/Users/dulguun/oak/Users/dulguun/230916_Genome-wide_Screen_Analysis/mageckflute")

gdata_2 = ReadRRA(file3)
head(gdata_2)

## Volcano plot
# LogFDR vs RRA score
# Score = Log2FC
gdata_2$LogFDR = -log10(gdata_2$FDR)
p1 = ScatterView(gdata_2, x = "Score", y = "LogFDR", label = "id", xlab = "Log2FC", ylab = "-log10(p-value)",
                 model = "volcano", top = 5, group_col = c("red", "darkgreen"))
p1
ggsave(paste0(OUTDIR,"/EC_total_vs_D0_volcano.png"), width = 7, height = 6)

# Log2FC vs p-value
p2 = VolcanoView(gdata_2, x = "Score", y = "FDR", Label = "id", group_col = c("red", "darkgreen"))
p2
ggsave(paste0(OUTDIR,"/EC_total_vs_D0_volcano_2.png"), width = 8, height = 6)
