knitr::opts_chunk$set(tidy=FALSE, cache=TRUE, dev="png", message=FALSE, error=FALSE, warning=TRUE)

if(!"MAGeCKFlute" %in% installed.packages()) BiocManager::install("MAGeCKFlute")
if(!"clusterProfiler" %in% installed.packages()) BiocManager::install("clusterProfiler")
if(!"ggplot2" %in% installed.packages()) BiocManager::install("ggplot2")

library(MAGeCKFlute)
library(clusterProfiler)
library(ggplot2)

## path to the gene summary file (required)
file1 = "/Users/cjmunger/Documents/230523_Pipeline/combined/mageck/Pos_vs_Neg/Pos_vs_Neg.gene_summary.txt"
                  
## path to the sgRNA summary file (optional)
file2 = "/Users/cjmunger/Documents/230523_Pipeline/combined/mageck/Pos_vs_Neg/Pos_vs_Neg.sgrna_summary.txt"

# Run FluteRRA with only gene summary file
FluteRRA(file1, proj="EC_Total_D0", organism="hsa", outdir = "./")

# Run FluteRRA with both gene summary file and sgRNA summary file
FluteRRA(file1, file2, proj="EC_Total_D0", organism="hsa", outdir = "./")

gdata = ReadRRA(file1)
head(gdata)


#LogFDR vs RRA score
gdata$LogFDR = -log10(gdata$FDR)
p1 = ScatterView(gdata, x = "Score", y = "LogFDR", label = "id", 
                 model = "volcano", top = 5)
print(p1)



#p-value vs RRA score
p2 = VolcanoView(gdata, x = "Score", y = "FDR", top=5,ylab = "-log10(FDR)", x_cutoff=0,Label = "id")
print(p2)


#Rankplot with names for top 5 strongest promoters/inhibitors of differentiation
gdata$Rank = rank(gdata$Score)
p1 = ScatterView(gdata, x = "Rank", y = "Score", label = "id", 
                 top = 5, auto_cut_y = TRUE, ylab = "Log2FC", 
                 groups = c("top", "bottom"))
print(p1)

#KEGG Pathway Enrichment Analysis
geneList= gdata$Score
names(geneList) = gdata$id
enrich_pos = EnrichAnalyzer(geneList = geneList[geneList>0.5], 
                            method = "HGT", type = "KEGG")
enrich_neg = EnrichAnalyzer(geneList = geneList[geneList< -0.5], 
                            method = "HGT", type = "KEGG")


EnrichedView(enrich_pos, mode = 1, top = 5, bottom = 0)
EnrichedView(enrich_neg, mode = 2, top = 0, bottom = 20)

print(enrich_neg$geneName)
print(enrich_pos$geneName)

