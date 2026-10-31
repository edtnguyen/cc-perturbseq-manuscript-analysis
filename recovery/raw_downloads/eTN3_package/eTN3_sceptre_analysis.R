# eTN3 Split-seq analysis

library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)
library(grid)
library(gridExtra)
library(reshape2)
library(ggrastr)


##############################################################################################################################

#devtools::install_github("katsevich-lab/sceptre")

library(sceptre)

# 1. import data



gene_matrix_path <- ("/Users/juliannguyen/Rdir/eTN3_EC_pilot/data/WT_combined/MOI8/DGE_filtered/count_matrix.mtx")
grna_matrix_path <- ("/Users/juliannguyen/Rdir/eTN3_EC_pilot/data/crispr_combined/MOI1/guide_RNAs_filtered/count_matrix.mtx")
all_genes_path <- ("/Users/juliannguyen/Rdir/eTN3_EC_pilot/data/WT_combined/MOI8/DGE_filtered/all_genes.csv")
all_grnas_path <- ("/Users/juliannguyen/Rdir/eTN3_EC_pilot/data/crispr_combined/MOI8/guide_RNAs_filtered/all_guides.csv")




all_grnas_path <- ("/Users/juliannguyen/Rdir/eTN5_DAhits/crispr_combined/sample_D3/guide_RNAs_filtered/all_guides.csv")



grna_group_data_frame <- data.frame(read.csv("/Users/juliannguyen/Rdir/eTN3_EC_pilot/params/eTN3_gRNA_group_table.csv"))


sceptre_object_highmoi <- import_data_from_parse(gene_mat_fp = gene_matrix_path,
                                                 grna_mat_fp = grna_matrix_path,
                                                 all_genes_fp = all_genes_path,
                                                 all_grnas_fp = all_grnas_path,
                                                 moi = "high",
                                                 grna_target_data_frame = grna_group_data_frame)



# 2-3. set analysis parameters, assign gRNAs

discovery_pairs_highmoi <- data.frame(read.csv("/Users/juliannguyen/Rdir/eTN3_EC_pilot/params/eTN3_gRNA_group_pairs.csv"))


positive_control_pairs_highmoi <- construct_positive_control_pairs(
   sceptre_object = sceptre_object_highmoi
 )
 
head(positive_control_pairs_highmoi)


## Step 2: Determine the set of response-gRNA group pairs to analyze
# set of response-gRNA group pairs to test for association
# each gRNA group is mapped to each gene
#response_grna_group_pairs <- generate_all_pairs(response_matrix,grna_group_data_frame)
# or 
#response_grna_group_pairs <- data.frame(read.csv("params/eTN3_gRNA_group_pairs.csv"))

sceptre_object_highmoi <- set_analysis_parameters(
  sceptre_object = sceptre_object_highmoi,
  discovery_pairs = discovery_pairs_highmoi,
#  positive_control_pairs = positive_control_pairs,
  side = "both", control_group = "complement")
# grna_integration_strategy = "singleton",
#grna_integration_strategy = "bonferroni",

print(sceptre_object_highmoi)

plot_grna_count_distributions(
  sceptre_object = sceptre_object_highmoi,
  n_grnas_to_plot = 9
)



plot_covariates(sceptre_object_highmoi, p_mito_threshold = 0.22)


sceptre_object_highmoi <- assign_grnas(
  sceptre_object = sceptre_object_highmoi,
  method = "mixture",
#  threshold = 10,
 # method = "thresholding",
# method = "maximum",
  parallel = TRUE
 )

#grnas_to_plot <- c("OCT4_g+22_single_emprl",
#                   "POU5F1_single_JW_lib1",
#                   "NT_00642_single_JW_lib1",
#                   "CDKN1A_single_JW_lib1",
#                   "NANOG+g21_single_emprl")

#plot(sceptre_object_highmoi, grnas_to_plot = grnas_to_plot)

plot(sceptre_object_highmoi)

# 4. QC 

sceptre_object_highmoi <- run_qc(
  sceptre_object = sceptre_object_highmoi,
  p_mito_threshold = 0.22
)

plot(sceptre_object_highmoi)


#5. Calibration check


sceptre_object_highmoi <- run_calibration_check(
  sceptre_object = sceptre_object_highmoi,
  parallel = FALSE
)

plot(sceptre_object_highmoi)

calibration_result <- get_result(
  sceptre_object = sceptre_object_highmoi,
  analysis = "run_calibration_check"
)
head(calibration_result)

#6 Discovery analysis

sceptre_object_highmoi <- run_discovery_analysis(
  sceptre_object = sceptre_object_highmoi,
  parallel = FALSE
)

plot(sceptre_object_highmoi)

discovery_result <- get_result(
  sceptre_object = sceptre_object_highmoi,
  analysis = "run_discovery_analysis"
)

discovery_result

write.csv(discovery_result, ""/Users/juliannguyen/Rdir/eTN3_EC_pilot/targetKD.csv", row.names = FALSE)

