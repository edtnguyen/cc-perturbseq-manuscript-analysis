library(ggplot2)
library(readxl)
library(dplyr)
library(writexl)
install.packages("writexl")

setwd("/Users/dulguun/oak/Users/dulguun/231023_Validation_Screen_Files")

# make a baseline dataframe
fileName = "d3pos_d3neg_mageckLFC"
data <- read.csv(paste(fileName, '.csv', sep = ''))
subset_data <- data %>% 
  select(Gene, Original.FC.Mageck, Average.U.Test.FC) %>%       # select desired columns
  filter(!(Gene %in% c("NON-TARGETING", "SAFE-TARGETING")))    # remove desired rows
  
# add a column indicating whether or not a gene is a genome-wide screen hit
gw_hit_list <- read.csv('original_genome-wide_results/gw_245_hits.txt', sep = '\t')
gw_hit_list$Gene <- toupper(gw_hit_list$Gene)
gw_hit_list <- gw_hit_list %>% mutate(Gene = case_when(
  Gene == "RBM14" ~ "RBM14-RBM4",
  TRUE ~ Gene  # Keep other values unchanged
))

subset_data$Gene %in% gw_hit_list$Gene %>% sum()

subset_data <- subset_data %>% mutate(genome_wide_screen_hit = Gene %in% gw_hit_list$Gene)

# add a column indicating whether or not a gene is a validation screen hit
validation_hit_list_path <- "/Users/dulguun/oak/Users/dulguun/240522_EC_Enhancer_Screen_Design/Genome-wide Validation Screen/validation_u_test_significant.xlsx"
validation_hit_list <- read_excel(validation_hit_list_path, sheet = 'all u-test significant')

validation_hit_list$Gene %in% subset_data$Gene %>% sum()

subset_data <- subset_data %>% mutate(validation_screen_hit = Gene %in% validation_hit_list$Gene)

# rename columns
subset_data <- subset_data %>% 
  rename(
    gene = Gene,
    genome_wide_screen_log2fc = Original.FC.Mageck,
    validation_screen_log2fc = Average.U.Test.FC
  )

# Export the data frame to a CSV file
write.csv(subset_data, "ec_screen_hits.csv", row.names = FALSE)


#to find row that don't match in gene names
#false_rows <- validation_hit_list %>% filter(!Gene %in% subset_data$Gene)
