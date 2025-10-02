# Library -----------------------------------------------------------------

library(readr)
library(tidyr)
library(dplyr)
library(Matrix)
library(tibble)
library(Rediscover)
library(maftools)
library(ComplexHeatmap)

# Load data and process data ----------------------------------------------

# Load the MAF
wgs_data <- read.delim("00_data/wgs_data/TCGA-BLCA-WGS-gt1percFreq.maf",
                       comment.char = "#", stringsAsFactors = FALSE)

#wgs_maf <- read.maf(maf = wgs_data)

## Remove the Silent mutations
#wgs_data <- wgs_data[wgs_data$Variant_Classification != "Silent", ]

# Remove non-coding mutations
wgs_data <- wgs_data[wgs_data$Variant_Classification %in%
                       c("Missense_Mutation", "RNA", "Splice_Region",
                         "Nonsense_Mutation", "In_Frame_Del", "Frame_Shift_Del",
                         "Splice_Site", "Nonstop_Mutation", "In_Frame_Ins",
                         "Translation_Start_Site", "Frame_Shift_Ins"),]

# Filter down the maf - get unique gene-sample combinations
mut_pres <- wgs_data %>% 
  select(Hugo_Symbol, Tumor_Sample_Barcode) %>% 
  distinct()  # Use distinct() instead of unique() for clarity

# Identify genes with at least 20 mutations (do this before creating the matrix)
genes_to_keep <- mut_pres %>% 
  dplyr::count(Hugo_Symbol, name = "n_mutations") %>% 
  filter(n_mutations >= 20) %>%
  arrange(desc(n_mutations))

# Filter mutation data to only include genes with >=20 mutations
mut_pres_filtered <- mut_pres %>%
  filter(Hugo_Symbol %in% genes_to_keep$Hugo_Symbol)

# Create binary matrix (samples as rows, genes as columns)
mut_matrix <- mut_pres_filtered %>% 
  mutate(present = 1) %>% 
  pivot_wider(names_from = Hugo_Symbol,
              values_from = present,
              values_fill = 0) %>%
  arrange(Tumor_Sample_Barcode)

sample_names <- mut_matrix$Tumor_Sample_Barcode

mut_matrix <- mut_matrix |> select(-Tumor_Sample_Barcode)
mut_matrix_m <- as.matrix(mut_matrix)
#rownames(mut_matrix_m) <- sample_names

gene_names <- colnames(mut_matrix)


# This is where you add the new lines to assign row and column names.
rownames(mut_matrix_m) <- sample_names
colnames(mut_matrix_m) <- gene_names


# Analysis ----------------------------------------------------------------

# Estimate mutation probabilities

PMA <- getPM(mut_matrix_m)

mut_ex <- getMutex(mut_matrix_m,PMA)

## Subset by threshold -----------------------------------------------------

# Plot --------------------------------------------------------------------

#somaticInteractions(maf = wgs_maf, top = 50, pvalue = c(1e-2, 2e-3))


Heatmap(as.matrix(mut_ex), name = "Mutex score")
