# Library -----------------------------------------------------------------

library(readr)
library(dplyr)
library(tidyr)
library(tidyverse)
library(tibble)

# Load data ---------------------------------------------------------------

# LOad the MAF
wgs_data <- read.delim("00_data/wgs_data/TCGA-BLCA-WGS-gt1percFreq.maf",
                       comment.char = "#", stringsAsFactors = FALSE)

# Remove the Silent mutations
wgs_data <- wgs_data[wgs_data$Variant_Classification != "Silent",]

# filter down the maf
mut_pres <- wgs_data |> select(Hugo_Symbol,Tumor_Sample_Barcode) |> unique()

# Create binary matrix
mut_matrix <- mut_pres |> 
  mutate(present = 1) |> 
  pivot_wider(names_from = Tumor_Sample_Barcode,
              values_from = present,
              values_fill = list(present = 0)) |> 
  arrange(Hugo_Symbol)

# Transpose the matrix
mut_matrix <- mut_matrix %>%
  column_to_rownames(var = "Hugo_Symbol") %>%  # move gene names to rownames
  t() %>%                                      # transpose
  as.data.frame() %>%                          # convert back to data frame
  rownames_to_column(var = "Tumor_Sample_Barcode")  # move sample names to a column


# Filter data -------------------------------------------------------------

# Make a list of genes to be filtered 
genes <- mut_pres |> group_by(Hugo_Symbol) |> summarise(n=n()) |> arrange(desc(n))
genes <- genes[genes$n >= 20, ]

# Filter the matrix 
filtered_matrix <- mut_matrix %>%
  select(Tumor_Sample_Barcode, all_of(genes$Hugo_Symbol))

filtered_matrix$Unknown <- NULL



# Check matrix frequency --------------------------------------------------

matrix_freq <- colSums(filtered_matrix[, -1]) |> as.data.frame()
matrix_freq <- rownames_to_column(matrix_freq,"Gene")
colnames(matrix_freq) <- c('Gene','freq')
matrix_freq$freq <- matrix_freq$freq/4.07
