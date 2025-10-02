# Library -----------------------------------------------------------------

library(readr)
library(dplyr)
library(vegan)
library(ggplot2)

# Load and process data ---------------------------------------------------------------

# Load the MAF
wgs_dat <- read.delim("00_data/wgs_data/TCGA-BLCA-WGS-gt1percFreq.maf",
                       comment.char = "#", stringsAsFactors = FALSE)

# Remove the Silent mutations
wgs_dat <- wgs_dat[wgs_dat$Variant_Classification != "Silent",]

# Trim the maf
tumor_mutations <- wgs_dat |> select(Hugo_Symbol,Tumor_Sample_Barcode)

# Remove rows with missing values
tumor_mutations <- tumor_mutations[complete.cases(tumor_mutations), ]

# Create a binary matrix: rows = samples, columns = genes
# Remove duplicate gene-sample combinations before pivoting
mutation_matrix <- tumor_mutations %>%
  distinct(Hugo_Symbol, Tumor_Sample_Barcode) %>%  # Remove duplicates
  mutate(present = 1) %>%
  pivot_wider(names_from = Hugo_Symbol, 
              values_from = present, 
              values_fill = 0,
              id_cols = Tumor_Sample_Barcode)

# Convert to matrix format
mutation_matrix <- mutation_matrix %>%
  column_to_rownames("Tumor_Sample_Barcode") %>%  # Use this instead of rownames()
  as.matrix()

# PCoA analysis -----------------------------------------------------------

# Calculate distance matrix (Jaccard for binary data)
dist_matrix <- vegdist(mutation_matrix, method = "jaccard")

# Perform PCoA
pcoa_result <- cmdscale(dist_matrix, eig = TRUE, k = 2)

# Create plotting dataframe
pcoa_df <- data.frame(
  Sample = rownames(mutation_matrix),
  PC1 = pcoa_result$points[, 1],
  PC2 = pcoa_result$points[, 2]
)

# Create PCoA plot --------------------------------------------------------

# Calculate variance explained
var_explained <- pcoa_result$eig[1:2] / sum(pcoa_result$eig) * 100

# Plot
ggplot(pcoa_df, aes(x = PC1, y = PC2)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
    title = "PCoA of Tumor Samples Based on Mutation Profiles"
  ) +
  theme_minimal()

# Filter down to just oncogenes -------------------------------------------

# Import and extract the gene names
oncogenes_table <- read_tsv('00_data/oncoKB_genelist/2025-05-20_oncoKB_list.tsv')
oncogenes <- oncogenes_table$Hugo_Symbol

# Filter the onco genes
wgs_genes <- wgs_dat$Hugo_Symbol
oncogenes_filtered <- intersect(oncogenes, wgs_genes)

# Filter the MAF to just oncogenes
wgs_dat_onco <- wgs_dat[wgs_dat$Hugo_Symbol %in% oncogenes_filtered,]

# Trim the maf
tumor_mutations <- wgs_dat_onco |> select(Hugo_Symbol,Tumor_Sample_Barcode)

# Remove rows with missing values
tumor_mutations <- tumor_mutations[complete.cases(tumor_mutations), ]

# Create a binary matrix: rows = samples, columns = genes
# Remove duplicate gene-sample combinations before pivoting
mutation_matrix <- tumor_mutations %>%
  distinct(Hugo_Symbol, Tumor_Sample_Barcode) %>%  # Remove duplicates
  mutate(present = 1) %>%
  pivot_wider(names_from = Hugo_Symbol, 
              values_from = present, 
              values_fill = 0,
              id_cols = Tumor_Sample_Barcode)

# Convert to matrix format
mutation_matrix <- mutation_matrix %>%
  column_to_rownames("Tumor_Sample_Barcode") %>%  # Use this instead of rownames()
  as.matrix()

## PCoA analysis -----------------------------------------------------------

# Calculate distance matrix (Jaccard for binary data)
dist_matrix <- vegdist(mutation_matrix, method = "jaccard")

# Perform PCoA
pcoa_result <- cmdscale(dist_matrix, eig = TRUE, k = 2)

# Create plotting dataframe
pcoa_df <- data.frame(
  Sample = rownames(mutation_matrix),
  PC1 = pcoa_result$points[, 1],
  PC2 = pcoa_result$points[, 2]
)

## Create PCoA plot --------------------------------------------------------

# Calculate variance explained
var_explained <- pcoa_result$eig[1:2] / sum(pcoa_result$eig) * 100

# Plot
ggplot(pcoa_df, aes(x = PC1, y = PC2)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
    title = "PCoA of Tumor Samples Based on FIltered Mutation Profiles"
  ) +
  theme_minimal()

