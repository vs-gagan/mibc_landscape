# Library -----------------------------------------------------------------

library(maftools)
library(tidyr)
library(dplyr)
library(rtracklayer)
library(GenomicRanges)
library(ggplot2)


# Load maf files ----------------------------------------------------------

wgs_data <- read.delim("00_data/wgs_data/TCGA-BLCA-WGS-gt1percFreq.maf",
                       comment.char = "#", stringsAsFactors = FALSE) |> 
  data.frame()

wxs_data <- read.delim("00_data/wxs_data/blca_tcga_pub_2017/data_mutations.txt",
                       comment.char = "#", stringsAsFactors = FALSE) |> 
  data.frame()

# Add "chr" prefix to the Chromosome column in wxs_data
wxs_data$Chromosome <- paste0("chr", wxs_data$Chromosome)

# filtered maf data -----------------------------------------------------------

# filter for silent mutations
wgs_data_filtered <- wgs_data |> filter(Variant_Classification != "Silent")
wxs_data_filtered <- wxs_data |> filter(Variant_Classification != "Silent")

# Convert position columns to numeric
wxs_data_filtered$Start_Position <- as.numeric(wxs_data_filtered$Start_Position)
wxs_data_filtered$End_Position <- as.numeric(wxs_data_filtered$End_Position)



# Load the Gencode data
gencode_gtf <- import('00_data/gencode_data/gencode.v48.annotation.gtf')

gene_loci_gtf <- gencode_gtf[gencode_gtf$type == "gene"]

gene_loci_data <- data.frame(
  seqname = as.character(seqnames(gene_loci_gtf)),
  start = start(gene_loci_gtf),
  end = end(gene_loci_gtf),
  strand = as.character(strand(gene_loci_gtf)),
  gene_id = mcols(gene_loci_gtf)$gene_id,
  gene_name = mcols(gene_loci_gtf)$gene_name
)
# remove redundant variables for memory's sake
remove(wgs_data,wxs_data,gencode_gtf,gene_loci_gtf)

# Convert gene loci data to GRanges
gene_loci_gr <- GRanges(
  seqnames = gene_loci_data$seqname,
  ranges = IRanges(start = gene_loci_data$start, end = gene_loci_data$end),
  strand = gene_loci_data$strand,
  gene_id = gene_loci_data$gene_id,
  gene_name = gene_loci_data$gene_name
)

# Convert WGS data to GRanges
wgs_data_gr <- GRanges(
  seqnames = wgs_data_filtered$Chromosome,
  ranges = IRanges(start = wgs_data_filtered$Start_Position,
                   end = wgs_data_filtered$End_Position),
  strand = wgs_data_filtered$Strand
)


# Filter out rows with NA values in Start_Position or End_Position
wxs_data_filtered <- wxs_data_filtered[!is.na(wxs_data_filtered$Start_Position)
                                       & !is.na(wxs_data_filtered$End_Position), ]

# Convert WXS data to GRanges
wxs_data_gr <- GRanges(
  seqnames = wxs_data_filtered$Chromosome,
  ranges = IRanges(start = wxs_data_filtered$Start_Position,
                   end = wxs_data_filtered$End_Position),
  strand = wxs_data_filtered$Strand
)


# Find overlaps
wgs_gene_overlap <- findOverlaps(wgs_data_gr, gene_loci_gr)
wxs_gene_overlap <- findOverlaps(wxs_data_gr, gene_loci_gr)


# Filter WGS data
wgs_data_filtered_genes <- wgs_data_filtered[queryHits(wgs_gene_overlap), ]
wgs_data_filtered_genes <- wgs_data_filtered_genes |> 
  select(Hugo_Symbol,Start_Position,End_Position,Tumor_Sample_Barcode,HGVSc,Gene)
wxs_data_filtered_genes <- wxs_data_filtered[queryHits(wxs_gene_overlap), ]
wxs_data_filtered_genes <- wxs_data_filtered_genes |> 
  select(Hugo_Symbol,Start_Position,End_Position,Tumor_Sample_Barcode,HGVSc,i_dbNSFP_Ensembl_geneid)

colnames(wxs_data_filtered_genes)[colnames(wxs_data_filtered_genes) == 
                                    "i_dbNSFP_Ensembl_geneid"] <- "Gene"

# Remove version numbers from gene_id
gene_loci_data <- gene_loci_data %>%
  mutate(gene_id_clean = sub("\\..*", "", gene_id))

# Prepare a lookup table for gene_id and gene coordinates
gene_coords <- gene_loci_data %>%
  select(gene_id_clean, start, end)
colnames(gene_coords) <- c("Gene","gene_start","gene_end")

# Join with wgs_data_filtered_genes & wxs_data_filtered_genes
wgs_data_filtered_genes <- wgs_data_filtered_genes %>%
  left_join(gene_coords, by = "Gene")

wxs_data_filtered_genes <- wxs_data_filtered_genes %>%
  left_join(gene_coords, by = "Gene")

# Calculate TMB per tumor sample (not per gene) ------------------------

# For WGS data - count mutations per tumor sample
wgs_tmb <- wgs_data_filtered_genes %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(
    mutations = n(),
    .groups = 'drop'
  ) %>%
  mutate(
    # For WGS, normalize by total coding sequence length or genome size
    # Using approximate coding sequence length: ~30 Mb
    coding_length_mb = 30,
    tmb = mutations / coding_length_mb
  )

# For WXS data - count mutations per tumor sample  
wxs_tmb <- wxs_data_filtered_genes %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(
    mutations = n(),
    .groups = 'drop'
  ) %>%
  mutate(
    # For WXS, normalize by exome size: ~30 Mb
    exome_length_mb = 30,
    tmb = mutations / exome_length_mb
  )

# TMB comparison between WGS and WXS --------------------------------------

# Rename columns for clarity
wgs_tmb_clean <- wgs_tmb
names(wgs_tmb_clean)[names(wgs_tmb_clean) == "tmb"] <- "wgs_tmb"
wgs_tmb_clean <- wgs_tmb_clean[, c("Tumor_Sample_Barcode", "wgs_tmb")]

wxs_tmb_clean <- wxs_tmb
names(wxs_tmb_clean)[names(wxs_tmb_clean) == "tmb"] <- "wxs_tmb"
wxs_tmb_clean <- wxs_tmb_clean[, c("Tumor_Sample_Barcode", "wxs_tmb")]

wgs_tmb_clean$Tumor_Sample_Barcode <- substr(wgs_tmb_clean$Tumor_Sample_Barcode, 1, 15)

# Join the TMB data for samples present in both datasets
tmb_comparison <- wgs_tmb_clean %>%
  inner_join(wxs_tmb_clean, by = "Tumor_Sample_Barcode")

# Add some summary statistics
tmb_comparison <- tmb_comparison %>%
  mutate(
    tmb_diff = wgs_tmb - wxs_tmb,
    tmb_ratio = wgs_tmb / wxs_tmb
  )

# Print summary
cat("Number of samples in WGS:", nrow(wgs_tmb), "\n")
cat("Number of samples in WXS:", nrow(wxs_tmb), "\n")
cat("Number of matched samples:", nrow(tmb_comparison), "\n")

# View first few rows
head(tmb_comparison)

# Plot the TMB comparison -------------------------------------------------

# 1. Scatter plot comparing WGS vs WXS TMB
ggplot(tmb_comparison, aes(x = wxs_tmb, y = wgs_tmb)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  labs(
    title = "TMB Comparison: WGS vs WXS",
    x = "WXS TMB (mutations/Mb)",
    y = "WGS TMB (mutations/Mb)",
    subtitle = paste("n =", nrow(tmb_comparison), "matched samples")
  ) +
  theme_minimal()

# 2. Box plot comparison
tmb_long <- tmb_comparison %>%
  select(Tumor_Sample_Barcode, wgs_tmb, wxs_tmb) %>%
  pivot_longer(cols = c(wgs_tmb, wxs_tmb), 
               names_to = "Method", 
               values_to = "TMB") %>%
  mutate(Method = ifelse(Method == "wgs_tmb", "WGS", "WXS"))

ggplot(tmb_long, aes(x = Method, y = TMB, fill = Method)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(
    title = "TMB Distribution: WGS vs WXS",
    x = "Sequencing Method",
    y = "Tumor Mutational Burden (mutations/Mb)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# 3. Correlation analysis
correlation <- cor(tmb_comparison$wgs_tmb, tmb_comparison$wxs_tmb, 
                   use = "complete.obs")
cat("\nCorrelation between WGS and WXS TMB:", round(correlation, 3), "\n")
