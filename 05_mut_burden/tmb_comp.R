# WGS vs WXS Coverage and VAF Comparison
# Focused analysis comparing coverage depth and variant allele frequency

# Library -----------------------------------------------------------------

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

# Standardize chromosome naming -------------------------------------------

wxs_data$Chromosome <- paste0("chr", wxs_data$Chromosome)

# Filter for non-silent mutations ----------------------------------------

wgs_data_filtered <- wgs_data %>% 
  filter(Variant_Classification != "Silent")

wxs_data_filtered <- wxs_data %>% 
  filter(Variant_Classification != "Silent")

# Convert position columns to numeric for WXS data
wxs_data_filtered$Start_Position <- as.numeric(wxs_data_filtered$Start_Position)
wxs_data_filtered$End_Position <- as.numeric(wxs_data_filtered$End_Position)

# Filter out rows with NA values
wgs_data_filtered <- wgs_data_filtered[!is.na(wgs_data_filtered$Start_Position) & 
                                         !is.na(wgs_data_filtered$End_Position), ]

wxs_data_filtered <- wxs_data_filtered[!is.na(wxs_data_filtered$Start_Position) & 
                                         !is.na(wxs_data_filtered$End_Position), ]

# Load Gencode data for gene overlaps ------------------------------------

gencode_gtf <- import('00_data/gencode_data/gencode.v48.annotation.gtf')
gene_loci_gtf <- gencode_gtf[gencode_gtf$type == "gene"]

gene_loci_gr <- GRanges(
  seqnames = as.character(seqnames(gene_loci_gtf)),
  ranges = IRanges(start = start(gene_loci_gtf), end = end(gene_loci_gtf)),
  strand = as.character(strand(gene_loci_gtf))
)

# Find overlaps with gene regions ----------------------------------------

wgs_data_gr <- GRanges(
  seqnames = wgs_data_filtered$Chromosome,
  ranges = IRanges(start = wgs_data_filtered$Start_Position,
                   end = wgs_data_filtered$End_Position),
  strand = wgs_data_filtered$Strand
)

wxs_data_gr <- GRanges(
  seqnames = wxs_data_filtered$Chromosome,
  ranges = IRanges(start = wxs_data_filtered$Start_Position,
                   end = wxs_data_filtered$End_Position),
  strand = wxs_data_filtered$Strand
)

wgs_gene_overlap <- findOverlaps(wgs_data_gr, gene_loci_gr)
wxs_gene_overlap <- findOverlaps(wxs_data_gr, gene_loci_gr)

# Filter to gene regions only --------------------------------------------

wgs_data_filtered_genes <- wgs_data_filtered[queryHits(wgs_gene_overlap), ] %>%
  dplyr::select(Tumor_Sample_Barcode, t_depth, t_ref_count, t_alt_count)

wxs_data_filtered_genes <- wxs_data_filtered[queryHits(wxs_gene_overlap), ] %>%
  dplyr::select(Tumor_Sample_Barcode, t_ref_count, t_alt_count)

# Convert coverage columns to numeric for WXS data
wxs_data_filtered_genes$t_ref_count <- as.numeric(wxs_data_filtered_genes$t_ref_count)
wxs_data_filtered_genes$t_alt_count <- as.numeric(wxs_data_filtered_genes$t_alt_count)

# Calculate t_depth for WXS data
wxs_data_filtered_genes$t_depth <- wxs_data_filtered_genes$t_ref_count + wxs_data_filtered_genes$t_alt_count

# Replace infinite or NaN values with NA
wxs_data_filtered_genes$t_depth[!is.finite(wxs_data_filtered_genes$t_depth)] <- NA

# Remove redundant variables for memory
remove(wgs_data, wxs_data, gencode_gtf, gene_loci_gtf, gene_loci_gr, 
       wgs_data_gr, wxs_data_gr, wgs_gene_overlap, wxs_gene_overlap,
       wgs_data_filtered, wxs_data_filtered)

# Coverage Analysis -------------------------------------------------------

# Tumor depth coverage analysis
wgs_coverage <- wgs_data_filtered_genes %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(
    median_coverage = median(t_depth, na.rm = TRUE),
    mean_coverage = mean(t_depth, na.rm = TRUE),
    min_coverage = min(t_depth, na.rm = TRUE),
    max_coverage = max(t_depth, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(method = "WGS")

wxs_coverage <- wxs_data_filtered_genes %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(
    median_coverage = median(t_depth, na.rm = TRUE),
    mean_coverage = mean(t_depth, na.rm = TRUE),
    min_coverage = min(t_depth, na.rm = TRUE),
    max_coverage = max(t_depth, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(method = "WXS")

coverage_comparison <- bind_rows(wgs_coverage, wxs_coverage)

# VAF Analysis ------------------------------------------------------------

# VAF (Variant Allele Frequency) analysis
wgs_vaf <- wgs_data_filtered_genes %>%
  mutate(vaf = t_alt_count / (t_ref_count + t_alt_count)) %>%
  filter(!is.na(vaf) & is.finite(vaf)) %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(
    median_vaf = median(vaf, na.rm = TRUE),
    mean_vaf = mean(vaf, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(method = "WGS")

wxs_vaf <- wxs_data_filtered_genes %>%
  mutate(vaf = t_alt_count / (t_ref_count + t_alt_count)) %>%
  filter(!is.na(vaf) & is.finite(vaf)) %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(
    median_vaf = median(vaf, na.rm = TRUE),
    mean_vaf = mean(vaf, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(method = "WXS")

vaf_comparison <- bind_rows(wgs_vaf, wxs_vaf)

# Visualization -----------------------------------------------------------

# Plot 3: Coverage comparison
p3 <- ggplot(coverage_comparison, aes(x = method, y = median_coverage, fill = method)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  labs(
    title = "Coverage Comparison: WGS vs WXS",
    x = "Sequencing Method",
    y = "Median Coverage Depth"
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("WGS" = "#1f77b4", "WXS" = "#ff7f0e"))

# Plot 4: VAF comparison
p4 <- ggplot(vaf_comparison, aes(x = method, y = median_vaf, fill = method)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  labs(
    title = "Variant Allele Frequency Comparison: WGS vs WXS",
    x = "Sequencing Method",
    y = "Median VAF"
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("WGS" = "#1f77b4", "WXS" = "#ff7f0e"))

# Display plots
print(p3)
print(p4)

# Summary Statistics ------------------------------------------------------

cat("Coverage Summary:\n")
coverage_summary <- coverage_comparison %>%
  group_by(method) %>%
  summarise(
    avg_median_coverage = round(mean(median_coverage, na.rm = TRUE), 1),
    sd_median_coverage = round(sd(median_coverage, na.rm = TRUE), 1),
    .groups = 'drop'
  )
print(coverage_summary)

cat("\nVAF Summary:\n")
vaf_summary <- vaf_comparison %>%
  group_by(method) %>%
  summarise(
    avg_median_vaf = round(mean(median_vaf, na.rm = TRUE), 3),
    sd_median_vaf = round(sd(median_vaf, na.rm = TRUE), 3),
    .groups = 'drop'
  )
print(vaf_summary)