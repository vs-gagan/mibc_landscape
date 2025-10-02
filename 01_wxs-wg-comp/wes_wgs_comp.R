# Load Packages -----------------------------------------------------------
library(maftools)
library(tidyverse)
library(data.table)

# Load and Filter the WGS data --------------------------------------------
# Import the raw MAF file
wgs_data = read.maf('00_data/wgs_data/TCGA-BLCA-WGS.maf')
# Generate gene summary
wgs_gene_summary <- getGeneSummary(wgs_data)


#Import the oncoKB gene list
oncoKB_genelist = read_tsv('00_data/oncoKB_genelist/2025-05-20_oncoKB_list.tsv')


# Filter the gene summary with the oncoKB gene list and genes with a frequeny of
#  10 or more
wgs_gene_filtered <-
  wgs_gene_summary[Hugo_Symbol %in% oncoKB_genelist$Hugo_Symbol] |>
  arrange(desc(total)) |>
  filter(total >= 10)

# # save the gene names
# write.table(
#   wgs_gene_filtered$Hugo_Symbol,
#   '01_wxs-wg-comp/outputs/filtered_genelist.txt',
#   row.names = F,
#   col.names = F,
#   quote = F
# )

# calculate the percentage of samples with mutations in the gene
wgs_gene_filtered$mut_perc = wgs_gene_filtered$total / 408 * 100

# Create a filtered MAF file with the filtered gene list
filtered_wgs_data = subsetMaf(wgs_data, genes = wgs_gene_filtered$Hugo_Symbol)

# # plot summaries for the filtered MAF
# plotmafSummary(filteredMAF)
# oncoplot(filteredMAF)

# Load and filter WXS data -------------------------------------------
# Import the WXS MAF file
wxs_data <- read.maf('00_data/wxs_data/blca_tcga_pub_2017/data_mutations.txt')

# Generate gene summary
wxs_gene_summary <- getGeneSummary(wxs_data)

#Import the oncoKB gene list
oncoKB_genelist = read_tsv('00_data/oncoKB_genelist/2025-05-20_oncoKB_list.tsv')

# Filter the gene summary with the oncoKB gene list and genes with a frequeny of
#  10 or more
wxs_gene_filtered <-
  wxs_gene_summary[Hugo_Symbol %in% oncoKB_genelist$Hugo_Symbol] |>
  arrange(desc(total)) |>
  filter(total >= 10)

# calculate the percentage of samples with mutations in the gene
wxs_gene_summary$mut_perc = wxs_gene_summary$total / 412 * 100

# Create a filtered MAF file with the filtered gene list
filtered_wxs_data = subsetMaf(wxs_data, genes = wxs_gene_filtered$Hugo_Symbol)

# #MAF summary
# plotmafSummary(filtered_wxs_data)

# Trim and Prep data -----------------------------------------------

# Convert and prepare WXS data
wxs_trimmed <- filtered_wxs_data@data |>
  mutate(source = "WXS") |>
  select(
    Hugo_Symbol,
    Start_Position,
    End_Position,
    Tumor_Sample_Barcode,
    HGVSc
  ) |>
  as.data.table()

# Remove transcript prefix from HGVSc
wxs_trimmed[, HGVSc := sub(".*:(c\\..+)", "\\1", HGVSc)]

# Convert and prepare WGS data
wgs_trimmed <- filtered_wgs_data@data |>
  mutate(source = "WGS") |>
  select(
    Hugo_Symbol,
    Start_Position,
    End_Position,
    Tumor_Sample_Barcode,
    HGVSc
  ) |>
  as.data.table()

# Standardize Tumor_Sample_Barcode in WGS data (use first 4 parts, remove letters from 4th)
wgs_trimmed[,
  c("p1", "p2", "p3", "p4", "p5", "p6", "p7") := tstrsplit(
    Tumor_Sample_Barcode,
    "-",
    fixed = TRUE
  )
]
wgs_trimmed[, p4 := sub("[A-Z]+$", "", p4)]
wgs_trimmed[, Tumor_Sample_Barcode := paste(p1, p2, p3, p4, sep = "-")]
wgs_trimmed[, c("p1", "p2", "p3", "p4", "p5", "p6", "p7") := NULL]

# Combine unique genes of interest from both datasets
goi_list <- unique(c(wxs_trimmed$Hugo_Symbol, wgs_trimmed$Hugo_Symbol))


# Create Mutation Count Table by Gene and Origin (WES/WGS) ----------------

# Step 1: Create mutation ID for uniqueness
wxs_trimmed$mut_id <- paste(
  wxs_trimmed$Tumor_Sample_Barcode,
  wxs_trimmed$Hugo_Symbol,
  wxs_trimmed$HGVSc,
  sep = "_"
)
wgs_trimmed$mut_id <- paste(
  wgs_trimmed$Tumor_Sample_Barcode,
  wgs_trimmed$Hugo_Symbol,
  wgs_trimmed$HGVSc,
  sep = "_"
)

# Step 2: Label the source of each mutation
wxs_trimmed$source <- "WES"
wgs_trimmed$source <- "WGS"

# Step 3: Combine both datasets
merged_mutations <- rbind(
  wxs_trimmed[, .(Hugo_Symbol, mut_id, source)],
  wgs_trimmed[, .(Hugo_Symbol, mut_id, source)]
)

# Step 4: Create a table of how many times each mutation appears
#  (once = unique, twice = shared)
mutation_counts <- merged_mutations[, .N, by = .(Hugo_Symbol, mut_id)]
shared_mutations <- mutation_counts[N == 2] 

# Step 5: Count mutations by gene and category
# Total per gene (combined WES + WGS)
mut_per_gene <- merged_mutations[, .(total_mutations = length(unique(mut_id))), by = Hugo_Symbol]
#setnames(mut_per_gene, "N", "total_mutations")

# Shared mutations per gene
shared_by_gene <- shared_mutations[, .N, by = Hugo_Symbol]
setnames(shared_by_gene, "N", "shared_mutations")

# WXS-only mutations
wxs_only <- setdiff(wxs_trimmed$mut_id, wgs_trimmed$mut_id)
wxs_only_data <- wxs_trimmed[mut_id %in% wxs_only]
wxs_only_by_gene <- wxs_only_data[, .N, by = Hugo_Symbol]
setnames(wxs_only_by_gene, "N", "wxs_unique")

# WGS-only mutations
wgs_only <- setdiff(wgs_trimmed$mut_id, wxs_trimmed$mut_id)
wgs_only_data <- wgs_trimmed[mut_id %in% wgs_only]
wgs_only_by_gene <- wgs_only_data[, .N, by = Hugo_Symbol]
setnames(wgs_only_by_gene, "N", "wgs_unique")

# Step 6: Merge all summaries into one data frame
mutation_summary <- Reduce(
  function(...) merge(..., all = TRUE, by = "Hugo_Symbol"),
  list(mut_per_gene, shared_by_gene, wxs_only_by_gene, wgs_only_by_gene)
)

mutation_summary <- mutation_summary[order(-mutation_summary$total_mutations), ]

# Step 7: Replace NA with 0
mutation_summary[is.na(mutation_summary)] <- 0


# stats -------------------------------------------------------------------

top_gene_stats <- mutation_summary |> arrange(desc(total_mutations)) |>
  head(n=100) |> select(Hugo_Symbol,shared_mutations,wxs_unique,wgs_unique)

top_stat_long <- top_gene_stats |> pivot_longer(
  cols = c(shared_mutations,wxs_unique,wgs_unique),
  names_to = "mutation_type",
  values_to = "mutation_number"
)

ggplot(top_stat_long,
       aes(x=mutation_type,y = mutation_number,fill = mutation_type))+
  geom_boxplot()+
  scale_y_log10()+
  labs(title = "Distribution of Mutation Types",
       x = "Mutation Type",
       y = "Number of Mutations",
       fill = "Mutation Type") +
  scale_fill_manual(
    values = c(
      "wxs_unique" = "skyblue",
      "shared_mutations" = "darkgreen",
      "wgs_unique" = "salmon"
    ))+
  theme_minimal()

top_stat_summary <- top_gene_stats |> summarise(
  shared_mutations_mean = mean(shared_mutations),
  shared_mutations_sd = sd(shared_mutations),
  shared_mutations_min = min(shared_mutations),
  shared_mutations_max = max(shared_mutations),
  wxs_unique_mean = mean(wxs_unique),
  wxs_unique_sd = sd(wxs_unique),
  wxs_unique_min = min(wxs_unique),
  wxs_unique_max = max(wxs_unique),
  wgs_unique_mean = mean(wgs_unique),
  wgs_unique_sd = sd(wgs_unique),
  wgs_unique_min = min(wgs_unique),
  wgs_unique_max = max(wgs_unique)
) |>
  pivot_longer(
    cols = everything(),
    names_to = c("mutation_type", "stat"),
    names_sep = "_(?=[^_]+$)"  # splits at the last underscore
  ) |>
  pivot_wider(
    names_from = stat,
    values_from = value
  )


# WXS vs WGS plot ---------------------------------------------------------
# libraries
library(ggplot2)
library(reshape2) # For melt

# Select top n genes by total mutations
top_genes <- mutation_summary[order(-total_mutations)][1:25, Hugo_Symbol]

# Reshape to long format
plot_data <- melt(
  mutation_summary,
  id.vars = "Hugo_Symbol",
  measure.vars = c("wxs_unique", "shared_mutations", "wgs_unique"),
  variable.name = "Mutation_Type",
  value.name = "Count"
)

# Filter to top genes only
plot_data_top <- plot_data[plot_data$Hugo_Symbol %in% top_genes, ]

# Set factor order to put shared mutations in the middle of the stack
plot_data_top$Mutation_Type <- factor(
  plot_data_top$Mutation_Type,
  levels = c("wxs_unique", "shared_mutations", "wgs_unique")
)

# Plot
ggplot(
  plot_data_top,
  aes(x = reorder(Hugo_Symbol, -Count), y = Count, fill = Mutation_Type)
) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = c(
      "wxs_unique" = "skyblue",
      "shared_mutations" = "darkgreen",
      "wgs_unique" = "salmon"
    )
  ) +
  labs(
    title = "Mutation Type Breakdown per Gene (Top 25)",
    x = "Gene",
    y = "Number of Mutations",
    fill = "Mutation Type"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Type of Mutations -------------------------------------------------------

## Retrim the dataframes --------------------------------------------------

# Convert and prepare WXS data
wxs_trimmed <- filtered_wxs_data@data |>
  mutate(source = "WXS") |>
  select(
    Hugo_Symbol,
    Start_Position,
    End_Position,
    Tumor_Sample_Barcode,
    HGVSc,
    Variant_Classification,
    Variant_Type,
    t_ref_count,
    t_alt_count
  ) |>
  as.data.table()

# Remove transcript prefix from HGVSc
wxs_trimmed[, HGVSc := sub(".*:(c\\..+)", "\\1", HGVSc)]

# Convert and prepare WGS data
wgs_trimmed <- filtered_wgs_data@data |>
  mutate(source = "WGS") |>
  select(
    Hugo_Symbol,
    Start_Position,
    End_Position,
    Tumor_Sample_Barcode,
    HGVSc,
    Variant_Classification,
    Variant_Type,
    t_ref_count,
    t_alt_count
  ) |>
  as.data.table()

# Standardize Tumor_Sample_Barcode in WGS data (use first 4 parts, remove letters from 4th)
wgs_trimmed[,
            c("p1", "p2", "p3", "p4", "p5", "p6", "p7") := tstrsplit(
              Tumor_Sample_Barcode,
              "-",
              fixed = TRUE
            )
]
wgs_trimmed[, p4 := sub("[A-Z]+$", "", p4)]
wgs_trimmed[, Tumor_Sample_Barcode := paste(p1, p2, p3, p4, sep = "-")]
wgs_trimmed[, c("p1", "p2", "p3", "p4", "p5", "p6", "p7") := NULL]

wxs_trimmed <- wxs_trimmed[Hugo_Symbol %in% top_genes]
wgs_trimmed <- wgs_trimmed[Hugo_Symbol %in% top_genes]

## Variant Type Distribution Analysis -------------------------------------

# Step 1: Create mutation categories like you did before
wxs_trimmed$mut_id <- paste(wxs_trimmed$Tumor_Sample_Barcode, 
                            wxs_trimmed$Hugo_Symbol, 
                            wxs_trimmed$HGVSc, sep = "_")
wgs_trimmed$mut_id <- paste(wgs_trimmed$Tumor_Sample_Barcode, 
                            wgs_trimmed$Hugo_Symbol, 
                            wgs_trimmed$HGVSc, sep = "_")

# Step 2: Identify mutation categories
wxs_only_ids <- setdiff(wxs_trimmed$mut_id, wgs_trimmed$mut_id)
wgs_only_ids <- setdiff(wgs_trimmed$mut_id, wxs_trimmed$mut_id)
shared_ids <- intersect(wxs_trimmed$mut_id, wgs_trimmed$mut_id)

# Step 3: Create summary tables for Variant_Type
variant_type_summary <- rbind(
  wxs_trimmed[mut_id %in% wxs_only_ids, .(Category = "WXS_Unique", Variant_Type)],
  wgs_trimmed[mut_id %in% wgs_only_ids, .(Category = "WGS_Unique", Variant_Type)],
  wxs_trimmed[mut_id %in% shared_ids, .(Category = "Shared", Variant_Type)]
)

# Step 4: Count and visualize
variant_counts <- variant_type_summary[, .N, by = .(Category, Variant_Type)]


### Plot ------------------------------------------------------------------

variant_counts$Category <- factor(variant_counts$Category, 
                                  levels = c("WXS_Unique", "Shared", "WGS_Unique"))



ggplot(variant_counts, aes(x = Variant_Type, y = N, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("WXS_Unique" = "skyblue", 
                               "WGS_Unique" = "salmon",
                               "Shared" = "darkgreen")) +
  labs(title = "Variant Type Distribution by Category",
       x = "Variant Type",
       y = "Number of Mutations",
       fill = "Category") +
  theme_minimal() +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = N), position = position_dodge(width = 0.9), vjust = -0.3, size = 3)



## Variant Classification Distribution Analysis ---------------------------

# Similar process but with Variant_Classification
classification_summary <- rbind(
  wxs_trimmed[mut_id %in% wxs_only_ids, .(Category = "WXS_Unique", Variant_Classification)],
  wgs_trimmed[mut_id %in% wgs_only_ids, .(Category = "WGS_Unique", Variant_Classification)],
  wxs_trimmed[mut_id %in% shared_ids, .(Category = "Shared", Variant_Classification)]
)

classification_counts <- classification_summary[, .N, by = .(Category, Variant_Classification)]

### Plot ------------------------------------------------------------------

classification_counts$Category <- factor(classification_counts$Category, 
                                         levels = c("WXS_Unique", "Shared", "WGS_Unique"))

ggplot(classification_counts, aes(x = Variant_Classification, y = N, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("WXS_Unique" = "skyblue", 
                               "WGS_Unique" = "salmon",
                               "Shared" = "darkgreen")) +
  scale_y_log10() +  # Log scale for y-axis
  labs(title = "Variant Classification Distribution by Category (Log Scale)",
       x = "Variant Classification",
       y = "Number of Mutations (Log Scale)",
       fill = "Category") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = N), position = position_dodge(width = 0.9), vjust = -0.3, size = 2.5)


