# Library -----------------------------------------------------------------
library(tidyverse)
library(dplyr)
library(data.table)
library(readr)
library(ggplot2)
library(ggrepel)

# Sort mutation sites by freq ---------------------------------------------

# Import maf data directly
wgs_data <- read.delim("00_data/wgs_data/TCGA-BLCA-WGS-gt1percFreq.maf",
                       comment.char = "#", stringsAsFactors = FALSE)
# Trim the WGS data frame
wgs_data_trimmed = wgs_data |>
  select(Chromosome,Start_Position,End_Position,Variant_Type)

# Add a position column
wgs_data_trimmed$pos = 
  paste(wgs_data_trimmed$Chromosome,':',wgs_data_trimmed$Start_Position,
        sep = '')

# Create a new dataframe with positions frequency
pos_freq = wgs_data_trimmed |> count(pos) |>  arrange(desc(n))

pos_detail <- wgs_data_trimmed |> select(Chromosome, Start_Position,pos) |>
  unique()

pos_data <- left_join(pos_freq,pos_detail,by = "pos")

# checking the positions against cCREdatabase -----------------------------
ccre_data <- read_tsv("00_data/cCRE_data/GRCh38-cCREs.bed", 
                  col_names = c("Chromosome", "Start_Position", "End_Position",
                                "cCRE_ID", "accession", "annotations"))

# Convert to data.tables
pos_data <- as.data.table(pos_data)
ccre_data <- as.data.table(ccre_data)

# Perform non-equi join
hotspot_data <- ccre_data[pos_data, 
                      on = .(Chromosome, 
                             Start_Position <= Start_Position, 
                             End_Position >= Start_Position), 
                      nomatch = 0L]

# hotspot defination ------------------------------------------------------

dist <- hotspot_data$n
hotspot_threshold <- quantile(dist, 0.99)
hotspot_data_unfiltered <- hotspot_data
hotspot_data <- hotspot_data[n>=hotspot_threshold,]

# plot --------------------------------------------------------------------

ggplot(head(hotspot_data,n=20), aes(x = reorder(pos, n), y = n)) +
  geom_col() +
  coord_flip() +
  theme_minimal()

# Calculate weighted statistics
weighted_analysis <- hotspot_data %>%
  group_by(annotations) %>%
  summarise(
    hotspot_count = n(),
    total_mutations = sum(n),
    weighted_mean = sum(n * n) / sum(n),  # Weight each hotspot by its own frequency
    max_frequency = max(n),
    top_25_percent = quantile(n, 0.75)
  ) %>%
  arrange(desc(total_mutations))


# Bubble plot - size represents total mutations in that annotation
annotation_bubble <- hotspot_data %>%
  group_by(annotations) %>%
  summarise(
    hotspot_count = n(),
    total_mutations = sum(n),
    avg_mutations = mean(n)
  )

ggplot(annotation_bubble, aes(x = hotspot_count, y = avg_mutations)) +
  geom_point(aes(size = total_mutations), alpha = 0.7) +
  geom_text_repel(aes(label = annotations), 
                  box.padding = 0.5, 
                  max.overlaps = Inf) +
  labs(title = "cCRE Annotations: Hotspot Count vs Average Mutations",
       subtitle = "Bubble size = Total mutations across all hotspots",
       x = "Number of Hotspots", 
       y = "Average Mutations per Hotspot",
       size = "Total Mutations") +
  theme_minimal()

# stats -------------------------------------------------------------------

hotspot_data_stats <- hotspot_data |> summarise(
  mean = mean(n),
  sd = sd(n)
)

# Genes associated with hotspots -------------------------------------------

# Load the Gencode data
library(rtracklayer)
gencode_gtf <- import('00_data/gencode_data/gencode.v48.annotation.gtf')

# Filter for genes only and convert to data.table
genes <- gencode_gtf[gencode_gtf$type == "gene"]
genes_dt <- as.data.table(genes) |> 
  select(seqnames, start, end, gene_id, gene_name, gene_type) 
  
# Rename columns by name
names(genes_dt)[names(genes_dt) == "seqnames"] <- "Chromosome"
names(genes_dt)[names(genes_dt) == "start"] <- "Gene_Start"
names(genes_dt)[names(genes_dt) == "end"] <- "Gene_End"

# Add TSS positions and prepare for association
genes_strand <- gencode_gtf[gencode_gtf$type == "gene"]
strand_info <- data.table(
  gene_id = genes_strand$gene_id,
  strand = as.character(strand(genes_strand))
)

# Merge strand info with genes_dt
genes_dt <- merge(genes_dt, strand_info, by = "gene_id")

# Calculate TSS position
genes_dt[, TSS := ifelse(strand == "+", Gene_Start, Gene_End)]

# Clean chromosome names to match hotspot data
genes_dt[, Chromosome := as.character(Chromosome)]

# Function to associate hotspots with genes based on annotation rules
associate_hotspots_to_genes <- function(hotspot_dt, genes_dt) {
  
  results <- list()
  
  for(i in 1:nrow(hotspot_dt)) {
    hotspot <- hotspot_dt[i]
    chr <- hotspot$Chromosome
    pos <- hotspot$Start_Position
    annotation <- hotspot$annotations
    
    # Filter genes on same chromosome
    chr_genes <- genes_dt[Chromosome == chr]
    if(nrow(chr_genes) == 0) next
    
    # Calculate distances
    chr_genes[, distance_to_tss := abs(TSS - pos)]
    chr_genes[, distance_to_gene_body := pmin(
      abs(Gene_Start - pos), 
      abs(Gene_End - pos),
      ifelse(pos >= Gene_Start & pos <= Gene_End, 0, Inf)
    )]
    
    # Apply annotation-specific rules
    associated_genes <- switch(annotation,
                               "PLS" = {
                                 # Promoter-like: within 200bp of TSS
                                 chr_genes[distance_to_tss <= 200]
                               },
                               "pELS" = {
                                 # Proximal enhancer: within 2kb of TSS
                                 chr_genes[distance_to_tss <= 2000]
                               },
                               "dELS" = {
                                 # Distal enhancer: nearest gene within 1MB
                                 min_dist <- min(chr_genes$distance_to_tss)
                                 if(min_dist <= 1000000) {
                                   chr_genes[distance_to_tss == min_dist]
                                 } else {
                                   data.table() # empty if too far
                                 }
                               },
                               "CA-CTCF" = ,
                               "CA" = ,
                               "CA-TF" = ,
                               "CA-H3K4me3" = {
                                 # Chromatin architecture: nearest 3 genes within 500kb
                                 chr_genes[distance_to_tss <= 500000][order(distance_to_tss)][1:min(3, .N)]
                               },
                               "TF" = {
                                 # TF binding: nearest gene within 100kb
                                 min_dist <- min(chr_genes$distance_to_tss)
                                 if(min_dist <= 100000) {
                                   chr_genes[distance_to_tss == min_dist]
                                 } else {
                                   data.table()
                                 }
                               },
                               # Default case
                               {
                                 min_dist <- min(chr_genes$distance_to_tss)
                                 chr_genes[distance_to_tss == min_dist]
                               }
    )
    
    # Add hotspot information to the associated genes
    if(nrow(associated_genes) > 0) {
      associated_genes[, ':='(
        hotspot_pos = hotspot$pos,
        hotspot_chr = hotspot$Chromosome,
        hotspot_position = hotspot$Start_Position,
        hotspot_freq = hotspot$n,
        cCRE_annotation = hotspot$annotations,
        cCRE_ID = hotspot$cCRE_ID,
        association_rule = annotation
      )]
      
      results[[i]] <- associated_genes
    }
  }
  
  # Combine results
  if(length(results) > 0) {
    return(rbindlist(results, fill = TRUE))
  } else {
    return(data.table())
  }
}

# Run the association
hotspot_gene_associations <- associate_hotspots_to_genes(hotspot_data, genes_dt)

# Filter relevant gene types
hotspot_gene_associations <- hotspot_gene_associations[
  gene_type %in% c("protein_coding", "lncRNA", "miRNA", 
                   "transcribed_processed_pseudogene", 
                   "transcribed_unprocessed_pseudogene", "misc_RNA")
]

# Summary by annotation type
annotation_summary <- hotspot_gene_associations[, .(
  unique_hotspots = uniqueN(hotspot_pos),
  unique_genes = uniqueN(gene_name),
  total_associations = .N,
  avg_distance = round(mean(distance_to_tss, na.rm = TRUE)),
  median_distance = round(median(distance_to_tss, na.rm = TRUE)),
  total_mutations = sum(hotspot_freq)
), by = cCRE_annotation][order(-total_mutations)]

print("Summary by annotation type:")
print(annotation_summary)


## Plots --------------------------------------------------------------------

# Unique Hotspots per cCRE
ggplot(annotation_summary, aes(x = reorder(cCRE_annotation, -unique_hotspots), 
                               y = unique_hotspots)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_text(aes(label = unique_hotspots), vjust = -0.5, size = 3.5) +
  labs(title = "Unique Hotspots per cCRE Annotation",
       x = "cCRE Annotation",
       y = "Number of Unique Hotspots") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# Genes per cCRE
ggplot(annotation_summary, aes(x = reorder(cCRE_annotation, -unique_genes), 
                               y = unique_genes)) +
  geom_col(fill = "darkorange", alpha = 0.8) +
  geom_text(aes(label = unique_genes), vjust = -0.5, size = 3.5) +
  labs(title = "Unique Genes per cCRE Annotation",
       x = "cCRE Annotation",
       y = "Number of Unique Genes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))



# Top genes by mutation burden
top_genes <- hotspot_gene_associations[, .(
  associated_hotspots = uniqueN(hotspot_pos),
  total_mutations = sum(hotspot_freq),
  avg_hotspot_freq = round(mean(hotspot_freq)),
  max_hotspot_freq = max(hotspot_freq),
  annotation_types = paste(unique(cCRE_annotation), collapse = ","),
  avg_distance = round(mean(distance_to_tss))
), by = .(gene_name, gene_type)][order(-total_mutations)]

protein_coding_genes <- top_genes[gene_type == "protein_coding"] |> 
  select(gene_name, associated_hotspots, total_mutations, max_hotspot_freq)

print("Top 20 genes by mutation burden:")
print(head(top_genes, 20))

write.table(protein_coding_genes$gene_name,
            "02_wgs-hotspots/outputs/prot_gene.txt",
            row.names = F,col.names = F,quote = F)
## Plot --------------------------------------------------------------------

top_genes_top <- top_genes[order(-total_mutations)] |> 
  head(n=25)


# Total Mutations per Protein-Coding Gene
ggplot(top_genes_top, aes(x = reorder(gene_name, -total_mutations), 
                                     y = total_mutations)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_text(aes(label = total_mutations), vjust = -0.5, size = 3) +
  labs(title = "Total Mutations per Protein-Coding Gene",
       x = "Gene Name", y = "Total Mutations") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))


protein_coding_genes_top <- protein_coding_genes[order(-total_mutations)] |> 
  head(n=25)


# Total Mutations per Protein-Coding Gene
ggplot(protein_coding_genes_top, aes(x = reorder(gene_name, -total_mutations), 
                                 y = total_mutations)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_text(aes(label = total_mutations), vjust = -0.5, size = 3) +
  labs(title = "Total Mutations per Protein-Coding Gene",
       x = "Gene Name", y = "Total Mutations") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))

# go analysis -------------------------------------------------------------

go_data <- read.delim("02_wgs-hotspots/GO_Biological_Process_2025_table.txt",
                      header = TRUE, sep = "\t")

go_filtered <- go_data |> select(Term,Combined.Score) |> 
  arrange(desc(Combined.Score)) |> head(n=50)

go_filtered$Term <- sub("\\s*\\(GO:[0-9]+\\)", "", go_filtered$Term)

ggplot(go_filtered, aes(x = reorder(Term, Combined.Score), y = Combined.Score)) +
  geom_col(fill = "skyblue", color = "black") + # Add color and border to bars
  labs(
    title = "Top 50 GO Biological Processes by Combined Score",
    subtitle = "Enrichment analysis of gene hotspots",
    x = "GO Biological Process Term",
    y = "Combined Score"
  ) +
  coord_flip() + # Flip the axes for better readability of long labels
  theme_minimal() + # Use a clean theme
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_blank() # Remove horizontal grid lines
  )

mibc_focused <- go_filtered[go_data$Term %in% c(
  "Cellular Response to Caffeine (GO:0071313)",
  "Response to Alkaloid (GO:0043279)", 
  "Termination of RNA Polymerase II Transcription (GO:0006369)",
  "Cellular Response to Purine-Containing Compound (GO:0071415)",
  "Negative Regulation of Transforming Growth Factor Beta Production (GO:0071635)",
  "Embryonic Digestive Tract Morphogenesis (GO:0048557)",
  "Neuronal Action Potential (GO:0019228)",
  "Fibroblast Growth Factor Receptor Signaling Pathway (GO:0008543)",
  "Regulation of Fibroblast Growth Factor Receptor Signaling Pathway (GO:0040036)",
  "Negative Regulation of Fibroblast Growth Factor Receptor Signaling Pathway (GO:0040037)",
  "Antibody-Dependent Cellular Cytotoxicity (GO:0001788)",
  "Glycerolipid Metabolic Process (GO:0046486)",
  "Action Potential (GO:0001508)",
  "Antigen Processing and Presentation of Peptide Antigen via MHC Class II (GO:0002495)",
  "Myelination (GO:0042552)"
), ]

mibc_focused <- na.omit(mibc_focused)

ggplot(mibc_focused, aes(x = reorder(Term, Combined.Score), y = Combined.Score)) +
  geom_col(fill = "skyblue", color = "black") + # Add color and border to bars
  labs(
    title = "Filtered GO Biolocial Processes",
    subtitle = "Enrichment analysis of gene hotspots",
    x = "GO Biological Process Term",
    y = "Combined Score"
  ) +
  coord_flip() + # Flip the axes for better readability of long labels
  theme_minimal() + # Use a clean theme
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_blank() # Remove horizontal grid lines
  )
