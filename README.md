
# MIBC Landscape Analysis

A comparative study of whole-genome sequencing (WGS) versus whole-exome sequencing (WXS) in characterising the mutational landscape of Muscle-Invasive Bladder Cancer (MIBC) using TCGA data.

**Key Finding:** WGS provides superior genomic characterisation of MIBC compared to WXS, offering enhanced structural variant detection, more accurate TMB measurements (7-fold difference), and identification of 133 regulatory hotspots missed by exome sequencing.

## Repository Structure

```
mibc_landscape/
├── 00_data/                    # Data directory (not tracked)
│   ├── wgs_data/              # WGS MAF files (408 TCGA-BLCA samples)
│   ├── wxs_data/              # WXS data from Robertson et al. 2017
│   ├── oncoKB_genelist/       # OncoKB cancer gene annotations
│   ├── cCRE_data/             # ENCODE cis-regulatory elements
│   └── gencode_data/          # GENCODE v48 gene annotations
├── 01_wxs-wg-comp/            # WXS vs WGS comparison analysis
├── 02_wgs-hotspots/           # Regulatory hotspot identification
├── 03_mut_exclu/              # Mutual exclusivity analysis
├── 04_mut_sig/                # Mutational signature analysis
├── 05_mut_burden/             # TMB and VAF comparison
├── 06_clustering/             # PCoA clustering analysis
└── scratch.R                  # Exploratory analysis script

```

## Requirements

### R Version

-   R >= 4.5.1
-   RStudio 2025.05.0 (Build 496)

### Required R Packages

```r
# Data manipulation
library(readr)        # v2.1.5
library(dplyr)
library(tidyr)
library(tidyverse)    # v2.0.0
library(tibble)
library(data.table)   # v1.17.8

# Genomic data
library(maftools)     # v2.25.10
library(rtracklayer)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

# Analysis tools
library(mutSignatures) # v2.1.5
library(Rediscover)
library(vegan)
library(Matrix)

# Visualisation
library(ggplot2)      # v3.5.2
library(ggrepel)
library(ComplexHeatmap)
library(NMF)

# Other
library(corpcor)
library(fs)

```

## Data Requirements

-   **408 TCGA-BLCA tumour samples** WGS data from the TCGA-BLCA project
-   **Matched WXS data** from Robertson et al. (2017)
-   **Reference genome:** GRCh38 (patch release 14)
-   **Gene annotations:** GENCODE release 48
-   **Regulatory elements:** ENCODE cCRE database
-   **Cancer gene list:** OncoKB (version 2025-05-20)

Place data files in the `00_data/` directory structure as shown above.

## Analysis Modules

### 1. WXS vs WGS Comparison (`01_wxs-wg-comp/`)

Compares mutation detection patterns between platforms. Filters to OncoKB oncogenes (≥10 occurrences), classifies mutations as platform-unique or shared, and analyses variant types.

**Key Results:** WGS shows superior sensitivity for insertions/deletions and exclusively detects dinucleotide/trinucleotide polymorphisms. TP53 showed 62.7% concordance between platforms.

### 2. Mutational Hotspot Analysis (`02_wgs-hotspots/`)

Identifies recurrent mutation sites (99th percentile, ≥33 mutations/position), maps to ENCODE cCRE elements, and associates with genes using distance-based rules.

**Key Results:** 133 hotspots identified including TERT promoter (230 mutations) and ADGRG6 enhancer (150 mutations). Distal enhancers most abundant (82 hotspots), promoters highest density (61.7 mutations/hotspot).

### 3. Mutual Exclusivity Analysis (`03_mut_exclu/`)

Creates binary mutation matrix for genes with ≥20 mutations, calculates mutex scores using Rediscover package.

**Results:** Mutex score heatmap generated. Analysis completed but not fully explored due to time constraints.

### 4. Mutational Signature Analysis (`04_mut_sig/`)

Extracts SNVs, annotates trinucleotide context, generates 96-category mutation matrix, performs signature deconvolution using SigProfiler/SIGNAL.

**Key Results:** Eight SBS signatures identified. SBS5 (smoking-related, median 1,156.73 mutations) and SBS3 (BRCA-deficiency, median 737.26 mutations) dominate. Note: Expected APOBEC signatures not prominent (requires investigation).

### 5. Tumour Mutational Burden (`05_mut_burden/`)

Calculates TMB per sample normalised by target size (~30 Mb), performs regression analysis, compares coverage depth and VAF distributions.

**Key Results:** WGS median 35 mutations/Mb vs WXS median 5 mutations/Mb (7-fold difference). WGS shows lower VAF (0.08-0.10) indicating superior subclonal detection. Clinical implication: WXS may misclassify patients for immunotherapy.

### 6. Clustering Analysis (`06_clustering/`)

Creates binary mutation matrix, calculates Jaccard distance, performs PCoA on all genes and oncogenes subset.

**Results:** PCoA plots showing sample relationships and variance explained by principal coordinates.

## Reproducing This Study

### Data Preparation

1.  Setup data directory:

```bash
mkdir -p 00_data/{wgs_data,wxs_data,oncoKB_genelist,cCRE_data,gencode_data}

```

2.  Obtain TCGA data (requires dbGaP authorisation)
    
3.  Download reference datasets (ENCODE cCRE, GENCODE v48, OncoKB, Robertson et al. 2017 WXS)
    

### Running the Analysis

Each analysis can be run independently:

```r
# Example: Run WXS vs WGS comparison
source("01_wxs-wg-comp/wes_wgs_comp.R")

# Run other modules similarly
source("02_wgs-hotspots/mut_freq_by_pos.R")
source("03_mut_exclu/wgs_mut_exc.R")
source("04_mut_sig/mut_sig.R")
source("05_mut_burden/mut_brdn.R")
source("06_clustering/clustering.R")

```

**Note:** Analyses are numbered logically but can be executed in any sequence.

## Output Files

Analysis results are organised in `outputs/` subdirectories within each module (not tracked by git):

### Module 1: WXS vs WGS Comparison

-   **Stacked bar chart:** Mutation type breakdown for top 25 genes (WXS-unique, shared, WGS-unique)
-   **Box plot:** Mutation count distribution by detection category
-   **Box plot:** Variant type distribution (SNP, DEL, INS, DNP, TNP) with log scale
-   **Box plot:** Variant classification distribution (missense, nonsense, frameshift, etc.) with log scale

### Module 2: Mutational Hotspot Analysis

-   **Bar chart:** Top 20 most frequently mutated genomic positions
-   **Bubble plot:** cCRE annotation metrics (hotspot count vs average mutations)
-   **Bar chart:** Unique hotspots per cCRE annotation category
-   **Bar chart:** Unique genes associated per cCRE annotation
-   **Bar chart:** Total mutations per protein-coding gene (top 25)
-   **Bar chart:** GO Biological Process enrichment results (top 50)
-   **Bar chart:** Filtered GO processes relevant to MIBC
-   **Text file:** `prot_gene.txt` - List of protein-coding genes associated with hotspots
-   **Text file:** GO enrichment analysis table

### Module 3: Mutual Exclusivity Analysis

-   **Heatmap:** Mutex scores showing gene pair relationships

### Module 4: Mutational Signature Analysis

-   **CSV file:** `mut_matrix.csv` - 96 trinucleotide mutation count matrix for all samples
-   **CSV file:** `avg.csv` - Average signature profile across cohort
-   **Box plot:** Mutational signature contributions (8 SBS signatures) with log scale

### Module 5: Tumour Mutational Burden

-   **Box plot:** TMB distribution comparison (WGS vs WXS)
-   **Scatter plot:** TMB correlation between platforms with linear regression fit
-   **Box plot:** Coverage depth comparison between platforms
-   **Box plot:** Variant allele frequency (VAF) comparison between platforms

### Module 6: Clustering Analysis

-   **PCoA plot:** Sample clustering using all genes
-   **PCoA plot:** Sample clustering using oncogenes only (OncoKB-filtered)

## Study Limitations (Technical)

1.  **MAF annotation bias:** Mutation annotation formats primarily designed for coding regions may influence non-coding variant interpretation
    
2.  **Computational requirements:** Substantial infrastructure needed for WGS data processing
    
3.  **Single cohort:** TCGA-BLCA dataset only; generalisability across populations not validated
    
4.  **Distance-based gene associations:** Regulatory hotspot-gene relationships based on genomic proximity may not reflect true regulatory interactions (requires ChIP-seq/DAP-seq validation)
    
5.  **Mutational signature methodology:** Absence of expected APOBEC signatures warrants further investigation
    

----------

**Note:** Detailed citations and acknowledgements are provided in the accompanying project report.

**Author:** Gagan Shetteppanavar Vishaya  
**Project Director:** Dr. Andrew Mason
