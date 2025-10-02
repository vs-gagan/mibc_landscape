# Loading libraries --------------------------------------------------------

library(readr)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.UCSC.hg19)
library(mutSignatures)
library(NMF)
library(corpcor)

# Load and filer data -----------------------------------------------------

# Load the MAF
wgs_data <- read.delim("00_data/wgs_data/TCGA-BLCA-WGS-gt1percFreq.maf",
                       comment.char = "#", stringsAsFactors = FALSE)

# Filter for SNVs only
wgs_data_snv <- wgs_data |> filter(Variant_Type == "SNP")

# load human genome
hg38 <- BSgenome.Hsapiens.UCSC.hg38

data("mutSigData")

wgs_data_snv_ctx <- attachContext(
  mutData       = wgs_data_snv,
  chr_colName   = "Chromosome",
  start_colName = "Start_Position",
  end_colName   = "End_Position",
  nucl_contextN = 3,   # trinucleotide
  BSGenomeDb    = hg38
)

wgs_data_snv_ctx <- removeMismatchMut(
  mutData         = wgs_data_snv_ctx,
  refMut_colName  = "Reference_Allele",
  context_colName = "context",
  refMut_format   = "N"
)

wgs_data_snv_ctx <- attachMutType(
  mutData       = wgs_data_snv_ctx,
  ref_colName   = "Reference_Allele",
  var_colName   = "Tumor_Seq_Allele2",  # use ALT allele column
  context_colName = "context"
)

maf_counts <- countMutTypes(
  mutTable         = wgs_data_snv_ctx,
  sample_colName   = "Tumor_Sample_Barcode",  # sample ID column in MAF
  mutType_colName  = "mutType"
)

count_df <- getCounts(maf_counts)



write.csv(count_df,'04_mut_sig/output/mut_matrix.csv')

# # Deconvoluting the signatures --------------------------------------------
# 
# count_df <- read.csv("04_mut_sig/output/mut_matrix.csv",row.names = 1, 
#                      check.names = FALSE)
# 
# mut_cnt<- as.mutation.counts(count_df)
# 
# cosmix <- getCosmicSignatures()
# cosmix.ovcar <- cosmix[c(13,4,2)]
# 
# pre_pro <- prelimProcessAssess(input = mut_cnt, approach = "counts")
# 
# params <- setMutClusterParams(num_processesToExtract = 3,
#                                   approach = "counts",
#                                   num_totIterations = 50,
#                                   num_parallelCores = 8, 
#                                   debug = FALSE,
#                                   algorithm = "alexa")
# 
# analysis <- decipherMutationalProcesses(input = mut_cnt,
#                                             params = params)
# sigs <- analysis$Results$signatures
# 
# mSign <- matchSignatures(mutSign = sigs, reference = cosmix)
# 
# deconv<- resolveMutSignatures(mut_cnt,sigs)
# 
# msigPlot(expos)
# 
# getSignatureIdentifiers(deconv$Results$count.result)

# avg ---------------------------------------------------------------------

count_df_avg <- count_df

count_df_avg$Average <- rowMeans(count_df, na.rm = TRUE)

count_df_avg <- count_df_avg |> select(Average)

write.csv(count_df_avg,'04_mut_sig/output/avg.csv')
