library(maftools)

big_maf <- read.delim("00_data/wgs_data/TCGA-BLCA-WGS-gt1percFreq.maf",
                      comment.char = "#", stringsAsFactors = FALSE)

barcode_list <- big_maf$Tumor_Sample_Barcode |> unique()

barcode_subset <- barcode_list[100:110]

test_maf <- big_maf[big_maf$Tumor_Sample_Barcode %in% barcode_subset,]

write.table(test_maf,
            file = "00_data/test.maf",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)