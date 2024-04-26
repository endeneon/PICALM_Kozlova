# run HOMER enrichment of iMG/non-iMG SNPs for Jubao, size +/-250
# Siwei 24 Apr 2024

# init #####
library(readxl)

# load iMG SNPs #####
df_raw <-
  read_excel("AD_ASoC_iMG_for_HOMER_TF_motif_enrichment.xlsx",
             sheet = 1)

df_writeout <-
  df_raw

write.table(df_writeout,
            file = "iMG_SNPs_24Apr2024.bed",
            quote = F, sep = "\t",
            row.names = F, col.names = F)

# load non-iMG SNPs #####
df_raw <-
  read_excel("AD_ASoC_iMG_for_HOMER_TF_motif_enrichment.xlsx",
             sheet = 2)

df_writeout <-
  df_raw

write.table(df_writeout,
            file = "non_iMG_SNPs_24Apr2024.bed",
            quote = F, sep = "\t",
            row.names = F, col.names = F)
