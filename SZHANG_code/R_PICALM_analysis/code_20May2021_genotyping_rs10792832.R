# Siwei 20 May 2021
# Genotyping to find het individuals on rs10792832 site

# init
library(readr)
library(readxl)
library(factoextra)

# load data
## load 60 MGS samples
MGS_60_ID <- read_excel("60_MGS_ID.xlsx")

## load in genotype data (each row is a sample)
## column 1-6 are not part of the genotyping data
## genotype data starts at column 7
df_genotype <- read_table2("MGS_AD_snps.ped", 
                           col_names = FALSE, 
                           col_types = cols(
                             .default = col_character(),
                             X3 = col_double(),
                             X4 = col_double(),
                             X5 = col_double(),
                             X6 = col_double()
                           ))

## load in SNP location data
## note the SNP positions are in GRCh37/hg19 !
df_SNP_location <- read_delim("MGS_AD_snps.map", 
                              "\t", escape_double = FALSE, col_names = FALSE, 
                              trim_ws = TRUE)

## extract SNP GT info
df_genotype_gt <- df_genotype[, 7:ncol(df_genotype)]

## separate GT info into two alleles
df_genotype_gt_allele_1 <- 
  data.frame(df_genotype_gt[, seq(1, ncol(df_genotype_gt), 2)], 
             stringsAsFactors = F)
df_genotype_gt_allele_2 <- 
  data.frame(df_genotype_gt[, seq(2, ncol(df_genotype_gt), 2)], 
             stringsAsFactors = F)

## assemble alleles
## use sapply with simplify = "array", will return a matrix
df_genotype_gt_assembled <- 
  sapply(1:5333, function(x) {
    paste(df_genotype_gt_allele_1[x, ],
          df_genotype_gt_allele_2[x, ],
          sep = ",")
  }, 
  simplify = "array")
## convert matrix to data frame
df_genotype_gt_assembled <- as.data.frame(df_genotype_gt_assembled)

## assign rownames and colnames as SNP names and sample names
rownames(df_genotype_gt_assembled) <- df_SNP_location$X2
colnames(df_genotype_gt_assembled) <- df_genotype$X1

## extract desired SNPs
## rs10792832
df_sub <- 
  df_genotype_gt_assembled[rownames(df_genotype_gt_assembled) %in% 'rs10792832', ]
df_sub_index <- unlist(df_sub)
## genotypes: "G,G" "A,G" "A,A" "0,0"
df_sub <- df_sub[, df_sub_index %in% "A,G"]

MGS_60_ID_het_rs10792832 <-
  MGS_60_ID[MGS_60_ID$ID %in% colnames(df_sub), ]

write.table(MGS_60_ID_het_rs10792832,
            file = "MGS_het_ID_rs10792832.txt",
            row.names = F, quote = F, sep = "\t")
