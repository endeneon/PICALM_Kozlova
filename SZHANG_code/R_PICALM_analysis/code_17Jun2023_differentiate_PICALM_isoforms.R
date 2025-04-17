# Siwei 17 Jun 2023
# identify which isoforms were affected in PICALM G/A/KD samples

# Re-analyse Alena's data using Kallisto pseudocounts

# init
{
  library(tximport)
  library(readxl)

  library(EnsDb.Hsapiens.v86)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)

  library(stringr)

  library(edgeR)
  library(DESeq2)

  library(sva)
  library(factoextra)

  library(ggplot2)
  library(ggrepel)

}

txdb <-
  makeTxDbFromGFF(file = "~/Data/Databases/Genomes/hg38/gencode.v35.annotation.gtf",
                  dataSource = "Gencode_v35",
                  organism = "Homo Sapiens",
                  taxonomyId = 9606)

k <- keys(x = txdb,
          keytype = "TXNAME")

###
kallisto_tsv_files <-
  dir(path = "../kallisto_output",
      pattern = ".*.tsv",
      recursive = T,
      full.names = T)
# kallisto_tsv_files <-
#   kallisto_tsv_files[str_detect(string = kallisto_tsv_files,
#                                 pattern = "CD04")]
kallisto_tsv_files <-
  kallisto_tsv_files[str_detect(string = kallisto_tsv_files,
                                pattern = "CD09")]

sample_names <-
  dir(path = "../kallisto_output",
      pattern = ".*.tsv",
      recursive = T,
      full.names = F)
# sample_names <-
#   sample_names[str_detect(string = sample_names,
#                           pattern = "CD04")]
sample_names <-
  sample_names[str_detect(string = sample_names,
                          pattern = "CD09")]
sample_names <-
  str_split(string = sample_names,
            pattern = "\\/",
            simplify = T)[, 1]
names(kallisto_tsv_files) <-
  sample_names


tx2gene <-
  AnnotationDbi::select(x = txdb,
                        keys = k,
                        columns = "GENEID",
                        keytype = "TXNAME")
PICALM_tx <-
  tx2gene[tx2gene$GENEID == "ENSG00000073921.18", ]

# direct quantify different transcripts rather than collapsing into genes ####

df_tx_raw <-
  tximport(files = kallisto_tsv_files,
           type = "kallisto",
           txIn = T,
           txOut = T,
           countsFromAbundance = "lengthScaledTPM",
           ignoreTxVersion = F,
           ignoreAfterBar = T)

sum(rowSums(df_tx_raw$abundance) == 0)

# import metadata  #####
df_metadata <-
  read_excel("../PICALM_RNA_samples_description.xlsx")
# df_metadata <-
#   df_metadata[df_metadata$Cell_line_ID == "CD04", ]
df_metadata <-
  df_metadata[df_metadata$Cell_line_ID == "CD09", ]

df_metadata$Sample_Name <-
  str_replace_all(string = df_metadata$Sample_Name,
                  pattern = "\\-",
                  replacement = "_")
df_metadata$line_clone_geno <-
  str_c(df_metadata$Cell_line_ID,
        df_metadata$Clone_ID,
        df_metadata$genotype_edited,
        sep = "_")

df_metadata <-
  df_metadata[match(colnames(df_tx_raw$abundance),
                    df_metadata$Sample_Name)
              ,]

# send data to deseq2 ####
dds <-
  DESeqDataSetFromTximport(txi = df_tx_raw,
                           colData = df_metadata,
                           design = ~ genotype_edited)
dds$genotype_edited <-
  relevel(dds$genotype_edited,
          ref = "unedited")

snow_param <-
  SnowParam(workers = 6,
            type = "FORK",
            progressbar = T)
register(snow_param,
         default = T)
registered()

dds <- DESeq(dds,
             BPPARAM = MulticoreParam(6))

resultsNames(dds)
# [1] "Intercept"                             "Cell_line_ID_CD09_vs_CD04"
# [3] "genotype_edited_CRISPRoff_vs_unedited" "genotype_edited_risk_vs_unedited"

# results_dds <-
#   results(object = dds,
#           contrast = c(0, 0, 1),
#           parallel = T,
#           BPPARAM = MulticoreParam(6))

results_dds <-
  results(object = dds,
          contrast = c("genotype_edited",
                       "risk",
                       "unedited"),
          parallel = T,
          BPPARAM = MulticoreParam(6))

# rownames(results_dds)

results_edited_vs_unedited <-
  results_dds[rownames(results_dds) %in% PICALM_tx$TXNAME, ]
results_edited_vs_unedited <-
  as.data.frame(results_edited_vs_unedited)
write.table(results_edited_vs_unedited,
            file = "PICALM_tx_risk_vs_nonrisk_CD09.tsv",
            quote = F, sep = "\t",
            row.names = T, col.names = T)

results_edited_vs_unedited_CD04 <-
  results_edited_vs_unedited
write.table(results_edited_vs_unedited_CD04,
            file = "PICALM_tx_risk_vs_nonrisk_CD04.tsv",
            quote = F, sep = "\t",
            row.names = T, col.names = T)

results_dds <-
  results(object = dds,
          contrast = c("genotype_edited",
                       "CRISPRoff",
                       "unedited"),
          parallel = T,
          BPPARAM = MulticoreParam(6))

##
results_CRISPRoff_vs_unedited <-
  results_dds[rownames(results_dds) %in% PICALM_tx$TXNAME, ]
results_CRISPRoff_vs_unedited <-
  as.data.frame(results_CRISPRoff_vs_unedited)
