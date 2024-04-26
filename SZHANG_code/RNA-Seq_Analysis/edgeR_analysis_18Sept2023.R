# Siwei 18 Sept 2023
# Analyse Alena's RNASeq results in-house

# init ####
{
  library(edgeR)
  library(readr)
  library(readxl)

  library(Rfast)
  library(factoextra)
  library(dplyr)
  library(stringr)

  library(ggplot2)
  library(RColorBrewer)
}

# load data ####
# ! Novogene "xls" files are actually tab-delimited text files
df_raw <-
  read_delim("gene_count.xls",
             delim = "\t", escape_double = FALSE,
             trim_ws = TRUE)
df_raw <-
  df_raw[!duplicated(df_raw$gene_name), ]
# keep gene names as a separate vector
gene_name_list <-
  df_raw$gene_name

df_4_DGE <-
  df_raw[, 1:17]
df_4_DGE$gene_id <- NULL
rownames(df_4_DGE) <- gene_name_list

# make metadata list
df_metadata_raw <-
  read_excel("RNA-seq_samples_description_Alena_9-18-23.xlsx")


df_metadata <-
  df_metadata_raw[order(df_metadata_raw$sample_name_long), ]
df_metadata$sample_name <-
  str_split(string = df_metadata$sample_name_long,
            pattern = "_CD",
            simplify = T)[, 1]
df_metadata$type <-
  c(rep_len("non-risk",
            length.out = 3),
    rep_len("non-risk-off",
            length.out = 2),
    rep_len("risk",
            length.out = 6),
    rep_len("non-risk",
            length.out = 3),
    rep_len("non-risk-off",
            length.out = 2))

# make df_4_DGE consisitent
df_4_DGE <-
  df_4_DGE[, match(x = df_metadata$sample_name,
                   table = colnames(df_4_DGE))]
rownames(df_4_DGE) <- gene_name_list

df_edgeR_DGE <-
  df_4_DGE[!rowSums(df_4_DGE, na.rm = T) == 0, ]
rownames(df_edgeR_DGE) <-
  gene_name_list[!rowSums(df_4_DGE, na.rm = T) == 0]

df_edgeR_DGE <-
  DGEList(counts = as.matrix(df_edgeR_DGE),
          samples = df_metadata$sample_name,
          group = df_metadata$type,
          genes = rownames(df_edgeR_DGE),
          remove.zeros = T)

gene_lengths <-
  df_raw$gene_length[df_raw$gene_name %in% rownames(df_edgeR_DGE)]

gene_lookup_table <-
  df_raw[, 18:ncol(df_raw)]
gene_lookup_table <-
  gene_lookup_table[gene_lookup_table$gene_name %in% rownames(df_edgeR_DGE), ]
fpkm_count_log <-
  rpkm(df_edgeR_DGE,
       gene.length = gene_lookup_table$gene_length,
       normalized.lib.sizes = T,
       log = T)
fpkm_count_log <-
  as.data.frame(fpkm_count_log)
saveRDS(fpkm_count_log,
        file = "~/Data/FASTQ/CellNet/Alena_PICALM_fpkm_log.RDs")

# cpm_gene_count <-
#   as.data.frame(cpm(df_edgeR_DGE))
cpm_gene_count <-
  as.data.frame(cpm(df_edgeR_DGE,
                    normalized.lib.sizes = T))
cpm_gene_count[rownames(cpm_gene_count) %in% "PICALM", ]

# write.table(cpm_gene_count,
#             file = "Alena_PICALM_samples_gene_count_in_cpm_19Sept2023.tsv",
#             quote = F, sep = "\t",
#             row.names = T, col.names = T)

# Plot PCA #####
# cpm_gene_count <-
#   as.data.frame(cpm(df_edgeR_DGE))
cpm_cutoff <- 1

df_edgeR_DGE <-
  df_edgeR_DGE[(rowSums(cpm_gene_count[, df_metadata$type %in%
                                         "non-risk"] > cpm_cutoff) >= 5) |
                 (rowSums(cpm_gene_count[, df_metadata$type %in%
                                           "risk"] > cpm_cutoff) >= 5) |
                 (rowSums(cpm_gene_count[, df_metadata$type %in%
                                           "non-risk-off"] > cpm_cutoff) >= 3)
               , ]
nrow(df_edgeR_DGE)

df_edgeR_DGE <-
  calcNormFactors(df_edgeR_DGE)




design_glm <-
  model.matrix(~ 0 +
                 Cell_line_ID +
                 # Clone_ID +
                 type,
               data = df_metadata)
rownames(design_glm) <- colnames(df_edgeR_DGE)

df_edgeR_DGE <-
  estimateDisp(df_edgeR_DGE,
               design = design_glm,
               robust = T)

plotBCV(df_edgeR_DGE)

df_edgeR_GLM <-
  glmQLFit(df_edgeR_DGE,
           design = design_glm,
           robust = T)
plotQLDisp(df_edgeR_GLM)

design_glm

## calc risk vs non-risk ####
df_edgeR_GLM_risk_v_nonrisk <-
  glmQLFTest(df_edgeR_GLM,
             coef = 4)
summary(decideTestsDGE(df_edgeR_GLM_risk_v_nonrisk))
df_writeout <-
  df_edgeR_GLM_risk_v_nonrisk$table
df_writeout$FDR <-
  p.adjust(p = df_writeout$PValue,
           method = "fdr")
df_writeout[rownames(df_writeout) %in% "PICALM", ]
write.table(df_writeout,
            file = "Alena_PICALM_edgeR_GLM_risk_v_nonrisk_20Sept2023.tsv",
            quote = F, sep = "\t",
            row.names = T, col.names = T)

df_edgeR_GLM_offrisk_v_nonrisk <-
  glmQLFTest(df_edgeR_GLM,
             coef = 3)
summary(decideTestsDGE(df_edgeR_GLM_offrisk_v_nonrisk))
df_writeout <-
  df_edgeR_GLM_offrisk_v_nonrisk$table
df_writeout$FDR <-
  p.adjust(p = df_writeout$PValue,
           method = "fdr")
write.table(df_writeout,
            file = "Alena_PICALM_edgeR_GLM_offrisk_v_nonrisk_20Sept2023.tsv",
            quote = F, sep = "\t",
            row.names = T, col.names = T)

df_edgeR_GLM_risk_v_offrisk <-
  glmQLFTest(df_edgeR_GLM,
             contrast = c(0, 0, -1, 1))
summary(decideTestsDGE(df_edgeR_GLM_risk_v_offrisk))
df_writeout <-
  df_edgeR_GLM_risk_v_offrisk$table
df_writeout$FDR <-
  p.adjust(p = df_writeout$PValue,
           method = "fdr")
write.table(df_writeout,
            file = "Alena_PICALM_edgeR_GLM_risk_v_offrisk_20Sept2023.tsv",
            quote = F, sep = "\t",
            row.names = T, col.names = T)

df_edgeR_GLM_09vs04 <-
  glmQLFTest(df_edgeR_GLM,
             contrast = c(-1, 1, 0, 0))
summary(decideTestsDGE(df_edgeR_GLM_09vs04))
df_writeout <-
  df_edgeR_GLM_09vs04$table
df_writeout$FDR <-
  p.adjust(p = df_writeout$PValue,
           method = "fdr")
write.table(df_writeout,
            file = "Alena_PICALM_edgeR_GLM_09vs04_20Sept2023.tsv",
            quote = F, sep = "\t",
            row.names = T, col.names = T)


## risk vs non-risk
df_metadata_sub <-
  df_metadata[df_metadata$type %in% c("non-risk",
                                      "risk"), ]

df_edgeR_DGE <-
  df_4_DGE[!rowSums(df_4_DGE, na.rm = T) == 0, ]
df_edgeR_DGE <-
  df_edgeR_DGE[, df_metadata$type %in% c("non-risk",
                                         "risk")]
rownames(df_edgeR_DGE) <-
  gene_name_list[!rowSums(df_4_DGE, na.rm = T) == 0]

df_edgeR_DGE <-
  DGEList(counts = as.matrix(df_edgeR_DGE),
          samples = df_metadata_sub$sample_name,
          group = df_metadata_sub$type,
          genes = rownames(df_edgeR_DGE),
          remove.zeros = T)

gene_lengths <-
  df_raw$gene_length[df_raw$gene_name %in% rownames(df_edgeR_DGE)]

cpm_gene_count <-
  as.data.frame(cpm(df_edgeR_DGE))
cpm_cutoff <- 1

df_edgeR_DGE <-
  df_edgeR_DGE[rowSums(cpm_gene_count[, df_metadata_sub$type %in% "non-risk"] > cpm_cutoff) >= 2 |
                 rowSums(cpm_gene_count[, df_metadata_sub$type %in% "risk"] > cpm_cutoff) >= 2, ]
rownames(df_edgeR_DGE)

df_edgeR_DGE <-
  calcNormFactors(df_edgeR_DGE)

design_glm <-
  model.matrix(~ 0 +
                 Cell_line_ID +
                 # Clone_ID +
                 type,
               data = df_metadata_sub)
rownames(design_glm) <- colnames(df_edgeR_DGE)

df_edgeR_DGE <-
  estimateDisp(df_edgeR_DGE,
               design = design_glm,
               robust = T)

plotBCV(df_edgeR_DGE)

df_edgeR_GLM <-
  glmQLFit(df_edgeR_DGE,
           design = design_glm,
           robust = T)
plotQLDisp(df_edgeR_GLM)

design_glm

df_edgeR_GLM <-
  glmQLFTest(df_edgeR_GLM,
             coef = 3)


# rpkm_df_edgeR_DGE <-
#   rpkm(df_edgeR_DGE,
#        gene.length = gene_lengths,
#        normalized.lib.sizes = T)
# rpkm_df_edgeR_DGE <-
#   top_frac(x = as.data.frame(rpkm_df_edgeR_DGE),
#            n = 0.5,
#            wt = "rowSums")


pca_rpkm <-
  cpm(df_edgeR_DGE,
      normalized.lib.sizes = T)
pca_rpkm[rownames(pca_rpkm) %in% "PICALM", ]

pca_rpkm <-
  prcomp(t(as.matrix(pca_rpkm)),
         center = T,
         scale. = T)

fviz_pca_ind(pca_rpkm,
             geom = c("point", "text"),
             pointsize = 2,
             addEllipses = T,
             habillage = df_edgeR_DGE$samples$group,
             palette = brewer.pal(n = 5,
                                  name = "Dark2"),
             invisible = "quali")

# use glmFit
df_edgeR_NBGLM <-
  glmFit(df_edgeR_DGE,
         design = design_glm)

## calc risk vs non-risk ####
df_edgeR_NBGLM_risk_v_nonrisk <-
  glmLRT(df_edgeR_GLM,
             coef = 4)
summary(decideTestsDGE(df_edgeR_NBGLM_risk_v_nonrisk))
df_writeout <-
  df_edgeR_NBGLM_risk_v_nonrisk$table
df_writeout$FDR <-
  p.adjust(p = df_writeout$PValue,
           method = "fdr")
df_writeout[rownames(df_writeout) %in% "PICALM", ]
