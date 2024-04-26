# Siwei 29 Sept 2023
# plot PCA of Alena's microglia with iPS-derived microglia + human
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

  library(sva)
}

set.seed(42)
# load data ####
# all data aligned by STAR and GENCODE v35
# Alena_batch_Nov2022 <-
#   read_delim("STAR_aligned/Alena_microglia_LPS_Nov2022.txt",
#              delim = "\t", escape_double = FALSE,
#              trim_ws = TRUE)
# Alena_batch_Nov2022 <-
#   Alena_batch_Nov2022[, 1:14]

Alena_batch_march <-
  read_delim("STAR_aligned/Alena_microglia_March_STAR.txt",
             delim = "\t", escape_double = FALSE,
             trim_ws = TRUE)
Alena_batch_August <-
  read_delim("STAR_aligned/Alena_microglia_August_STAR.txt",
             delim = "\t", escape_double = FALSE,
             trim_ws = TRUE)
#
# iPS_microglia <-
#   read_delim("STAR_aligned/iPS_microglia_NG.txt",
#              delim = "\t", escape_double = FALSE,
#              trim_ws = TRUE)
#
# human_microglia <-
#   read_delim("STAR_aligned/raw_counts_human_microglia_phs.txt",
#              delim = "\t", escape_double = FALSE,
#              trim_ws = TRUE)

load("ENSG_gene_index.RData")


#####
# Alena_Nov2022_meta <-
#   read_xlsx("STAR_aligned/RNA-seq_samples_Alena_Whitney_Nov2022_description_Alena.xlsx")
# Alena_Nov2022_meta <-
#   Alena_Nov2022_meta[1:16, ]
# Alena_Nov2022_meta$genotype_edited <-
#   Alena_Nov2022_meta$`PICALM genotype`

Alena_march_meta <-
  read_xlsx("STAR_aligned/Alena_microglia_Mar_metadata.xlsx")
Alena_march_meta$genotype_edited[Alena_march_meta$genotype_edited %in% "unedited"] <-
  "non-risk"
Alena_march_meta$genotype_edited[Alena_march_meta$genotype_edited %in% "CRISPRoff"] <-
  "non-risk-CRISPRoff"
Alena_August_meta <-
  read_xlsx("RNA-seq_samples_description_Alena_9-18-23.xlsx")
Alena_August_meta$PICALM_genotype[13:16] <- "non-risk-CRISPRoff"
#


#####

df_master_raw <-
  merge(x = Alena_batch_march,
        y = Alena_batch_August,
        by = "Geneid")

df_master_raw$Geneid <-
  str_split(df_master_raw$Geneid,
            pattern = "\\.",
            simplify = T)[, 1]
colnames(df_master_raw)


df_master_raw <-
  df_master_raw[!duplicated(df_master_raw$Geneid), ]
df_master_raw <-
  merge(x = df_master_raw,
        y = ENSG_anno_gene_indexed,
        by = "Geneid")
df_master_raw$Geneid <- NULL
df_master_raw <-
  df_master_raw[!duplicated(df_master_raw$Gene_Symbol), ]
gene_list <-
  df_master_raw$Gene_Symbol
rownames(df_master_raw) <-
  df_master_raw$Gene_Symbol
df_master_raw$Gene_Symbol <- NULL

#####
colnames(df_master_raw)
colnames(df_master_raw) <-
  str_replace_all(string = colnames(df_master_raw),
                  pattern = "_",
                  replacement = "-")
colnames(df_master_raw)
colnames(df_master_raw) <-
  str_replace_all(string = colnames(df_master_raw),
                  pattern = "-$",
                  replacement = "")
colnames(df_master_raw)
colnames(df_master_raw) <-
  str_replace_all(string = colnames(df_master_raw),
                  pattern = "^A1-",
                  replacement = "CD09-A1-")
colnames(df_master_raw) <-
  str_replace_all(string = colnames(df_master_raw),
                  pattern = "^A1off-",
                  replacement = "CD09-A1off-")
colnames(df_master_raw) <-
  str_replace_all(string = colnames(df_master_raw),
                  pattern = "^A3-",
                  replacement = "CD09-A3-")
colnames(df_master_raw) <-
  str_replace_all(string = colnames(df_master_raw),
                  pattern = "^A5-",
                  replacement = "CD04-A5-")
colnames(df_master_raw) <-
  str_replace_all(string = colnames(df_master_raw),
                  pattern = "^G3-",
                  replacement = "CD04-G3-")
colnames(df_master_raw) <-
  str_replace_all(string = colnames(df_master_raw),
                  pattern = "^G3off-",
                  replacement = "CD04-G3off-")

colnames(df_master_raw)[1:18] <-
  str_c(colnames(df_master_raw)[1:18],
        "_march")
# colnames(df_master_raw)[(19 + 13):(34 + 13)] <-
#   str_remove(string = colnames(df_master_raw)[(19 + 13):(34 + 13)],
#              pattern = )
colnames(df_master_raw)[19:ncol(df_master_raw)] <-
  str_c(colnames(df_master_raw)[19:ncol(df_master_raw)],
        "_August")

colnames(df_master_raw)


##### Jump
df_metadata_raw <-
  data.frame(sample_names = c(Alena_Nov2022_meta$`Sample Name`,
                              Alena_march_meta$Sample_Name,
                              Alena_August_meta$Original_name,
                              iPS_microglia_meta$Run,
                              human_microglia_meta$Run),
             sample_type = c(Alena_Nov2022_meta$genotype_edited,
                             Alena_march_meta$genotype_edited,
                             Alena_August_meta$PICALM_genotype,
                             iPS_microglia_meta$cell_type_group,
                             human_microglia_meta$cell_type_group))
df_metadata_raw$batch <-
  c(rep_len("Alena_Nov2022",
            length.out = 16),
    rep_len("Alena_March",
            length.out = 18),
    rep_len("Alena_August",
            length.out = 16),
    rep_len("iPSC_microglia",
            length.out = 8),
    rep_len("human_microglia",
            length.out = 4))

# df_master_raw_backup <- df_master_raw
# df_metadata_raw_backup <- df_metadata_raw
# df_master_raw <- df_master_raw_backup
# df_metadata_raw <- df_metadata_raw_backup

colnames(df_master_raw)[(1 + 13):(18 + 13)] <-
  str_c(colnames(df_master_raw)[(1 + 13):(18 + 13)],
        "_march")
# colnames(df_master_raw)[(19 + 13):(34 + 13)] <-
#   str_remove(string = colnames(df_master_raw)[(19 + 13):(34 + 13)],
#              pattern = )
colnames(df_master_raw)[(19 + 13):(34 + 13)] <-
  str_c(colnames(df_master_raw)[(19 + 13):(34 + 13)],
        "_August")

df_metadata_raw$sample_names[(1 + 16):(18 + 16)] <-
  str_c(df_metadata_raw$sample_names[(1 + 16):(18 + 16)],
        "_march")
df_metadata_raw$sample_names[((19 + 16)):(34 + 16)] <-
  str_c(df_metadata_raw$sample_names[(19 + 16):(34 + 16)],
        "_August")

df_master_4_DGE <-
  df_master_raw[, colnames(df_master_raw) %in% df_metadata_raw$sample_names]
colnames(df_master_4_DGE)
df_metadata_4_DGE <-
  df_metadata_raw[df_metadata_raw$sample_names %in% colnames(df_master_4_DGE), ]

df_master_4_DGE_combat <-
  ComBat_seq(counts = as.matrix(df_master_4_DGE),
             batch = as.factor(df_metadata_4_DGE$batch))

# df_master_4_DGE_combat <- df_master_4_DGE

df_edgeR_DGE <-
  DGEList(counts = as.matrix(df_master_4_DGE_combat),
          samples = df_metadata_4_DGE$sample_name,
          group = df_metadata_4_DGE$sample_type,
          genes = rownames(df_master_4_DGE_combat),
          remove.zeros = T)

cpm_cutoff <- 1
cpm_gene_count <-
  as.data.frame(cpm(df_edgeR_DGE,
                    normalized.lib.sizes = T))

df_edgeR_DGE <-
  df_edgeR_DGE[(rowSums(cpm_gene_count[, df_metadata_4_DGE$batch %in%
                                         "Alena_Match"] > cpm_cutoff) >= 12) |
                 (rowSums(cpm_gene_count[, df_metadata_4_DGE$batch %in%
                                           "Alena_August"] > cpm_cutoff) >= 10) |
                 (rowSums(cpm_gene_count[, df_metadata_4_DGE$batch %in%
                                           "iPSC_microglia"] > cpm_cutoff) >= 5) |
                 (rowSums(cpm_gene_count[, df_metadata_4_DGE$batch %in%
                                           "human_microglia"] > cpm_cutoff) >= 2)
               , ]

cpm_2_plot <-
  cpm(df_edgeR_DGE,
      normalized.lib.sizes = T,
      log = T)

pca_cpm <-
  prcomp(t(as.matrix(cpm_2_plot)),
         center = T,
         scale. = T)

fviz_pca_ind(pca_cpm,
             geom = c("point"),
             pointsize = 2,
             addEllipses = F,
             habillage = df_metadata_4_DGE$sample_type,
             palette = brewer.pal(n = 11,
                                  name = "Paired"),
             invisible = "quali") +
  ggtitle(label = "All 3 batches of Alena's Microglia pooled, with SVA/ComBat",
          subtitle = paste("Note: exvivo and invitro-6h samples were generated by single-end",
                           "RNA-seq and have inherent discreptancies from paired-end ones,",
                           "they may appear more distant from the centre without ComBat normalisation.",
                           sep = "\n"))

save.image("batch_1_2_3_w_iPS_MG_human_MG.RData")

####
df_metadata_raw <-
  data.frame(sample_names = c(Alena_Nov2022_meta$`Sample Name`,
                              Alena_march_meta$Sample_Name,
                              Alena_August_meta$Original_name),
             sample_type = c(Alena_Nov2022_meta$genotype_edited,
                             Alena_march_meta$genotype_edited,
                             Alena_August_meta$PICALM_genotype))
df_metadata_raw$batch <-
  c(rep_len("Alena_Nov2022",
            length.out = 16),
    rep_len("Alena_March",
            length.out = 18),
    rep_len("Alena_August",
            length.out = 16))
df_metadata_raw

df_metadata_raw$stimuli <-
  c(Alena_Nov2022_meta$stimuli,
    rep_len(NA,
            length.out = length(Alena_march_meta$Sample_Name)),
    rep_len(NA,
            length.out = length(Alena_August_meta$Original_name)))

df_metadata_raw$sample_names[(1 + 16):(18 + 16)] <-
  str_c(df_metadata_raw$sample_names[(1 + 16):(18 + 16)],
        "_march")
df_metadata_raw$sample_names[((19 + 16)):(34 + 16)] <-
  str_c(df_metadata_raw$sample_names[(19 + 16):(34 + 16)],
        "_August")

df_metadata_raw$stimuli[is.na(df_metadata_raw$stimuli)] <- "unstimulated"
df_metadata_raw <-
  df_metadata_raw[df_metadata_raw$stimuli %in% "unstimulated", ]

df_master_4_DGE <-
  df_master_raw[, colnames(df_master_raw) %in% df_metadata_raw$sample_names]
colnames(df_master_4_DGE)
df_metadata_4_DGE <-
  df_metadata_raw[df_metadata_raw$sample_names %in% colnames(df_master_4_DGE), ]

df_master_4_DGE_combat <-
  ComBat_seq(counts = as.matrix(df_master_4_DGE),
             batch = as.factor(df_metadata_4_DGE$batch),
             group = as.factor(df_metadata_4_DGE$sample_type))

# df_master_4_DGE_combat <- df_master_4_DGE

df_edgeR_DGE <-
  DGEList(counts = as.matrix(df_master_4_DGE_combat),
          samples = df_metadata_4_DGE$sample_name,
          group = df_metadata_4_DGE$sample_type,
          genes = rownames(df_master_4_DGE_combat),
          remove.zeros = T)

df_edgeR_DGE <-
  calcNormFactors(df_edgeR_DGE)

cpm_cutoff <- 1
cpm_gene_count <-
  as.data.frame(cpm(df_edgeR_DGE,
                    normalized.lib.sizes = T))

df_edgeR_DGE <-
  df_edgeR_DGE[(rowSums(cpm_gene_count[, df_metadata_4_DGE$batch %in%
                                         "Alena_March"] > cpm_cutoff) >= 5) |
                 (rowSums(cpm_gene_count[, df_metadata_4_DGE$batch %in%
                                         "Alena_March"] > cpm_cutoff) >= 5) |
                 (rowSums(cpm_gene_count[, df_metadata_4_DGE$batch %in%
                                           "Alena_August"] > cpm_cutoff) >= 5)
               , ]

cpm_2_plot <-
  cpm(df_edgeR_DGE,
      normalized.lib.sizes = T,
      log = T)

pca_cpm <-
  prcomp(t(as.matrix(cpm_2_plot)),
         center = T,
         scale. = T)

fviz_pca_ind(pca_cpm,
             geom = c("point"),
             repel = T,
             pointsize = 2,
             addEllipses = F,
             habillage = df_metadata_4_DGE$sample_type,
             palette = brewer.pal(n = 6,
                                  name = "Dark2"),
             invisible = "quali") +
  # geom_text(guide = "none") +
  ggtitle(label = "All 3 batches of Alena's Microglia pooled, \nwith SVA/ComBat")

design_matrix <-
  model.matrix(~ 0 +
                 batch +
                 sample_type,
               data = df_metadata_4_DGE)


df_edgeR_DGE <-
  estimateDisp(df_edgeR_DGE,
               design = design_matrix,
               robust = T)

df_edgeR_QLM <-
  glmQLFit(df_edgeR_DGE,
           design = design_matrix,
           robust = T)

# risk vs nonrisk
df_edgeR_output <-
  glmQLFTest(df_edgeR_QLM,
             coef = 5)
summary(decideTestsDGE(df_edgeR_output))

# non-risk-off vs nonrisk
df_edgeR_output <-
  glmQLFTest(df_edgeR_QLM,
             coef = 4)
summary(decideTestsDGE(df_edgeR_output))


exp_table <-
  df_edgeR_output$table
exp_table$FDR <-
  p.adjust(exp_table$PValue,
           method = "fdr")
exp_table[rownames(exp_table) %in% "PICALM", ]

cpm_PICALM <-
  cpm_gene_count[rownames(cpm_gene_count) %in% "PICALM", ]
write.table(cpm_PICALM,
            file = "cpm_PICALM_all_post_ComBat.txt",
            quote = F, sep = "\t",
            row.names = T, col.names = T)
write.table(df_metadata_4_DGE,
            file = "PICALM_all_batches_unstimulated_only_metadata.txt",
            quote = F, sep = "\t",
            row.names = T, col.names = T)


##### what if CD04-G3-march and CD04-A5-march were labelled wrong? #####
# (interchanged identities)
# try switch them and see the results
####
df_metadata_raw <-
  data.frame(sample_names = c(Alena_march_meta$Sample_Name,
                              Alena_August_meta$Original_name),
             sample_type = c(Alena_march_meta$genotype_edited,
                             Alena_August_meta$PICALM_genotype))

df_metadata_raw$batch <-
  c(rep_len("Alena_March",
            length.out = 18),
    rep_len("Alena_August",
            length.out = 16))
df_metadata_raw


df_metadata_raw$sample_names[1:18] <-
  str_c(df_metadata_raw$sample_names[1:18],
        "_march")
df_metadata_raw$sample_names[19:nrow(df_metadata_raw)] <-
  str_c(df_metadata_raw$sample_names[19:nrow(df_metadata_raw)],
        "_August")

# renumber rownames
rownames(df_metadata_raw) <-
  1:nrow(df_metadata_raw)

### switch G3 and A5 march identities here
df_metadata_raw$sample_type[1:3] <- "risk"
df_metadata_raw$sample_type[4:6] <- "non-risk"

df_master_4_DGE <-
  df_master_raw[, colnames(df_master_raw) %in% df_metadata_raw$sample_names]
colnames(df_master_4_DGE)
df_metadata_4_DGE <-
  df_metadata_raw[df_metadata_raw$sample_names %in% colnames(df_master_4_DGE), ]
df_metadata_4_DGE <-
  df_metadata_4_DGE[match(x = colnames(df_master_4_DGE),
                          table = df_metadata_4_DGE$sample_names), ]

# inspect consistency
colnames(df_master_4_DGE)
df_metadata_4_DGE$sample_names
all(colnames(df_master_4_DGE) == df_metadata_4_DGE$sample_names)

df_master_4_DGE_combat <-
  ComBat_seq(counts = as.matrix(df_master_4_DGE),
             batch = as.factor(df_metadata_4_DGE$batch),
             group = as.factor(df_metadata_4_DGE$sample_type))

# df_master_4_DGE_combat <- df_master_4_DGE

df_edgeR_DGE <-
  DGEList(counts = as.matrix(df_master_4_DGE_combat),
          samples = df_metadata_4_DGE$sample_name,
          group = df_metadata_4_DGE$sample_type,
          genes = rownames(df_master_4_DGE_combat),
          remove.zeros = T)

df_edgeR_DGE <-
  calcNormFactors(df_edgeR_DGE)

cpm_cutoff <- 1
cpm_gene_count <-
  as.data.frame(cpm(df_edgeR_DGE,
                    normalized.lib.sizes = T))

df_edgeR_DGE <-
  df_edgeR_DGE[(rowSums(cpm_gene_count[, df_metadata_4_DGE$sample_type %in%
                                           "risk"] > cpm_cutoff) >= 8) |
                 (rowSums(cpm_gene_count[, df_metadata_4_DGE$sample_type %in%
                                           "non-risk"] > cpm_cutoff) >= 8) |
                 (rowSums(cpm_gene_count[, df_metadata_4_DGE$sample_type %in%
                                           "non-risk-CRISPRoff"] > cpm_cutoff) >= 7)
               , ]

nrow(df_edgeR_DGE)

cpm_2_plot <-
  cpm(df_edgeR_DGE,
      normalized.lib.sizes = T,
      log = T)

pca_cpm <-
  prcomp(t(as.matrix(cpm_2_plot)),
         center = T,
         scale. = T)

df_metadata_4_plot <-
  df_metadata_4_DGE[match(x = rownames(pca_cpm$x),
                          table = df_metadata_4_DGE$sample_names)
                      , ]
df_metadata_4_DGE$sample_names
rownames(pca_cpm$x)
# df_metadata_4_DGE[match(x = df_metadata_4_DGE$sample_names,
#                         table = rownames(pca_cpm$x)), ]$sample_names
# df_metadata_4_DGE[match(x = rownames(pca_cpm$x),
#                         table = df_metadata_4_DGE$sample_names), ]$sample_names

fviz_pca_ind(pca_cpm,
             geom = c("point", "text"),
             repel = T,
             pointsize = 2,
             addEllipses = F,
             habillage = factor(df_metadata_4_plot$sample_type),
             # palette = brewer.pal(n = 3,
             #                      name = "Dark2"),
             palette = "aaas",
             invisible = "quali") +
  # geom_text(guide = "none") +
  ggtitle(label = "All 3 batches of Alena's Microglia pooled, \nwith SVA/ComBat")

fviz_pca_ind(pca_cpm,
             geom = c("point", "text"),
             repel = T,
             pointsize = 2,
             addEllipses = F,
             habillage = factor(df_metadata_4_plot$sample_type),
             # palette = brewer.pal(n = 3,
             #                      name = "Dark2"),
             palette = "aaas",
             invisible = "quali") +
  # geom_text(guide = "none") +
  ggtitle(label = "Batch 2 and 3 of Alena's Microglia pooled, \nwith SVA/ComBat")


design_matrix <-
  model.matrix(~ 0 +
                 batch +
                 sample_type,
               data = df_metadata_4_DGE)


df_edgeR_DGE <-
  estimateDisp(df_edgeR_DGE,
               design = design_matrix,
               robust = T)

df_edgeR_QLM <-
  glmQLFit(df_edgeR_DGE,
           design = design_matrix,
           robust = T)

# risk vs nonrisk
df_edgeR_output <-
  glmQLFTest(df_edgeR_QLM,
             coef = 4)
summary(decideTestsDGE(df_edgeR_output))

# # non-risk-off vs nonrisk
# df_edgeR_output <-
#   glmQLFTest(df_edgeR_QLM,
#              coef = 3)
# summary(decideTestsDGE(df_edgeR_output))


exp_table <-
  df_edgeR_output$table
exp_table$FDR <-
  p.adjust(exp_table$PValue,
           method = "fdr")
exp_table[rownames(exp_table) %in% "PICALM", ]

cpm_PICALM <-
  cpm_gene_count[rownames(cpm_gene_count) %in% "PICALM", ]
write.table(cpm_PICALM,
            file = "cpm_PICALM_all_post_ComBat.txt",
            quote = F, sep = "\t",
            row.names = T, col.names = T)
write.table(df_metadata_4_DGE,
            file = "PICALM_all_batches_unstimulated_only_metadata.txt",
            quote = F, sep = "\t",
            row.names = T, col.names = T)

save.image("Alena_batch_02_03_only_03Oct2023.RData")
