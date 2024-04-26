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

{
  ### switch G3 and A5 march identities here
  df_metadata_raw$sample_type[1:3] <- "risk"
  df_metadata_raw$sample_type[4:6] <- "non-risk"

  ### reassign sample identities
  df_metadata_raw$sample_names[1] <- "CD04-A5-1_march"
  df_metadata_raw$sample_names[2] <- "CD04-A5-2_march"
  df_metadata_raw$sample_names[3] <- "CD04-A5-3_march"

  df_metadata_raw$sample_names[4] <- "CD04-G3-1_march"
  df_metadata_raw$sample_names[5] <- "CD04-G3-2_march"
  df_metadata_raw$sample_names[6] <- "CD04-G3-3_march"

}

df_master_4_DGE <-
  df_master_raw[, colnames(df_master_raw) %in% df_metadata_raw$sample_names]
colnames(df_master_4_DGE)

### also need to reassign the master data frame
# colnames(df_master_4_DGE)
{
  colnames(df_master_4_DGE)[1:6] <-
    c("CD04-G3-1_march",
      "CD04-G3-2_march",
      "CD04-G3-3_march",
      "CD04-A5-1_march",
      "CD04-A5-2_march",
      "CD04-A5-3_march")
}

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

# pca_cpm$x

df_metadata_4_plot <-
  df_metadata_4_DGE[match(x = rownames(pca_cpm$x),
                          table = df_metadata_4_DGE$sample_names)
                      , ]
df_metadata_4_plot$sample_names
rownames(pca_cpm$x)
all(colnames(df_master_4_DGE) == df_metadata_4_plot$sample_names)
# df_metadata_4_DGE[match(x = df_metadata_4_DGE$sample_names,
#                         table = rownames(pca_cpm$x)), ]$sample_names
# df_metadata_4_DGE[match(x = rownames(pca_cpm$x),
#                         table = df_metadata_4_DGE$sample_names), ]$sample_names

# ! use prcomp$x, it is the data set projected !on! the PC

# pca_cpm$rotation[, 1]
# rownames(pca_cpm$rotation[1, ])
# pca_cpm$x[, 1]
# rownames(pca_cpm$x[1, 1])

pca_cpm$x <-
  pca_cpm$x[!(df_metadata_4_plot$sample_names %in%
              c("CD04-G3-1_march",
                "CD04-G3-2_August",
                "CD04-A5-2_August"))
            , ]

df_metadata_4_plot <-
  df_metadata_4_plot[!df_metadata_4_plot$sample_names %in%
                      c("CD04-G3-1_march",
                        "CD04-G3-2_August",
                        "CD04-A5-2_August"), ]

pca_cpm$x <-
  pca_cpm$x[!((df_metadata_4_plot$sample_type %in%
                 "non-risk-CRISPRoff") &
                (str_detect(string = df_metadata_4_plot$sample_names,
                            pattern = "march",
                            negate = T)))
            , ]

rownames(pca_cpm$x)

df_metadata_4_plot <-
  df_metadata_4_plot[!((df_metadata_4_plot$sample_type %in%
                         "non-risk-CRISPRoff") &
                        (str_detect(string = df_metadata_4_plot$sample_names,
                                    pattern = "march",
                                    negate = T)))
                    , ]

### remove all non-risk-CRISPRoff samples #####
pca_cpm$x <-
  pca_cpm$x[!((df_metadata_4_plot$sample_type %in%
                 "non-risk-CRISPRoff"))
            , ]

df_metadata_4_plot <-
  df_metadata_4_plot[!((df_metadata_4_plot$sample_type %in%
                          "non-risk-CRISPRoff"))
                     , ]

all(rownames(pca_cpm$x) == df_metadata_4_plot$sample_names)

fviz_pca_ind(pca_cpm,
             geom = c("point", "text"),
             repel = T,
             pointsize = 2,
             addEllipses = T,
             habillage = factor(df_metadata_4_plot$sample_type),
             ellipse.level = 0.95,
             # palette = brewer.pal(n = 3,
             #                      name = "Dark2"),
             # palette = "aaas",
             palette = c("#3B4992",
                         "#008B45"),
             invisible = "quali") +
  # geom_text(guide = "none") +
  ggtitle(label = "Ellipse.level = 0.95 (default)")

# plot without CD04-A5-1_march
# pca_cpm$x <-
#   pca_cpm$x[str_detect(string = df_metadata_4_plot$sample_names,
#                             pattern = "CD04-A5-1_march",
#                             negate = T)
#             , ]
#
# rownames(pca_cpm$x)
#
# df_metadata_4_plot <-
#   df_metadata_4_plot[str_detect(string = df_metadata_4_plot$sample_names,
#                                  pattern = "CD04-A5-1_march",
#                                  negate = T)
#                      , ]
#
# all(rownames(pca_cpm$x) == df_metadata_4_plot$sample_names)
#
# fviz_pca_ind(pca_cpm,
#              geom = c("point", "text"),
#              repel = T,
#              pointsize = 2,
#              addEllipses = F,
#              habillage = factor(df_metadata_4_plot$sample_type),
#              # palette = brewer.pal(n = 3,
#              #                      name = "Dark2"),
#              palette = "aaas",
#              invisible = "quali") +
#   # geom_text(guide = "none") +
#   ggtitle(label = "Batches 2 and 3 of Alena's Microglia, \nwith SVA/ComBat")

### Re-calculate DGEList using the current gene and sample list #####
df_master_4_DGE_combat_sub <-
  df_master_4_DGE_combat[rownames(df_master_4_DGE_combat) %in%
                           rownames(pca_cpm$rotation),
                         colnames(df_master_4_DGE_combat) %in%
                           df_metadata_4_plot$sample_names]

# df_master_4_DGE_combat <- df_master_4_DGE

df_edgeR_DGE <-
  DGEList(counts = as.matrix(df_master_4_DGE_combat_sub),
          samples = df_metadata_4_plot$sample_name,
          group = df_metadata_4_plot$sample_type,
          genes = rownames(df_master_4_DGE_combat_sub),
          remove.zeros = T)

df_edgeR_DGE <-
  calcNormFactors(df_edgeR_DGE)

# cpm_cutoff <- 1
# cpm_gene_count <-
#   as.data.frame(cpm(df_edgeR_DGE,
#                     normalized.lib.sizes = T))

# try add line as a blocking factor
df_metadata_4_plot$cell_line <-
  str_split(string = df_metadata_4_plot$sample_names,
            pattern = "-",
            simplify = T)[, 1]

design_matrix <-
  model.matrix(~ 0 +
                 batch +
                 cell_line +
                 sample_type,
               data = df_metadata_4_plot)

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

write.table(exp_table,
            file = "risk_vs_nonrisk_exp_table_06Oct2023.txt",
            row.names = T, col.names = T,
            sep = "\t", quote = F)

df_metadata_4_plot_writeout <-
  df_metadata_4_plot
df_metadata_4_plot_writeout$cell_line <- NULL

df_metadata_4_plot_writeout$PICALM_CPM <-
  unlist(cpm(df_edgeR_DGE)[rownames(cpm(df_edgeR_DGE)) %in% "PICALM", ])
write.table(df_metadata_4_plot_writeout,
            file = "sample_table_w_PICALM_cpm_outliers_removed_match_exp_table_06Oct2023.txt",
            row.names = T, col.names = T,
            sep = "\t", quote = F)


write.table(exp_table,
            file = "nonrisk-CRISPRoff_vs_nonrisk_exp_table_09Oct2023.txt",
            row.names = T, col.names = T,
            sep = "\t", quote = F)

# ! results of replicate Oct 03 results and remove off-risk August group ####
# ! best results so far
save.image("Alena_batch_02_03_only_rep_03Oct2023_best_results_so_far_05Oct2023.RData")

## Do not use the code below! ######

### try use 2 groups only
### Re-calculate DGEList using the current gene and sample list #####
df_metadata_4_plot <-
  df_metadata_4_plot[str_detect(string = df_metadata_4_plot$sample_names,
                                pattern = "CD04-A5-1_march",
                                negate = T)
                     , ]
df_metadata_4_plot <-
  df_metadata_4_plot[!(df_metadata_4_plot$sample_type %in%
                        "non-risk-CRISPRoff")
                     , ]


df_master_4_DGE_combat_sub <-
  df_master_4_DGE_combat[rownames(df_master_4_DGE_combat) %in%
                           rownames(pca_cpm$rotation),
                         colnames(df_master_4_DGE_combat) %in%
                           df_metadata_4_plot$sample_names]

# df_master_4_DGE_combat <- df_master_4_DGE

df_edgeR_DGE <-
  DGEList(counts = as.matrix(df_master_4_DGE_combat_sub),
          samples = df_metadata_4_plot$sample_name,
          group = df_metadata_4_plot$sample_type,
          genes = rownames(df_master_4_DGE_combat_sub),
          remove.zeros = T)

# df_edgeR_DGE <-
#   calcNormFactors(df_edgeR_DGE)

# cpm_cutoff <- 1
# cpm_gene_count <-
#   as.data.frame(cpm(df_edgeR_DGE,
#                     normalized.lib.sizes = T))

# try add line as a blocking factor
df_metadata_4_plot$cell_line <-
  str_split(string = df_metadata_4_plot$sample_names,
            pattern = "-",
            simplify = T)[, 1]

design_matrix <-
  model.matrix(~ 0 +
                 batch +
                 cell_line +
                 sample_type,
               data = df_metadata_4_plot)

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

# non-risk-off vs nonrisk
df_edgeR_output <-
  glmQLFTest(df_edgeR_QLM,
             coef = 3)
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

save.image("Alena_batch_02_03_only_03Oct2023.RData")
