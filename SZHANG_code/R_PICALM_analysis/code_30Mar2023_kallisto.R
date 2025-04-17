# Siwei 30 Mar 2023

# Re-analyse Alena's data using Kallisto pseudocounts

# init
library(tximport)
library(readxl)

library(EnsDb.Hsapiens.v86)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

library(stringr)

library(edgeR)
library(sva)
library(factoextra)

library(ggplot2)
library(ggrepel)


###
kallisto_tsv_files <-
  dir(path = "../kallisto_output",
      pattern = ".*.tsv",
      recursive = T,
      full.names = T)

sample_names <-
  dir(path = "../kallisto_output",
      pattern = ".*.tsv",
      recursive = T,
      full.names = F)
sample_names <-
  str_split(string = sample_names,
            pattern = "\\/",
            simplify = T)[, 1]
names(kallisto_tsv_files) <-
  sample_names

# txdb <- EnsDb.Hsapiens.v86

# to keep transcript ID consistency (since it has .xx).
# use txdb created from Gencode v35 gtf file

txdb <-
  makeTxDbFromGFF(file = "~/Data/Databases/Genomes/hg38/gencode.v35.annotation.gtf",
                  dataSource = "Gencode_v35",
                  organism = "Homo Sapiens",
                  taxonomyId = 9606)
keytypes(txdb)
# [1] "CDSID"    "CDSNAME"  "EXONID"   "EXONNAME" "GENEID"   "TXID"     "TXNAME"
columns(txdb)
# [1] "CDSCHROM"   "CDSEND"     "CDSID"      "CDSNAME"    "CDSPHASE"   "CDSSTART"
# [7] "CDSSTRAND"  "EXONCHROM"  "EXONEND"    "EXONID"     "EXONNAME"   "EXONRANK"
# [13] "EXONSTART"  "EXONSTRAND" "GENEID"     "TXCHROM"    "TXEND"      "TXID"
# [19] "TXNAME"     "TXSTART"    "TXSTRAND"   "TXTYPE"

k <- keys(x = txdb,
          keytype = "TXNAME")
tx2gene <-
  AnnotationDbi::select(x = txdb,
                        keys = k,
                        columns = "GENEID",
                        keytype = "TXNAME")


# We can set ignoreTxVersion=T to ignore the transcript version after .
df_raw <-
  tximport(files = kallisto_tsv_files,
           type = "kallisto",
           tx2gene = tx2gene,
           ignoreAfterBar = T)

ref_db_gene_symbol <- EnsDb.Hsapiens.v86
keytypes(ref_db_gene_symbol)
ref_k <- keys(x = ref_db_gene_symbol,
              keytype = "SYMBOL")
gene2symbol <-
  AnnotationDbi::select(ref_db_gene_symbol,
                        keys = ref_k,
                        columns = "GENEID",
                        keytype = "SYMBOL")
gene2symbol <-
  gene2symbol[!duplicated(gene2symbol$SYMBOL), ]
gene2symbol <-
  gene2symbol[!duplicated(gene2symbol$GENEID), ]

### Get exp value from kallisto imports
df_4_DGE <-
  as.data.frame(df_raw$counts)

df_4_DGE$Geneid <-
  rownames(df_4_DGE)
df_4_DGE$Geneid <-
  str_split(string = df_4_DGE$Geneid,
            pattern = "\\.",
            simplify = T)[, 1]
df_4_DGE <-
  df_4_DGE[!duplicated(df_4_DGE$Geneid), ]

df_4_DGE <-
  merge(x = df_4_DGE,
        y = gene2symbol,
        by.x = "Geneid",
        by.y = "GENEID")

gene_names_list <- df_4_DGE$SYMBOL
df_4_DGE$Geneid <- NULL
df_4_DGE$SYMBOL <- NULL
rownames(df_4_DGE) <- gene_names_list

df_metadata <-
  read_excel("../PICALM_RNA_samples_description.xlsx")
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
  df_metadata[match(colnames(df_4_DGE),
                    df_metadata$Sample_Name), ]

df_metadata$Cell_line_ID <- as.factor(df_metadata$Cell_line_ID)
df_metadata$Clone_ID <- as.factor(df_metadata$Clone_ID)
df_metadata$biological_replicates <- as.factor(df_metadata$biological_replicates)
df_metadata$PICALM_isRisk <- as.factor(df_metadata$PICALM_isRisk)

df_metadata$genotype_edited[df_metadata$genotype_edited == "unedited"] <-
  "non_risk"
df_metadata$genotype_edited <- factor(df_metadata$genotype_edited,
                                      levels = c("non_risk",
                                                 "risk",
                                                 "CRISPRoff"))
df_metadata$rs10792832_genotype <- as.factor(df_metadata$rs10792832_genotype)

# df_4_DGE_backup <- df_4_DGE
# df_metadata_backup <- df_metadata

df_4_DGE <- df_4_DGE_backup
df_metadata <- df_metadata_backup

df_DGE <-
  DGEList(counts = df_4_DGE,
          samples = colnames(df_4_DGE),
          group = df_metadata$genotype_edited,
          genes = rownames(df_4_DGE),
          remove.zeros = T)

cpm_gene_count <-
  as.data.frame(cpm(df_DGE))
sum(rowMeans(cpm_gene_count) > 0.5)

cpm_cutoff <- 0.5 # 16301 genes
df_DGE <-
  df_DGE[rowSums(cpm_gene_count[, df_metadata$genotype_edited %in% "non_risk"] > cpm_cutoff) > 4 |
           rowSums(cpm_gene_count[, df_metadata$genotype_edited %in% "CRISPRoff"] > cpm_cutoff) > 4 |
           rowSums(cpm_gene_count[, df_metadata$genotype_edited %in% "risk"] > cpm_cutoff) > 4 , ]

df_DGE <- calcNormFactors(df_DGE)

# ! NOT REGRESSING OUT PICALM_isRisk) (rs10792832_genotype)
# ! for PCA plot only !
df_design_matrix <-
  model.matrix(~ 0 +
                 genotype_edited +
                 Cell_line_ID,
               data = df_metadata)


df_DGE <-
  estimateDisp(df_DGE,
               design = df_design_matrix,
               robust = T)
plotBCV(df_DGE)

result_QLM <-
  glmQLFit(df_DGE,
           design = df_design_matrix,
           robust = T)
plotQLDisp(result_QLM)

result_gene_list <-
  glmQLFTest(result_QLM,
             coef = 1:3)
summary(decideTestsDGE(result_gene_list))
result_gene_list <-
  topTags(result_gene_list,
          n = 2000)
result_gene_list <-
  result_gene_list$table$genes


df_PCA_to_plot <-
  prcomp(x = t(cpm_gene_count[rownames(cpm_gene_count) %in% result_gene_list, ]),
         center = T,
         scale. = T)


fviz_pca_ind(df_PCA_to_plot,
             repel = T,
             col.ind = df_metadata$genotype_edited,
             palette = "aaas",
             invisible = "quali",
             title = "PICALM_rs10792832_2000_genes_kallisto")


## calculate unedited vs risk
## note the changes in matrix.model
## disable combat_seq

df_4_DGE <- df_4_DGE_backup

# df_4_DGE <-
#   ComBat_seq(counts = as.matrix(df_4_DGE),
#              batch = df_metadata$Cell_line_ID,
#              group = df_metadata$genotype_edited)

df_DGE <-
  DGEList(counts = df_4_DGE,
          samples = colnames(df_4_DGE),
          group = df_metadata$genotype_edited,
          genes = rownames(df_4_DGE),
          remove.zeros = T)

cpm_gene_count <-
  as.data.frame(cpm(df_DGE))
sum(rowMeans(cpm_gene_count) > 0.5)

cpm_cutoff <- 0.5 # 16301 genes
df_DGE <-
  df_DGE[rowSums(cpm_gene_count[, df_metadata$genotype_edited %in% "unedited"] > cpm_cutoff) > 4 |
           rowSums(cpm_gene_count[, df_metadata$genotype_edited %in% "CRISPRoff"] > cpm_cutoff) > 4 |
           rowSums(cpm_gene_count[, df_metadata$genotype_edited %in% "risk"] > cpm_cutoff) > 4 , ]

df_DGE <- calcNormFactors(df_DGE)

# ! NOT REGRESSING OUT PICALM_isRisk) (rs10792832_genotype)


df_design_matrix <-
  model.matrix(~ genotype_edited +
                 Cell_line_ID,
               data = df_metadata)


df_DGE <-
  estimateDisp(df_DGE,
               design = df_design_matrix,
               robust = T)
plotBCV(df_DGE)

result_QLM <-
  glmQLFit(df_DGE,
           design = df_design_matrix,
           robust = T)
plotQLDisp(result_QLM)
result_gene_list <-
  glmQLFTest(result_QLM,
             coef = 2) # non_risk:risk:CRISPRoff
summary(decideTestsDGE(result_gene_list))
topTags(result_gene_list)

df_to_plot <- result_gene_list$table
df_to_plot$FDR <-
  p.adjust(df_to_plot$PValue,
           method = "fdr")

df_to_plot$y <-
  0 - log10(df_to_plot$PValue)

df_to_plot$gene_label <- NA
df_to_plot$gene_label[rownames(df_to_plot) %in% "PICALM"] <- "PICALM"

ggplot(data = df_to_plot,
       aes(x = logFC,
           y = y,
           label = gene_label,
           colour = ifelse(FDR < 0.05,
                           yes = "FDR < 0.05",
                           no = "FDR > 0.05"))) +
  geom_point(size = 0.5) +
  scale_colour_manual(values = c("red", "black")) +
  guides(colour = guide_legend(title = "")) +
  labs(title = paste("rs10792832, PICALM;\n",
                     "Risk vs non-risk, 11887 genes;\n",
                     "CPM > 0.5 in each group",
                     sep = "")) +
  xlab("logFC") +
  ylab("-log10Pvalue") +
  # xlim(c(0, 1)) +
  # ylim(c(0, 10)) +
  geom_text_repel(min.segment.length = 0,
                  force = 100,
                  colour = "red") +
  theme_classic() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 16))

## CRISPRoff vs non-risk
result_gene_list <-
  glmQLFTest(result_QLM,
             coef = 3) # non_risk:risk:CRISPRoff
summary(decideTestsDGE(result_gene_list))
topTags(result_gene_list)

df_to_plot <- result_gene_list$table
df_to_plot$FDR <-
  p.adjust(df_to_plot$PValue,
           method = "fdr")

df_to_plot$gene_label <- NA
df_to_plot$gene_label[rownames(df_to_plot) %in% "PICALM"] <- "PICALM"

df_to_plot$y <-
  0 - log10(df_to_plot$PValue)

ggplot(data = df_to_plot,
       aes(x = logFC,
           y = y,
           label = gene_label,
           colour = ifelse(FDR < 0.05,
                           yes = "FDR < 0.05",
                           no = "FDR > 0.05"))) +
  geom_point(size = 0.5) +
  scale_colour_manual(values = c("red", "black")) +
  guides(colour = guide_legend(title = "")) +
  labs(title = paste("rs10792832, PICALM;\n",
                     "CRISPRoff vs non-risk, 11887 genes;\n",
                     "CPM > 0.5 in each group",
                     sep = "")) +
  xlab("logFC") +
  ylab("-log10Pvalue") +
  geom_text_repel(min.segment.length = 0,
                  force = 100,
                  colour = "red") +
  # xlim(c(0, 1)) +
  # ylim(c(0, 10)) +
  theme_classic() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 16))


### Analyze with lines separated

## Line 04

df_4_DGE <- df_4_DGE_backup
df_metadata <- df_metadata_backup

df_4_DGE <- df_4_DGE[, df_metadata$Cell_line_ID %in% "CD04"]
df_metadata <- df_metadata[df_metadata$Cell_line_ID %in% "CD04", ]
# df_metadata$

df_DGE <-
  DGEList(counts = df_4_DGE,
          samples = colnames(df_4_DGE),
          group = df_metadata$genotype_edited,
          genes = rownames(df_4_DGE),
          remove.zeros = T)

cpm_gene_count <-
  as.data.frame(cpm(df_DGE))
sum(rowMeans(cpm_gene_count) > 0.5)

cpm_gene_count[rownames(cpm_gene_count) %in% "PICALM", ]

mean(unlist(cpm_gene_count[rownames(cpm_gene_count) %in% "PICALM", ])[1:3])

cpm_cutoff <- 0.5 # 16301 genes
df_DGE <-
  df_DGE[rowSums(cpm_gene_count[, df_metadata$genotype_edited %in% "unedited"] > cpm_cutoff) > 2 |
           rowSums(cpm_gene_count[, df_metadata$genotype_edited %in% "CRISPRoff"] > cpm_cutoff) > 2 |
           rowSums(cpm_gene_count[, df_metadata$genotype_edited %in% "risk"] > cpm_cutoff) > 2 , ]

df_DGE <- calcNormFactors(df_DGE)


df_design_matrix <-
  model.matrix(~ genotype_edited,
               data = df_metadata)

df_DGE <-
  estimateDisp(df_DGE,
               design = df_design_matrix,
               robust = T)
plotBCV(df_DGE)

result_QLM <-
  glmQLFit(df_DGE,
           design = df_design_matrix,
           robust = T)
plotQLDisp(result_QLM)

result_gene_list <-
  glmQLFTest(result_QLM,
             coef = 2) # non_risk:CRISPRoff:risk
summary(decideTestsDGE(result_gene_list))
topTags(result_gene_list)

df_to_plot <- result_gene_list$table
df_to_plot$FDR <-
  p.adjust(df_to_plot$PValue,
           method = "fdr")
# write.table(df_to_plot,
#             file = "PICALM_risk_vs_non_risk_21Mar2023.tsv",
#             quote = F, sep = "\t",
#             row.names = T, col.names = T)
# write.table(df_to_plot,
#             file = "PICALM_CD04_risk_vs_non_risk_28Mar2023.tsv",
#             quote = F, sep = "\t",
#             row.names = T, col.names = T)

df_to_plot$y <-
  0 - log10(df_to_plot$PValue)

df_to_plot$gene_label <- NA
df_to_plot$gene_label[rownames(df_to_plot) %in% "PICALM"] <- "PICALM"

ggplot(data = df_to_plot,
       aes(x = logFC,
           y = y,
           label = gene_label,
           colour = ifelse(FDR < 0.05,
                           yes = "FDR < 0.05",
                           no = "FDR > 0.05"))) +
  geom_point(size = 0.5) +
  scale_colour_manual(values = c("red", "black")) +
  guides(colour = guide_legend(title = "")) +
  labs(title = paste("rs10792832, PICALM;\n",
                     "CRISPRoff vs non-risk, 11611 genes;\n",
                     "CPM > 0.5 in each group",
                     sep = "")) +
  xlab("logFC") +
  ylab("-log10Pvalue") +
  # xlim(c(0, 1)) +
  # ylim(c(0, 10)) +
  geom_text_repel(min.segment.length = 0,
                  force = 100,
                  colour = "red") +
  theme_classic() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 16))


## Risk vs non-risk
result_gene_list <-
  glmQLFTest(result_QLM,
             coef = 3) # non_risk:CRISPRoff:risk
summary(decideTestsDGE(result_gene_list))
topTags(result_gene_list)

df_to_plot <- result_gene_list$table
df_to_plot$FDR <-
  p.adjust(df_to_plot$PValue,
           method = "fdr")


#
# write.table(df_to_plot,
#             file = "PICALM_CRISPRoff_vs_non_risk_21Mar2023.tsv",
#             quote = F, sep = "\t",
#             row.names = T, col.names = T)

# write.table(df_to_plot,
#             file = "PICALM_CD04_CRISPRoff_vs_non_risk_28Mar2023.tsv",
#             quote = F, sep = "\t",
#             row.names = T, col.names = T)
df_to_plot <- result_gene_list$table
df_to_plot$FDR <-
  p.adjust(df_to_plot$PValue,
           method = "fdr")

df_to_plot$gene_label <- NA
df_to_plot$gene_label[rownames(df_to_plot) %in% "PICALM"] <- "PICALM"

df_to_plot$y <-
  0 - log10(df_to_plot$PValue)

ggplot(data = df_to_plot,
       aes(x = logFC,
           y = y,
           label = gene_label,
           colour = ifelse(FDR < 0.05,
                           yes = "FDR < 0.05",
                           no = "FDR > 0.05"))) +
  geom_point(size = 0.5) +
  scale_colour_manual(values = c("red", "black")) +
  guides(colour = guide_legend(title = "")) +
  labs(title = paste("rs10792832, PICALM;\n",
                     "Risk vs non-risk, 11611 genes;\n",
                     "CPM > 0.5 in each group",
                     sep = "")) +
  xlab("logFC") +
  ylab("-log10Pvalue") +
  geom_text_repel(min.segment.length = 0,
                  force = 100,
                  colour = "red") +
  # xlim(c(0, 1)) +
  # ylim(c(0, 10)) +
  theme_classic() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 16))


## Line 09

df_4_DGE <- df_4_DGE_backup
df_metadata <- df_metadata_backup

df_4_DGE <- df_4_DGE[, df_metadata$Cell_line_ID %in% "CD09"]
df_metadata <- df_metadata[df_metadata$Cell_line_ID %in% "CD09", ]

df_DGE <-
  DGEList(counts = df_4_DGE,
          samples = colnames(df_4_DGE),
          group = df_metadata$genotype_edited,
          genes = rownames(df_4_DGE),
          remove.zeros = T)

cpm_gene_count <-
  as.data.frame(cpm(df_DGE))
sum(rowMeans(cpm_gene_count) > 0.5)

cpm_cutoff <- 0.5 # 16301 genes
df_DGE <-
  df_DGE[rowSums(cpm_gene_count[, df_metadata$genotype_edited %in% "unedited"] > cpm_cutoff) > 2 |
           rowSums(cpm_gene_count[, df_metadata$genotype_edited %in% "CRISPRoff"] > cpm_cutoff) > 2 |
           rowSums(cpm_gene_count[, df_metadata$genotype_edited %in% "risk"] > cpm_cutoff) > 2 , ]

df_DGE <- calcNormFactors(df_DGE)


df_design_matrix <-
  model.matrix(~ genotype_edited,
               data = df_metadata)

df_DGE <-
  estimateDisp(df_DGE,
               design = df_design_matrix,
               robust = T)
plotBCV(df_DGE)

result_QLM <-
  glmQLFit(df_DGE,
           design = df_design_matrix,
           robust = T)
plotQLDisp(result_QLM)

result_gene_list <-
  glmQLFTest(result_QLM,
             coef = 2) # non_risk:CRISPRoff:risk
summary(decideTestsDGE(result_gene_list))
topTags(result_gene_list)

df_to_plot <- result_gene_list$table
df_to_plot$FDR <-
  p.adjust(df_to_plot$PValue,
           method = "fdr")
# write.table(df_to_plot,
#             file = "PICALM_risk_vs_non_risk_21Mar2023.tsv",
#             quote = F, sep = "\t",
#             row.names = T, col.names = T)
# write.table(df_to_plot,
#             file = "PICALM_CD09_risk_vs_non_risk_28Mar2023.tsv",
#             quote = F, sep = "\t",
#             row.names = T, col.names = T)

df_to_plot$y <-
  0 - log10(df_to_plot$PValue)

df_to_plot$gene_label <- NA
df_to_plot$gene_label[rownames(df_to_plot) %in% "PICALM"] <- "PICALM"

ggplot(data = df_to_plot,
       aes(x = logFC,
           y = y,
           label = gene_label,
           colour = ifelse(FDR < 0.05,
                           yes = "FDR < 0.05",
                           no = "FDR > 0.05"))) +
  geom_point(size = 0.5) +
  scale_colour_manual(values = c("red", "black")) +
  guides(colour = guide_legend(title = "")) +
  labs(title = paste("rs10792832, PICALM;\n",
                     "Risk vs non-risk, 12257 genes;\n",
                     "CPM > 0.5 in each group",
                     sep = "")) +
  xlab("logFC") +
  ylab("-log10Pvalue") +
  # xlim(c(0, 1)) +
  # ylim(c(0, 10)) +
  geom_text_repel(min.segment.length = 0,
                  force = 100,
                  colour = "black") +
  theme_classic() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 16))


## Risk vs non-risk
result_gene_list <-
  glmQLFTest(result_QLM,
             coef = 3) # non_risk:CRISPRoff:non-risk
summary(decideTestsDGE(result_gene_list))
topTags(result_gene_list)

df_to_plot <- result_gene_list$table
df_to_plot$FDR <-
  p.adjust(df_to_plot$PValue,
           method = "fdr")


df_to_plot$y <-
  0 - log10(df_to_plot$PValue)

df_to_plot$gene_label <- NA
df_to_plot$gene_label[rownames(df_to_plot) %in% "PICALM"] <- "PICALM"

ggplot(data = df_to_plot,
       aes(x = logFC,
           y = y,
           label = gene_label,
           colour = ifelse(FDR < 0.05,
                           yes = "FDR < 0.05",
                           no = "FDR > 0.05"))) +
  geom_point(size = 0.5) +
  scale_colour_manual(values = c("red", "black")) +
  guides(colour = guide_legend(title = "")) +
  labs(title = paste("rs10792832, PICALM;\n",
                     "CRISPRoff vs non-risk, 12257 genes;\n",
                     "CPM > 0.5 in each group",
                     sep = "")) +
  xlab("logFC") +
  ylab("-log10Pvalue") +
  # xlim(c(0, 1)) +
  # ylim(c(0, 10)) +
  geom_text_repel(min.segment.length = 0,
                  force = 100,
                  colour = "red") +
  theme_classic() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 16))
