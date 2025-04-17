# Siwei 20 Mar 2023
# Use tximport and EnsDB for gene name translation

# init
library(readr)
library(readxl)

library(stringr)

library(tximport)
library(EnsDb.Hsapiens.v86)
library(AnnotationDbi)

library(edgeR)
library(variancePartition)
library(factoextra)

library(sva)

library(ggplot2)
library(ggrepel)

# library(clusterProfiler)
# library(enrichplot)
# library(DOSE)

# library(org.Hs.eg.db)

# load files ####
db <- EnsDb.Hsapiens.v86

# db <- org.Hs.eg.db

columns(db)
# > columns(db)
# [1] "ENTREZID"            "EXONID"              "EXONIDX"             "EXONSEQEND"
# [5] "EXONSEQSTART"        "GENEBIOTYPE"         "GENEID"              "GENENAME"
# [9] "GENESEQEND"          "GENESEQSTART"        "INTERPROACCESSION"   "ISCIRCULAR"
# [13] "PROTDOMEND"          "PROTDOMSTART"        "PROTEINDOMAINID"     "PROTEINDOMAINSOURCE"
# [17] "PROTEINID"           "PROTEINSEQUENCE"     "SEQCOORDSYSTEM"      "SEQLENGTH"
# [21] "SEQNAME"             "SEQSTRAND"           "SYMBOL"              "TXBIOTYPE"
# [25] "TXCDSSEQEND"         "TXCDSSEQSTART"       "TXID"                "TXNAME"
# [29] "TXSEQEND"            "TXSEQSTART"          "UNIPROTDB"           "UNIPROTID"
# [33] "UNIPROTMAPPINGTYPE"

keytypes(db)
# > keytypes(db)
# [1] "ENTREZID"            "EXONID"              "GENEBIOTYPE"         "GENEID"
# [5] "GENENAME"            "PROTDOMID"           "PROTEINDOMAINID"     "PROTEINDOMAINSOURCE"
# [9] "PROTEINID"           "SEQNAME"             "SEQSTRAND"           "SYMBOL"
# [13] "TXBIOTYPE"           "TXID"                "TXNAME"              "UNIPROTID"

keys(db, keytype = "GENEID")

# gene2symbol <-
#   select(db,
#          keys = keys(db, keytype = "GENEID"),
#          columns = c("GENEID", "SYMBOL"),
#          keytype = "SYMBOL")

df_raw <-
  read_delim("../STAR_output/output/ReadsPerGene_STAR.txt",
             delim = "\t", escape_double = FALSE,
             trim_ws = TRUE)
colnames(df_raw) <- str_remove(colnames(df_raw), "_$")

df_4_DGE <- df_raw
df_4_DGE$Geneid <-
  str_split(string = df_4_DGE$Geneid,
            pattern = "\\.",
            simplify = T)[, 1]
df_4_DGE <-
  df_4_DGE[!duplicated(df_4_DGE$Geneid), ]

# keys require a list object as input
# Explicitely use AnnotationDbi::select (otherwise it confuses with
# clusterProfiler::select)
gene2symbol <-
  AnnotationDbi::select(db,
         keys = list(GeneIdFilter(df_4_DGE$Geneid)),
         columns = c("GENEID", "SYMBOL"))
gene2symbol <-
  gene2symbol[!duplicated(gene2symbol$SYMBOL), ]

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
df_metadata$genotype_edited <- factor(df_metadata$genotype_edited,
                                      levels = c("unedited",
                                                 "risk",
                                                 "CRISPRoff"))
df_metadata$rs10792832_genotype <- as.factor(df_metadata$rs10792832_genotype)


## for plotting PCA only !!
## try combat_seq from sva
# df_4_DGE_backup <- df_4_DGE
df_4_DGE <- df_4_DGE_backup

# df_4_DGE <-
#   ComBat_seq(counts = as.matrix(df_4_DGE),
#              batch = df_metadata$Cell_line_ID,
#              group = df_metadata$genotype_edited)

df_DGE <-
  DGEList(counts = as.matrix(df_4_DGE),
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

# ! NOT REGRESSING OUT PICALM_isRisk) (rs10792832_genotype)
# ! for PCA plot only !
df_design_matrix <-
  model.matrix(~ 0 +
                 genotype_edited +
                 Cell_line_ID +
                 RIN +
                 conc,
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
             title = "PICALM_rs10792832_2000_genes")


## calculate unedited vs risk
## note the changes in matrix.model
df_4_DGE <- df_4_DGE_backup

df_4_DGE <-
  ComBat_seq(counts = as.matrix(df_4_DGE),
             batch = df_metadata$Cell_line_ID,
             group = df_metadata$genotype_edited)

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

# ! NOT REGRESSING OUT PICALM_isRisk) (rs10792832_genotype)
# ! for PCA plot only !
df_design_matrix <-
  model.matrix(~ genotype_edited +
                 # Cell_line_ID +
                 RIN +
                 conc,
               data = df_metadata)
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
             coef = 2) # unedited:risk:CRISPRoff
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
                     "Risk vs non-risk, 16039 genes;\n",
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
             coef = 3) # unedited:risk:CRISPRoff
summary(decideTestsDGE(result_gene_list))
topTags(result_gene_list)

df_to_plot <- result_gene_list$table
df_to_plot$FDR <-
  p.adjust(df_to_plot$PValue,
           method = "fdr")

df_to_plot$gene_label <- NA
df_to_plot$gene_label[rownames(df_to_plot) %in% "PICALM"] <- "PICALM"
#
# write.table(df_to_plot,
#             file = "PICALM_CRISPRoff_vs_non_risk_21Mar2023.tsv",
#             quote = F, sep = "\t",
#             row.names = T, col.names = T)

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
                     "CRISPRoff vs non-risk, 16039 genes;\n",
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

## CRISPRoff vs risk
result_gene_list <-
  glmQLFTest(result_QLM,
             contrast = c(0, -1, 1, 0)) # unedited:risk:CRISPRoff
summary(decideTestsDGE(result_gene_list))
topTags(result_gene_list)

df_to_plot <- result_gene_list$table
df_to_plot$FDR <-
  p.adjust(df_to_plot$PValue,
           method = "fdr")
# write.table(df_to_plot,
#             file = "PICALM_CRISPRoff_vs_risk_21Mar2023.tsv",
#             quote = F, sep = "\t",
#             row.names = T, col.names = T)

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
                     "CRISPRoff vs risk, 16039 genes;\n",
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


### Run GSEA analysis with clusterProfileR

# CRISPRoff vs non-risk
result_gene_list <-
  glmQLFTest(result_QLM,
             coef = 3) # unedited:risk:CRISPRoff
summary(decideTestsDGE(result_gene_list))
topTags(result_gene_list)

df_to_plot <- result_gene_list$table
df_to_plot$FDR <-
  p.adjust(df_to_plot$PValue,
           method = "fdr")

input_gene_list <-
  df_to_plot$FDR
names(input_gene_list) <- rownames(df_to_plot)
input_gene_list <-
  na.omit(input_gene_list)
input_gene_list <-
  sort(input_gene_list,
       decreasing = T)

# check keytypes
keytypes(org.Hs.eg.db)
# > keytypes(org.Hs.eg.db)
# [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS"
# [6] "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"
# [11] "GENETYPE"     "GO"           "GOALL"        "IPI"          "MAP"
# [16] "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"
# [21] "PMID"         "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"
# [26] "UNIPROT"

gse_output <-
  gseGO(geneList = input_gene_list,
        ont = "BP",
        keyType = "SYMBOL",
        nPermSimple = 10000,
        minGSSize = 5,
        maxGSSize = 1000,
        pvalueCutoff = 0.05,
        pAdjustMethod = "none",
        seed = 42,
        OrgDb = "org.Hs.eg.db",
        verbose = T)
dotplot(gse_output,
        showCategory = 10,
        split = ".sign") +
  facet_grid(. ~ .sign) +
  theme(axis.text.y.left = element_text(size = 8))


# risk vs non-risk
result_gene_list <-
  glmQLFTest(result_QLM,
             coef = 2) # unedited:risk:CRISPRoff
summary(decideTestsDGE(result_gene_list))
topTags(result_gene_list)

df_to_plot <- result_gene_list$table
df_to_plot$FDR <-
  p.adjust(df_to_plot$PValue,
           method = "fdr")

input_gene_list <-
  df_to_plot$FDR
names(input_gene_list) <- rownames(df_to_plot)
input_gene_list <-
  na.omit(input_gene_list)
input_gene_list <-
  sort(input_gene_list,
       decreasing = T)

gse_output <-
  gseGO(geneList = input_gene_list,
        ont = "BP",
        keyType = "SYMBOL",
        nPermSimple = 10000,
        minGSSize = 5,
        maxGSSize = 1000,
        pvalueCutoff = 0.05,
        pAdjustMethod = "none",
        seed = 42,
        OrgDb = "org.Hs.eg.db",
        verbose = T)

dotplot(gse_output,
        showCategory = 10,
        split = ".sign") +
  facet_grid(. ~ .sign) +
  theme(axis.text.y.left = element_text(size = 8))

# CRISPRoff vs risk
result_gene_list <-
  glmQLFTest(result_QLM,
             contrast = c(0, -1, 1, 0)) # unedited:risk:CRISPRoff
summary(decideTestsDGE(result_gene_list))
topTags(result_gene_list)

df_to_plot <- result_gene_list$table
df_to_plot$FDR <-
  p.adjust(df_to_plot$PValue,
           method = "fdr")

input_gene_list <-
  df_to_plot$FDR
names(input_gene_list) <- rownames(df_to_plot)
input_gene_list <-
  na.omit(input_gene_list)
input_gene_list <-
  sort(input_gene_list,
       decreasing = T)

gse_output <-
  gseGO(geneList = input_gene_list,
        ont = "BP",
        keyType = "SYMBOL",
        nPermSimple = 10000,
        minGSSize = 5,
        maxGSSize = 1000,
        pvalueCutoff = 0.05,
        pAdjustMethod = "none",
        seed = 42,
        OrgDb = "org.Hs.eg.db",
        verbose = T)

dotplot(gse_output,
        showCategory = 10,
        split = ".sign") +
  facet_grid(. ~ .sign) +
  theme(axis.text.y.left = element_text(size = 8))
