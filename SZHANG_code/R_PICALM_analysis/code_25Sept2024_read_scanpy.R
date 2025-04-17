# Siwei 25 Sept 2024
# clustering was done by scanpy, seems good

# init ####
{
  library(Seurat)
  library(Signac)

  # library(remotes)

  # library(SeuratDisk)
  library(anndata)

  library(edgeR)

  library(future)

  library(stringr)

  library(harmony)
  library(MAST)
  library(SingleCellExperiment)

  library(sctransform)
  library(glmGamPoi)

  library(future)
  library(dplyr)
  library(data.table)
  library(reshape2)

  library(readr)

  library(ggplot2)
  library(gplots)
  library(RColorBrewer)
  library(ggdark)
  library(viridis)

  library(lsr)

  library(factoextra)
}

plan("multisession", workers = 2)
options(mc.cores = 32)
set.seed(42)
options(future.globals.maxSize = 229496729600)

#####
# raw_h5ad_new <-
#   anndata::read_h5ad("AD_Homeostatis/ad_raw_annotated.h5ad")
# raw_h5ad_new
# names(raw_h5ad_new)
# raw_obs <- raw_h5ad_new$obs
#
#


MG_h5ad_scanpy_raw <-
  anndata::read_h5ad("AD_Homeostatis/ad_MG_raw_w_cluster.h5ad")
names(MG_h5ad_scanpy_raw)
MG_h5ad_scanpy_cpm <-
  anndata::read_h5ad("AD_Homeostatis/ad_MG_cpm_w_cluster.h5ad")

MG_sctransform_reference <-
  readRDS("AD_Homeostatis/GSE254205_sctransform_harmony_seurat.RDs")


MG_h5ad_scanpy_raw
MG_var <- MG_h5ad_scanpy_raw$var
MG_obs <- MG_h5ad_scanpy_raw$obs
MG_h5ad_scanpy_raw$obs_names
colnames(MG_obs)

df_metadata <-
  MG_obs
# MG_h5ad_scanpy_raw$obsm
# unique(df_metadata$leiden)
unique(df_metadata$cluster)

MG_h5ad_scanpy_cpm$n_vars
head(MG_h5ad_scanpy_cpm$var_names)

View(head(MG_h5ad_scanpy_raw$X))

# Assay_seurat_from_scanpy_df <-
#   CreateAssay5Object(counts = t(MG_h5ad_scanpy_raw$X),
#                      data = t(MG_h5ad_scanpy_cpm$X),
#                      project = "AD_GSE254205_MG_cpm_as_data)",
#                      meta.data = MG_h5ad_scanpy_raw$obs)
Assay_seurat_from_scanpy_df <-
  CreateAssay5Object(counts = t(MG_h5ad_scanpy_raw$X),
                     data = t(MG_h5ad_scanpy_cpm$X))
Seurat_from_scanpy_df <-
  CreateSeuratObject(counts = t(MG_h5ad_scanpy_raw$X),
                     # data = t(MG_h5ad_scanpy_cpm$X),
                     project = "AD_GSE254205_MG_cpm_as_data)",
                     meta.data = MG_h5ad_scanpy_raw$obs)
Seurat_from_scanpy_df@assays$RNA@layers$data <-
  Assay_seurat_from_scanpy_df$data

Seurat_from_scanpy_df$cluster
Idents(Seurat_from_scanpy_df) <- "cluster"
exp_by_cluster <-
  AggregateExpression(Seurat_from_scanpy_df,
                      assays = "RNA",
                      group.by = "cluster",
                      normalization.method = "RC",
                      scale.factor = 1e+6,
                      return.seurat = T,
                      verbose = T)

high_var_genes <-
  read_csv("AD_Homeostatis/high_var_gene_list.csv")
high_var_genes <-
  high_var_genes[high_var_genes$highly_variable == 1, ]

exp_by_cluster_hivar <-
  as.data.frame(exp_by_cluster@assays$RNA@layers$data)
rownames(exp_by_cluster_hivar) <-
  rownames(exp_by_cluster)
exp_by_cluster_hivar <-
  exp_by_cluster_hivar[rownames(exp_by_cluster_hivar) %in% high_var_genes$index, ]

# get gene exp from selected genes only
selected_gene_group <-
  c("P2RY12",
    "P2RY13",
    # "CX3CR1",
    "TMEM119",
    "CD9",
    "ITGAX",
    "CLEC7A",
    "CD63",
    "SPP1",
    "LPL",
    "TREM2",
    "APOE",
    "NAMPT",
    "ACSL1",
    "DPYD",
    "CD163",
    "PICALM")

(rownames(Seurat_from_scanpy_df) %in% c("P2RY12",
                                           "P2RY13",
                                           "CX3CR1",
                                           "TMEM119",
                                           "CD9",
                                           "ITGAX",
                                           "CLEC7A",
                                           "CD63",
                                           "SPP1",
                                           "LPL",
                                           "TREM2",
                                           "APOE",
                                           "NAMPT",
                                           "ACSL1",
                                           "DPYD",
                                           "CD163"))

c("P2RY12",
  "P2RY13",
  # "CX3CR1",
  "TMEM119",
  "CD9",
  "ITGAX",
  "CLEC7A",
  "CD63",
  "SPP1",
  "LPL",
  "TREM2",
  "APOE",
  "NAMPT",
  "ACSL1",
  "DPYD",
  "CD163") %in%
  rownames(Seurat_from_scanpy_df)

selected_exp_by_cluster <-
  AggregateExpression(Seurat_from_scanpy_df,
                      features = selected_gene_group,
                      assays = "RNA",
                      group.by = "cluster",
                      normalization.method = "RC",
                      scale.factor = 1e+6,
                      return.seurat = T,
                      verbose = T)
shared_cpm_scRNA <-
  as.matrix(selected_exp_by_cluster@assays$RNA@layers$data)
rownames(shared_cpm_scRNA) <-
  rownames(selected_exp_by_cluster)
  # unlist(selected_exp_by_cluster@assays$RNA@layers$data@Dim[[1]])
shared_cpm_scRNA <-
  shared_cpm_scRNA[order(rownames(shared_cpm_scRNA)), ]



# load Alena's data
Alena_metadata <-
  readxl::read_xlsx("../Alena_RNASeq_23Aug2023/STAR_aligned/Alena_microglia_Mar_metadata.xlsx",
             sheet = 1)
df_Alena_raw <-
  read.table("../Alena_RNASeq_23Aug2023/STAR_aligned/Alena_microglia_August_STAR.txt",
             header = T, sep = "\t")
load("../Alena_RNASeq_23Aug2023/ENSG_gene_index.RData")

df_4_DGE <- df_Alena_raw
df_4_DGE$Geneid <-
  str_split(string = df_4_DGE$Geneid,
            pattern = "\\.",
            simplify = T)[, 1]
df_4_DGE <-
  df_4_DGE[!duplicated(df_4_DGE$Geneid), ]

df_4_DGE <-
  merge(x = df_4_DGE,
        y = ENSG_anno_gene_indexed,
        by = "Geneid")
df_4_DGE <-
  df_4_DGE[!duplicated(df_4_DGE$Geneid), ]
df_4_DGE <-
  df_4_DGE[!duplicated(df_4_DGE$Gene_Symbol), ]
gene_names = df_4_DGE$Gene_Symbol
df_4_DGE$Geneid <- NULL
df_4_DGE$Gene_Symbol <- NULL

rownames(df_4_DGE) <-
  gene_names

Alena_metadata$off <-
  c(rep_len(x = "",
            length.out = 12),
    rep_len(x = "off",
            length.out = 6))
Alena_metadata$sample_col <-
  str_c(Alena_metadata$Clone_ID,
        Alena_metadata$off,
        "_",
        Alena_metadata$biological_replicates,
        "_")
Alena_metadata <-
  Alena_metadata[order(Alena_metadata$sample_col), ]

df_4_DGE <-
  df_4_DGE[, order(colnames(df_4_DGE))]
ncol(df_4_DGE)

Alena_metadata <-
  Alena_metadata[Alena_metadata$sample_col %in% colnames(df_4_DGE), ]
all(Alena_metadata$sample_col == colnames(df_4_DGE))

df_4_DGE_combat <-
  df_4_DGE[rowSums(df_4_DGE) > 0, ]
df_4_DGE_combat <-
  sva::ComBat_seq(as.matrix(df_4_DGE_combat),
                  batch = Alena_metadata$Cell_line_ID)

DGE_alena <-
  DGEList(counts = df_4_DGE_combat,
          genes = rownames(df_4_DGE_combat),
          samples = colnames(df_4_DGE_combat))

DGE_alena <-
  calcNormFactors(DGE_alena)
CPM_alena <-
  as.data.frame(edgeR::cpm(DGE_alena))
selected_gene_group %in%
  rownames(CPM_alena)

intersecting_gene_list <-
  high_var_genes$index[high_var_genes$index %in% rownames(CPM_alena)]

rownames(exp_by_cluster_hivar)

# df_intersect_hivar <-
#   exp_by_cluster_hivar[rownames(exp_by_cluster_hivar %in% intersecting_gene_list), ]
df_4_PCA <-
  cbind(exp_by_cluster_hivar[rownames(exp_by_cluster_hivar) %in% intersecting_gene_list, ],
        CPM_alena[rownames(CPM_alena) %in% intersecting_gene_list, ])

df_PCA_to_plot <-
  prcomp(x = t(log1p(df_4_PCA)),
         center = T,
         scale. = T)

fviz_pca_ind(df_PCA_to_plot,
             repel = T,
             col.ind = c(rep_len("human", length.out = 13),
                         Alena_metadata$rs10792832_genotype),
             palette = "aaas",
             invisible = "quali",
             title = "")


# shared_cpm_Alena <-
  # CPM_alena[rownames(CPM_alena) %in% selected_gene_group, ]
# shared_cpm_Alena <-
#   shared_cpm_Alena[order(rownames(shared_cpm_Alena)), ]
shared_cpm_Alena <-
  shared_cpm_Alena[match(x = selected_gene_group,
                            table = rownames(shared_cpm_Alena)) , ]
heatmap.2(x = as.matrix(shared_cpm_Alena[, c(1,2,3,12,13,14,
                                             6:11)]),
          scale = "row",
          col = "viridis",
          trace = "none",
          Rowv = F,
          Colv = F)

heatmap.2(x = as.matrix(shared_cpm_Alena[, c(1,2,3,12,14,
                                             6:7, 10:11)]),
          scale = "row",
          col = "viridis",
          Rowv = F,
          Colv = F)

cpm_risk_vs_nonrisk <-
  cbind(rowMeans(as.matrix(shared_cpm_Alena[, c(1,2,3)])),
        rowMeans(as.matrix(shared_cpm_Alena[, c(12,14)])),
        rowMeans(as.matrix(shared_cpm_Alena[, c(6:7)])),
        rowMeans(as.matrix(shared_cpm_Alena[, c(10:11)])))
colnames(cpm_risk_vs_nonrisk) <-
  c("non_risk_A1",
    "non_risk_G3",
    "risk_A3",
    "risk_A5")
# cpm_risk_vs_nonrisk <-
#   cpm_risk_vs_nonrisk[rownames(cpm_risk_vs_nonrisk) %in% selected_gene_group, ]
# cpm_risk_vs_nonrisk <-
#   cpm_risk_vs_nonrisk[selected_gene_group %in% rownames(cpm_risk_vs_nonrisk) , ]
cpm_risk_vs_nonrisk <-
  cpm_risk_vs_nonrisk[match(x = selected_gene_group,
                            table = rownames(cpm_risk_vs_nonrisk)) , ]
heatmap.2(x = as.matrix(cpm_risk_vs_nonrisk),
          scale = "row",
          col = "viridis",
          Rowv = F,
          Colv = F,
          cexCol = 1,
          trace = "none")

col_heads <-
  c(5, 1, 3, 2, 8, 0, 7, 6, 4)
heatmap.2(x = as.matrix(shared_cpm_scRNA[match(x = selected_gene_group,
                                               table = rownames(shared_cpm_scRNA)),
                                         col_heads + 1]),
          scale = "row",
          rowsep = c(3, 11, 15),
          colsep = c(5, 7),
          Colv = F,
          col = "viridis",
          labCol = col_heads,
          trace = "none",
          density.info = "none",
          Rowv = F)

DimPlot(Seurat_from_scanpy_df)

shared_cpm_scRNA_table <-
  shared_cpm_scRNA[match(x = selected_gene_group,
                         table = rownames(shared_cpm_scRNA)),
                   col_heads + 1]


