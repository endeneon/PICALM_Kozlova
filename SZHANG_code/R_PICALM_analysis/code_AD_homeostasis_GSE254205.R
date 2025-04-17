# Siwei 26 Jul 2024
# Import scRNA-seq data of GSE254025
# Check the expression of PICALM in risk vs non-risk cells
# also compare Alena's iMG snRNA-seq with GSE254025's HOMEO, LDAM, DAM populations
# run hierachical clustering

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
}

# install seuratDisk
# library(cli)
# library(crayon)
# library(hdf5r)
# library(Matrix)
# library(R6)
# library(rlang)
# library(Seurat)
# library(SeuratObject)
# library(stringi)
# library(withr)

# remotes::install_github("mojaveazure/seurat-disk")
# install.packages('anndata')

plan("multisession", workers = 6)
options(mc.cores = 32)
set.seed(42)
options(future.globals.maxSize = 229496729600)

# test loading of h5AD files from GSE254205 ####
# SeuratDisk::Convert(source = "AD_Homeostatis/GSE254205_ad_raw.h5ad",
#                     dest = "AD_Homeostatis/GSE_254205.h5seurat",
#                     verbose = T)
# GSE254205 <-
#   LoadH5Seurat(file = "AD_Homeostatis/GSE_254205.h5seurat")
raw_h5ad <-
  anndata::read_h5ad("AD_Homeostatis/GSE254205_ad_raw.h5ad")
# # View(raw_h5ad)
names(raw_h5ad)
# [1] ".__enclos_env__"         "raw"                     "uns"
# [4] "shape"                   "varp"                    "varm"
# [7] "var_names"               "var"                     "n_vars"
# [10] "obsp"                    "obsm"                    "obs_names"
# [13] "obs"                     "n_obs"                   "isbacked"
# [16] "is_view"                 "T"                       "layers"
# [19] "filename"                "X"                       ".get_py_object"
# [22] ".set_py_object"          "print"                   "write_loom"
# [25] "write_h5ad"              "write_csvs"              "transpose"
# [28] "to_df"                   "strings_to_categoricals" "rename_categories"
# [31] "copy"                    "concatenate"             "chunked_X"
# [34] "chunk_X"                 "uns_keys"                "varm_keys"
# [37] "var_names_make_unique"   "var_keys"                "obsm_keys"
# [40] "obs_names_make_unique"   "obs_keys"                "initialize"
raw_h5ad$obs # metadata
raw_h5ad$var_names # gene names
raw_h5ad$shape
# raw_df <-
#   CreateSeuratObject(counts = t(raw_h5ad$X),
#                      project = "AD_GSE254205",
#                      meta.data = raw_h5ad$obs)
# saveRDS(raw_df,
#         file = "AD_Homeostatis/Seurat_GSE254205.RDs")

seurat_raw <-
  readRDS(file = "AD_Homeostatis/Seurat_GSE254205.RDs")
rownames(seurat_raw)[str_detect(string = rownames(seurat_raw),
                                pattern = "^MT-")]

seurat_GSE254205_list <-
  SplitObject(seurat_raw,
              split.by = "sample")
seurat_GSE254205_list <-
  lapply(X = seurat_GSE254205_list,
         FUN = PercentageFeatureSet,
         pattern = "^MT-",
         col.name = "percent.mt")
seurat_GSE254205_list <-
  lapply(X = seurat_GSE254205_list,
         FUN = SCTransform,
         method = "glmGamPoi",
         return.only.var.genes = F,
         vars.to.regress = "percent.mt",
         verbose = T)
saveRDS(seurat_GSE254205_list,
        file = "AD_Homeostatis/GSE254205_sctransform_list.RDs")

seurat_GSE254205_list <-
  readRDS("AD_Homeostatis/GSE254205_sctransform_list.RDs")

var.features <-
  SelectIntegrationFeatures(seurat_GSE254205_list,
                            nfeatures = 3000,
                            fvf.nfeatures = 5000,
                            verbose = T)

seurat_GSE254205_sct <-
  merge(x = seurat_GSE254205_list[[1]],
        y = seurat_GSE254205_list[2:length(seurat_GSE254205_list)],
        merge.data = T)
VariableFeatures(seurat_GSE254205_sct) <-
  var.features

seurat_GSE254205_sct <-
  seurat_GSE254205_sct %>%
  RunPCA(verbose = T,
         seed.use = 42) %>%
  RunHarmony(assay.use = "SCT",
             group.by.vars = "sample",
             plot_convergence = T,
             reduction.save = "harmony",
             verbose = T) %>%
  RunUMAP(reduction = "harmony",
          dims = 1:50,
          verbose = T) %>%
  FindNeighbors(reduction = "harmony",
                dims = 1:50,
                verbose = T)

seurat_GSE254205_sct_clusters <-
  seurat_GSE254205_sct %>%
  FindClusters(resolution = 0.02,
               verbose = T)

DimPlot(seurat_GSE254205_sct_clusters,
        # cols = brewer.pal(n = 12,
        #                   name = "Set3"),
        label = T) #+
  dark_theme_classic()

rownames(seurat_GSE254205_sct_clusters)
FeaturePlot(seurat_GSE254205_sct_clusters,
            features = "CD74")

saveRDS(seurat_GSE254205_sct_clusters,
        file = "AD_Homeostatis/GSE254205_sctransform_harmony_seurat_10_clusters.RDs")
# seurat_GSE254205_sct_clusters <-
#   readRDS("AD_Homeostatis/GSE254205_sctransform_harmony_seurat_10_clusters.RDs")
# q(save = "no")

# isolate microglia (cluster 4) #####
microglia_GSE254205 <-
  subset(seurat_GSE254205_sct_clusters,
         subset = seurat_clusters == 4)

saveRDS(microglia_GSE254205,
        file = "AD_Homeostatis/GSE254205_sctransform_harmony_seurat_microglia_only.RDs")

microglia_GSE254205 <-
  microglia_GSE254205 %>%
  RunPCA(verbose = T,
         seed.use = 42) %>%
  RunHarmony(assay.use = "SCT",
             group.by.vars = "sample",
             plot_convergence = T,
             reduction.save = "harmony",
             verbose = T) %>%
  RunUMAP(reduction = "harmony",
          dims = 1:50,
          verbose = T) %>%
  FindNeighbors(reduction = "harmony",
                dims = 1:50,
                verbose = T)

microglia_GSE254205_clusters <-
  microglia_GSE254205 %>%
  FindClusters(resolution = 1,
               verbose = T)

microglia_GSE254205_clusters <-
  subset(microglia_GSE254205_clusters,
         subset = (seurat_clusters %in% 0:7))

microglia_GSE254205_clusters_final <-
  microglia_GSE254205_clusters %>%
  RunPCA(verbose = T,
         seed.use = 42) %>%
  RunHarmony(assay.use = "SCT",
             group.by.vars = "sample",
             plot_convergence = T,
             reduction.save = "harmony",
             verbose = T) %>%
  RunUMAP(reduction = "harmony",
          dims = 1:50,
          verbose = T) %>%
  FindNeighbors(reduction = "harmony",
                dims = 1:50,
                verbose = T)

microglia_GSE254205_clusters_final_10 <-
  microglia_GSE254205_clusters_final %>%
  FindClusters(resolution = 1.85,
               random.seed = 42,
               future.seed = T,
               method = 4,
               verbose = T)
DimPlot(microglia_GSE254205_clusters_final_10,
        # cols = brewer.pal(n = 12,
        #                   name = "Set3"),
        label = T) #+

# microglia_GSE254205_clusters_final_10 <-
#   subset(microglia_GSE254205_clusters_final_10,
#          subset = (seurat_clusters %in% 0:7))
# microglia_GSE254205_clusters_final_10 <-
#   microglia_GSE254205_clusters_final_10 %>%
#   FindClusters(resolution = 1.5,
#                verbose = T)
# microglia_GSE254205_clusters_final_10 <-
#   microglia_GSE254205_clusters_final_10 %>%
#   FindClusters(resolution = 2.0,
#                verbose = T)

#
saveRDS(microglia_GSE254205_clusters_final_10,
        file = "AD_Homeostatis/GSE254205_sctransform_harmony_seurat_microglia_only.RDs")
microglia_GSE254205_clusters_final_10 <-
  readRDS("AD_Homeostatis/GSE254205_sctransform_harmony_seurat_microglia_only.RDs")


microglia_GSE254205_gene_exp <-
  AggregateExpression(microglia_GSE254205_clusters_final_10,
                      assays = "SCT",
                      group.by = "seurat_clusters",
                      normalization.method = "RC")
microglia_GSE254205_gene_exp <-
  as.data.frame(microglia_GSE254205_gene_exp)

microglia_GSE254205_gene_exp$SCT.g10 <-
  (microglia_GSE254205_gene_exp$SCT.g10 +
     microglia_GSE254205_gene_exp$SCT.g11) / 2
microglia_GSE254205_gene_exp$SCT.g11 <- NULL

DEG_microglia_clusters <-
  DGEList(counts = as.matrix(microglia_GSE254205_gene_exp),
          genes = rownames(microglia_GSE254205_gene_exp),
          samples = colnames(microglia_GSE254205_gene_exp),
          remove.zeros = T)
DEG_microglia_clusters <-
  calcNormFactors(DEG_microglia_clusters)

cpm_microglia_clusters <-
  as.data.frame(edgeR::cpm(DEG_microglia_clusters))

microglia_GSE254205_gene_exp <-
  cpm_microglia_clusters[rownames(cpm_microglia_clusters) %in% c("P2RY12",
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
                                                                 "CD163"), ]

# microglia_GSE254205_clusters_final_10_lib_size <-
# heatmap_gene_list <-
#   rownames(microglia_GSE254205_gene_exp)
# microglia_GSE254205_gene_exp$SCT.g9 <- NULL

microglia_GSE254205_gene_exp_plot <-
  microglia_GSE254205_gene_exp

# microglia_GSE254205_gene_exp_plot <-
#   microglia_GSE254205_gene_exp[, (c(0,1,2,5,9,
#                                    3,6,
#                                    4,7,8,10) + 1)]
microglia_GSE254205_gene_exp_plot <-
  microglia_GSE254205_gene_exp_plot[c("P2RY12", "P2RY13", "CX3CR1", "TMEM119",
                                        "CD9", "ITGAX", "CLEC7A", "CD63", "SPP1",
                                        "LPL", "TREM2", "APOE",
                                        "NAMPT", "ACSL1", "DPYD", "CD163"), ]

heatmap.2(as.matrix(log(microglia_GSE254205_gene_exp_plot)),
          col = "viridis",
          scale = "row",
          trace = "none",
          Rowv = F,
          # Colv = F,
          sepcolor = "black", sepwidth = 0.5)

heatmap.2(as.matrix(log(microglia_GSE254205_gene_exp_plot[, c(7,0,5,6,1,
                                                              2,4,
                                                              3,8,9,
                                                              10) + 1])),
          col = "viridis",
          scale = "row",
          trace = "none",
          rowsep = c(4, 12),
          colsep = c(5, 7, 10),
          Rowv = F,
          Colv = F,
          sepcolor = "white",
          sepwidth = c(0.1, 0.1))
save.image(file = "Homeo_DAM_LDAM_MACRO_microglia.RData")

df_raw_Alena_iMG <-
  read.table("PICALM_risk_vs_non_risk_28Mar2023.tsv",
             header = T,
             sep = "\t")
df_raw_Alena_iMG_genes <-
  df_raw_Alena_iMG[rownames(df_raw_Alena_iMG) %in% c("P2RY12", "P2RY13", "CX3CR1", "TMEM119",
                                                     "CD9", "ITGAX", "CLEC7A", "CD63", "SPP1",
                                                     "LPL", "TREM2", "APOE",
                                                     "NAMPT", "ACSL1", "DPYD", "CD163"), ]
df_raw_Alena_iMG_genes$gene_symbol <-
  rownames(df_raw_Alena_iMG_genes)

df_CPM_log2_scRNAseq <-
  log2(microglia_GSE254205_gene_exp_plot[, c(7,0,5,6,1,
                                            2,4,
                                            3,8,9,
                                            10) + 1])
df_CPM_log2_scRNAseq <-
  df_CPM_log2_scRNAseq[rownames(df_CPM_log2_scRNAseq) %in% rownames(df_raw_Alena_iMG_genes), ]
df_CPM_log2_scRNAseq$gene_symbol <-
  rownames(df_CPM_log2_scRNAseq)

df_2_plot <-
  merge(x = df_CPM_log2_scRNAseq,
        y = df_raw_Alena_iMG_genes,
        by = "gene_symbol")
df_2_plot$FDR <- NULL
df_2_plot$PValue <- NULL
df_2_plot$F <- NULL
df_2_plot$logFC <- NULL

df_2_plot_gene_list <-
  df_2_plot$gene_symbol
df_2_plot$gene_symbol <- NULL
rownames(df_2_plot) <-
  df_2_plot_gene_list

df_2_plot <-
  df_2_plot[c("P2RY12", "P2RY13", "TMEM119",
              "CD9", "ITGAX", "CLEC7A", "CD63", "SPP1",
              "LPL", "TREM2", "APOE",
              "NAMPT", "ACSL1", "DPYD", "CD163"), ]

colnames(df_2_plot)[12] <- "iMG"
colnames(df_2_plot)[1:11] <-
  c("g0", "g1", "g2", "g5", "g9",
    "g3", "g6",
    "g4", "g7", "g8",
    "g10")

heatmap.2(as.matrix(df_2_plot),
          col = "viridis",
          scale = "row",
          trace = "none",
          rowsep = c(3, 11),
          colsep = c(5, 7, 10, 11),
          Rowv = F,
          Colv = F,
          tracecol = NULL,
          sepcolor = "white",
          sepwidth = c(0.1, 0.1))

save.image("plot_Alena_HOMEO_DAM_LDAM_MACRO_13Aug2024.RData")


# try to read the h5ad generated from the raw data #####
clustered_h5ad <-
  anndata::read_h5ad("AD_Homeostatis/ad_red_annotated.h5ad")

names(clustered_h5ad)
clustered_h5ad$obs # metadata
# raw_h5ad$var_names # gene names
# raw_h5ad$shape

Idents(microglia_GSE254205_clusters_final_10)

DimPlot(microglia_GSE254205_clusters_final_10,
        alpha = 0.8, pt.size = 1)
