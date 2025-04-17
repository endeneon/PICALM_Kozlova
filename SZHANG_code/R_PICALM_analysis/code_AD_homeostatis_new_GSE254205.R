# Siwei 19 Sept 2024
# Import scRNA-seq data of GSE254025
# found another source of raw_h5ad, claimed to be annotated
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

plan("multisession", workers = 2)
options(mc.cores = 32)
set.seed(42)
options(future.globals.maxSize = 229496729600)

raw_h5ad_new <-
  anndata::read_h5ad("AD_Homeostatis/ad_raw_annotated.h5ad")
names(raw_h5ad_new)

df_metadata_raw <-
  raw_h5ad_new$obs
unique(df_metadata_raw$cell_type)

raw_df <-
  CreateSeuratObject(counts = t(raw_h5ad_new$X),
                     project = "AD_GSE254205_new_source)",
                     meta.data = raw_h5ad_new$obs)
microglia_byannotation_GSE254205 <-
  subset(raw_df,
         subset = cell_type == "Microglia")

microglia_byannotation_GSE254205_list <-
  SplitObject(microglia_byannotation_GSE254205,
              split.by = "sample")
microglia_byannotation_GSE254205_list <-
  lapply(X = microglia_byannotation_GSE254205_list,
         FUN = SCTransform,
         method = "glmGamPoi",
         return.only.var.genes = F,
         vars.to.regress = c("pct_counts_mt",
                             "pct_counts_rb"),
         verbose = T)

var.features <-
  SelectIntegrationFeatures(microglia_byannotation_GSE254205_list,
                            nfeatures = 2000,
                            fvf.nfeatures = 5000,
                            verbose = T)

sct_microglia_byannotation_GSE254205 <-
  merge(x = microglia_byannotation_GSE254205_list[[1]],
        y = microglia_byannotation_GSE254205_list[2:length(microglia_byannotation_GSE254205_list)],
        merge.data = T)
VariableFeatures(sct_microglia_byannotation_GSE254205) <-
  var.features

sct_microglia_byannotation_GSE254205 <-
  sct_microglia_byannotation_GSE254205 %>%
  RunPCA(verbose = T,
         seed.use = 42,
         features = var.features,
         npcs = 50) %>%
  RunHarmony(assay.use = "SCT",
             group.by.vars = "sample",
             plot_convergence = T,
             reduction.save = "harmony",
             verbose = T) %>%
  RunUMAP(reduction = "harmony",
          dims = 1:20,
          verbose = T) %>%
  FindNeighbors(reduction = "harmony",
                dims = 1:20,
                verbose = T)

sct_microglia_clusters <-
  sct_microglia_byannotation_GSE254205 %>%
  FindClusters(resolution = 0.45,
               random.seed = 42,
               algorithm = 4,
               verbose = T)

DimPlot(sct_microglia_clusters,
        label = T)
sct_microglia_clusters$state <-
  "HOMEO"
sct_microglia_clusters$state[sct_microglia_clusters$seurat_clusters %in% c(4, 9)] <-
  "DAM"
sct_microglia_clusters$state[sct_microglia_clusters$seurat_clusters %in% c(5, 8, 10)] <-
  "LDAM"

Idents(sct_microglia_clusters) <- "state"
DimPlot(sct_microglia_clusters,
        label = T)



DoHeatmap(sct_microglia_clusters,
          features = c("P2RY12",
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
                       "CD163"),
          group.by = "seurat_clusters")

sct_microglia_clusters$uuid <-
  names(sct_microglia_clusters$orig.ident)

# Homeo_sig_scores <-
#   AggregateExpression(sct_microglia_clusters,
#                       assays = "RNA",
#                       features = c("P2RY12",
#                                    "P2RY13",
#                                    "CX3CR1"),
#                       group.by = "uuid",
#                       normalization.method = "LogNormalize",
#                       scale.factor = 1e6,
#                       return.seurat = T)
# Homeo_sig_scores <-
#   as.data.frame(Homeo_sig_scores@assays$RNA$data)
# sct_microglia_clusters$Homeo_sig_scores_2 <-
#   10^(0 - (rowSums(Homeo_sig_scores, na.rm = T)))
DefaultAssay(sct_microglia_clusters) <- "RNA"
Homeo_exp_genes <-
  as.data.frame(as.matrix(microglia_byannotation_GSE254205@assays$RNA$counts[rownames(microglia_byannotation_GSE254205) %in%
                                                   c("P2RY12",
                                                     "P2RY13",
                                                     "CX3CR1",
                                                     "TMEM119"), ]))
gross_counts_per_cell <-
  as.data.frame(as.matrix(colSums2(microglia_byannotation_GSE254205@assays$RNA$counts)))
Homeo_sig_scores <-
  colMeans2(as.matrix(Homeo_exp_genes)) / gross_counts_per_cell$V1
sct_microglia_clusters$Homeo_sig_scores <-
  Homeo_sig_scores

FeaturePlot(sct_microglia_clusters,
            features = "Homeo_sig_scores",
            cols = c("white",
                     "darkblue"),
            alpha = 0.7)

DAM_exp_genes <-
  as.data.frame(as.matrix(microglia_byannotation_GSE254205@assays$RNA$counts[rownames(microglia_byannotation_GSE254205) %in%
                                                                               c("CD9",
                                                                                 "ITGAX",
                                                                                 "CLEC7A",
                                                                                 "CD63",
                                                                                 "SPP1",
                                                                                 "LPL",
                                                                                 "TREM2",
                                                                                 "APOE"), ]))
DAM_sig_scores <-
  colMeans2(as.matrix(DAM_exp_genes)) / gross_counts_per_cell$V1
sct_microglia_clusters$DAM_sig_scores <-
  DAM_sig_scores

FeaturePlot(sct_microglia_clusters,
            features = "DAM_sig_scores",
            cols = c("white",
                     "darkblue"),
            alpha = 0.7)

LDAM_exp_genes <-
  as.data.frame(as.matrix(microglia_byannotation_GSE254205@assays$RNA$counts[rownames(microglia_byannotation_GSE254205) %in%
                                                                               c("NAMPT",
                                                                                 "ACSL1",
                                                                                 "DPYD",
                                                                                 "CD163"), ]))
LDAM_sig_scores <-
  colMeans2(as.matrix(LDAM_exp_genes)) / gross_counts_per_cell$V1
sct_microglia_clusters$LDAM_sig_scores <-
  LDAM_sig_scores

FeaturePlot(sct_microglia_clusters,
            features = "LDAM_sig_scores",
            cols = c("white",
                     "darkblue"),
            alpha = 0.7)




seurat_2_plot <-
  AggregateExpression(sct_microglia_clusters,
                      assays = "SCT",
                      features = c("P2RY12",
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
                                   "CD163"),
                      group.by = "seurat_clusters",
                      normalization.method = "RC",
                      scale.factor = 1e6,
                      return.seurat = T)
df_2_plot <-
  seurat_2_plot@assays$SCT@layers$scale.data
df_2_plot <-
  as.data.frame(df_2_plot)
rownames(df_2_plot) <-
  seurat_2_plot@assays$SCT@features[[1]]
colnames(df_2_plot) <-
  seurat_2_plot@assays$SCT@cells[[1]]
# df_2_plot$Gene <-
#   rownames(df_2_plot)

df_4_heatmap <-
  melt(df_2_plot)


heatmap.2(as.matrix(df_2_plot[, c(1,6,3,2,7,
                                  9,4,
                                  8,10,5)]),
          col = "viridis",
          scale = "row",
          trace = "none",
          Rowv = F,
          Colv = F,
          colsep = c(5, 7, 10),
          sepcolor = "white",
          sepwidth = c(0.1, 0.1))
# seurat_2_plot@assays$SCT@features[[1]]
#
# saveRDS(raw_df,
#         file = "AD_Homeostatis/Seurat_GSE254205.RDs")

