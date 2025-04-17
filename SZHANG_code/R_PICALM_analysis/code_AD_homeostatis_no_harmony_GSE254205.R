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

microglia_byannotation_GSE254205 <-
  microglia_byannotation_GSE254205 %>%
  NormalizeData(normalization.method = "LogNormalize",
                scale.factor = 10000,
                verbose = T) %>%
  FindVariableFeatures(selection.method = "vst",
                       nfeatures = 2000,
                       verbose = T) %>%
  ScaleData(features = rownames(microglia_byannotation_GSE254205),
            vars.to.regress = c("pct_counts_mt",
                                "pct_counts_rb",
                                "orig.ident"))

microglia_byannotation_GSE254205 <-
  microglia_byannotation_GSE254205 %>%
  RunPCA(verbose = T,
         seed.use = 42,
         features = VariableFeatures(microglia_byannotation_GSE254205),
         npcs = 50) %>%
  RunUMAP(reduction = "pca",
          dims = 1:20,
          verbose = T) %>%
  FindNeighbors(reduction = "pca",
                dims = 1:20,
                verbose = T)

sct_microglia_w_cluster <-
  microglia_byannotation_GSE254205 %>%
  FindClusters(resolution = 0.16,
               random.seed = 42,
               algorithm = 4,
               verbose = T)
DimPlot(sct_microglia_w_cluster,
        label = T)

sct_microglia_w_cluster$uuid <-
  names(sct_microglia_w_cluster$orig.ident)
length(sct_microglia_w_cluster$orig.ident) # 6659

HOMEO_avgExp <-
  AverageExpression(sct_microglia_w_cluster,
                    features = c("P2RY12",
                                 "P2RY13",
                                 "CX3CR1",
                                 "TMEM119"),
                    group.by = "uuid",
                    layer = "scale.data",
                    verbose = T)
HOMEO_avgExp <-
  as.data.frame(as.matrix(HOMEO_avgExp$RNA))

HOMEO_sig_scores <-
  colMeans2(as.matrix(HOMEO_avgExp))
names(HOMEO_sig_scores) <- NULL
sct_microglia_w_cluster$HOMEO_sig_scores <-
  HOMEO_sig_scores

FeaturePlot(sct_microglia_w_cluster,
            features = "HOMEO_sig_scores",
            cols = viridis(n = 50))

LDAM_avgExp <-
  AverageExpression(sct_microglia_w_cluster,
                    features = c("NAMPT",
                                 "ACSL1",
                                 "DPYD",
                                 "CD163"),
                    group.by = "uuid",
                    layer = "scale.data",
                    verbose = T)
LDAM_avgExp <-
  as.data.frame(as.matrix(LDAM_avgExp$RNA))

LDAM_sig_scores <-
  colMeans2(as.matrix(LDAM_avgExp))
names(LDAM_sig_scores) <- NULL
sct_microglia_w_cluster$LDAM_sig_scores <-
  LDAM_sig_scores

# FeaturePlot(sct_microglia_w_cluster,
#             features = "LDAM_sig_scores",
#             cols = c("white",
#                      "darkblue"),
#             alpha = 0.7)

FeaturePlot(sct_microglia_w_cluster,
            features = "LDAM_sig_scores",
            cols = viridis(n = 50),
            alpha = 0.7)
