# 09 Nov 2023 Siwei
# sample GABA, nmglut, and npglut neurons
# Use 2000 cells per type each

# init ####
{
  library(Seurat)
  library(Signac)
  library(readr)
  library(future)
  library(parallel)
  library(ggplot2)
  library(RColorBrewer)
  library(stringr)
}

plan("multisession", workers = 2)
options(expressions = 20000)
options(future.globals.maxSize = 207374182400)
options(future.seed = T)
set.seed(42)

# load data #####
load("~/backuped_space/Siwei_misc_R_projects/Alena_RNASeq_23Aug2023/ENSG_gene_index.RData")
# Velmeshev 2023
neuron_subsetted_scARC_Velmeshev_2023 <-
  readRDS("scARC_Velmeshev_neuron_only_2023.RDs")
# 3 cell types of 2000 each
cell_3_types_2K_each <-
  readRDS("cell_3_types_from_76_lines_2K_each_4_Velmeshev.RDs")

## preprocess 3 cell types #####
cell_3_types_2K_each[["GABA"]]$subtype <- "GABA_iPS"
cell_3_types_2K_each[["nmglut"]]$subtype <- "nmglut_iPS"
cell_3_types_2K_each[["npglut"]]$subtype <- "npglut_iPS"


## preprocess Velmeshev #####
gene_list_Velmeshev_2023 <-
  merge(x = data.frame(Geneid = rownames(neuron_subsetted_scARC_Velmeshev_2023)),
        y = ENSG_anno_gene_indexed,
        by = "Geneid")

## use Velmeshev 2023 as the reference #####
# integration
DefaultAssay(neuron_subsetted_scARC_Velmeshev_2023)
Integrate_raw_df <-
  vector(mode = "list",
         length = 4L)
names(Integrate_raw_df) <-
  c("GABA", "nmglut", "npglut", "Velmeshev_2023")
Integrate_raw_df[1:3] <-
  cell_3_types_2K_each

dim(neuron_subsetted_scARC_Velmeshev_2023)
Integrate_raw_df[["Velmeshev_2023"]] <-
  neuron_subsetted_scARC_Velmeshev_2023

Integrate_raw_df[["Velmeshev_2023"]] <-
  Integrate_raw_df[["Velmeshev_2023"]][rownames(Integrate_raw_df[["Velmeshev_2023"]]) %in%
                                         gene_list_Velmeshev_2023$Geneid[!duplicated(gene_list_Velmeshev_2023$Gene_Symbol)], ]
dim(Integrate_raw_df[["Velmeshev_2023"]])

gene_list_Velmeshev_2023 <-
  gene_list_Velmeshev_2023[!duplicated(gene_list_Velmeshev_2023$Gene_Symbol), ]
Gene_Symbol_Velmeshev_2023 <-
  gene_list_Velmeshev_2023$Gene_Symbol[rownames(Integrate_raw_df[["Velmeshev_2023"]]) %in%
                                         gene_list_Velmeshev_2023$Geneid[!duplicated(gene_list_Velmeshev_2023$Gene_Symbol)]]

Integrate_raw_df[["Velmeshev_2023"]]@assays$RNA@counts@Dimnames[[1]] <-
  Gene_Symbol_Velmeshev_2023
Integrate_raw_df[["Velmeshev_2023"]]@assays$RNA@data@Dimnames[[1]] <-
  Gene_Symbol_Velmeshev_2023
rownames(Integrate_raw_df[["Velmeshev_2023"]])
dim(Integrate_raw_df[["Velmeshev_2023"]])

DefaultAssay(Integrate_raw_df[[1]]) <- "RNA"
DefaultAssay(Integrate_raw_df[[2]]) <- "RNA"
DefaultAssay(Integrate_raw_df[[3]]) <- "RNA"

## get a vector of shared gene lists across all 4 subsets #####
shared_gene_list <-
  intersect(rownames(Integrate_raw_df[[1]]),
            rownames(Integrate_raw_df[[2]]))
shared_gene_list <-
  intersect(shared_gene_list,
            rownames(Integrate_raw_df[[3]]))
shared_gene_list <-
  intersect(shared_gene_list,
            rownames(Integrate_raw_df[[4]]))

## SCT on each of the list items ####
Integrate_raw_df <-
  lapply(X = Integrate_raw_df,
         FUN = function(x) {
           x <- ScaleData(x)
           x <- SCTransform(x)
           x <- FindVariableFeatures(object = x,
                                     selection.method = "vst",
                                     nfeatures = 3000)
         })

# The expected workflow for integratinge assays produced by SCTransform is
# SelectIntegrationFeatures -> PrepSCTIntegration -> FindIntegrationAnchors.
anchors_4_integration <-
  SelectIntegrationFeatures(object.list = Integrate_raw_df,
                            nfeatures = 3000,
                            verbose = T)
Integrate_df <-
  PrepSCTIntegration(object.list = Integrate_raw_df,
                     anchor.features = anchors_4_integration,
                     verbose = T)
Integrate_df <-
  FindIntegrationAnchors(object.list = Integrate_df,
                         dims = 1:30,
                         anchor.features = anchors_4_integration,
                         normalization.method = "SCT")

## run integration
Integrate_df <-
  IntegrateData(anchorset = Integrate_df,
                # dims = 1:30,
                # normalization.method = "SCT",
                verbose = T)

## pre-process data on the integrated df
Integrate_df_4_plot <-
  ScaleData(Integrate_df,
            verbose = T)
Integrate_df_4_plot <-
  RunPCA(Integrate_df_4_plot,
         verbose = T)
Integrate_df_4_plot <-
  RunUMAP(Integrate_df_4_plot,
          dims = 1:30,
          verbose = T)

# get clusters
Integrate_df_4_plot <-
  FindNeighbors(Integrate_df_4_plot,
                dims = 1:30,
                verbose = T)
Integrate_df_4_plot <-
  FindClusters(Integrate_df_4_plot,
               # dims = 1:30,
               verbose = T)

DimPlot(Integrate_df_4_plot,
        group.by = "subtype",
        label = T,
        repel = T)

Integrate_df_4_plot$development_stage[Integrate_df_4_plot$subtype %in%
                                        "GABA_iPS"] <-
  "GABA_iPS"
Integrate_df_4_plot$development_stage[Integrate_df_4_plot$subtype %in%
                                        "nmglut_iPS"] <-
  "nmglut_iPS"
Integrate_df_4_plot$development_stage[Integrate_df_4_plot$subtype %in%
                                        "npglut_iPS"] <-
  "npglut_iPS"

DimPlot(Integrate_df_4_plot,
        group.by = "development_stage",
        cols = c(rev(brewer.pal(n = 10,
                                name = "Set3")),
                 brewer.pal(n = 8,
                            name = "Dark2")),
        label = T,
        repel = T)
FeaturePlot(Integrate_df_4_plot,
        features = "SLC17A6")
FeaturePlot(Integrate_df_4_plot,
            features = "SLC17A7")
FeaturePlot(Integrate_df_4_plot,
            features = "GAD2")
FeaturePlot(Integrate_df_4_plot,
            features = "GAD1")
FeaturePlot(Integrate_df_4_plot,
            features = "MAP2")
FeaturePlot(Integrate_df_4_plot,
            features = "VIM")
FeaturePlot(Integrate_df_4_plot,
            features = "SOX2")
FeaturePlot(Integrate_df_4_plot,
            features = "NANOG")
DefaultAssay(Integrate_df_4_plot) <- "SCT"

save.image("Use_2000_cells_Vesmeshev_2023.RData")
# Integrate_df_4_plot$

unique(Integrate_df_4_plot$cell_type)
