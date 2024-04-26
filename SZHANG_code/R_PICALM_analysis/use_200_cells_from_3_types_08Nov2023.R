# 08 Nov 2023 Siwei
# sample GABA, nmglut, and npglut neurons


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

# load data
load("~/backuped_space/Siwei_misc_R_projects/Human_cortex_100_lines_comparison/pre_merge_w_HE_ident.RData")

# load new 200 cells/each cell object, 08 Nov 2023
cell_3_types_w_GSE120046 <-
  readRDS("cell_3_types_from_76_lines_200_each_4_GSE120046.RDs")

cell_3_types_wo_GSE120046_200 <-
  cell_3_types_w_GSE120046[1:3]

# merge data, set metadata
cell_3_types_wo_GSE120046_200[["GABA"]]$subtype <- "GABA_iPS"
cell_3_types_wo_GSE120046_200[["nmglut"]]$subtype <- "nmglut_iPS"
cell_3_types_wo_GSE120046_200[["npglut"]]$subtype <- "npglut_iPS"

HE_ident <-
  pre_merge_GSE120046$orig.ident
HE_ident[str_detect(string = HE_ident,
                    pattern = "^HE.*",
                    negate = T)] <-
  NA
sum(is.na(HE_ident))
pre_merge_GSE120046$HE.ident.merged <-
  str_remove(string = HE_ident,
             pattern = "[0-9]$")
pre_merge_GSE120046$ref_identity <-
  "GSE120046"
cell_3_types_wo_GSE120046_200[["GABA"]]$ref_identity <- "GABA_iPS"
cell_3_types_wo_GSE120046_200[["nmglut"]]$ref_identity <- "nmglut_iPS"
cell_3_types_wo_GSE120046_200[["npglut"]]$ref_identity <- "npglut_iPS"

save(list = c("pre_merge_GSE120046",
              "cell_3_types_wo_GSE120046_200"),
     file = "premerge_raw_df_08Nov2023.RData")
load("premerge_raw_df_08Nov2023.RData")

# integration #####
Integrate_raw_df <-
  vector(mode = "list",
         length = 4L)
names(Integrate_raw_df) <-
  c("GABA", "nmglut", "npglut", "GSE120046")
Integrate_raw_df[1:3] <-
  cell_3_types_wo_GSE120046_200
Integrate_raw_df[["GSE120046"]] <-
  pre_merge_GSE120046

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
Integrate_df <-
  IntegrateData(anchorset = Integrate_df,
                dims = 1:30,
                normalization.method = "SCT",
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

## assign identities for plotting
Integrate_df_4_plot$cell.origin <-
  "GSE120046"
Integrate_df_4_plot$cell.origin[Integrate_df_4_plot$subtype %in%
                                  c("GABA_iPS")] <-
  "GABA_iPS"
Integrate_df_4_plot$cell.origin[Integrate_df_4_plot$subtype %in%
                                  c("nmglut_iPS")] <-
  "nmglut_iPS"
Integrate_df_4_plot$cell.origin[Integrate_df_4_plot$subtype %in%
                                  c("npglut_iPS")] <-
  "npglut_iPS"


# Integrate_df_4_plot$cell.origin[Integrate_df_4_plot$subtype %in%
#                                   c("GABA_iPS",
#                                     "nmglut_iPS",
#                                     "npglut_iPS")] <-
#   "scARC"

DimPlot(Integrate_df_4_plot,
        group.by = "cell.origin",
        label = T,
        repel = T)

DimPlot(Integrate_df_4_plot,
        group.by = "HE.ident.merged",
        label = T,
        repel = T)

DimPlot(Integrate_df_4_plot,
        group.by = "type",
        label = T,
        repel = T)

scARC_Velmeshev_2023 <-
  readRDS("~/Data/Databases/scARC_prenatal_postnatal_human_cort_dev_Velmeshev_2023/local.rds")

# test Velmeshev et al 2023 #####
# > unique(scARC_Velmeshev_2023$cell_type)
# [1] oligodendrocyte precursor cell native cell
# [3] neural cell                    oligodendrocyte
# [5] microglial cell                astrocyte
DimPlot(scARC_Velmeshev_2023,
        group.by = "cell_type",
        label = T,
        repel = T)
DimPlot(scARC_Velmeshev_2023,
        group.by = "development_stage",
        label = T,
        repel = T)
FeaturePlot(scARC_Velmeshev_2023,
            features = "ENSG00000128683")
# DefaultAssay(scARC_Velmeshev_2023) <- "SCT"
# scARC_Velmeshev_2023$


Idents(scARC_Velmeshev_2023) <-
  "cell_type"
neuron_subsetted_scARC_Velmeshev_2023 <-
  subset(x = scARC_Velmeshev_2023,
         idents = 'neural cell',
         downsample = 10000)
DimPlot(neuron_subsetted_scARC_Velmeshev_2023,
        group.by = "development_stage",
        label = F,
        repel = F)
# scARC_Velmeshev_2023$
# unique(neuron_subsetted_scARC_Velmeshev_2023$disease)
# neuron_subsetted_scARC_Velmeshev_2023 <-
#   subset(x = neuron_subsetted_scARC_Velmeshev_2023,
#          cells = str_detect(string = neuron_subsetted_scARC_Velmeshev_2023$development_stage,
#                             pattern = "LMP|month\\-old|newborn|1\\-year\\-old"))
neuron_subsetted_scARC_Velmeshev_2023 <-
  neuron_subsetted_scARC_Velmeshev_2023[ ,
                                        str_detect(string = neuron_subsetted_scARC_Velmeshev_2023$development_stage,
                                                           pattern = "LMP|month\\-old|newborn|1\\-year\\-old")]
neuron_subsetted_scARC_Velmeshev_2023 <-
  neuron_subsetted_scARC_Velmeshev_2023[ ,
                                         str_detect(string = neuron_subsetted_scARC_Velmeshev_2023$development_stage,
                                                    pattern = "21\\-year\\-old",
                                                    negate = T)]


DimPlot(neuron_subsetted_scARC_Velmeshev_2023,
        group.by = "development_stage",
        cols = c(brewer.pal(n = 8,
                            name = "Set1"),
                 brewer.pal(n = 8,
                            name = "Dark2")),
        label = F,
        repel = F)

saveRDS(neuron_subsetted_scARC_Velmeshev_2023,
        file = "scARC_Velmeshev_neuron_only_2023.RDs")

load("~/backuped_space/Siwei_misc_R_projects/Alena_RNASeq_23Aug2023/ENSG_gene_index.RData")
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
  cell_3_types_wo_GSE120046_200

##
dim(neuron_subsetted_scARC_Velmeshev_2023)
Integrate_raw_df[["Velmeshev_2023"]] <-
  neuron_subsetted_scARC_Velmeshev_2023
# dim(Integrate_raw_df[["Velmeshev_2023"]])
# length(rownames(Integrate_raw_df[["Velmeshev_2023"]]))
# length(gene_list_Velmeshev_2023$Geneid[!duplicated(gene_list_Velmeshev_2023$Gene_Symbol)])
# length(rownames(Integrate_raw_df[["Velmeshev_2023"]]) %in%
#           gene_list_Velmeshev_2023$Geneid[!duplicated(gene_list_Velmeshev_2023$Gene_Symbol)])
# sum(rownames(Integrate_raw_df[["Velmeshev_2023"]]) %in%
#          gene_list_Velmeshev_2023$Geneid[!duplicated(gene_list_Velmeshev_2023$Gene_Symbol)])
# dim(gene_list_Velmeshev_2023$Geneid[!duplicated(gene_list_Velmeshev_2023$Gene_Symbol), ])

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


# SCTransform(ScaleData(Integrate_raw_df[["Velmeshev_2023"]]))
# Integrate_raw_df[["Velmeshev_2023"]]@assays$RNA@data@Dimnames[[1]] <-
#   Gene_Symbol_Velmeshev_2023
# Integrate_raw_df[["Velmeshev_2023"]]@assays$SCT@counts@Dimnames[[1]] <-
#   Gene_Symbol_Velmeshev_2023
# Integrate_raw_df[["Velmeshev_2023"]]@assays$SCT@data@Dimnames[[1]] <-
#   Gene_Symbol_Velmeshev_2023

DefaultAssay(Integrate_raw_df[[1]]) <- "RNA"
DefaultAssay(Integrate_raw_df[[2]]) <- "RNA"
DefaultAssay(Integrate_raw_df[[3]]) <- "RNA"

## get a vector of shared gene lists across all 4 subsets
shared_gene_list <-
  intersect(rownames(Integrate_raw_df[[1]]),
            rownames(Integrate_raw_df[[2]]))
shared_gene_list <-
  intersect(shared_gene_list,
            rownames(Integrate_raw_df[[3]]))
shared_gene_list <-
  intersect(shared_gene_list,
            rownames(Integrate_raw_df[[4]]))


# Integrate_raw_df_backup <- Integrate_raw_df

Integrate_raw_df[[1]] <-
  Integrate_raw_df[[1]][rownames(Integrate_raw_df[[1]]) %in%
                          shared_gene_list, ]
Integrate_raw_df[[2]] <-
  Integrate_raw_df[[2]][rownames(Integrate_raw_df[[2]]) %in%
                          shared_gene_list, ]
Integrate_raw_df[[3]] <-
  Integrate_raw_df[[3]][rownames(Integrate_raw_df[[3]]) %in%
                          shared_gene_list, ]
Integrate_raw_df[[4]] <-
  Integrate_raw_df[[4]][rownames(Integrate_raw_df[[4]]) %in%
                          shared_gene_list, ]



Integrate_raw_df[[1]]@assays$integrated <- NULL
Integrate_raw_df[[1]]@assays$SCT <- NULL
Integrate_raw_df[[2]]@assays$integrated <- NULL
Integrate_raw_df[[2]]@assays$SCT <- NULL
Integrate_raw_df[[3]]@assays$integrated <- NULL
Integrate_raw_df[[3]]@assays$SCT <- NULL


Integrate_raw_df <-
  lapply(X = Integrate_raw_df,
         FUN = function(x) {
           x <- ScaleData(x)
           x <- SCTransform(x)
           x <- FindVariableFeatures(object = x,
                                     selection.method = "vst",
                                     nfeatures = 3000)
         })

rownames(Integrate_raw_df[["Velmeshev_2023"]])

# Integrate_raw_df_backup <- Integrate_raw_df
DefaultAssay(Integrate_raw_df[[1]])
DefaultAssay(Integrate_raw_df[[2]])
DefaultAssay(Integrate_raw_df[[3]])
DefaultAssay(Integrate_raw_df[[4]])

# Integrate_raw_df[["GABA"]]@assays$integrated <- NULL
# Integrate_raw_df[["nmglut"]]@assays$integrated <- NULL
# Integrate_raw_df[["npglut"]]@assays$integrated <- NULL
# Integrate_raw_df[["Velmeshev_2023"]]@assays$integrated <- NULL

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

# deal with Error: Attempting to add a different number of cells and/or features
rownames(Integrate_df@object.list[[1]])
rownames(Integrate_df@object.list[[2]])
rownames(Integrate_df@object.list[[3]])
rownames(Integrate_df@object.list[[4]])

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

## assign identities for plotting
Integrate_df_4_plot$cell.origin <-
  "Velmeshev_2023"
Integrate_df_4_plot$cell.origin[Integrate_df_4_plot$subtype %in%
                                  c("GABA_iPS")] <-
  "GABA_iPS"
Integrate_df_4_plot$cell.origin[Integrate_df_4_plot$subtype %in%
                                  c("nmglut_iPS")] <-
  "nmglut_iPS"
Integrate_df_4_plot$cell.origin[Integrate_df_4_plot$subtype %in%
                                  c("npglut_iPS")] <-
  "npglut_iPS"


###
