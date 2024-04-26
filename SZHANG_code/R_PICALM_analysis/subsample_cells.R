# 06 Nov 2023 Siwei
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

# load data #####
load("~/NVME/scARC_Duan_018/018-029_combined_analysis/018-029_RNA_integrated_labeled_with_harmony.RData")

# subset cells from Lexi's integrated_labeled object
Assays(integrated_labeled)
DefaultAssay(integrated_labeled) <- "RNA"
Idents(integrated_labeled) <- "cell.type"

raw_76_lines_list <-
  vector(mode = "list",
         length = 3L)
names(raw_76_lines_list) <-
  c("GABA", "nmglut", "npglut")

downsampled_integrated <-
  subset(x = integrated_labeled,
         downsample = 100)

raw_76_lines_list[["GABA"]] <-
  subset(x = downsampled_integrated,
         idents = 'GABA')
raw_76_lines_list[["nmglut"]] <-
  subset(x = downsampled_integrated,
         idents = 'nmglut')
raw_76_lines_list[["npglut"]] <-
  subset(x = downsampled_integrated,
         idents = 'npglut')

saveRDS(raw_76_lines_list,
        file = "raw_76_lines_list_100_cell_each.RDS")
rm(integrated_labeled)
rm(downsampled_integrated)

# load premerge dataset
pre_merge_GSE120046 <-
  readRDS(file = "pre_merge_GSE120046_cleaned.RDs")
unique(pre_merge_GSE120046$subtype)

# merge
raw_76_lines_list[["GABA"]]$subtype <- "GABA_iPS"
raw_76_lines_list[["nmglut"]]$subtype <- "nmglut_iPS"
raw_76_lines_list[["npglut"]]$subtype <- "npglut_iPS"

merged_GSE120046_76_lines <-
  merge(x = pre_merge_GSE120046,
        y = raw_76_lines_list[1:3],
        merge.data = T)
VariableFeatures_3000 <-
  FindVariableFeatures(object = pre_merge_GSE120046,
                       selection.method = "vst",
                       nfeatures = 3000)

merged_GSE120046_76_lines <-
  ScaleData(merged_GSE120046_76_lines)
merged_GSE120046_76_lines <-
  SCTransform(merged_GSE120046_76_lines,
              residual.features = VariableFeatures_3000@assays$RNA@var.features)

merged_GSE120046_76_lines <-
  RunPCA(merged_GSE120046_76_lines,
         verbose = T)
# merged_GSE120046_76_lines <-
#   RunTSNE(merged_GSE120046_76_lines,
#           reduction = "pca")
merged_GSE120046_76_lines <-
  RunUMAP(merged_GSE120046_76_lines,
          dims = 1:30,
          verbose = T)

merged_GSE120046_76_lines <-
  FindNeighbors(merged_GSE120046_76_lines,
                dims = 1:30,
                verbose = T)
merged_GSE120046_76_lines <-
  FindClusters(merged_GSE120046_76_lines,
               verbose = T)

DimPlot(merged_GSE120046_76_lines,
        group.by = "subtype",
        label = T,
        repel = T)



Integrate_raw_df <-
  vector(mode = "list",
         length = 4L)
Integrate_raw_df[1:3] <-
  raw_76_lines_list[1:3]
names(Integrate_raw_df)

names(Integrate_raw_df) <-
  c("GABA", "nmglut", "npglut", "reference")
Integrate_raw_df[["reference"]] <-
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

# save.image("pre_merge_GSE120046_76_lines.RData")
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
                         dims = 1:20,
                         anchor.features = anchors_4_integration,
                         normalization.method = "SCT")
Integrate_df <-
  IntegrateData(anchorset = Integrate_df,
                dims = 1:30,
                normalization.method = "SCT",
                verbose = T)

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


Integrate_df_4_plot$cell.origin[Integrate_df_4_plot$subtype %in%
                                  c("GABA_iPS",
                                    "nmglut_iPS",
                                    "npglut_iPS")] <-
  "scARC"

DimPlot(Integrate_df_4_plot,
        group.by = "cell.origin",
        label = T,
        repel = T)

DimPlot(Integrate_df_4_plot,
        group.by = "type",
        label = T,
        repel = T)

DimPlot(Integrate_df_4_plot,
        group.by = "subtype",
        label = T,
        repel = T)

DimPlot(Integrate_df_4_plot,
        group.by = "HE.ident",
        label = T,
        repel = T)

DimPlot(Integrate_df_4_plot,
        group.by = "HE.ident.merged",
        label = T,
        repel = T)

Integrate_df_4_plot$HE.ident <-
  Integrate_df_4_plot$orig.ident

Integrate_df_4_plot$HE.ident[str_detect(string = Integrate_df_4_plot$orig.ident,
                                        pattern = "^HE.*",
                                        negate = T)] <-
  NA

Integrate_df_4_plot$HE.ident.merged <-
  str_remove(string = Integrate_df_4_plot$HE.ident,
             pattern = "[0-9]$")

FeaturePlot(Integrate_df_4_plot,
            features = "SLC17A6")

save.image("pre_merge_w_HE_ident.RData")

