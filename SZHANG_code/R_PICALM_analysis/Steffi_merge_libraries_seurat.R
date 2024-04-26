# Siwei 21 Sept 2023
# Merge Steffi's 4 libraries by Seurat

# put all h5 files under ./h5_input/ like below:
# > h5list
# [1] "h5_input/CD_IgG_filtered_feature_bc_matrix.h5"
# [2] "h5_input/CD_PDL1_filtered_feature_bc_matrix.h5"
# [3] "h5_input/MCD_IgG_filtered_feature_bc_matrix.h5"
# [4] "h5_input/MCD_PDL1_filtered_feature_bc_matrix.h5"
# and run the whole script in RStudio
# Will NOT generate images in the folder if run in R Console only !
# make sure all libraries have been properly installed

# init ####
library(Seurat)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(stringr)
library(readr)

library(pals)
library(future)
library(dplyr)

plan("multisession", workers = 2)
options(expressions = 20000)
options(future.globals.maxSize = 207374182400)


# read h5 and get sample names ####
h5list <-
  sort(list.files(path = "h5_input",
                  pattern = "filtered_feature_bc_matrix\\.h5",
                  recursive = T,
                  include.dirs = F,
                  full.names = T))

objlist <-
  vector(mode = "list",
         length = length(h5list))
names(objlist) <-
  str_remove_all(string = str_split(string = h5list,
                                    pattern = '/',
                                    simplify = T)[, 2],
                 pattern = "_filtered_feature_bc_matrix\\.h5")

for (i in 1:length(objlist)) {
  print(h5list[i])
  h5file <-
    Read10X_h5(filename = h5list[i])
  obj <-
    CreateSeuratObject(counts = h5file,
                       project = names(objlist)[i])
  obj$sample.origin <-
    names(objlist)[i]
  objlist[[i]] <- obj
  rm(list = c("h5file",
              "obj"))
}

# QC ####
# merge without correcting anything first to see distribution
rownames(objlist[[1]])[str_detect(string = rownames(objlist[[1]]),
                                  pattern = "mt")]

for (i in 1:length(objlist)) {
  objlist[[i]][["percent.mt"]] <-
    PercentageFeatureSet(objlist[[i]],
                         pattern = "^mt-")
}

for (i in 1:length(objlist)) {
  if (i == 1) {
    rough_merged <-
      objlist[[i]]
  } else {
    rough_merged <-
      merge(rough_merged,
            objlist[[i]])
  }
}

# check nfeature, ncount, pct mt distribution
VlnPlot(rough_merged,
        features = c("nFeature_RNA",
                     "nCount_RNA",
                     "percent.mt"),
        ncol = 3,
        pt.size = 0)
sum(rough_merged$percent.mt > 15) # 544
sum(rough_merged$nFeature_RNA > 7500) # 28
sum(rough_merged$nCount_RNA > 40000) # 255
sum(rough_merged$nFeature_RNA < 300) # 355
sum(rough_merged$nCount_RNA < 500) # 0

QCed_lst <-
  vector(mode = "list",
         length = length(objlist))
for (i in 1:length(objlist)) {
  QCed_lst[[i]] <-
    subset(objlist[[i]],
           subset = nFeature_RNA > 300 &
             nFeature_RNA < 7500 &
             nCount_RNA > 500 &
             nCount_RNA < 40000 &
             percent.mt < 20)
}

rm(rough_merged)
gc()

# normalize ####
transformed_lst <-
  vector(mode = "list",
         length(QCed_lst))

for (i in 1:length(QCed_lst)){
  print(paste0("now at ", i))
  # proceed with normalization
  obj <-
    PercentageFeatureSet(QCed_lst[[i]],
                         pattern = c("^mt-"),
                         col.name = "percent.mt")

  obj <-
    FindVariableFeatures(obj,
                         nfeatures = 8000)
  transformed_lst[[i]] <-
    SCTransform(obj,
                vars.to.regress = "percent.mt",
                method = "glmGamPoi",
                return.only.var.genes = F,
                variable.features.n = 8000,
                seed.use = 42,
                verbose = T)
}

# integration using reciprocal PCA ####
# find anchors
features <-
  SelectIntegrationFeatures(object.list = transformed_lst,
                            nfeatures = 3000,
                            fvf.nfeatures = 3000)
transformed_lst <-
  PrepSCTIntegration(transformed_lst,
                     anchor.features = features)
transformed_lst <-
  lapply(X = transformed_lst,
         FUN = function(x) {
           x <-
             ScaleData(x,
                       features = features)
           x <-
             RunPCA(x,
                    features = features)
         })

anchors <-
  FindIntegrationAnchors(object.list = transformed_lst,
                         anchor.features = features,
                         reference = c(1, 2, 3),
                         reduction = "rpca",
                         normalization.method = "SCT",
                         scale = F,
                         dims = 1:50,
                         verbose = T)
#Found 5701 anchors

# integrate
integrated <-
  IntegrateData(anchorset = anchors,
                verbose = T)

# QC ####
integrated[["percent.mt"]] <-
  PercentageFeatureSet(integrated,
                       assay = "RNA",
                       pattern = "^mt-")
DefaultAssay(integrated) <- "RNA"
Idents(integrated) <- "orig.ident"
VlnPlot(integrated,
        features = c("nFeature_RNA",
                     "nCount_RNA",
                     "percent.mt"),
        ncol = 3,
        fill.by = "feature",
        pt.size = 0)

# clustering ####
unique(integrated$sample.origin)

DefaultAssay(integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
integrated <- ScaleData(integrated)
integrated <-
  RunPCA(integrated,
         verbose = T,
         seed.use = 42)
integrated <-
  RunUMAP(integrated,
          reduction = "pca",
          dims = 1:30,
          seed.use = 42)
integrated <-
  FindNeighbors(integrated,
                reduction = "pca",
                dims = 1:30)
integrated <-
  FindClusters(integrated,
               resolution = 0.1, # adjust here for number of clusters, larger is more
               random.seed = 42)


## make plots ####
DimPlot(integrated,
        label = T,
        repel = T,
        group.by = "seurat_clusters") +
  # NoLegend() +
  ggtitle("Steffi data") +
  theme(text = element_text(size = 10),
        axis.text = element_text(size = 10))

lib_colors <-
  DiscretePalette(n = length(unique(Idents(integrated))),
                  "alphabet")
p <-
  DimPlot(integrated,
          label = F,
          cols = lib_colors,
          group.by = "sample.origin") +
  ggtitle("by origin") +
  theme(text = element_text(size = 10),
        axis.text = element_text(size = 10))
p$layers[[1]]$aes_params$alpha <- 0.2 # adjust here for effects
p

FeaturePlot(integrated,
            features = c("Ccr2",
                         "Aqp9",
                         "Gzma",
                         "Bcl11b"),
            ncol = 2)

save(transformed_lst,
     file = "Steffi_4_libraries_integrated_Seurat.RData")
