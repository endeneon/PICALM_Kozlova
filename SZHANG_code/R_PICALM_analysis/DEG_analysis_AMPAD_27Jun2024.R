# Siwei 27 Jun 2024
# Import scRNA-seq data of previous microglia to identify potential samples
# try different DEG analysis methods


# init ####
{
  library(Seurat)
  library(Signac)

  library(edgeR)
  library(DESeq2)
  library(MAST)

  library(future)

  library(stringr)

  library(harmony)

  library(readr)

  library(dplyr)

  library(ggplot2)
  library(RColorBrewer)
}

plan("multisession", workers = 6)
set.seed(42)
options(future.globals.maxSize = 229496729600)

subsetted_185_samples_Seurat <-
        readRDS(file = "Seurat_sun_PFC_cells_185_genotypes.RDs")

## try to use NormalizeData workflow than SCTransform
DefaultAssay(subsetted_185_samples_Seurat) <- "RNA"

# > colnames(subsetted_185_samples_Seurat@meta.data)
# [1] "orig.ident"             "nCount_RNA"             "nFeature_RNA"
# [4] "subject"                "brainRegion"            "batch"
# [7] "barcode"                "percent.mt"             "age_death.x"
# [10] "msex.x"                 "pmi.x"                  "ADdiag3types"
# [13] "percent.rp"             "seurat_clusters"        "individualID"
# [16] "individualIdSource"     "species"                "sex"
# [19] "projid"                 "Study"                  "msex.y"
# [22] "educ"                   "race"                   "spanish"
# [25] "apoe_genotype"          "age_at_visit_max"       "age_first_ad_dx"
# [28] "age_death.y"            "cts_mmse30_first_ad_dx" "cts_mmse30_lv"
# [31] "pmi.y"                  "braaksc"                "ceradsc"
# [34] "cogdx"                  "dcfdx_lv"               "rs561655_G_A"
# [37] "rs3851179_T_C"
subsetted_185_samples_Seurat_logNormalise_umap <-
  subsetted_185_samples_Seurat %>%
  NormalizeData(assay = "RNA",
                normalization.method = "LogNormalize",
                verbose = T) %>%
  FindVariableFeatures(selection.method = "vst",
                       binning.method = "equal_frequency",
                       verbose = T) %>%
  ScaleData(vars.to.regress = c("percent.mt",
                                "percent.rp",
                                "msex.x"),
            model.use = "poisson",
            use.umi = T,
            do.scale = T,
            do.center = T,
            verbose = T) %>%
  RunPCA(seed.use = 42,
         verbose = T) %>%
  FindNeighbors(reduction = "pca",
                verbose = T) %>%
  FindClusters(resolution = 0.5,
               algorithm = 1,
               random.seed = 42,
               verbose = T) %>%
  RunTSNE(reduction = "pca",
          seed.use = 42,
          tsne.method = "Rtsne",
          dims = 2:30)

## umap-learn did not work ####
subsetted_185_samples_Seurat_logNormalise_umap <-
  RunUMAP(subsetted_185_samples_Seurat_logNormalise_umap,
          dims = 1:30,
          # umap.method = "umap-learn",
          seed.use = 42,
          verbose = T)

DimPlot(subsetted_185_samples_Seurat_logNormalise_umap,
        group.by = "rs561655_G_A")


####
Idents(subsetted_185_samples_Seurat_logNormalise_umap) <- "rs561655_G_A"
subsetted_185_samples_GG_vs_AA <-
  FindMarkers(subsetted_185_samples_Seurat_logNormalise_umap,
              # assay = "RNA",
              # slot = "scale.data",
              ident.1 = 2,
              ident.2 = 0,
              latent.vars = "ADdiag3types",
              logfc.threshold = 0,
              test.use = "LR",
              # fc.results = "data",
              # fc.slot = "data",
              # mean.fxn = "data",
              base = 2,
              verbose = T)

subsetted_185_samples_GG_vs_AA[rownames(subsetted_185_samples_GG_vs_AA) == "PICALM", ]


unique(subsetted_185_samples_Seurat_logNormalise_umap$ADdiag3types)
Idents(subsetted_185_samples_Seurat_logNormalise_umap) <- "ADdiag3types"
subsetted_185_samples_GG_vs_AA <-
  FindMarkers(subsetted_185_samples_Seurat_logNormalise_umap,
              # assay = "RNA",
              # slot = "scale.data",
              ident.1 = "lateAD",
              ident.2 = "nonAD",
              # latent.vars = "rs561655_G_A",
              logfc.threshold = 0,
              test.use = "wilcox",
              # fc.results = "data",
              # fc.slot = "data",
              # mean.fxn = "data",
              base = 2,
              verbose = T)
subsetted_185_samples_GG_vs_AA[rownames(subsetted_185_samples_GG_vs_AA) == "PICALM", ]

# run a pseudobulk for fast data first ####, use PFC only
Idents(subsetted_185_samples_Seurat_logNormalise_umap) <- "brainRegion"
subsetted_185_PFC <-
  subset(subsetted_185_samples_Seurat_logNormalise_umap,
         subset = (brainRegion == "PFC"))




agg_185_samples <-
  AggregateExpression(subsetted_185_PFC,
                      assays = "RNA",
                      group.by = "projid",
                      verbose = T)
agg_185_samples <-
  as.matrix(agg_185_samples$RNA)
sum(agg_185_samples == 0)
sum(apply(agg_185_samples, 2, function(x) {
  return(is.integer(x))
}))


agg_185_metadata <-
  subsetted_185_PFC@meta.data
agg_185_metadata <-
  agg_185_metadata[!(duplicated(agg_185_metadata$projid)), ]

# run a quick edgeR/voom analysis


agg_185_metadata$batch <-
  as.factor(agg_185_metadata$batch)
agg_185_metadata$msex.x <-
  as.factor(agg_185_metadata$msex.x)
agg_185_metadata$ADdiag3types <-
  factor(agg_185_metadata$ADdiag3types,
         levels = c("nonAD",
                    "earlyAD",
                    "lateAD"))
agg_185_metadata$apoe_genotype <-
  as.factor(agg_185_metadata$apoe_genotype)
agg_185_metadata$braaksc <-
  as.factor(agg_185_metadata$braaksc)
agg_185_metadata$ceradsc <-
  as.factor(agg_185_metadata$ceradsc)
agg_185_metadata$cogdx <-
  as.factor(agg_185_metadata$cogdx)
agg_185_metadata$dcfdx_lv <-
  as.factor(agg_185_metadata$dcfdx_lv)
agg_185_metadata$rs561655_G_A <-
  factor(agg_185_metadata$rs561655_G_A,
         levels = c(0, 1, 2))

agg_185_metadata <-
  agg_185_metadata[order(agg_185_metadata$projid), ]
head(agg_185_metadata$projid)

agg_185_samples <-
  agg_185_samples[, order(colnames(agg_185_samples))]
head(colnames(agg_185_samples))

summary(agg_185_samples)
sum(agg_185_samples < 0)

DEG_185_agg <-
  DGEList(counts = agg_185_samples,
          # genes = rownames(agg_185_samples),
          # samples = colnames(agg_185_samples),
          group = agg_185_metadata$rs561655_G_A,
          remove.zeros = T)

table(filterByExpr(DEG_185_agg))

DEG_185_agg <-
  DEG_185_agg[filterByExpr(DEG_185_agg),
              ,
              keep.lib.sizes = F]
sum(rownames(DEG_185_agg) %in% "PICALM")


DEG_185_agg <-
  calcNormFactors(DEG_185_agg)
DEG_185_agg$samples
plotMDS(DEG_185_agg)

mx_design_185 <-
  model.matrix(~ rs561655_G_A +
                 # dcfdx_lv +
                 # cogdx +
                 # ceradsc +
                 # braaksc +
                 # apoe_genotype +
                 ADdiag3types +
                 # batch +
                 msex.x,
               data = agg_185_metadata)



DEG_185_agg <-
  estimateDisp(DEG_185_agg)
DEG_185_agg$common.dispersion
plotBCV(DEG_185_agg)

fit_185_agg <-
  glmQLFit(DEG_185_agg,
           design = mx_design_185,
           robust = T)
plotQLDisp(fit_185_agg)

View(mx_design_185)
summary(fit_185_agg)
summary(DEG_185_agg)

sum(fit_185_agg$counts < 0)
str(fit_185_agg$counts)
sum(apply(fit_185_agg$counts, 2, function(x) {
  return(is.integer(x))
}))

fit_185_F <-
  edgeR::glmQLFTest(fit_185_agg,
                    contrast = c(-1, 0, 0, 1, 0, 0))

fit_185_df <-
  fit_185_F$table
fit_185_df[rownames(fit_185_df) %in% "PICALM", ]

cpm_PICALM <-
  as.data.frame(edgeR::cpm(fit_185_F))
cpm_PICALM <-
  cpm_PICALM[rownames(cpm_PICALM) %in% "PICALM", ]

df_2_plot <-
  data.frame(cpm = unlist(cpm_PICALM),
             genotype = factor(agg_185_metadata$rs561655_G_A),
             diag = factor(agg_185_metadata$ADdiag3types))

ggplot(df_2_plot,
       aes(x = genotype,
           y = cpm,
           colour = diag)) +
  geom_boxplot(outlier.size = 0.5) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2),
             size = 0.5) +
  scale_colour_manual(values = brewer.pal(n = 3,
                                          name = "Dark2")) +
  ylim(0, 11) +
  theme_classic()

ggplot(df_2_plot,
       aes(x = diag,
           y = cpm,
           colour = diag)) +
  geom_boxplot(outlier.size = 0.5) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2),
             size = 0.5) +
  scale_colour_manual(values = brewer.pal(n = 3,
                                          name = "Dark2")) +
  # ylim(0, 11) +
  theme_classic()


t_cpm <-
  as.data.frame(t(cpm_PICALM))
all(rownames(t_cpm) == rownames(agg_185_metadata))

t_cpm$barcode <-
  rownames(t_cpm)
t_cpm <-
  t_cpm[order(t_cpm$barcode), ]


agg_185_metadata_reordered <-
  agg_185_metadata[order(agg_185_metadata$barcode), ]
all(rownames(t_cpm) == agg_185_metadata_reordered$barcode)

df_plot_4_jubao <-
  merge(x = t_cpm,
        y = agg_185_metadata_reordered,
        by = "barcode")
df_plot_4_jubao$barcode <- NULL
df_plot_4_jubao$rs10792832_geno <-
  "AG"
df_plot_4_jubao$rs10792832_geno[df_plot_4_jubao$rs561655_G_A == 0] <- "AA"
df_plot_4_jubao$rs10792832_geno[df_plot_4_jubao$rs561655_G_A == 2] <- "GG"

ggplot(df_plot_4_jubao,
       aes(x = rs10792832_geno,
           y = PICALM,
           colour = ADdiag3types)) +
  geom_boxplot(outlier.size = 0.5) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2),
             size = 0.5) +
  scale_colour_manual(values = brewer.pal(n = 3,
                                          name = "Dark2")) +
  # ylim(0, 11) +
  theme_classic()

write.table(df_plot_4_jubao,
            file = "PICALM_cpm_indiv_w_metadata_28Jun2024.tsv",
            quote = F, sep = "\t",
            row.names = F, col.names = T)


## calculate genotype vs exp
df_plot_4_jubao$rs10792832_geno <-
  factor(df_plot_4_jubao$rs10792832_geno,
         levels = c("AA",
                    "AG",
                    "GG"))

ggplot(df_plot_4_jubao,
       aes(x = rs10792832_geno,
           y = PICALM,
           group = rs10792832_geno)) +
  geom_boxplot(outlier.size = 1,
               width = 0.5) +
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               dotsize = 0.5,
               binwidth = 0.0002,stackratio = 0.1,
               position = position_jitter(width = 0.2)) +
  # geom_point(position = "jitter",
  #            size = 0.2) +
  theme_classic()

cpm_target <-
  as.data.frame(edgeR::cpm(fit_185_F))
cpm_target <-
  cpm_target[rownames(cpm_target) %in% "APOE", ]
rownames(cpm_target)[str_detect(string = rownames(cpm_target),
                                pattern = "APOE")]
