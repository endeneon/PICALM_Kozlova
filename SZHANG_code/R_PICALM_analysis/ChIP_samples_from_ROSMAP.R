# Siwei 21 Jun 2024
# Import scRNA-seq data of previous microglia to identify potential samples

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
}

plan("multisession", workers = 6)
set.seed(42)
options(future.globals.maxSize = 229496729600)

## read back ####
processed_sun_et_al <-
  readRDS("Seurat_sun_PFC_cells_harmony.RDs")
unique(processed_sun_et_al)

# readin AMP-AD Chip imputed projids and compare, hg19 ####

AMPAD_ROSMAP_GWAS_samples <-
  read_delim("ROSMAP_genotypes_hg19/AMP-AD_ROSMAP_Rush-Broad_AffymetrixGenechip6_Imputed.csv",
             delim = "\t", escape_double = FALSE,
             trim_ws = TRUE)
AMPAD_ROSMAP_GWAS_samples <-
  unlist(AMPAD_ROSMAP_GWAS_samples)
AMPAD_ROSMAP_GWAS_samples <-
  c("rsid", "REF", "Alt",
    AMPAD_ROSMAP_GWAS_samples)

## hg19 REFs: rs10792832=A, rs3851179=T, rs561655=G; ####
## hg19 Alts: rs10792832=G, rs3851179=C, rs561655=A; ####
## dosage = 0*(0/0) + 1*(0/1) + 2*(1/1)
AMPAD_rs10792832_LD_rsids <-
  read_table("ROSMAP_genotypes_hg19/chr11_561655.txt",
             col_names = FALSE)
colnames(AMPAD_rs10792832_LD_rsids) <-
  AMPAD_ROSMAP_GWAS_samples
Chip_Imputed_samples <-
  AMPAD_rs10792832_LD_rsids[,
                            4:ncol(AMPAD_rs10792832_LD_rsids)]
rownames(Chip_Imputed_samples) <-
  str_c(AMPAD_rs10792832_LD_rsids$rsid,
        AMPAD_rs10792832_LD_rsids$REF,
        AMPAD_rs10792832_LD_rsids$Alt,
        sep = "_")
Chip_Imputed_samples <-
  as.data.frame(t(Chip_Imputed_samples))
Chip_Imputed_samples <-
  Chip_Imputed_samples[-c(nrow(Chip_Imputed_samples)), ]

ROS_ChIP_imputed_samples <-
  Chip_Imputed_samples[str_detect(string = rownames(Chip_Imputed_samples),
                                  pattern = "^ROS.*")
                       ,]
ROS_ChIP_imputed_samples <-
  ROS_ChIP_imputed_samples[order(rownames(ROS_ChIP_imputed_samples)), ]

Seurat_samples <-
  processed_sun_et_al$projid
Seurat_samples <-
  str_c("ROS",
        Seurat_samples)

Seurat_samples <-
  unique(Seurat_samples)
Seurat_samples <-
  Seurat_samples[order(Seurat_samples)]

sum(rownames(ROS_ChIP_imputed_samples) %in% Seurat_samples)
# 185 shared samples

Seurat_samples_in_genotype_mx <-
  ROS_ChIP_imputed_samples[(rownames(ROS_ChIP_imputed_samples) %in% Seurat_samples), ]
Seurat_samples_in_genotype_mx <-
  Seurat_samples_in_genotype_mx[order(rownames(Seurat_samples_in_genotype_mx)), ]
# Seurat_samples_in_genotype_mx <-
#   as.numeric(str_remove(rownames(Seurat_samples_in_genotype_mx),
#                         pattern = "ROS"))

## subset Seurat object to get these 185 samples only
Seurat_sun_et_al$projid <-
  str_c("ROS",
        Seurat_sun_et_al$projid)
subsetted_185_samples_Seurat <-
  subset(Seurat_sun_et_al,
         subset = (projid %in% rownames(Seurat_samples_in_genotype_mx)))

# saveRDS(subsetted_185_samples_Seurat,
#         file = "Seurat_sun_PFC_cells_185_genotypes.RDs")

# rs3851179 REF=T in hg19 !!
# rs3851179 T/T=0/0, linked to rs10792832 A/A (0/0) ####
temp_metadata <-
  subsetted_185_samples_Seurat@meta.data
Seurat_samples_in_genotype_mx$projid <-
  rownames(Seurat_samples_in_genotype_mx)
temp_metadata <-
  merge(x = temp_metadata,
        y = Seurat_samples_in_genotype_mx,
        by = "projid")
subsetted_185_samples_Seurat$rs561655_G_A <-
  temp_metadata$rs561655_G_A
subsetted_185_samples_Seurat$rs3851179_T_C <-
  temp_metadata$rs3851179_T_C
rm(temp_metadata)

saveRDS(subsetted_185_samples_Seurat,
        file = "Seurat_sun_PFC_cells_185_genotypes.RDs")


## it seems rs561655 is more in LD w/ rs10792832 (integers)
unique(subsetted_185_samples_Seurat$rs561655_G_A)
Idents(subsetted_185_samples_Seurat) <- "rs561655_G_A"
subsetted_185_samples_GG_vs_AA <-
    FindMarkers(subsetted_185_samples_Seurat,
                assay = "SCT",
                slot = "data",
                ident.1 = 2,
                ident.2 = 0,
                # latent.vars = "rs561655_G_A",
                logfc.threshold = 0,
                test.use = "DESeq2",
                # fc.results = "data",
                # fc.slot = "data",
                # mean.fxn = "data",
                base = 2,
                verbose = T)



data("maits",
     package = 'MAST')
MAST_maits <-
  maits
class

# # just test the dosage model, no zlm hurdle design matrix applied
# Idents(subsetted_185_samples_Seurat) <- "orig.ident"
# PICALM_rs10792832 <-
#   FindMarkers(subsetted_185_samples_Seurat,
#               # assay = "SCT",
#               # slot = "data",
#               ident.1 = "immune",
#               latent.vars = "rs561655_G_A",
#               logfc.threshold = 0,
#               test.use = "MAST",
#               # fc.results = "data",
#               # fc.slot = "data",
#               # mean.fxn = "data",
#               base = 2,
#               verbose = T)

# Build the SummarizedExperiment() for MAST
View(head(subsetted_185_samples_Seurat@assays$SCT@scale.data))
subsetted_185_samples_Seurat <-
  NormalizeData(subsetted_185_samples_Seurat,
                normalization.method = "CLR",
                margin = 2,
                verbose = T,
                assay = "RNA")
