# Siwei 21 Jun 2024
# Import scRNA-seq data of previous microglia to identify potential samples

# init ####
{
  library(Seurat)
  library(Signac)

  library(edgeR)

  library(future)
  library(parallel)

  library(stringr)

  library(harmony)

  library(readr)
}

plan("multisession", workers = 6)
set.seed(42)
options(future.globals.maxSize = 229496729600)


## load immune cells from "Human microglial state dynamics in Alzheimer’s disease progression" ####
Immune_cells <-
  readRDS("human_microglia_state_dynamics/Immune_cells.rds")

sun_et_al_cells <-
  readRDS("human_microglia_state_dynamics/ROSMAP.ImmuneCells.6regions.snRNAseq.counts.rds")
sun_et_al_cells_raw_meta <-
  readRDS("human_microglia_state_dynamics/ROSMAP.ImmuneCells.6regions.snRNAseq.meta.rds")

Idents(Immune_cells)
DimPlot(Immune_cells,
        alpha = 0.3)

microglia_P2RY12 <-
  Immune_cells[, Immune_cells$cell_type_high_resolution %in% c('Mic P2RY12')]

DimPlot(microglia_P2RY12)
Idents(microglia_P2RY12) <- "projid"
#
# MG_P2RY12_split <-
#   SplitObject(microglia_P2RY12,
#               split.by = "projid")

# collapse to pseudobulk by projid

pseudo_MG_P2RY12 <-
  AggregateExpression(microglia_P2RY12,
                      assays = "RNA",
                      return.seurat = F,
                      group.by = "projid",
                      verbose = T)

pseudo_MG_P2RY12 <-
  as.data.frame(pseudo_MG_P2RY12)

# sun_et_al_meta <-
#   readRDS("sun_et_al/personal.broadinstitute.org/cboix/sun_victor_et_al_data/ROSMAP.ImmuneCells.6regions.snRNAseq.meta.rds")
# sun_et_al_meta <-
#   as.data.frame(sun_et_al_meta)


# load ROSMAP lookup matrices ####
lookup_indivID_2_sunID <-
  read.csv(file = "human_microglia_state_dynamics/metadata/ROSMAP_meta_syn52293430/MIT_ROSMAP_Multiomics_individual_metadata.csv")
lookup_indivID_2_sunID <-
  lookup_indivID_2_sunID[, c(1:4, 19)]
lookup_indivID_2_sunID <-
  lookup_indivID_2_sunID[!duplicated(lookup_indivID_2_sunID$individualID), ]
lookup_indivID_2_sunID<-
  lookup_indivID_2_sunID[!duplicated(lookup_indivID_2_sunID$subject), ]

lookup_indivID_2_projID <-
  read.csv(file = "human_microglia_state_dynamics/metadata/ROSMAP_meta_all/ROSMAP_clinical.csv")
lookup_indivID_2_projID <-
  lookup_indivID_2_projID[!(duplicated(lookup_indivID_2_projID$projid)), ]
lookup_indivID_2_projID <-
  lookup_indivID_2_projID[!(duplicated(lookup_indivID_2_projID$individualID)), ]


lookup_indivID_projID_sunID <-
  merge(lookup_indivID_2_sunID,
        lookup_indivID_2_projID,
        by = "individualID")

# merge tables


## create Seurat object and merge/subset ##
length(unique(lookup_indivID_projID_sunID$individualID))
length(unique(lookup_indivID_projID_sunID$projid))
length(unique(lookup_indivID_projID_sunID$subject))
# length(unique(lookup_indivID_projID_sunID$))

dedup_lookup_indivID_projID_sunID <-
  lookup_indivID_projID_sunID[!(duplicated(lookup_indivID_projID_sunID$subject)), ]
dedup_lookup_indivID_projID_sunID <-
  lookup_indivID_projID_sunID[!(duplicated(lookup_indivID_projID_sunID$individualID)), ]
dedup_lookup_indivID_projID_sunID <-
  lookup_indivID_projID_sunID[!(duplicated(lookup_indivID_projID_sunID$projid)), ]

sum(duplicated(dedup_lookup_indivID_projID_sunID$individualID))
sum(duplicated(dedup_lookup_indivID_projID_sunID$projid))
sum(duplicated(dedup_lookup_indivID_projID_sunID$subject))

colnames(sun_et_al_cells) <-
  str_split(string = colnames(sun_et_al_cells),
            pattern = "\\.",
            simplify = T)[, 2]


sun_et_al_meta <-
  readRDS("sun_et_al/personal.broadinstitute.org/cboix/sun_victor_et_al_data/ROSMAP.ImmuneCells.6regions.snRNAseq.meta.rds")
sun_et_al_meta <-
  sun_et_al_meta[!duplicated(sun_et_al_meta$barcode), ]



sun_et_al_cell_meta <-
  merge(x = sun_et_al_meta,
        y = dedup_lookup_indivID_projID_sunID,
        by.x = "subject",
        by.y = "subject")
sun_et_al_cell_meta <-
  sun_et_al_cell_meta[!duplicated(sun_et_al_cell_meta$barcode), ]


sun_et_al_cells_filtered <-
  sun_et_al_cells[,
                  !duplicated(colnames(sun_et_al_cells))]
ncol(sun_et_al_cells_filtered)

sun_et_al_cells_filtered <-
  sun_et_al_cells_filtered[,
                  colnames(sun_et_al_cells_filtered) %in% sun_et_al_cell_meta$barcode]
ncol(sun_et_al_cells_filtered)
sun_et_al_cells_filtered <-
  sun_et_al_cells_filtered[,
                           order(colnames(sun_et_al_cells_filtered))]
rownames(sun_et_al_cells_filtered) ## use this one as count mx, ncol=171998 !!####
# mtx_sun_count <-
#   sun_et_al_cell_meta[, order(colnames(sun_et_al_cell_meta))]
# mtx_sun_count <-
#   mtx_sun_count[, !(duplicated(colnames(mtx_sun_count)))]


meta_sun_sorted <-
  sun_et_al_cell_meta[order(sun_et_al_cell_meta$barcode), ]
nrow(meta_sun_sorted) ## use this one as meta df, nrow=171998 !!####

rownames(meta_sun_sorted) <-
  meta_sun_sorted$barcode


Seurat_sun_et_al <-
  CreateSeuratObject(counts = sun_et_al_cells_filtered,
                     meta.data = meta_sun_sorted,
                     project = "MG")

saveRDS(Seurat_sun_et_al,
        file = "Seurat_sun_171998_cells.RDs")

save.image("workspace_21Jun2024.RData")



## recover memory by removing all other objects ####
gc()

{
  # Seurat_sun_et_al <-
  #   readRDS("Seurat_sun_171998_cells.RDs")

  ## use DLPFC 1st, then add hippocampus; ####
  # unique(Seurat_sun_et_al$brainRegion)
  # [1] "PFC"               "EntorhinalCortex"  "Thalamus"          "AngularGyrus"      "MidtemporalCortex"
  # [6] "Hippocampus"
  Seurat_sun_et_al <-
    Seurat_sun_et_al[, Seurat_sun_et_al$brainRegion == "PFC"]

  Seurat_sun_et_al <-
    Seurat_sun_et_al %>%
    SCTransform(seed.use = 42,
                vars.to.regress = c("percent.mt",
                                    "percent.rp"),
                verbose = T) %>%
    RunPCA(npcs = 50,
           verbose = T)

  # unique(Seurat_sun_et_al$subject)
  Seurat_sun_et_al$subject <-
    as.factor(Seurat_sun_et_al$subject)
  Seurat_sun_et_al$msex.x <-
    as.factor(Seurat_sun_et_al$msex.x)
  Seurat_sun_et_al$pmi.x <-
    as.factor(Seurat_sun_et_al$pmi.x)
  # Seurat_sun_et_al$brainRegion <-
  #   as.factor(Seurat_sun_et_al$brainRegion)

  # processed_sun_et_al <-
  #   IntegrateLayers(Seurat_sun_et_al,
  #                   method = HarmonyIntegration,
  #                   orig.reduction = "pca",
  #                   new.reduction = "harmony",
  #                   # normalzation.method = "SCT",
  #                   verbose = T)
  processed_sun_et_al <-
    RunHarmony(Seurat_sun_et_al,
               group.by.vars = c("subject",
                                 "msex.x",
                                 "pmi.x"),
               reduction.use = "pca",
               assay.use = "SCT",
               reduction.save = "harmony",
               # dims.use = 1:40,
               ncores = 40,
               plot_convergence = T,
               verbose = T)

  processed_sun_et_al <-
    processed_sun_et_al %>%
    FindNeighbors(dims = 1:30,
                  reduction = "harmony") %>%
    FindClusters(resolution = 0.5) %>%
    RunUMAP(dims = 1:30,
            seed.use = 42,
            reduction = "harmony",
            reduction.name = "umap.harmony")

}

DimPlot(processed_sun_et_al)
saveRDS(processed_sun_et_al,
        file = "Seurat_sun_PFC_cells_harmony.RDs")

## read back ####
processed_sun_et_al <-
  readRDS("Seurat_sun_PFC_cells_harmony.RDs")
unique(processed_sun_et_al)

# readin AMP-AD Chip imputed projids and compare, hg19 ####
Chip_Imputed_samples <-
  read.csv(file = "ROSMAP_genotypes_hg19/AMP-AD_ROSMAP_Rush-Broad_AffymetrixGenechip6_Imputed.csv")

Seurat_samples <-
  processed_sun_et_al$projid
Seurat_samples <-
  str_c("ROS",
        Seurat_samples)

Seurat_samples <-
  unique(Seurat_samples)

sum(Chip_Imputed_samples$gwas_id %in% Seurat_samples)
# 185 shared samples
Seurat_samples_in_genotype_mx <-
  Chip_Imputed_samples[Chip_Imputed_samples$gwas_id %in% Seurat_samples, ]
Seurat_samples_in_genotype_mx <-
  as.numeric(str_remove(Seurat_samples_in_genotype_mx,
                        pattern = "ROS"))

## subset Seurat object to get these 185 samples only
subsetted_185_samples_Seurat <-
  subset(Seurat_sun_et_al,
         subset = (projid %in% Seurat_samples_in_genotype_mx))

saveRDS(subsetted_185_samples_Seurat,
        file = "Seurat_sun_PFC_cells_185_genotypes.RDs")

# rs3851179 REF=C in hg19 !!
# rs3851179 C/C=0/0, linked to rs10792832 G/G (1/1) ####
# inverted between hg19 and hg38

df_rs3851179_raw <-
  read_delim("ROSMAP_genotypes_hg19/imputed_rs3851799_genotypes.txt",
             delim = "\t", escape_double = FALSE,
             trim_ws = TRUE)

df_rs3851179_transposed_hg19_CC_00 <-
  data.frame(individualID = colnames(df_rs3851179_raw),
             genotype = unlist(df_rs3851179_raw[1, ]))
df_rs3851179_transposed_hg19_CC_00 <-
  df_rs3851179_transposed_hg19_CC_00[-c(1:9), ]
df_rs3851179_transposed_hg19_CC_00$individualID[1]

df_rs3851179_transposed_hg19_CC_00$CEL <-
  str_split(string = df_rs3851179_transposed_hg19_CC_00$individualID,
            pattern = "\\.",
            simplify = T)[, 3]
df_rs3851179_transposed_hg19_CC_00$ROS_No <-
  str_split(string = df_rs3851179_transposed_hg19_CC_00$CEL,
            pattern = "_ROS",
            simplify = T)[, 2]

## ! deduplicate ####
df_rs3851179_transposed_hg19_CC_00 <-
  df_rs3851179_transposed_hg19_CC_00[!duplicated(df_rs3851179_transposed_hg19_CC_00$ROS_No), ]

sum(df_rs3851179_transposed_hg19_CC_00$ROS_No %in%
      processed_sun_et_al$projid) # 131

## subset 131 indiv from the intersection ####
## subset Seurat object to get these 135 samples only
subsetted_131_samples_Seurat <-
  subset(Seurat_sun_et_al,
         subset = (projid %in% df_rs3851179_transposed_hg19_CC_00$ROS_No))

length(Seurat_sun_et_al$projid %in% df_rs3851179_transposed_hg19_CC_00$ROS_No)
length(df_rs3851179_transposed_hg19_CC_00$ROS_No)

unique(df_rs3851179_transposed_hg19_CC_00$genotype)


PFC_metadata_info <-
  Seurat_sun_PFC_cells_harmony@meta.data
PFC_metadata_info <-
  PFC_metadata_info[!duplicated(PFC_metadata_info$individualID), ]

saveRDS(PFC_metadata_info,
        file = "Sun_PFC_metadata_info_14Jan2025.RDs")


PFC_174346_info <-
  Seurat_sun_174346_cells@meta.data
PFC_174346_info <-
  PFC_174346_info[!duplicated(PFC_174346_info$subject), ]
