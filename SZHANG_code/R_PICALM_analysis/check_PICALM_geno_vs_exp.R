# Siwei 13 Jan 2025

# lookup projID

{
  library(stringr)
  library(Seurat)

  library(parallel)
  library(future)


  library(glmGamPoi)

  library(edgeR)

  library(data.table)

  library(readr)

  plan("multisession", workers = 3)
  # options(mc.cores = 32)
  set.seed(42)
  options(future.globals.maxSize = 429496729600)

  setwd("~/backuped_space/Siwei_misc_R_projects/R_Alena_PICALM")
}

Immune_cells <-
  readRDS("~/Data/FASTQ/sage_synapse/syn52368912/Immune_cells.rds")
Immune_cells_meta <-
  Immune_cells@meta.data
Immune_cells_meta$barcode <-
  rownames(Immune_cells_meta)
unique(str_split(Immune_cells_meta$barcode,
                 pattern = '-',
                 simplify = T)[, 3])
length(unique(Immune_cells_meta$projid))
Immune_cells_meta$trimmed_barcode <-
  str_c(str_split(string = Immune_cells_meta$barcode,
                  pattern = '-',
                  simplify = T)[, 1],
        str_split(string = Immune_cells_meta$barcode,
                  pattern = '-',
                  simplify = T)[, 2],
        sep = '-')

# Immune_cells_cluster_0 <-
#   Immune_cells_meta[str_detect(string = Immune_cells_meta$barcode,
#                                pattern = "\\-0$"), ]
# unique(Immune_cells_cluster_0$projid)



Immune_lookup_table <-
  Immune_cells_meta[!duplicated(Immune_cells_meta$projid), ]

Sun_meta <-
  readRDS("sun_et_al/personal.broadinstitute.org/cboix/sun_victor_et_al_data/ROSMAP.ImmuneCells.6regions.snRNAseq.meta.rds")
length(unique(Sun_meta$subject))


combined_Immune_cells_meta_Sun_meta <-
  merge(x = Immune_cells_meta,
        y = Sun_meta,
        by.x = "trimmed_barcode",
        by.y = "barcode")
length(unique(combined_Immune_cells_meta_Sun_meta$projid))

Sun_immune_lookup_table <-
  combined_Immune_cells_meta_Sun_meta[!duplicated(combined_Immune_cells_meta_Sun_meta$projid), ]
Sun_immune_lookup_table <-
  Sun_immune_lookup_table[, c(2, 5, 15)]

Immune_cells_meta_reassign <-
  merge(x = Immune_cells_meta,
        y = Sun_immune_lookup_table,
        by.x = "projid",
        by.y = "projid",
        all.x = T)

Immune_cells$ADdiag3types <-
  Immune_cells_meta_reassign$ADdiag3types

Immune_cells_noNA <-
  Immune_cells[, !is.na(Immune_cells$ADdiag3types)]
Immune_cells_noNA <-
  Immune_cells_noNA[, Immune_cells_noNA$cell_type_high_resolution == 'Mic P2RY12']

Immune_cells_noNA <-
  SCTransform(Immune_cells_noNA,
              vst.flavor = "v2",
              verbose = T)

Idents(Immune_cells_noNA) <- "projid"

Immune_pseudobulk <-
  PseudobulkExpression(Immune_cells_noNA,
                       assays = "RNA",
                       return.seurat = F,
                       group.by = "ident",
                       normalization.method = "LogNormalize")
Immune_pseudobulk <-
  Immune_pseudobulk$RNA
Immune_pseudobulk_PICALM <-
  data.frame(expression = Immune_pseudobulk[rownames(Immune_pseudobulk) == "PICALM", ])
Immune_pseudobulk_PICALM$g_projid <-
  rownames(Immune_pseudobulk_PICALM)

Immune_lookup_table$g_projid <-
  str_c("g",
        Immune_lookup_table$projid,
        sep = "")
Sun_immune_lookup_table$g_projid <-
  str_c("g",
        Sun_immune_lookup_table$projid,
        sep = "")

Immune_pseudobulk_PICALM <-
  merge(x = Immune_pseudobulk_PICALM,
        y = Sun_immune_lookup_table,
        by = "g_projid")

wilcox.test(x = Immune_pseudobulk_PICALM$expression[Immune_pseudobulk_PICALM$ADdiag3types == "nonAD"],
            y = Immune_pseudobulk_PICALM$expression[Immune_pseudobulk_PICALM$ADdiag3types == "lateAD"],
            paired = F,
            alternative = "t")

t.test(x = Immune_pseudobulk_PICALM$expression[Immune_pseudobulk_PICALM$ADdiag3types == "nonAD"],
            y = Immune_pseudobulk_PICALM$expression[Immune_pseudobulk_PICALM$ADdiag3types == "lateAD"],
            paired = F, var.equal = T,
            alternative = "t")

wilcox.test(x = Immune_pseudobulk_PICALM$expression[Immune_pseudobulk_PICALM$ADdiag3types == "nonAD"],
            y = Immune_pseudobulk_PICALM$expression[Immune_pseudobulk_PICALM$ADdiag3types == "earlyAD"],
            paired = F,
            alternative = "t")

mean(Immune_pseudobulk_PICALM$expression[Immune_pseudobulk_PICALM$ADdiag3types == "nonAD"])
mean(Immune_pseudobulk_PICALM$expression[Immune_pseudobulk_PICALM$ADdiag3types == "lateAD"])
mean(Immune_pseudobulk_PICALM$expression[Immune_pseudobulk_PICALM$ADdiag3types == "earlyAD"])


# write a function to loop through all data files

calc_PICALM_type <-
  function(x) {

    # Immune_cells <-
    #   readRDS(x)
    Immune_cells <-
      readRDS(all_RDS_files[5])
    Immune_cells_meta <-
      Immune_cells@meta.data
    Immune_cells_meta$barcode <-
      rownames(Immune_cells_meta)
    print(unique(str_split(Immune_cells_meta$barcode,
                           pattern = '-',
                           simplify = T)[, 3]))
    print(length(unique(Immune_cells_meta$projid)))
    Immune_cells_meta$trimmed_barcode <-
      str_c(str_split(string = Immune_cells_meta$barcode,
                      pattern = '-',
                      simplify = T)[, 1],
            str_split(string = Immune_cells_meta$barcode,
                      pattern = '-',
                      simplify = T)[, 2],
            sep = '-')

    # Immune_cells_cluster_0 <-
    #   Immune_cells_meta[str_detect(string = Immune_cells_meta$barcode,
    #                                pattern = "\\-0$"), ]
    # unique(Immune_cells_cluster_0$projid)



    Immune_lookup_table <-
      Immune_cells_meta[!duplicated(Immune_cells_meta$projid), ]

    Sun_meta <-
      readRDS("sun_et_al/personal.broadinstitute.org/cboix/sun_victor_et_al_data/ROSMAP.ImmuneCells.6regions.snRNAseq.meta.rds")
    print(length(unique(Sun_meta$subject)))

    meta_indiv_index_from_immune <-
      readRDS("~/backuped_space/Siwei_misc_R_projects/R_Alena_PICALM/meta_indiv_index_from_immune.RDs")

    # combined_Immune_cells_meta_Sun_meta <-
    #   merge(x = Immune_cells_meta,
    #         y = Sun_meta,
    #         by.x = "trimmed_barcode",
    #         by.y = "barcode")
    # print(length(unique(combined_Immune_cells_meta_Sun_meta$projid)))

    Sun_immune_lookup_table <-
      combined_Immune_cells_meta_Sun_meta[!duplicated(combined_Immune_cells_meta_Sun_meta$projid), ]
    Sun_immune_lookup_table <-
      Sun_immune_lookup_table[, c(2, 5, 15)]

    Immune_cells_meta_reassign <-
      merge(x = Immune_cells_meta,
            y = Sun_immune_lookup_table,
            by.x = "projid",
            by.y = "projid",
            all.x = T)

    Immune_cells$ADdiag3types <-
      Immune_cells_meta_reassign$ADdiag3types

    print('=== 3 === ')

    Immune_cells <-
      UpdateSeuratObject(Immune_cells)

    Immune_cells_noNA <-
      Immune_cells[, !is.na(Immune_cells$ADdiag3types)]
    # Immune_cells_noNA <-
    #   Immune_cells_noNA[, Immune_cells_noNA$cell_type_high_resolution == 'Mic P2RY12']

    print('=== assigned ADdiag3types === ')

    Immune_cells_noNA <-
      NormalizeData(Immune_cells_noNA,
                    assay = "RNA",
                    normalization.method = "LogNormalize")

    # Immune_cells_noNA <-
    #   SCTransform(Immune_cells_noNA,
    #               vst.flavor = "v2",
    #               verbose = T)

    Idents(Immune_cells_noNA) <- "projid"

    Immune_pseudobulk <-
      AggregateExpression(Immune_cells_noNA,
                           assays = "RNA",
                           return.seurat = F,
                           group.by = "ident",
                           normalization.method = "LogNormalize")
    Immune_pseudobulk <-
      Immune_pseudobulk$RNA

    Immune_pseudobulk_DGE <-
      DGEList(counts = as.matrix(Immune_pseudobulk),
              sample = colnames(Immune_pseudobulk),
              genes = rownames(Immune_pseudobulk))
    Immune_pseudobulk_logCPM <-
      as.data.frame(cpm(Immune_pseudobulk_DGE,
                        log = T))

    Immune_pseudobulk_logCPM <-
      data.frame(expression = Immune_pseudobulk_logCPM[rownames(Immune_pseudobulk_logCPM) == "PICALM", ])
    Immune_pseudobulk_logCPM$g_projid <-
      rownames(Immune_pseudobulk_logCPM)

    Immune_lookup_table$g_projid <-
      str_c("g",
            Immune_lookup_table$projid,
            sep = "")
    Sun_immune_lookup_table$g_projid <-
      str_c("g",
            Sun_immune_lookup_table$projid,
            sep = "")

    Immune_pseudobulk_logCPM <-
      merge(x = Immune_pseudobulk_logCPM,
            y = Sun_immune_lookup_table,
            by = "g_projid")

    print(Immune_pseudobulk_logCPM)

    return(list(df = Immune_pseudobulk_logCPM
                # df2 = Immune_pseudobulk
                # nonADvslateAD = wilcox.test(x = Immune_pseudobulk_PICALM$expression[Immune_pseudobulk_PICALM$ADdiag3types == "nonAD"],
                #             y = Immune_pseudobulk_PICALM$expression[Immune_pseudobulk_PICALM$ADdiag3types == "lateAD"],
                #             paired = F,
                #             alternative = "t"),
                #
                # nonADvsearlyAD = wilcox.test(x = Immune_pseudobulk_PICALM$expression[Immune_pseudobulk_PICALM$ADdiag3types == "nonAD"],
                #             y = Immune_pseudobulk_PICALM$expression[Immune_pseudobulk_PICALM$ADdiag3types == "earlyAD"],
                #             paired = F,
                #             alternative = "t"),

                # mean_nonAD = mean(Immune_pseudobulk_PICALM$expression[Immune_pseudobulk_PICALM$ADdiag3types == "nonAD"]),
                # mean_lateAD = mean(Immune_pseudobulk_PICALM$expression[Immune_pseudobulk_PICALM$ADdiag3types == "lateAD"]),
                # mean_earlyAD = mean(Immune_pseudobulk_PICALM$expression[Immune_pseudobulk_PICALM$ADdiag3types == "earlyAD"])))
    ))
    #
  }

#
all_RDS_files <-
  list.files("~/Data/FASTQ/sage_synapse/syn52368912",
             pattern = ".*.rds",
             full.names = T,
             recursive = F)
names(all_RDS_files) <-
  str_split(string = all_RDS_files,
            pattern = '/',
            simplify = T)[, 8]
# names(all_RDS_files)


PICALM_results <-
  vector(mode = "list",
         length = length(all_RDS_files))
names(PICALM_results) <-
  names(all_RDS_files)

for (i in 1:(length(all_RDS_files) - 1)) {
  print(all_RDS_files[i])

  try({
    PICALM_results[[i]] <-
      calc_PICALM_type(x = all_RDS_files[i])
  })

  # PICALM_results
}

saveRDS(PICALM_results,
        file = "PICALM_results_cell_types_new.RDs")

#
# Immune_lookup_table$ROS_id <- NA
#
# Immune_lookup_table$ROS_id[Immune_lookup_table$bare_id < 10000000] <-
#   str_c("ROS0",
#         Immune_lookup_table$bare_id[Immune_lookup_table$bare_id < 10000000])
# Immune_lookup_table$ROS_id[Immune_lookup_table$bare_id > 10000000] <-
#   str_c("ROS",
#         Immune_lookup_table$bare_id[Immune_lookup_table$bare_id > 10000000])
#
# subsetted_185_samples_Seurat <-
#   readRDS("Seurat_sun_PFC_cells_185_genotypes.RDs")
#
# sum(Immune_lookup_table$ROS_id %in% unique(subsetted_185_samples_Seurat$projid))
#
# subsetted_185_samples_Seurat$bare_projid <-
#   as.numeric(str_remove_all(string = subsetted_185_samples_Seurat$projid,
#                             pattern = "ROS"))
# sum(is.na(subsetted_185_samples_Seurat$bare_projid))
#
# Immune_cells$AD_type <- NA
# Immune_cells$AD_type[Immune_cells$projid %in% subsetted_185_samples_Seurat$bare_projid]
#
# length(Immune_cells$projid)
# length(Immune_cells$projid %in% subsetted_185_samples_Seurat$bare_projid[subsetted_185_samples_Seurat$ADdiag3types == "nonAD"])
#
# Immune_cells$AD_type[Immune_cells$projid %in% unique(subsetted_185_samples_Seurat$bare_projid[subsetted_185_samples_Seurat$ADdiag3types == "nonAD"])] <- "nonAD"
# Immune_cells$AD_type[Immune_cells$projid %in% unique(subsetted_185_samples_Seurat$bare_projid[subsetted_185_samples_Seurat$ADdiag3types == "earlyAD"])] <- "earlyAD"
# Immune_cells$AD_type[Immune_cells$projid %in% unique(subsetted_185_samples_Seurat$bare_projid[subsetted_185_samples_Seurat$ADdiag3types == "lateAD"])] <- "lateAD"
#
# Immune_Microglia <-
#   Immune_cells[, Immune_cells$cell_type_high_resolution == 'Mic P2RY12']
#
# unique(Immune_Microglia$AD_type)
# Immune_Microglia <-
#   Immune_Microglia[, !(is.na(Immune_Microglia$AD_type))]
#
#
# Immune_Microglia <-
#   SCTransform(Immune_Microglia,
#               vst.flavor = "v2",
#               verbose = T,
#               seed.use = 42)
#
# # Idents(Immune_Microglia) <- "AD_type"
# # Idents(Immune_Microglia) <- "projid"
# #
# # Immune_Microglia <-
# #   Immune_Microglia %>%
# #   NormalizeData(normalization.method = "LogNormalize",
# #                 assay = "RNA") %>%
# #   FindVariableFeatures() %>%
# #   ScaleData()
#
# Idents(Immune_Microglia) <- "projid"
#
# # VlnPlot(Immune_Microglia,
# #         features = c("PICALM"))
#
#
# Immune_pseudobulk <-
#   PseudobulkExpression(Immune_Microglia,
#                        assays = "SCT",
#                        return.seurat = F,
#                        group.by = "ident",
#                        normalization.method = "LogNormalize")
#
# pseudobulk_group <-
#   Immune_Microglia@meta.data
# pseudobulk_group <-
#   pseudobulk_group[!duplicated(pseudobulk_group$projid), ]
# pseudobulk_group$groupid <-
#   str_c("g",
#         pseudobulk_group$projid,
#         sep = "")
# pseudobulk_group <-
#   pseudobulk_group[!(is.na(pseudobulk_group$AD_type)), ]
#
#
# pseudobulk_exp_ <-
#   as.data.frame(Immune_pseudobulk$SCT[rownames(Immune_pseudobulk$SCT) == "PICALM", ])
# pseudobulk_exp_$group_name <-
#   rownames(pseudobulk_exp_)
#
# pseudobulk_exp_ <-
#   merge(x = pseudobulk_exp_,
#         y = pseudobulk_group,
#         by.x = "group_name",
#         by.y = "groupid")
# colnames(pseudobulk_exp_)[2] <- "PICALM_exp"
#
# mean(pseudobulk_exp_$PICALM_exp[pseudobulk_exp_$AD_type == "nonAD"])
# mean(pseudobulk_exp_$PICALM_exp[pseudobulk_exp_$AD_type == "lateAD"])
# mean(pseudobulk_exp_$PICALM_exp[pseudobulk_exp_$AD_type == "earlyAD"])
#
# sum(pseudobulk_exp_$AD_type == "nonAD")
# sum(pseudobulk_exp_$AD_type == "lateAD")
# sum(pseudobulk_exp_$AD_type == "earlyAD")
#
# wilcox.test(x = pseudobulk_exp_$PICALM_exp[pseudobulk_exp_$AD_type == "nonAD"],
#             y = pseudobulk_exp_$PICALM_exp[pseudobulk_exp_$AD_type == "lateAD"],
#             paired = F,
#             alternative = "t")
#
# wilcox.test(x = pseudobulk_exp_$PICALM_exp[pseudobulk_exp_$AD_type == "nonAD"],
#             y= pseudobulk_exp_$PICALM_exp[pseudobulk_exp_$AD_type == "earlyAD"],
#             paired = F,
#             alternative = "t")
#
#
# # BiocManager::install('glmGamPoi')
#
# sun_ROSMAP_meta <-
#   readRDS("sun_et_al/personal.broadinstitute.org/cboix/sun_victor_et_al_data/ROSMAP.ImmuneCells.6regions.snRNAseq.meta.rds")
# sun_ROSMAP_count <-
#   readRDS("sun_et_al/personal.broadinstitute.org/cboix/sun_victor_et_al_data/ROSMAP.ImmuneCells.6regions.snRNAseq.counts.rds")

Sun_185_projid <-
  Seurat_sun_PFC_cells_185_genotypes@meta.data
Sun_185_projid <-
  Sun_185_projid[!duplicated(Sun_185_projid$subject), ]
Sun_185_projid <-
  Sun_185_projid[Sun_185_projid$Study == "ROS", ]
Sun_185_projid$bare_projid <-
  as.numeric(str_remove_all(Sun_185_projid$projid,
                            pattern = "^ROS"))


AMP_AD_ROSMAP_samples <-
  read.table("ROSMAP_genotypes_hg19/AMP-AD_ROSMAP_Rush-Broad_AffymetrixGenechip6_Imputed.fam",
             header = T, sep = " ")
rs3851179_genotype_raw <-
  read_delim("ROSMAP_genotypes_hg19/grepped_rs3851179.txt",
             delim = " ", escape_double = FALSE,
             col_names = FALSE, trim_ws = TRUE)
rs3851179_genotype_raw <-
  as.data.frame(t(rs3851179_genotype_raw))
  # read.table("ROSMAP_genotypes_hg19/bim_sets/rs3851179.txt",
  #            sep = "\t", header = T)
rs3851179_genotype <-
  rs3851179_genotype_raw[rs3851179_genotype_raw$V2 != './.', ]
rs3851179_genotype$ROSMAP <-
  str_split(rs3851179_genotype$V1,
            pattern = "_",
            simplify = T)[, 2]


rs3851179_genotype_Sun <-
  rs3851179_genotype[rs3851179_genotype$ROSMAP %in% Sun_185_projid$projid, ]
table(rs3851179_genotype_Sun$V2)
