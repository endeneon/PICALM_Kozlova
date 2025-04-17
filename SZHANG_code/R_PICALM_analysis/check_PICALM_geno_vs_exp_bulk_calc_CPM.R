# Siwei 13 Jan 2025

# lookup projID

{
  library(stringr)
  library(Seurat)

  library(parallel)
  library(future)


  library(glmGamPoi)

  library(data.table)
  library(limma)
  library(edgeR)

  library(ggplot2)
  library(RColorBrewer)
}


plan("multisession", workers = 3)
# options(mc.cores = 32)
set.seed(42)
options(future.globals.maxSize = 429496729600)

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

wilcox.test(x = Immune_pseudobulk_PICALM$expression[Immune_pseudobulk_PICALM$ADdiag3types == "nonAD"],
            y = Immune_pseudobulk_PICALM$expression[Immune_pseudobulk_PICALM$ADdiag3types == "earlyAD"],
            paired = F,
            alternative = "t")

mean(Immune_pseudobulk_PICALM$expression[Immune_pseudobulk_PICALM$ADdiag3types == "nonAD"])
mean(Immune_pseudobulk_PICALM$expression[Immune_pseudobulk_PICALM$ADdiag3types == "lateAD"])
mean(Immune_pseudobulk_PICALM$expression[Immune_pseudobulk_PICALM$ADdiag3types == "earlyAD"])


# filter out low-expression genes ####
FilterGenes <-
  function(object, min.value=1, min.cells = 0, genes = NULL) {
    parameters.to.store <- as.list(environment(), all = TRUE)[names(formals("FilterGenes"))]
    object <- Seurat:::SetCalcParams(object = object, calculation = "FilterGenes", ... = parameters.to.store)
    genes.use <- rownames(object@data)

    if (!is.null(genes)) {
      genes.use <- intersect(genes.use, genes)
      object@data <- object@data[genes.use, ]
      return(object)
    } else if (min.cells > 0) {
      num.cells <- Matrix::rowSums(object@data > min.value)
      genes.use <- names(num.cells[which(num.cells >= min.cells)])
      object@data <- object@data[genes.use, ]
      return(object)
    } else {
      return(object)
    }
  }

# write a function to loop through all data files ####

calc_PICALM_type <-
  function(x) {

    Immune_cells_raw <-
      readRDS(x)
    # Immune_cells_raw <-
    #   readRDS(all_RDS_files[5])
    metadata.raw <-
      Immune_cells_raw@meta.data

    Immune_cells <-
      CreateSeuratObject(counts = Immune_cells_raw@assays$RNA@counts,
                         assay = "RNA",
                         meta.data = metadata.raw,
                         min.cells = 500,
                         min.features = 500)


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


    combined_Immune_cells_meta_Sun_meta <-
      merge(x = Immune_cells_meta,
            y = Sun_meta,
            by.x = "trimmed_barcode",
            by.y = "barcode")
    print(length(unique(combined_Immune_cells_meta_Sun_meta$projid)))

    Sun_immune_lookup_table <-
      combined_Immune_cells_meta_Sun_meta[!duplicated(combined_Immune_cells_meta_Sun_meta$projid), ]
    Sun_immune_lookup_table <-
      Sun_immune_lookup_table[, c("projid",
                                  "subject",
                                  "ADdiag3types")]
    # Sun_immune_lookup_table <-
    #   Sun_immune_lookup_table[, c(2, 5, 15)]

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

    # Immune_cells_noNA <-
    #   SCTransform(Immune_cells_noNA,
    #               vst.flavor = "v2",
    #               verbose = T)

    Immune_cells_noNA <-
      NormalizeData(Immune_cells_noNA,
                    normalization.method = "LogNormalize",
                    scale.factor = 10000)

    # Immune_cells_noNA <-
    #   FilterGenes(object = Immune_cells_noNA,
    #               min.value = 1,
    #               min.cells = 100)

    Idents(Immune_cells_noNA) <- "projid"

    # Immune_pseudobulk <-
    #   PseudobulkExpression(Immune_cells_noNA,
    #                        assays = "RNA",
    #                        return.seurat = F,
    #                        group.by = "ident",
    #                        normalization.method = "LogNormalize")
    Immune_pseudobulk <-
      AggregateExpression(Immune_cells_noNA,
                          assays = "RNA",
                          return.seurat = T,
                          group.by = c("projid",
                                       "ADdiag3types"))

    # Immune_pseudobulk$projid.ADdiag3types <-
    #   paste(Immune_pseudobulk$projid,
    #         Immune_pseudobulk$ADdiag3types,
    #         sep = "-")
    #
    Idents(Immune_pseudobulk) <- "ADdiag3types"

    bulk.mono.de.lateAD <-
      FindMarkers(object = Immune_pseudobulk,
                  ident.1 = "lateAD",
                  ident.2 = "nonAD",
                  test.use = "DESeq2")

    # bulk.mono.de.lateAD[rownames(bulk.mono.de.lateAD) == "PICALM", ]

    bulk.mono.de.earlyAD <-
      FindMarkers(object = Immune_pseudobulk,
                  ident.1 = "earlyAD",
                  ident.2 = "nonAD",
                  test.use = "DESeq2")

    return(list(bulk.mono.de.earlyAD,
                bulk.mono.de.lateAD))

    # bulk.mono.de.earlyAD[rownames(bulk.mono.de.earlyAD) == "PICALM", ]

    # Immune_pseudobulk <-
    #   Immune_pseudobulk$RNA
    # Immune_pseudobulk_PICALM <-
    #   data.frame(expression = Immune_pseudobulk[rownames(Immune_pseudobulk) == "PICALM", ])
    # Immune_pseudobulk_PICALM$g_projid <-
    #   rownames(Immune_pseudobulk_PICALM)
    #
    # Immune_lookup_table$g_projid <-
    #   str_c("g",
    #         Immune_lookup_table$projid,
    #         sep = "")
    # Sun_immune_lookup_table$g_projid <-
    #   str_c("g",
    #         Sun_immune_lookup_table$projid,
    #         sep = "")
    #
    # Immune_pseudobulk_PICALM <-
    #   merge(x = Immune_pseudobulk_PICALM,
    #         y = Sun_immune_lookup_table,
    #         by = "g_projid")

    # return(list(nonADvslateAD = wilcox.test(x = Immune_pseudobulk_PICALM$expression[Immune_pseudobulk_PICALM$ADdiag3types == "nonAD"],
    #                         y = Immune_pseudobulk_PICALM$expression[Immune_pseudobulk_PICALM$ADdiag3types == "lateAD"],
    #                         paired = F,
    #                         alternative = "t"),
    #
    #             nonADvsearlyAD = wilcox.test(x = Immune_pseudobulk_PICALM$expression[Immune_pseudobulk_PICALM$ADdiag3types == "nonAD"],
    #                         y = Immune_pseudobulk_PICALM$expression[Immune_pseudobulk_PICALM$ADdiag3types == "earlyAD"],
    #                         paired = F,
    #                         alternative = "t"),
    #
    #             mean_nonAD = mean(Immune_pseudobulk_PICALM$expression[Immune_pseudobulk_PICALM$ADdiag3types == "nonAD"]),
    #             mean_lateAD = mean(Immune_pseudobulk_PICALM$expression[Immune_pseudobulk_PICALM$ADdiag3types == "lateAD"]),
    #             mean_earlyAD = mean(Immune_pseudobulk_PICALM$expression[Immune_pseudobulk_PICALM$ADdiag3types == "earlyAD"])))
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

for (i in 1:length(all_RDS_files)) {
  print(all_RDS_files[i])

  try({
    PICALM_results[[i]] <-
      calc_PICALM_type(x = all_RDS_files[i])
  })

  # PICALM_results
}

saveRDS(PICALM_results,
        file = "PICALM_results_cell_types_1_2.RDs")


# non-function version for one-off calculation ####

calc_p_and_plot <-
  function(x) {
    Immune_cells_raw <-
      readRDS(all_RDS_files[5])
    # print(x)
    # Immune_cells_raw <-
    #   Immune_cells_raw[, Immune_cells_raw$cell_type_high_resolution == 'Mic P2RY12']
    #   readRDS(all_RDS_files[5])
    # metadata.raw <-
    #   Immune_cells_raw@meta.data
    #
    # Immune_cells <-
    #   CreateSeuratObject(counts = Immune_cells_raw@assays$RNA@counts,
    #                      assay = "RNA",
    #                      meta.data = metadata.raw,
    #                      min.cells = 1,
    #                      min.features = 1)
    Immune_cells <-
      Immune_cells_raw


    Immune_cells_meta <-
      Immune_cells@meta.data

    Sun_meta_w_projid <-
      readRDS("Sun_PFC_metadata_info_14Jan2025.RDs")
    # Immune_cells_meta$barcode <-
    #   rownames(Immune_cells_meta)
    # print(unique(str_split(Immune_cells_meta$barcode,
    #                        pattern = '-',
    #                        simplify = T)[, 3]))
    # print(length(unique(Immune_cells_meta$projid)))
    # Immune_cells_meta$trimmed_barcode <-
    #   str_c(str_split(string = Immune_cells_meta$barcode,
    #                   pattern = '-',
    #                   simplify = T)[, 1],
    #         str_split(string = Immune_cells_meta$barcode,
    #                   pattern = '-',
    #                   simplify = T)[, 2],
    #         sep = '-')

    # Immune_cells_cluster_0 <-
    #   Immune_cells_meta[str_detect(string = Immune_cells_meta$barcode,
    #                                pattern = "\\-0$"), ]
    # unique(Immune_cells_cluster_0$projid)



    # Immune_lookup_table <-
    #   Immune_cells_meta[!duplicated(Immune_cells_meta$projid), ]

    # Sun_meta <-
    #   readRDS("sun_et_al/personal.broadinstitute.org/cboix/sun_victor_et_al_data/ROSMAP.ImmuneCells.6regions.snRNAseq.meta.rds")
    # print(length(unique(Sun_meta$subject)))
    # Sun_meta <-
    #   Sun_meta[Sun_meta$brainRegion == "PFC", ]
    #

    # combined_Immune_cells_meta_Sun_meta <-
    #   merge(x = Immune_cells_meta,
    #         y = Sun_meta,
    #         by.x = "trimmed_barcode",
    #         by.y = "barcode")
    # print(length(unique(combined_Immune_cells_meta_Sun_meta$projid)))
    #
    # Sun_immune_lookup_table <-
    #   combined_Immune_cells_meta_Sun_meta[!duplicated(combined_Immune_cells_meta_Sun_meta$projid), ]
    # Sun_immune_lookup_table <-
    #   Sun_immune_lookup_table[, c("projid",
    #                               "subject",
    #                               "ADdiag3types")]
    # # Sun_immune_lookup_table <-
    #   Sun_immune_lookup_table[, c(2, 5, 15)]
    # Sun_immune_lookup_table <-
    #   readRDS("Sun_immune_lookup_table.RDs")

    Immune_cells_meta_reassign <-
      merge(x = Immune_cells_meta,
            y = Sun_meta_w_projid,
            by.x = "projid",
            by.y = "projid",
            all.x = T)
    Immune_cells_meta_reassign <-
      merge(x = Immune_cells_meta_reassign,
            y = rs10792832_ROSMAP,
            by = "projid",
            all.x = T)

    Immune_cells <-
      Immune_cells[, order(Immune_cells$projid)]
    Immune_cells_meta_reassign <-
      Immune_cells_meta_reassign[order(Immune_cells_meta_reassign$projid), ]

    #
    Immune_cells$ADdiag3types <-
      Immune_cells_meta_reassign$ADdiag3types
    Immune_cells$IndividualID <-
      Immune_cells_meta_reassign$individualID
    Immune_cells$rs10792832_genotype <-
      Immune_cells_meta_reassign$value

    # Immune_cells <-
    #   Immune_cells[, Immune_cells$projid %in% Sun_samples]


    print('=== 3 === ')

    Immune_cells <-
      UpdateSeuratObject(Immune_cells)

    Immune_cells_noNA <-
      Immune_cells[, !is.na(Immune_cells$ADdiag3types)]
    Immune_cells_noNA <-
      Immune_cells[, !is.na(Immune_cells$IndividualID)]
    Immune_cells_noNA <-
      Immune_cells_noNA[, Immune_cells_noNA$cell_type_high_resolution == 'Mic P2RY12']

    print('=== assigned ADdiag3types === ')

    # Immune_cells_noNA <-
    #   SCTransform(Immune_cells_noNA,
    #               vst.flavor = "v2",
    #               verbose = T)

    # Immune_cells_noNA <-
    #   NormalizeData(Immune_cells_noNA,
    #                 normalization.method = "LogNormalize",
    #                 scale.factor = 10000)

    # Immune_cells_noNA <-
    #   FilterGenes(object = Immune_cells_noNA,
    #               min.value = 1,
    #               min.cells = 100)

    Idents(Immune_cells_noNA) <- "projid"

    # Immune_pseudobulk <-
    #   PseudobulkExpression(Immune_cells_noNA,
    #                        assays = "RNA",
    #                        return.seurat = F,
    #                        group.by = "ident",
    #                        normalization.method = "LogNormalize")
    # Immune_pseudobulk <-
    #   AggregateExpression(Immune_cells_noNA,
    #                       assays = "RNA",
    #                       return.seurat = T,
    #                       group.by = c("projid"))
    #
    # mtx_Immune_pseudobulk <-
    #   as.matrix(Immune_pseudobulk@assays$RNA@layers$counts)
    # Immune_cells_noNA <-
    #   NormalizeData(Immune_cells_noNA)

    Immune_pseudobulk <-
      AggregateExpression(Immune_cells_noNA,
                           assays = "RNA",
                           return.seurat = T,
                           group.by = "projid",
                           normalization.method = "LogNormalize")

    mtx_Immune_pseudobulk <-
      as.matrix(Immune_pseudobulk@assays$RNA@layers$counts)




    rownames(mtx_Immune_pseudobulk) <-
      rownames(Immune_pseudobulk$RNA)
    colnames(mtx_Immune_pseudobulk) <-
      colnames(Immune_pseudobulk$RNA)

    meta_Immune <-
      Immune_cells_noNA@meta.data
    meta_Immune <-
      meta_Immune[!duplicated(meta_Immune$projid), ]

    mtx_immune_index <-
      data.frame(projid = colnames(mtx_Immune_pseudobulk))

    meta_Immune$g_projid <-
      paste0('g',
             meta_Immune$projid)
    mtx_immune_index <-
      merge(mtx_immune_index,
            meta_Immune,
            by.x = "projid",
            by.y = "g_projid")
    mtx_immune_index <-
      mtx_immune_index[order(mtx_immune_index$projid), ]

    mtx_Immune_pseudobulk <-
      mtx_Immune_pseudobulk[, order(colnames(mtx_Immune_pseudobulk))]


    mtx_cpm_Immune_pseudobulk <-
      DGEList(counts = mtx_Immune_pseudobulk,
              genes = rownames(mtx_Immune_pseudobulk),
              samples = colnames(mtx_Immune_pseudobulk),
              group = mtx_immune_index$ADdiag3types)





    cpm_output <-
      cpm(mtx_cpm_Immune_pseudobulk, log = T)
    cpm_output[rownames(cpm_output) == "PICALM", ]
    df_2_plot <-
      data.frame(exp = cpm_output[rownames(cpm_output) == "PICALM", ])
    df_2_plot$projid <-
      rownames(df_2_plot)
    # df_2_plot$ADdiag3types <-
    #   mtx_immune_index$ADdiag3types
    df_2_plot <-
      merge(df_2_plot,
            mtx_immune_index,
            by.x = "projid",
            by.y = "projid")
    df_2_plot <-
      df_2_plot[!is.na(df_2_plot$rs10792832_genotype), ]

    # scatter_plot <-
      ggplot(df_2_plot,
             aes(x = ADdiag3types,
                 y = exp)) +
      geom_boxplot(outliers = F) +
      # geom_dotplot(binaxis = "y",
      #              stackdir = "center", binwidth = 0.2,
      #              dotsize = 0.3,
      #              width = 0.3) +
      geom_jitter(width = 0.2) +
      theme_classic() #+
      # ggtitle(basename(x))

    earlyAD_vs_nonAD <-
      wilcox.test(x = df_2_plot$exp[df_2_plot$ADdiag3types == "earlyAD"],
             y = df_2_plot$exp[df_2_plot$ADdiag3types == "nonAD"],
             # var.equal = T,
             paired = F, alternative = "t")
    lateAD_vs_nonAD <-
      t.test(x = df_2_plot$exp[df_2_plot$ADdiag3types == "lateAD"],
             y = df_2_plot$exp[df_2_plot$ADdiag3types == "nonAD"],
             var.equal = T,
             paired = F, alternative = "t")

    return(list(df = df_2_plot,
                scatter_plot = scatter_plot,
                earlyAD_vs_nonAD = earlyAD_vs_nonAD,
                lateAD_vs_nonAD = lateAD_vs_nonAD))
  }

write.table(df_2_plot,
            file = "Sun_184_samples_w_diag_14Jan2025.tsv",
            row.names = F, col.names = T,
            quote = F, sep = "\t")



#
# calc ####
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

# Seurat_sun_PFC_cells_185_genotypes <-
#   readRDS("Seurat_sun_PFC_cells_185_genotypes.RDs")

rs10792832_ROSMAP <-
  read.delim("ROSMAP_genotypes_hg19/ROSMAP_sun_rs10792832_genotype.txt")
rs10792832_ROSMAP$projid <-
  as.numeric(str_remove(string = rs10792832_ROSMAP$ROSprojid,
                        pattern = "ROS"))
rs10792832_ROSMAP <-
  rs10792832_ROSMAP[!is.na(rs10792832_ROSMAP$projid), ]

Sun_samples <-
  unique(Seurat_sun_PFC_cells_185_genotypes$projid)
Sun_samples <-
  as.numeric(str_remove(Sun_samples,
                        pattern = "ROS"))

PICALM_pseudo_results <-
  vector(mode = "list",
         length = length(all_RDS_files))
names(PICALM_pseudo_results) <-
  names(all_RDS_files)

for (i in 1:(length(all_RDS_files) - 1)) {
  print(all_RDS_files[i])

  try({
    PICALM_pseudo_results[[i]] <-
      calc_p_and_plot(x = all_RDS_files[i])
  })
}



microglia_result <-
  calc_p_and_plot(x = all_RDS_files[5])


saveRDS(meta_Immune,
        file = "meta_indiv_index_from_immune.RDs")
saveRDS(Sun_immune_lookup_table,
        file = "Sun_immune_lookup_table.RDs")





View(Seurat_sun_PFC_cells_185_genotypes@meta.data)


print(PICALM_pseudo_results[[1]][["scatter_plot"]])
PICALM_pseudo_results[[1]][["earlyAD_vs_nonAD"]]
PICALM_pseudo_results[[1]][["lateAD_vs_nonAD"]]

print(PICALM_pseudo_results[[2]][["scatter_plot"]])
PICALM_pseudo_results[[2]][["earlyAD_vs_nonAD"]]
PICALM_pseudo_results[[2]][["lateAD_vs_nonAD"]]




print(PICALM_pseudo_results[[3]][["scatter_plot"]])
print(PICALM_pseudo_results[[4]][["scatter_plot"]])
print(PICALM_pseudo_results[[5]][["scatter_plot"]])
PICALM_pseudo_results[[5]][["earlyAD_vs_nonAD"]]
PICALM_pseudo_results[[5]][["lateAD_vs_nonAD"]]


print(PICALM_pseudo_results[[6]][["scatter_plot"]])
print(PICALM_pseudo_results[[7]][["scatter_plot"]])
print(PICALM_pseudo_results[[8]][["scatter_plot"]])
