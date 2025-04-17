
DEG_185_agg <-
  DGEList(counts = agg_185_samples,
          # genes = rownames(agg_185_samples),
          # samples = colnames(agg_185_samples),
          group = agg_185_metadata$rs561655_G_A,
          remove.zeros = T)

# table(filterByExpr(DEG_185_agg))

# DEG_185_agg <-
#   DEG_185_agg[filterByExpr(DEG_185_agg),
#               ,
#               keep.lib.sizes = F]
# sum(rownames(DEG_185_agg) %in% "PICALM")


DEG_185_agg <-
  calcNormFactors(DEG_185_agg)
cpm_target <-
  as.data.frame(edgeR::cpm(DEG_185_agg))

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
# rownames(cpm_target)[str_detect(string = rownames(cpm_target),
#                                 pattern = "APOE")]

## APOE
t_cpm <-
  as.data.frame(t(cpm_target))
# all(rownames(t_cpm) == rownames(agg_185_metadata))

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
           y = APOE,
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

save.image("workspace_4_mast_28Jun2024.RData")
q(save = F)
