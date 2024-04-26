# Siwei 21 Oct 2023
# Transform raw input, calculate P val from 95% CI, and plot SE+P

# init #####
library(ggplot2)
library(RColorBrewer)
library(readr)
library(stringr)
library(dplyr)

# load data #####
# df_raw <-
#   read_table("sum_2_plot_b4_transform.txt",
#              col_names = FALSE)

df_raw <-
  read_table("R_GWAS_pre_transform_23Oct2023.txt",
             col_names = FALSE)

df_2_plot <- df_raw

colnames(df_2_plot) <-
  c("Disease", "SNP_category",
    "lnOR_Estimate", "lnCI95low", "lnCI95High")

df_2_plot$Disease <-
  str_split(string = df_2_plot$Disease,
            pattern = '/',
            simplify = T)[, 2]

df_2_plot$SNP_category <-
  str_split(string = df_2_plot$SNP_category,
            pattern = "\\.",
            simplify = T)[, 1]
# df_2_plot <-
#   data.frame(Disease = df_raw$Disease,
#              SNP_category = df_raw$SNP_category,
#              log2OR = log2(exp(df_raw$lnOR_Estimate)),
#              `-log10P` = 0 - log10(exp((exp(df_raw$lnCI95High) - exp(df_raw$lnCI95low)) /
#                                               (2 * 1.96) *
#                                               (0 - 0.717) -
#                                               0.416 *
#                                               ((exp(df_raw$lnCI95High) - exp(df_raw$lnCI95low)) /
#                                                               (2 * 1.96)) ^ 2)))

z <- df_2_plot$lnOR_Estimate / ((df_2_plot$lnCI95High - df_2_plot$lnCI95low) / (2 * 1.96))
# p  <- exp(z * (-0.717) - 0.416 * z ^ 2)

df_2_plot_original <-
  data.frame(Disease = df_2_plot$Disease,
             SNP_category = df_2_plot$SNP_category,
             log2OR = log2(exp(df_2_plot$lnOR_Estimate)),
             `neglog10P` = 0 - log10(exp(z * (-0.717) - 0.416 * z ^ 2)))

df_2_plot_bubble <- df_2_plot_original
## remove undesired rows (diseases) first
# [1] "Alz_sumstats_Jansenetal_2019_Dec"
# [2] "daner_ADHD_2019"
# [3] "EUR.CD.GWAS"
# [4] "EUR.UC.GWAS"
# [5] "GIANT_Height_2014"
# [6] "Intelligence_Jansenetal_2020"
# [7] "iPSYCH-PGC_ASD_2019_hg19_wo_blockinfo_rsid"
# [8] "MDD_jamapsy_Giannakopoulou_2021"
# [9] "meta_result_EPICInterAct_T2D_GWAS_summary_4_ldsc_corrected"
# [10] "Neurotisim_Nagel_2018"
# [11] "Parkinson_nallsEtAl2019_excluding23andMe_allVariants_hg19_wo_blkinfo_rsid"
# [12] "PGC3_SCZ_w3_2021_w_blockinfo_no_rs"
# [13] "PGC_Alc_dep_EUR_2018"
# [14] "PGC_BIP_2021"
# [15] "PGC_UKB_depression_500199"
# [16] "PTSD_all_freeze2"
# [17] "T2D_GWAS_corrected"

df_2_plot_bubble <-
  df_2_plot_bubble[!(df_2_plot_bubble$Disease %in%
                        c("T2D_GWAS_corrected",
                          "MDD_jamapsy_Giannakopoulou_2021")), ]

## make Disease names consistent and easily readable
### make an 1:1 lookup table
replace_disease_lookup_table <-
  c("Alz_sumstats_Jansenetal_2019_Dec" = 'Alzheimers (2019)',
    "daner_ADHD_2019" = 'ADHD (2019)',
    "EUR.CD.GWAS" = 'Crohn\'s Disease (2015)',
    "EUR.UC.GWAS" = 'Ulcerative Colitis (2015)',
    "GIANT_Height_2014" = 'Height (2014)',
    "Intelligence_Jansenetal_2020" = 'Intelligence (2020)',
    "iPSYCH-PGC_ASD_2019_hg19_wo_blockinfo_rsid" = 'Asperger\'s Syndrome (2019)',
    "meta_result_EPICInterAct_T2D_GWAS_summary_4_ldsc_corrected" = 'Type 2 Diabetes (2012)',
    "Neurotisim_Nagel_2018" = 'Neuroticism (2018)',
    "Parkinson_nallsEtAl2019_excluding23andMe_allVariants_hg19_wo_blkinfo_rsid" = 'Parkinson\'s Disease (2019)',
    "PGC3_SCZ_w3_2021_w_blockinfo_no_rs" = 'Schizophrenia PGCw3 (2022)',
    "PGC_Alc_dep_EUR_2018" = 'Alcohol Dependence PGC3(2018)',
    "PGC_BIP_2021" = 'Bipolar PGC3 (2021)',
    "PGC_UKB_depression_500199" = 'MDD UKB (2018)',
    "PTSD_all_freeze2" = 'PTSD PGC2 (2019)')

df_2_plot_bubble$Disease <-
  str_replace_all(df_2_plot_bubble$Disease,
                  pattern = replace_disease_lookup_table)
df_2_plot_bubble$Disease <-
  factor(df_2_plot_bubble$Disease,
         levels = unique(df_2_plot_bubble$Disease)[c(11, 13, 14, 9,
                                                     7, 2, 6, 12, 15,
                                                     1, 10, 5, 3, 4, 8)])
# reorder $Disease by current index

# df_2_plot_bubble <-
#   df_2_plot_bubble %>%
#   arrange(Disease) %>%
#   mutate(Disease = factor(df_2_plot_bubble$Disease,
#                           levels = unique(df_2_plot_bubble$Disease)[c(11, 13, 14, 9,
#                                                                       7, 2, 6, 12, 15,
#                                                                       1, 10, 5, 3, 4, 8)]))



ggplot(df_2_plot_bubble,
       aes(x = Disease,
           y = SNP_category,
           size = neglog10P,
           fill = log2OR)) +
  scale_fill_gradient2(low = "darkblue",
                       mid = "white",
                       high = "darkred") +
  scale_radius() +
  geom_point(shape = 21) +
  scale_y_discrete(limits = rev) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 315,
                                   vjust = 0,
                                   hjust = 0))


df_2_plot_bubble2 <-
  df_2_plot_bubble[str_detect(string = df_2_plot_bubble$SNP_category,
                              pattern = "_[016]hr")
                   , ]
df_2_plot_bubble2 <-
  df_2_plot_bubble2[str_detect(string = df_2_plot_bubble2$SNP_category,
                              pattern = "unique",
                              negate = T)
                   , ]
ggplot(df_2_plot_bubble2,
       aes(x = Disease,
           y = SNP_category,
           size = neglog10P,
           fill = log2OR)) +
  scale_fill_gradient2(low = "darkblue",
                       mid = "white",
                       high = "darkred",
                       limits = c(-10, 10)) +
  scale_radius() +
  geom_point(shape = 21,
             stroke = 0.5) +
  scale_y_discrete(limits = rev) +
  labs(size = '-log10Pval') +
  ylab('SNP Category') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 315,
                                   vjust = 0,
                                   hjust = 0))


df_2_plot_bubble2 <-
  df_2_plot_bubble[str_detect(string = df_2_plot_bubble$SNP_category,
                              pattern = "_[016]hr")
                   , ]
df_2_plot_bubble2 <-
  df_2_plot_bubble2[str_detect(string = df_2_plot_bubble2$SNP_category,
                               pattern = "unique",
                               negate = F)
                    , ]
ggplot(df_2_plot_bubble2,
       aes(x = Disease,
           y = SNP_category,
           size = neglog10P,
           fill = log2OR)) +
  scale_fill_gradient2(low = "darkblue",
                       mid = "white",
                       high = "darkred",
                       limits = c(-10, 10)) +
  scale_radius() +
  geom_point(shape = 21) +
  scale_y_discrete(limits = rev) +
  labs(size = '-log10Pval') +
  ylab('SNP Category') +
  theme_bw() +
  # theme_classic() +
  theme(axis.text.x = element_text(angle = 315,
                                   vjust = 0,
                                   hjust = 0))
