# Siwei 26 Jan 2025
# plot flow cyto
# need to downsample each subset to 5800

# init ####
{
  library(readxl)

  library(stringr)
  library(ggplot2)

  library(scales)
  library(reshape2)

  library(RColorBrewer)
  library(ggpubr)
  library(ggridges)

  library(dplyr)
  library(data.table)

  library(DescTools)
  library(multcomp)

  library(gridExtra)

  library(ggallin)
}

set.seed(42)


csv_list <-
  list.files(path = "Batch_2_of_new_plots/Alena-iMG-sub_population_export",
             pattern = ".*.csv",
             full.names = T)

raw_data_list <-
  vector(mode = "list",
         length = length(csv_list))

for (i in 1:length(csv_list)) {
  df_input <-
    read.csv(csv_list[i])
  # FITC_val <-
  #   df_input$FITC.A
  FITC_min <-
    min(df_input$FITC.A,
        df_input = T)
  print(paste(min(df_input$FITC.A),
              max(df_input$FITC.A)))
  # df_input$FITC.A <-
  #   df_input$FITC.A +
  #   FITC_min + 0.1
  raw_data_list[[i]] <-
    df_input
  # raw_data_list[[i]] <-
  #   df_input[df_input$FITC.A > 0, ]
}

# hist(raw_data_list[[1]]$FITC.A,breaks = 100,
#      freq = F)
# hist(raw_data_list[[2]]$FITC.A,breaks = 100,
#      freq = F)

names(raw_data_list) <-
  basename(path = csv_list)

# sum(raw_data_list[[3]]$FITC.A > 0)



subsetted_list <-
  vector(mode = "list",
         length = length(csv_list))
names(subsetted_list) <-
  names(raw_data_list)

for (i in 1:length(raw_data_list)) {
  # subsetted_list[[i]] <-
  #   raw_data_list[[i]][sample(x = nrow(raw_data_list[[i]]),
  #                             size = 5600,
  #                             replace = F), ]
  subsetted_list[[i]] <-
    raw_data_list[[i]]
  subsetted_list[[i]]$sample_name <-
    names(raw_data_list)[i]
}



df_2_plot <-
  do.call(what = "rbind",
          args = subsetted_list)
min(df_2_plot$FITC.A)
max(df_2_plot$FITC.A)

# reverselog_trans <-
#   function(base = exp(1)) {
#   trans <- function(x) -log(x, base)
#   inv <- function(x) base ^(-x)
#   trans_new(paste0("reverselog-", format(base)), trans, inv,
#             log_breaks(base = base),
#             domain = c(1e-100, Inf))
# }

# ggplot(df_2_plot,
#        aes(x = FITC.A + min(df_2_plot$FITC.A) + 1,
#            y = sample_name,
#            fill = sample_name)) +
#   geom_density_ridges2() +
#   scale_x_log10() +
#   # scale_x_continuous(transform = 'log1p') +
#   scale_fill_manual(values = brewer.pal(n = 4,
#                                         name = "Dark2")) +
#   theme_classic()
df_2_plot$cell_line <-
  str_split(string = df_2_plot$sample_name,
            pattern = "_",
            simplify = T)[, 3]
df_2_plot$genotype <-
  str_split(string = df_2_plot$sample_name,
            pattern = "_",
            simplify = T)[, 4]

df_2_plot$samples <-
  str_c(df_2_plot$cell_line,
        df_2_plot$genotype,
        sep = '-')

df_2_plot$samples <-
  factor(df_2_plot$samples,
         levels = rev(c(
                        "CD04-non-risk",
                        "CD04-risk",
                        "CD09-non-risk",
                        "CD09-risk")))

ggplot(df_2_plot,
       aes(x = FITC.A,
           y = samples,
           fill = samples)) +
  geom_density_ridges2(na.rm = T,
                       panel_scaling = F,
                       scale = 0.9) +
  # geom_density() +
  scale_x_continuous(transform = pseudolog10_trans,,
                     breaks = c(0, 10, 100, 1000, 10000, 100000),
                     labels = scientific(x = c(0, 10, 100, 1000, 10000, 100000)),
                     limits = c(0, 2e5),
                     expand = c(0, 0)) +
  # scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x))) +
  # scale_y_continuous(expand = expansion(mult = c(0, 0.01)),
  #                    limits = c(0, 1)) +
  # scale_x_continuous(labels = scientific) +
  # scale_fill_manual(values = brewer.pal(n = 4,
  #                                       name = "Dark2")) +
  xlab("FITC-A intensity") +
  ylab("") +

  scale_fill_manual(values = rev(c("darkblue",
                                   "darkred",
                                   "darkblue",
                                   "darkred"))) +
  annotation_logticks(sides = "b", scaled = T) +
  geom_vline(xintercept = 1000,
             linetype = 2) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 315,
                                   hjust = 0),
        axis.text = element_text(size = 12,
                                 colour = "black"),
        legend.position = "none")



  # facet_wrap(samples ~ .,
  #            strip.position = "y",
  #            ncol = 2)


ggplot(df_2_plot,
       aes(x = FITC.A,
           y = samples,
           # group = samples,
           fill = samples)) +
  geom_density_ridges2(na.rm = T,
                       panel_scaling = F,
                       scale = 0.9) +
  # geom_density() +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  # scale_y_continuous(expand = expansion(mult = c(0, 0.01)),
  #                    limits = c(0, 1)) +
  # scale_x_continuous(labels = scientific) +
  scale_fill_manual(values = brewer.pal(n = 4,
                                        name = "Dark2")) +
  annotation_logticks(sides = "b") +
  theme_classic() +
  ggtitle("Sub-sampled to 5600 cells/each")

# plot statistics ####

df_raw <-
  read_excel(path = "Batch_2_of_new_plots/Bodipy_sorting_data_only.xlsx")

## CD04 ####
df_2_plot <-
  df_raw[df_raw$Line == "CD04", ]
df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk"))

colnames(df_2_plot) <-
  make.names(colnames(df_2_plot))

df_2_plot$Differnetiations <-
  as.factor(df_2_plot$Differnetiations)

lm_model <-
  lmerTest::lmer(Bodipy..total.cells ~
                   Genotype +
                   (1|Differnetiations),
                 data = df_2_plot)
summary(lm_model)

ggerrorplot(df_2_plot,
            x = "Genotype",
            y = "Bodipy..total.cells",
            color = "Genotype",
            # shape = "Clones",
            # group = "Genotype",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkred",
                                 "darkblue"),
                      guide = "none") +
  stat_compare_means(label = "p.signif",
                     label.y.npc = 0.95,
                     # label.x.npc = 0,
                     # method = 't.test',
                     method = 'wilcox.test',
                     hide.ns = F,
                     ref.group = "non-risk",
                     paired = F) +
  geom_jitter(aes(colour = Genotype),
              width = 0.1,
              size = 1) +
  # scale_shape_manual(values = c(1, 2)) +
  stat_summary(geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black") +
  labs(x = "",
       y = "Bodypy..total.cells") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0,
                                max(df_raw$`Bodipy+/total cells`,
                                    na.rm = T) * 1.1)) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 12,
                                 colour = "black"),
        legend.position = "none",
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 1)) +
  ggtitle("CD04")

## CD09 ####
df_2_plot <-
  df_raw[df_raw$Line == "CD09", ]
df_2_plot$Genotype <-
  factor(df_2_plot$Genotype,
         levels = c("risk",
                    "non-risk"))
colnames(df_2_plot) <-
  make.names(colnames(df_2_plot))

df_2_plot$Differnetiations <-
  as.factor(df_2_plot$Differnetiations)

lm_model <-
  lmerTest::lmer(Bodipy..total.cells ~
                   Genotype +
                   (1|Differnetiations),
                 data = df_2_plot)
summary(lm_model)

ggerrorplot(df_2_plot,
            x = "Genotype",
            y = "Bodipy..total.cells",
            color = "Genotype",
            # shape = "Clones",
            # group = "Genotype",
            width = 0.2,
            # facet.by = "Gene",
            # ncol = 4,
            error.plot = "errorbar") +
  scale_colour_manual(values = c("darkred",
                                 "darkblue"),
                      guide = "none") +
  stat_compare_means(label = "p.signif",
                     label.y.npc = 0.95,
                     # label.x.npc = 0,
                     # method = 't.test',
                     method = 'wilcox.test',
                     hide.ns = F,
                     ref.group = "non-risk",
                     paired = F) +
  geom_jitter(aes(colour = Genotype),
              width = 0.1,
              size = 1) +
  # scale_shape_manual(values = c(1, 2)) +
  stat_summary(geom = "crossbar",
               fun = "mean",
               width = 0.5,
               linewidth = 0.2,
               colour = "black") +
  labs(x = "",
       y = "Bodypy..total.cells") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0,
                                max(df_raw$`Bodipy+/total cells`,
                                    na.rm = T) * 1.1)) +
  # stat_compare_means(label = "p.signif")
  theme_classic() +
  theme(axis.text = element_text(size = 12,
                                 colour = "black"),
        legend.position = "none",
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 1)) +
  ggtitle("CD09")
