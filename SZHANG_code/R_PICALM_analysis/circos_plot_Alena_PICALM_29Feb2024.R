# Siwei 29 Feb 2024
# Make circos plots for Alena's PICALM paper

# init ####
library(readxl)
library(circlize)

library(RColorBrewer)
library(colorRamps)

library(viridis)
library(colorspace)
library(colorRamps)

library(stringr)
library(ggplot2)
# library(RColorBrewer)

## make line-dot plot #####
df_pathways_2_genes <-
  read_excel("IPA_pathway_for_graphing.xlsx",
             sheet = 3)
df_gene_log2FC_raw <-
  read_excel("IPA_pathway_for_graphing.xlsx",
             sheet = 4)

df_pathways_2_genes <-
  df_pathways_2_genes[1:4, ]
df_pathways_2_genes <-
  df_pathways_2_genes[order(df_pathways_2_genes$`-log_pVal`), ]
df_pathways_2_genes$Ingenuity_Canonical_Pathways <-
  c("P_Macro",
    "P_Act_SREBF",
    "P_Chl",
    "P_MHC")
df_extended <-
  lapply(1:nrow(df_pathways_2_genes),
        # margin = 1,
        FUN = function(x) {
          data.frame(Pathways = df_pathways_2_genes[x, 1],
                     `-logPval` = df_pathways_2_genes[x, 2],
                     Genes = unlist(str_split(df_pathways_2_genes[x, 5],
                                              pattern = ",",
                                              simplify = F)))
        })
df_extended <-
  do.call(what = rbind,
          args = df_extended)
# df_extended <-
#   df_extended[order(df_extended$X.log_pVal)]

unique_gene_list <-
  unique(df_extended$Genes)
df_gene_log2FC <-
  df_gene_log2FC_raw[df_gene_log2FC_raw$Gene %in% unique_gene_list, ]
df_gene_log2FC <-
  df_gene_log2FC[order(df_gene_log2FC$logFC), ]

df_2_plot <-
  merge(x = df_extended,
        y = df_gene_log2FC,
        by.x = "Genes",
        by.y = "Gene",
        all.x = T)

# df_2_plot <-
#   df_2_plot[, c(1,2,4,3)]

group_order <-
  c(df_gene_log2FC$Gene[order(df_gene_log2FC$logFC)],
    df_pathways_2_genes$Ingenuity_Canonical_Pathways[order(df_pathways_2_genes$`-log_pVal`)])
# group_order_vector <-
#   group_order
group_order_vector <-
  c(rep_len("Genes",
            length.out = (length(group_order) - 4)),
    rep_len("Pathways",
            length.out = 4))
names(group_order_vector) <-
  group_order
group_order_vector <-
  factor(group_order_vector,
         levels = unique(group_order_vector))

df_2_plot$width_genes <-
  360 / 3 * 2 /
  (length(group_order_vector) - 4)
names(table(df_2_plot$Genes))
names(table(df_2_plot$Genes))[table(df_2_plot$Genes) == 2]
df_2_plot$width_genes[df_2_plot$Genes %in%
                        names(table(df_2_plot$Genes))[table(df_2_plot$Genes) == 2]] <-
  360 / 3 * 1 /
  (length(group_order_vector) - 4)

df_2_plot$width_pathways <- 5
  # 360 / 3  / 4

df_2_plot <-
  df_2_plot[, c(1,2,5,6,3,4)]

map2color <-
  function(x, palette, limits=NULL) {
  if (is.null(limits)) limits = range(x)
  palette[findInterval(x,
                       seq(limits[1],
                           limits[2],
                           length.out = length(palette) + 1),
                       all.inside = TRUE)]
  }

# map2color(x = sort(unique(df_2_plot$logFC)),
#           palette = rev(colorspace::heat_hcl(n = 100,
#                                              h = c(0, -100),
#                                              l = c(75, 40),
#                                              c = c(40, 80),
#                                              power = 1,
#                                              fixup = T)))
# map2color(x = sort(unique(df_2_plot$X.log_pVal)),
#           palette = rev(magenta2green(n = 100)))
range(df_2_plot$logFC)
log(range(df_2_plot$logFC)[2]) /
  log(range(df_2_plot$logFC)[2] - range(df_2_plot$logFC)[1])

grid_colour <-
  c(map2color(x = sort(unique(df_2_plot$logFC)),
              palette = colorRampPalette(colors = c("darkblue",
                                                    "white",
                                                    "darkred"),
                                         bias = 1)(100)),
    map2color(x = sort(unique(df_2_plot$X.log_pVal)),
              palette = viridis(n = 100)))
# grid_colour <-
#   c(map2color(x = sort(unique(df_2_plot$logFC)),
#               palette = rev(colorspace::heat_hcl(n = 100,
#                                                  h = c(0, -100),
#                                                  l = c(75, 40),
#                                                  c = c(40, 80),
#                                                  power = 1,
#                                                  fixup = T))),
#     map2color(x = sort(unique(df_2_plot$X.log_pVal)),
#               palette = viridis(n = 100)))
names(grid_colour) <-
  group_order
#
# cut(x = unique(df_2_plot$X.log_pVal),
#     breaks = 10)
# range(unique(df_2_plot$X.log_pVal))
#
# map_range <-
#   function(x) {
#     (x - range(x, na.rm = T)[1]) /
#           (range(x, na.rm = T)[2] - range(x, na.rm = T)[1])
#   }
# map_range(unique(df_2_plot$X.log_pVal))
# map_range(unique(df_2_plot$logFC))
#
#
#
# diff(unique(df_2_plot$X.log_pVal))
# colorRampPalette(colors = c("blue",
#                             "red"))
# grid_col_vector <-
#   c(rep_len(x = (length(group_order_vector) / 3 / 4),
#             length.out = 4),
#     rep_len(x = (length(group_order_vector) / 3 * 2 /
#                    (length(group_order_vector) - 4)),
#             length.out = length(group_order_vector) - 4))
# grid_col_vector <-
#   c(rep_len(x = (360 / 3 / 4),
#             length.out = 4),
#     rep_len(x = (360 / 3 * 2 /
#                    (length(group_order_vector) - 4)),
#             length.out = length(group_order_vector) - 4))
# names(grid_col_vector) <-
#   names(group_order_vector)

# circos.info()
circos.clear()
# par(mfrow = c(1,2))
circos.par(start.degree = 0)
chordDiagramFromDataFrame(df = df_2_plot[, 1:4],
                          annotationTrack = c("grid"),
                          group = group_order_vector,
                          grid.col = grid_colour,
                          annotationTrackHeight = 0.05,
                          preAllocateTracks = list(track.height = 0.05))
## customise sector labels
circos.track(track.index = 1,
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter,
                           CELL_META$ylim[1],
                           CELL_META$sector.index,
                           facing = "clockwise",
                           niceFacing = T,
                           adj = c(0, 0),
                           cex = 0.7)
             },
             bg.border = NA)

# chordDiagramFromDataFrame(df = df_2_plot[, 1:4],
#                           annotationTrack = c("name",
#                                               "grid"),
#                           group = group_order_vector,
#                           grid.col = grid_colour)
circos.clear()
# circos.info()
# df_2_plot <-
#   df_2_plot[, c(1,2,4,3,5,6)]

## plot colour scales ####
df_color_scale <-
  data.frame(x = 1:100,
             y = seq(from = range(df_2_plot$logFC)[1],
                     to = range(df_2_plot$logFC)[2],
                     length.out = 100),
             fill = colorRampPalette(colors = c("darkblue",
                                                "white",
                                                "darkred"),
                                     bias = 1)(100))

# ggplot(df_color_scale,
#        aes(x = x)) +
#   geom_bar(fill = df_color_scale$fill) +
#   theme_minimal() +
#   theme(legend.position = "none")

seq(from = range(df_2_plot$logFC)[1],
    to = range(df_2_plot$logFC)[2],
    length.out = 100)
# [1] -2.132579  1.993261

ggplot(df_color_scale) +
  geom_tile(aes(x = 1,
                y = y),
            fill = df_color_scale$fill) +
  ylab("logFC") +
  coord_flip() +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10,
                                   colour = "black"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())


df_color_scale <-
  data.frame(x = 1:100,
             y = seq(from = range(df_2_plot$X.log_pVal)[1],
                     to = range(df_2_plot$X.log_pVal)[2],
                     length.out = 100),
             fill = viridis(n = 100))

ggplot(df_color_scale) +
  geom_tile(aes(x = 1,
                y = y),
            fill = df_color_scale$fill) +
  ylab("logFC") +
  coord_flip() +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10,
                                   colour = "black"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())
