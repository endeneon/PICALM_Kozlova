# Siwei 27 Sept 2023
# Make Jubao's table into bubble plot


# init
{
  library(readxl)
  library(ggplot2)
  library(RColorBrewer)
  library(stringr)
  library(scales)
}

df_raw <-
  read_excel("MAST_case_control_diff_Gene_to_R.xlsx")
df_raw <-
  df_raw[-(24:33), ]

# split into ups and downs, since downs require their OR set to negative

df_raw_positive <-
  as.data.frame(cbind(df_raw$Term,
                      df_raw[, str_detect(string = colnames(df_raw),
                                          pattern = "_up_")]))
df_raw_negative <-
  as.data.frame(cbind(df_raw$Term,
                      df_raw[, str_detect(string = colnames(df_raw),
                                          pattern = "_down_")]))
nrow(df_raw_positive)
nrow(df_raw_negative)

func_table_transform <-
  function(data) {
    for (i in seq(2, ncol(data) - 1, 2)) {
      print(i)
      if (i == 2) {
        return_df <-
          data[, 1:3]
        return_df$cell_type_time <-
          str_c(str_split(string = colnames(return_df)[2],
                          pattern = "_do|_up",
                          simplify = T)[, 1],
                "hr")
        colnames(return_df) <-
          c("Term", "FDR", "OR", "cell_type_time")
      } else {
        temp_df <-
          cbind(data[, 1],
                data[, i:(i + 1)])
        temp_df$cell_type_time <-
          str_c(str_split(string = colnames(temp_df)[2],
                          pattern = "_do|_up",
                          simplify = T)[, 1],
                "hr")
        print(colnames(temp_df))
        colnames(temp_df) <-
          c("Term", "FDR", "OR", "cell_type_time")
        print(colnames(temp_df))
        print(colnames(return_df))
        return_df <-
          as.data.frame(rbind(return_df,
                              temp_df))
      }
    }
    return(return_df)
  }

df_raw_pos_reorg <-
  func_table_transform(data = df_raw_positive)
# df_raw_pos_reorg$cell_type_time <-
df_raw_pos_reorg$direction <- "positive"

df_raw_neg_reorg <-
  func_table_transform(data = df_raw_negative)
df_raw_neg_reorg$OR <-
  0 - df_raw_neg_reorg$OR
df_raw_neg_reorg$direction <- "negative"

df_assembled <-
  as.data.frame(rbind(df_raw_pos_reorg,
                      df_raw_neg_reorg))

df_assembled$OR[is.na(df_assembled$OR)] <- 0
df_assembled$FDR[is.na(df_assembled$FDR)] <- 0


df_assembled_up <-
  df_assembled[df_assembled$direction %in% "positive", ]
df_assembled_down <-
  df_assembled[df_assembled$direction %in% "negative", ]

df_assembled_plot <-
  df_assembled_up
df_assembled_plot$FDR <-
  df_assembled_up$FDR +
  df_assembled_down$FDR
df_assembled_plot$OR <-
  df_assembled_up$OR +
  df_assembled_down$OR

df_assembled_plot$FDR[df_assembled_plot$FDR == 0] <- 1

df_assembled_plot$Term <-
  str_split(string = df_assembled_plot$Term,
            pattern = ' \\(',
            simplify = T)[, 1]
df_assembled_plot$Term <-
  str_replace_all(string = df_assembled_plot$Term,
                  pattern = "dependent\\ cotranslational\\ ",
                  replacement = "cotranslational\\ \n")
df_assembled_plot$Term <-
  str_replace_all(string = df_assembled_plot$Term,
                  pattern = ",\\ ",
                  replacement = ",\\ \n")

ggplot(data = df_assembled_plot,
       aes(x = Term,
           y = cell_type_time,
           size = abs(OR),
           alpha = 0 - log10(FDR),
           fill = factor(ifelse(test = OR > 0,
                                yes = "Up-regulated",
                                no = "Down-regulated")))) +
  geom_point(shape = 21) +
  scale_size_area() +
  scale_alpha_continuous() +
  # scale_fill_discrete() +
  scale_fill_manual(values = c("orangered2",
                               "darkblue"),
                    breaks = c("Up-regulated",
                               "Down-regulated")) +
  labs(alpha = '-log10FDR',
       size = "Odds Ratio",
       fill = "Direction") +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0),
        legend.position = "bottom")

ggplot(data = df_assembled_plot,
       aes(x = factor(Term,
                      levels = Term[1:30]),
           y = cell_type_time,
           size = abs(OR),
           alpha = 0 - log10(FDR),
           fill = factor(ifelse(test = OR > 0,
                                yes = "Up-regulated",
                                no = "Down-regulated")))) +
  geom_point(shape = 21) +
  scale_size_area() +
  scale_alpha_continuous(breaks = c(0, 10, 20)) +
  # scale_fill_discrete() +
  scale_fill_manual(values = c("orangered2",
                               "darkblue"),
                    breaks = c("Up-regulated",
                               "Down-regulated")) +
  labs(alpha = '-log10FDR',
       size = "Odds Ratio",
       fill = "Direction",
       x = "GO Terms",
       y = "Cell type + Time") +
  scale_y_discrete(limits = rev) +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "right")

###
df_plot_alt <-
  df_assembled_plot

df_plot_alt$logFDR <-
  0 - log10(df_plot_alt$FDR)
df_plot_alt$logFDR[df_plot_alt$OR < 0] <-
  0 - df_plot_alt$logFDR[df_plot_alt$OR < 0]

ggplot(data = df_plot_alt,
       aes(x = factor(Term,
                      levels = Term[1:30]),
           y = cell_type_time,
           size = abs(OR),
           # alpha = 0 - log10(FDR),
           fill = logFDR)) +
  geom_point(shape = 21) +
  scale_size_area() +
  scale_fill_gradientn(colours = c("darkblue",
                                   "white",
                                   "darkred"),
                       values = rescale(c(min(df_plot_alt$logFDR),
                                          0,
                                          max(df_plot_alt$logFDR))),
                       limits = c(min(df_plot_alt$logFDR),
                                  max(df_plot_alt$logFDR)
                                  )) +
  # scale_alpha_continuous(breaks = c(0, 10, 20)) +
  # scale_fill_gra(values = c("orangered",
  #                                  "white",
  #                                  "darkblue"),
  #                       breaks = c(min(df_plot_alt$logFDR),
  #                                  0,
  #                                  max(df_plot_alt$logFDR)))
  # scale_fill_discrete() +
  # scale_fill_manual(values = c("orangered2",
  #                              "darkblue"),
  #                   breaks = c("Up-regulated",
  #                              "Down-regulated")) +
  labs(fill = "FDR",
       size = "Odds Ratio",
       # fill = "Direction",
       x = "GO Terms",
       y = "Cell type + Time") +
  guides(fill = "none") +
  scale_y_discrete(limits = rev) +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "right")


df_plot_alt_sub <-
  df_plot_alt[df_plot_alt$logFDR < 0, ]
df_plot_alt_sub$logFDR <-
  0 - df_plot_alt_sub$logFDR

ggplot(data = df_plot_alt_sub,
       aes(x = factor(Term),
           y = cell_type_time,
           size = abs(OR),
           # alpha = 0 - log10(FDR),
           fill = logFDR)) +
  geom_point(shape = 21) +
  scale_size_area() +
  scale_fill_gradientn(colours = c("white",
                                   "darkred"),
                       # values = rescale(c(max(df_plot_alt_sub$logFDR),
                       #                    0)),
                       limits = c(0, max(df_plot_alt_sub$logFDR))) +
  # scale_alpha_continuous(breaks = c(0, 10, 20)) +
  # scale_fill_gra(values = c("orangered",
  #                                  "white",
  #                                  "darkblue"),
  #                       breaks = c(min(df_plot_alt$logFDR),
  #                                  0,
  #                                  max(df_plot_alt$logFDR)))
  # scale_fill_discrete() +
  # scale_fill_manual(values = c("orangered2",
  #                              "darkblue"),
  #                   breaks = c("Up-regulated",
#                              "Down-regulated")) +
labs(fill = "-log10FDR",
     size = "Odds Ratio",
     # fill = "Direction",
     x = "GO Terms",
     y = "Cell type + Time") +
  guides(size = "none") +
  scale_y_discrete(limits = rev) +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "right")

df_plot_alt_sub <-
  df_plot_alt[df_plot_alt$logFDR > 0, ]
# df_plot_alt_sub$logFDR <-
#   0 - df_plot_alt_sub$logFDR

ggplot(data = df_plot_alt_sub,
       aes(x = factor(Term),
           y = cell_type_time,
           size = abs(OR),
           # alpha = 0 - log10(FDR),
           fill = logFDR)) +
  geom_point(shape = 21) +
  scale_size_area() +
  scale_fill_gradientn(colours = c("white",
                                   "darkred"),
                       # values = rescale(c(max(df_plot_alt_sub$logFDR),
                       #                    0)),
                       limits = c(0, max(df_plot_alt_sub$logFDR))) +
  # scale_alpha_continuous(breaks = c(0, 10, 20)) +
  # scale_fill_gra(values = c("orangered",
  #                                  "white",
  #                                  "darkblue"),
  #                       breaks = c(min(df_plot_alt$logFDR),
  #                                  0,
  #                                  max(df_plot_alt$logFDR)))
  # scale_fill_discrete() +
  # scale_fill_manual(values = c("orangered2",
  #                              "darkblue"),
  #                   breaks = c("Up-regulated",
#                              "Down-regulated")) +
labs(fill = "-log10FDR",
     size = "Odds Ratio",
     # fill = "Direction",
     x = "GO Terms",
     y = "Cell type + Time") +
  guides(size = "none") +
  scale_y_discrete(limits = rev) +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "right")
