# Siwei 01 Mar 2024
# plot insert sumstats of 5 cell types

# init ####
library(readr)
library(ggplot2)

library(RColorBrewer)

library(stringr)
library(dplyr)

# load a file list of all fragment files
frag_file <-
  dir(path = "insert_sumstats",
      pattern = "*.txt",
      full.names = T,
      recursive = T)

# select 8 files from the cell type
cell_type_to_plot <- "MG"

file_list_to_readin <-
  frag_file[str_detect(string = frag_file,
                       pattern = cell_type_to_plot)]
file_list_to_readin
file_list_to_readin <-
  file_list_to_readin[str_detect(string = file_list_to_readin,
                                 pattern = "bwa_barcode_WASPed",
                                 negate = F)]
file_list_to_readin <-
  file_list_to_readin[order(file_list_to_readin)]
file_list_to_readin <-
  file_list_to_readin[str_detect(string = file_list_to_readin,
                                 pattern = "10x",
                                 negate = T)]
file_list_to_readin <-
  file_list_to_readin[str_detect(string = file_list_to_readin,
                                 pattern = "MG_",
                                 negate = T)]
# file_list_to_readin <-
#   file_list_to_readin[1:8]
file_list_to_readin
file_list_to_readin <-
  file_list_to_readin[sample(x = (length(file_list_to_readin) - 3),
                             size = 8)]

frag_list_2_bind <-
  vector(mode = "list",
         length = length(file_list_to_readin))
frag_list_2_bind <-
  lapply(X = 1:length(file_list_to_readin),
         FUN = function(x) {
           frag_df_raw <-
             read_delim(file_list_to_readin[x],
                        delim = "\t", escape_double = FALSE,
                        col_names = FALSE, trim_ws = TRUE)
           frag_df_raw$Sample <-
             file_list_to_readin[x]
           return(frag_df_raw)
         })

df_2_plot <-
  do.call(what = "rbind",
          args = frag_list_2_bind)
colnames(df_2_plot)[1:5] <-
  c("frag_size", "counts", "pos_strand", "neg_strand", "outward")
df_2_plot$Sample <-
  str_split(string = df_2_plot$Sample,
            pattern = "\\/",
            simplify = T)[, 2]
df_2_plot$Sample <-
  str_split(string = df_2_plot$Sample,
            pattern = "\\_",
            simplify = T)[, 1]

df_2_test <-
  df_2_plot %>%
  group_by(Sample) %>%
  summarise(counts = counts,
            frac_counts = counts / sum(counts),
            frag_size = frag_size)
# df_2_test$frag_size <-
#   factor(df_2_test$frag_size,
#          levels = sort(unique(df_2_test$frag_size)))
df_2_test$Sample <-
  factor(df_2_test$Sample,
         levels = unique(df_2_test$Sample))

# ggplot(df_2_test,
#        aes(x = frag_size,
#            y = log(frac_counts),
#            group = Sample,
#            colour = Sample)) +
#   geom_line() +
#   # geom_smooth(formula = log(y) ~ x,
#   #             method = "loess",
#   #             se = F) +
#   scale_colour_manual(values = brewer.pal(n = 8,
#                                           name = "Dark2")) +
#   ylim(-8, -4) +
#   theme_classic()

ggplot(df_2_test,
       aes(x = frag_size,
           y = counts,
           group = Sample,
           colour = Sample)) +
  geom_line(alpha = 1) +
  # geom_smooth(formula = log(y) ~ x,
  #             method = "loess",
  #             se = F) +
  xlab("Insert size") +
  ylab("Reads per sample") +
  xlim(0, 800) +
  coord_cartesian(expand = F) +
  scale_colour_manual(values = brewer.pal(n = 8,
                                          name = "Dark2"),
                      guide = "none") +
  # guides()
  # ylim(8, 15) +
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black")) +
  ggtitle("Microglia")


## DN ####

# select 8 files from the cell type
cell_type_to_plot <- "DN"

file_list_to_readin <-
  frag_file[str_detect(string = frag_file,
                       pattern = cell_type_to_plot)]
file_list_to_readin <-
  file_list_to_readin[str_detect(string = file_list_to_readin,
                                 pattern = "bwa_barcode_WASPed",
                                 negate = F)]
file_list_to_readin <-
  file_list_to_readin[order(file_list_to_readin)]
file_list_to_readin
file_list_to_readin <-
  file_list_to_readin[str_detect(string = file_list_to_readin,
                                 pattern = "10x",
                                 negate = T)]
# file_list_to_readin <-
#   file_list_to_readin[str_detect(string = file_list_to_readin,
#                                  pattern = "MG_",
#                                  negate = T)]
# file_list_to_readin <-
#   file_list_to_readin[1:8]
file_list_to_readin <-
  file_list_to_readin[sample(x = (length(file_list_to_readin) - 3),
                             size = 8)]

frag_list_2_bind <-
  vector(mode = "list",
         length = length(file_list_to_readin))
frag_list_2_bind <-
  lapply(X = 1:length(file_list_to_readin),
         FUN = function(x) {
           frag_df_raw <-
             read_delim(file_list_to_readin[x],
                        delim = "\t", escape_double = FALSE,
                        col_names = FALSE, trim_ws = TRUE)
           frag_df_raw$Sample <-
             file_list_to_readin[x]
           return(frag_df_raw)
         })

df_2_plot <-
  do.call(what = "rbind",
          args = frag_list_2_bind)
colnames(df_2_plot)[1:5] <-
  c("frag_size", "counts", "pos_strand", "neg_strand", "outward")
df_2_plot$Sample <-
  str_split(string = df_2_plot$Sample,
            pattern = "\\/",
            simplify = T)[, 2]
df_2_plot$Sample <-
  str_split(string = df_2_plot$Sample,
            pattern = "\\_",
            simplify = T)[, 1]

df_2_test <-
  df_2_plot %>%
  group_by(Sample) %>%
  summarise(counts = counts,
            frac_counts = counts / sum(counts),
            frag_size = frag_size)
# df_2_test$frag_size <-
#   factor(df_2_test$frag_size,
#          levels = sort(unique(df_2_test$frag_size)))
df_2_test$Sample <-
  factor(df_2_test$Sample,
         levels = unique(df_2_test$Sample))

ggplot(df_2_test,
       aes(x = frag_size,
           y = counts,
           group = Sample,
           colour = Sample)) +
  geom_line(alpha = 1) +
  # geom_smooth(formula = log(y) ~ x,
  #             method = "loess",
  #             se = F) +
  xlab("Insert size") +
  ylab("Reads per sample") +
  xlim(0, 800) +
  coord_cartesian(expand = F) +
  scale_colour_manual(values = brewer.pal(n = 8,
                                          name = "Dark2"),
                      guide = "none") +
  # guides()
  # ylim(8, 15) +
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black")) +
  ggtitle("Dopaminergic")

## GABA ####

# select 8 files from the cell type
cell_type_to_plot <- "GA"

file_list_to_readin <-
  frag_file[str_detect(string = frag_file,
                       pattern = cell_type_to_plot)]
file_list_to_readin <-
  file_list_to_readin[str_detect(string = file_list_to_readin,
                                 pattern = "barcode_WASPed",
                                 negate = F)]
file_list_to_readin <-
  file_list_to_readin[order(file_list_to_readin)]
file_list_to_readin <-
  file_list_to_readin[str_detect(string = file_list_to_readin,
                                 pattern = "10x",
                                 negate = T)]
# file_list_to_readin <-
#   file_list_to_readin[str_detect(string = file_list_to_readin,
#                                  pattern = "bowtie2",
#                                  negate = T)]
file_list_to_readin
# file_list_to_readin <-
#   file_list_to_readin[str_detect(string = file_list_to_readin,
#                                  pattern = "MG_",
#                                  negate = T)]
file_list_to_readin <-
  file_list_to_readin[1:8]
# file_list_to_readin <-
#   file_list_to_readin[sample(x = (length(file_list_to_readin) - 3),
#                              size = 8)]

frag_list_2_bind <-
  vector(mode = "list",
         length = length(file_list_to_readin))
frag_list_2_bind <-
  lapply(X = 1:length(file_list_to_readin),
         FUN = function(x) {
           frag_df_raw <-
             read_delim(file_list_to_readin[x],
                        delim = "\t", escape_double = FALSE,
                        col_names = FALSE, trim_ws = TRUE)
           frag_df_raw$Sample <-
             file_list_to_readin[x]
           return(frag_df_raw)
         })

df_2_plot <-
  do.call(what = "rbind",
          args = frag_list_2_bind)
colnames(df_2_plot)[1:5] <-
  c("frag_size", "counts", "pos_strand", "neg_strand", "outward")
df_2_plot$Sample <-
  str_split(string = df_2_plot$Sample,
            pattern = "\\/",
            simplify = T)[, 2]
df_2_plot$Sample <-
  str_split(string = df_2_plot$Sample,
            pattern = "\\_",
            simplify = T)[, 1]

df_2_test <-
  df_2_plot %>%
  group_by(Sample) %>%
  summarise(counts = counts,
            frac_counts = counts / sum(counts),
            frag_size = frag_size)
# df_2_test$frag_size <-
#   factor(df_2_test$frag_size,
#          levels = sort(unique(df_2_test$frag_size)))
df_2_test$Sample <-
  factor(df_2_test$Sample,
         levels = unique(df_2_test$Sample))

ggplot(df_2_test,
       aes(x = frag_size,
           y = counts,
           group = Sample,
           colour = Sample)) +
  geom_line(alpha = 1) +
  # geom_smooth(formula = log(y) ~ x,
  #             method = "loess",
  #             se = F) +
  xlab("Insert size") +
  ylab("Reads per sample") +
  xlim(0, 800) +
  coord_cartesian(expand = F) +
  scale_colour_manual(values = brewer.pal(n = 8,
                                          name = "Dark2"),
                      guide = "none") +
  # guides()
  # ylim(8, 15) +
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black")) +
  ggtitle("GABAergic")

## NGN2 ####
# select 8 files from the cell type
cell_type_to_plot <- c("Glut")

file_list_to_readin <-
  frag_file[str_detect(string = frag_file,
                       pattern = "Glut|R21")]
file_list_to_readin <-
  file_list_to_readin[str_detect(string = file_list_to_readin,
                                 pattern = "_barcode_WASPed",
                                 negate = F)]
file_list_to_readin <-
  file_list_to_readin[order(file_list_to_readin)]
# file_list_to_readin <-
#   file_list_to_readin[str_detect(string = file_list_to_readin,
#                                  pattern = "10x",
#                                  negate = T)]
# file_list_to_readin <-
#   file_list_to_readin[str_detect(string = file_list_to_readin,
#                                  pattern = "bowtie2",
#                                  negate = T)]
file_list_to_readin
# file_list_to_readin <-
#   file_list_to_readin[str_detect(string = file_list_to_readin,
#                                  pattern = "MG_",
#                                  negate = T)]
file_list_to_readin <-
  file_list_to_readin[c(1:7, 9)]
# file_list_to_readin <-
#   file_list_to_readin[sample(x = (length(file_list_to_readin) - 3),
#                              size = 8)]

frag_list_2_bind <-
  vector(mode = "list",
         length = length(file_list_to_readin))
frag_list_2_bind <-
  lapply(X = 1:length(file_list_to_readin),
         FUN = function(x) {
           frag_df_raw <-
             read_delim(file_list_to_readin[x],
                        delim = "\t", escape_double = FALSE,
                        col_names = FALSE, trim_ws = TRUE)
           frag_df_raw$Sample <-
             file_list_to_readin[x]
           return(frag_df_raw)
         })

df_2_plot <-
  do.call(what = "rbind",
          args = frag_list_2_bind)
colnames(df_2_plot)[1:5] <-
  c("frag_size", "counts", "pos_strand", "neg_strand", "outward")
df_2_plot$Sample <-
  str_split(string = df_2_plot$Sample,
            pattern = "\\/",
            simplify = T)[, 2]
df_2_plot$Sample <-
  str_split(string = df_2_plot$Sample,
            pattern = "\\_S",
            simplify = T)[, 1]
df_2_plot$Sample <-
  str_replace_all(string = df_2_plot$Sample,
            pattern = "Glut_rapid_neuron20_",
            replacement = "NGN2-")
df_2_test <-
  df_2_plot %>%
  group_by(Sample) %>%
  summarise(counts = counts,
            frac_counts = counts / sum(counts),
            frag_size = frag_size)
# df_2_test$frag_size <-
#   factor(df_2_test$frag_size,
#          levels = sort(unique(df_2_test$frag_size)))
df_2_test$Sample <-
  factor(df_2_test$Sample,
         levels = unique(df_2_test$Sample))

ggplot(df_2_test,
       aes(x = frag_size,
           y = counts,
           group = Sample,
           colour = Sample)) +
  geom_line(alpha = 1) +
  # geom_smooth(formula = log(y) ~ x,
  #             method = "loess",
  #             se = F) +
  xlab("Insert size") +
  ylab("Reads per sample") +
  xlim(0, 800) +
  coord_cartesian(expand = F) +
  scale_colour_manual(values = brewer.pal(n = 8,
                                          name = "Dark2"),
                      guide = "none") +
  # guides()
  # ylim(8, 15) +
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black")) +
  ggtitle("NGN2-Glut")

## Ast ####
# select 8 files from the cell type
cell_type_to_plot <- "Ast"

file_list_to_readin <-
  frag_file[str_detect(string = frag_file,
                       pattern = cell_type_to_plot)]
# file_list_to_readin <-
#   file_list_to_readin[str_detect(string = file_list_to_readin,
#                                  pattern = "bwa_barcode_WASPed",
#                                  negate = T)]
file_list_to_readin <-
  file_list_to_readin[order(file_list_to_readin)]
# file_list_to_readin <-
#   file_list_to_readin[str_detect(string = file_list_to_readin,
#                                  pattern = "10x",
#                                  negate = T)]
# file_list_to_readin <-
#   file_list_to_readin[str_detect(string = file_list_to_readin,
#                                  pattern = "bowtie2",
#                                  negate = T)]
file_list_to_readin
# file_list_to_readin <-
#   file_list_to_readin[str_detect(string = file_list_to_readin,
#                                  pattern = "MG_",
#                                  negate = T)]
file_list_to_readin <-
  file_list_to_readin[c(1:7, 9)]
# file_list_to_readin <-
#   file_list_to_readin[sample(x = (length(file_list_to_readin) - 3),
#                              size = 8)]

frag_list_2_bind <-
  vector(mode = "list",
         length = length(file_list_to_readin))
frag_list_2_bind <-
  lapply(X = 1:length(file_list_to_readin),
         FUN = function(x) {
           frag_df_raw <-
             read_delim(file_list_to_readin[x],
                        delim = "\t", escape_double = FALSE,
                        col_names = FALSE, trim_ws = TRUE)
           frag_df_raw$Sample <-
             file_list_to_readin[x]
           return(frag_df_raw)
         })

df_2_plot <-
  do.call(what = "rbind",
          args = frag_list_2_bind)
colnames(df_2_plot)[1:5] <-
  c("frag_size", "counts", "pos_strand", "neg_strand", "outward")
df_2_plot$Sample <-
  str_split(string = df_2_plot$Sample,
            pattern = "\\/",
            simplify = T)[, 2]
df_2_plot$Sample <-
  str_split(string = df_2_plot$Sample,
            pattern = "\\_S",
            simplify = T)[, 1]
# df_2_plot$Sample <-
#   str_replace_all(string = df_2_plot$Sample,
#                   pattern = "Glut_rapid_neuron20_",
#                   replacement = "NGN2-")
df_2_test <-
  df_2_plot %>%
  group_by(Sample) %>%
  summarise(counts = counts,
            frac_counts = counts / sum(counts),
            frag_size = frag_size)
# df_2_test$frag_size <-
#   factor(df_2_test$frag_size,
#          levels = sort(unique(df_2_test$frag_size)))
df_2_test$Sample <-
  factor(df_2_test$Sample,
         levels = unique(df_2_test$Sample))

ggplot(df_2_test,
       aes(x = frag_size,
           y = counts,
           group = Sample,
           colour = Sample)) +
  geom_line(alpha = 1) +
  # geom_smooth(formula = log(y) ~ x,
  #             method = "loess",
  #             se = F) +
  xlab("Insert size") +
  ylab("Reads per sample") +
  xlim(0, 800) +
  coord_cartesian(expand = F) +
  scale_colour_manual(values = brewer.pal(n = 8,
                                          name = "Dark2"),
                      guide = "none") +
  # guides()
  # ylim(8, 15) +
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black")) +
  ggtitle("Astrocytes")
# ggplot(df_2_test,
#        aes(x = frag_size,
#            y = log(frac_counts),
#            group = Sample,
#            colour = Sample)) +
#   geom_line() +
#   # geom_smooth(formula = log(y) ~ x,
#   #             method = "loess",
#   #             se = F) +
#   scale_colour_manual(values = brewer.pal(n = 8,
#                                           name = "Dark2")) +
#   ylim(-8, -4) +
#   theme_classic()

ggplot(df_2_test,
       aes(x = frag_size,
           y = counts,
           group = Sample,
           colour = Sample)) +
  geom_line(alpha = 1) +
  # geom_smooth(formula = log(y) ~ x,
  #             method = "loess",
  #             se = F) +
  xlab("Insert size") +
  # scale_colour_manual(values = brewer.pal(n = 8,
  #                                         name = "Dark2")) +
  # ylim(8, 15) +
  theme_classic() +
  ggtitle("Astrocytes, bowtie2-aligned")
