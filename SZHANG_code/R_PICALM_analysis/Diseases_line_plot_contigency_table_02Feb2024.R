# Siwei 02 Feb 2024
# Plot Jubao's Diseases bar plot

# init ####
library(readxl)

library(reshape2)
library(plyr)

library(ggplot2)
library(scales)

library(RColorBrewer)

# load data #####
df_raw <-
  read_excel("Jubao_Diseases_bar_plot_v2.xlsx")

df_list_by_cell_type <-
  split(x = df_raw,
        f = df_raw$Cell_Type)

df_list_pval_by_cell_type <-
  vector(mode = "list",
         length = length(df_list_by_cell_type))

## write a function to return the plotting values #####
return_pval <-
  function(df_each_type) {
    # df_each_type <-
    #   df_list_by_cell_type[[1]]
    ### record cell type
    cell_type_record <-
      unlist(unique(df_each_type$Cell_Type))
    dynamic_pattern_record <-
      df_each_type$Dynamic_Pattern
    df_each_type <-
      df_each_type[, 2:9]

    df_results_2_return <-
      as.data.frame(matrix(0,
                           nrow = (nrow(df_each_type) - 1),
                           ncol = (ncol(df_each_type) - 1)))
    colnames(df_results_2_return) <-
      colnames(df_each_type)[1:(length(colnames(df_each_type)) - 1)]
    rownames(df_results_2_return) <-
      dynamic_pattern_record[1:(length(dynamic_pattern_record) - 1)]

    ## fill column-wise first
    for (row_number in 1:nrow(df_results_2_return)) {
      for (col_number in 1:ncol(df_results_2_return)) {
        # print(paste(row_number,
        #             col_number))
        matrix_2_test <-
          matrix(c(unlist(df_each_type[row_number, col_number]),
                   unlist(df_each_type[(nrow(df_results_2_return) + 1), col_number]),
                   unlist(df_each_type[row_number, (ncol(df_results_2_return) + 1)]),
                   unlist(df_each_type[(nrow(df_results_2_return) + 1), (ncol(df_results_2_return) + 1)])),
                 nrow = 2,
                 byrow = T)
        df_results_2_return[row_number, col_number] <-
          0 - log2(fisher.test(x = matrix_2_test)$p.value)
        if ((row_number == 1) & (col_number == 1)) {
          print(matrix_2_test)
          print(paste("-log2P = ",
                      0 - log2(fisher.test(x = matrix_2_test)$p.value)))
        }
      }
    }
    df_results_2_return$Dynamic_Pattern <-
      rownames(df_results_2_return)
    df_results_2_return$Cell_Type <-
      cell_type_record
  return(df_results_2_return)
  }

for (m in 1:length(df_list_by_cell_type)) {
  print(m)
  df_list_pval_by_cell_type[[m]] <-
    return_pval(df_list_by_cell_type[[m]])
}

df_cum_2_plot <-
  ldply(df_list_pval_by_cell_type,
        data.frame)
colnames(df_cum_2_plot)

df_cum_2_plot_melt <-
  reshape2::melt(data = df_cum_2_plot,
                 id = c(8, 9))
colnames(df_cum_2_plot_melt)[3] <- "Diseases"
colnames(df_cum_2_plot_melt)[4] <- '-log2Pval'

df_cum_2_plot_melt$Dynamic_Pattern <-
  factor(df_cum_2_plot_melt$Dynamic_Pattern,
         levels = c("up-up",
                    "up-flat",
                    "up-down",
                    "flat-up",
                    "down-up",
                    "flat-down",
                    "down-down",
                    "down-flat"))

df_cum_2_plot_melt$Diseases <-
  factor(df_cum_2_plot_melt$Diseases,
         levels = unique(df_cum_2_plot_melt$Diseases))

ggplot(df_cum_2_plot_melt,
       aes(x = factor(Dynamic_Pattern),
           y = `-log2Pval`,
           # group = Diseases,
           colour = factor(Diseases))) +
  # geom_line() +
  geom_line(aes(group = Diseases)) +
  scale_colour_manual(values = brewer.pal(n = 8,
                                          name = "Dark2"),
                      name = "Diseases") +
  facet_wrap(~ Cell_Type,
             nrow = 1) +
  labs(y = "-log2 P value",
       x = "Dynamic Patterns") +
  # ylim(c(0, 50)) +
  # guides(fill = "Dynamic patterns") +
  # scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0))

0-log2(fisher.test(matrix(c(6,1,1683,1376), nrow = 2, byrow = T))$p.value)
