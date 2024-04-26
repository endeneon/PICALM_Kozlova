# Siwei 16 Oct 2023

# init ####
library(ggplot2)
library(readr)
library(stringr)

df_raw <-
  read_table("SCZw3_summary_2_transform.txt",
             col_names = FALSE)

df_2_plot <-
  data.frame(annotation = str_split(string = df_raw$X1,
                                    pattern = '\\.',
                                    simplify = T)[, 1])

## from Min Qiao's python code
df_2_plot$log2_fold_enrichment <-
  round(log2(exp(df_raw$X2)),
        digits = 4)

df_2_plot$CI_low <-
  round(log2(exp(df_raw$X3)),
        digits = 4)
df_2_plot$CI_hi <-
  round(log2(exp(df_raw$X4)),
        digits = 4)

ggplot(data = df_2_plot,
       aes(y = log2_fold_enrichment,
           x = annotation,
           ymin = CI_low,
           ymax = CI_hi)) +
  geom_point(colour = "darkred") +
  geom_linerange() +
  scale_x_discrete(limit = rev) +
  coord_flip() +
  theme_classic()
