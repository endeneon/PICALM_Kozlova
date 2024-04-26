# Siwei 22 Feb 2024
# plot Jubao's lysosome genes

# init ####
library(ggplot2)
library(RColorBrewer)

library(scales)

library(readxl)

df_raw <-
  read_excel(path = "plot_lysosome_genes.xlsx")

df_2_plot <-
  df_raw
df_2_plot$`-log10P` <-
  0 - log10(df_2_plot$`p-value`)

df_2_plot$Gene <-
  factor(df_2_plot$Gene,
         levels = df_2_plot$Gene[1:27])

unique(df_2_plot$Category)
df_2_plot$Category <-
  factor(df_2_plot$Category,
         levels = unique(df_2_plot$Category)[c(1,2,4,3)])

ggplot(df_2_plot,
       aes(x = Category,
           y = Gene,
           size = `-log10P`,
           fill = log2FC)) +
  geom_point(shape = 21) +
  scale_fill_gradient2(low = "darkblue",
                       mid = "white",
                       high = "darkred",
                       midpoint = 0) +
  # legend(position = "below") +
  scale_x_discrete(limits = rev) +
  ylab("") +
  theme_classic() +
  coord_flip() +
  theme(legend.position = "bottom",
        plot.margin = margin(0.2, 0.5, 0, 0.2,
                               unit = "in"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0),
        axis.text = element_text(colour = "black"))
