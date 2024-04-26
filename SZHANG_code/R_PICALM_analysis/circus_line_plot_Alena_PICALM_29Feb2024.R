# Siwei 29 Feb 2024
# Make line-dot and circus plots for Alena's PICALM paper

# init ####
library(readxl)
library(ggplot2)
library(RColorBrewer)

## make line-dot plot #####
df_raw <-
  read_excel("IPA_pathway_for_graphing.xlsx",
             sheet = 2)

df_2_plot <-
  df_raw
df_2_plot$Ingenuity_Canonical_Pathways <-
  factor(df_2_plot$Ingenuity_Canonical_Pathways,
         levels = unique(df_2_plot$Ingenuity_Canonical_Pathways))

ggplot(df_2_plot,
       aes(x = Ingenuity_Canonical_Pathways,
           y = `-log_pVal`,
           fill = zScore)) +
  geom_col() +
  scale_fill_gradient2(low = "darkblue",
                       mid = "white",
                       high = "darkred",
                       limits = c(0 - max(df_2_plot$zScore),
                                  max(df_2_plot$zScore))) +
  xlab("Ingenuity Canonical Pathways") +
  ylab("Z-score") +
  scale_x_discrete(limits = rev) +
  guides(size = guide_legend('-logP value')) +
  coord_flip() +
  theme_classic() +
  theme(axis.text = element_text(colour = "black"),
        axis.text.x = element_text(size = 12))



ggplot(df_2_plot,
       aes(x = Ingenuity_Canonical_Pathways,
           y = zScore,
           ymax = zScore,
           size = `-log_pVal`)) +
  geom_linerange(ymin = 0,
                 size = 1) +
  geom_point(colour = "darkred") +
  scale_radius(limits = c(3, 7)) +
  scale_x_discrete(limits = rev) +
  xlab("Ingenuity Canonical Pathways") +
  ylab("Z-score") +
  guides(size = guide_legend('-logP value')) +
  coord_flip() +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12))
