# 14 Jan 2024
# make Fig 3A/B

# init ####
library(readxl)
library(ggplot2)
library(RColorBrewer)

library(reshape2)

# load data ####
df_raw_Fig3a <-
  read_excel("Fig_3a.xlsx")

df_to_plot <-
  df_raw_Fig3a
colnames(df_to_plot) <-
  c("Gene", "CW20107", "KOLF2.2J", "CD-14")

df_to_plot <-
  melt(df_to_plot)
colnames(df_to_plot) <-
  c("Gene",
    "Cell_Line",
    "Post_Sorting_Viability")

df_to_plot$Gene <-
  factor(df_to_plot$Gene,
         levels = c("CUL1", "SETD1A", "HERC1", "Negative Control",
                    "TRIO", "AKAP11", "ASH1L", "XPO7",
                    "CACNA1G",
                    "KDM6B", "CHD8", "SMARCC2",
                    "GRIA3", "GABRA1",
                    "SCN2A", "DLL1", "HCN4", "GRIN2A",
                    "ARID1B", "ANKRD11",
                    "RB1CC1", "SP4", "KMT2C", "SHANK3"))
df_to_plot$Post_Sorting_Viability <-
  df_to_plot$Post_Sorting_Viability / 100

ggplot(df_to_plot,
       aes(x = Gene,
           y = Post_Sorting_Viability,
           fill = Cell_Line)) +
  geom_bar(position = "dodge",
           stat = "identity") +
  scale_fill_manual(values = brewer.pal(n = 3,
                                        name = "Dark2")) +
  xlab("Gene") +
  ylab("Post Sorting Viability") +
  scale_y_continuous(labels = scales::percent) +
  labs(fill = "Cell Line") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0))


# Fig.3b ####
df_raw_Fig3b <-
  read_excel("Fig_3a.xlsx",
             sheet = 2)

df_to_plot <-
  df_raw_Fig3b
colnames(df_to_plot) <-
  c("Gene", "CW20107", "KOLF2.2J", "CD-14")

df_to_plot <-
  melt(df_to_plot)
colnames(df_to_plot) <-
  c("Gene",
    "Cell_Line",
    "Editing_Efficiency")

# df_to_plot$Gene <-
#   factor(df_to_plot$Gene,
#          levels = c("CUL1", "SETD1A", "HERC1",
#                     # "Negative Control",
#                     "TRIO", "AKAP11", "ASH1L", "XPO7",
#                     "CACNA1G",
#                     "KDM6B", "CHD8", "SMARCC2",
#                     "GRIA3", "GABRA1",
#                     "SCN2A", "DLL1", "HCN4", "GRIN2A",
#                     "ARID1B", "ANKRD11",
#                     "RB1CC1", "SP4", "KMT2C", "SHANK3"))
df_to_plot$Gene <-
  factor(df_to_plot$Gene,
         levels = c("HERC1", "ANKRD11", "KMT2C", "GABRA1",
                    "AKAP11", "SHANK3", "SETD1A", "CUL1",
                    "TRIO", "GRIA3", "SP4", "CACNA1G",
                    "CHD8", "ASH1L", "GRIN2A", "ARID1B",
                    "KDM6B", "DLL1", "SCN2A", "SMARCC2",
                    "RB1CC1", "HCN4", "XPO7"))
# df_to_plot$Post_Sorting_Viability <-
#   df_to_plot$Post_Sorting_Viability / 100

ggplot(df_to_plot,
       aes(x = Gene,
           y = Editing_Efficiency,
           fill = Cell_Line)) +
  geom_bar(position = "dodge",
           stat = "identity") +
  scale_fill_manual(values = brewer.pal(n = 3,
                                        name = "Dark2")) +
  xlab("Gene") +
  ylab("Editing Efficiency") +
  scale_y_continuous(labels = scales::percent) +
  labs(fill = "Cell Line") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0))
