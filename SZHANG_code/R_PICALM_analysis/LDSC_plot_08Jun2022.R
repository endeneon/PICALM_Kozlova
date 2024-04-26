# Siwei 08 Jun 2022
# plot h2 and enrichment for LDSC results

# init

library(readxl)
library(ggplot2)
library(RColorBrewer)


# load data
LDSC_MG_AST_h2_enrichment <- 
  read_excel("LDSC_MG_AST_h2_enrichment_w_Alz.xlsx", 
             col_names = FALSE)


df_to_plot <- as.data.frame(t(LDSC_MG_AST_h2_enrichment))
colnames(df_to_plot) <- unlist(df_to_plot[1, ])
df_to_plot <- df_to_plot[-1, ]

df_to_plot$cell_type <-
  rep_len(x = c("Ast", "MGlia", "NGN2-Glut"),
          length.out = 21)
df_to_plot$cell_type <- 
  factor(df_to_plot$cell_type,
         levels = c("Ast", "MGlia", "NGN2-Glut"))


df_to_plot$h2 <- as.numeric(df_to_plot$h2)
df_to_plot$enrichment <- as.numeric(df_to_plot$enrichment)



ggplot(df_to_plot, 
       aes(y = cell_type,
           x = Diseases,
#           shape = ifelse(h2 > 0, 21, 23),
           size = abs(enrichment),
           fill = ifelse(enrichment > 0, "B", "R"),
           alpha = abs(h2))) +
  geom_point(shape = ifelse(df_to_plot$h2 > 0, 21, 23)) +
  scale_y_discrete(limits = rev(levels(df_to_plot$Diseases))) +
  scale_fill_manual(values = c("darkred", "darkblue")) +
  scale_alpha_continuous(range = c(0.2, 1)) +
  labs(size = "Enrichment",
       alpha = "h2") +
  guides(fill = "none") +
  xlab("Diseases") +
  ylab("Peak origin") +
  theme_bw()

# remove all h2 < 0
df_to_plot_v2 <- df_to_plot

df_to_plot_v2$enrichment[df_to_plot_v2$h2 < 0] <- NA

ggplot(df_to_plot_v2, 
       aes(y = cell_type,
           x = Diseases,
           #           shape = ifelse(h2 > 0, 21, 23),
           size = abs(enrichment),
           fill = ifelse(enrichment > 0, "B", "R"),
           alpha = abs(h2))) +
  geom_point(shape = ifelse(df_to_plot_v2$h2 > 0, 21, 23)) +
  scale_y_discrete(limits = rev(levels(df_to_plot$Diseases))) +
  scale_fill_manual(values = c("darkred", "darkblue")) +
  scale_alpha_continuous(range = c(0.2, 1)) +
  labs(size = "Enrichment",
       alpha = "h2") +
  guides(fill = "none") +
  xlab("Diseases") +
  ylab("Peak origin") +
  theme_bw()
