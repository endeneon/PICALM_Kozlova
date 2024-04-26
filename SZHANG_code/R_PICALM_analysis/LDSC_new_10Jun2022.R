# Siwei 08 Jun 2022
# plot h2 and enrichment for LDSC results

# init

library(readxl)
library(readr)

library(ggplot2)
library(RColorBrewer)


# load data
combined_output <- 
  read_delim("combined_output.txt", 
             delim = "\t", escape_double = FALSE, 
             col_names = FALSE, trim_ws = TRUE)

df_to_plot <- combined_output


colnames(df_to_plot) <- c("Enrichment", "Enrichment_pVal")


df_to_plot$cell_type <-
  rep_len(x = c("Ast", "MGlia", "GABA", "NGN2-Glut"),
          length.out = nrow(df_to_plot))
df_to_plot$cell_type <- 
  factor(df_to_plot$cell_type,
         levels = c("Ast", "MGlia", "GABA", "NGN2-Glut"))

df_to_plot$Enrichment <- as.numeric(df_to_plot$Enrichment)
df_to_plot$Enrichment_pVal <- as.numeric(df_to_plot$Enrichment_pVal)

df_to_plot$Diseases <-
  rep(c("ADHD", "Alz", "Bipolar", "SCZw3", "IBD", 
        "MDD", "ASD", "T2D", "UC"),
      each = 4)
df_to_plot$Diseases <-
  factor(df_to_plot$Diseases,
         levels = c("ADHD", "Alz", "ASD", "Bipolar", "MDD", "SCZw3", 
                    "IBD", "T2D", "UC"))

# df_to_plot_backup <- df_to_plot
df_to_plot$Enrichment[df_to_plot$Enrichment < 0] <- NA

ggplot(df_to_plot, 
       aes(y = cell_type,
           x = Diseases,
           #           shape = ifelse(h2 > 0, 21, 23),
           size = Enrichment,
           fill = (0 - log10(Enrichment_pVal)))) +
  geom_point(shape = 21) +
  # scale_y_discrete(limits = rev(levels(df_to_plot$cell_type))) +
  scale_fill_gradient(low = "white", high = "darkred") +
  # scale_alpha_continuous(range = c(0.2, 1)) +
  labs(size = "Enrichment",
       fill = "Enrichment -log10P value") +
  # guides(fill = "none") +
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
  scale_fill_manual(values = c("white", "darkred")) +
  scale_alpha_continuous(range = c(0.2, 1)) +
  labs(size = "Enrichment",
       alpha = "h2") +
  guides(fill = "none") +
  xlab("Diseases") +
  ylab("Peak origin") +
  theme_bw()
