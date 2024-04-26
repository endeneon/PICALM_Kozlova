# make disease/gene group plots for LDSC enrichment
# peaks enriched from bulk ATAC-seq Ast, MG, GA, NGN2
# focus on the enrichment of Als and Alz
# Siwei 26 Jun 2023

# init
library(ggplot2)
library(readr)
library(RColorBrewer)
library(stringr)

# load data #####
raw_df <-
  read_delim("combined_output_Ast_MG_GA_NGN2_26Jun2023.txt",
                delim = "\t", escape_double = FALSE,
                trim_ws = TRUE)

raw_df <-
  read_delim("combined_output_Ast_MG_GA_NGN2_GA_27Jun2023.txt",
             delim = "\t", escape_double = FALSE,
             trim_ws = TRUE)

raw_df <-
  read_delim("combined_output_Ast_MG_GA_NGN2_GA_30Jun2023.txt",
             delim = "\t", escape_double = FALSE,
             trim_ws = TRUE)

raw_df <-
  read_delim("combined_output_26Feb2024.txt",
             delim = "\t", escape_double = FALSE,
             trim_ws = TRUE)

# assign corresponding column names #####
df_data <- raw_df
colnames(df_data)[1] <- "Disease"

df_data$Disease <-
  str_split(string = df_data$Disease,
            pattern = "\\.",
            simplify = T)[, 1]

df_data$Cell_type <-
  rep_len(c("Ast", "MG", "GA", "NGN2-Glut", "DN"),
          length.out = nrow(df_data))
df_data$Cell_type <-
  factor(df_data$Cell_type,
         levels = c("Ast", "MG", "GA", "NGN2-Glut", "DN"))
df_data$`-log10P` <-
  0 - log10(df_data$Enrichment_p)

# make plot in Xiaotong's style #####
df_to_plot <-
  df_data[df_data$Cell_type %in% "MG", ]

ggplot(df_to_plot,
       aes(x = Disease,
           y = (Enrichment),
           ymin = (Enrichment - Enrichment_std_error),
           ymax = (Enrichment + Enrichment_std_error),
           size = `-log10P`)) +
  geom_linerange(linewidth = 1) +
  geom_point(shape = 21,
              fill = "darkred",
              color = "transparent") +
  ylab("Enrichment (fold)") +
  scale_radius() +
  # ylim(0, 4) +
  scale_x_discrete(limits = rev) +
  coord_flip() +
  theme_bw() +
  ggtitle("MG peaks, q=0.05")

df_to_plot <-
  df_data[df_data$Cell_type %in% "Ast", ]

ggplot(df_to_plot,
       aes(x = Disease,
           y = (Enrichment),
           ymin = (Enrichment - Enrichment_std_error),
           ymax = (Enrichment + Enrichment_std_error),
           size = `-log10P`)) +
  geom_linerange(linewidth = 1) +
  geom_point(shape = 21,
             fill = "darkred",
             color = "transparent") +
  ylab("Enrichment (fold)") +
  scale_radius() +
  # ylim(0, 4) +
  scale_x_discrete(limits = rev) +
  coord_flip() +
  theme_bw() +
  ggtitle("Ast peaks, q=0.05")

df_to_plot <-
  df_data[df_data$Cell_type %in% "GA", ]

ggplot(df_to_plot,
       aes(x = Disease,
           y = (Enrichment),
           ymin = (Enrichment - Enrichment_std_error),
           ymax = (Enrichment + Enrichment_std_error),
           size = `-log10P`)) +
  geom_linerange(linewidth = 1) +
  geom_point(shape = 21,
             fill = "darkred",
             color = "transparent") +
  ylab("Enrichment (fold)") +
  scale_radius() +
  # ylim(0, 4) +
  scale_x_discrete(limits = rev) +
  coord_flip() +
  theme_bw() +
  ggtitle("GA peaks, q=0.05")

df_to_plot <-
  df_data[df_data$Cell_type %in% "DN", ]

ggplot(df_to_plot,
       aes(x = Disease,
           y = (Enrichment),
           ymin = (Enrichment - Enrichment_std_error),
           ymax = (Enrichment + Enrichment_std_error),
           size = `-log10P`)) +
  geom_linerange(linewidth = 1) +
  geom_point(shape = 21,
             fill = "darkred",
             color = "transparent") +
  ylab("Enrichment (fold)") +
  scale_radius() +
  # ylim(0, 4) +
  scale_x_discrete(limits = rev) +
  coord_flip() +
  theme_bw() +
  ggtitle("DN peaks, q=0.05")

df_to_plot <-
  df_data[df_data$Cell_type %in% "NGN2-Glut", ]

ggplot(df_to_plot,
       aes(x = Disease,
           y = (Enrichment),
           ymin = (Enrichment - Enrichment_std_error),
           ymax = (Enrichment + Enrichment_std_error),
           size = `-log10P`)) +
  geom_linerange(linewidth = 1) +
  geom_point(shape = 21,
             fill = "darkred",
             color = "transparent") +
  ylab("Enrichment (fold)") +
  scale_radius() +
  # ylim(0, 4) +
  scale_x_discrete(limits = rev) +
  coord_flip() +
  theme_bw() +
  ggtitle("NGN2-Glut peaks, q=0.05")

# make plot in conventional style #####


df_to_plot <-
  df_data
# df_to_plot$`-logP` <-
#   0 - log10(df_to_plot$Enrichment_p)
# make the plot #####

df_to_plot <-
  df_to_plot[str_detect(string = df_to_plot$Disease,
                        pattern = "Als_15*",
                        negate = T), ]
df_to_plot <-
  df_to_plot[str_detect(string = df_to_plot$Disease,
                        pattern = "bentham*",
                        negate = T), ]
df_to_plot$Disease <-
  str_replace_all(string = df_to_plot$Disease,
                  pattern = "Mc_",
                  replacement = "Mc-")
# df_to_plot$Disease <-
#   str_replace_all(string = df_to_plot$Disease,
#                   pattern = "Mc_",
#                   replacement = "Mc-")
df_to_plot$Disease <-
  str_split(string = df_to_plot$Disease,
            pattern = "_",
            simplify = T)[, 1]
df_to_plot$Disease <-
  str_replace_all(string = df_to_plot$Disease,
                  pattern = "daner-SCZ3",
                  replacement = "Schizophrenia")
df_to_plot$Disease <-
  str_replace_all(string = df_to_plot$Disease,
                  pattern = "jointGWASMc-",
                  replacement = "")
df_to_plot$Disease <-
  str_replace_all(string = df_to_plot$Disease,
                  pattern = "Total",
                  replacement = "Total-Cholestrol")
df_to_plot$Disease <-
  str_replace_all(string = df_to_plot$Disease,
                  pattern = "MDD2018",
                  replacement = "MDD")
df_to_plot$Disease <-
  str_replace_all(string = df_to_plot$Disease,
                  pattern = "MDD",
                  replacement = "Depression")
df_to_plot$Disease <-
  str_replace_all(string = df_to_plot$Disease,
                  pattern = "PGC-ASD-2019",
                  replacement = "Autism")
df_to_plot$Disease <-
  str_replace_all(string = df_to_plot$Disease,
                  pattern = "neuroticism",
                  replacement = "Neuroticism")
df_to_plot$Disease <-
  str_replace_all(string = df_to_plot$Disease,
                  pattern = "Alz",
                  replacement = "Alzheimer\\'s")
df_to_plot$Disease <-
  str_replace_all(string = df_to_plot$Disease,
                  pattern = "Parkinson",
                  replacement = "Parkinson\\'s")

unique(df_to_plot$Disease)

write.table(df_to_plot,
            file = "df_plot_sLDSC_w_h2.txt",
            row.names = F,
            quote = F, sep = "\t")


###
df_to_plot$Category <- NULL
write.table(df_to_plot,
            file = "df_plot_sLDSC.txt",
            row.names = F,
            quote = F, sep = "\t")


df_to_plot$Disease <-
  factor(df_to_plot$Disease,
         levels = c("Schizophrenia",
                    "Bipolar",
                    "Autism",
                    "Depression",
                    "ADHD",
                    "Neuroticism",
                    "Intelligence",
                    "PTSD",
                    "Alzheimer's",
                    "Parkinson's",
                    "HDL",
                    "LDL",
                    "Total-Cholestrol",
                    "T2D"))



ggplot(df_to_plot,
       aes(x = Disease,
           y = Cell_type,
           fill = Enrichment,
           size = `-log10P`)) +
  geom_point(shape = 21) +
  scale_fill_gradient(low = "white",
                      high = "darkred") +
  # scale_size_continuous() +
  # scale_radius() +
  xlab("Diseases") +
  labs(size = "-log10P") +
  scale_y_discrete(limits = rev) +
  ylab("Peak types") +
  # scale_alpha_continuous(df_to_plot$`-logP`,
  #                        range = c(0, 1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 315,
                                   vjust = 0,
                                   hjust = 0),
        axis.text = element_text(size = 12)) +
  ggtitle("s-LDSC enrichment")

df_to_plot <-
  df_to_plot[!(df_to_plot$Cell_type %in% "GA"), ]

ggplot(df_to_plot,
       aes(x = Disease,
           y = Cell_type,
           fill = Enrichment,
           size = `-log10P`)) +
  geom_point(shape = 21) +
  scale_fill_gradient(low = "white",
                      high = "darkred") +
  # scale_size_continuous() +
  # scale_radius() +
  xlab("Diseases") +
  labs(size = "-log10P") +
  scale_y_discrete(limits = rev) +
  ylab("Peak types") +
  # scale_alpha_continuous(df_to_plot$`-logP`,
  #                        range = c(0, 1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 315,
                                   vjust = 0,
                                   hjust = 0),
        axis.text = element_text(size = 12)) +
  ggtitle("s-LDSC enrichment")
