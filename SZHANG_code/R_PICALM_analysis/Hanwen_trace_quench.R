# Siwei 13 Sept 2023
# make trace plot for Hanwen
# compensate for background quenching effect

# init ####
library(readxl)
library(baseline)
library(gWidgets2)
library(tidyr)

library(ggplot2)
library(RColorBrewer)

# load data ####
Stimulation_trace_signal <-
  read_excel("Stimulation_traces-7hr-013123.xlsx",
             sheet = "Cells", col_names = FALSE)
Stimulation_trace_bkgrd <-
  read_excel("Stimulation_traces-7hr-013123.xlsx",
             sheet = "Background", col_names = FALSE)

Stimulation_trace_signal <-
  Stimulation_trace_signal[, -1]
# Stimulation_trace_bkgrd <-
#   Stimulation_trace_bkgrd[, -1]

# Stimulation_trace_signal_t <-
#   as.matrix(t(Stimulation_trace_signal))
# Stimulation_trace_signal_t <-
#   baseline(spectra = Stimulation_trace_signal_t,
#            method = "irls")

# corrected_stim <-
#   getCorrected(Stimulation_trace_signal_t)
# plot(t(corrected_stim)[, 1],
#      xlim = c(1, 15))

####

Stimulation_trace_signal_plateau <-
  as.matrix(Stimulation_trace_signal[1:14, ])

add_timepoint <- vector(mode = "list",
                        length = 9L)

for (i in 1:length(add_timepoint)) {
  add_timepoint[[i]] <- Stimulation_trace_signal_plateau
}


Stim_output <- Stimulation_trace_signal_plateau
for (i in 1:length(add_timepoint)) {
  Stim_output <-
    rbind(Stim_output,
          Stimulation_trace_signal_plateau)
}

random_stim_output <-
  data.frame(x = rnorm(n = nrow(Stim_output),
                       mean = colMeans(Stim_output)[1],
                       sd = 10),
             y = rnorm(n = nrow(Stim_output),
                       mean = colMeans(Stim_output)[2],
                       sd = 10),
             z = rnorm(n = nrow(Stim_output),
                       mean = colMeans(Stim_output)[3],
                       sd = 10))
colnames(random_stim_output) <-
  colnames(Stimulation_trace_signal)

Stimulation_trace_signal_all <-
  rbind(random_stim_output,
        Stimulation_trace_signal)
# Stimulation_trace_signal_all <-
#   Stimulation_trace_signal_all[, -1]

Stimulation_trace_signal_t <-
  as.matrix(t(Stimulation_trace_signal_all))
Stimulation_trace_signal_t <-
  baseline(spectra = Stimulation_trace_signal_t,
           method = "irls")

corrected_stim <-
  getCorrected(Stimulation_trace_signal_t)
corrected_stim <-
  as.data.frame(t(corrected_stim))
corrected_stim <-
  corrected_stim - min(corrected_stim, na.rm = T)

plot_stim <-
  as.data.frame(cbind(1:nrow(corrected_stim),
                      corrected_stim))
colnames(plot_stim) <-
  c("x", "cell_1", "cell_2", "cell_3")

plot_stim_long <-
  gather(plot_stim,
         key = "intensity",
         value = "value",
         c("cell_1", "cell_2", "cell_3"))

ggplot(plot_stim_long,
       aes(x = x,
           y = value,
           color = factor(intensity))) +
  geom_smooth(stat = "identity") +
  facet_wrap(~ intensity,
             ncol = 1) +
  xlim(1, 1000) +
  xlab("frames") +
  ylab("Baseline-corrected intensity level") +
  theme_linedraw() +
  scale_colour_brewer(palette = "Dark2",
                      name = "intensity")

save.image(file = "Hanwen_plot_baseline_3_cells.RData")


corrected_stim_single <-
  rowMeans(corrected_stim,
           na.rm = T)
plot_stim <-
  as.data.frame(cbind(1:length(corrected_stim_single),
                      corrected_stim_single))
colnames(plot_stim) <-
  c("x", "intensity")

plot_stim$x <-
  plot_stim$x - 140

ggplot(plot_stim,
       aes(x = x / 360,
           y = intensity)) +
  geom_smooth(stat = "identity",
              colour = "darkblue") +
  xlim(1, 1000) +
  xlab("frames") +
  ylab("Baseline-corrected intensity level") +
  theme_linedraw()

ggplot(plot_stim,
       aes(x = x / 360,
           y = intensity)) +
  geom_smooth(stat = "identity",
              colour = "darkred") +
  scale_x_continuous(breaks = c(0:7)) +
  # xlim(1, 1000) +
  xlab("hours post KCl stimulation") +
  ylab("Baseline-corrected \nintensity level") +
  theme_linedraw() +
  theme(axis.text = element_text(size = 12))

save.image(file = "Hanwen_plot_baseline_3_cells.RData")
