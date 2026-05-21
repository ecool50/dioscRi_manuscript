# N_ref sensitivity figure (response-only).
# Two panels: VAE reconstruction error (MSE) and pipeline runtime vs K.
# AUC vs K is in manuscript Supp Fig S5F and is not duplicated here.
#
# Input:  results/nref/nref_sensitivity_results.csv
# Output: figures/supp_fig_s7_nref_sensitivity.pdf

suppressPackageStartupMessages({
  library(here)
  library(ggplot2)
  library(patchwork)
})

revision_root <- here::here("revision")
nref <- read.csv(file.path(revision_root, "results/nref/nref_sensitivity_results.csv"))
nref$time_minutes <- nref$time_seconds / 60
nref_levels <- sort(unique(nref$n_ref))

p_mse <- ggplot(nref, aes(x = n_ref, y = recon_mse)) +
  geom_line(colour = "#E69F00", linewidth = 1) +
  geom_point(colour = "#E69F00", size = 3) +
  scale_x_continuous(breaks = nref_levels) +
  labs(x = "Number of reference samples", y = "MSE",
       title = "VAE reconstruction error") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

p_time <- ggplot(nref, aes(x = n_ref, y = time_minutes)) +
  geom_line(colour = "#4A90D9", linewidth = 1) +
  geom_point(colour = "#4A90D9", size = 3) +
  scale_x_continuous(breaks = nref_levels) +
  labs(x = "Number of reference samples", y = "Time (minutes)",
       title = "Pipeline runtime") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

p_s7 <- p_mse | p_time

out <- file.path(revision_root, "figures/supp_fig_s7_nref_sensitivity.pdf")
ggsave(out, p_s7, width = 10, height = 4)
cat("Saved:", out, "\n")
