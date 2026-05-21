# Supp Fig S8: stratified tier comparison, BioHEART + CMV, MSE + AUC.
# Four panels: A=BH MSE, B=BH AUC, C=CMV MSE, D=CMV AUC.
#
# Inputs:  results/nref/nref_stratified_results.csv
#          results/nref/cmv_nref_stratified_results.csv
#          results/ablation/no_norm_comparison_results.csv
#          results/ablation/basic_norm_comparison_results.csv
# Output:  figures/supp_fig_s9_s10_combined.pdf

suppressPackageStartupMessages({ library(here); library(ggplot2); library(patchwork) })

revision_root <- here::here("revision")
strat     <- read.csv(file.path(revision_root, "results/nref/nref_stratified_results.csv"))
cmv_strat <- read.csv(file.path(revision_root, "results/nref/cmv_nref_stratified_results.csv"))
no_norm   <- read.csv(file.path(revision_root, "results/ablation/no_norm_comparison_results.csv"))
basic_norm<- read.csv(file.path(revision_root, "results/ablation/basic_norm_comparison_results.csv"))
no_norm_auc    <- no_norm$point_auc
basic_norm_auc <- basic_norm$point_auc
no_norm_col    <- '#D62728'
basic_norm_col <- '#9467BD'

tier_cols <- c('Top 10%' = '#009E73', 'Bottom 10%' = '#D55E00', 'Random' = '#56B4E9')
strat$tier     <- factor(strat$tier,     levels = c('Top 10%', 'Random', 'Bottom 10%'))
cmv_strat$tier <- factor(cmv_strat$tier, levels = c('Top 10%', 'Random', 'Bottom 10%'))

make_box <- function(data, y_col, y_lab, title, ylims = NULL) {
  p <- ggplot(data, aes(x = tier, y = .data[[y_col]], fill = tier)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.4) +
    scale_fill_manual(values = tier_cols) +
    labs(x = NULL, y = y_lab, title = title) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')
  if (!is.null(ylims)) p <- p + scale_y_continuous(limits = ylims)
  p
}

p_bh_mse <- make_box(strat, 'recon_mse', 'MSE',
                     'BioHEART: VAE reconstruction error') +
  theme(plot.title = element_text(size = 11))

p_bh_auc <- (make_box(strat, 'auc_both', 'AUC',
                       'BioHEART: AUC (both features)', c(0.3, 1.0)) +
             theme(plot.title = element_text(size = 11))) +
  geom_hline(yintercept = no_norm_auc, linetype = 'dashed',
             colour = no_norm_col, linewidth = 0.7) +
  geom_hline(yintercept = basic_norm_auc, linetype = 'dotted',
             colour = basic_norm_col, linewidth = 0.7) +
  annotate('label', x = 0.55, y = 0.97,
           label = sprintf('No norm = %.2f', no_norm_auc),
           hjust = 0, vjust = 0.5, size = 2.8, colour = no_norm_col,
           fill = 'white', label.size = 0.3,
           label.padding = unit(0.12, 'lines')) +
  annotate('label', x = 0.55, y = 0.90,
           label = sprintf('Basic norm = %.2f', basic_norm_auc),
           hjust = 0, vjust = 0.5, size = 2.8, colour = basic_norm_col,
           fill = 'white', label.size = 0.3,
           label.padding = unit(0.12, 'lines'))

p_cmv_mse <- make_box(cmv_strat, 'recon_mse', 'MSE',
                      'CMV: VAE reconstruction error') +
  theme(plot.title = element_text(size = 11))

p_cmv_auc <- make_box(cmv_strat, 'auc_both', 'AUC',
                      'CMV: AUC (both features)', c(0.3, 1.0)) +
  theme(plot.title = element_text(size = 11))

p_s9_s10 <- (p_bh_mse | p_bh_auc) / (p_cmv_mse | p_cmv_auc) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = 'bold', size = 14))

out <- file.path(revision_root, "figures/supp_fig_s9_s10_combined.pdf")
ggsave(out, p_s9_s10, width = 9, height = 8)
cat("Saved:", out, "\n")
