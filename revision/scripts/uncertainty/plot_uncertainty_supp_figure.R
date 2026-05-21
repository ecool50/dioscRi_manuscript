# Supplementary Figure: Uncertainty Quantification Across All Datasets
# All four datasets use the same method: bootstrap CI on test predictions
# (2000 stratified resamples). Forest plot + bootstrap distribution.
# Run after the four bootstrap CI scripts have completed.

rm(list = ls())
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

bioheart_root <- Sys.getenv("BIOHEART_ROOT", unset = "~/Documents/Academic/PhD/bioheart_analysis")
revision_root <- here::here("revision")
base_dir <- bioheart_root  # alias: upstream inputs (raw data, sce objects, etc.)
dataset_order <- c('BioHEART-CT', 'Breast Cancer', 'COVID-19', 'CMV SDY519')

dataset_cols <- c('BioHEART-CT'   = '#8B0000',
                  'Breast Cancer' = '#006400',
                  'COVID-19'      = '#FF8C00',
                  'CMV SDY519'    = '#4A4A4A')

theme_set(theme_classic(base_size = 14))

# --- Load bootstrap distributions (all four datasets, same method) ---

bioheart_boot <- read.csv(file.path(revision_root, 'results/uncertainty/bioheart_bootstrap_distribution.csv'))
bioheart_boot$dataset <- 'BioHEART-CT'

breast_boot <- read.csv(file.path(revision_root, 'results/uncertainty/breast_cancer_bootstrap_distribution.csv'))
breast_boot$dataset <- 'Breast Cancer'

covid_boot <- read.csv(file.path(revision_root, 'results/uncertainty/covid_bootstrap_distribution.csv'))
covid_boot$dataset <- 'COVID-19'

cmv_boot <- read.csv(file.path(revision_root, 'results/uncertainty/cmv_bootstrap_distribution.csv'))
cmv_boot$dataset <- 'CMV SDY519'

# --- Load point estimates + CIs (from bootstrap CI scripts) ---

bioheart_ci <- read.csv(file.path(revision_root, 'results/uncertainty/bioheart_bootstrap_ci.csv'))
bioheart_ci$dataset <- 'BioHEART-CT'

breast_ci <- read.csv(file.path(revision_root, 'results/uncertainty/breast_cancer_bootstrap_ci.csv'))
breast_ci$dataset <- 'Breast Cancer'

covid_ci <- read.csv(file.path(revision_root, 'results/uncertainty/covid_bootstrap_ci.csv'))
covid_ci$dataset <- 'COVID-19'

cmv_ci <- read.csv(file.path(revision_root, 'results/uncertainty/cmv_bootstrap_ci.csv'))
cmv_ci$dataset <- 'CMV SDY519'

summary_df <- bind_rows(bioheart_ci, breast_ci, covid_ci, cmv_ci) %>%
  dplyr::select(dataset, point_auc, ci_lower, ci_upper)

summary_df$dataset <- factor(summary_df$dataset, levels = dataset_order)
summary_df <- summary_df %>% arrange(dataset)

summary_df$ci_text <- sprintf("%.3f [%.3f, %.3f]",
                               summary_df$point_auc,
                               summary_df$ci_lower,
                               summary_df$ci_upper)

# --- Panel A: Forest plot (point AUC + 95% CI) ---

# Reverse for top-to-bottom display order
summary_df$dataset_y <- factor(summary_df$dataset, levels = rev(dataset_order))

p_forest <- ggplot(summary_df, aes(x = point_auc, y = dataset_y, colour = dataset)) +
  geom_point(size = 4) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                 height = 0.25, linewidth = 0.8) +
  geom_text(aes(label = ci_text, x = ci_upper + 0.02),
            hjust = 0, size = 3.5, colour = 'black') +
  geom_vline(xintercept = 0.5, linetype = 'dashed', colour = 'grey60') +
  scale_colour_manual(values = dataset_cols) +
  scale_x_continuous(limits = c(0.4, 1.15), breaks = seq(0.4, 1.0, 0.1)) +
  labs(x = 'AUC', y = NULL,
       title = 'A) Point AUC with 95% bootstrap CI (n = 2000 stratified resamples)') +
  theme(plot.title = element_text(face = 'bold', size = 13),
        legend.position = 'none')

# --- Panel B: Bootstrap distributions (violin) ---

all_boot <- bind_rows(bioheart_boot, breast_boot, covid_boot, cmv_boot)
all_boot$dataset <- factor(all_boot$dataset, levels = dataset_order)

p_violin <- ggplot(all_boot, aes(x = dataset, y = auc, fill = dataset)) +
  geom_violin(alpha = 0.7, scale = 'width', trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = 'white', alpha = 0.6) +
  scale_fill_manual(values = dataset_cols) +
  scale_y_continuous(limits = c(0.4, 1.05), breaks = seq(0.4, 1.0, 0.1)) +
  labs(x = NULL, y = 'AUC',
       title = 'B) Bootstrap AUC distribution') +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title = element_text(face = 'bold', size = 13),
        legend.position = 'none')

# --- Combine ---

p_combined <- p_forest / p_violin + plot_layout(heights = c(1, 1.3))

ggsave(file.path(revision_root, 'figures/supp_uncertainty_quantification.pdf'),
       p_combined, width = 9, height = 6.5)

cat("Figure saved to figures/supp_uncertainty_quantification.pdf\n")

# --- Print summary table ---
cat("\n--- Summary Table (uniform: bootstrap CI on test predictions) ---\n")
print(summary_df %>%
        dplyr::select(dataset, point_auc, ci_lower, ci_upper) %>%
        mutate(across(where(is.numeric), ~round(., 3))))
