# Supplementary Figure 9: Uncertainty quantification across all four datasets,
# extended to include CellCNN and DeepCNN baselines alongside dioscRi.
# Panel A: 2x2 grid of forest sub-panels, one per dataset, each showing the
#          three methods with point AUC + 95% bootstrap CI.
# Panel B: violin + box plot of bootstrap AUC distributions, 4 dataset columns
#          with 3 dodged method violins per column.

rm(list = ls())
suppressPackageStartupMessages({
  library(here)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(patchwork)
  library(pROC)
})

revision_root <- here::here("revision")
DLC_root <- Sys.getenv("DEEPLEARNING_CYTOF_ROOT",
  unset = "~/Documents/Academic/PhD/DeepLearning_CyTOF")

dataset_order <- c("BioHEART-CT", "Breast Cancer", "COVID-19", "CMV SDY519")
method_order  <- c("dioscRi", "CellCNN", "DeepCNN")

method_cols <- c(dioscRi = "#c0392b",
                 CellCNN = "#f39c12",
                 DeepCNN = "#34495e")

theme_set(theme_classic(base_size = 13))

# --- Bootstrap helper ---
bootstrap_auc_samples <- function(labels, scores, n = 2000, seed = 1994) {
  set.seed(seed)
  pos <- which(labels == 1)
  neg <- which(labels == 0)
  aucs <- numeric(n)
  for (i in seq_len(n)) {
    idx <- c(sample(pos, length(pos), replace = TRUE),
             sample(neg, length(neg), replace = TRUE))
    suppressMessages(
      aucs[i] <- as.numeric(pROC::auc(labels[idx], scores[idx]))
    )
  }
  list(point = suppressMessages(as.numeric(pROC::auc(labels, scores))),
       lo    = unname(quantile(aucs, 0.025)),
       hi    = unname(quantile(aucs, 0.975)),
       aucs  = aucs)
}

load_pred <- function(path) {
  d <- read.csv(path)
  if (all(c("label", "prediction_score") %in% names(d)))
    return(list(labs = d$label, preds = d$prediction_score))
  list(labs = d$truth, preds = d$pred)
}

# --- Baseline prediction files for all 4 datasets ---
baseline_files <- list(
  list(dataset = "BioHEART-CT",   method = "CellCNN",
       path = file.path(DLC_root, "Bioheart_combined/data/bioheart_cell_cnn_aucs.csv")),
  list(dataset = "BioHEART-CT",   method = "DeepCNN",
       path = file.path(DLC_root, "Bioheart_combined/data/bioheart_deep_cytof_aucs.csv")),
  list(dataset = "Breast Cancer", method = "CellCNN",
       path = file.path(DLC_root, "Breast_Cancer_Wagner_2019/data/breast_cancer_cell_cnn_aucs.csv")),
  list(dataset = "Breast Cancer", method = "DeepCNN",
       path = file.path(DLC_root, "Breast_Cancer_Wagner_2019/data/breast_cancer_deep_cytof_aucs.csv")),
  list(dataset = "COVID-19",      method = "CellCNN",
       path = file.path(DLC_root, "COVID_19_PBMC_Mathew_2020/data/bcell_cell_cnn_aucs.csv")),
  list(dataset = "COVID-19",      method = "DeepCNN",
       path = file.path(DLC_root, "COVID_19_PBMC_Mathew_2020/data/bcell_deep_cytof_aucs.csv")),
  list(dataset = "CMV SDY519",    method = "CellCNN",
       path = file.path(DLC_root, "DeepLearningCyTOF/cmv_cell_cnn_aucs.csv")),
  list(dataset = "CMV SDY519",    method = "DeepCNN",
       path = file.path(DLC_root, "DeepLearningCyTOF/cmv_deep_cytof_aucs.csv"))
)

baseline_boots <- lapply(baseline_files, function(e) {
  pr <- load_pred(e$path)
  b  <- bootstrap_auc_samples(pr$labs, pr$preds)
  list(dataset = e$dataset, method = e$method,
       point = b$point, lo = b$lo, hi = b$hi, aucs = b$aucs)
})

# --- dioscRi distributions (already saved) ---
dio_files <- list(
  list(dataset = "BioHEART-CT",
       dist = "results/uncertainty/bioheart_bootstrap_distribution.csv",
       ci   = "results/uncertainty/bioheart_bootstrap_ci.csv"),
  list(dataset = "Breast Cancer",
       dist = "results/uncertainty/breast_cancer_bootstrap_distribution.csv",
       ci   = "results/uncertainty/breast_cancer_bootstrap_ci.csv"),
  list(dataset = "COVID-19",
       dist = "results/uncertainty/covid_bootstrap_distribution.csv",
       ci   = "results/uncertainty/covid_bootstrap_ci.csv"),
  list(dataset = "CMV SDY519",
       dist = "results/uncertainty/cmv_bootstrap_distribution.csv",
       ci   = "results/uncertainty/cmv_bootstrap_ci.csv")
)

dio_boots <- lapply(dio_files, function(e) {
  dist <- read.csv(file.path(revision_root, e$dist))
  ci   <- read.csv(file.path(revision_root, e$ci))
  auc_col <- if ("auc" %in% names(dist)) "auc" else
             names(dist)[sapply(dist, is.numeric)][1]
  list(dataset = e$dataset, method = "dioscRi",
       point = ci$point_auc[1], lo = ci$ci_lower[1], hi = ci$ci_upper[1],
       aucs  = dist[[auc_col]])
})

all_boots <- c(dio_boots, baseline_boots)

# --- Summary table ---
summary_df <- do.call(rbind, lapply(all_boots, function(b)
  tibble(dataset = b$dataset, method = b$method,
         point = b$point, lo = b$lo, hi = b$hi)))
summary_df$dataset <- factor(summary_df$dataset, levels = dataset_order)
summary_df$method  <- factor(summary_df$method,  levels = method_order)
summary_df$ci_text <- sprintf("%.3f [%.3f, %.3f]",
                              summary_df$point, summary_df$lo, summary_df$hi)

cat("\n--- Summary (point AUC and 95% bootstrap CI) ---\n")
print(summary_df %>% arrange(dataset, method) %>%
        mutate(across(where(is.numeric), ~round(., 3))))

# --- Helper: per-dataset mini forest panel ---
# Label format: "point [lo, hi]" to the right of the upper CI bound.
# Extra horizontal room reserved past 1.0 so labels never clip.
make_forest <- function(ds_name) {
  d <- summary_df %>% filter(dataset == ds_name)
  d$method_y <- factor(d$method, levels = rev(method_order))

  ggplot(d, aes(x = point, y = method_y, colour = method)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey60") +
    geom_errorbar(aes(xmin = lo, xmax = hi),
                  width = 0.22, linewidth = 0.7, orientation = "y") +
    geom_point(size = 3.0) +
    geom_text(aes(label = ci_text, x = hi + 0.015),
              hjust = 0, size = 2.4, colour = "grey15") +
    scale_colour_manual(values = method_cols) +
    scale_x_continuous(limits = c(0.30, 1.55),
                       breaks = seq(0.4, 1.0, 0.2),
                       expand = expansion(mult = 0)) +
    scale_y_discrete(limits = rev(method_order),
                     expand = expansion(add = c(0.5, 0.5))) +
    labs(title = ds_name, x = "AUC", y = NULL) +
    theme(legend.position = "none",
          plot.title   = element_text(face = "bold", size = 12, hjust = 0.5),
          axis.text.y  = element_text(size = 10, face = "bold"),
          axis.title.x = element_text(size = 11),
          plot.margin  = margin(8, 12, 8, 8),
          panel.grid.major.y = element_blank())
}

panel_a <- (make_forest("BioHEART-CT")   | make_forest("Breast Cancer")) /
           (make_forest("COVID-19")      | make_forest("CMV SDY519"))   +
           plot_annotation(
             title = "A) Point AUC with 95% bootstrap CI (n = 2000 stratified resamples)",
             theme = theme(plot.title = element_text(face = "bold", size = 13))
           )

# --- Panel B: violin plot (datasets x methods) ---
dist_df <- do.call(rbind, lapply(all_boots, function(b)
  tibble(dataset = b$dataset, method = b$method, auc = b$aucs)))
dist_df$dataset <- factor(dist_df$dataset, levels = dataset_order)
dist_df$method  <- factor(dist_df$method,  levels = method_order)

panel_b <- ggplot(dist_df, aes(x = dataset, y = auc,
                                fill = method, colour = method)) +
  geom_violin(alpha = 0.40, linewidth = 0.4, trim = FALSE,
              position = position_dodge(width = 0.8)) +
  geom_boxplot(width = 0.10, fill = "white", outlier.shape = NA,
               linewidth = 0.4,
               position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = method_cols) +
  scale_colour_manual(values = method_cols) +
  scale_y_continuous(limits = c(0.30, 1.05),
                     breaks = seq(0.4, 1.0, 0.1)) +
  labs(title = "B) Bootstrap AUC distribution",
       x = NULL, y = "AUC", fill = NULL, colour = NULL) +
  guides(fill = guide_legend(override.aes = list(alpha = 0.7)),
         colour = "none") +
  theme(plot.title  = element_text(face = "bold", size = 13),
        axis.text.x = element_text(size = 12, face = "bold",
                                   angle = 30, hjust = 1),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", colour = NA),
        panel.grid.minor.y = element_blank())

p_combined <- (wrap_elements(panel_a) / panel_b) +
              plot_layout(heights = c(1.0, 1.0))

ggsave(file.path(revision_root, "figures/supp_uncertainty_quantification.pdf"),
       p_combined, width = 9, height = 9)

cat("\nFigure saved to figures/supp_uncertainty_quantification.pdf\n")
