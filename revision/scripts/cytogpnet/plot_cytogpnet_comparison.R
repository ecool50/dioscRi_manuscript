# Response Figure R1: cell-to-patient method comparison on BioHEART-CT.
# Two-panel layout matching Supp Fig 9 of the manuscript:
#   A) Forest plot with point AUC + 95% stratified bootstrap CI per method.
#   B) Violin + box plot of the bootstrap AUC distributions per method.
# All four methods (dioscRi, cytoGPNet, DeepCNN, CellCNN) have per-sample
# test predictions, so each gets the same 2000-resample stratified bootstrap.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(pROC)
})

revision_root <- here::here("revision")
# CellCNN / DeepCNN per-sample predictions on BioHEART-CT are shipped with
# the repo under revision/results/baselines/BioHEART_combined/.
baselines_root <- file.path(revision_root, "results", "baselines")

# ---- bootstrap AUC distribution ----
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

# ---- cytoGPNet (revision results) ----
cyto <- read.csv(file.path(revision_root,
  "results/cytogpnet/bioheart/model_output/test_predictions1.csv"))
cyto_b <- bootstrap_auc_samples(cyto$label, cyto$prediction_score)

# ---- CellCNN per-sample test predictions ----
cc <- read.csv(file.path(baselines_root,
  "BioHEART_combined/bioheart_cell_cnn_aucs.csv"))
cc_b <- bootstrap_auc_samples(cc$truth, cc$pred)

# ---- DeepCNN (a.k.a. DL CyTOF) per-sample test predictions ----
dc <- read.csv(file.path(baselines_root,
  "BioHEART_combined/bioheart_deep_cytof_aucs.csv"))
dc_b <- bootstrap_auc_samples(dc$truth, dc$pred)

# ---- dioscRi: existing 2000-resample bootstrap distribution + point AUC ----
dio_dist <- read.csv(file.path(revision_root,
  "results/uncertainty/bioheart_bootstrap_distribution.csv"))
dio_col <- {
  if ("bioheart" %in% names(dio_dist)) "bioheart"
  else if ("auc" %in% names(dio_dist)) "auc"
  else names(dio_dist)[sapply(dio_dist, is.numeric)][1]
}
dio_aucs <- dio_dist[[dio_col]]
dio_b <- list(point = 0.825,
              lo    = unname(quantile(dio_aucs, 0.025)),
              hi    = unname(quantile(dio_aucs, 0.975)),
              aucs  = dio_aucs)

# ---- combine ----
boot <- list(dioscRi = dio_b, cytoGPNet = cyto_b, DeepCNN = dc_b, CellCNN = cc_b)

methods_top_to_bottom <- c("dioscRi", "cytoGPNet", "DeepCNN", "CellCNN")
methods_left_to_right <- methods_top_to_bottom

col <- c(dioscRi   = "#c0392b",
        cytoGPNet  = "#27ae60",
        DeepCNN    = "#34495e",
        CellCNN    = "#f39c12")

res <- tibble::tibble(
  method = methods_top_to_bottom,
  auc    = vapply(boot[methods_top_to_bottom], `[[`, numeric(1), "point"),
  lo     = vapply(boot[methods_top_to_bottom], `[[`, numeric(1), "lo"),
  hi     = vapply(boot[methods_top_to_bottom], `[[`, numeric(1), "hi")
)
res$method  <- factor(res$method, levels = rev(methods_top_to_bottom))   # bottom->top
res$label   <- sprintf("%.3f [%.3f, %.3f]", res$auc, res$lo, res$hi)
res$label_x <- res$hi + 0.015

cat("\n=== AUCs for cytoGPNet comparison figure ===\n")
print(res)

# ---- Panel A: forest plot ----
panel_a <- ggplot(res, aes(x = auc, y = method, colour = method)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey60") +
  geom_errorbar(aes(xmin = lo, xmax = hi),
                width = 0.25, linewidth = 0.7, orientation = "y") +
  geom_point(size = 3.4) +
  geom_text(aes(x = label_x, label = label),
            hjust = 0, size = 3.5, colour = "grey15") +
  scale_colour_manual(values = col) +
  scale_x_continuous(limits = c(0.20, 1.18),
                     breaks = seq(0.2, 1.0, 0.1),
                     expand = expansion(mult = 0)) +
  scale_y_discrete(limits = rev(methods_top_to_bottom)) +
  labs(title = "A) Point AUC with 95% bootstrap CI (n = 2000 stratified resamples)",
       x = "AUC", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none",
        plot.title  = element_text(face = "bold", size = 13),
        axis.text.y = element_text(size = 12, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank())

# ---- Panel B: violin + box for all four methods ----
dist_df <- bind_rows(
  tibble::tibble(method = "dioscRi",   auc = dio_b$aucs),
  tibble::tibble(method = "cytoGPNet", auc = cyto_b$aucs),
  tibble::tibble(method = "DeepCNN",   auc = dc_b$aucs),
  tibble::tibble(method = "CellCNN",   auc = cc_b$aucs)
) |>
  mutate(method = factor(method, levels = methods_left_to_right))

panel_b <- ggplot(dist_df, aes(x = method, y = auc,
                               fill = method, colour = method)) +
  geom_violin(alpha = 0.45, linewidth = 0.5, trim = FALSE) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA,
               linewidth = 0.5) +
  scale_fill_manual(values = col) +
  scale_colour_manual(values = col) +
  scale_y_continuous(limits = c(0.15, 1.05),
                     breaks = seq(0.2, 1.0, 0.1)) +
  scale_x_discrete(limits = methods_left_to_right) +
  labs(title = "B) Bootstrap AUC distribution",
       x = NULL, y = "AUC") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none",
        plot.title  = element_text(face = "bold", size = 13),
        axis.text.x = element_text(size = 12, face = "bold",
                                   angle = 30, hjust = 1),
        panel.grid.minor.y = element_blank())

fig <- panel_a / panel_b + plot_layout(heights = c(1, 1.3))

out <- file.path(revision_root, "figures/supp_fig_cytogpnet_comparison.pdf")
ggsave(out, fig, width = 9, height = 6.5)
cat("\nSaved:", out, "\n")

write.csv(res, file.path(revision_root, "results/cytogpnet/method_comparison_bioheart.csv"),
          row.names = FALSE)
