# Supplementary figure: cell-to-patient method comparison on BioHEART-CT.
# Compares cytoGPNet against dioscRi, CellCNN, and DeepCNN. dioscRi and
# cytoGPNet have stratified bootstrap CIs computed from their test predictions;
# CellCNN/DeepCNN are point estimates from Figure 4 of the main manuscript.
#
# Inputs:
#   results/cytogpnet/bioheart/model_output/test_predictions1.csv
#   results/uncertainty/bioheart_bootstrap_distribution.csv
#
# Output:
#   figures/supp_fig_cytogpnet_comparison.pdf

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ggplot2)
  library(pROC)
})

revision_root <- here::here("revision")

# ---- AUC + stratified bootstrap CI from predictions ----
bootstrap_auc <- function(labels, scores, n = 2000, seed = 1994) {
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
  point <- suppressMessages(as.numeric(pROC::auc(labels, scores)))
  list(point = point,
       lo    = unname(quantile(aucs, 0.025)),
       hi    = unname(quantile(aucs, 0.975)))
}

# ---- cytoGPNet: from test_predictions1.csv ----
cyto <- read.csv(file.path(revision_root,
  "results/cytogpnet/bioheart/model_output/test_predictions1.csv"))
cyto_ci <- bootstrap_auc(cyto$label, cyto$prediction_score)

# ---- dioscRi: reuse the existing bootstrap distribution (already 2000 resamples) ----
dio_dist <- read.csv(file.path(revision_root,
  "results/uncertainty/bioheart_bootstrap_distribution.csv"))
dio_col <- {
  if ("bioheart" %in% names(dio_dist)) {
    "bioheart"
  } else if ("auc" %in% names(dio_dist)) {
    "auc"
  } else {
    names(dio_dist)[sapply(dio_dist, is.numeric)][1]
  }
}
dio_aucs <- dio_dist[[dio_col]]
dio_ci   <- list(point = 0.825,
                 lo    = unname(quantile(dio_aucs, 0.025)),
                 hi    = unname(quantile(dio_aucs, 0.975)))

# ---- CellCNN / DeepCNN: point estimates from main-manuscript Figure 4 ----
res <- tibble::tibble(
  method = c("dioscRi", "cytoGPNet", "CellCNN", "DeepCNN"),
  auc    = c(dio_ci$point,  cyto_ci$point, 0.60, 0.62),
  lo     = c(dio_ci$lo,     cyto_ci$lo,   NA,   NA),
  hi     = c(dio_ci$hi,     cyto_ci$hi,   NA,   NA),
  has_ci = c(TRUE, TRUE, FALSE, FALSE)
)

# Order top-to-bottom: dioscRi, cytoGPNet, DeepCNN, CellCNN.
# scale_y_discrete(limits = ...) places the first listed level at the bottom.
method_order_bottom_to_top <- c("CellCNN", "DeepCNN", "cytoGPNet", "dioscRi")
res$method <- factor(res$method, levels = method_order_bottom_to_top)

# Label x-position = right end of CI (where present) or the point itself, plus a margin
res$label    <- ifelse(res$has_ci,
                       sprintf("%.3f [%.3f, %.3f]", res$auc, res$lo, res$hi),
                       sprintf("%.3f (point estimate)", res$auc))
res$label_x  <- ifelse(res$has_ci, res$hi, res$auc) + 0.015

cat("\n=== AUCs for cytoGPNet comparison figure ===\n")
print(res)

# ---- Forest plot ----
p <- ggplot(res, aes(x = auc, y = method, colour = method)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(data = subset(res, has_ci),
                 aes(xmin = lo, xmax = hi), height = 0.18, linewidth = 0.7) +
  geom_point(size = 3.5) +
  geom_text(aes(x = label_x, label = label),
            hjust = 0, size = 3.4, colour = "black") +
  scale_colour_manual(values = c(dioscRi   = "#c0392b",
                                  cytoGPNet = "#27ae60",
                                  CellCNN   = "#f39c12",
                                  DeepCNN   = "#34495e")) +
  scale_x_continuous(limits = c(0.35, 1.25), breaks = seq(0.4, 1.0, 0.1)) +
  scale_y_discrete(limits = method_order_bottom_to_top) +
  labs(title = "Cell-to-patient method comparison on BioHEART-CT",
       subtitle = "Point AUC with 95% stratified bootstrap CI (n = 2000 resamples) for dioscRi and cytoGPNet;\nCellCNN and DeepCNN are point estimates from Figure 4 (no resampling).",
       x = "AUC", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey30", size = 9, lineheight = 1.1),
        panel.grid.major.y = element_blank())

out <- file.path(revision_root, "figures/supp_fig_cytogpnet_comparison.pdf")
ggsave(out, p, width = 7.5, height = 3.6)
cat("\nSaved:", out, "\n")

# ---- Persist the comparison table alongside ----
write.csv(res, file.path(revision_root, "results/cytogpnet/method_comparison_bioheart.csv"),
          row.names = FALSE)
