# BioHEART-CT Ablation Study
# Systematically removes/replaces pipeline components to isolate contributions.
#
# Ablations:
#   1. Full pipeline (baseline)
#   2. No normalization (raw data)
#   3. No hierarchy (flat groups)
#   4. Elastic net (replaces overlapping group lasso)
#   5. Proportions only
#   6. Means only

rm(list = ls())
library(dioscRi)
library(ggplot2)
library(FuseSOM)
library(SingleCellExperiment)
library(tidyverse)
library(data.table)
library(data.tree)
library(caret)
library(pROC)
library(glmnet)
tensorflow::set_random_seed(1994)

# --- Config ---
nclust <- 11
bioheart_root <- Sys.getenv("BIOHEART_ROOT", unset = "~/Documents/Academic/PhD/bioheart_analysis")
revision_root <- here::here("revision")
base_dir <- bioheart_root  # alias: upstream inputs (raw data, sce objects, etc.)
useMarkers <- c('HLA-DR', 'CD3', 'CD4', 'CD8a', 'CD25', 'CD127', 'FoxP3', 'CD27',
                'KLRG1', 'CD56', 'CD45RO', 'CD45RA', 'CD192-CCR2', 'CD194-CCR4',
                'CD196-CCR6',
                'CD39', 'CD38', 'Ki67', 'CD183-CXCR3', 'CCR7', 'CD19', 'CD20',
                'IgD', 'CD14', 'CD304', 'CD141', 'CD1c-PE')
cluster_col <- "clusters_11"

# --- Load data ---
cat("Loading data...\n")
load(file.path(base_dir, 'data/sce_dat/study_4_MMD_VAE_updated.RData'))  # sce_norm
load(file.path(base_dir, 'data/sce_dat/study_3_MMD_VAE_updated.RData'))  # sce_3_norm

# --- Clinical data ---
clinicaldata <- colData(sce_norm) %>%
  as.data.frame() %>%
  dplyr::select(sample_id, Gensini) %>%
  distinct()
clinicaldata$Gensini_bin <- factor(if_else(clinicaldata$Gensini > 0, 1, 0))

clinicaldata_3 <- colData(sce_3_norm) %>%
  as.data.frame() %>%
  dplyr::select(sample_id, Gensini_bin) %>%
  distinct()
clinicaldata_3$Gensini_bin <- as.factor(clinicaldata_3$Gensini_bin)

# --- Helper: run prediction pipeline with options ---
run_ablation <- function(sce_train, sce_test, cd_train, cd_test,
                         assay_name = 'norm',
                         use_hierarchy = TRUE,
                         model_type = 'grLasso',  # 'grLasso' or 'elasticnet'
                         feature_type = 'both',    # 'both', 'prop', 'mean'
                         cluster_col_train = 'clusters_11',
                         cluster_col_test = 'clusters') {

  # Use existing clusters_11 if available and using norm assay (matches original pipeline),
  # otherwise re-cluster (e.g., for raw data ablation)
  if (assay_name == 'norm' && 'clusters_11' %in% names(colData(sce_train))) {
    sce_train$clusters_abl <- sce_train$clusters_11
    cat("    Using existing clusters_11\n")
  } else {
    sce_train <- runFuseSOM(sce_train, numClusters = nclust, assay = assay_name,
                            verbose = FALSE, clusterCol = 'clusters_abl')
    cat("    Re-clustering on", assay_name, "assay\n")
  }

  # LDA on VAE embeddings
  train_x <- reducedDim(sce_train, type = "VAE") %>% as.data.frame()
  train_x$cellTypes <- as.factor(sce_train$clusters_abl)
  test_x <- reducedDim(sce_test, type = "VAE") %>% as.data.frame()
  lda_fit <- MASS::lda(cellTypes ~ ., data = train_x)
  test_clusters <- predict(lda_fit, test_x)$class
  sce_test$clusters_abl <- as.factor(test_clusters)

  # Features
  prop_train <- computeFeatures(sce = sce_train, featureType = 'prop',
                                cellTypeCol = 'clusters_abl', sampleCol = 'sample_id',
                                logit = TRUE, useMarkers = useMarkers)
  means_train <- computeFeatures(sce = sce_train, featureType = 'mean',
                                 cellTypeCol = 'clusters_abl', sampleCol = 'sample_id',
                                 logit = TRUE, useMarkers = useMarkers)
  colnames(means_train) <- as.factor(colnames(means_train))

  prop_test <- computeFeatures(sce = sce_test, featureType = 'prop',
                               cellTypeCol = 'clusters_abl', sampleCol = 'sample_id',
                               logit = TRUE, useMarkers = useMarkers)
  means_test <- computeFeatures(sce = sce_test, featureType = 'mean',
                                cellTypeCol = 'clusters_abl', sampleCol = 'sample_id',
                                logit = TRUE, useMarkers = useMarkers)
  colnames(prop_test) <- as.factor(colnames(prop_test))
  colnames(means_test) <- as.factor(colnames(means_test))

  # Response
  row_train <- rownames(prop_train)
  condition_train <- factor(cd_train[match(row_train, cd_train$sample_id), "Gensini_bin"])
  y_train <- as.numeric(levels(condition_train))[condition_train]

  row_test <- rownames(prop_test)
  condition_test <- factor(cd_test[match(row_test, cd_test$sample_id), "Gensini_bin"])
  y_test <- as.numeric(levels(condition_test))[condition_test]

  # Select features based on type
  if (feature_type == 'prop') {
    X_train <- prop_train
    X_test <- prop_test
  } else if (feature_type == 'mean') {
    X_train <- means_train
    X_test <- means_test
  } else {
    X_train <- cbind(prop_train, means_train)
    X_test <- cbind(prop_test, means_test)
  }

  # Scale
  scaleVals <- preProcess(X_train, method = c('scale'))
  X_train <- predict(scaleVals, X_train) %>% as.matrix()
  X_test <- predict(scaleVals, X_test) %>% as.matrix()

  if (model_type == 'elasticnet') {
    # Elastic net via glmnet
    cv_fit <- cv.glmnet(X_train, y_train, family = 'binomial', alpha = 0.5, nfolds = 5)
    preds <- predict(cv_fit, X_test, s = 'lambda.min', type = 'response')[, 1]
    auc_val <- as.numeric(auc(response = factor(y_test), predictor = preds, quiet = TRUE))
    preds_df <- data.frame(Truth = factor(y_test), Predicted = preds)
  } else {
    # Group lasso
    if (use_hierarchy && feature_type != 'mean') {
      trees <- generateTree(features = prop_train, method = "ward")
      tree <- trees$tree
      if (feature_type == 'prop') {
        groups <- generateGroups(tree = tree, nClust = nclust,
                                proportions = prop_train, means = means_train,
                                type = "prop")
      } else {
        groups <- generateGroups(tree = tree, nClust = nclust,
                                proportions = prop_train, means = means_train)
      }
    } else {
      # Flat groups: each feature is its own group
      groups <- as.list(seq(1, ncol(X_train), 1))
    }
    groups <- lapply(groups, sort)

    fit <- fitModel(xTrain = X_train, yTrain = y_train,
                    groups = groups, penalty = 'grLasso', seed = 1994,
                    modelSummary = FALSE)
    test_auc <- plotAUC(fit = fit$fit, xTest = X_test, yTest = y_test, title = "")
    auc_val <- as.numeric(auc(test_auc$preds$Truth, test_auc$preds$Predicted, quiet = TRUE))
    preds_df <- test_auc$preds
  }

  # Bootstrap CI (stratified, 2000 resamples)
  set.seed(1994)
  roc_obj <- roc(preds_df$Truth, preds_df$Predicted, quiet = TRUE)
  ci_result <- ci.auc(roc_obj, method = "bootstrap", boot.n = 2000,
                      boot.stratified = TRUE, conf.level = 0.95)

  return(list(auc = auc_val, ci_lower = as.numeric(ci_result[1]),
              ci_upper = as.numeric(ci_result[3])))
}

# --- Run all combinations ---
grid <- expand.grid(
  assay_name = c('norm', 'raw'),
  use_hierarchy = c(TRUE, FALSE),
  model_type = c('grLasso', 'elasticnet'),
  feature_type = c('both', 'prop', 'mean'),
  stringsAsFactors = FALSE
)

# Label each combination
grid$label <- paste0(
  if_else(grid$assay_name == 'norm', 'Norm', 'Raw'), ' + ',
  if_else(grid$use_hierarchy, 'Hierarchy', 'Flat'), ' + ',
  if_else(grid$model_type == 'grLasso', 'GroupLasso', 'ElasticNet'), ' + ',
  grid$feature_type
)

results <- data.frame(
  ablation = character(), assay = character(), hierarchy = character(),
  model = character(), features = character(), auc = numeric(),
  ci_lower = numeric(), ci_upper = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(grid)) {
  g <- grid[i, ]
  cat(sprintf("Running %d/%d: %s...\n", i, nrow(grid), g$label))

  res <- tryCatch(
    run_ablation(sce_norm, sce_3_norm, clinicaldata, clinicaldata_3,
                 assay_name = g$assay_name, use_hierarchy = g$use_hierarchy,
                 model_type = g$model_type, feature_type = g$feature_type),
    error = function(e) { cat("  ERROR:", e$message, "\n"); list(auc = NA, ci_lower = NA, ci_upper = NA) }
  )

  results <- rbind(results, data.frame(
    ablation = g$label, assay = g$assay_name, hierarchy = g$use_hierarchy,
    model = g$model_type, features = g$feature_type,
    auc = res$auc, ci_lower = res$ci_lower, ci_upper = res$ci_upper
  ))
  cat(sprintf("  AUC = %.4f [%.4f, %.4f]\n", res$auc, res$ci_lower, res$ci_upper))
}

# Mark the full pipeline
results$is_full <- results$assay == 'norm' & results$hierarchy == TRUE &
                   results$model == 'grLasso' & results$features == 'both'
# --- Summary ---
cat("\n--- BioHEART-CT Ablation Results (all 24 combinations) ---\n")
print(results %>% arrange(desc(auc)))

write.csv(results, file.path(revision_root, 'results/ablation/bioheart_ablation_results.csv'), row.names = FALSE)
cat("\nSaved to data/bioheart_ablation_results.csv\n")

# --- Plot ---
library(patchwork)

# Order by AUC
results <- results %>% arrange(auc)
row_order <- results$ablation

# --- Plot ---
# Faceted by normalization, short labels, CI as text

# Short label: Hierarchy/Flat + GroupLasso/ElasticNet + Both/Prop/Mean
results$short_label <- paste0(
  if_else(results$hierarchy == TRUE, 'Hier', 'Flat'), ' + ',
  if_else(results$model == 'grLasso', 'GrLasso', 'ENet'), ' + ',
  results$features
)
results$norm_label <- if_else(results$assay == 'norm', 'MMD-VAE Normalized', 'Raw (no normalization)')

# AUC + CI text label
results$ci_label <- sprintf("%.3f [%.3f, %.3f]", results$auc, results$ci_lower, results$ci_upper)

# Order by AUC within each facet
results <- results %>%
  group_by(norm_label) %>%
  arrange(auc, .by_group = TRUE) %>%
  ungroup()
results$short_label <- factor(results$short_label,
                               levels = unique(results$short_label))

results$highlight <- case_when(
  results$is_full ~ 'Full pipeline',
  results$assay == 'norm' ~ 'Normalized',
  TRUE ~ 'Raw'
)

p <- ggplot(results, aes(x = short_label, y = auc, fill = highlight)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = ci_label), hjust = -0.05, size = 2.8) +
  geom_hline(yintercept = 0.5, linetype = 'dashed', colour = 'grey60') +
  scale_fill_manual(values = c('Full pipeline' = '#8B0000',
                                'Normalized' = '#4A90D9',
                                'Raw' = '#B0B0B0')) +
  scale_y_continuous(limits = c(0, 1.35), breaks = seq(0, 1, 0.2)) +
  coord_flip() +
  facet_wrap(~ norm_label, scales = 'free_y', ncol = 1) +
  labs(x = NULL, y = 'AUC', fill = NULL,
       title = 'BioHEART-CT Ablation Study',
       subtitle = 'AUC with 95% bootstrap confidence intervals') +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = 'bold'),
        plot.subtitle = element_text(size = 11, colour = 'grey40'),
        legend.position = 'none',
        strip.text = element_text(face = 'bold', size = 12))

ggsave(file.path(revision_root, 'figures/supp_ablation_bioheart.pdf'),
       p, width = 10, height = 10)
cat("Figure saved to figures/supp_ablation_bioheart.pdf\n")
