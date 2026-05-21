# BioHEART-CT Downsampling Sensitivity [R2-4]
# Varies cells per sample (1K, 2K, 5K, 10K, 20K), retrains the full pipeline
# at each level, and reports ARI and AUC.

# Pin Python environment
Sys.setenv(RETICULATE_PYTHON = '~/.virtualenvs/r-reticulate/bin/python')
library(reticulate)

rm(list = ls())
library(dioscRi)
library(ggplot2)
library(FuseSOM)
library(SingleCellExperiment)
library(tidyverse)
library(data.table)
library(data.tree)
library(tensorflow)
library(keras3)
library(caret)
library(pROC)
library(mclust)  # for adjustedRandIndex

stopifnot(packageVersion('keras3') == '1.2.0')
tf <- import('tensorflow')
stopifnot(tf$version$VERSION == '2.16.2')
cat('Versions OK: keras3 1.2.0, TF 2.16.2\n')

tensorflow::set_random_seed(1994)

# --- Config ---
cell_levels <- c(1000, 2000, 5000, 10000, 20000)
nclust <- 11
nRefSamples <- 2
bioheart_root <- Sys.getenv("BIOHEART_ROOT", unset = "~/Documents/Academic/PhD/bioheart_analysis")
revision_root <- here::here("revision")
base_dir <- bioheart_root  # alias: upstream inputs (raw data, sce objects, etc.)
useMarkers <- c('HLA_DR', 'CD3', 'CD4', 'CD8a', 'CD25', 'CD127', 'FoxP3', 'CD27',
                'KLRG1', 'CD56', 'CD45RO', 'CD45RA', 'CD192_CCR2', 'CD194_CCR4',
                'CD196_CCR6', 'CD39', 'CD38', 'Ki67', 'CD183_CXCR3', 'CCR7',
                'CD19', 'CD20', 'IgD', 'CD14', 'CD304', 'CD141', 'CD1c_PE')

# --- Load 20K data (superset of 10K) ---
cat("Loading BioHEART-CT 20K data...\n")
df_full <- fread(file.path(base_dir, 'data/raw_data/bioheart_ct_cytof_data_b4_mg_20K.csv'),
                 nThread = 7) %>% as.data.frame()

# Also load validation cohort (study 3) from raw CSV - fixed, not downsampled
cat("Loading validation cohort...\n")
df_val_full <- fread(file.path(base_dir, 'misc/Study_3_2019/all_batches_processed_mg_10K.csv'),
                     nThread = 7) %>% as.data.frame()

# Validation clinical data
clinicaldata_3 <- df_val_full %>%
  dplyr::select(sample_id, Gensini_bin) %>%
  distinct()
clinicaldata_3$Gensini_bin <- as.factor(clinicaldata_3$Gensini_bin)

# Discovery clinical data
clinicaldata <- df_full %>%
  dplyr::select(sample_id, Gensini) %>%
  distinct()
clinicaldata$Gensini_bin <- factor(if_else(clinicaldata$Gensini > 0, 1, 0))

# --- Run pipeline at each downsampling level ---
results <- data.frame(cells_per_sample = integer(), ari = numeric(), auc = numeric())

for (n_cells in cell_levels) {
  cat(sprintf("\n=== %dK cells per sample ===\n", n_cells / 1000))
  set.seed(1994)

  # Subsample: shuffle then take up to n_cells per sample
  set.seed(1994)
  df_sub <- df_full[sample(nrow(df_full)), ] %>%
    group_by(sample_id) %>%
    slice_head(n = n_cells) %>%
    ungroup() %>%
    as.data.frame()

  cat(sprintf("  Total cells: %d\n", nrow(df_sub)))

  # Preprocess: arcsinh + min-max
  df_sub[, useMarkers] <- cyCombine::transform_asinh(df_sub[, useMarkers],
                                                      markers = useMarkers, derand = FALSE)
  preProcValues <- preProcess(df_sub[, useMarkers], method = c("range"))
  df_sub[, useMarkers] <- predict(preProcValues, df_sub[, useMarkers])

  # Reference sample selection
  cat("  Selecting reference samples...\n")
  res_ind <- computeReferenceSample(df_sub, useMarkers, sampleCol = 'sample_id', N = nRefSamples)

  # Train VAE
  cat("  Training VAE...\n")
  train_dat <- df_sub[which(df_sub$sample_id %in% res_ind$topNSamples), ]
  val_dat <- df_sub[which(df_sub$sample_id %in% res_ind$bottomNSamples), ]

  vae_model <- trainVAEModel(trainData = train_dat[, useMarkers],
                             valData = val_dat,
                             batchSize = 32L, lambda = 0.01,
                             useMarkers = useMarkers, epochs = 100)

  # Decode discovery
  cat("  Decoding discovery...\n")
  train_decoded <- decodeSamples(newSamples = as.matrix(df_sub[, useMarkers]),
                                  vae = vae_model$vae)
  df_norm <- train_decoded[[1]]
  colnames(df_norm) <- useMarkers

  # Preprocess validation with same scaling as discovery
  cat("  Preprocessing validation...\n")
  df_val <- df_val_full
  df_val[, useMarkers] <- cyCombine::transform_asinh(df_val[, useMarkers],
                                                      markers = useMarkers, derand = FALSE)
  df_val[, useMarkers] <- predict(preProcValues, df_val[, useMarkers])

  # Decode validation
  cat("  Decoding validation...\n")
  val_decoded <- decodeSamples(newSamples = as.matrix(df_val[, useMarkers]),
                                vae = vae_model$vae)
  val_norm <- val_decoded[[1]]
  colnames(val_norm) <- useMarkers

  # Build discovery SCE
  sce_disc <- SingleCellExperiment(
    assays = list(norm = t(df_norm)),
    colData = df_sub[, c('sample_id', 'Gensini', 'mg_cell_type_distinct')]
  )
  reducedDims(sce_disc) <- list(VAE = train_decoded$encoded %>% as.matrix() %>% as.data.frame())

  # Cluster discovery
  cat("  Clustering...\n")
  sce_disc <- runFuseSOM(sce_disc, numClusters = nclust, assay = 'norm',
                         verbose = FALSE, clusterCol = 'clusters')

  # ARI: compare unsupervised clusters to manual gates
  ari <- adjustedRandIndex(sce_disc$clusters, sce_disc$mg_cell_type_distinct)
  cat(sprintf("  ARI = %.4f\n", ari))

  # Build validation SCE
  sce_val <- SingleCellExperiment(
    assays = list(norm = t(val_norm)),
    colData = df_val[, c('sample_id', 'Gensini_bin')]
  )
  reducedDims(sce_val) <- list(VAE = val_decoded$encoded %>% as.matrix() %>% as.data.frame())

  # LDA cell type classification
  cat("  Fitting LDA...\n")
  train_x <- reducedDim(sce_disc, type = "VAE") %>% as.data.frame()
  train_x$cellTypes <- as.factor(sce_disc$clusters)
  test_x <- reducedDim(sce_val, type = "VAE") %>% as.data.frame()
  lda_fit <- MASS::lda(cellTypes ~ ., data = train_x)
  test_clusters <- predict(lda_fit, test_x)$class
  sce_val$clusters <- as.factor(test_clusters)

  # Features - discovery
  prop_train <- computeFeatures(sce = sce_disc, featureType = 'prop',
                                cellTypeCol = 'clusters', sampleCol = 'sample_id',
                                logit = TRUE, useMarkers = useMarkers)
  means_train <- computeFeatures(sce = sce_disc, featureType = 'mean',
                                 cellTypeCol = 'clusters', sampleCol = 'sample_id',
                                 logit = TRUE, useMarkers = useMarkers)
  colnames(means_train) <- as.factor(colnames(means_train))

  row_names <- rownames(prop_train)
  condition <- factor(clinicaldata[match(row_names, clinicaldata$sample_id), "Gensini_bin"])
  y_train <- as.numeric(levels(condition))[condition]

  # Features - validation
  prop_test <- computeFeatures(sce = sce_val, featureType = 'prop',
                               cellTypeCol = 'clusters', sampleCol = 'sample_id',
                               logit = TRUE, useMarkers = useMarkers)
  means_test <- computeFeatures(sce = sce_val, featureType = 'mean',
                                cellTypeCol = 'clusters', sampleCol = 'sample_id',
                                logit = TRUE, useMarkers = useMarkers)
  colnames(prop_test) <- as.factor(colnames(prop_test))
  colnames(means_test) <- as.factor(colnames(means_test))

  row_names_test <- rownames(prop_test)
  condition_test <- factor(clinicaldata_3[match(row_names_test, clinicaldata_3$sample_id), "Gensini_bin"])
  y_test <- as.numeric(levels(condition_test))[condition_test]

  # Tree and groups
  trees <- generateTree(features = prop_train, method = "ward")
  tree <- trees$tree
  groups <- generateGroups(tree = tree, nClust = nclust,
                           proportions = prop_train, means = means_train)
  groups <- lapply(groups, sort)

  # Combine and scale
  X_train <- cbind(prop_train, means_train)
  scaleVals <- preProcess(X_train, method = c('scale'))
  X_train <- predict(scaleVals, X_train) %>% as.matrix()
  X_test <- cbind(prop_test, means_test) %>% as.matrix()
  X_test <- predict(scaleVals, X_test)

  # Fit and predict
  cat("  Fitting model...\n")
  fit <- fitModel(xTrain = X_train, yTrain = y_train,
                  groups = groups, penalty = 'grLasso', seed = 1994,
                  modelSummary = FALSE)
  test_auc <- plotAUC(fit = fit$fit, xTest = X_test, yTest = y_test, title = "")
  auc_val <- as.numeric(auc(test_auc$preds$Truth, test_auc$preds$Predicted))

  cat(sprintf("  AUC = %.4f\n", auc_val))

  results <- rbind(results, data.frame(cells_per_sample = n_cells, ari = ari, auc = auc_val))

  keras3::clear_session()
  gc()
}

# --- Summary ---
cat("\n--- Downsampling Sensitivity Results ---\n")
print(results)

write.csv(results, file.path(revision_root, 'results/downsampling/bioheart_downsampling_results.csv'), row.names = FALSE)

# --- Plot ---
results$cells_label <- paste0(results$cells_per_sample / 1000, 'K')
results$cells_label <- factor(results$cells_label, levels = paste0(cell_levels / 1000, 'K'))

# Highlight the 10K bar (manuscript default)
results$highlight <- if_else(results$cells_per_sample == 10000, 'Manuscript default', 'Other')

p <- ggplot(results, aes(x = cells_label, y = auc, fill = highlight)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.3f", auc)), vjust = -0.5, size = 4) +
  geom_hline(yintercept = 0.5, linetype = 'dashed', colour = 'grey60') +
  scale_fill_manual(values = c('Manuscript default' = '#8B0000', 'Other' = '#B0B0B0')) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1.0, 0.2)) +
  labs(x = 'Cells per sample', y = 'AUC',
       title = 'Effect of downsampling on predictive performance') +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = 'bold'),
        legend.position = 'none')

ggsave(file.path(revision_root, 'figures/supp_downsampling_sensitivity.pdf'),
       p, width = 6, height = 5)
cat("\nFigure saved to figures/supp_downsampling_sensitivity.pdf\n")
