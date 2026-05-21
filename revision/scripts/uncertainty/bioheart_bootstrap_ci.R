# BioHEART-CT - Bootstrap CIs on Validation Cohort
# The train/test split is fixed (discovery=train, validation=test).
# We run the full pipeline once, then bootstrap the validation predictions
# to get 95% CIs on the AUC.
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
tensorflow::set_random_seed(1994)

# --- Config ---
nBootstrap <- 2000
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

# --- Load normalised data ---
cat("Loading data...\n")
load(file.path(base_dir, 'data/sce_dat/study_4_MMD_VAE_updated.RData'))  # sce_norm (discovery)
load(file.path(base_dir, 'data/sce_dat/study_3_MMD_VAE_updated.RData'))  # sce_3_norm (validation)

# --- Cluster discovery data ---
cat("Clustering discovery cohort...\n")
sce_norm$cluster_col <- NULL
sce_norm <- runFuseSOM(sce_norm, numClusters = nclust, assay = 'norm',
                       verbose = FALSE, clusterCol = 'clusters_norm')

# --- LDA cell type classification ---
cat("Fitting LDA...\n")
train_x <- reducedDim(sce_norm, type = "VAE") %>% as.data.frame()
train_x$cellTypes <- as.factor(sce_norm$clusters_11)
test_x <- reducedDim(sce_3_norm, type = "VAE") %>% as.data.frame()
lda_fit <- MASS::lda(cellTypes ~ ., data = train_x)
df_3_clusters <- predict(lda_fit, test_x)$class
sce_3_norm$clusters <- as.factor(df_3_clusters)

# --- Clinical data ---
clinicaldata <- colData(sce_norm) %>%
  as.data.frame() %>%
  dplyr::select(sample_id, Gensini, Gender, Age) %>%
  distinct()
clinicaldata$Gensini_bin <- factor(if_else(clinicaldata$Gensini > 0, 1, 0))

clinicaldata_3 <- colData(sce_3_norm) %>%
  as.data.frame() %>%
  dplyr::select(sample_id, Gensini_bin, Gender, Age) %>%
  distinct()
clinicaldata_3$Gensini_bin <- as.factor(clinicaldata_3$Gensini_bin)

# --- Features: discovery ---
cat("Computing features...\n")
data_logit <- computeFeatures(sce = sce_norm, featureType = 'prop',
                              cellTypeCol = cluster_col, sampleCol = 'sample_id',
                              logit = TRUE, useMarkers = useMarkers)
row_names <- rownames(data_logit)
condition <- factor(clinicaldata[match(row_names, clinicaldata$sample_id), "Gensini_bin"])

markerMeanCellType <- computeFeatures(sce = sce_norm, featureType = 'mean',
                              cellTypeCol = cluster_col, sampleCol = 'sample_id',
                              logit = TRUE, useMarkers = useMarkers)
colnames(markerMeanCellType) <- as.factor(colnames(markerMeanCellType))

# --- Features: validation ---
data_3_logit <- computeFeatures(sce = sce_3_norm, featureType = 'prop',
                              cellTypeCol = 'clusters', sampleCol = 'sample_id',
                              logit = TRUE, useMarkers = useMarkers)
row_names_3 <- rownames(data_3_logit)
condition_3 <- factor(clinicaldata_3[match(row_names_3, clinicaldata_3$sample_id), "Gensini_bin"])
colnames(data_3_logit) <- as.factor(colnames(data_3_logit))

markerMeanCellType_3 <- computeFeatures(sce = sce_3_norm, featureType = 'mean',
                              cellTypeCol = 'clusters', sampleCol = 'sample_id',
                              logit = TRUE, useMarkers = useMarkers)
colnames(markerMeanCellType_3) <- as.factor(colnames(markerMeanCellType_3))

# --- Tree and groups ---
nclust <- length(unique(sce_norm$clusters_norm))
trees <- generateTree(features = data_logit, method = "ward")
tree <- trees$tree
groups <- generateGroups(tree = tree, nClust = nclust,
                         proportions = data_logit, means = markerMeanCellType)

# --- Fit model ---
cat("Fitting model...\n")
X_train <- cbind(data_logit, markerMeanCellType)
scaleVals <- preProcess(X_train, method = c('scale'))
X_train <- predict(scaleVals, X_train) %>% as.matrix()
y_train <- as.numeric(levels(condition))[condition]

X_test <- cbind(data_3_logit, markerMeanCellType_3) %>% as.matrix()
X_test <- predict(scaleVals, X_test)
y_test <- as.numeric(levels(condition_3))[condition_3]

groups <- lapply(groups, sort)
fit <- fitModel(xTrain = X_train, yTrain = y_train, modelSummary = FALSE,
                groups = groups, penalty = "grLasso", seed = 1994)

# --- Get predictions ---
test_auc_res <- plotAUC(fit = fit$fit, xTest = X_test, yTest = y_test, title = "")
point_auc <- auc(test_auc_res$preds$Truth, test_auc_res$preds$Predicted)
cat(sprintf("Point AUC: %.4f\n", point_auc))

# --- Bootstrap CIs via pROC (stratified, 2000 resamples) ---
# Stratified bootstrap preserves case/control ratio per resample (Robin et al., 2011)
# 2000 resamples recommended for CIs (Carpenter & Bithell, 2000)
cat(sprintf("Running stratified bootstrap (%d resamples)...\n", nBootstrap))
preds_df <- test_auc_res$preds

set.seed(1994)
roc_obj <- roc(preds_df$Truth, preds_df$Predicted, quiet = TRUE)
ci_result <- ci.auc(roc_obj, method = "bootstrap", boot.n = nBootstrap,
                    boot.stratified = TRUE, conf.level = 0.95)

# Also extract the full bootstrap distribution for plotting
set.seed(1994)
boot_aucs <- replicate(nBootstrap, {
  pos_idx <- which(preds_df$Truth == 1)
  neg_idx <- which(preds_df$Truth == 0)
  boot_pos <- sample(pos_idx, length(pos_idx), replace = TRUE)
  boot_neg <- sample(neg_idx, length(neg_idx), replace = TRUE)
  boot_preds <- preds_df[c(boot_pos, boot_neg), ]
  as.numeric(auc(boot_preds$Truth, boot_preds$Predicted, quiet = TRUE))
})

# --- Summary ---
cat(sprintf("\n--- BioHEART-CT Bootstrap Results (n=%d) ---\n", nBootstrap))
cat(sprintf("Point AUC: %.4f\n", point_auc))
cat(sprintf("95%% CI (pROC stratified bootstrap): [%.4f, %.4f]\n",
            ci_result[1], ci_result[3]))
cat(sprintf("Bootstrap mean: %.4f, SD: %.4f\n", mean(boot_aucs), sd(boot_aucs)))

# Save
results <- data.frame(
  point_auc = as.numeric(point_auc),
  mean_auc = mean(boot_aucs),
  sd_auc = sd(boot_aucs),
  ci_lower = as.numeric(ci_result[1]),
  ci_upper = as.numeric(ci_result[3])
)
write.csv(results, file.path(revision_root, 'results/uncertainty/bioheart_bootstrap_ci.csv'), row.names = FALSE)

# Save full bootstrap distribution for plotting
boot_df <- data.frame(bootstrap_id = 1:length(boot_aucs), auc = boot_aucs)
# write.csv(boot_df, file.path(revision_root, 'results/uncertainty/bioheart_bootstrap_distribution.csv'), row.names = FALSE)

cat("Results saved.\n")
