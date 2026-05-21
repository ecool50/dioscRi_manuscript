# CMV SDY519 - Bootstrap CIs on Fixed Test Set
# The split follows the DeepCNN paper: SDY519 held out for testing.
# We run the full pipeline once, then bootstrap the test predictions.

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

# --- Config ---
nBootstrap <- 2000
nclust <- 11
bioheart_root <- Sys.getenv("BIOHEART_ROOT", unset = "~/Documents/Academic/PhD/bioheart_analysis")
revision_root <- here::here("revision")
base_dir <- bioheart_root  # alias: upstream inputs (raw data, sce objects, etc.)
sce_dir <- file.path(base_dir, 'data/dioscRi_analysis_data/sce_objects/CMV_Study_SDY519')

split_and_extract <- function(string, n) {
  strsplit(string, "\\.")[[1]][[2]]
}

useMarkers <- c("TCRGD","IGD","HLADR","CD94","CD85J","CD8","CD56",
                "CD45RA","CD4","CD38","CD33","CD3","CD28","CD27",
                "CD25","CD24","CD20","CD19","CD161","CD16","CD14","CD127","CCR7")

# --- Load normalised data ---
cat("Loading data...\n")
load(file.path(sce_dir, 'cmv_train_MMD_VAE_sce.RData'))  # sce_train
load(file.path(sce_dir, 'cmv_test_MMD_VAE_sce.RData'))    # sce_test

# --- Cluster training data ---
cat("Clustering...\n")
sce_train$clusters <- NULL
sce_test$clusters <- NULL
sce_train <- runFuseSOM(sce_train, numClusters = nclust, assay = 'norm',
                        verbose = FALSE)

# --- LDA classifier ---
cat("Fitting LDA...\n")
train_x <- reducedDim(sce_train, type = "VAE") %>% as.data.frame()
train_x$cellTypes <- as.factor(sce_train$clusters)
test_x <- reducedDim(sce_test, type = "VAE") %>% as.data.frame()
lda_fit <- MASS::lda(cellTypes ~ ., data = train_x)
test_clusters <- predict(lda_fit, test_x)$class
sce_test$clusters <- as.factor(test_clusters)

# --- Clinical data ---
clinicaldata <- colData(sce_train) %>%
  as.data.frame() %>%
  dplyr::select(sample_id, CMV_Ab) %>%
  distinct()
clinicaldata$condition <- factor(if_else(clinicaldata$CMV_Ab == "True", 1, 0))

clinicaldata_test <- colData(sce_test) %>%
  as.data.frame() %>%
  dplyr::select(sample_id, CMV_Ab) %>%
  distinct()
clinicaldata_test$condition <- factor(if_else(clinicaldata_test$CMV_Ab == "True", 1, 0))

# --- Features: training ---
cat("Computing features...\n")
data_logit <- computeFeatures(sce = sce_train, featureType = 'prop',
                              cellTypeCol = 'clusters', sampleCol = 'sample_id',
                              logit = TRUE, useMarkers = useMarkers)
row_names <- rownames(data_logit)
condition <- factor(clinicaldata[match(row_names, clinicaldata$sample_id), "condition"])

markerMeanCellType <- computeFeatures(sce = sce_train, featureType = 'mean',
                              cellTypeCol = 'clusters', sampleCol = 'sample_id',
                              logit = TRUE, useMarkers = useMarkers)

# --- Features: testing ---
data_test_logit <- computeFeatures(sce = sce_test, featureType = 'prop',
                              cellTypeCol = 'clusters', sampleCol = 'sample_id',
                              logit = TRUE, useMarkers = useMarkers)
row_names_test <- rownames(data_test_logit)
condition_test <- factor(clinicaldata_test[match(row_names_test, clinicaldata_test$sample_id),
                                           "condition"])

markerMeanCellType_test <- computeFeatures(sce = sce_test, featureType = 'mean',
                              cellTypeCol = 'clusters', sampleCol = 'sample_id',
                              logit = TRUE, useMarkers = useMarkers)

# --- Tree and groups ---
nclust <- length(unique(test_clusters))
tree <- generateTree(features = data_logit, method = "ward")
tree <- tree$tree
groups <- generateGroups(tree = tree, nClust = nclust,
                         proportions = data_logit, means = markerMeanCellType)

# --- Fit model ---
cat("Fitting model...\n")
X_train <- cbind(data_logit, markerMeanCellType)
scaleVals <- preProcess(X_train, method = c('scale'))
X_train <- predict(scaleVals, X_train) %>% as.matrix()
y_train <- as.numeric(levels(condition))[condition]

X_test <- cbind(data_test_logit, markerMeanCellType_test) %>% as.matrix()
X_test <- predict(scaleVals, X_test)
y_test <- as.numeric(levels(condition_test))[condition_test]

fit <- fitModel(xTrain = X_train, yTrain = y_train, groups = groups,
                penalty = "grLasso", seed = 1994, modelSummary = FALSE)

# --- Get predictions ---
test_auc_res <- plotAUC(fit = fit$fit, xTest = X_test, yTest = y_test, title = "")
point_auc <- auc(test_auc_res$preds$Truth, test_auc_res$preds$Predicted)
cat(sprintf("Point AUC: %.4f\n", point_auc))

# --- Bootstrap CIs via pROC (stratified, 2000 resamples) ---
cat(sprintf("Running stratified bootstrap (%d resamples)...\n", nBootstrap))
preds_df <- test_auc_res$preds

set.seed(1994)
roc_obj <- roc(preds_df$Truth, preds_df$Predicted, quiet = TRUE)
ci_result <- ci.auc(roc_obj, method = "bootstrap", boot.n = nBootstrap,
                    boot.stratified = TRUE, conf.level = 0.95)

# Full bootstrap distribution for plotting
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
cat(sprintf("\n--- CMV SDY519 Bootstrap Results (n=%d) ---\n", nBootstrap))
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
write.csv(results, file.path(revision_root, 'results/uncertainty/cmv_bootstrap_ci.csv'), row.names = FALSE)

boot_df <- data.frame(bootstrap_id = 1:length(boot_aucs), auc = boot_aucs)
write.csv(boot_df, file.path(revision_root, 'results/uncertainty/cmv_bootstrap_distribution.csv'), row.names = FALSE)

cat("Results saved.\n")
