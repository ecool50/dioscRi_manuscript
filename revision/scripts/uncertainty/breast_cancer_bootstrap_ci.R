# Breast Cancer - Bootstrap CIs on Test Set
# Uses the pre-normalized SCEs (canonical train/test split + MMD-VAE).
# Parallels bioheart_bootstrap_ci.R for uniform reporting across all four datasets.
rm(list = ls())
suppressPackageStartupMessages({
  library(dioscRi)
  library(FuseSOM)
  library(SingleCellExperiment)
  library(tidyverse)
  library(data.tree)
  library(caret)
  library(pROC)
})
set.seed(1994)

# --- Config ---
nBootstrap <- 2000
nclust <- 11
bioheart_root <- Sys.getenv("BIOHEART_ROOT", unset = "~/Documents/Academic/PhD/bioheart_analysis")
revision_root <- here::here("revision")
base_dir <- bioheart_root  # alias: upstream inputs (raw data, sce objects, etc.)
data_dir <- '~/Documents/Academic/PhD/DeepLearning_CyTOF/Breast_Cancer_Wagner_2019/data'

useMarkers <- c("AKT","AR","BCL2","CA9","CD24","CD44","CD49f","cMET","cMYC",
                "CyclinB1","ECadherin","EGFR","EpCAM","ERa","EZH2","H3K27me3",
                "HER2","HLADR","K14","K5","K7","K8K18","Ki67","p53","panK",
                "PRB","PTEN","SMA","Survivin","Vimentin")

# --- Load normalised data ---
cat("Loading pre-normalized SCEs...\n")
load(file.path(data_dir, 'breast_cancer_train_MMD_VAE_sce.RData'))   # sce_norm
load(file.path(data_dir, 'breast_cancer_test_MMD_VAE_sce.RData'))    # sce_norm_test

sce_norm$sample_id <- sce_norm$file_name
sce_norm_test$sample_id <- sce_norm_test$file_name

# --- Cluster training data ---
cat("Clustering training cohort...\n")
sce_norm <- runFuseSOM(sce_norm, numClusters = nclust, assay = 'norm',
                       verbose = FALSE, clusterCol = 'clusters')

# --- LDA cell type classification ---
cat("Fitting LDA...\n")
train_x <- reducedDim(sce_norm, type = "VAE") %>% as.data.frame()
train_x$cellTypes <- as.factor(sce_norm$clusters)
test_x <- reducedDim(sce_norm_test, type = "VAE") %>% as.data.frame()
lda_fit <- MASS::lda(cellTypes ~ ., data = train_x)
test_clusters <- predict(lda_fit, test_x)$class
sce_norm_test$clusters <- as.factor(test_clusters)

# --- Clinical data ---
clinicaldata <- colData(sce_norm) %>% as.data.frame() %>%
  dplyr::select(sample_id, condition) %>% distinct()
clinicaldata$y <- factor(if_else(clinicaldata$condition == "T", 1, 0))

clinicaldata_test <- colData(sce_norm_test) %>% as.data.frame() %>%
  dplyr::select(sample_id, condition) %>% distinct()
clinicaldata_test$y <- factor(if_else(clinicaldata_test$condition == "T", 1, 0))

# --- Features ---
cat("Computing features...\n")
prop_train <- computeFeatures(sce = sce_norm, featureType = 'prop',
                              cellTypeCol = 'clusters', sampleCol = 'sample_id',
                              logit = TRUE, useMarkers = useMarkers)
means_train <- computeFeatures(sce = sce_norm, featureType = 'mean',
                               cellTypeCol = 'clusters', sampleCol = 'sample_id',
                               logit = TRUE, useMarkers = useMarkers)
colnames(means_train) <- as.factor(colnames(means_train))

prop_test <- computeFeatures(sce = sce_norm_test, featureType = 'prop',
                             cellTypeCol = 'clusters', sampleCol = 'sample_id',
                             logit = TRUE, useMarkers = useMarkers)
means_test <- computeFeatures(sce = sce_norm_test, featureType = 'mean',
                              cellTypeCol = 'clusters', sampleCol = 'sample_id',
                              logit = TRUE, useMarkers = useMarkers)
colnames(prop_test) <- as.factor(colnames(prop_test))
colnames(means_test) <- as.factor(colnames(means_test))

# --- Tree and groups ---
trees <- generateTree(features = prop_train, method = "ward")
tree <- trees$tree
groups <- generateGroups(tree = tree, nClust = nclust,
                         proportions = prop_train, means = means_train)
groups <- lapply(groups, sort)

# --- Fit model ---
cat("Fitting group lasso...\n")
row_train <- rownames(prop_train)
condition <- factor(clinicaldata[match(row_train, clinicaldata$sample_id), "y"])
y_train <- as.numeric(levels(condition))[condition]

row_test <- rownames(prop_test)
condition_test <- factor(clinicaldata_test[match(row_test, clinicaldata_test$sample_id), "y"])
y_test <- as.numeric(levels(condition_test))[condition_test]

X_train <- cbind(prop_train, means_train)
scaleVals <- preProcess(X_train, method = c('scale'))
X_train <- predict(scaleVals, X_train) %>% as.matrix()

X_test <- cbind(prop_test, means_test) %>% as.matrix()
X_test <- predict(scaleVals, X_test)

fit <- fitModel(xTrain = X_train, yTrain = y_train, modelSummary = FALSE,
                groups = groups, penalty = "grLasso", seed = 1994)

test_auc_res <- plotAUC(fit = fit$fit, xTest = X_test, yTest = y_test, title = "")
point_auc <- as.numeric(auc(test_auc_res$preds$Truth, test_auc_res$preds$Predicted))
cat(sprintf("Point AUC: %.4f\n", point_auc))

# --- Bootstrap CIs (stratified, 2000 resamples) ---
cat(sprintf("Running stratified bootstrap (%d resamples)...\n", nBootstrap))
preds_df <- test_auc_res$preds

set.seed(1994)
roc_obj <- roc(preds_df$Truth, preds_df$Predicted, quiet = TRUE)
ci_result <- ci.auc(roc_obj, method = "bootstrap", boot.n = nBootstrap,
                    boot.stratified = TRUE, conf.level = 0.95)

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
cat(sprintf("\n--- Breast Cancer Bootstrap Results (n=%d) ---\n", nBootstrap))
cat(sprintf("Point AUC: %.4f\n", point_auc))
cat(sprintf("95%% CI (pROC stratified bootstrap): [%.4f, %.4f]\n",
            ci_result[1], ci_result[3]))
cat(sprintf("Bootstrap mean: %.4f, SD: %.4f\n", mean(boot_aucs), sd(boot_aucs)))

# Save
results <- data.frame(
  point_auc = point_auc,
  mean_auc = mean(boot_aucs),
  sd_auc = sd(boot_aucs),
  ci_lower = as.numeric(ci_result[1]),
  ci_upper = as.numeric(ci_result[3])
)
write.csv(results, file.path(revision_root, 'results/uncertainty/breast_cancer_bootstrap_ci.csv'),
          row.names = FALSE)

boot_df <- data.frame(bootstrap_id = seq_along(boot_aucs), auc = boot_aucs)
write.csv(boot_df, file.path(revision_root, 'results/uncertainty/breast_cancer_bootstrap_distribution.csv'),
          row.names = FALSE)

cat("Saved.\n")
