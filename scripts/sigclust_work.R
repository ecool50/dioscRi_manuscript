library(multiview)

dat <- data_ratio

# temp <- .ConsensusTest(t(dat), n_min = 10)
# .PlotRes(temp$Tree, datName = "", use_labs = T)
# 
clusters <- varClusterMod(dat, clus.summary = "max.pw",
                         corr.method = "pearson", corr.use = "complete.obs",
                         clus.method = "complete", percentile = 0.5)



nClust <- length(unique(clusters))
meta_features <- data.frame(matrix(nrow = nrow(dat), ncol = nClust))

count <- 1
for (i in unique(clusters)) {
  cluster_features <- which(clusters == i)
  # print(cluster_features)
  meta_features[, count] <- rowMeans(dat[, cluster_features, drop = FALSE])
  count <- count + 1
}



result <- ClassifyR::crossValidate(meta_features, condition, classifier = c("elasticNetGLM"),
                                   nFeatures = NULL, nFolds = 5, nRepeats = 10, nCores = 4,
                                   performanceType = 'AUC', selectionOptimisation = "none")
ROCplot(result)


cv_dat_glm <- KFoldCustom(train_x = meta_features, train_y = condition, nRepeats = 10,
                          method = 'glm', alpha = 1)
cv_dat_svm <- KFoldCustom(train_x = meta_features, train_y = condition, nRepeats = 10,
                          method = 'svm')



pred_5_fs <- prediction(preds, condition)
perf_5_fs <- ROCR::performance(pred_5_fs, "tpr", "fpr")


plt_dat_5 = data.frame(
  FPR = perf_5_fs@x.values[[1]],
  TPR = perf_5_fs@y.values[[1]]
)

auc_perf <- ROCR::performance(pred_5_fs, measure = "auc")
auc <- auc_perf@y.values[[1]]

# generate AUC plot
ggplot(plt_dat_5, aes(x = FPR, y = TPR)) +
  geom_line(colour = "blue") +
  labs(x = perf_5_fs@x.name, y = perf_5_fs@y.name) +
  geom_abline(slope = 1, intercept = 0) + theme_bw() +
  theme(
    plot.title = element_text(color="Black", size=16, face="bold", hjust = 0.5),
    plot.subtitle = element_text(color = "red", size = 16, hjust = 0.5),
    axis.title.x = element_text(color="Black", size=16),
    axis.title.y = element_text(color="Black", size=16),
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 16),
    strip.text.x = element_text(size = 16),
    legend.title = element_text(size=16), #change legend title font size
    legend.text = element_text(size=16) #change legend text font size)
  )  + ggtitle(paste("Study 5 - FuseSOM AUC =", round(auc, 2)))
