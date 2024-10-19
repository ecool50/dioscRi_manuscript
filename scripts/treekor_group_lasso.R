# Run treekor
exprs <- t(assay(sce_norm, "counts"))
clusters <- colData(sce_norm)$clusters
classes <- colData(sce_norm)$gensini_bin
samples <- colData(sce_norm)$sample_id


clust_tree <- getClusterTree(exprs,
                             clusters,
                             hierarchy_method="average", 
                             scale_exprs = T)

hc <- clust_tree$clust_tree

res <- MLGL:::preliminaryStep(hc)
var <- res$var
group <- res$group
x <- data_logit[, var] %>%
  as.matrix()
x_test <- data_3_logit[, var] %>%
  as.matrix()

fit <- cv.grpreg(X = x, y  = y, family = "binomial", type.meaure = 'auc',
                 penalty = 'grLasso', alpha = 0.01,
                 group = group, 
                 nfolds = 5, seed = 1994)
plot(fit)

predict(fit, type="group", lambda=fit$lambda.min)

pred_test <- predict(fit, X = x_test,
                     type = "response", which = 'lambda.min')

pred_fs <- prediction(pred_test, condition_3)
perf_fs <- ROCR::performance(pred_fs, "tpr", "fpr")


plt_dat = data.frame(
  FPR = perf_fs@x.values[[1]],
  TPR = perf_fs@y.values[[1]]
)

auc_perf <- ROCR::performance(pred_fs, measure = "auc")
auc <- auc_perf@y.values[[1]]

# generate AUC plot
p_auc <- ggplot(plt_dat, aes(x = FPR, y = TPR)) +
  geom_line(colour = "blue") +
  labs(x = perf_fs@x.name, y = perf_fs@y.name) +
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
  )  + ggtitle(paste("Test Set - AUC =", round(auc, 2)))

p_auc


coefs <- coef(fit, s = 'lambda.min') %>%
  as.matrix() %>%
  as.data.frame()
coefs$feature <- rownames(coefs)

coefs <- coefs[coefs$V1 != 0, ]

MP_gLasso(cv_object = fit, group = group, lambda.type = "min", sort.type = "max")
