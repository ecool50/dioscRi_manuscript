exprs <- t(assay(sce_norm, "counts"))
clusters <- colData(sce_norm)$clusters
samples <- colData(sce_norm)$sample_id
classes <- colData(sce_norm)$gensini_bin

clust_tree <- getClusterTreeCustom(exprs,
                                   clusters,
                                   hierarchy_method="ward", 
                                   scale_exprs = T)

tested_tree <- testTree(clust_tree$clust_tree,
                        clusters=clusters,
                        samples=samples,
                        classes=classes)

hc <- clust_tree$clust_tree

res <- getHCGroups(hc)
var <- res$varGroup
group <- res$indexGroup

var2 <- var + max(var)
group_2 <- group + max(group)


d_mean <- dist(scale(t(t_mean)))
hc_mean <- hclust(d_mean, method = 'ward')
res_mean <- getHCGroups(hc_mean)
var_mean <- res_mean$varGroup
group_mean <- res_mean$indexGroup
var_mean <- var_mean + max(var2)
group_mean <- group_mean + max(group_2)


var_final <- c(var, var2, var_mean)
group_final <- c(group, group_2, group_mean)

groups <- findSequenceIndices(group_final)
var_groups<- list()
for(i in 1:length(groups)){
  var_groups[[i]] <- var_final[groups[[i]]]
}

X_train <- cbind(data, data_logit, t_mean) %>% 
  scale()

y_train <- if_else(condition == 1, 1, 0) 

set.seed(999)
cvfit <- cv.grpregOverlap(X_train, y_train, var_groups, penalty = 'grLasso', 
                          family = 'binomial',  alpha = 0.05)
# plot(cvfit)

fit <- grpregOverlap(X_train, y, var_groups, nfolds = 10,
                     family='binomial', alpha = 0.05,
                     returnX.latent = T, returnOverlap = FALSE, 
                     lambda = cvfit$lambda.min, penalty = 'grLasso')

predict(fit, type="group")

X_test_study_5 <- cbind(data_5, data_5_logit, t_5_mean)
X_test_study_5 <- scale(X_test_study_5, center=attr(X_train, "scaled:center"), 
                scale=attr(X_train, "scaled:scale"))


X_test_study_3 <- cbind(data_3, data_3_logit, t_3_mean)
X_test_study_3 <- scale(X_test_study_3, center=attr(X_train, "scaled:center"), 
                        scale=attr(X_train, "scaled:scale"))

pred_study_3 <- predict(fit, X = X_test_study_3,
                        type = "response")

pred_study_3 <- prediction(pred_study_3, condition_3)
pred_study_3_perf <- ROCR::performance(pred_study_3, "tpr", "fpr")


plt_dat_3 = data.frame(
  FPR = pred_study_3_perf@x.values[[1]],
  TPR = pred_study_3_perf@y.values[[1]]
)

auc_perf_study_3 <- ROCR::performance(pred_study_3, measure = "auc")
auc_study_3 <- auc_perf_study_3@y.values[[1]]

# generate AUC plot
p_auc_study_3 <- ggplot(plt_dat_3, aes(x = FPR, y = TPR)) +
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
  )  + ggtitle(paste("Overlap Group - Study 3 - AUC =", round(auc_study_3, 2)))


pred_study_5 <- predict(fit, X = X_test_study_5,
                        type = "response")

pred_study_5 <- prediction(pred_study_5, condition_5)
pred_study_5_perf <- ROCR::performance(pred_study_5, "tpr", "fpr")


plt_dat_5 = data.frame(
  FPR = pred_study_5_perf@x.values[[1]],
  TPR = pred_study_5_perf@y.values[[1]]
)

auc_perf_study_5 <- ROCR::performance(pred_study_5, measure = "auc")
auc_study_5 <- auc_perf_study_5@y.values[[1]]

# generate AUC plot
p_auc_study_5 <- ggplot(plt_dat_5, aes(x = FPR, y = TPR)) +
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
  )  + ggtitle(paste("Overlap Group - Study 5 - AUC =", round(auc_study_5, 2)))


cowplot::plot_grid(p_auc_study_3, p_auc_study_5)

coefs <- coef(cvfit, s = 'lambda.min') %>%
  as.matrix() %>%
  as.data.frame()
coefs$feature <- rownames(coefs)
coefs <- coefs[coefs$V1 != 0, ]


res_mp <- MP_gLasso(fit = fit, group = fit$grp.vec, lambda.type = "min", sort.type = "max")
girafe(code = print(res_mp$plot), width_svg = 8, height_svg = 8)

