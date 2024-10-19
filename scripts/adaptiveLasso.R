exprs <- t(assay(sce_norm, "counts"))
clusters <- colData(sce_norm)$clusters
samples <- colData(sce_norm)$sample_id


clust_tree <- getClusterTreeCustom(exprs,
                                   clusters,
                                   hierarchy_method="average", 
                                   scale_exprs = T)

hc <- clust_tree$clust_tree

res <- MLGL:::preliminaryStep(hc)
var <- res$var
group <- res$group


## Proportions
x_1 <- data[, var] %>%
  as.matrix() 
x_1_test <- data_3[, var] %>%
  as.matrix()

## Logit Proportions
group_2 <- group + max(group)
x_2 <- data_logit[, var] %>%
  as.matrix()
x_2_test <- data_3_logit[, var] %>%
  as.matrix()

# Marker mean
d <- dist(scale(t(t_mean)))
hc3 <- hclust(d, method = 'average')
res3 <- MLGL:::preliminaryStep(hc3)
var_3 <- res3$var
group_3 <- res3$group
group_3 <- group_3 + max(group_2)
x_3 <- t_mean[, var_3] %>%
  as.matrix() 
x_3_test <- t_3_mean[, var_3] %>%
  as.matrix()

# consolidate the groups
group_final <- c(group, group_2, group_3)

# setup training and testing data
x <- cbind(x_1, x_2, x_3)
x_test <- cbind(x_1_test, x_2_test, x_3_test) 

meta_features <- computeClusters(groups = groups, dat = x) %>%
  as.matrix()

meta_features_test <- computeClusters(groups = groups, dat = x_test) %>%
  as.matrix()

y <- if_else(condition == 1, 1, 0)
## Ridge Regression to create the Adaptive Weights Vector
set.seed(999)
cv.ridge <- cv.glmnet(meta_features, y, family='binomial', alpha=0, parallel=TRUE, 
                      standardize=TRUE, relax = F)

best_ridge_coef <- as.numeric(coef(cv.ridge, s = cv.ridge$lambda.min))[-1]
## Adaptive Lasso
set.seed(999)
cv.lasso <- cv.glmnet(meta_features, y, family='binomial', alpha=0.5, parallel=TRUE, 
                      standardize=TRUE, relax = T, type.measure = 'auc',
                      penalty.factor=1/(abs(best_ridge_coef)))
plot(cv.lasso)
coef(cv.lasso, s='lambda.min')
coef <- coef(cv.lasso, s='lambda.min')
selected_attributes <- (coef@i[-1]+1)
pred_test <- predict(cv.lasso, newx = meta_features_test, 
                     type = "response", s = 'lambda.min')

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


# Extract the coefficients from all three models
# coefs_adl <- coef(cv.lasso, s = 'lambda.min') %>%
#   as.matrix() %>%
#   as.data.frame()

# Combine the coefficients
# coefs_adl <- data.frame(rownames(coefs_adl), coefs_adl)
# colnames(coefs_adl) <- c('Feature', 'AdaptiveLasso')
# coefs_adl <- coefs_adl[coefs_adl$AdaptiveLasso != 0,]
# coefs_adl <- coefs_adl[-1, ] %>%
#   as.data.frame()

# Combine prediction and truth
# result_test <- data.frame("Truth"=condition_3,
#                           "Prob"=pred_test[, 1])
## ---- results='asis',message=FALSE,fig.width = 5-------------------------
# Use decision tree to find the cell subsets that are associated the AML
# 
# 
# X <- as.data.frame(dat_test[, coefs_adl$Feature, drop=F])
# P = pred_test[, 1]
# colnames(X) <- gsub("[^[:alnum:]]", "_", colnames(X))
# X <- cbind.data.frame(P = P, X)
# 
# fit <- rpart(P ~ ., method = "anova", data = X, control = list(minsplit = 5))
# rpart.plot(fit)


# Assume 'list_of_lists' is your original list of lists, where each sublist contains cluster labels for each data point at that level
# num_levels <- length(h_res$cutree_list)
# num_data_points <- length(h_res$cutree_list[[1]]) # Assuming all levels have the same number of data points

# Initialize an empty matrix
# cluster_matrix <- matrix(nrow = num_levels, ncol = num_data_points)

# Fill the matrix with cluster labels
# for (i in 1:num_levels) {
#   cluster_matrix[i, ] <- h_res$cutree_list[[i]]
# }

# Convert to a dataframe if necessary
# cluster_df <- as.data.frame(cluster_matrix)

# Optional: Naming the rows and columns
# rownames(cluster_df) <- paste0("Level_", 1:num_levels)
# colnames(cluster_df) <- names(h_res$cutree_list[[1]])

# cluster_df <- t(cluster_df) %>%
#   as.data.frame()
# cluster_df$label <- rownames(cluster_df)

# temp <- cluster_df[, c('Level_2','Level_3', 'label')]
# labs <- temp[temp$Level_3 == 1, 'label']


# new_dat <- data %>% 
#   dplyr::select(-c(labs))
# sin(new_dat/2)^2 %>%
#   as.matrix() %>%
#   boxplot()




