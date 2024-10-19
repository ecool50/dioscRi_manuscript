library(pvclust)
set.seed(1994)
dat <- cbind(data, data_logit)
dat_test <- cbind(data_3, data_3_logit)

res.tree <- .ConsensusTest(dat, n_min = 10)
clusters <- .MultiCutTRee(res.tree$Tree)

# data_temp <- data_ratio[, temp_clusters$clusters[[2]]]

meta_features <- data.frame(matrix(nrow = nrow(data), ncol = length(unique(temp_clusters))))

for(i in 1:length(unique(clusters))){
  cur_clusters <- clusters[[i]]
  meta_features[, i] <- rowMeans(dat[, cur_clusters, drop = FALSE])
  # print(val)
}
meta_features_test <- data.frame(matrix(nrow = nrow(dat_test), ncol = length(unique(temp_clusters))))

for(i in 1:length(unique(clusters))){
  cur_clusters <- clusters[[i]]
  meta_features_test[, i] <- rowMeans(dat_test[, cur_clusters, drop = FALSE])
  # print(val)
}

h_res <- runHOPACH(t(meta_features))
#
res_hopach <- computeHOPACHClusters(h_res = h_res, dat = meta_features, group = TRUE)
meta_features <- res_hopach$meta_features %>%
  as.matrix()
groups <- findSequenceIndices(res_hopach$groups)

meta_features_test <- computeHOPACHClusters(h_res = h_res, dat = meta_features_test, group = FALSE) %>%
  as.matrix()

x <- as.matrix(meta_features) # Removes class
y <- condition # Only class

## Ridge Regression to create the Adaptive Weights Vector
set.seed(999)
cv.ridge <- cv.glmnet(x, y, family='binomial', alpha=0, parallel=TRUE, standardize=TRUE)
w3 <- 1/abs(matrix(coef(cv.ridge, s=cv.ridge$lambda.min)
                   [, 1][2:(ncol(x)+1)] ))^2 ## Using gamma = 1
w3[w3[,1] == Inf] <- 999999999 ## Replacing values estimated as Infinite for 999999999

## Adaptive Lasso
set.seed(999)
cv.lasso <- cv.glmnet(x, y, family='binomial', alpha=1, parallel=TRUE, standardize=TRUE, 
                      type.measure='auc', penalty.factor=w3)
plot(cv.lasso)
# plot(cv.lasso$glmnet.fit, xvar="lambda", label=TRUE)
# abline(v = log(cv.lasso$lambda.min))
# abline(v = log(cv.lasso$lambda.1se))
coef(cv.lasso, s='lambda.min')
coef <- coef(cv.lasso, s='lambda.min')
selected_attributes <- (coef@i[-1]+1) ## Considering the structure of the data frame dataF as shown earlier


pred_test <- predict(cv.lasso, newx = as.matrix(meta_features_test), 
                     type = "response", s = 'lambda.min')

pred_fs <- prediction(pred_test[, 1], condition_3)
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

# Combine prediction and truth
result_test <- data.frame("Truth"=condition_3,
                          "Prob"=pred_test[, 1])
## ---- results='asis',message=FALSE,fig.width = 5-------------------------
# Use decision tree to find the cell subsets that are associated the AML


x <- as.data.frame(meta_features_test[, selected_attributes, drop=F])
P = pred_test[, 1]
colnames(x) <- gsub("[^[:alnum:]]", "_", colnames(x))
x <- cbind.data.frame(P = P, x)

fit <- rpart(P ~ ., method = "anova", data = x, control = list(minsplit = 5))
rpart.plot(fit)