dat <- cbind(data, data_logit)
colnames(dat) <- janitor::make_clean_names(colnames(dat))
dat <- as.matrix(dat)
size <- computeGridSize(dat)
dat_test <- cbind(data_3, data_3_logit)
colnames(dat_test) <- janitor::make_clean_names(colnames(dat_test))
dat_test <- as.matrix(dat_test)

h_res <- runHOPACHCustom(t(dat), kmax = size, K = size)

res_hopach <- computeHOPACHClusters(h_res = h_res, dat = dat, group = TRUE)
dat <- res_hopach$meta_features %>%
  as.matrix()
groups <- findSequenceIndices(res_hopach$groups)

dat_test <- computeHOPACHClusters(h_res = h_res, dat = dat_test, group = FALSE) %>%
  as.matrix()

# groups <- rep(1:4, c(ncol(data), ncol(data_logit), ncol(data_ratio), ncol(t)))
# groups <- findSequenceIndices(groups)


set.seed(1234)
pl_fit1 <- prioritylasso(X = dat, Y = condition, family = "binomial", type.measure = "auc", 
                    blocks = groups, nfolds = 10,
                     standardize = FALSE, lambda.type = "lambda.min", block1.penalization = T)
pl_fit1$min.cvm


coeff1 <- pl_fit1$coefficients
coeff1 <- coeff1[coeff1 != 0] %>%
  na.omit()
print(round(coeff1, 4))

library(pROC)

pl1_score <- dat[ , names(coeff1), drop=F] %*% coeff1
temp <- predict(pl_fit1, newdata = dat, type = 'response')
pl1_roc <- roc(factor(condition), predict(pl_fit1, newdata = dat, type = 'response')[,1])

auc(pl1_roc)

# plot.roc(pl1_roc, grid=0.1)

# dat_test <- cbind(data_3, data_3_logit, data_ratio_3, t_3) %>%
#   as.matrix()

# colnames(dat_test) <- janitor::make_clean_names(colnames(dat_test))

val1_score <- dat_test[ , names(coeff1), drop=F] %*% coeff1
val1_roc <- roc(factor(condition_3), c(val1_score))
auc(val1_roc)

# plot.roc(pl1_roc, grid=0.1, print.auc = TRUE)
plot.roc(val1_roc, col="red")

legend("topleft", legend=paste0("Test set AUC: ", round(auc(val1_roc), 2)),
       col="red", lwd=2)


# Combine prediction and truth
result_test <- data.frame("Truth"=condition_3,
                          "Prob"=val1_score[, 1])

# Plot the prediction
stripchart(result_test$Prob~result_test$Truth, jitter = 0.1,
           vertical = TRUE, method = "jitter", pch = 20,
           xlab="Truth",ylab="Predicted Prob of Gensini")

## ---- results='asis',message=FALSE,fig.width = 5-------------------------
# Use decision tree to find the cell subsets that are associated the AML.
TG_test <- treeGate(P = val1_score[, 1],
                    x= as.data.frame(dat_test[, names(coeff1), drop=F ]))

