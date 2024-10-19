  library(grpreg)
  library(rpart)
  library(rpart.plot)

  dat <- cbind(data, data_logit)
  dat_test <- cbind(data_3, data_3_logit)
  colnames(dat) <- janitor::make_clean_names(colnames(dat))
  
  h_res <- runHOPACHCustom(data = t(dat), dissimilarity_metric = 'cor')
  
  
  res <- computeHOPACHClusters(h_res, dat = dat, group = TRUE)
  dat <- res$meta_features
  groups <- res$groups
  
  dat_test <- computeHOPACHClusters(h_res = h_res, dat = dat_test, group = F)
  
  
  fit.cv <- cv.grpreg(as.matrix(dat), as.numeric(condition), 
                      group = groups, 
                penalty="grLasso", family = "binomial", alpha = 1, seed = 1995,
                nfolds = 5, trace = TRUE, returnY = TRUE)
  
  # coef(fit.cv)
  plot(fit.cv)
  sig_vars <- predict(fit.cv, type="vars", lambda=fit.cv$lambda.min) %>%
    names()
  
  
  preds <- predict(fit.cv, X=as.matrix(res$meta_features), type = 'response')
  
  cal_1 <- beta_calibration(preds, condition, parameters="am")
  preds <- beta_predict(preds, cal_1)
  
  pred_fs <- prediction(preds, condition)
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
    )  + ggtitle(paste("Training Set - AUC =", round(auc, 2)))
  
  p_auc


preds_test <- predict(fit.cv, X=as.matrix(dat_test), type = 'response')
cal_2 <- beta_calibration(preds_test, condition_3, parameters="am")
preds_test <- beta_predict(preds_test, cal_2)

pred_fs <- prediction(preds_test, condition_3)
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
                          "Prob"=preds_test)
## ---- results='asis',message=FALSE,fig.width = 5-------------------------
# Use decision tree to find the cell subsets that are associated the AML


x <- as.data.frame(dat_test[, sig_vars, drop=F])
P = val1_score[, 1]
colnames(x) <- gsub("[^[:alnum:]]", "_", colnames(x))
x <- cbind.data.frame(P = P, x)

fit <- rpart(P ~ ., method = "anova", data = x, control = list(minsplit = 5))
rpart.plot(fit)
 