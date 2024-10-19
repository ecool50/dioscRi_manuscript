train_x <- cbind(data, data_logit, t_mean) %>%
  scale()
test_x_study_3 <- cbind(data_3, data_3_logit, t_3_mean)
test_x_study_3 <- scale(test_x_study_3, center=attr(train_x, "scaled:center"), 
                        scale=attr(train_x, "scaled:scale"))

test_x_study_5 <- cbind(data_5, data_5_logit, t_5_mean)
test_x_study_5 <- scale(test_x_study_5, center=attr(train_x, "scaled:center"), 
                        scale=attr(train_x, "scaled:scale"))

train_y <- if_else(condition == 1, 1, 0) 

## Fit a base Logistic regression model
set.seed(999)
cv.ridge <- cv.glmnet(train_x, train_y, family='binomial', alpha=0, parallel=TRUE, 
                      standardize=F, relax = F)

best_ridge_coef <- as.numeric(coef(cv.ridge, s = cv.ridge$lambda.min))[-1]
## Adaptive Lasso
set.seed(999)
lr_base <- cv.glmnet(train_x, train_y, family='binomial', alpha=0.05, parallel=TRUE, 
                      standardize=F, relax = F, type.measure = 'auc',
                      penalty.factor=1/(abs(best_ridge_coef)))

plot(lr_base)

pred_study_3 <- predict(lr_base, newx = X_test_study_3, 
                        s = lr_base$lambda.min, type = "response")

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
  )  + ggtitle(paste("Adaptive Lasso - Study 3 - AUC =", round(auc_study_3, 2)))


pred_study_5 <- predict(lr_base, newx = X_test_study_5,
                        type = "response", s = lr_base$lambda.min)

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
  )  + ggtitle(paste("Adaptive Lasso - Study 5 - AUC =", round(auc_study_5, 2)))


cowplot::plot_grid(p_auc_study_3, p_auc_study_5)