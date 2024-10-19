## Data = considering that we have a data frame named dataF, with its first column being the class
dat <- cbind(data_logit)
colnames(dat) <- janitor::make_clean_names(colnames(dat))

dat_test <- cbind(data_3_logit)
colnames(dat_test) <- janitor::make_clean_names(colnames(dat_test))
dat_test <- as.matrix(dat_test)

hc <- mhclust(t(dat), verb = T)
hc <- fixNonMonotHca(hc)
dend <- as.dendrogram(hc)
dend %>% 
  set("labels_cex", 1) %>%
  plot(type = "rectangle", ylab = "Height")

res <- MLGL:::preliminaryStep(hc)

# order group (for gglasso)
group <- res$group
var <- res$var
weights <- res$weight
# ord <- order(group)
# groupord <- group[ord]
# # order var according to group
# varord <- var[ord]
# 
# # transform group to have consecutive numbers (for gglasso)
# groupb <- cumsum(!duplicated(groupord))
group <- paste0('Group ', as.character(group)) %>% 
  as.factor()
# heights <- findSequenceIndices(res$heights)

X <- as.matrix(dat[, var])

x <- as.matrix(dat[, var]) # Removes class
x_test <- as.matrix(dat_test[,var])
y <- condition %>%
  as.numeric()# Only class


# sds <- colSds(x)
# x <- scale(x, center = T, scale = sds)
# set.seed(999)
fit <- cv.grpreg(X = x, y  = y, family = "binomial", type.meaure = 'auc',
                 penalty = 'grLasso',
                 group.multiplier = as.numeric(sqrt(table(group))),
                 group = group, 
                 nfolds = 10, seed = 1994)
plot(fit)

predict(fit, type="group", lambda=fit$lambda.min)

# x_test <- scale(x_test, center = T, scale = sds)
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

MP_gLasso(cv_object = fit, group = group, lambda.type = "min", sort.type = "mean")
