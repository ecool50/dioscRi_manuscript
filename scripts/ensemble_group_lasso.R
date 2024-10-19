# Run treekor
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

# order group (for gglasso)
ord <- order(group)
groupord <- group[ord]
# order var according to group
varord <- var[ord]

# transform group to have consecutive numbers (for gglasso)
groupb <- cumsum(!duplicated(groupord))

## Proportions
train_x_prop <- data[, varord] %>%
  as.matrix() 
test_prop_study_3 <- data_3[, varord] %>%
  as.matrix()
test_prop_study_5 <- data_5[, varord] %>%
  as.matrix()

## Logit Proportions
group_2 <- group + max(groupb)
train_x_logit <- data_logit[, varord] %>%
  as.matrix()
test_logit_study_3 <- data_3_logit[, varord] %>%
  as.matrix()
test_logit_study_5 <- data_5_logit[, varord] %>%
  as.matrix()

# Marker mean
d <- dist(scale(t(t_mean)))
hc3 <- hclust(d, method = 'average')
res3 <- MLGL:::preliminaryStep(hc3)
var_3 <- res3$var
group_3 <- res3$group
group_3 <- group_3 + max(group_2)
train_x_mean <- t_mean[, var_3] %>%
  as.matrix() 
test_mean_study_3 <- t_3_mean[, var_3] %>%
  as.matrix()
test_mean_study_5 <- t_5_mean[, var_3] %>%
  as.matrix()

# consolidate the groups
group_final <- c(group, group_2, group_3)

# setup training and testing datasets
train_x <- cbind(train_x_prop, train_x_logit, train_x_mean)
test_x_study_3 <- cbind(test_prop_study_3, test_logit_study_3, test_mean_study_3) 
test_x_study_5 <- cbind(test_prop_study_5, test_logit_study_5, test_mean_study_5) 
y <- condition %>%
  as.numeric()

train_x <- scale(train_x)

test_x_study_3 <- scale(test_x_study_3, center=attr(train_x, "scaled:center"), 
                scale=attr(train_x, "scaled:scale"))

# fit model and evaluate
fit <- cv.grpreg(X = train_x, y  = y, family = "binomial", type.meaure = 'auc',
                 penalty = 'grSCAD', alpha = 1,
                 group = group_final,
                 nfolds = 10, seed = 1994)
plot(fit)

predict(fit, type="group", lambda=fit$lambda.min)

## Predict on study 3
pred_study_3 <- predict(fit, X = test_x_study_3,
                     type = "response", which = 'lambda.min')

pred_study_3 <- prediction(pred_study_3, condition_3)
pred_study_3_perf <- ROCR::performance(pred_study_3, "tpr", "fpr")


plt_dat_3 = data.frame(
  FPR = pred_study_3_perf@x.values[[1]],
  TPR = pred_study_3_perf@y.values[[1]]
)

auc_perf_study_3 <- ROCR::performance(pred_study_3, measure = "auc")
auc_study_3 <- auc_perf_study_3@y.values[[1]]

# generate AUC plot
p_auc_study_3 <- ggplot(plt_dat_5, aes(x = FPR, y = TPR)) +
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
  )  + ggtitle(paste("Study 3 - AUC =", round(auc_study_3, 2)))




## Predict on study 5
pred_study_5 <- predict(fit, X = test_x_study_5,
                        type = "response", which = 'lambda.min')

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
  )  + ggtitle(paste("Study 5 - AUC =", round(auc_study_5, 2)))


## combine both plots
cowplot::plot_grid(p_auc_study_3, p_auc_study_5)




coefs <- coef(fit, s = 'lambda.min') %>%
  as.matrix() %>%
  as.data.frame()
coefs$feature <- rownames(coefs)

coefs <- coefs[coefs$V1 != 0, ]

res_mp <- MP_gLasso(cv_object = fit, group = group_final, lambda.type = "min", sort.type = "max")
girafe(code = print(res_mp$plot), width_svg = 8, height_svg = 8)


res_dat <- res_mp$data
res_dat$beta_name <- janitor::make_clean_names(res_dat$beta_name)
res_dat$beta_name <- clean_strings(res_dat$beta_name)

# get the mean groups
res_dat_mean <- res_dat[res_dat$beta_name %like% 'mean', ]

# get the prop groups
res_dat_prop <- res_dat[res_dat$beta_name %like% 'prop', ]

# get the logit groups
res_dat_logit <- res_dat[res_dat$beta_name %like% 'logit', ]

## Highlight all the single groups from mean
h_mean_single <- res_dat_mean[table(res_dat_mean$groups) == 1, ]$beta_name
hc3$labels <- janitor::make_clean_names(hc3$labels)
dend_mean <- as.dendrogram(hc3)
# create a dendrogram with colored labels:
dend_mean <- color_labels(dend_mean, labels = h_mean_single ,
                          col = rainbow(n = length(h_mean_single)))
# ploting it
# plot(dend_mean)

## Highlight all the single groups for prop
h_prop_single <- res_dat_prop[table(res_dat_prop$groups) == 1, ]$beta_name
dend_prop <- hc %>%
  as.dendrogram()
dend_prop <- color_labels(dend_prop, labels = gsub('_prop', '', h_prop_single),
                          col = rainbow(n = length(h_prop_single)))
# plot(dend_prop)


## Highlight all the single groups for logit
h_logit_single <- res_dat_logit[table(res_dat_logit$groups) == 1, ]$beta_name
dend_logit <- hc %>%
  as.dendrogram()
dend_logit <- color_labels(dend_logit, labels = gsub('_logit', '', h_logit_single),
                          col = rainbow(n = length(h_logit_single)))
# plot(dend_logit)
# 

x_sub_48 <- train_x[, which(group_final %in% 48)] %>% as.data.frame()
  
x_sub_51 <- train_x[, which(group_final %in% 51)] %>% as.data.frame()

x_sub <- data %>%
  dplyr::select(cluster_29, cluster_2)
colnames(x_sub) <- c('cluster_29_prop', "cluster_2_prop")
x_sub <- sin(x_sub/2)^2

x_sub$Condition <- condition

p1 <- x_sub %>% 
  ggplot(aes(cluster_29_prop, cluster_2_prop)) +
  geom_point(aes(color = Condition))

p2 <- x_sub %>%
  melt() %>%
  ggplot(aes(variable, value)) +
  geom_boxplot(aes(color = Condition))

## Fit a logistic model on groups 48 and 51
fit_48 <- cv.glmnet(x = as.matrix(x_sub_48), y = y, type.measure = 'auc', alpha = 1,
                    family = 'binomial')
plot(fit_48)
coef(fit_48)

fit_51 <- cv.glmnet(x = as.matrix(x_sub_51), y = y, type.measure = 'auc', alpha = 1,
                    family = 'binomial')
coef(fit_51)

