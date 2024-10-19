#Generate clustering result, the input matrix has rows as samples and columns as genes

set.seed(1994)
dat <- x_train %>%
  dplyr::group_by(cellTypes) %>%
  slice_sample(n = 500)

val_set <- vfold_cv(dat, v = 5, repeats = 5)

lr_mod <- 
  multinom_reg(penalty = tune(), mixture = 1) %>% 
  set_engine("nnet")


lr_recipe <- 
  recipe(cellTypes ~ ., data = dat) %>% 
  step_zv(all_predictors()) %>% 
  step_normalize(all_predictors())

lr_workflow <- 
  workflow() %>% 
  add_model(lr_mod) %>% 
  add_recipe(lr_recipe)

lr_reg_grid <- tibble(penalty = 10^seq(-4, -1, length.out = 30))

lr_reg_grid %>% top_n(-5) # lowest penalty values

set.seed(1994)
lr_res <- 
  lr_workflow %>% 
  tune_grid(val_set,
            grid = lr_reg_grid,
            control = control_grid(save_pred = TRUE),
            metrics = metric_set(roc_auc))


lr_plot <- 
  lr_res %>% 
  collect_metrics() %>% 
  ggplot(aes(x = penalty, y = mean)) + 
  geom_point() + 
  geom_line() + 
  ylab("Area under the ROC Curve") +
  scale_x_log10(labels = scales::label_number())

lr_plot 

lr_best <-
  lr_res %>% 
  select_best()

lr_res %>% 
  collect_predictions(parameters = lr_best) %>% 
  roc_auc(condition, .pred_0)

lr_auc <- 
  lr_res %>% 
  collect_predictions(parameters = lr_best) %>% 
  roc_curve(condition, .pred_0) %>% 
  mutate(model = "Logistic Regression")

autoplot(lr_auc)


top_models <-
  lr_res %>% 
  gte(n = 15) %>% 
  arrange(penalty) 
top_models



