dat <- cbind(dat, condition)

### Random Forest
rf_mod <- 
  rand_forest(trees = 1000) %>% 
  set_engine("ranger") %>% 
  set_mode("classification")

set.seed(345)
folds <- vfold_cv(dat, v = 5)

rf_wf <- 
  workflow() %>%
  add_model(rf_mod) %>%
  add_formula(condition ~ .)

set.seed(456)
rf_fit_rs <- 
  rf_wf %>% 
  fit_resamples(folds)

collect_metrics(rf_fit_rs)

### Deep Learning

biv_rec <- 
  recipe(condition ~ ., data = dat) %>%
  step_normalize(all_predictors())

nnet_spec <- 
  mlp(epochs = 1000, hidden_units = 10, penalty = 0.01, learn_rate = 0.01) %>% 
  set_engine("brulee", validation = 0) %>% 
  set_mode("classification")

nnet_wflow <- 
  biv_rec %>% 
  workflow(nnet_spec)

set.seed(345)
folds <- vfold_cv(dat, v = 5)

nn_wf <- 
  workflow() %>%
  add_model(nnet_spec) %>%
  add_formula(condition ~ .)

set.seed(456)
nn_fit_rs <- 
  nn_wf %>% 
  fit_resamples(folds)

collect_metrics(nn_fit_rs)

### Model Tunning
tune_spec <- 
  decision_tree(
    cost_complexity = tune(),
    tree_depth = tune()
  ) %>% 
  set_engine("rpart") %>% 
  set_mode("classification")

tree_grid <- grid_regular(cost_complexity(),
                          tree_depth(),
                          levels = 5)

set.seed(345)

tree_wf <- workflow() %>%
  add_model(tune_spec) %>%
  add_formula(condition ~ .)

tree_res <- 
  tree_wf %>% 
  tune_grid(
    resamples = folds,
    grid = tree_grid
  )

tree_res %>% 
  collect_metrics()
