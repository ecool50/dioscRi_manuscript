dat <- cbind(meta_features_fs, condition)

data_split <- initial_split(dat, strata = "condition")
data_train <- training(data_split)
data_test <- testing(data_split)

# dat_upsampled <- recipe(~., dat) %>%
#   step_upsample(condition, over_ratio = 0.5) %>%
#   prep() %>%
#   bake(new_data = NULL)

set.seed(1994)
val_set <- vfold_cv(dat, v = 5, repeats = 1, strata = condition)

xgb_spec <- boost_tree(
  trees = 1000,
  tree_depth = tune(), min_n = tune(),
  loss_reduction = tune(),                     ## first three: model complexity
  sample_size = tune(), mtry = tune(),         ## randomness
  learn_rate = tune()                          ## step size
) %>%
  set_engine("xgboost") %>%
  set_mode("classification")

xgb_spec


xgb_grid <- grid_latin_hypercube(
  tree_depth(),
  min_n(),
  loss_reduction(),
  sample_size = sample_prop(),
  finalize(mtry(), dat),
  learn_rate(),
  size = 2
)

xgb_wf <- workflow() %>%
  add_formula(condition ~ .) %>%
  add_model(xgb_spec)

xgb_wf

doParallel::registerDoParallel()

set.seed(234)
xgb_res <- tune_grid(
  xgb_wf,
  resamples = val_set,
  grid = xgb_grid,
  control = control_grid(save_pred = TRUE)
)

xgb_res
