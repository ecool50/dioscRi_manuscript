library(tidymodels)
library(probably)
library(discrim)

tidymodels_prefer()
theme_set(theme_bw())
options(pillar.advice = FALSE, pillar.min_title_chars = Inf)

dat <- cbind(data, condition)
set.seed(8928)
split <- initial_split(dat, strata = condition)
dat_tr <- training(split)
dat_te <- testing(split)

dat_rs <- vfold_cv(dat_tr, strata = condition)


bayes_wflow <-
  workflow() %>%
  add_formula(condition ~ .) %>%
  add_model(parsnip::naive_Bayes(mode = "classification"))


cls_met <- metric_set(roc_auc, brier_class)
# We'll save the out-of-sample predictions to visualize them. 
ctrl <- control_resamples(save_pred = TRUE)

bayes_res <-
  bayes_wflow %>%
  fit_resamples(dat_rs, metrics = cls_met, control = ctrl)

collect_metrics(bayes_res)

collect_predictions(bayes_res) %>%
  ggplot(aes(.pred_0)) +
  geom_histogram(col = "white", bins = 40) +
  facet_wrap(~ condition, ncol = 1) +
  geom_rug(col = "blue", alpha = 1 / 2) + 
  labs(x = "Probability Estimate of 0")

cal_plot_breaks(bayes_res)

cal_plot_logistic(bayes_res)

logit_val <- cal_validate_logistic(bayes_res, metrics = cls_met, save_pred = TRUE)
collect_metrics(logit_val)


collect_predictions(logit_val) %>%
  filter(.type == "calibrated") %>%
  cal_plot_windowed(truth = condition, estimate = .pred_0, step_size = 0.025) +
  ggtitle("Logistic calibration via GAM")


set.seed(1212)
iso_val <- cal_validate_isotonic_boot(bayes_res, metrics = cls_met, 
                                      save_pred = TRUE, times = 25)
collect_metrics(iso_val)

collect_predictions(iso_val) %>%
  filter(.type == "calibrated") %>%
  cal_plot_windowed(truth = condition, estimate = .pred_0, step_size = 0.025) +
  ggtitle("Isotonic regression calibration")


beta_val <- cal_validate_beta(bayes_res, metrics = cls_met, save_pred = TRUE)
collect_metrics(beta_val)


collect_predictions(beta_val) %>%
  filter(.type == "calibrated") %>%
  cal_plot_windowed(truth = condition, estimate = .pred_0, step_size = 0.025) +
  ggtitle("Beta calibration")


dat_cal <- cal_estimate_beta(bayes_res)
bayes_fit <- bayes_wflow %>% fit(data = dat_tr)

dat_test_pred <- augment(bayes_fit, new_data = dat_te)
dat_test_pred %>% cls_met(condition, .pred_0)

dat_test_cal_pred <-
  dat_test_pred %>%
  cal_apply(dat_cal)
dat_test_cal_pred %>% dplyr::select(condition, starts_with(".pred_"))


dat_test_cal_pred %>% cls_met(condition, .pred_0)

dat_test_cal_pred %>%
  cal_plot_windowed(truth = condition, estimate = .pred_0, step_size = 0.025)
