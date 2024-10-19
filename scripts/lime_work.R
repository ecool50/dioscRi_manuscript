suppressPackageStartupMessages({ library (readxl)
  library (keras);library (lime);library (tidyquant)
  library (rsample);library (recipes);library (yardstick)
  library (corrr);library(knitr);library (DT) })
# Explaining a model based on tabular data
smp_size <- floor(0.75 * nrow(data))

## set the seed to make your partition reproducible
set.seed(123)
train_ind <- sample(seq_len(nrow(meta_features_fs)), size = smp_size)

train <- meta_features_fs[train_ind, ]
test <- meta_features_fs[-train_ind, ]

condition_num <- as.numeric(levels(condition))[condition]

tensorflow::set_random_seed(1994)
model_fs <- getModel(data = train)

es_cb <- callback_early_stopping(patience = 15, monitor = 'val_loss',
                                 restore_best_weights = TRUE)

model_fs %>% keras::fit(
  as.matrix(train), 
  to_categorical(as.numeric(condition[train_ind]) - 1),
  epochs = 20,
  batch_size = 16,
  validation_split = 0.1,
  callbacks = es_cb,
  verbose = 1
)

keras::evaluate(model_fs, as.matrix(test),
                to_categorical(as.numeric(condition[-train_ind]) - 1),
                batch_size = 16, verbose = 0, sample_weight = NULL)[[2]]


model_type.keras.engine.sequential.Sequential <- function(x, ...) {
  "classification"}

predict_model.keras.engine.sequential.Sequential <- function (x, newdata, type, ...) {
  pred <- predict(object = x, x = as.matrix(newdata))
  data.frame (Positive = pred[, 1], Negative = pred[, 2]) }

predict_model (x       = model_fs, 
               newdata = test, 
               type    = 'raw')

explainer <- lime::lime(
  x              = train, 
  model          = model_fs, 
  bin_continuous = FALSE)


explanation <- lime::explain (
  test[1:6, ], # Just to show first 10 cases
  explainer    = explainer, 
  n_permutations = 5000,
  feature_select = 'lasso_path',
  n_labels     = 1, # explaining a `single class`(Polarity)
  n_features   = 15, # returns top four features critical to each case
  kernel_width = sqrt(ncol(data)) * 0.75)

plot_features (explanation) +
  labs (title = "LIME: Feature Importance Visualization",
        subtitle = "Hold Out (Test) Set, First 10 Cases Shown")

plot_explanations (explanation) +
  labs (title = "LIME Feature Importance Heatmap",
        subtitle = "Hold Out (Test) Set, First 10 Cases Shown")

corrr_analysis <- test %>%
  mutate (Gensini = condition_num[-train_ind]) %>%
  correlate () %>%
  focus (Gensini) %>%
  dplyr::rename (feature = term) %>%
  arrange (abs(Gensini)) %>%
  mutate (feature = as_factor(feature))


corrr_analysis %>%
  
  ggplot (aes (x = Gensini, y = fct_reorder(feature, desc(Gensini)))) +
  geom_point () +
  # Positive Correlations - Contribute to Gensini--------------------------------------------
geom_segment (aes(xend = 0, yend = feature), 
              color = palette_light()[[2]], 
              data = corrr_analysis %>% filter(Gensini > 0)) +
  geom_point (color = palette_light()[[2]], 
              data = corrr_analysis %>% filter(Gensini > 0)) +
  # Negative Correlations - Prevent Gensini--------------------------------------------------
geom_segment (aes(xend = 0, yend = feature), 
              color = palette_light()[[1]], 
              data = corrr_analysis %>% filter(Gensini < 0)) +
  geom_point (color = palette_light()[[1]], 
              data = corrr_analysis %>% filter(Gensini < 0)) +
  # Vertical lines-------------------------------------------------------------------------
geom_vline (xintercept = 0, color = palette_light()[[5]], size = 1, linetype = 2) +
  geom_vline (xintercept = -0.25, color = palette_light()[[5]], size = 1, linetype = 2) +
  geom_vline (xintercept = 0.25, color = palette_light()[[5]], size = 1, linetype = 2) +
  # Aesthetics-----------------------------------------------------------------------------
theme_tq () +
  labs (title = "Gensini Correlation Analysis",
        subtitle = "Positive Correlations vs. Negative Correlations", 
        y = "Feature Importance")
