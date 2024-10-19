computePred <- function(train_x, test_x = NULL, train_y, test_y = NULL,
         weights_path = NULL){
  
  model <- keras_model_sequential(input_shape = ncol(train_x)) %>%
    layer_dense(256, activation = "relu") %>%
    layer_batch_normalization() %>%
    layer_dense(256, activation = "relu") %>%
    layer_dropout(0.3) %>%
    layer_batch_normalization() %>%
    layer_dense(256, activation = "relu") %>%
    layer_dropout(0.3) %>%
    layer_batch_normalization() %>%
    layer_dense(2, activation = "sigmoid")
  
  metrics <- list(
    # metric_false_negatives(name = "fn"),
    # metric_false_positives(name = "fp"),
    metric_auc(name = "AUC")
    # metric_precision(name = "precision"),
    # metric_recall(name = "recall")
  )
  model %>% compile(
    optimizer = optimizer_adam(weight_decay = 0.001),
    loss = "binary_crossentropy",
    metrics = metrics
  )
  # class_weight <- list("0" = weight_for_0,
  #                      "1" = weight_for_1)
  es_cb <- keras::callback_early_stopping(patience = 5)
  checkpoint <- keras::callback_model_checkpoint(filepath =
                                                   paste0(weights_path,'.hdf5'),
                                                 save_best_only = T)
  
  # train_features <- as.matrix(train_df[feature_names])
  # train_targets <- as.matrix(train_df$Class)
  # validation_data <- list(
  #   as.matrix(val_df[feature_names]),
  #   as.matrix(val_df$Class))
  
  model %>%
    fit(as.matrix(train_x), to_categorical(train_y),
        # validation_data = validation_data,
        # class_weight = class_weight,
        validation_split = 0.5,
        batch_size = 60, epochs = 100,
        callbacks = list(checkpoint, es_cb),
        verbose = 2)
}
