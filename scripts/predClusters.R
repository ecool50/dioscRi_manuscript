predClusters <- function(train_x, train_y, nClust,
                         weights_path = NULL){
  
  tensorflow::set_random_seed(seed = 1994)
  l2_penalty <- 1e-04
  
  class_weights <- function(labels, case = 1, mu = 0.15) {
    # do bin count
    weights <- table(labels)
    if (case == 1){
      # sklearn.utils.class_weight.compute_class_weight approach
      # http://scikit-learn.org/stable/modules/generated/sklearn.utils.class_weight.compute_class_weight.html
      weights <- sum(weights) / (length(weights) * weights)
    } else if (case == 2) {
      weights <- log(mu * sum(weights) / weights)
      weights <- ifelse(weights < 1, 1, weights)
    } else if (case == 3) {
      weights <- ceiling(max(weights) / weights)
    } else {
      weights <- weights / sum(weights)
    }
    # create and return list
    setNames(as.list(weights), names(weights))
  }
  
  train_y <- factor(train_y)
  weights <- class_weights(as.numeric(train_y) - 1, case = 1)
  train_y <- to_categorical(as.numeric(train_y) - 1)
  # print(train_y)
  # print(weights)
  
  model <- keras_model_sequential() %>%
    layer_dense(units = 128, input_shape = ncol(train_x), 
                activation = 'mish', 
                kernel_regularizer = regularizer_l2(l2_penalty)) %>%
    layer_dense(units = 64, activation = 'mish', 
                kernel_regularizer = regularizer_l2(l2_penalty)) %>%
    layer_dense(units = 32, activation = 'mish',
                kernel_regularizer = regularizer_l2(l2_penalty)) %>%
    layer_dense(units = nClust, activation = 'softmax')
  
  lr_schedule <- learning_rate_schedule_exponential_decay(
    0.001, 
    decay_steps = 100000, 
    decay_rate = 0.96, 
    staircase = TRUE) 
  opt <- optimizer_adamax(learning_rate = lr_schedule)
  
  model %>%
    compile(optimizer = opt,
            loss = 'categorical_crossentropy',
            metrics = "AUC")
  
  # es_cb <- callback_early_stopping(patience = 15, monitor = 'val_loss')
  es_cb <- keras::callback_reduce_lr_on_plateau(
    monitor = "val_auc",
    min_delta = 1e-04,  # No improvement threshold for early stopping based on moving average
    patience = 15,   # Number of epochs to wait before stopping if no improvement
    verbose = 1
  )
  
  es_cb_1 <- keras::callback_early_stopping(
    min_delta = 1e-04,
    patience = 15, 
    verbose = 1
  )
  
  # Adjusted Number of Epochs and Batch Size
  model %>%
    fit(x = as.matrix(train_x),
        y = train_y, 
        epochs = 100,  # Increased number of epochs
        batch_size = 16,  # Adjusted batch size
        validation_split = 0.1, 
        callbacks = list(es_cb, es_cb_1),
        class_weights = weights)
  
  return(model)
  
  
}
