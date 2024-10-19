trainDAE <- function(train_dat, input_size = 27, encoding_dim = 25, batch_size = 32,
                     noise_factor = 0.1,l2_reg = 0,  
                     dropout_rate = 0.0, epochs = 100, useMarkers){
  
  train_dat <- as.matrix(train_dat[, useMarkers])
  input_size <- ncol(train_dat)
  tensorflow::set_random_seed(seed = 1994)
  input_layer <- layer_input(shape = c(input_size))

  train_dat_noisy <- train_dat + noise_factor * tf$random$normal(shape = tf$shape(train_dat))
  train_dat_noisy <- tf$clip_by_value(train_dat_noisy, 0, 1)
  
  # Encoder
  encoded_layer <- layer_dense(input_layer, units = encoding_dim, activation = 'relu',
                               kernel_regularizer = regularizer_l2(l2_reg)) %>%
    layer_dropout(rate = dropout_rate)
  
  encoded_layer <- layer_dense(encoded_layer, units = encoding_dim, activation = 'relu',
                              kernel_regularizer = regularizer_l2(l2_reg)) %>%
    layer_dropout(rate = dropout_rate)
  
  decoded_layer <- layer_dense(encoded_layer, units = input_size, activation = 'sigmoid',
                              kernel_regularizer = regularizer_l2(l2_reg)) %>%
    layer_dropout(rate = dropout_rate)

  
  # Autoencoder model
  autoencoder_model <- keras_model(inputs = input_layer, outputs = decoded_layer)
  
  # Encoder model
  encoder_model <- keras_model(inputs = input_layer, outputs = encoded_layer)
  
  
  lr_schedule <- learning_rate_schedule_exponential_decay(
    1e-3,
    decay_steps = 100000, 
    decay_rate = 0.95, 
    staircase = FALSE) 
  opt <- optimizer_rmsprop(learning_rate = lr_schedule,
                           momentum = 0.9, centered = TRUE)
  
  # es_cb <- callback_early_stopping(patience = 15, monitor = 'val_loss')
  es_cb <- callback_reduce_lr_on_plateau(
    monitor = "val_loss",
    min_delta = 1e-04,  # No improvement threshold for early stopping based on moving average
    patience = 15,   # Number of epochs to wait before stopping if no improvement
    verbose = 1
  )
  
  es_cb_1 <- callback_early_stopping(
    min_delta = 1e-04,
    patience = 5, 
    verbose = 1
  )
  
  autoencoder_model %>% compile(optimizer = opt, loss = 'mse')
  autoencoder_model %>% fit(x = train_dat_noisy,
                            y = train_dat,
                            validation_split = 0.1,
                            shuffle = TRUE,
                            batch_size = batch_size,
                            callbacks = list(es_cb_1),
                            epochs = epochs)
  
  return(autoencoder_model)
}