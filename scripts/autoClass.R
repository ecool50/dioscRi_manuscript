robust_kernel <- function(alpha, sigma = 0.22) {
  1 / (op_sqrt(2 * pi * sigma)) * 
    op_exp(- op_square(alpha) / (2 * sigma * sigma))
}

correntropy_loss <- function(y_true, y_pred) {
  loss_corr <- -op_sum(robust_kernel(y_pred - y_true))
  loss_corr
}

autoClass <- function(train_dat, input_size = 27, encoding_dim = 25, 
                      l2_reg = 0,  dropout_rate = 0.0, 
                      classifier_weight = 0.9, epochs = 300,
                      full_dat, weights_path = NULL, nclust = 9) {
  
  tensorflow::set_random_seed(seed = 1994)
  
  train_dat <- train_dat %>% as.matrix()
  # batch_dat <- batch_dat %>% as.matrix()
  
  # noise_filter <- ruta::noise_gaussian()
  # train_dat_noisy <- ruta::apply_filter(noise_filter, train_dat)
  
  input_layer <- layer_input(shape = c(input_size))
  mid_layer <- input_layer
  
  # Define layer sizes based on encoding dimension
  # size1 <- max(encoding_dim * 2, input_size)
  size1 <- encoding_dim * 1.5
  size2 <- encoding_dim
  
  # Encoder
  # mid_layer <- layer_dense(mid_layer, units = size1, activation = 'mish',
  #                          kernel_regularizer = regularizer_l2_l2(l2_l2_reg))
  encode_layer <- layer_dense(mid_layer, units = size1, activation = 'relu',
                           kernel_regularizer = regularizer_l2(l2_reg)) %>%
    layer_dropout(rate = dropout_rate)
  
  
  bottle_neck_layer <- layer_dense(encode_layer, units = size2, activation = 'relu',
                           kernel_regularizer = regularizer_l2(l2_reg))
  # mid_layer <- layer_dropout(mid_layer, dropout_rate)
  
  decode_layer <- bottle_neck_layer
  
  # Decoder
  decode_layer <- layer_dense(decode_layer, units = size1, activation = 'relu',
                           kernel_regularizer = regularizer_l2(l2_reg)) %>%
    layer_dropout(rate = dropout_rate)
  
  output_layer <- layer_dense(decode_layer, units = input_size, 
                              activation = 'softplus',
                              name = 'reconstruction', 
                              kernel_regularizer = regularizer_l2(l2_reg))
  
  neck_layer <- encode_layer
  classifier_layer <- layer_dense(neck_layer, units = nclust+1, activation = 'softmax',
                                  name = 'classification', 
                                  kernel_regularizer = regularizer_l2(l2_reg))
  
  if(classifier_weight == 0){
    model <- keras_model(inputs = input_layer, outputs = output_layer)
  } else{
    model <- keras_model(inputs = input_layer, outputs = c(output_layer, classifier_layer))
  }
  
  
  es_cb <- callback_early_stopping(patience = 15, monitor = 'val_loss')
  
  initial_learning_rate <- 0.001
  lr_schedule <- learning_rate_schedule_exponential_decay(
    initial_learning_rate = initial_learning_rate,
    decay_steps = 100000,
    decay_rate = 0.96
  )
  opt <- optimizer_adam(learning_rate = lr_schedule, 
                         beta_1 = 0.975, 
                         beta_2 = 0.999)
  
  if(classifier_weight == 0){
    model %>% compile(optimizer = opt, loss = correntropy_loss)
    model %>% fit(x = train_dat,
                         y = train_dat,
                         validation_split = 0.1,
                         shuffle = TRUE,
                         batch_size = 16,
                         epochs = 100,
                         callbacks = list(es_cb))
    
    predictions <- predict(model, as.matrix(full_dat))
    norm_data <- predictions %>% as.data.frame()
  } else{
    clusters <- kmeans(train_dat, centers = nclust)$cluster
    # clusters <- FuseSOM::runFuseSOM(train_dat, numClusters = nclust)$clusters
    # clusters <- gsub("cluster_", "", clusters)
    clusters <- to_categorical(clusters)
    model %>% keras::compile(optimizer = opt, 
                      loss = list('reconstruction' = correntropy_loss,
                                  'classification' = correntropy_loss),
                      loss_weights = list('reconstruction' = 1-classifier_weight,
                                          'classification' = classifier_weight))
    model %>% keras::fit(x = train_dat,
                         y = list('reconstruction' = train_dat,
                                  'classification' = clusters),
                         validation_split = 0.1,
                         batch_size = 16,
                         shuffle = TRUE,
                         epochs = epochs,
                         callbacks = list(es_cb))
    
    
  }
  

  return(list(autoencoder_model = model))
}

autoEncoderModel <- function(train_dat, input_size = 27, encoding_dim = 25, 
                             l2_reg = 0,  dropout_rate = 0.0, epochs = 100){
  
  train_dat <- as.matrix(train_dat)
  input_size <- ncol(train_dat)
  tensorflow::set_random_seed(seed = 1994)
  input_layer <- layer_input(shape = c(input_size))
  mid_layer <- input_layer
  
  # Define layer sizes based on encoding dimension
  size1 <- input_size - 4
  size2 <- size1 - 4
  size3 <- size2 + 7
  if(size1 > input_size){
    size1 <- closest_power_of_2_less(input_size)
  }
  
  # Encoder
  encoded_layer <- layer_dense(mid_layer, units = size1, activation = 'relu',
                               kernel_regularizer = regularizer_l2(l2_reg)) %>%
    layer_dropout(rate = dropout_rate)
  
  decode_layer <- encoded_layer
  
  decode_layer <- layer_dense(decode_layer, units = size2, activation = 'relu',
                               kernel_regularizer = regularizer_l2(l2_reg)) %>%
    layer_dropout(rate = dropout_rate)
  
  decode_layer <- layer_dense(decode_layer, units = size1, activation = 'relu',
                              kernel_regularizer = regularizer_l2(l2_reg)) %>%
    layer_dropout(rate = dropout_rate)
  
  # Decoder
  output_layer <- layer_dense(decode_layer, units = input_size, 
                              activation = 'sigmoid', 
                              kernel_regularizer = regularizer_l2(l2_reg),
                              bias_regularizer = regularizer_l2(l2_reg),
                              name = 'reconstruction')
  
  # Autoencoder model
  autoencoder_model <- keras_model(inputs = input_layer, outputs = output_layer)
  
  # Encoder model
  encoder_model <- keras_model(inputs = input_layer, outputs = encoded_layer)
  
  
  lr_schedule <- learning_rate_schedule_exponential_decay(
    0.001, 
    decay_steps = 100000, 
    decay_rate = 0.96, 
    staircase = FALSE) 
  opt <- optimizer_adam(learning_rate = lr_schedule)
  
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
  
  autoencoder_model %>% compile(optimizer = opt, loss = correntropy_loss)
  autoencoder_model %>% fit(x = train_dat,
                                     y = train_dat,
                                     validation_split = 0.1,
                                     shuffle = TRUE,
                                     batch_size = 32,
                                     callbacks = list(es_cb, es_cb_1),
                                     epochs = epochs)
    
  return(list(autoencoder_model = autoencoder_model, encoder_model = encoder_model))
}
