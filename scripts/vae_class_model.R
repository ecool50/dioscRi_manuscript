compute_kernel <- function(x, y) {
  x_size <- tf$shape(x)[1]
  y_size <- tf$shape(y)[1]
  dim <- tf$shape(x)[2]
  
  tiled_x <- tf$tile(tf$reshape(x, tf$stack(list(x_size, 1L, dim))), tf$stack(list(1L, y_size, 1L)))
  tiled_y <- tf$tile(tf$reshape(y, tf$stack(list(1L, y_size, dim))), tf$stack(list(x_size, 1L, 1L)))
  
  return(tf$exp(-tf$reduce_mean(tf$square(tiled_x - tiled_y), axis = 2L) / tf$cast(dim, tf$float32)))
}

compute_mmd <- function(x, y, sigma_sqr = 1.0) {
  x_kernel <- compute_kernel(x, x)
  y_kernel <- compute_kernel(y, y)
  xy_kernel <- compute_kernel(x, y)
  
  return(tf$reduce_mean(x_kernel) + tf$reduce_mean(y_kernel) - 2 * tf$reduce_mean(xy_kernel))
}

train_vae_class_model <- function(train_dat, useMarkers, epochs = 80, 
                                  latent_dim = 2L,
                            original_dim = 27L, batch_size = 16) {
  
  tensorflow::set_random_seed(seed = 1994)
  # Normalize and prepare the training data
  x_train <- train_dat[, useMarkers] %>%
    as.matrix()
  
  original_dim <- ncol(train_dat)
  intermediate_dim <- original_dim - 4
  intermediate_dim_2 <- intermediate_dim - 4
  latent_dim <- intermediate_dim_2 - 3
  
  # Define the encoder
  
  encoder_inputs <- layer_input(shape = original_dim)
  
  x <- encoder_inputs %>%
    layer_dense(intermediate_dim, activation = "relu") %>%
    layer_dense(intermediate_dim_2, activation = "relu") 
  
  z_mean  <- x %>% layer_dense(latent_dim, name = "z_mean", 
                               activation = 'gelu')
  encoder <- keras_model(encoder_inputs, z_mean,
                         name = "encoder")
  
  decoder_outputs <- z_mean %>%
    layer_dense(intermediate_dim_2, activation = "relu") %>%
    layer_dense(intermediate_dim, activation = "relu") %>%
    layer_dense(original_dim, activation = "sigmoid")
  
  decoder <- keras_model(z_mean, decoder_outputs,
                         name = "decoder")
  
  
  model_vae <- new_model_class(
    classname = "VAE",
    
    initialize = function(encoder, decoder, ...) {
      super$initialize(...)
      self$encoder <- encoder
      self$decoder <- decoder
      self$sampler <- layer_sampler()
      self$total_loss_tracker <-
        metric_mean(name = "total_loss")
      self$reconstruction_loss_tracker <-
        metric_mean(name = "reconstruction_loss")
      self$mmd_loss_tracker <-
        metric_mean(name = "mmd_loss")
    },
    
    metrics = mark_active(function() {
      list(
        self$total_loss_tracker,
        self$reconstruction_loss_tracker,
        self$mmd_loss_tracker
      )
    }),
    
    train_step = function(data) {
      x <- data
      with(tf$GradientTape() %as% tape, {
        
        z_mean <- self$encoder(x)
        
        
        reconstruction <- decoder(z_mean)
        reconstruction_loss <-
          loss_binary_crossentropy(x, reconstruction) %>%
          sum(axis = c(1)) %>%
          mean()
        
        true_samples <- tf$random$normal(shape = tf$shape(z_mean))
        mmd_loss <- compute_mmd(true_samples, z_mean)
        total_loss <- reconstruction_loss + 0.1*mmd_loss
      })
      
      grads <- tape$gradient(total_loss, self$trainable_weights)
      self$optimizer$apply_gradients(zip_lists(grads, self$trainable_weights))
      
      self$total_loss_tracker$update_state(total_loss)
      self$reconstruction_loss_tracker$update_state(reconstruction_loss)
      self$mmd_loss_tracker$update_state(mmd_loss)
      
      list(total_loss = self$total_loss_tracker$result(),
           reconstruction_loss = self$reconstruction_loss_tracker$result(),
           mmd_loss = self$mmd_loss_tracker$result()
           )
    }
  )
  
  
  vae <- model_vae(encoder, decoder)
  lr_schedule <- learning_rate_schedule_exponential_decay(
    1e-3,
    decay_steps = 100000, 
    decay_rate = 0.95, 
    staircase = FALSE) 
  es_cb <- callback_early_stopping(
    min_delta = 1e-04, 
    monitor = 'total_loss',
    mode = 'min',
    patience = 5, 
    verbose = 1, 
    restore_best_weights = TRUE
  )
  
  # First stage training with Adam optimizer
  opt <- optimizer_rmsprop(learning_rate = lr_schedule, momentum = 0.9, centered = TRUE)
  vae %>% compile(optimizer = opt)
  vae %>% keras3::fit(x_train,
                      epochs = epochs,
                      shuffle = TRUE)
  
  return(list(vae = vae, encoder = encoder))
}


# Function to decode new samples
decode_samples <- function(new_samples, vae, latent_dim = 8L, batch_size = 16) {
  tensorflow::set_random_seed(seed = 1994)
  z_mean <- predict(vae$encoder, new_samples)

  # Step 3: Decode the Latent Representation
  decoded_samples <- predict(vae$decoder, z_mean) %>%
    as.data.frame()
  return(list(decoded = decoded_samples, encoded = z_mean))
}