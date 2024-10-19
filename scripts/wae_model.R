mmd_penalty <- function(pz, qz, z_dim = 16L, sigma_z = 1.0) {
  # Calculates the unbiased U-statistic of the MMD between the prior and the aggregated posterior.
  # The kernel function used is the inverse multiquadratic (IMQ).
  
  batch_size <- tf$shape(qz)[1]
  
  # ||x - y||**2 = ||x||**2 + ||y||**2 - <x, y> is used to calculate distances
  norms_pz <- tf$reduce_sum(tf$square(pz), axis = 1L, keepdims = TRUE)
  dotprods_pz <- tf$matmul(pz, pz, transpose_b = TRUE)
  distances_pz <- norms_pz + tf$transpose(norms_pz) - 2 * dotprods_pz
  
  norms_qz <- tf$reduce_sum(tf$square(qz), axis = 1L, keepdims = TRUE)
  dotprods_qz <- tf$matmul(qz, qz, transpose_b = TRUE)
  distances_qz <- norms_qz + tf$transpose(norms_qz) - 2 * dotprods_qz
  
  dotprods <- tf$matmul(qz, pz, transpose_b = TRUE)
  distances <- norms_qz + tf$transpose(norms_pz) - 2 * dotprods
  
  cbase <- tf$constant(2 * z_dim * sigma_z, dtype = tf$float32)
  stat <- tf$constant(0, dtype = tf$float32)
  nf <- tf$cast(batch_size, dtype = tf$float32)
  
  # Looking at "various scales" using the property that the sum of two positive definite kernels is still
  # a positive definite kernel.
  scales <- c(0.1, 0.2, 0.5, 1, 2, 5, 10)
  for (scale in scales) {
    C <- cbase * scale
    res1 <- C / (C + distances_qz)
    res1 <- res1 + C / (C + distances_pz)
    res1 <- tf$multiply(res1, 1 - tf$eye(batch_size, dtype = tf$float32))
    res1 <- tf$reduce_sum(res1) / (nf * nf - nf)
    
    res2 <- C / (C + distances)
    res2 <- tf$reduce_sum(res2) * 2 / (nf * nf)
    
    stat <- stat + (res1 - res2)
  }
  
  stat
}


robust_kernel <- function(alpha, sigma = 1) {
  tf_2pi = tf$constant(sqrt(2*pi), dtype=tf$float32)
   tf$exp(tf_2pi*(-(alpha * alpha) / (2 * sigma * sigma)))
}

correntropy_loss <- function(y_true, y_pred, sigma= 1) {
  loss_corr <- -tf$reduce_mean(robust_kernel(y_pred - y_true, sigma = sigma))
  loss_corr
}


compute_kernel <- function(x, y) {
  x_size <- tf$shape(x)[1]
  y_size <- tf$shape(y)[1]
  dim <- tf$shape(x)[2]
  
  tiled_x <- tf$tile(
    tf$reshape(x, tf$stack(list(x_size, 1L, dim))),
    tf$stack(list(1L, y_size, 1L))
  )
  
  tiled_y <- tf$tile(
    tf$reshape(y, tf$stack(list(1L, y_size, dim))),
    tf$stack(list(x_size, 1L, 1L))
  )
  
  return(tf$exp(-tf$reduce_mean(tf$square(tiled_x - tiled_y), axis = 2L) / tf$cast(dim, tf$float32)))
}

compute_mmd <- function(x, y, sigma_sqr = 1.0) {
  x_kernel <- compute_kernel(x, x)
  y_kernel <- compute_kernel(y, y)
  xy_kernel <- compute_kernel(x, y)
  
  return(tf$reduce_mean(x_kernel) + tf$reduce_mean(y_kernel) - 2 * tf$reduce_mean(xy_kernel))
}



train_wae_mmd_model <- function(train_dat, useMarkers, epochs = 100, latent_dim = 2L, 
                                batch_size = 16, sigma = 1) {
  
  tensorflow::set_random_seed(seed = 1994)
  # Normalize and prepare the training data
  x_train <- train_dat[, useMarkers] %>%
    as.matrix()
  
  original_dim <- ncol(train_dat)
  
  # Updated dimensions for additional layers
  intermediate_dim <- 24
  # intermediate_dim_2 <- intermediate_dim - 4
  latent_dim <- 16
  
  # Define the encoder with more layers
  encoder_inputs <- layer_input(shape = original_dim)
  
  encoder_outputs <- encoder_inputs %>%
    layer_dense(intermediate_dim, activation = "relu") %>%
    # layer_dense(intermediate_dim_2, activation = "relu") %>%
    layer_dense(units = latent_dim, activation = "linear") # Latent representation
  
  encoder <- keras_model(encoder_inputs, encoder_outputs,
                         name = "encoder")
  
  # Define the decoder with more layers
  latent_inputs <- layer_input(shape = c(latent_dim))
  
  decoder_outputs <- encoder_outputs %>%
    # layer_dense(intermediate_dim_2, activation = "relu") %>%
    layer_dense(intermediate_dim, activation = "relu") %>%
    layer_dense(original_dim, activation = "sigmoid") # Output layer
  
  decoder <- keras_model(encoder_outputs, decoder_outputs,
                         name = "decoder")
  
  # Custom VAE model class definition remains the same
  model_vae <- new_model_class(
    classname = "VAE",
    
    initialize = function(encoder, decoder, ...) {
      super$initialize(...)
      self$encoder <- encoder
      self$decoder <- decoder
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
        
        z_mean %<-% self$encoder(x)
        true_samples <- tf$random$normal(shape = c(tf$cast(batch_size, tf$int64), 
                                                   tf$cast(latent_dim, tf$int64)))
        
        reconstruction <- decoder(z_mean)
        
        reconstruction_loss <- correntropy_loss(x, reconstruction)
        
        mmd_loss <- mmd_penalty(true_samples, z_mean)
        
        total_loss <- reconstruction_loss + 0.1*mmd_loss
      })
      
      grads <- tape$gradient(total_loss, self$trainable_weights)
      self$optimizer$apply_gradients(zip_lists(grads, self$trainable_weights))
      
      self$total_loss_tracker$update_state(total_loss)
      self$reconstruction_loss_tracker$update_state(reconstruction_loss)
      self$mmd_loss_tracker$update_state(mmd_loss)
      
      list(total_loss = self$total_loss_tracker$result(),
           reconstruction_loss = self$reconstruction_loss_tracker$result(),
           mmd_loss = self$mmd_loss_tracker$result())
    }
  )
  
  lr_schedule <- learning_rate_schedule_exponential_decay(
    1e-3,
    decay_steps = 100000, 
    decay_rate = 0.95, 
    staircase = FALSE) 
  
  # First stage training with Adam optimizer
  opt <- optimizer_rmsprop(learning_rate = lr_schedule,
                           momentum = 0.9, centered = TRUE)
  
  vae <- model_vae(encoder, decoder)
  vae %>% compile(optimizer = opt)
  
  # Early stopping callback for the first stage
  es_cb <- callback_early_stopping(
    min_delta = 1e-04, 
    monitor = 'total_loss',
    mode = 'min',
    patience = 5, 
    verbose = 1, 
    restore_best_weights = TRUE
  )
  
  # Train the model with Adam optimizer
  vae %>% keras3::fit(
    x_train,
    batch_size = batch_size,
    epochs = epochs,
    # callbacks = list(es_cb),
    shuffle = TRUE
  )
  
  return(list(vae = vae, encoder = encoder))
}


# Function to decode new samples
decode_samples_wae <- function(new_samples, vae, batch_size = 16) {
  tensorflow::set_random_seed(seed = 1994)
  encoded <- predict(vae$encoder, new_samples)
  z_mean <- encoded
  
  # Step 3: Decode the Latent Representation
  decoded_samples <- predict(vae$decoder, z_mean, batch_size = batch_size) %>%
    as.data.frame()
  return(list(decoded = decoded_samples, encoded = z_mean))
}