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

posTANH <- function(x){
  op_tanh(x) + 1
}

liSHT <- function(x){
  x * op_tanh(x)
}

train_vae_kl_model <- function(train_dat, useMarkers, epochs = 80, 
                               latent_dim = 2L, lambda = 0.1, val_dat,
                               original_dim = 27L, batch_size = 16) {
  
  tensorflow::set_random_seed(seed = 1994)
  
  # Normalize and prepare the training data
  x_train <- train_dat[, useMarkers] %>%
    as.matrix()
  x_val <- val_dat[, useMarkers] %>%
    as.matrix()
  
  original_dim <- ncol(train_dat)
  intermediate_dim <- original_dim - 4
  intermediate_dim_2 <- intermediate_dim - 4
  latent_dim <- intermediate_dim_2 - 3
  
  # Define the encoder
  
  encoder_inputs <- layer_input(shape = original_dim)
  
  x <- encoder_inputs %>%
    layer_dense(intermediate_dim, activation = "mish") %>%
    layer_dense(intermediate_dim_2, activation = "mish")
  
  z_mean    <- x %>% layer_dense(latent_dim, name = "z_mean")
  z_log_var <- x %>% layer_dense(latent_dim, name = "z_log_var")
  encoder <- keras_model(encoder_inputs, list(z_mean, z_log_var),
                         name = "encoder")
  
  layer_sampler <- new_layer_class(
    classname = "Sampler",
    call = function(z_mean, z_log_var) {
      epsilon <- tf$random$normal(shape = tf$shape(z_mean))
      z_mean + exp(0.5 * z_log_var) * epsilon }
  )
  
  latent_inputs <- layer_input(shape = c(latent_dim))
  
  decoder_outputs <- latent_inputs %>%
    layer_dense(intermediate_dim_2, activation = "mish") %>%
    layer_dense(intermediate_dim, activation = "mish") %>%
    layer_dense(original_dim, activation = "sigmoid")
  
  decoder <- keras_model(latent_inputs, decoder_outputs,
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
      self$kl_loss_tracker <-
        metric_mean(name = "kl_loss")
    },
    
    metrics = mark_active(function() {
      list(
        self$total_loss_tracker,
        self$reconstruction_loss_tracker,
        self$kl_loss_tracker
      )
    }),
    
    train_step = function(data) {
      with(tf$GradientTape() %as% tape, {
        data <- data[[1]]
        c(z_mean, z_log_var) %<-% self$encoder(data)
        z <- self$sampler(z_mean, z_log_var)
        
        reconstruction <- decoder(z)
        # reconstruction_loss <- tf$reduce_mean(tf$square(data - reconstruction))
        reconstruction_loss <- loss_binary_crossentropy(data, reconstruction) %>%
        sum(axis = c(1)) %>%
        mean()
        # true_samples <- tf$random$normal(shape = tf$shape(z_mean))
        kl_loss <- -0.5 * (1 + z_log_var - z_mean^2 - exp(z_log_var))
        total_loss <- reconstruction_loss + lambda*kl_loss
      })
      
      grads <- tape$gradient(total_loss, self$trainable_weights)
      self$optimizer$apply_gradients(zip_lists(grads, self$trainable_weights))
      
      self$total_loss_tracker$update_state(total_loss)
      self$reconstruction_loss_tracker$update_state(reconstruction_loss)
      self$kl_loss_tracker$update_state(kl_loss)
      
      list(total_loss = self$total_loss_tracker$result(),
           reconstruction_loss = self$reconstruction_loss_tracker$result(),
           kl_loss = self$kl_loss_tracker$result())
    },
    
    # Custom validation step
    test_step = function(data) {
      data <- data[[1]]
      c(z_mean, z_log_var) %<-% self$encoder(data)
      z <- self$sampler(z_mean, z_log_var)
      
      reconstruction <- decoder(z)
      # reconstruction_loss <- tf$reduce_mean(tf$square(data - reconstruction))
      reconstruction_loss <- loss_binary_crossentropy(data, reconstruction) %>%
        sum(axis = c(1)) %>%
        mean()
      # true_samples <- tf$random$normal(shape = tf$shape(z_mean))
      kl_loss <- -0.5 * (1 + z_log_var - z_mean^2 - exp(z_log_var))
      total_loss <- reconstruction_loss + lambda*kl_loss
      
      self$total_loss_tracker$update_state(total_loss)
      self$reconstruction_loss_tracker$update_state(reconstruction_loss)
      self$kl_loss_tracker$update_state(kl_loss)
      
      list(
        total_loss = self$total_loss_tracker$result(),
        reconstruction_loss = self$reconstruction_loss_tracker$result(),
        kl_loss = self$kl_loss_tracker$result()
      )
    }
  )
  
  
  vae <- model_vae(encoder, decoder)
  # Define learning rate schedule
  lr_schedule <- learning_rate_schedule_exponential_decay(
    5e-4,
    decay_steps = 100000, 
    decay_rate = 0.95, 
    staircase = FALSE
  ) 
  
  # Early stopping callback
  es_cb <- callback_early_stopping(
    min_delta = 1e-04, 
    monitor = 'val_total_loss',
    mode = 'min',
    patience = 25, 
    verbose = 1, 
    restore_best_weights = TRUE
  )
  
  # Compile the model
  opt <- optimizer_adam(learning_rate = lr_schedule,weight_decay = 0.1)
  
  vae %>% compile(optimizer = opt)
  
  # Fit the model, with validation data
  vae %>% fit(x_train, 
              x_train,  # For VAE, input is the same as output
              batch_size = batch_size,
              epochs = epochs,
              validation_data = list(x_val, x_val),
              shuffle = TRUE,
              callbacks = list(es_cb))
  
  return(list(vae = vae, encoder = encoder))
}


# Function to decode new samples
decode_samples_kl <- function(new_samples, vae) {
  tensorflow::set_random_seed(seed = 1994)
  encoded <- predict(vae$encoder, new_samples)
  
  z_mean <- encoded[[1]]
  z_log_var <- encoded[[2]]
  
  epsilon <- tf$random$normal(shape = tf$shape(z_mean))
  z <- z_mean + exp(0.5 * z_log_var) * epsilon
  
  # Step 3: Decode the Latent Representation
  decoded_samples <- predict(vae$decoder, z_mean) %>%
    as.data.frame()
  return(list(decoded = decoded_samples, encoded = z_mean))
}