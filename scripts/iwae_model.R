train_iwae_model <- function(train_dat, useMarkers, epochs = 80, 
                             latent_dim = 2L, lambda = 0.1, val_dat,
                             original_dim = 27L, batch_size = 16, K = 5) {
  
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
    layer_dense(intermediate_dim, activation = "relu") %>%
    layer_dense(intermediate_dim_2, activation = "relu")
  
  z_mean    <- x %>% layer_dense(latent_dim, name = "z_mean")
  z_log_var <- x %>% layer_dense(latent_dim, name = "z_log_var")
  encoder <- keras_model(encoder_inputs, list(z_mean, z_log_var),
                         name = "encoder")
  
  layer_sampler <- new_layer_class(
    classname = "Sampler",
    call = function(z_mean, z_log_var) {
      epsilon <- tf$random$normal(shape = tf$shape(z_mean))
      z_mean + exp(0.5 * z_log_var) * epsilon
    }
  )
  
  # Define the decoder
  latent_inputs <- layer_input(shape = c(latent_dim))
  
  decoder_outputs <- latent_inputs %>%
    layer_dense(intermediate_dim_2, activation = "relu") %>%
    layer_dense(intermediate_dim, activation = "relu") %>%
    layer_dense(original_dim, activation = "sigmoid")
  
  decoder <- keras_model(latent_inputs, decoder_outputs, name = "decoder")
  
  # IWAE Loss Function
  model_iwae <- new_model_class(
    classname = "IWAE",
    
    initialize = function(encoder, decoder, K, ...) {
      super$initialize(...)
      self$encoder <- encoder
      self$decoder <- decoder
      self$sampler <- layer_sampler()
      self$total_loss_tracker <- metric_mean(name = "total_loss")
      self$reconstruction_loss_tracker <- metric_mean(name = "reconstruction_loss")
      self$iwae_loss_tracker <- metric_mean(name = "iwae_loss")
      self$K <- K  # Number of importance samples
    },
    
    metrics = mark_active(function() {
      list(
        self$total_loss_tracker,
        self$reconstruction_loss_tracker,
        self$iwae_loss_tracker
      )
    }),
    
    # Training Step
    train_step = function(data) {
      with(tf$GradientTape() %as% tape, {
        data <- data[[1]]
        c(z_mean, z_log_var) %<-% self$encoder(data)
        
        # Sampling K latents
        total_loss <- 0
        log_p_x_given_z_all <- 0
        
        for (i in seq_len(self$K)) {
          z <- self$sampler(z_mean, z_log_var)
          reconstruction <- decoder(z)
          
          reconstruction_loss <- loss_binary_crossentropy(data, reconstruction) %>%
            sum(axis = c(1)) %>%
            mean()
          
          # KL Divergence
          kl_loss <- -0.5 * (1 + z_log_var - z_mean^2 - exp(z_log_var))
          total_loss <- total_loss + reconstruction_loss + lambda * kl_loss
        }
        
        # Averaging over K samples
        total_loss <- total_loss / self$K
      })
      
      grads <- tape$gradient(total_loss, self$trainable_weights)
      self$optimizer$apply_gradients(zip_lists(grads, self$trainable_weights))
      
      self$total_loss_tracker$update_state(total_loss)
      self$reconstruction_loss_tracker$update_state(reconstruction_loss)
      self$iwae_loss_tracker$update_state(total_loss)
      
      list(total_loss = self$total_loss_tracker$result(),
           reconstruction_loss = self$reconstruction_loss_tracker$result(),
           iwae_loss = self$iwae_loss_tracker$result())
    },
    
    # Validation Step
    test_step = function(data) {
      data <- data[[1]]
      c(z_mean, z_log_var) %<-% self$encoder(data)
      
      total_loss <- 0
      
      for (i in seq_len(self$K)) {
        z <- self$sampler(z_mean, z_log_var)
        reconstruction <- decoder(z)
        
        reconstruction_loss <- loss_binary_crossentropy(data, reconstruction) %>%
          sum(axis = c(1)) %>%
          mean()
        
        kl_loss <- -0.5 * (1 + z_log_var - z_mean^2 - exp(z_log_var))
        total_loss <- total_loss + reconstruction_loss + lambda * kl_loss
      }
      
      total_loss <- total_loss / self$K
      
      self$total_loss_tracker$update_state(total_loss)
      self$reconstruction_loss_tracker$update_state(reconstruction_loss)
      self$iwae_loss_tracker$update_state(total_loss)
      
      list(
        total_loss = self$total_loss_tracker$result(),
        reconstruction_loss = self$reconstruction_loss_tracker$result(),
        iwae_loss = self$iwae_loss_tracker$result()
      )
    }
  )
  
  # Instantiate the IWAE model
  iwae <- model_iwae(encoder, decoder, K)
  
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
    patience = 30, 
    verbose = 1, 
    restore_best_weights = TRUE
  )
  
  # Compile the model
  opt <- optimizer_rmsprop(learning_rate = lr_schedule,
                           momentum = 0.9, centered = TRUE)
  
  iwae %>% compile(optimizer = opt)
  
  # Fit the model, with validation data
  iwae %>% fit(x_train, 
               x_train,  # For IWAE, input is the same as output
               batch_size = batch_size,
               epochs = epochs,
               validation_data = list(x_val, x_val),
               shuffle = TRUE,
               callbacks = list(es_cb))
  
  return(list(iwae = iwae, encoder = encoder))
}

# Function to decode new samples
decode_samples_iwae <- function(new_samples, vae) {
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