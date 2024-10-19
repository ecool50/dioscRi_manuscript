# MSE loss
MSE_loss <- function(Y, X) {
  ret <- tf$square(X - Y)
  ret <- tf$reduce_sum(ret)
  return(ret)
}

# SE loss
SE_loss <- function(Y, X) {
  ret <- tf$square(X - Y)
  ret <- tf$reduce_sum(ret)
  return(ret)
}

# Gaussian Cross-Entropy Loss
Gaussian_CE_loss <- function(Y, X, beta, sigma = 0.5) { 
  Dim  <- tf$shape(X)[2]
  const1 <- tf$cast(-((1 + beta) / beta), tf$float32)
  const2 <- tf$cast(1 / ((2 * pi * (sigma^2))^(beta * Dim / 2)), tf$float32)
  SE <- SE_loss(Y, X)
  term1 <- tf$cast(tf$exp(-(beta / (2 * (sigma^2))) * SE),tf$float32)
  loss <- tf$reduce_sum(const1 * (const2 * term1 - 1))
  return(loss)
}

# Beta Loss Function
beta_loss_function <- function(recon_x, x, mu, logvar, beta = 0.01) {
  
  if (beta > 0) {
    # If beta is nonzero, use the beta entropy
    BBCE <- Gaussian_CE_loss(recon_x, x, beta)
  } else {
    # If beta is zero, use mean squared error
    BBCE <- MSE_loss(recon_x, x, beta)
  }
  
  # Compute KL divergence
  KLD <- -0.5 * tf$reduce_sum(1 + logvar - tf$square(mu) - tf$exp(logvar))
  
  return(BBCE + KLD)
}


train_rvae_model <- function(train_dat, useMarkers, epochs = 80, latent_dim = 2L, cell_types,
                            original_dim = 27L, batch_size = 16) {
  
  tensorflow::set_random_seed(seed = 1994)
  # Normalize and prepare the training data
  x_train <- train_dat[, useMarkers] %>%
    as.matrix()
  
  intermediate_dim <- original_dim - 4
  intermediate_dim_2 <- intermediate_dim - 4
  
  # Define the encoder
  
  encoder_inputs <- layer_input(shape = original_dim)
  
  x <- encoder_inputs %>%
    layer_dense(intermediate_dim, activation = "relu") %>%
    layer_dense(intermediate_dim_2, activation = "relu") 
  
  z_mean  <- x %>% layer_dense(latent_dim, name = "z_mean")
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
    layer_dense(intermediate_dim_2, activation = "relu") %>%
    layer_dense(intermediate_dim, activation = "relu") %>%
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
    },
    
    metrics = mark_active(function() {
      list(
        self$total_loss_tracker
      )
    }),
    
    train_step = function(data) {
      x <- data
      with(tf$GradientTape() %as% tape, {
        
        c(z_mean, z_log_var) %<-% self$encoder(x)
        z <- self$sampler(z_mean, z_log_var)
        
        
        reconstruction <- decoder(z)
        
        total_loss <- beta_loss_function(recon_x = reconstruction, x = x, mu = z_mean, logvar = z_log_var)
      })
      
      grads <- tape$gradient(total_loss, self$trainable_weights)
      self$optimizer$apply_gradients(zip_lists(grads, self$trainable_weights))
      
      self$total_loss_tracker$update_state(total_loss)
      
      list(total_loss = self$total_loss_tracker$result())
    }
  )
  
  
  vae <- model_vae(encoder, decoder)
  vae %>% compile(optimizer = optimizer_rmsprop())
  vae %>% keras3::fit(x_train,
                      epochs = epochs,
                      shuffle = TRUE)
  
  return(list(vae = vae, encoder = encoder))
}


# Function to decode new samples
decode_samples <- function(new_samples, vae, latent_dim = 8L, batch_size = 16) {
  tensorflow::set_random_seed(seed = 1994)
  encoded <- predict(vae$encoder, new_samples)
  z_mean <- encoded[[1]]
  z_log_var <- encoded[[2]]
  
  # Use the sampler to sample from the latent space
  epsilon <-  tf$random$normal(shape = c(nrow(z_mean), latent_dim))
  z <- z_mean + exp(0.5 * z_log_var) * epsilon
  
  # Step 3: Decode the Latent Representation
  decoded_samples <- predict(vae$decoder, z, batch_size = batch_size) %>%
    as.data.frame()
  return(list(decoded = decoded_samples, encoded = z_mean))
}