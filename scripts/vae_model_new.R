mmd_penalty <- function(pz, qz, batch_size = 32L, 
                        sigma_z = 1., z_dim = 16L){
  
  # This method calculates the unbiased U-statistic estimator of
  #         the MMD with the IMQ kernel. It's taken from
  #         https://github.com/tolstikhin/wae/blob/master/wae.py#L233
  # 
  #         Here the property that the sum of positive definite kernels is 
  #         still a p.d. kernel is used. Various kernels calculated at different
  #         scales are summed together in order to "simultaneously look at various
  #   scales [https://github.com/tolstikhin/wae/issues/2].
  
  
  norms_pz <- tf$reduce_sum(tf$square(pz), axis = 1L, keepdims=TRUE)
  dotprods_pz <- tf$matmul(pz,qz, transpose_b=TRUE)
  distances_pz <- norms_pz + tf$transpose(norms_pz) - 2. * dotprods_pz
  
  norms_qz <- tf$reduce_sum(tf$square(qz), axis = 1L, keepdims=TRUE)
  dotprods_qz <- tf$matmul(qz,qz, transpose_b=TRUE)
  distances_qz <- norms_qz + tf$transpose(norms_qz) - 2. * dotprods_qz
  
  dotprods <- tf$matmul(qz, pz, transpose_b=TRUE)
  distances <- norms_qz + tf$transpose(norms_pz) - 2. * dotprods
  
  cbase <- tf$constant(2. * z_dim * sigma_z)
  stat = tf$constant(0.)
  nf = tf$cast(batch_size, dtype=tf$float32)
  
  scales <- c(0.1, 0.2, 0.5, 1., 2., 5., 10.)
  for (scale in scales) {
    C <- cbase * scale
    res1 <- C / (C + distances_qz)
    res1 <- res1 + (C / (C + distances_pz))
    res1 <- tf$multiply(res1, 1. - tf$eye(batch_size, dtype = tf$float32))
    res1 <- tf$reduce_sum(res1) / (nf * nf - nf)
    
    res2 <- C / (C + distances)
    res2 <- tf$reduce_sum(res2) * 2. / (nf * nf)
    
    stat <- stat + (res1 - res2)
  }
  
  stat
}

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

robust_kernel <- function(alpha, sigma = 1) {
  tf_2pi = tf$constant(sqrt(2*pi), dtype=tf$float32)
  tf$exp(tf_2pi*(-(alpha * alpha) / (2 * sigma * sigma)))
}

correntropy_loss <- function(y_true, y_pred, sigma = 1) {
  loss_corr <- -tf$reduce_mean(robust_kernel(y_pred - y_true, sigma = sigma))
  loss_corr
}


# Custom loss function combining reconstruction and MMD loss
vae_mmd_loss <- function(y_true, y_pred, encoder, 
                         lambda_mmd = 1, sigma = 1.0) {
  # Reconstruction loss (MSE)
  recon_loss <- loss_mean_squared_error(y_true, y_pred)
  
  # MMD loss
  z_real <- encoder(y_true)  # Real latent representation
  z_fake <- tf$random$normal(shape = tf$shape(z_real)) * sigma  # Fake latent samples
  latent_dim <- dim(z_real)[2]
  mmd_loss <- compute_mmd(z_real, z_fake)
  # Total loss
  recon_loss + (lambda_mmd * mmd_loss)
}

# Define the custom loss function
custom_loss <- function(train_z, train_xr, train_x) {
  # Get batch size and latent dimension
  batch_size <- tf$shape(train_z)[1]
  latent_dim <- dim(train_z)[2]
  
  # Sample from a standard normal distribution (random noise)
  true_samples <- tf$random$normal(shape = c(batch_size, latent_dim), 
                                   mean = 0.0, stddev = 1.0)
  
  # Compute MMD loss
  loss_mmd <- compute_mmd(true_samples, train_z)
  
  # Compute reconstruction loss (mean squared error)
  loss_nll <- tf$reduce_mean(tf$square(train_xr - train_x))
  
  # Combine the losses (MMD loss + reconstruction loss)
  loss <- loss_nll + loss_mmd
  
  return(loss)
}


train_vae_model <- function(train_dat, useMarkers, epochs = 80, 
                            latent_dim = 2L, lambda = 0.1,
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
    layer_dense(intermediate_dim_2, activation = "relu") %>%
    layer_dense(latent_dim, activation = "gelu")
  
  # z_mean  <- x %>% layer_dense(latent_dim, name = "z_mean")
  encoder <- keras_model(encoder_inputs, x,
                         name = "encoder")
  
  decoder_outputs <- x %>%
    layer_dense(intermediate_dim_2, activation = "relu") %>%
    layer_dense(intermediate_dim, activation = "relu") %>%
    layer_dense(original_dim, activation = "sigmoid")
  
  decoder <- keras_model(x, decoder_outputs,
                         name = "decoder")
  
  # Full model (Encoder -> Decoder)
  train_z <- encoder(encoder_inputs)
  train_xr = decoder(train_z)
  vae <- keras_model(inputs = encoder_inputs, 
                     outputs = train_xr, name = "VAE")
  
  loss <- custom_loss(train_z, train_xr, encoder_inputs)
  vae$add_loss(loss)
  
  # Compile the model with the custom loss function
  lr_schedule <- learning_rate_schedule_exponential_decay(
    1e-3,
    decay_steps = 100000, 
    decay_rate = 0.95, 
    staircase = FALSE) 
  
  es_cb <- callback_early_stopping(
    min_delta = 1e-04, 
    monitor = 'val_loss',
    mode = 'auto',
    patience = 15, 
    verbose = 1, 
    restore_best_weights = TRUE
  )
  
  # First stage training with Adam optimizer
  opt <- optimizer_rmsprop(learning_rate = lr_schedule,
                           momentum = 0.9, centered = TRUE)
  
  vae %>% compile(
    optimizer = opt
  )
  
  # Fit the model
  vae %>% fit(
    x_train, x_train,
    epochs = epochs,
    callbacks = es_cb,
    batch_size = batch_size,  # Make sure batch_size is used correctly
    validation_split = 0.1
  )
  
  return(list(encoder = encoder, decoder = decoder, vae = vae))
}


# Function to decode new samples
decode_samples <- function(new_samples, model) {
  tensorflow::set_random_seed(seed = 1994)
  z_mean <- predict(model$encoder, new_samples)
  
  # Step 3: Decode the Latent Representation
  decoded_samples <- predict(model$decoder, z_mean) %>%
    as.data.frame()
  return(list(decoded = decoded_samples, encoded = z_mean))
}