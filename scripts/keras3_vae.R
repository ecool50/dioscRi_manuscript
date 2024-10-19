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

original_dim <- 27L
latent_dim <- 2L
intermediate_dim <- 24L
intermediate_dim_2 <- 21
batch_size<- 128
n_cell_types <- length(unique(train_dat$mg_cell_type_distinct))

encoder_inputs <- layer_input(shape = original_dim)

x <- encoder_inputs %>%
  layer_dense(intermediate_dim, activation = "relu") %>%
  layer_dense(intermediate_dim_2, activation = "relu")

z_mean  <- x %>% layer_dense(latent_dim, name = "z_mean")
z_log_var <- x %>% layer_dense(latent_dim, name = "z_log_var")
encoder <- keras_model(encoder_inputs, list(z_mean, z_log_var),
                       name = "encoder")

encoder


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

decoder

classifier_outputs <- latent_inputs %>%
  layer_dense(units = n_cell_types, 
              activation = "softmax", name = "cell_type_pred", 
              kernel_regularizer = regularizer_l1_l2(0))

classifier <- keras_model(latent_inputs, classifier_outputs, name = "classifier")


model_vae <- new_model_class(
  classname = "VAE",
  
  initialize = function(encoder, decoder,classifier_outputs, ...) {
    super$initialize(...)
    self$encoder <- encoder
    self$decoder <- decoder
    self$sampler <- layer_sampler()
    self$classifier <- classifier
    self$total_loss_tracker <-
      metric_mean(name = "total_loss")
    self$reconstruction_loss_tracker <-
      metric_mean(name = "reconstruction_loss")
    self$kl_loss_tracker <-
      metric_mean(name = "kl_loss")
    self$mmd_loss_tracker <- 
      metric_mean(name = "mmd_loss")
    self$class_loss_tracker <- metric_mean(name = "class_loss")
  },
  
  metrics = mark_active(function() {
    list(
      self$total_loss_tracker,
      self$reconstruction_loss_tracker,
      self$kl_loss_tracker,
      self$mmd_loss_tracker,
      self$class_loss_tracker
    )
  }),
  
  train_step = function(data) {
    x <- data[[1]]
    y <- data[[2]]
    with(tf$GradientTape() %as% tape, {
      
      c(z_mean, z_log_var) %<-% self$encoder(x)
      z <- self$sampler(z_mean, z_log_var)
      
      
      reconstruction <- decoder(z)
      cell_type_pred <- self$classifier(z_mean)
      reconstruction_loss <-
        loss_binary_crossentropy(x, reconstruction) %>%
        sum(axis = c(1)) %>%
        mean()
      
      kl_loss <- -0.5 * (1 + z_log_var - z_mean^2 - exp(z_log_var))
      mmd_loss <- compute_mmd(z, z_mean)
      class_loss <- loss_categorical_focal_crossentropy(y, cell_type_pred)
      
      total_loss <- reconstruction_loss + mean(kl_loss) + mean(mmd_loss) + mean(class_loss)
    })
    
    grads <- tape$gradient(total_loss, self$trainable_weights)
    self$optimizer$apply_gradients(zip_lists(grads, self$trainable_weights))
    
    self$total_loss_tracker$update_state(total_loss)
    self$reconstruction_loss_tracker$update_state(reconstruction_loss)
    self$kl_loss_tracker$update_state(kl_loss)
    self$mmd_loss_tracker$update_state(mmd_loss)
    self$class_loss_tracker$update_state(class_loss)
    
    list(total_loss = self$total_loss_tracker$result(),
         reconstruction_loss = self$reconstruction_loss_tracker$result(),
         kl_loss = self$kl_loss_tracker$result(),
         mmd_loss = self$mmd_loss_tracker$result(),
         class_loss = self$class_loss_tracker$result())
  }
)

mnist <- dataset_mnist()

## normalize so the range is (0,1)
# x_train <- df_asinh[, useMarkers]/max(df_asinh[, useMarkers])
# x_test <- df_3[, useMarkers]/max(df_3[, useMarkers])
cell_types <- train_dat$mg_cell_type_distinct
cell_types <- to_categorical(as.numeric(cell_types %>% as.factor()) - 1)

vae <- model_vae(encoder, decoder, classifier_outputs)
vae %>% compile(optimizer = optimizer_adam())
vae %>% fit(train_dat[, useMarkers] %>% as.matrix(), 
            y = cell_types,
            epochs = 20,
            shuffle = TRUE)


library(ggplot2)
x_test_encoded <- predict(encoder, train_dat[, useMarkers] %>% as.matrix(), batch_size = batch_size)

x_test_encoded <- x_test_encoded[[1]] %>%
  as.data.frame() %>%
  dplyr::mutate(clusters = train_dat$mg_cell_type_distinct)

x_test_encoded %>%
  mutate(class = as.factor(clusters)) %>%
  ggplot(aes(x = V1, y = V2, colour = class)) + 
  geom_point() +
  scattermore::geom_scattermore()+
  theme_classic(base_size = 14) 


# cl <- makePSOCKcluster(6)
# registerDoParallel(cl)
# 
# fit.control <- trainControl(method = "repeatedcv", number = 5, repeats = 3, allowParallel = TRUE)
# 
# set.seed(1994)  
# fit <- train(clusters ~ ., data = x_test_encoded, method = "multinom", trControl = fit.control, trace = TRUE)
# 
# stopCluster(cl)
# fit

