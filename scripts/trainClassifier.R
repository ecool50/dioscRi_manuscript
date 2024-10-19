 trainClassifier <- function(train, layer_sizes = c(90, 90, 45, 45),
                    activation = 'tanh', l2_penalty = 1e-4,
                    weights_path = 'None', train_label){
   
   nClusts <- length(unique(train_label))
   input_layer <- layer_input(shape = ncol(train))
   hidden_1 <- layer_dense(object = input_layer, units = layer_sizes[[1]],
                           activation = activation, regularizer_l2(l = l2_penalty)) 
   hidden_2 <- layer_dense(object = hidden_1, units = layer_sizes[[2]],
                           activation = activation, regularizer_l2(l = l2_penalty))
   hidden_3 <- layer_dense(object = hidden_2, units = layer_sizes[[3]],
                           activation = activation, regularizer_l2(l = l2_penalty))
   hidden_4 <- layer_dense(object = hidden_3, units = layer_sizes[[4]],
                           activation = activation, regularizer_l2(l = l2_penalty))
   output_layer <- layer_dense(object = hidden_4, units = nClusts+1, 
                               activation = 'softmax')

   # encoder <- keras_model(inputs = input_layer, outputs = output_layer)
   model <- keras_model(inputs = input_layer, outputs = output_layer)
   
   lr_schedule <- learning_rate_schedule_exponential_decay(
     1e-3,
     decay_steps = 100000, 
     decay_rate = 0.95, 
     staircase = FALSE) 
   
   opt <- optimizer_rmsprop(learning_rate = lr_schedule,
                            momentum = 0.9, centered = TRUE)
   
   model %>% 
     compile(optimizer = opt, 
             loss = 'categorical_crossentropy',
             metrics = keras3::metric_auc(multi_label = TRUE))
   
   es_cb <- keras3::callback_early_stopping(patience = 15, monitor = 'val_loss', 
                                           mode = 'auto')
   # checkpoint <- keras::callback_model_checkpoint(filepath = paste0(weights_path,'.hdf5'), 
   #                                                save_best_only = T)
   
   model %>%
     fit(x = as.matrix(train),
         y = to_categorical(as.numeric(train_label)), 
         epochs = 10,
         batch_size = 32, # try also 128 and 256
         validation_split = 0.1,
         callbacks = es_cb)
   
   return(model)
 }
