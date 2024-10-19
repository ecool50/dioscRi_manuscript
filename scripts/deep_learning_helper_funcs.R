fitCelltypeModel <- function(train_x, train_y){
  # set.seed(1994)
  nnet_model <- SKM::deep_learning(
    y = train_y,
    x =  as.matrix(train_x),
    epochs_number = 100,
    layers = list(list(neurons_number = 27, 
                       neurons_proportion = NULL, 
                       activation = "relu")),
    optimizer = 'adamax', seed = 1994, 
    loss_function = 'categorical_crossentropy', 
    batch_size = 16, learning_rate = 0.005,
    verbose = TRUE
  )
  
  return(nnet_model)
}

computeConcordance <- function(truth, predicted){
  ARI = aricode::ARI(truth, predicted)
  NMI = aricode::NMI(truth, predicted)
  
  return(list(ARI = ARI, NMI = NMI))
}
