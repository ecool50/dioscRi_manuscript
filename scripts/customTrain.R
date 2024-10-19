customTrain <- function(train_x, test_x, ratio = 0.75,
                        nFolds = 5, nRepeats = 5, model = 'lgr'){
  
  # generate seed values
  seeds <- sample(1:1e6, nRepeats, replace=FALSE)
  
  # generate model parameters
  h2o.init(max_mem_size = "20g")  # connect to H2O instance
  # train <- train_x %>% as.data.frame() 
  # train <- as.h2o(train)
  
  # test <- test_x
  # test <- as.h2o(test)
  
  y <- 'gensini_bin'
  x <- setdiff(names(data), y)
  
  train_res <- c()
  test_res <- c()
  
  
  for(i in 1:length(seeds)){
    # data_split <- h2o.splitFrame(data = as.h2o(data), ratios = ratio, 
    #                              seed = seeds[[i]])
    
    if(model == 'gbm'){
      my_model <- h2o.gbm(x = x,
                        y = y,
                        training_frame = data_split[[1]],
                        distribution = "bernoulli",
                        ntrees = 50,
                        max_depth = 3,
                        min_rows = 2,
                        learn_rate = 0.2,
                        nfolds = nFolds,
                        keep_cross_validation_predictions = F,
                        seed = seeds[[i]])
    }
    
    
    if(model == 'lgr'){
      my_model <- h2o.glm(x = x, # Vector of predictor variable names
                          y = y, # Name of response/dependent variable
                          training_frame = train_x, # Training data
                          seed = seeds[[i]],        # Seed for random numbers
                          family = "binomial",   # Outcome variable
                          lambda_search = TRUE,  # Optimum regularisation lambda
                          alpha = 0,           # Elastic net regularisation
                          nfolds = nFolds, 
                          max_iterations = 1000,
                          balance_classes = T,# N-fold cross validation
      )
    }
    if(model == 'deep'){
      # Build and train the model:
      my_model <- h2o.deeplearning(x = x,
                                   y = y,
                                   distribution = "bernoulli",
                                   hidden = c(1),
                                   reproducible = FALSE,
                                   activation = "tanh",
                                   single_node_mode = FALSE,
                                   balance_classes = T,
                                   force_load_balance = FALSE,
                                   seed = seeds[[i]],
                                   score_training_samples = 0,
                                   score_validation_samples = 0,
                                   training_frame = train_x,
                                   stopping_rounds = 0, 
                                   nfolds = nFolds)
    }
    
    
    train_res[[i]] <- h2o.auc(my_model, xval = TRUE)
    test_res[[i]] <- h2o.performance(my_model, newdata = test_x) %>% h2o.auc()
    
    # perf <- h2o.performance(my_model, data_split[[2]])
    # plot(perf, type = "roc")
  }
  
  return(list(train = mean(unlist(train_res)), 
              test = mean(unlist(test_res))))
}