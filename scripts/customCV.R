KFoldCustomOSCC <- function(features, k = 5, nRepeats = 1, method = "svm",
                            metaData, nFeatures = 2, response, sampleCol){
  # tensorflow::set_random_seed(1994)
  
  res <- c()
  res_test <- c()
  truth_list <- c()
  predictions_list <- c()
  cur_best <- 0
  cur_best_features <- c()
  best_overall <- c()
  
  count = 1
  features_new <- features
  features_new$sample_id <- rownames(features) %>% as.numeric()
  features_new <- merge(features_new, metaData)
  # seeds <- sample(1:1e6, nRepeats, replace=FALSE)
  for(i in 1:nRepeats){
    # tensorflow::set_random_seed(seed = seeds[[i]])
    res_cur <- c()
    res_cur_test <- c()
    cv <- SKM::cv_kfold_strata(features_new[, response], k = k)
    
    for(j in 1:k){
      idx <- cv[[j]]$training
      idx_val <- cv[[j]]$testing
      cur_x_train <- features_new[idx, colnames(features)]
      
      
      cur_y_train <- features_new[idx, response]
      testRes <- colTest(cur_x_train,
                         condition = cur_y_train, type = 'ttest')
      
      # sigFeatures <- nestedcv::glmnet_filter(y = cur_y_train, x = cur_x_train,
      #                                         nfilter = NULL,
      #                                       type = "names")
      sigFeatures <- testRes$cluster
      cur_x_train <- cur_x_train[, sigFeatures]
      colnames(cur_x_train) <- janitor::make_clean_names(colnames(cur_x_train))
      
      cur_x_val <- features_new[idx_val, sigFeatures]
      colnames(cur_x_val) <- janitor::make_clean_names(cur_x_val)

      cur_y_val <- features_new[idx_val, response]
      
      if(method == 'glm'){
        model <-  glmnet::cv.glmnet(as.matrix(cur_x_train), 
                                    cur_y_train, family = 'binomial',
                                    alpha = 0.5)
        t <- rownames(coef(model, s = 'lambda.min'))[coef(model, s = 'lambda.min')[,1]!= 0]
        t <- t[t != '(Intercept)']
        
      } else if(method == 'rf') {
        model <-  SKM::random_forest(
          as.matrix(cur_x_train),
          factor(cur_y_train),
          verbose = F,
          seed = 1994
        )
        
        t <- sigFeatures
      } else if(method == "svm"){
        model <-  SKM::support_vector_machine(
          as.matrix(cur_x_train),
          factor(cur_y_train),
          verbose = F,
          seed = 1994)
        t <- sigFeatures
        
      } else {
        model <- getModel(data = cur_x_train)
        es_cb <- callback_early_stopping(patience = 30, monitor = 'val_loss')
        
        model %>% fit(
          as.matrix(cur_x_train), 
          to_categorical(as.numeric(cur_y_train) - 1),
          epochs = 100,
          batch_size = 16,
          validation_split = 0.1,
          callbacks = es_cb,
          verbose = 0
        )
        
        t <- sigFeatures
      }
      
      # make predictions
      if(method == 'glm'){
        predictions <- predict(model, as.matrix(cur_x_val), type = 'response',
                               s = 'lambda.min')
        predictions_list[[count]] <- predictions[, 1]
        truth_list[[count]] <- as.numeric(cur_y_val)
        
        roc_pred_val <- ROCR::prediction(predictions[, 1],
                                         cur_y_val)
        
        auc_perf <- ROCR::performance(roc_pred_val, measure = "auc")
        auc_j <- auc_perf@y.values[[1]]
      } else if(method == 'nn'){
        auc_j <- evaluate(model, as.matrix(cur_x_val), 
                                            to_categorical(as.numeric(cur_y_val) - 1), 
                                            batch_size = 16, verbose = 0, 
                                            sample_weight = NULL)[[2]]
      }
      else{
        predictions <- predict(model, as.matrix(cur_x_val))
        predictions_list[[count]] <- predictions$probabilities[, 2]
        truth_list[[count]] <- as.numeric(cur_y_val)
        
        roc_pred_val <- ROCR::prediction(predictions$probabilities[, 2],
                                         cur_y_val)
        
        auc_perf <- ROCR::performance(roc_pred_val, measure = "auc")
        auc_j <- auc_perf@y.values[[1]]
      }
      
      res_cur[[j]] <- auc_j
      
      if(auc_j > cur_best){
        cur_best <- auc_j
        cur_best_features <- t
        best_overall <- colnames(cur_x_train)
        message(paste('The best AUC score is: ', round(cur_best, 2)))
      }
    }
    res[[i]] <- mean(unlist(res_cur))
    
  }
  message(paste('The mean CV AUC is:', round(mean(unlist(res)), 2)))
  return(cur_best_features)
}


runKFeatures <- function(data, featureRange = c(2,5,10,15,20,25), clinicalData,
                         sampleCol, response, k = 5, nRepeats = 20, method = "svm"){
  
  # define a list to hold the plots
  myplots <- vector('list', length = length(featureRange)*2)
  
  for (i in length(featureRange)) {
    # generate top n features
    nFeatures <-featureRange[[i]]
    t <- KFoldCustomOSCC(features = data, k = k, nRepeats = nRepeats,
                         metaData = clinicalData,
                         nFeatures = nFeatures, sampleCol = sampleCol,
                         response = response, method = method)
    # t <- sort(table(t), decreasing=TRUE)[1:nFeatures] %>% names()
    
  #   # Generate classification data
  #   x <- rownames(data)
  #   colnames(data) <- janitor::make_clean_names(colnames(data))
  #   
  #   data_final <- data[, t]
  #   colnames(data_final) <- janitor::make_clean_names(colnames(data_final))
  #   
  #   data_final$response <- factor(clinicalData[match(x,clinicalData[, sampleCol]), 
  #                                              response])
  #   data_final <- na.omit(data_final)
  #   
  #   # build classification model
  #   n <- ncol(data_final)
  #   train_x <- data_final[, 1:n-1] %>%
  #     as.matrix()
  #   train_y <- data_final$response
  #   # train_y <- factor(if_else(train_y == 'N', 0, 1))
  #   
  #   model_svm <- support_vector_machine(
  #     as.matrix(train_x),
  #     factor(train_y),
  #     kernel = "linear",
  #     verbose = F,
  #     seed = 1994)
  #   
  #   # do prediction
  #   predictions <- predict(model_svm, as.matrix(train_x))
  #   
  #   roc.pred.test <- ROCR::prediction(predictions$probabilities[, 2], 
  #                                     train_y)
  #   perf <- ROCR::performance(roc.pred.test, "tpr", "fpr")
  #   
  #   
  #   plt_dat = data.frame(
  #     FPR = perf@x.values[[1]],
  #     TPR = perf@y.values[[1]]
  #   )
  #   
  #   auc_perf <- ROCR::performance(roc.pred.test, measure = "auc")
  #   auc <- auc_perf@y.values[[1]]
  #   
  #   # generate AUC plot
  #   p.auc <- ggplot(plt_dat, aes(x = FPR, y = TPR)) +
  #     geom_line(colour = "blue") +
  #     labs(x = perf@x.name, y = perf@y.name) +
  #     geom_abline(slope = 1, intercept = 0) + theme_bw() + 
  #     theme(
  #       plot.title = element_text(color="Black", size=16, face="bold", hjust = 0.5),
  #       plot.subtitle = element_text(color = "red", size = 16, hjust = 0.5),
  #       axis.title.x = element_text(color="Black", size=16),
  #       axis.title.y = element_text(color="Black", size=16),
  #       axis.text.y = element_text(size = 16),
  #       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 16),
  #       strip.text.x = element_text(size = 16),
  #       legend.title = element_text(size=16), #change legend title font size
  #       legend.text = element_text(size=16) #change legend text font size)
  #     )  + ggtitle(paste("AUC for top ", nFeatures, 'features: ' ,
  #                        round(auc, digits = 2)))
  #   # generate model weights
  #   sv <- model_svm[["fitted_model"]][["SV"]]
  #   coef <- model_svm[["fitted_model"]][["coefs"]]
  #   
  #   w = t(sv) %*% coef %>%
  #     as.data.frame()
  #   b = -model_svm[["fitted_model"]][["rho"]]
  #   
  #   w$Feature <- rownames(w)
  #   colnames(w)[[1]] <- 'Weight'
  #   
  #   # plot model weights
  #   w$Feature <- gsub('fox_p3', 'foxp3', w$Feature)
  #   w$Feature <- gsub('pan_ck', 'panck', w$Feature)
  #   w$Feature <- stri_replace_last(w$Feature, fixed = "_", ": ") 
  #   w$Feature <- gsub('_', ' ', w$Feature)
  #   
  #   p.weights <- ggplot(w, aes(x=reorder(Feature, Weight), y=Weight, fill=Weight))+
  #     geom_bar(stat="identity") +
  #     scale_fill_gradient2(low="darkblue", high="darkred") +
  #     ggtitle("Plot of feature weights for SVM model")+
  #     labs(x = "Feature",
  #          y = "Weight") +
  #     theme(axis.title = element_text(size = 15, face = "bold"),
  #           axis.text = element_text(size = 15),
  #           plot.title = element_text(size = 16, hjust = 0.5, face = "bold")) + 
  #     coord_flip()
  #   
  #   # myplots[[2*i]] <- p.weights
  #   p <- cowplot::plot_grid(p.auc, p.weights)
  }
  
  return(t)
}

kFoldCV <- function(train_x, train_y, p = 0.8, k = 5, nRepeats = 1){
  # set.seed(1994)
  
  seeds <- sample(1:1e6, nRepeats, replace=FALSE)
  
  res <- c()
  for(i in 1:nRepeats){
    set.seed(seeds[[i]])
    res_cur <- c()
    # cv <- modelr::crossv_kfold(train_x, k = k, )
    cv <- SKM::cv_kfold_strata(train_y, k = k)
    
    message(paste('Now running iteration:', i))
    for(j in 1:k){
      message('Now running for fold: ', j)
      # idx <- cv$train[[j]]$idx
      # idx_val <- cv$test[[j]]$idx
      idx <- cv[[j]]$training
      idx_val <- cv[[j]]$testing
      cur_x_train <- train_x[idx, ]
      cur_y_train <- train_y[idx]
      
      cur_x_val <- train_x[idx_val, ]
      cur_y_val <- train_y[idx_val]
      
      # fit the model 
      res_j <- computePred(train_x = cur_x_train, train_y = cur_y_train,
                           test_x = cur_x_val, test_y = cur_y_val)
      
      res_cur[[j]] <- res_j
    }
    res[[i]] <- mean(unlist(res_cur))
  }
  
  
  return(unlist(res))
}

KFoldCustom <- function(train_x, train_y, p = 0.8, k = 5, nRepeats = 1,
                        method = 'glm', alpha = 0.5){
  set.seed(1994)
  seeds <- sample(1:1e6, nRepeats, replace=FALSE)
  
  res <- c()
  truth_list <- c()
  predictions_list <- c()
  count = 1
  for(i in 1:nRepeats){
    set.seed(seeds[[i]])
    res_cur <- c()
    # cv <- modelr::crossv_kfold(train_x, k = k, )
    cv <- SKM::cv_kfold_strata(train_y, k = k)
    
    # message(paste('Now running iteration:', i))
    for(j in 1:k){
      idx <- cv[[j]]$training
      idx_val <- cv[[j]]$testing
      cur_x_train <- train_x[idx, ]
      cur_y_train <- train_y[idx]
      
      cur_x_val <- train_x[idx_val, ]
      cur_y_val <- train_y[idx_val]
      
      if(method == 'glm'){
        model <-  glmnet::cv.glmnet(as.matrix(cur_x_train), 
                                    cur_y_train, family = 'binomial',
                                    alpha = alpha, nfolds = 5, type.measure = 'auc')
      } else if(method == 'rf') {
        model <-  SKM::random_forest(
          as.matrix(cur_x_train),
          factor(cur_y_train),
          verbose = F,
          seed = 1994
        )
      } else {
        model <-  SKM::support_vector_machine(
          as.matrix(cur_x_train),
          factor(cur_y_train),
          verbose = F,
          seed = 1994)
      }
      
      # make predictions
      if(method == 'glm'){
        predictions <- predict(model, as.matrix(cur_x_val), type = 'response', 
                               s = 'lambda.min')
        predictions_list[[count]] <- predictions[, 1]
        truth_list[[count]] <- as.numeric(cur_y_val)
        
        roc_pred_val <- ROCR::prediction(predictions[, 1],
                                         cur_y_val)
      } else{
        predictions <- predict(model, as.matrix(cur_x_val))
        predictions_list[[count]] <- predictions$probabilities[, 2]
        truth_list[[count]] <- as.numeric(cur_y_val)
        
        roc_pred_val <- ROCR::prediction(predictions$probabilities[, 2],
                                         cur_y_val)
      }
      auc_perf <- ROCR::performance(roc_pred_val, measure = "auc")
      auc_j <- auc_perf@y.values[[1]]
      
      # res[[count]] <- auc_j
      res_cur[[j]] <- auc_j
      
      count = count + 1
    }
    res[[i]] <- round(mean(unlist(res_cur)), 2)
  }
  
  
  message(paste('The mean CV AUC is:', round(mean(unlist(res)), 2)))
  cv_dat <- list(predictions = predictions_list,
                 labels = truth_list)
  
  return(list(cv_res = unlist(res), cv_dat = cv_dat))
}


KFoldSvmMaj <- function(train_x, train_y, p = 0.8, k = 5, nRepeats = 1,
                        test_x, test_y, model_weights){
  seeds <- sample(1:1e6, nRepeats, replace=FALSE)
  
  res <- c()
  res_test <- c()
  for(i in 1:nRepeats){
    set.seed(seeds[[i]])
    res_cur <- c()
    res_cur_test <- c()
    # cv <- modelr::crossv_kfold(train_x, k = k, )
    cv <- SKM::cv_kfold_strata(train_y, k = k)
    
    message(paste('Now running iteration:', i))
    for(j in 1:k){
      message('Now running for fold: ', j)
      # idx <- cv$train[[j]]$idx
      # idx_val <- cv$test[[j]]$idx
      idx <- cv[[j]]$training
      idx_val <- cv[[j]]$testing
      cur_x_train <- train_x[idx, ]
      cur_y_train <- train_y[idx]
      
      cur_x_val <- train_x[idx_val, ]
      cur_y_val <- train_y[idx_val]
      
      optim_lambda <- svmmajcrossval(
        cur_x_train, cur_y_train,
        mc.cores = 30,
        ngroup = 10,
        return.model = TRUE,
        verbose = T,
        weights.obs = model_weights[idx],
        spline.knots = 8, spline.degree = 3
        
      )$param.opt[[1]]
      
      model <- svmmaj(
        cur_x_train, cur_y_train, spline.knots = 8, spline.degree = 3,
        lambda = optim_lambda, weights.obs = model_weights[idx])
      
      # make predictions
      predictions <- predict(model, cur_x_val, cur_y_val, show.plot = F)
      auc_j <- SVMMaj:::auc(predictions, cur_y_val)
      
      res_cur[[j]] <- auc_j
      
    }
    res[[i]] <- mean(unlist(res_cur))
  }
  
  predictions_test <- predict(model, test_x, test_y, show.plot = F)
  
  auc_test <- SVMMaj:::auc(predictions_test, test_y)
  
  
  
  message(paste('The mean CV AUC is:', round(mean(unlist(res)), 2)))
  message(paste("The test AUC is:", round(auc_test,2)))
  
  
  return(unlist(res))
}


mean_normalize <- function(x) {
  (x - mean(x)) / (max(x) - min(x))
}


kFOldCustomNested <- function(features, nFolds = 5, clinicalData, sampleCol, 
                              response, nRepeats = 10, featureRange = 10,
                              method = 'nn'){
  res <- rep(NA, nFolds*nRepeats)
  top_features <- c()
  count <- 1
  tensorflow::set_random_seed(1994)
  seeds <- sample(1:1e6, nRepeats, replace=FALSE)
  for(i in 1:nRepeats){
    set.seed(seeds[[i]])
    for(i in 1:nFolds){
      cv <- SKM::cv_kfold_strata(clinicalData[, response], k = nFolds)
      idx <- cv[[i]]$training
      idx_val <- cv[[i]]$testing
      cur_x_train <- features[idx, ]
      cur_y_train <- clinicalData[, response][idx]
      cur_x_val <- features[idx_val, ]
      cur_y_val <- clinicalData[, response][idx_val]
      clinicaldata_cur <- clinicalData[idx, ]
      
      # print(table(cur_y_train))
      
      cur_y_train_num <- as.numeric(levels(cur_y_train))[cur_y_train]
      # cur_x_train <- cbind(cur_x_train, cur_y_train_num)
      
      # res_efs <- ensemble_fs(data = cur_x_train, 
      #                        classnumber = ncol(cur_x_train), 
      #                        runs = 10, cor_threshold = 0.75,
      #                        selection = c(FALSE, FALSE, TRUE, FALSE, FALSE, 
      #                                     FALSE, FALSE,
      #                                     FALSE))
      # 
      # res_efs <- t(res_efs) %>% as.data.frame()
      # res_efs$feature <- rownames(res_efs)
      # # order based on correlaton
      # res_efs <- res_efs[order(res_efs$S_cor, decreasing = TRUE), ]
      # 
      # cur_features <- res_efs$feature[1:featureRange]
      
      # cur_features <- nestedcv::correl_filter(y = as.numeric(cur_y_train), 
      #                                         x = as.matrix(cur_x_train),
      #                                         type = "names")
      
      # res_g <- feature.selection(as.matrix(cur_x_train), cur_y_train_num,
      #                            d = featureRange)
      
      # idx <- res_g$integrated_selected_feature_index
      # dat <- cur_x_train %>%
      #   min_max_scaling()
      # cur_y_train_num <- ifelse(cur_y_train == 0, -1, 1)
      # cur_y_val_num <- ifelse(cur_y_val == 0, -1, 1)
      # 
      # hc <- bootstrapHclust(dat,frac=1,B=100, nCore = 7)
      # groupWeight <- computeGroupSizeWeight(hc)
      # 
      # res_mgml <- MLGL(as.matrix(dat), cur_y_train_num, loss = 'logit',
      #                  verbose = TRUE, hc = hc, weightSizeGroup = groupWeight)
      # 
      # out <- HMT(res_mgml,as.matrix(dat), cur_y_train_num, control="FWER",alpha=0.1)
      # summary(out)
      
      # dat_new <- dat[, out$var] %>%
      #   as.data.frame()
      
      cur_features <- nestedcv::anova_filter(y = as.numeric(cur_y_train),
                                              x = as.matrix(cur_x_train),
                                              type = "names")
      # if(ncol(dat_new) < 2){
      #   next
      # }
      # message("Anova success")
      # 
      # cur_features_ttest <- nestedcv::ttest_filter(y = cur_y_train, x = cur_x_train,
      #                                              type = "names")
      # message("T-test success")
      
      # cur_features_wilcox <- nestedcv::wilcoxon_filter(y = cur_y_train, x = cur_x_train,
      #                                              type = "names", p_cutoff = 0.1)
      # message("Wilcoxon success")
      
      # cur_features_glmnet <- nestedcv::glmnet_filter(y = cur_y_train, x = cur_x_train,
      #                                                type = "names")
      # message("GLMNET success")
      # 
      # cur_features <- union(cur_features_glmnet, 
      #                           union(cur_features_anova, cur_features_ttest))
      
      # cur_features <- runKFeatures(cur_x_train, clinicalData = clinicaldata,
      #                              featureRange = featureRange,k = nFolds, nRepeats = 1,
      #                              sampleCol = sampleCol, response = response, method = 'glm')
      cur_features <- cur_features[!is.na(cur_features)]
      if(length(cur_features) == 0){
        next
      }
      # print(cur_features)
      cur_x_train <- cur_x_train[, cur_features]
      
      cur_x_val <- cur_x_val[, cur_features]
      # model <-  glmnet::cv.glmnet(as.matrix(cur_x_train[, cur_features]), 
      #                             cur_y_train, family = 'binomial',
      #                             alpha = 1, nfolds = 5, type.measure = 'auc')
      
      # predictions <- predict(model, as.matrix(cur_x_val[, cur_features]), type = 'response', 
      #                        s = 'lambda.min')
      # # predictions_list[[count]] <- predictions[, 1]
      # 
      # roc_pred_val <- ROCR::prediction(predictions[, 1],
      #                                  cur_y_val)
      # 
      # auc_perf <- ROCR::performance(roc_pred_val, measure = "auc")
      
      if(method == 'nn'){
        model <- getModel(data = as.matrix(cur_x_train))
        es_cb <- callback_early_stopping(patience = 15, monitor = 'val_loss', 
                                         mode = "min", 
                                         restore_best_weights = TRUE)
        
        model %>% fit(
          as.matrix(cur_x_train[, cur_features]), 
          to_categorical(as.numeric(cur_y_train) - 1),
          epochs = 100,
          batch_size = 10,
          validation_split = 0.1,
          callbacks = es_cb,
          verbose = 0
        )
        
        validation_score <- evaluate(model, as.matrix(cur_x_val), 
                                            to_categorical(as.numeric(cur_y_val) - 1), 
                                            batch_size = 1, verbose = 0, sample_weight = NULL)[[2]]
      } else if(method == 'glm'){
        # model <-  glmnet::cv.glmnet(as.matrix(cur_x_train[, cur_features]), 
        #                             cur_y_train, family = 'binomial',
        #                             alpha = 1, nfolds = 5, 
        #                             type.measure = 'auc')
        # 
        # predictions <- predict(model, as.matrix(cur_x_val[, cur_features]), 
        #                        type = 'response', 
        #                        s = 'lambda.min')
        # roc_pred_val <- ROCR::prediction(predictions,
        #                                  cur_y_val)
        # 
        # auc_perf <- ROCR::performance(roc_pred_val, measure = "auc")
        # validation_score <- auc_perf@y.values[[1]]
        
        validation_score <- optimizeLR(cur_x_train[,cur_features], cur_x_val[, cur_features], 
                                       cur_y_train, cur_y_val)
      } else if(method == 'svm') {
        model <-  SKM::support_vector_machine(
          as.matrix(cur_x_train[, cur_features]),
          factor(cur_y_train),
          verbose = F,
          seed = 1994)
        predictions <- predict(model, as.matrix(cur_x_val[, cur_features]))
        
        roc_pred_val <- ROCR::prediction(predictions$probabilities[, 2], cur_y_val)
        
        auc_perf <- ROCR::performance(roc_pred_val, measure = "auc")
        validation_score <- auc_perf@y.values[[1]]
      }
      
    
      
      res[[count]] <- validation_score
      # message(paste("Current score:", validation_score))
      top_features <- c(top_features, cur_features)

      
      
      count <- count + 1
      
    }
  }
  
  
  return(list(scores = res, features = top_features))
}

getModel <- function(data, l1_l2_reg = 0){
  num_features <- ncol(data)
  # units <- 2 ^ floor(log(num_features, base = 2))
  units <- (num_features %/% 2) * 2
  model <- keras_model_sequential() %>%
    layer_dense(units = num_features, activation = 'relu', 
                input_shape = c(num_features),
                kernel_regularizer = regularizer_l1_l2(l1_l2_reg)) %>%
    layer_dropout(rate = 0) %>%
    # layer_dense(units = 100, activation = 'relu',
    #             kernel_regularizer = regularizer_l1_l2(l1_l2_reg)) %>%
    #  layer_dropout(rate = 0.0) %>%
    # layer_dense(units = 50, activation = 'relu',
    #             kernel_regularizer = regularizer_l1_l2(l1_l2_reg)) %>%
    #  layer_dropout(rate = 0.0) %>%
    # layer_dense(units = 25, activation = 'relu',
    #             kernel_regularizer = regularizer_l1_l2(l1_l2_reg)) %>%
    # layer_dropout(rate = 0.0) %>%
    layer_dense(units = 2, activation = 'sigmoid', 
                input_shape = c(num_features))
  
  initial_learning_rate <- 0.01
  lr_schedule <- learning_rate_schedule_exponential_decay(
    initial_learning_rate = initial_learning_rate,
    decay_steps = 100000,
    decay_rate = 0.96
  )
  opt <- optimizer_nadam(learning_rate = lr_schedule, 
                         beta_1 = 0.975, 
                         beta_2 = 0.999)
  
  # Compile the model
  model %>% compile(
    loss = 'binary_crossentropy',
    optimizer = opt,
    metrics = metric_auc()
  )
  
  return(model)
}

get_nparam <- function(mod, numvar) {
  
  coef(mod, s = with(mod, min(lambda[df == numvar])))
  
}


KFoldCustomNN <- function(train_x, train_y, k = 5, nRepeats = 1){
  set.seed(1994)
  seeds <- sample(1:1e6, nRepeats, replace=FALSE)
  count <- 1
  truth_list <- c()
  predictions_list <- c()
  res <- rep(k*nRepeats)
  for(i in 1:nRepeats){
    set.seed(seeds[[i]])
    res_cur <- c()
    cv <- SKM::cv_kfold_strata(train_y, k = k)
    
    message(paste('Now running iteration:', i))
    for(j in 1:k){
      idx <- cv[[j]]$training
      idx_val <- cv[[j]]$testing
      cur_x_train <- train_x[idx, ]
      cur_y_train <- train_y[idx]
      
      cur_x_val <- train_x[idx_val, ]
      cur_y_val <- train_y[idx_val]
      
      # instantiate the nn model
      cur_model <- getModel(data = cur_x_train)
      
      # train the model
      es_cb <- callback_early_stopping(patience = 30, monitor = 'val_loss', 
                                       restore_best_weights = TRUE)
      
      cur_model %>% fit(
        as.matrix(cur_x_train), 
        to_categorical(as.numeric(cur_y_train) - 1),
        epochs = 100,
        batch_size = 5,
        validation_split = 0.1,
        callbacks = es_cb,
        verbose = 0
      )
    
      # make predictions
      auc_j <- evaluate(cur_model, as.matrix(cur_x_val),
                                          to_categorical(as.numeric(cur_y_val) - 1),
                                          batch_size = 5, verbose = 0, sample_weight = NULL)[[2]]
      
      # predictions <- cur_model %>% predict(as.matrix(cur_x_val))
      # roc_pred_val <- ROCR::prediction(predictions, to_categorical(as.numeric(cur_y_val) - 1))

      # auc_perf <- ROCR::performance(roc_pred_val, measure = "auc")
      message(auc_j)
      # auc_j <- auc_perf@y.values[[1]]
      res[[count]] <- auc_j
      
      # predictions_list[[count]] <- predictions[, 1]
      # truth_list[[count]] <- as.numeric(cur_y_val)
      
      count <- count + 1
    }
    # res[[i]] <- round(mean(unlist(res_cur)), 2)
  }
  
  
  message(paste('The mean CV AUC is:', round(mean(unlist(res)), 2)))
  cv_dat <- list(predictions = predictions_list,
                 labels = truth_list)
  
  return(list(cv_res = unlist(res), cv_dat = cv_dat))
}

KFoldCustomMLGL <- function(train_x, train_y, k = 5, nRepeats = 1, nFolds = 5){
  set.seed(1994)
  seeds <- sample(1:1e6, nRepeats, replace=FALSE)
  count <- 1
  truth_list <- c()
  predictions_list <- c()
  res <- rep(k*nRepeats)
  for(i in 1:nRepeats){
    set.seed(seeds[[i]])
    res_cur <- c()
    cv <- SKM::cv_kfold_strata(train_y, k = k)
    
    message(paste('Now running iteration:', i))
    for(j in 1:k){
      idx <- cv[[j]]$training
      idx_val <- cv[[j]]$testing
      cur_x_train <- train_x[idx, ]
      cur_y_train <- train_y[idx]
      
      cur_x_val <- train_x[idx_val, ]
      cur_y_val <- train_y[idx_val]
      
      dat <- cur_x_train
      cur_y_train_num <- ifelse(cur_y_train == 0, -1, 1)
      cur_y_val_num <- ifelse(cur_y_val == 0, -1, 1)
      
      hc <- bootstrapHclust(dat,frac=1,B=100, nCore = 7)
      groupWeight <- computeGroupSizeWeight(hc)
      
      res_mgml <- MLGL(as.matrix(dat), cur_y_train_num, loss = 'logit', hc = hc, 
                  weightSizeGroup = groupWeight, verbose = TRUE)
      
      out <- HMT(res_mgml,as.matrix(dat), cur_y_train_num, control="FWER",alpha=0.1)
      summary(out)
      
      cur_x_train <- dat[, out$var] %>%
        as.matrix()
      cur_x_val <- cur_x_val[, out$var] %>%
        as.matrix()
      if(length(out$var) < 2){
        next
      }
      
      set.seed(999)
      message("Now running Adaptive Lasso")
      cv.ridge <- cv.glmnet(x = cur_x_train, y = cur_y_train, family='binomial', nfolds = nFolds,
                            alpha=0, parallel=TRUE, standardize=TRUE)
      best_ridge_coef <- as.numeric(coef(cv.ridge, s = cv.ridge$lambda.min))[-1]
      
      ## Adaptive Lasso
      set.seed(999)
      cv.lasso <- cv.glmnet(x = cur_x_train, y = cur_y_train, family='binomial', alpha=1, 
                            parallel=TRUE, standardize=TRUE,  nfolds = nFolds,
                            type.measure='auc', penalty.factor=1/(abs(best_ridge_coef)))
      
      message("Now running Priority Lasso")
      
      
      pred_test <- predict(cv.lasso, newx = as.matrix(cur_x_val), 
                           type = "response", s = 'lambda.min')
      validation_score <- my.auc(pred_test, cur_y_val)
      
      res[[count]] <- validation_score
      
      message(paste('The current AUC is:', round(validation_score, 2)))
      
      count <- count + 1
    }
    # res[[i]] <- round(mean(unlist(res_cur)), 2)
  }
  
  
  message(paste('The mean CV AUC is:', round(mean(unlist(res)), 2)))
  # cv_dat <- list(predictions = predictions_list,
  #                labels = truth_list)
  
  # return(list(cv_res = unlist(res), cv_dat = cv_dat))
}


optimizeLR <- function(dat_train, dat_test, y_train, y_test){
  
  dat_train <- cbind(dat_train, y_train)
  dat_test <- cbind(dat_test, y_test)
  
  set.seed(1994)
  val_set <- vfold_cv(dat_train, v = 5, repeats = 5)
  
  lr_mod <- 
    logistic_reg(penalty = tune(), mixture = 1) %>% 
    set_engine("glmnet")
  
  lr_recipe <- 
    recipe(y_train ~ ., data = dat_train) %>% 
    step_zv(all_predictors()) %>% 
    step_normalize(all_predictors())
  
  lr_workflow <- 
    workflow() %>% 
    add_model(lr_mod) %>% 
    add_recipe(lr_recipe)
  
  lr_reg_grid <- tibble(penalty = 10^seq(-4, -1, length.out = 30))
  
  lr_res <- 
    lr_workflow %>% 
    tune_grid(val_set,
              grid = lr_reg_grid,
              control = control_grid(save_pred = TRUE),
              metrics = metric_set(roc_auc))
  lr_best <-
    lr_res %>% 
    select_best()
  
  top_models <-
    lr_res %>% 
    show_best("roc_auc", n = 5) %>% 
    arrange(.metric) 
  # print(top_models)
  
  penalty_val <- lr_best %>%
    pull(penalty)
  
  last_lr_mod <- 
    logistic_reg(penalty = penalty_val, mixture = 1) %>% 
    set_engine("glmnet") %>% 
    set_mode("classification") %>%
    # Fit the model
    parsnip::fit(y_train~., data = dat_train)
  
  # Prediction Probabilities
  pred_proba_lr <- predict(last_lr_mod,
                           new_data = dat_test,
                           type = "prob")
  
  
  test_lr_results <- 
    bind_cols(y_test, pred_proba_lr)
  colnames(test_lr_results)[[1]] <- 'condition'
  
  auc <- test_lr_results %>%
    roc_auc(truth = condition, .pred_0)
  
  return(auc$.estimate)
}

KFoldCustomAdaptive <- function(train_x, train_y, k = 5, nRepeats = 1, nfolds = 5){
  set.seed(1994)
  seeds <- sample(1:1e6, nRepeats, replace=FALSE)
  count <- 1
  truth_list <- c()
  predictions_list <- c()
  res <- rep(k*nRepeats)
  for(i in 1:nRepeats){
    set.seed(seeds[[i]])
    res_cur <- c()
    cv <- SKM::cv_kfold_strata(train_y, k = k)
    
    message(paste('Now running iteration:', i))
    for(j in 1:k){
      idx <- cv[[j]]$training
      idx_val <- cv[[j]]$testing
      cur_x_train <- train_x[idx, ]
      cur_y_train <- train_y[idx]
      
      cur_x_val <- train_x[idx_val, ]
      cur_y_val <- train_y[idx_val]
      
      cur_x_train <- scale(cur_x_train)
      cur_x_val <- scale(cur_x_val, center=attr(cur_x_train, "scaled:center"), 
                        scale=attr(cur_x_train, "scaled:scale"))
      
      cur_y_train <- if_else(cur_y_train == 1, 1, 0) 
      
      set.seed(999)
      cv.ridge <- cv.glmnet(cur_x_train, cur_y_train, family='binomial', alpha=0, 
                            parallel=TRUE, nfolds = nfolds,
                            standardize=F, relax = F)
      
      best_ridge_coef <- as.numeric(coef(cv.ridge, s = cv.ridge$lambda.min))[-1]
      ## Adaptive Lasso
      set.seed(999)
      fit.cv <- cv.glmnet(cur_x_train, cur_y_train, family='binomial', 
                           alpha=1, parallel=TRUE, nfolds = nfolds,
                           standardize=F, relax = F, type.measure = 'auc',
                           penalty.factor=1/(abs(best_ridge_coef)))
      
      preds <- predict(fit.cv, newx=as.matrix(cur_x_val), 
                       type = 'response', 
                       s = fit.cv$lambda.min)
      
      pred_fs <- prediction(preds, cur_y_val)
      perf_fs <- ROCR::performance(pred_fs, "tpr", "fpr")
      
      
      auc_perf <- ROCR::performance(pred_fs, measure = "auc")
      validation_score <- auc_perf@y.values[[1]]
      
      res[[count]] <- validation_score
      
      message(paste('The current AUC is:', round(validation_score, 2)))
      
      count <- count + 1
    }
    # res[[i]] <- round(mean(unlist(res_cur)), 2)
  }
  
  
  message(paste('The mean CV AUC is:', round(mean(unlist(res)), 2)))
  # cv_dat <- list(predictions = predictions_list,
  #                labels = truth_list)
  
  # return(list(cv_res = unlist(res), cv_dat = cv_dat))
}


KFoldCustomPLasso <- function(train_x, train_y, k = 5, nRepeats = 1, 
                              nFolds = 10, n = 20){
  set.seed(1994)
  seeds <- sample(1:1e6, nRepeats, replace=FALSE)
  count <- 1
  truth_list <- c()
  predictions_list <- c()
  res <- rep(k*nRepeats)
  for(i in 1:nRepeats){
    set.seed(seeds[[i]])
    res_cur <- c()
    cv <- SKM::cv_kfold_strata(train_y, k = k)
    
    message(paste('Now running iteration:', i))
    for(j in 1:k){
      idx <- cv[[j]]$training
      idx_val <- cv[[j]]$testing
      cur_x_train <- train_x[idx, ]
      cur_y_train <- train_y[idx]
      
      cur_x_val <- train_x[idx_val, ]
      cur_y_val <- train_y[idx_val]
      
      # dat <- cur_x_train %>%
      #   scale()
      cur_y_train_num <- ifelse(cur_y_train == 0, -1, 1)
      cur_y_val_num <- ifelse(cur_y_val == 0, -1, 1)
      
      idx <- nestedcv::collinear(cur_x_train, rsq_cutoff = 0.75)
      if(length(idx) > 0){
        cur_x_train <- cur_x_train[, -c(idx)]
        cur_x_val <- cur_x_val[, -c(idx)]
      }
      
      if(ncol(cur_x_train) < 2){
        next
      }
      
      h_res <- runHOPACH(t(cur_x_train), kmax = )
      
      res_hopach <- computeHOPACHClusters(h_res = h_res, dat = cur_x_train, group = TRUE)
      cur_x_train <- res_hopach$meta_features %>%
        as.matrix()
      groups <- findSequenceIndices(res_hopach$groups)
      
      cur_x_val <- computeHOPACHClusters(h_res = h_res, dat = cur_x_val, group = FALSE) %>%
        as.matrix()
      
      require(doMC)
      registerDoMC(cores = 4)
      
      pl_fit1 <- prioritylasso(X = cur_x_train, Y = cur_y_train, family = "binomial", type.measure = "auc", 
                               blocks = groups, 
                               nfolds = nFolds, parallel =TRUE,
                               standardize = FALSE, lambda.type = "lambda.min", block1.penalization = T)
      
      # pl_fit1$min.cvm
      
      
      coeff1 <- pl_fit1$coefficients
      coeff1 <- coeff1[coeff1 != 0] %>%
        na.omit()
      # print(round(coeff1, 4))
      
      pl1_score <- cur_x_val[ , names(coeff1), drop=F] %*% coeff1
      pl1_roc <- roc(factor(cur_y_val), pl1_score[,1])
      
      validation_score <- auc(pl1_roc)
      
      res[[count]] <- validation_score
      
      message(paste('The current AUC is:', round(validation_score, 2)))
      
      count <- count + 1
    }
    # res[[i]] <- round(mean(unlist(res_cur)), 2)
  }
  
  
  message(paste('The mean CV AUC is:', round(mean(unlist(res)), 2)))
  # cv_dat <- list(predictions = predictions_list,
  #                labels = truth_list)
  
  # return(list(cv_res = unlist(res), cv_dat = cv_dat))
}


KFoldCustomEnsemble <- function(dat_prop, dat_logit, dat_mean, train_y,
                                k = 5, nRepeats = 1, alpha = 0.5, seed = 1994,
                              nFolds = 10, linkage = 'ward', numMarkers = 27){
  set.seed(1994)
  seeds <- sample(1:1e6, nRepeats, replace=FALSE)
  count <- 1
  truth_list <- c()
  predictions_list <- c()
  res <- rep(k*nRepeats)
  for(i in 1:nRepeats){
    set.seed(seeds[[i]])
    
    # shuffle <- sample(1:nrow(train_x))
    # train_x <- train_x[shuffle, ]
    # train_y <- train_y[shuffle]
    # 
    res_cur <- c()
    cv <- SKM::cv_kfold_strata(train_y, k = k)
    
    message(paste('Now running iteration:', i))
    for(j in 1:k){
      idx <- cv[[j]]$training
      idx_val <- cv[[j]]$testing
      
      cur_x_train_prop <- dat_prop[idx, ]
      cur_x_train_logit <- dat_logit[idx, ]
      cur_x_train_mean <- dat_mean[idx, ]
      cur_y_train <- train_y[idx]
      
      cur_x_val_prop <- dat_prop[idx_val, ]
      cur_x_val_logit <- dat_logit[idx_val, ]
      cur_x_val_mean <- dat_mean[idx_val, ]
      cur_y_val <- train_y[idx_val]
      
      print(table(cur_y_val))
      
      d_prop <- dist(t(cur_x_train_prop))
      
      hc_prop <- hclust(d_prop, method = linkage)
      hc_prop <- treekoR:::findChildren(ggtree(as.phylo(hc_prop), 
                                               ladderize = F, layout = 'dendrogram'))
      # subs_prop <- hc_prop$data$clusters
      # prop_dat <- hc_prop$data
      # prop_dat$group <- 1:length(subs_prop)
      
      d_logit <- dist(t(cur_x_train_logit))
      hc_logit <- hclust(d_logit, method = linkage)
      hc_logit <- treekoR:::findChildren(ggtree(as.phylo(hc_logit),
                                                ladderize = F, layout = 'dendrogram'))
      subs_logit <- hc_logit$data$clusters
      logit_dat <- hc_logit$data
      logit_dat$group <- 1:length(subs_logit)
      
      
      var_groups<- list()
      n <- length(subs_logit)
      # n2 <- 2*n
      nclust <- ncol(dat_logit)
      for(i in 1:length(subs_logit)){
        var_groups[[i]] <- which(colnames(dat_logit) %in% subs_logit[[i]])
        # var_groups[[i+n]] <- which(colnames(data_logit) %in% subs_logit[[i]]) + nclust
      }
      
      last_max <- ncol(data_logit)
      nclust <- last_max
      off_set <- length(var_groups)
      for (i in 1:nclust) {
        index <- off_set + i  # Calculate the new index in the list
        begin <- last_max + 1        # Start the next sequence right after the last max
        end <- begin + numMarkers - 1            # Each sublist has exactly 27 items
        
        # Add the new sequence to the res list at the new index
        var_groups[[index]] <- seq(begin, end)
        
        # Update last_max for the next iteration
        last_max <- end
      }
      
      cur_x_train <- cbind(cur_x_train_logit, cur_x_train_mean)
      scaleVals <- preProcess(cur_x_train, method = c('range'))
      cur_x_train <- predict(scaleVals, cur_x_train) %>%
        as.matrix()
      
      cur_x_val <- cbind(cur_x_val_logit, cur_x_val_mean) 
      
      cur_x_val <- predict(scaleVals, cur_x_val) %>%
        as.matrix()
      
      cur_y_train <- if_else(cur_y_train == 1, 1, 0) 
      
      set.seed(999)
      cvfit <- cv.grpregOverlap(cur_x_train, cur_y_train, var_groups, 
                                penalty = 'cMCP', family = 'binomial',  
                                alpha = alpha, seed = seed)
      
      fit <- grpregOverlap(cur_x_train, cur_y_train, var_groups, 
                           family='binomial', alpha = alpha,
                           returnX.latent = T, returnOverlap = FALSE, 
                           lambda = cvfit$lambda.min, penalty = 'cMCP')
      
      pred_test <- predict(fit, X = cur_x_val,
                           type = "response")
      
      
      validation_score <- my.auc(pred_test, cur_y_val)
      
      res[[count]] <- validation_score
      
      message(paste('The current AUC is:', round(validation_score, 2)))
      
      predictions_list[[count]] <- pred_test
      truth_list[[count]] <- as.numeric(cur_y_val)
      
      count <- count + 1
    }
    # res[[i]] <- round(mean(unlist(res_cur)), 2)
  }
  
  mean_cv <- unlist(res) %>% 
    na.omit() %>% 
    mean()
  message(paste('The mean CV AUC is:', round(mean_cv, 2)))
  cv_dat <- list(predictions = predictions_list,
                 labels = truth_list)
  
  return(list(cv_res = unlist(res), cv_dat = cv_dat))
}



KFoldCustomEnsembleAlpha <- function(dat_prop, dat_logit, dat_mean, train_y,
                                k = 5, nRepeats = 1, seed = 1994, root = NULL,
                                nFolds = 10, linkage = 'average', numMarkers = 27,
                                alpha_search = seq(0.01, 0.1, by = 0.01), penalty = 'grMCP') {
  set.seed(seed)
  seeds <- sample(1:1e6, nRepeats, replace=FALSE)
  best_alpha <- NA
  best_mean_cv <- 0
  
  for (alpha in alpha_search) {
    message(paste("Testing alpha:", alpha))
    res <- vector("list", k * nRepeats)
    count <- 1
    predictions_list <- list()
    truth_list <- list()
    
    for (i in 1:nRepeats) {
      set.seed(seeds[i])
      cv <- SKM::cv_kfold_strata(train_y, k = k)
      
      for (j in 1:k) {
        idx <- cv[[j]]$training
        idx_val <- cv[[j]]$testing
        
        cur_x_train_logit <- dat_logit[idx, ]
        cur_x_train_mean <- dat_mean[idx, ]
        cur_y_train <- train_y[idx]
        
        cur_x_val_logit <- dat_logit[idx_val, ]
        cur_x_val_mean <- dat_mean[idx_val, ]
        cur_y_val <- train_y[idx_val]
        
        if(!is.null(root)){
          hc_logit <- treekoR:::findChildren(ggtree(as.phylo(root), 
                                                    ladderize = F, layout = 'dendrogram'))
        } else{
          d_logit <- dist(t(cur_x_train_logit))
          hc_logit <- hclust(d_logit, method = linkage)
          hc_logit <- treekoR:::findChildren(ggtree(as.phylo(hc_logit),
                                                    ladderize = F, layout = 'dendrogram'))
        }
        
        subs_logit <- hc_logit$data$clusters
        logit_dat <- hc_logit$data
        logit_dat$group <- 1:length(subs_logit)
        
        
        var_groups<- list()
        n <- length(subs_logit)
        nclust <- ncol(dat_logit)
        for(i in 1:length(subs_logit)){
          var_groups[[i]] <- which(make_clean_names(colnames(dat_logit)) %in% make_clean_names(subs_logit[[i]]))
        }
        
        last_max <- ncol(dat_logit)
        nclust <- last_max
        off_set <- length(var_groups)
        for (i in 1:nclust) {
          index <- off_set + i  # Calculate the new index in the list
          begin <- last_max + 1        # Start the next sequence right after the last max
          end <- begin

          # Add the new sequence to the res list at the new index
          var_groups[[index]] <- seq(begin, end)

          # Update last_max for the next iteration
          last_max <- end
        }
        # var_groups[off_set + 1] <- list(seq(last_max + 1, ncol(dat_mean) + last_max))
        
        # Combine and scale training data
        cur_x_train <- cbind(cur_x_train_logit, cur_x_train_mean)
        scaleVals <- preProcess(cur_x_train, method = c('range'))
        cur_x_train <- predict(scaleVals, cur_x_train) %>% as.matrix()
        
        # Prepare and scale validation data
        cur_x_val <- cbind(cur_x_val_logit, cur_x_val_mean)
        cur_x_val <- predict(scaleVals, cur_x_val) %>% as.matrix()
        
        cur_y_train <- if_else(cur_y_train == 1, 1, 0) 
        
        # grp_weights <- compute_weights(cur_x_train, cur_y_train, 
        #                                groups = var_groups)
        
        # Fit model using grpregOverlap
        cvfit <- cv.grpregOverlap(cur_x_train, cur_y_train, var_groups, 
                                  nfolds = 10,
                                  penalty = penalty, family = 'binomial',  
                                  alpha = alpha, seed = seed)
        
        fit <- grpregOverlap(cur_x_train, cur_y_train, var_groups, 
                             family='binomial',
                             returnX.latent = T, returnOverlap = FALSE, 
                             lambda = cvfit$lambda.min, penalty = penalty)
        
        pred_test <- predict(fit, X = cur_x_val, type = "response")
        
        validation_score <- my.auc(pred_test, cur_y_val)
        message(paste('The current AUC is:', round(validation_score, 2)))
        
        res[[count]] <- validation_score
        predictions_list[[count]] <- pred_test
        truth_list[[count]] <- as.numeric(cur_y_val)
        
        count <- count + 1
      }
    }
    
    mean_cv <- unlist(res) %>% na.omit() %>% mean()
    message(paste('The mean CV AUC for alpha', alpha, 'is:', round(mean_cv, 2)))
    
    # Track the best alpha
    if (mean_cv > best_mean_cv) {
      best_mean_cv <- mean_cv
      best_alpha <- alpha
    }
  }
  
  message(paste("Best alpha:", best_alpha, "with mean CV AUC:", round(best_mean_cv, 2)))
  
  cv_dat <- list(predictions = predictions_list, labels = truth_list)
  
  return(list(best_alpha = best_alpha, cv_res = unlist(res), cv_dat = cv_dat, best_mean_cv = best_mean_cv))
}

selectApha <- function(X_train, y_train, groups, penalty = "cMCP", weights,
                      alpha_search = seq(0.1, 1, by = 0.1), seed = 999){
  best_aic <- Inf
  aics <- numeric(length = length(alpha_search))
  for(i in 1:length(alpha_search)){
    alpha <- alpha_search[i]
    message(paste("Testing alpha:", alpha))
    
    cvfit <- cv.grpregOverlap(X_train, y_train, groups, 
                              nfolds = 5, seed = seed,
                              penalty = penalty, family = 'binomial',  
                              alpha = alpha)
    
    fit <- grpregOverlap(X_train, y_train, groups,
                         family='binomial', alpha = alpha,
                         returnX.latent = T, returnOverlap = FALSE, 
                         lambda = cvfit$lambda.min, penalty = penalty)
    cur_aic <- BIC(fit) %>% round(2)
    aics[i] <- cur_aic
  if(cur_aic < best_aic){
    best_alpha <- alpha
    message(paste("current best AIC is:", cur_aic))
    best_aic <- cur_aic
    best_fit <- fit
  }
      
    
  }
  
  return(list(best_fit = best_fit, best_aic = best_aic, best_alpha = best_alpha,
              aics = aics))
  
}



GlassoCV <- function(dat_logit, dat_mean, train_y, root = NULL,
                     k = 5, nRepeats = 1, alpha = 0.02, seed = 1994,
                     nFolds = 10, linkage = 'ward', penalty = 'grMCP'){
    set.seed(1994)
    seeds <- sample(1:1e6, nRepeats, replace=FALSE)
    count <- 1
    truth_list <- c()
    predictions_list <- c()
    res <- rep(k*nRepeats)
    for(i in 1:nRepeats){
      set.seed(seeds[[i]])
      
      # shuffle <- sample(1:nrow(train_x))
      # train_x <- train_x[shuffle, ]
      # train_y <- train_y[shuffle]
      # 
      res_cur <- c()
      cv <- SKM::cv_kfold_strata(train_y, k = k)
      
      message(paste('Now running iteration:', i))
      for(j in 1:k){
        idx <- cv[[j]]$training
        idx_val <- cv[[j]]$testing
        
        cur_x_train_logit <- dat_logit[idx, ]
        cur_x_train_mean <- dat_mean[idx, ]
        cur_y_train <- train_y[idx]
        
        cur_x_val_logit <- dat_logit[idx_val, ]
        cur_x_val_mean <- dat_mean[idx_val, ]
        cur_y_val <- train_y[idx_val]
        
        if(!is.null(root)){
          hc_logit <- treekoR:::findChildren(ggtree(as.phylo(root), 
                                                    ladderize = F, layout = 'dendrogram'))
        } else{
          d_logit <- coop::pcor(cur_x_train_logit) %>%
            cor2dist() %>%
            as.dist()
          hc_logit <- hclust(d_logit, method = linkage)
          hc_logit <- treekoR:::findChildren(ggtree(as.phylo(hc_logit),
                                                    ladderize = F, layout = 'dendrogram'))
        }
        
        subs_logit <- hc_logit$data$clusters
        logit_dat <- hc_logit$data
        logit_dat$group <- 1:length(subs_logit)
        
        
        var_groups<- list()
        n <- length(subs_logit)
        nclust <- ncol(dat_logit)
        for(i in 1:length(subs_logit)){
          var_groups[[i]] <- which(make_clean_names(colnames(dat_logit)) %in% make_clean_names(subs_logit[[i]]))
        }
        
        last_max <- ncol(dat_logit)
        nclust <- last_max
        off_set <- length(var_groups)
        for (i in 1:nclust) {
          index <- off_set + i  # Calculate the new index in the list
          begin <- last_max + 1        # Start the next sequence right after the last max
          end <- begin
          
          # Add the new sequence to the res list at the new index
          var_groups[[index]] <- seq(begin, end)
          
          # Update last_max for the next iteration
          last_max <- end
        }
        # var_groups[off_set + 1] <- list(seq(last_max + 1, ncol(dat_mean) + last_max))
        
        # Combine and scale training data
        cur_x_train <- cbind(cur_x_train_logit, cur_x_train_mean)
        scaleVals <- preProcess(cur_x_train, method = c('range'))
        cur_x_train <- predict(scaleVals, cur_x_train) %>% as.matrix()
        
        # Prepare and scale validation data
        cur_x_val <- cbind(cur_x_val_logit, cur_x_val_mean)
        cur_x_val <- predict(scaleVals, cur_x_val) %>% as.matrix()
        
        cur_y_train <- if_else(cur_y_train == 1, 1, 0)
        
        alpha_search = seq(0.01, 0.1, 0.01)
        res_aic_cur <- selectApha(cur_x_train, cur_y_train, alpha_search = alpha_search,
                              groups = var_groups, penalty = penalty)
        
        alpha <- alpha_search[computeElbow(res_aic_cur$aics)]
        
        set.seed(999)
        cvfit <- cv.grpregOverlap(cur_x_train, cur_y_train, var_groups, 
                                  penalty = penalty, family = 'binomial',  
                                  alpha = alpha, seed = seed)
        
        fit <- grpregOverlap(cur_x_train, cur_y_train, var_groups, 
                             family='binomial', alpha = alpha,
                             returnX.latent = T, returnOverlap = FALSE, 
                             lambda = cvfit$lambda.min, penalty = penalty)
        
        pred_test <- predict(fit, X = cur_x_val,
                             type = "response")
        
        
        validation_score <- my.auc(pred_test, cur_y_val)
        
        res[[count]] <- validation_score
        
        message(paste('The current AUC is:', round(validation_score, 2)))
        
        predictions_list[[count]] <- pred_test
        truth_list[[count]] <- as.numeric(cur_y_val)
        
        count <- count + 1
      }
      # res[[i]] <- round(mean(unlist(res_cur)), 2)
    }
    
    mean_cv <- unlist(res) %>% 
      na.omit() %>% 
      mean()
    message(paste('The mean CV AUC is:', round(mean_cv, 2)))
    cv_dat <- list(predictions = predictions_list,
                   labels = truth_list)
    
    return(list(cv_res = unlist(res), cv_dat = cv_dat))
  }