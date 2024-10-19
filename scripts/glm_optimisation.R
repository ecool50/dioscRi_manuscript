


OLGLasso <- list(type = "Classification",
              library = "grpregOverlap",
              loop = NULL) 

## The parameters Element
prm <- data.frame(parameter = "alpha",
                  class = 'numeric',
                  label = "Alpha")

OLGLasso$parameters <- prm

# The grid element
OLGLassoGrid <- function(x, y, len = NULL, search = "grid") {
 
  if(search == "grid") {
    out <- data.frame(alpha = seq(.05, 1, length = len))
  } else {
    ## For random search, define ranges for the parameters then
    ## generate random values for them
    out <- data.frame(alpha = seq(.05, 1, length = len))
  }
  out
}

OLGLasso$grid <- OLGLassoGrid

# The fit Element
OLGLassoFit <- function(x, y, wts, param, lev, last, weights, classProbs, group, ...) { 
  grpregOverlap(x, y, groups, nfolds = 10,
                family='binomial', alpha = param$alpha,
                returnX.latent = F, returnOverlap = FALSE,
                penalty = 'grLasso')
}

OLGLasso$fit <- OLGLassoFit

# The predict Element
OLGLassoPred <- function(modelFit, newdata, preProc = NULL, submodels = NULL){
  grpregOverlap:::predict.grpregOverlap(modelFit, X = newdata,
                                        type = "response", which = 'lambda.min')
}
  
OLGLasso$predict <- OLGLassoPred

# The prob Element
OLGLasso$prob <- OLGLassoPred

# The sort Element
OLGLassoSort <- function(x) x[order(x$alpha),]
OLGLasso$sort <- OLGLassoSort

#  The levels Element


fitControl <- trainControl(method = "repeatedcv",
                           ## 10-fold CV...
                           number = 10,
                           ## repeated ten times
                           repeats = 3)

temp <- cbind(train_x, condition)
set.seed(825)
Laplacian <- caret::train(condition ~ ., data = temp, 
                   method = OLGLasso,
                   tuneLength = 8,
                   trControl = fitControl)
Laplacian