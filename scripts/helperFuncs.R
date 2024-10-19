#' This function calculates the proportion of cells belonging to 
#' each cell type in each sample, applicable to scRNA-seq and 
#' spatial proteomics. It also takes an optional logit argument 
#' that specifies whether the proportions should be logit-transformed.
#' @noRd
helper_proportion_raw <- function(data, logit = TRUE) {
  df <- data@meta.data
  df <- table(df$sample, data$celltype)
  df <- df / rowSums(df)
  
  # if logit transformation is needed,
  # need to do do the following to avoid infinite or NA
  if (logit) {
    try(df[df == 0] <- 0.001, silent = TRUE)
    try(df[df == 1] <- 0.999, silent = TRUE)
    df <- gtools::logit(df)
  }
  
  
  df <- as.data.frame(df)
  
  colnames(df) <- c("sample", "celltype", "proportion")
  
  
  df <- df |>
    tidyr::pivot_wider(names_from = "celltype", values_from = "proportion")
  df <- as.data.frame(df)
  rownames(df) <- df$sample
  df <- df[, -1]
  
  df <- df[unique(data$sample), ]
  
  return(df)
}

getProp <- function(cells, feature = "clusters", 
                    imageID = "sample_id", logit = TRUE) {
  if (is.data.frame(cells)) {
    df <- cells[, c(imageID, feature)]
  }
  
  if (is(cells, "SingleCellExperiment") | is(cells, "SpatialExperiment")) {
    df <- as.data.frame(SummarizedExperiment::colData(cells))[, c(imageID, feature)]
  }
  
  if (is(cells, "SegmentedCells")) {
    cellSummary <- cellSummary(cells, bind = TRUE)
    df <- as.data.frame(cellSummary[, c(imageID, feature)])
  }
  
  
  tab <- table(df[, imageID], df[, feature])
  tab <- sweep(tab, 1, rowSums(tab), "/")
  
  
  if (logit) {
    try(tab[tab == 0] <- 0.001, silent = TRUE)
    try(tab[tab == 1] <- 0.999, silent = TRUE)
    tab <- boot::logit(tab)
  }
  as.data.frame.matrix(tab)
}

sparsity <- function(input,
                     count.na.zero = FALSE) {
  sparsity <- dplyr::rename(
    tibble::enframe(
      apply(input, 2, function(i) {
        Matrix::nnzero(i, na.counted = count.na.zero)/nrow(input)
      })
    ),
    "VARIABLE"    = "name",
    "PERC_SPARSE" = "value"
  )
  return(sparsity)
}

ComputeReferenceSample <- function(data, markers, sample_col = 'sample_id', N = 3){
  cov_mats <- list()
  num_samples <- length(unique(data[[sample_col]]))
  samples <- unique(data[[sample_col]])
  norms <- matrix(0, num_samples, num_samples)
  
  message("Now precomputing covariance matrices")
  for (i in 1:num_samples) {
    dat_i <- data[data[[sample_col]] == samples[[i]], markers]
    cov_mats[[i]] <- Rfast::cova(as.matrix(dat_i))
  }
  
  message("Now computing pairwise Frobenius norms")
  for(i in 1:num_samples){
    cov_i <- cov_mats[[i]]
    
    for(j in 1:num_samples){
      cov_j <- cov_mats[[j]]
      
      covDiff <- cov_i - cov_j
      norms[i, j] <- Rfast::Norm(covDiff, type = 'F')
      norms[j, i] <- norms[i, j]
    }
  }
  
  avg_norms <- colMeans(na.omit(norms))
  sorted_indices <- order(avg_norms, decreasing = FALSE)
  top_N_samples <- samples[sorted_indices[1:N]]
  bottom_N_samples <- samples[sorted_indices[(length(samples)-N+1):length(samples)]]
  
  refSampleInd <- samples[which.min(avg_norms)]
  medianSampleInd <- samples[which(avg_norms == median(avg_norms))]
  valSampleInd <- samples[which.max(avg_norms)]
  
  return(list('Norms' = norms, 
              'avgNorms' = avg_norms,
              'topNSamples' = top_N_samples,
              'bottomNSamples' = bottom_N_samples,
              'medianSampleInd' = medianSampleInd,
              'refSampleInd' = refSampleInd,
              'valSampleInd' = valSampleInd))
}


processBatch <- function(batch_dat, markers=NULL, cutoff=0.25){
  
  
  train_spar <- sparsity(t(batch_dat))
  spar.val <- quantile(train_spar$PERC_SPARSE, cutoff)
  
  # filter on the median
  filtered_cells <- train_spar[which(train_spar$PERC_SPARSE > spar.val), ]
  
  high_spar <- batch_dat[filtered_cells$VARIABLE,]
  low_spar <- batch_dat[!(row.names(batch_dat) %in% filtered_cells$VARIABLE),] %>%
    as.matrix()
  p <- coop::sparsity(low_spar)
  
  # use h2o to denoise the data
  zero_noise <- ruta::noise_zeros(p = p)
  low_spar_corrupted <- ruta::apply_filter(zero_noise, low_spar)
  l2_penalty_ae = 1e-2
  # 
  encoding_dim <- 25
  input_dim <- ncol(batch_dat)
  
  # de_train <- h2o.deeplearning(
  #   x = seq_along(low_spar_corrupted),
  #   training_frame = low_spar_corrupted,
  #   autoencoder = TRUE,
  #   hidden = c(input_dim,encoding_dim,input_dim),
  #   activation = 'Tanh',
  #   sparse = TRUE,
  #   seed = 1994,
  #   validation_frame = as.h2o(low_spar)
  # )
  # 
  # denoised_batch <- h2o.deepfeatures(de_train, as.h2o(batch_dat), layer = 3) %>%
  #   as.data.frame()
  
  # define the keras model
  model_keras <- keras_model_sequential()
  model_keras %>%
    
    # First hidden layer
    layer_dense(
      units              = encoding_dim,
      activation         = "relu",
      input_shape        = input_dim) %>%
    # Dropout to prevent overfitting
    layer_dropout(rate = 0.1) %>%
    # Second hidden layer
    layer_dense(
      units              = encoding_dim,
      activation         = "relu") %>%
    # Dropout to prevent overfitting
    layer_dropout(rate = 0.1) %>%
    # Output layer
    layer_dense(
      units              = input_dim,
      activation         = "linear") %>%
    # Compile ANN
    compile(
      optimizer = 'rmsprop',
      loss      = 'kl_divergence',
      metrics = 'mae'
    )
  # fit the model
  model_keras %>% fit(
    x                = low_spar_corrupted,
    y                = low_spar_corrupted,
    validation_data = list(low_spar, low_spar),
    batch_size       = 128,
    epochs           = 10,
    validation_split = 0.2
  )
  
  # predict on the full dataset
  denoised_batch <- predict(object = model_keras, x = as.matrix(batch_dat)) %>% 
    as.data.frame()
  return(denoised_batch)
}

robust_scale <- function(data) {
  median_data <- apply(data, 2, median)
  IQR_data <- apply(data, 2, IQR)
  
  scaled_data <- sweep(data, 2, median_data, `-`) / IQR_data
  return(scaled_data)
}

scaling_UV <- function(inputMat) {        
  for (i in 1:ncol(inputMat)) {
    colmean <- colMeans(inputMat, na.rm=TRUE)
    colsd   <- apply( inputMat, 2, function(x) stats::sd(x, na.rm=TRUE) )
    return( t( apply( inputMat, 1, function(x) (x-colmean)/colsd ) ) )  #apply transposes the result when working on rows
  }
}

scaling_mean <- function(inputMat) {        
  for (i in 1:ncol(inputMat)) {
    colmean <- colMeans(inputMat, na.rm=TRUE)
    return( t( apply( inputMat, 1, function(x) (x-colmean)/colmean ) ) )  #apply transposes the result when working on rows
  }
}

# Function to check if column is bimodal
is_bimodal <- function(column) {
  test <- dip.test(column)
  # Here, p-value < 0.05 is used as a criterion, but you can adjust this threshold.
  return(test$p.value < 0.05)
}

min_max_scaling <- function(data_matrix, min_vals = NULL, max_vals = NULL) {
  if(is.null(min_vals) && is.null(max_vals)){
    min_vals <- apply(data_matrix, 2, min) # find minimum value for each column
    max_vals <- apply(data_matrix, 2, max) # find maximum value for each column
    message("True")
  }

  scaled_matrix <- sweep(data_matrix, 2, min_vals, "-") # subtract min from each element
  range_vals <- max_vals - min_vals
  scaled_matrix <- sweep(scaled_matrix, 2, range_vals, "/") # divide by range
  
  return(list(scaled = scaled_matrix, min_vals = min_vals, max_vals = max_vals))
}

# Function to perform within-sample standardization
standardize_data <- function(data) {
  means <- colMeans(data)  # Calculate column means
  sds <- apply(data, 2, sd)  # Calculate column standard deviations
  
  # Standardize each column using the calculated means and standard deviations
  standardized_data <- scale(data, center = means, scale = sds)
  
  return(standardized_data)
}


runFeatureSelection <- function(data, clinicalData, sample_col, response, condition,
                                featureRange = 5, nFolds = 5, method = 'glm', 
                                nRepeats = 10){
  colnames(data) <- janitor::make_clean_names(colnames(data))
  
  topKFeatures <- runKFeatures(data = data, clinicalData = clinicalData, 
                               featureRange = featureRange,k = nFolds,
                               sampleCol = sample_col, response = response)
  data_final <- data %>%
    dplyr::select(unique(topKFeatures))
  
  res <- KFoldCustom(train_x = data_final, train_y = condition, 
                     nRepeats = nRepeats,
                     method = method)
  
  return(round(mean(res$cv_res), 2))
  
}

calculate_ratios <- function(data) {
  num_clusters <- ncol(data)
  cluster_pairs <- combn(num_clusters, 2)
  
  # Name the columns based on cluster combinations
  column_names <- apply(cluster_pairs, 2, function(x) paste(colnames(data)[x], collapse = "_over_"))
  
  calculate_ratio <- function(pair) {
    cluster1 <- data[, pair[1]]
    cluster2 <- data[, pair[2]]
    return(ifelse(cluster2 == 0, 0, cluster1 / (cluster1 + cluster2)))
  }
  
  ratios <- apply(cluster_pairs, 2, calculate_ratio)
  colnames(ratios) <- column_names
  ratios <- log2(ratios + 1)
  return(ratios)
}

computeMetaFeatures <- function(data, nClust=9, type = 'ratio'){
  # Calculate the correlation matrix
  cor_matrix <- cor(data)
  
  # Hierarchical clustering on the correlation matrix
  hc <- hclust(as.dist(cor2dist(cor_matrix)), method = 'ward.D2')  # Use '1 - correlation' to cluster highly correlated features
  
  # Cut the tree to form clusters - choose 'k' based on the dendrogram
  clusters <- cutree(hc, nClust)
  
  data <- as.matrix(data)
  
  # Create meta-features
  meta_features <- data.frame(matrix(nrow = nrow(data), ncol = nClust))
  
  for (i in 1:nClust) {
    cluster_features <- which(clusters == i)
    if(type == 'ratio'){
      meta_features[, i] <- rowMedians(data[, cluster_features, drop = FALSE])
    } else{
      meta_features[, i] <- rowMeans(data[, cluster_features, drop = FALSE])
    }
    
  }
  
  # Rename the meta features
  colnames(meta_features) <- paste("meta_feature", 1:nClust, sep = "_")
  
  rownames(meta_features) <- rownames(data)
  return(list(Features = meta_features, clusters = clusters))
}




# Adjusting the threshold dynamically
correlationThreshold <- function(correlationmatrix, percentile = 80) {
  threshold <- quantile(abs(correlationmatrix[lower.tri(correlationmatrix)]), 
                        probs = percentile / 100)
  return(threshold)
}

performClustering <- function(graph, method = "louvain") {
  if (method == "louvain") {
    return(cluster_louvain(graph))
  } else if (method == "walktrap") {
    return(cluster_walktrap(graph))
  } else if (method == 'greedy'){
    return(cluster_fast_greedy(graph))
  }
  # Additional clustering methods can be added here
}

# Advanced statistical measures for meta-feature computation
computeMetaFeature <- function(data, cluster, type = "mean") {
  if (type == "mean") {
    return(rowMeans(data[, cluster, drop = FALSE]))
  } else if (type == "variance") {
    return(apply(data[, cluster, drop = FALSE], 1, var))
  }
  # More statistical measures can be added here
}

# Graph-based feature selection
selectImportantFeatures <- function(graph, topN = 10) {
  centralityScores <- betweenness(graph)
  importantFeatures <- order(centralityScores, decreasing = TRUE)[1:topN]
  return(importantFeatures)
}


computeMetaFeaturesGraph <- function(data, type = 'ratio', 
                                     clusteringMethod = "louvain", percentile = 80) {
  data <- as.matrix(data)
  set.seed(1994)
  # correlationmatrix <- cor(data)
  # distancematrix <- cor2dist(correlationmatrix)
  # 
  # DM1 <- as.matrix(distancematrix)
  
  s1 <- Spectrum::CNN_kernel(data)
  s2 <- Spectrum::rbfkernel_b(data)
  s3 <- Spectrum::ng_kernel(data)
  
  s_final <- list(s1,s2,s3)
  
  DM1 <- Spectrum::integrate_similarity_matrices(s_final, diffusion_iters = 6)
  
  threshold <- correlationThreshold(DM1, percentile)
  DM1[DM1 < threshold] = 0
  
  G1 <- graph.adjacency(DM1, mode = "undirected", weighted = TRUE, diag = FALSE)
  clusterlouvain <- performClustering(G1, method = clusteringMethod)
  
  clusters <- clusterlouvain$membership
  nClust <- length(unique(clusters))
  
  meta_features <- data.frame(matrix(nrow = nrow(data), ncol = nClust))
  
  for (i in 1:nClust) {
    cluster_features <- which(clusters == i)
    meta_features[, i] <- rowMeans(data[, cluster_features, drop = FALSE])
  }
  
  colnames(meta_features) <- paste("meta_feature", 1:nClust, sep = "_")
  rownames(meta_features) <- rownames(data)
  
  importantFeatures <- selectImportantFeatures(G1)
  
  return(list(Features = meta_features, clusters = clusters, 
              ImportantFeatures = importantFeatures, graph = G1))
}


computeMetaFeaturesHC <- function(data, nClust=9, type = 'ratio', percentile = 0.80){
  
  s1 <- Spectrum::CNN_kernel(data)
  s2 <- Spectrum::rbfkernel_b(data)
  s3 <- Spectrum::ng_kernel(data)
  
  s_final <- list(s2,s3)
  
  cormat <- Spectrum::integrate_similarity_matrices(s_final, diffusion_iters = 6)
  
  # threshold <- correlationThreshold(DM1, percentile)
  # DM1[DM1 < threshold] = 0
  
  corr.min <- quantile(abs(cormat[lower.tri(cormat)]), 
                       probs = percentile)
  distmat  <- as.dist(1 - cormat)
  clustree <- hclust(distmat, method = 'ward.D2')
  clusters <- cutree(clustree, h = 1 - corr.min)
  # 
  # hc <- hclust(as.dist(DM1), method = 'ward.D2')
  
  # clusters <- dynamicTreeCut::cutreeDynamic(hc, minClusterSize = 1)
  
  nClust <- length(unique(clusters))

  # Create meta-features
  meta_features <- data.frame(matrix(nrow = nrow(data), ncol = nClust))
  
  for (i in 1:nClust) {
    cluster_features <- which(clusters == i)
    if(type == 'ratio'){
      meta_features[, i] <- rowMeans(data[, cluster_features, drop = FALSE])
    } else{
      meta_features[, i] <- rowMeans(data[, cluster_features, drop = FALSE])
    }
    
  }
  
  # Rename the meta features
  colnames(meta_features) <- paste("meta_feature", 1:nClust, sep = "_")
  
  rownames(meta_features) <- rownames(data)
  return(list(Features = meta_features, clusters = clusters))
}


computeMetaFeaturesSPECTRUM <- function(data){
  
  data <- data %>%
    as.matrix()
  # Run SIMLR
  spectrum_res <- Spectrum(data, method = 2, tunekernel = TRUE,
                           kerneltype = 'stsc', diffusion_iters = 20)
  clusters <- spectrum_res$assignments
  
  nClust <- length(unique(clusters))
  
  # Create meta-features
  meta_features <- data.frame(matrix(nrow = nrow(data), ncol = nClust))
  
  for (i in 1:nClust) {
    cluster_features <- which(clusters == i)
    meta_features[, i] <- rowMeans(data[, cluster_features, drop = FALSE])
    meta_features[, i] <- rowMeans(data[, cluster_features, drop = FALSE])
    
  }
  
  # Rename the meta features
  colnames(meta_features) <- paste("meta_feature", 1:nClust, sep = "_")
  
  rownames(meta_features) <- rownames(data)
  return(list(Features = meta_features, clusters = clusters))
}

predictMetaFeatures <- function(data, clusters, type = 'ratio'){
  nClust <- length(unique(clusters))
  message(nClust)
  data <- as.matrix(data)
  # Create meta-features
  meta_features <- data.frame(matrix(nrow = nrow(data), ncol = nClust))
  
  for (i in 1:nClust) {
    cluster_features <- which(clusters == i)
    if(type == 'ratio'){
      meta_features[, i] <- rowMeans(data[, cluster_features, drop = FALSE])
    } else{
      meta_features[, i] <- rowMeans(data[, cluster_features, drop = FALSE])
    }
    
  }
  
  colnames(meta_features) <- paste("meta_feature", 1:nClust, sep = "_")
  
  rownames(meta_features) <- rownames(data)
  
  return(meta_features)
}


computeHOPACHClusters <- function(h_res, dat, group = TRUE){
  meta_features_list <- list()
  groups <- c()
  group_list <- list()
  for(i in 2:length(h_res$cutree_list)){
    # print(i)
    cur_clusters <- h_res$cutree_list[[i]]
    nClust <- length(unique(cur_clusters))
    # print(nClust)
    meta_features <- data.frame(matrix(nrow = nrow(dat), ncol = nClust))
    for (j in 1:nClust) {
      cluster_features <- which(cur_clusters == j)
      meta_features[, j] <- rowMeans(dat[, cluster_features, drop = FALSE])
    }
    colnames(meta_features) <- paste(paste0("level_",i), paste0('cluster_',1:nClust), sep = "_")
    
    # meta_features <- dat[, cur_clusters]
    # print(dim(meta_features))
    
    meta_features_list[[i-1]] <- meta_features %>% as.data.frame()
    
    if(group == TRUE){
      groups <- c(groups, rep(paste0("level_",i), nClust))
    }
    
  }
  
  meta_features <- do.call(cbind, Filter(Negate(is.null), meta_features_list))
  
  if(group == TRUE){
    return(list(meta_features = meta_features, groups = groups, meta_list = meta_features_list))
  }
  
  return(meta_features)
}


findSequenceIndices <- function(vec) {
  # Find the differences between adjacent elements
  diff_vec <- c(1, diff(as.integer(factor(vec))))
  
  # Identify the start and end points of each sequence
  starts <- which(diff_vec != 0)
  ends <- c(starts[-1] - 1, length(vec))
  
  # Create the sequences
  sequences <- mapply(seq, from = starts, to = ends, SIMPLIFY = FALSE)
  
  return(sequences)
}

fitPLasso <- function(data, condition){
  # run Hopach
  h_res <- runHOPACH(t(data))
  
  # Get meta features from hopach result
  res_hopach <- computeHOPACHClusters(h_res = h_res, dat = data, group = TRUE)
  data <- res_hopach$meta_features %>%
    as.matrix()
  # define the groups
  groups <- findSequenceIndices(res_hopach$groups)
  
  # fit the priority lasso model
  set.seed(1234)
  pl_fit1 <- prioritylasso(X = data, Y = condition, family = "binomial", type.measure = "auc", 
                           blocks = groups, nfolds = 10,
                           standardize = FALSE, 
                           lambda.type = "lambda.min", block1.penalization = T)
  
  return(list(fit = pl_fit1, h_res = h_res))
}


corFilter <- function(data, threshold = 0.75){
  idx <- nestedcv::collinear(data, rsq_cutoff = threshold)
  if(length(idx) > 0){
    data <- data[, -c(idx)]
  }
  return(data)
}

runHOPACHCustom <- function (data, K = 10, kmax = 5, dissimilarity_metric = "cor") 
{

  dist <- as.dist(cor(t(data)))
  clustresult <- hopach(data, K = K, dmat = dist, kmax = kmax, 
                        clusters = 'none', verbose = TRUE)
  final_labels <- strsplit(as.character(format(clustresult$final$labels, 
                                               scientific = FALSE)), "")
  level_list <- list()
  for (i in seq_len(nchar(clustresult$final$labels[1]))) {
    if (i != 1) {
      level_list[[i]] <- paste(level_list[[i - 1]], unlist(lapply(final_labels, 
                                                                  "[[", i)), sep = "")
    }
    else {
      level_list[[i]] <- unlist(lapply(final_labels, "[[", 
                                       i))
    }
  }
  cutree_list <- list()
  cutree_list[[1]] <- rep(1, length(rownames(data)))
  names(cutree_list[[1]]) <- rownames(data)
  cutree_list[2:(length(level_list) + 1)] <- lapply(level_list, 
                                                    function(x) {
                                                      x <- as.numeric(as.factor(x))
                                                      names(x) <- rownames(data)
                                                      return(x)
                                                    })
  cutree_list[[length(cutree_list) + 1]] <- seq_len(length(clustresult$final$labels))
  names(cutree_list[[length(cutree_list)]]) <- rownames(data)
  return(list(cutree_list = cutree_list))
}

runEnsemble <- function(train_x, train_y, test_x){
  size <- computeGridSize(train_x)
  h_res <- runHOPACH(t(train_x), kmax = size, K = size)
  #
  res_hopach <- computeHOPACHClusters(h_res = h_res, dat = train_x, group = TRUE)
  train_x <- res_hopach$meta_features %>%
    as.matrix()
  groups <- findSequenceIndices(res_hopach$groups)
  
  test_x <- computeHOPACHClusters(h_res = h_res, 
                                  dat = test_x, group = FALSE) %>%
    as.matrix()
  
  x <- as.matrix(train_x) # Removes class
  y <- train_y # Only class
  
  ## Ridge Regression to create the Adaptive Weights Vector
  set.seed(999)
  message('Now running Adaptive Lasso')
  cv.ridge <- cv.glmnet(x, y, family='binomial', alpha=0, parallel=TRUE, standardize=TRUE)
  w3 <- 1/abs(matrix(coef(cv.ridge, s=cv.ridge$lambda.min)
                     [, 1][2:(ncol(x)+1)] ))^2 ## Using gamma = 1
  w3[w3[,1] == Inf] <- 999999999 ## Replacing values estimated as Infinite for 999999999
  
  ## Adaptive Lasso
  set.seed(999)
  cv.lasso <- cv.glmnet(x, y, family='binomial', alpha=1, parallel=TRUE, standardize=TRUE, 
                        type.measure='auc', penalty.factor=w3)
  plot(cv.lasso)
  # plot(cv.lasso$glmnet.fit, xvar="lambda", label=TRUE)
  # abline(v = log(cv.lasso$lambda.min))
  # abline(v = log(cv.lasso$lambda.1se))
  coef(cv.lasso, s='lambda.min')
  coef <- coef(cv.lasso, s='lambda.min')
  selected_attributes <- (coef@i[-1]+1) ## Considering the structure of the data frame dataF as shown earlier
  
  message('Now running Priority Lasso')
  pl_fit <- prioritylasso(X = x, Y = y, family = "binomial", type.measure = "auc", 
                          blocks = groups, nfolds = 10,
                          standardize = TRUE, lambda.type = "lambda.min", block1.penalization = T)
  
  message('Now running Adaptive Lasso')
  ipf_fit <- cvr.ipflasso(x, y,family = 'binomial', type.measure="auc",
                          standardize=TRUE, alpha = 1, pf = rep(1,length(groups)),
                          blocks=groups,nfolds=10,ncv=10)
  
  pred_test <- predict(cv.lasso, newx = as.matrix(test_x), 
                       type = "response", s = 'lambda.min')
  
  pred_test_pl <- predict(pl_fit, newdata = test_x, type = "response")
  
  pred_test_ipf <- ipflasso.predict(ipf_fit, test_x)
  
  
  preds <- data.frame(adl = pred_test[, 1], pl = pred_test_pl, 
                      ipf = pred_test_ipf$probabilitiestest)
  preds$avg <- (preds$adl + preds$pl + preds$ipf)/3
  
  return(preds)
  
}

getHCGroups <- function (hc) {
  if (!inherits(hc, "hclust")) {
    stop("hc must be an hclust object.")
  }
  nr <- nrow(hc$merge)
  ind <- 1:(nr + 1) # Initial indices of observations
  gr <- 1:(nr + 1) # Initial group assignment, each observation in its own group
  ht <- numeric(nr + 1) # Initialize the heights vector with appropriate length
  
  # Correct initialization of heights for individual observations
  # Assuming the lowest height (i.e., 0) for individual observations as they have not been merged yet
  ht[] <- 0 
  
  for (i in 1:nr) {
    indpos <- which(hc$merge[i, ] > 0)
    pos <- hc$merge[i, indpos]
    neg <- abs(hc$merge[i, ][hc$merge[i, ] < 0])
    count <- 1
    while (length(pos) > 0) {
      pos <- hc$merge[pos, ]
      neg <- c(neg, abs(pos[pos < 0]))
      pos <- pos[pos > 0]
      count <- count + 1
    }
    ind <- c(ind, sort(neg))
    gr <- c(gr, rep(i + nr + 1, length(neg)))
    ht <- c(ht, rep(hc$height[i], length(neg)))
  }
  return(list(varGroup = ind, indexGroup = gr, heights = ht))
}



MP_gLasso <- function(fit, group, lambda.type = "min", sort.type = "mean", max.shown = 20, intercept=TRUE) {
  
  ## extracting value ##
  ## extracting value - coefficients and lambda ##
  
  coef <- fit$beta.latent %>% 
    as.matrix() %>%
    as.data.frame()
  lambda <- fit$lambda
  
  # print(head(coef))
  
  if (intercept) {
    coef_name <- dimnames(coef)[[1]][-1]
    coef <- coef[-1, ]
  } else {coef_name <- dimnames(coef)[[1]]}
  
  # print(coef_name)
  
  if (mode(group) == "numeric"){
    group <- as.character(group)
  }
  
  idx <- (coef != 0)
  idx2 <- group %in% group[idx]
  
  beta <- coef[idx2]
  beta_name <- coef_name[idx2]
  
  group_name <- group[idx2]
  group_size <- as.integer(table(group_name))
  
  gname <- unique(group_name)
  gid <- as.numeric(factor(group_name))
  
  len <- length(gname)
  
  if(len == 0){
    stop('No variable with non-zero coefficient')
  }
  if(len> max.shown){ 
    if(sort.type == 'mean'){
      avgbeta_abs <- tapply(abs(beta), group_name, mean)  
      oo <- order(avgbeta_abs, decreasing = TRUE)
      gname = rownames(avgbeta_abs[oo][1:max.shown])
      idx2 = group %in% gname
      beta <- coef[idx2]
      beta_name <- coef_name[idx2]
      
      group_name <- group[idx2]
      group_size <- as.integer(table(group_name))
      
      gname <- unique(group_name)
      gid <- as.numeric(factor(group_name))
      
      len <- length(gname)
    }else{
      max.beta_abs <- tapply(abs(beta), group_name, max)  
      oo <- order(max.beta_abs, decreasing = TRUE)
      gname = rownames(max.beta_abs[oo][1:max.shown])
      idx2 = group %in% gname
      beta <- coef[idx2]
      beta_name <- coef_name[idx2]
      
      group_name <- group[idx2]
      group_size <- as.integer(table(group_name))
      
      gname <- unique(group_name)
      gid <- as.numeric(factor(group_name))
      
      len <- length(gname)
    }
  }
  
  if (sort.type == "mean") {
    avgbeta_abs <- tapply(abs(beta), group_name, mean)  
    oo <- order(avgbeta_abs, decreasing = TRUE)
    scale <- tapply(abs(beta), group_name, mean) / tapply(abs(beta), group_name, max)
    
    data <- data.frame(beta=beta, group_name=factor(gid), beta_name = beta_name, groups = group_name)
    data$sign <- ifelse(data$beta > 0, "beta \nwith positive sign", "beta \nwith negative sign")
    
    p2 <- data %>%
      dplyr::group_by(group_name) %>%
      summarise(avg_beta = mean(abs(beta)), gsize=n()) %>%
      mutate(group_name = fct_reorder(group_name, avg_beta, .desc=TRUE),
             group_size = fct_reorder(factor(gsize), gsize, .desc=FALSE)) %>%
      ggplot() +
      geom_bar(aes(x=group_name, y=avg_beta, fill=group_size), col="black",
               stat="identity", width = 1) +
      theme_light() +
      scale_fill_grey(start=1, end=0.5) +
      scale_x_discrete(breaks=factor(unique(gid)), labels=gname) +
      coord_polar()
  } else {
    max.beta_abs <- tapply(abs(beta), group_name, max)  
    oo <- order(max.beta_abs, decreasing = TRUE)
    scale <- rep(1, length(unique(group_name)))
    
    data <- data.frame(beta=beta, group_name=factor(gid), beta_name = beta_name, group = group_name)
    data$sign <- ifelse(data$beta > 0, "beta \nwith positive sign", "beta \nwith negative sign")
    
    p2 <- data %>%
      dplyr::group_by(group_name) %>%
      dplyr::summarise(max_beta = max(abs(beta)), gsize=n()) %>%
      dplyr::mutate(group_name = fct_reorder(group_name, max_beta, .desc=TRUE),
             group_size = fct_reorder(factor(gsize), gsize, .desc=FALSE)) %>%
      ggplot() +
      geom_bar(aes(x=group_name, y=max_beta, fill=group_size), col="black",
               stat="identity", width = 1) +
      theme_light() +
      scale_fill_grey(start=1, end=0.5) +
      scale_x_discrete(breaks=factor(unique(gid)), labels=gname) +
      coord_polar()
  }
  
  position <- list(length = len)
  for (i in 1:len) {
    position[[i]] <- jitter(rep(i, group_size[oo][i]))
  }
  
  beta_scale <- list(length = len)
  for (i in 1:len) {
    beta_scale[[i]] <- split(abs(beta), group_name)[[i]] * scale[i]
  }
  
  beta_sign <- list(length = len)
  for (i in 1:len) {
    beta_sign[[i]] <- split(data$sign, group_name)[[i]]
  }
  
  beta_origin <- list(length = len)
  for (i in 1:len) {
    beta_origin[[i]] <- split(beta, group_name)[[i]]
  }
  
  indexx <- list(length = len)
  for (i in 1:len){
    indexx[[i]] <- split(seq(sum(group_size)), group_name)[[i]]
  }
  
  b_name <- list(length = len)
  for (i in 1:len) {
    b_name[[i]] <- split(beta_name, group_name)[[i]]
  }
  
  sample_points <- list(length= len)
  for (i in 1:len) {
    sample_points[[i]] <- data.frame(position=position[[i]], beta_scale=beta_scale[[oo[i]]],
                                     beta_sign=beta_sign[[oo[i]]], index=indexx[[oo[i]]], coef=beta_origin[[oo[i]]],
                                     name=b_name[[oo[i]]])
  }
  
  for (i in 1:len) {
    sample_points[[i]]$tt <- paste0("Name : ", sample_points[[i]]$name,
                                    "\nCoefficient : ", round(sample_points[[i]]$coef,4))
  }
  
  p3 <- p2
  for(i in 1:len){
    loop_input <- paste0('geom_point_interactive(data=sample_points[[i]], aes(x=position, y=beta_scale,
                       shape=beta_sign, colour=beta_sign, tooltip=tt, data_id=index))')
    p3 <- p3 + eval(parse(text = loop_input))
  }
  
  p4 <- p3 +
    scale_shape_manual(name = 'sign of beta', values=c(1, 17)) +
    scale_colour_manual(name = 'sign of beta', values = c('black', 'black', 'black')) +
    labs(title = "Group Lasso",
         subtitle = paste0("Lambda : ", round(lambda, 3), ", Lambda type : ",lambda.type,"\nSort type : ", sort.type),
         x="", y="", fill="Group size") +
    theme(plot.title = element_text(hjust = 0.5,size=13, face="bold"), 
          legend.key = element_rect(colour = "black"),
          plot.subtitle = element_text(size=11, face="bold")) +
    guides(fill = guide_legend(ncol=2)) +
    scale_y_continuous(labels = NULL)
  
  # girafe(code = print(p4), width_svg = 8, height_svg = 8)
  
  return(list(data = data, plot = p4))
}

bootstrapHclustCustom <- function(X, frac = 1, B = 50, method = "ward.D2", nCore = NULL) {
  t1 <- proc.time()
  n <- nrow(X)
  
  if (frac <= 0 | frac > 1) {
    stop("frac must be between 0 and 1.")
  }
  
  
  nInd <- floor(n * frac)
  d <- 0
  for (i in 1:B)
  {
    ind <- sample(n, nInd, replace = TRUE)
    d <- d + parDist(t(X[ind, ]), threads = nCore)
  }
  
  d <- d / B
  
  hc <- fastcluster::hclust(d, method = ifelse(is.character(method), method, "ward.D2"))
  t2 <- proc.time()
  tcah <- t2 - t1
  
  hc$tcah <- as.numeric((t2 - t1)[3])
  
  return(hc)
}


getClusterTreeCustom <- function(exprs,
                           clusters, hierarchy_method = 'average',
                           scale_exprs=TRUE) {
  clust_med_dt <- as.data.table(exprs)
  clust_med_dt[, cluster_id := clusters]
  # data table containing median
  res <- clust_med_dt[, lapply(.SD, median, na.rm=TRUE), by=cluster_id]
  res2 <- res[,.SD, .SDcols = !c('cluster_id')]
  rownames(res2) <- res[["cluster_id"]]
  
  if (scale_exprs) {
    res_unscaled <- res2
    res2[, (colnames(res2)) := lapply(.SD, scale), .SDcols=colnames(res2)]
  } else {
    res_unscaled <- res2
  }
  
  clust_dist <- as.dist(1 - cor(t(res2)))
  hc_dend <- hclust(clust_dist, method=hierarchy_method)
  hc_dend$labels <- as.character(res[["cluster_id"]])
  return(list(
    median_freq = res_unscaled,
    clust_tree = hc_dend
  ))
}


computeClusters <- function(groups, dat, group = TRUE){
  meta_features <- data.frame(matrix(nrow = nrow(dat), ncol = length(groups)))
  for(i in 1:length(groups)){
    cur_clusters <- groups[[i]]
    nClust <- length(unique(cur_clusters))
    meta_features[, i] <- rowMeans(dat[, cur_clusters, drop = FALSE])
    
    
  }
  return(meta_features)
}



# Assuming DipValue is a function that computes the dip statistic for a numeric vector
DipValue <- function(x) {
  # Placeholder for actual Dip Value computation
  # This should be replaced with actual computation or function call
  return(dip(x))  # Simplified example, replace with real computation
}

computeDipValues <- function(D){
  dim <- ncol(D)  # Number of dimensions/columns in D
  dipValues <- numeric(dim)  # To store dip values for each dimension
  
  # Compute Dip Value for each dimension
  for (i in 1:dim) {
    dipValues[i] <- DipValue(D[, i])
  }
  
  return(dipValues)
}

DipScale <- function(D, dipValues) {
  dim <- ncol(D)
  # Rescale values of each axis
  for (i in 1:dim) {
    D[, i] <- (D[, i] / max(D[, i])) * dipValues[i]
  }
  
  return(D)  # Return modified data
}

# Function to perform overlapping group Lasso regularization on input data
# Arguments:
#   - data: Numeric matrix of features for the proportion data
#   - data_logit: Numeric matrix of features for the logit-transformed proportion data
#   - data_mean: Numeric matrix of features for the gene means per cell type
# Returns:
#   - List containing the fitted model and scaling parameters

runOverlapLasso <- function(data, data_logit, data_mean, nclust = 30, numMarkers = 27,
                            alpha = 0.5, linkage = 'average', seed = 1994,
                            clinical = F, clinicalData = NULL, root = F){
  
  # Compute pairwise Euclidean distances for proportion data
  d_prop <- parallelDist::parallelDist(t(data), method = 'euclidean')
  # d_prop <- computeDistances(data)
  # Hierarchical clustering of primary data
  hc_prop <- hclust(d_prop, method = linkage)
  # Extract clusters from the dendrogram
  hc_prop <- treekoR:::findChildren(ggtree(as.phylo(hc_prop), 
                                           ladderize = F, layout = 'dendrogram'))
  # Extract clusters and assign group numbers
  subs_prop <- hc_prop$data$clusters
  prop_dat <- hc_prop$data
  prop_dat$group <- 1:length(subs_prop)
  
  # Compute pairwise Euclidean distances for logit-transformed proportion data
  if(root){
    hc_logit <- treekoR:::findChildren(ggtree(as.phylo(hc_logit),
                                              ladderize = F, layout = 'dendrogram'))
  } else {
    d_logit <- parallelDist::parallelDist(t(data_logit), method = 'euclidean')
    hc_logit <- hclust(d_logit, method = linkage)
    hc_logit <- treekoR:::findChildren(ggtree(as.phylo(hc_logit),
                                              ladderize = F, layout = 'dendrogram'))
  }
  
  # Extract clusters and assign group numbers
  subs_logit <- hc_logit$data$clusters
  logit_dat <- hc_logit$data
  logit_dat$group <- 1:length(subs_logit)
  
  # Compute pairwise Euclidean distances for gene means per cell type data
  d_mean <- parallelDist::parallelDist(t(data_mean), method = 'euclidean')
  # d_mean <- computeDistances(data_mean)
  # Hierarchical clustering of gene means per cell type data
  hc_mean <- hclust(d_mean, method = linkage)
  # Extract clusters from the dendrogram
  hc_mean <- treekoR:::findChildren(ggtree(as.phylo(hc_mean),
                                           ladderize = F, layout = 'dendrogram'))
  # Extract clusters and assign group numbers
  subs_mean <- hc_mean$data$clusters
  mean_dat <- hc_mean$data
  mean_dat$group <- 1:length(subs_mean)
  
  # Prepare variable groups for overlapping group Lasso
  var_groups<- list()
  n <- length(subs_logit)
  # n2 <- 2*n
  for(i in 1:length(subs_logit)){
    var_groups[[i]] <- which(colnames(data_logit) %in% subs_logit[[i]])
    # var_groups[[i+n]] <- which(colnames(data_logit) %in% subs_logit[[i]]) + nclust
  }
  
  # for(i in 1:length(subs_mean)){
  #   var_groups[[i+n]] <- which(colnames(data_mean) %in% subs_mean[[i]]) + nclust
  # }
  
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
  
  # var_groups[off_set + 1] <- list(seq(last_max + 1, ncol(markerMeanCellType) + last_max))
  
  if(clinical){
    n2 <- length(var_groups)
    idx <- ncol(data_logit) + ncol(data_mean)
    var_groups[[n2+1]] <- c(idx+1, idx+2)
  }
  
  # print(var_groups)
  # Combine all input data matrices
  if(clinical){
    print(idx)
    X_train <- cbind(data_logit, data_mean, clinicalData) 
    # Scale the combined data
    scaleVals <- preProcess(X_train[, 1:idx], method = c('range'))
    X_train[, 1:idx] <- predict(scaleVals, X_train[, 1:idx]) 
    X_train <- X_train %>%
      as.matrix()
  } else{
    X_train <- cbind(data_logit, data_mean) 
    # Scale the combined data
    scaleVals <- preProcess(X_train, method = c('range'))
    X_train <- predict(scaleVals, X_train) %>%
      as.matrix()
  }
  # colnames(X_train) <- janitor::make_clean_names(colnames(X_train))
  
  # Extract response variable
  y_train <- condition %>%
    as.numeric()
  
  # Cross-validation for overlapping group Lasso
  set.seed(999)
  cvfit <- cv.grpregOverlap(X_train, y_train, var_groups, penalty = 'cMCP', 
                            family = 'binomial',  alpha = alpha, 
                            nfolds = 10, seed = seed)
  # Plot cross-validation results
  plot(cvfit)
  
  # Fit overlapping group Lasso model
  fit <- grpregOverlap(X_train, y_train, var_groups,
                       family='binomial', alpha = alpha,
                       returnX.latent = T, returnOverlap = FALSE, 
                       lambda = cvfit$lambda.min, penalty = 'cMCP')
  
  # print(colnames(X_train))
  # Return fitted model and scaling parameters
  return(list(fit = fit, scaling = scaleVals, data = X_train,
              hc_prop = hc_prop, hc_logit = hc_logit, hc_mean = hc_mean))
}


# Function to plot a dendrogram with annotated coefficients for overlapping group Lasso
# Arguments:
#   - hc: Hierarchical clustering object
#   - data: Data frame containing coefficient information
#   - type: Type of coefficients to plot ('prop' for proportional, default is 'prop')
#   - heatmap: Boolean indicating whether to include a heatmap (default is FALSE)
# Returns:
#   - ggplot object representing the dendrogram with annotated coefficients

plotOverlapTreeOld <- function(hc, data, coefs, type = 'logit', 
                            heatmap = FALSE, offset = NULL){
  # Prepare coefficient data
  res_dat <- data
  # res_dat$beta_name <- janitor::make_clean_names(res_dat$beta_name)
  # print((res_dat$beta_name))
  # Extract groups from hierarchical clustering object
  hc_dat <- hc$data

  # hc_dat$label <- janitor::make_clean_names(hc_dat$label)
  # print((hc_dat$label))
  hc_dat$group <- hc_dat$node
  res_dat_new <- res_dat[res_dat$beta_name %like% type, ]
  res_dat_new$group <- as.numeric(res_dat_new$group)
  
  # if we are dealing with mean expression. Offset the 
  # points to fix naming
  if(type == 'mean'){
    res_dat_new$group <- res_dat_new$group - offset
  }
  
  res_dat_new <- base::merge(hc_dat, res_dat_new, all = T) %>%
    dplyr::select(-c('beta_name', 'group_name', 'sign'))
  res_dat_new$beta <- if_else(is.na(res_dat_new$beta), 0, res_dat_new$beta)
  res_dat_new_temp <- res_dat_new %>%
    dplyr::select(group, beta) %>%
    abs() %>% # make the internal beta values positive
   dplyr::group_by(group) %>%
    dplyr::summarise_at(vars(beta), list(mean))
  
  res_dat_new <- base::merge(res_dat_new %>%
                               dplyr::select(-beta), res_dat_new_temp)
  res_dat_new <- res_dat_new[!duplicated(res_dat_new), ]
  
  # Prepare coefficient names
  coefs_prop <- coefs[coefs$feature %like% type, ] 
  colnames(coefs_prop) <- c('beta_leaf', 'label')
  
  # Merge coefficient data with hierarchical clustering data
  res_dat_new <- base::merge(res_dat_new, coefs_prop, all = T)
  
  # Adjust coefficients for internal and leaf nodes
  res_dat_new$beta <- if_else(res_dat_new$isTip == T, 0, res_dat_new$beta)
  res_dat_new$beta_leaf <- if_else(is.na(res_dat_new$beta_leaf), 0, res_dat_new$beta_leaf)
  res_dat_new$beta_final <- res_dat_new$beta + res_dat_new$beta_leaf
  res_dat_new$beta_internal <- res_dat_new$beta
  
  # Categorize coefficients based on their values
  res_dat_new <- res_dat_new %>%
    mutate(beta = case_when(
      beta_final < 0 ~ "Leaf < 0",
      beta_final == 0 ~ "Leaf/Parent = 0",
      beta_final > 0 ~ "Leaf > 0",
    ))
  res_dat_new$beta <- if_else(res_dat_new$beta_final > 0 & res_dat_new$isTip == F, 
                              "Parent > 0", res_dat_new$beta)
  
  # Update hierarchical clustering data with adjusted coefficient information
  res_dat_new$label <- str_replace(res_dat_new$label, paste0('_',type), '')
  res_dat_new$label <- str_replace(res_dat_new$label, '_', " ")
  hc$data <- res_dat_new
  # print(res_dat_new)
  
  
  # Create ggplot object for dendrogram visualization
  p <- ggtree(hc$data, layout = 'dendrogram', hang = 0) +
    geom_tiplab(as_ylab = T, geom = "text", size = 24, color = 'black') +
    geom_nodepoint(aes(subset = beta_internal != 0, color = beta_internal), size = 5) +
    geom_nodepoint(aes(color = beta_internal), size = 0) +
    scale_color_gradient2(mid = "grey", high = muted("purple"), mid = "white",
                          midpoint=0, name = "Beta Internal") +
    new_scale_color() +
    geom_tippoint(aes(subset = beta_leaf != 0, color = beta_leaf), size = 5) + 
    scale_color_gradient2(low=muted("green"), high=muted("orange"), mid = "white",
                          midpoint=0, name = "Beta Leaf") +
    geom_tippoint(aes(color = beta_leaf), size = 0)
  
  # Add heatmap if specified
  if(heatmap){
    p <- p +
      gheatmap(p, data = t(dat_prop), offset = 0.005, hjust = 0,
               low = muted("blue"), high = muted("darkred"), color = 'white',
               legend_title="Proportions") +
      theme(legend.title = element_text(color = "black", size = 10),
            legend.text = element_text(color = "black", size = 10))
  }
  
  return(p)
}


# Define the function
plotOverlapTree <- function(model_fit, training_data, tree_data, nodes_to_remove) {
  
  # Function to sort each sublist
  sort_sublist <- function(x) {
    return(sort(x))
  }
  
  # Predict the variables from the model
  predicted_vars <- predict(model_fit, type = "groups")
  
  # Extract coefficients and filter the relevant features
  coefficients <- coef(model_fit) %>%
    as.matrix() %>%
    as.data.frame() %>%
    dplyr::mutate(feature = rownames(.)) %>%
    # dplyr::filter(V1 != 0) %>%
    dplyr::filter(feature %like% '_logit') %>%
    dplyr::filter(!feature %like% 'Intercept') 
  
  colnames(coefficients)[[1]] <- 'coefficient_value'
  
  # Initialize lists for groups, clusters, and beta values
  group_list <- list()
  clusters_list <- list()
  beta_values <- list()
  
  # Loop through each prediction
  for(i in 1:length(predicted_vars)){
    current_features <- model_fit$group[[predicted_vars[[i]]]]
    
    group_name <- paste0('group_', i)
    cluster_names <- colnames(training_data)[current_features]
    
    if(length(current_features) > 1){
      beta_value <- coefficients[colnames(training_data)[current_features], 1] %>% abs() %>% mean
    } else {
      beta_value <- coefficients[colnames(training_data)[current_features], 1]
    }
    
    group_list[[i]] <- group_name
    clusters_list[[i]] <- cluster_names
    beta_values[[i]] <- beta_value
  }
  
  # Combine groups, clusters, and beta values into a data frame
  result_data <- cbind(group_list, clusters_list, beta_values) %>%
    as.data.frame() %>%
    dplyr::filter(clusters_list %like% "logit")
  result_data$beta_values <- unlist(result_data$beta_values)
  
  # Modify tree data
  tree_cluster_data <- tree_data$data
  tree_cluster_data$clusters_list <- lapply(tree_cluster_data$clusters, function(x) gsub("_", " ", x))
  tree_cluster_data$clusters_list <- lapply(tree_cluster_data$clusters_list, function(x) gsub(" logit", "", x))
  
  # Function to remove the specified nodes from each list
  remove_nodes <- function(x, nodes) {
    return(setdiff(x, nodes))
  }
  
  # Apply the node removal and sorting to clusters
  tree_cluster_data$clusters_list <- lapply(tree_cluster_data$clusters_list, remove_nodes, nodes = nodes_to_remove)
  tree_cluster_data$clusters_list <- lapply(tree_cluster_data$clusters_list, sort_sublist)
  
  # Apply the modifications to the result clusters
  result_data$clusters_list <- lapply(result_data$clusters_list, function(x) gsub("_logit", "", x))
  result_data$clusters_list <- lapply(result_data$clusters_list, function(x) gsub("_", " ", x))
  result_data$clusters_list <- lapply(result_data$clusters_list, sort_sublist)
  
  # if we have an empty model
  if(nrow(result_data) < 1){
    message("No non-zero proportions.. Will generate null tree")
    result_data <- tree_cluster_data
    result_data$beta_internal <- rep(0, nrow(result_data))
    result_data$beta_leaf <- rep(0, nrow(result_data))
  } else{
    # Merge results with tree data and compute additional metrics
    result_data <- result_data %>%
      dplyr::full_join(tree_cluster_data, by = 'clusters_list') %>%
      dplyr::mutate(beta_internal = if_else(isTip == FALSE, beta_values, 0)) %>%
      dplyr::mutate(beta_leaf = if_else(isTip == TRUE, beta_values, 0)) %>%
      mutate_at(vars(beta_values, beta_leaf, beta_internal), ~replace_na(., 0))
  }
  
  
  result_data$label <- str_replace(result_data$label, paste0('_','logit'), '')
  result_data$label <- str_replace(result_data$label, '_', " ")
  
  # Create ggplot object for dendrogram visualization
  plot <- ggtree(result_data, layout = 'dendrogram', hang = 0) +
    geom_tiplab(as_ylab = T, geom = "text", size = 24, color = 'black') +
    geom_nodepoint(aes(subset = beta_internal != 0, color = beta_internal), size = 10) +
    geom_nodepoint(aes(color = beta_internal), size = 0) +
    scale_color_gradient2(mid = "grey", high = muted("purple"),
                          midpoint = 0, name = "Beta Internal",
                          guide = guide_colorbar(order = 2)) +
    new_scale_color() +
    geom_tippoint(aes(subset = beta_leaf != 0, color = beta_leaf), size = 10) + 
    scale_color_gradient2(low = muted("green"), high = muted("orange"), mid = "grey",
                          midpoint = 0, name = "Beta Leaf",
                          guide = guide_colorbar(order = 3)) +
    geom_tippoint(aes(color = beta_leaf), size = 0)
  
  
  return(plot)
}
# Example usage of the function
# p <- process_tree_data(fit, X_train, tree, nodes_to_delete)


cyCombineNormalize <- function(df,
                      markers = NULL,
                      norm_method = "scale",
                      ties.method = "average") {
  
  # Remove case-sensitivity
  norm_method <- norm_method %>% stringr::str_to_lower()
  ties.method <- ties.method %>% stringr::str_to_lower()
  
  # Error check
  if (norm_method == "rank" && ties.method %!in% c("average", "first", "last", "random", "max", "min")) {
    stop("When using norm_method = 'rank', please use an available ties.method (average, first, last, random, max, or min).")
  }
  
  # Messaging
  if (norm_method == "rank") {message("Ranking expression data..")
  } else if (norm_method == "scale") {message("Scaling expression data..")
  } else if (norm_method == "qnorm") {
    # message("Quantile normalizing expression data..")
    # Run quantile normalization
    df_normed <- cyCombine:::quantile_norm(df, markers = markers)
    return(df_normed)
  } else stop("Please use either 'scale', 'rank', or 'qnorm' as normalization method." )
  
  # Scale or rank at marker positions individually for every batch
  df_normed <- df %>%
    dplyr::group_by(.data$sample_id) %>%
    purrr::when(
      norm_method == "rank"  ~ dplyr::mutate(
        ., dplyr::across(dplyr::all_of(markers),
                         .fns = ~ {
                           if(sum(.x) == 0) stop("A marker is 0 for an entire batch. Please remove this marker.")
                           rank(.x, ties.method = ties.method) / length(.x)})),
      norm_method == "scale" ~ dplyr::mutate(., dplyr::across(dplyr::all_of(markers),
                                                              .fns = ~{
                                                                if(sum(.x) == 0) stop("A marker is 0 for an entire batch. Please remove this marker.")
                                                                scale(.x)}))
    ) %>%
    dplyr::ungroup()
  return(df_normed)
}



generateHeatmapMatrix <- function(data, markers = NULL, clusters = NULL, threshold = 2,
                          clusterMarkers = FALSE, fontSize = 14) {
  
  # do some house keeping
  if (is.null(clusters)) {
    stop("Please provide a vector of cluster labels")
  }
  
  if (!(is(data, "matrix") || is(data, "data.frame"))) {
    stop("Make sure you are passing in a dataframe or matrix")
  }
  
  
  # again if no markers are given, make sure all the columns are numeric
  if (is.null(markers)) {
    numNumeric <- sum(apply(data, 2, function(x) is.numeric(x)))
    if (numNumeric != ncol(data)) {
      stop(
        "If markers of interest are not provided, ",
        "make sure the data contains all numeric columns"
      )
    }
    message("No markers provided, will be using all columns as markers")
    markers <- colnames(data)
  }
  
  features <- data[, markers]
  # do some wrangling to get it in the proper format
  featuresHeatmap <- aggregate(
    . ~ as.character(clusters),
    features[, markers],
    mean
  )
  rownames(featuresHeatmap) <- featuresHeatmap[, 1]
  featuresHeatmap <- featuresHeatmap[, -1]
  
  # compute the marker expression
  featuresHeatmap <- sweep(featuresHeatmap, 2, colMeans(featuresHeatmap), "-")
  featuresHeatmap <- sweep(
    featuresHeatmap, 2, apply(featuresHeatmap, 2, sd), "/"
  )
  featuresHeatmap[featuresHeatmap > threshold] <- threshold
  featuresHeatmap[featuresHeatmap < -threshold] <- -threshold
  
  # compute the heatmap annotations
  annotationRow <- data.frame(Clusters = rownames(featuresHeatmap))
  
  rn <- rownames(featuresHeatmap)
  featuresHeatmap <- as.matrix(featuresHeatmap)
  rownames(featuresHeatmap) <- rn
  rownames(annotationRow) <- rownames(featuresHeatmap)
  
  gapRows <- which(!duplicated(substr(rownames(featuresHeatmap), 1, 2)))[-1] - 1
  
  # pHeat <- ggplotify::as.ggplot(pheatmap(featuresHeatmap,
  #                                        gaps_row = gapRows,
  #                                        annotation_row = annotationRow, annotation_legend = FALSE,
  #                                        cluster_cols = clusterMarkers, cluster_rows = FALSE,
  #                                        fontsize = fontSize
  # ))
  # return(pHeat)
  
  return(featuresHeatmap)
}


closest_power_of_2_less <- function(number) {
  power_of_2 <- 1
  while (number >= 2) {
    number <- number / 2
    power_of_2 <- power_of_2 * 2
  }
  return(power_of_2)
}


computeDistances <- function(data){
  pear <- stats::cor((data), method = "pearson")
  cosi <- coop::cosine(data)
  spear <- stats::cor((data), method = "spearman")
  
  dMat <- fuse(cor2dist(pear), cor2dist(cosi), cor2dist(spear))
  
  return(dMat)
}

plotPredProps <- function(data, truth, predicted, method){
  
  conf_matrix <- tibble("Truth" = truth,
                        "Predicted" = predicted) %>%
    table()
  
  data <- as.data.frame(melt(conf_matrix))
  names(data) <- c("Truth", "Predicted", "Count")
  
  # Calculate the total support for each actual class
  total_support <- aggregate(Count ~ Truth, data = data, sum)
  
  # Merge total support back to the original data for proportions
  data <- merge(data, total_support, by = "Truth", all.x = TRUE)
  
  names(data)[names(data) == "Count.x"] <- "Count"
  names(data)[names(data) == "Count.y"] <- "TotalSupport"
  
  # Calculate proportion
  data$Proportion <- data$Count / data$TotalSupport
  
  # Creating the plot
 p.prop <- ggplot(data, aes(x = Truth, y = Proportion, fill = Predicted)) +
    geom_bar(stat = "identity") +
    labs(title = paste0(method, " - Class Support and Prediction Proportions"),
         x = "Truth",
         y = "Proportion of Total Support",
         fill = "Predicted") +
    theme_minimal() + 
    theme_bw() +
    theme(axis.text.x = element_text(size = 12),
          plot.title = element_text(size = 15, hjust = 0.5),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12)
    ) + scale_color_d3(palette = 'category20') + 
    scale_fill_startrek()
 
 return(p.prop)
}

customNorm <- function(x){
  x <- abs(x)
  x <- log(x + sqrt(x^2 + 1))
  
  return(x)
}

# Function to compute MAD
mad_scaled <- function(x) {
  mad <- median(abs(x - median(x)))
  mad_scaled <- 1.4826 * mad
  return(mad_scaled)
}


robust_scaler <- function(x) {
  Q1 <- apply(x, 2, quantile, 0.25)
  Q3 <- apply(x, 2, quantile, 0.75)
  IQR <- Q3 - Q1
  scaled <- sweep(x, 2, Q1, FUN = "-")
  scaled <- sweep(scaled, 2, IQR, FUN = "/")
  return(scaled)
}

quantile_transformer <- function(x) {
  apply(x, 2, function(col) {
    qnorm(rank(col, ties.method = "random") / (length(col) + 1))
  })
}

compute_weights <- function(x_train, y_train, groups){
  
  cv.ridge <- cv.grpregOverlap(x_train, y_train, family='binomial', alpha=1e-16, group = groups,
                               nfolds = 10, seed = 999, penalty = 'cMCP')
  
  coefs <- as.numeric(coef(cv.ridge, s = cv.ridge$lambda.min))[-1]
  
  # Calculate group weights
  grp_weights <- numeric(length = length(groups))
  # Assuming 'index' is a list of lists, each containing indices of coefs
  for (i in seq_along(groups)) {
    # Extract coefficients for the current group using indices from the sublist
    current_coefs <- coefs[groups[[i]]]
    
    val <- 1 / sqrt(sum(abs(current_coefs^2)))
    # Calculate the group weight as the inverse of the square root of the sum of squares of coefficients in this group
    if(val == Inf){
      grp_weights[i] <- 1e-16
    } else{
      grp_weights[i] <- val
    }
    
  }
  
  return(grp_weights)
}

robustscale_fit <- function(data, dim=2, center=TRUE, scale=TRUE, preserveScale = TRUE) {
  medians <- NULL
  mads <- NULL
  
  if(center) {
    medians <- apply(data, dim, median, na.rm = TRUE)
  }
  
  if(scale) {
    if (is.null(medians)) {
      adjusted_data <- data
    } else {
      adjusted_data <- sweep(data, dim, medians, "-")
    }
    
    mads <- apply(adjusted_data, dim, mad, na.rm = TRUE)
    
    if(preserveScale) {
      mads <- mads / mean(mads, na.rm = TRUE)
    }
  }
  
  return(list(medians=medians, mads=mads))
}

robustscale_transform <- function(data, params, dim=2) {
  if (!is.null(params$medians)) {
    data <- sweep(data, dim, params$medians, "-")
  }
  if (!is.null(params$mads)) {
    data <- sweep(data, dim, params$mads, "/")
  }
  return(data)
}


computeElbow <- function(vals) {
  diffs <- diff(vals)
  # diffs <- diffs[-1]
  optKb <- which.max(abs(diffs)) + 1
  return(optKb)
}

process_data <- function(df, sample_ids, useMarkers, celltypes) {
  df_sub <- df[df$sample_id %in% sample_ids, ] %>%
    dplyr::filter(CellTypes %in% celltypes)
  clusters <- df_sub$CellTypes
  ### Calculate the median expression
  expr_median <- data.frame(df_sub[, useMarkers], cell_clustering = clusters) %>%
    group_by(cell_clustering) %>% summarize_all(list(~median(.)))
  
  df_med = as.data.frame(dplyr::select(expr_median, -cell_clustering))
  # df_med = scale(df_med)
  rownames(df_med) = paste('c',expr_median$cell_clustering,sep ='')
  ## Sort the cell clusters with hierarchical clustering
  d <- as.dist(1 - cor(t(df_med))) %>% 
    replace(is.na(.), 0)
  cluster_rows <- as.dendrogram(hclust(d, method = "ward.D2"))
  #cluster_rows  <- cluster_rows %>% set("nodes_pch", 19) 
  d <- as.dist(1 - cor((df_med))) %>% 
    replace(is.na(.), 0)
  hc <- hclust(d, method = "ward.D2")
  lineage_markers_ord <- hc$labels[hc$ord]
  lineage_markers_ord <- unlist(lapply(strsplit(lineage_markers_ord,'_'),function(x)x[1]))
  # colnames(exprs2) = unlist(lapply(strsplit(colnames(exprs2),'_'),function(x)x[1]))
  Cells = data.frame(df_sub[, useMarkers], clusters = clusters, 
                     group = df_sub$gensini_bin,sample = df_sub$sample_id, 
                     subject  = df_sub$sample_id) 
  dataCell = data.frame(cells = as.numeric(table(Cells$clusters)), row.names = names(table(Cells$clusters)))
  phylo <- as.phylo(cluster_rows)
  #gTree <- ggtree(phylo, layout = "circular", branch.length = 'none')%>%
  gTree <- ggtree(phylo, branch.length = 'none') 
  
  clustOrd <- as.character(na.omit(gTree$data$label[order(gTree$data$y, decreasing = TRUE)]))
  return(list(df_sub = df_sub, gTree = gTree, clustOrd = clustOrd, clusters = clusters))
}

plot_clustering_distr_wrapper <- function(expr, cell_clustering,clustOrd, 
                                          pval = NULL, title = ""){
  # Calculate the median expression
  cell_clustering <- factor(cell_clustering)
  expr_median <- data.frame(expr, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% summarize_all(funs(median))
  # Sort the cell clusters with hierarchical clustering
  #  d <- dist(expr_median[, colnames(expr)], method = "euclidean")
  #  cluster_rows <- hclust(d, method = "average")
  # Calculate cluster frequencies
  freq_clust <- table(cell_clustering)
  nm = names(freq_clust)
  freq_clust <- round(as.numeric(freq_clust)/sum(freq_clust)*100, 2)
  names(freq_clust) = nm
  # cell_clustering <- factor(cell_clustering
  #   labels = paste0(levels(cell_clustering), " (", freq_clust, "%)"))
  
  cell_clustering <- factor(cell_clustering,
                            labels = paste0(clustOrd, " (f = ", freq_clust[clustOrd], "%)"))
  
  # ### Data organized per cluster
  ggd <- melt(data.frame(cluster = cell_clustering, expr),
              id.vars = "cluster", value.name = "expression",
              variable.name = "antigen")
  ggd$antigen <- factor(ggd$antigen, levels = colnames(expr))
  ggd$reference <- "no"
  ### The reference data
  ggd_bg <- ggd
  ggd_bg$cluster <- "reference"
  ggd_bg$reference <- "yes"
  ggd_plot <- rbind(ggd, ggd_bg)
  ggd_plot$cluster <- factor(ggd_plot$cluster,
                             # levels = c(levels(cell_clustering)[rev(cluster_rows$order)], "reference"))
                             levels = c(rev(levels(cell_clustering)), "reference"))
  ggplot() +
    geom_density_ridges(data = ggd_plot, aes(x = expression, y = cluster,
                                             color = reference, fill = reference), alpha = 0.3) +
    facet_wrap( ~ antigen, scales = "free_x", nrow = 2) +
    labs(title = title) +
    theme_ridges() +
    theme(axis.text = element_text(size = 7),
          strip.text = element_text(size = 7), 
          legend.position = "none",
          plot.title = element_text(size = 20, hjust = 0.5))
}



fitVAES <- function(train_dat, test_dat, useMarkers, kRange = c(9,10,11),
                    dropout_rate = 0,
                    batch_size = 32L, latent_dim = 16L, epochs = 50){
  
  train_dat_norms <- matrix(nrow = nrow(train_dat), ncol = length(useMarkers), 
                            data = 0)
  test_dat_norms <- matrix(nrow = nrow(test_dat), ncol = length(useMarkers), 
                            data = 0)
  
  train_dat_encoded <- matrix(nrow = nrow(train_dat), ncol = latent_dim, 
                            data = 0)
  test_dat_encoded <- matrix(nrow = nrow(test_dat), ncol = latent_dim, 
                           data = 0)
  
  message("Now finding most stable sample")
  res_ind <- ComputeReferenceSample(train_dat, useMarkers, N = 2)
  
  topNsamples <- res_ind$topNSamples
  train_dat_sub <- train_dat[which(train_dat$sample_id %in% res_ind$topNSamples), ]
  for(i in 1:length(kRange)){
    
    # sample <- topNsamples[[i]]
    k <- kRange[[i]]
    message(paste("Now running for cluster:", k))
    cell_types <- kmeans(train_dat_sub[, useMarkers], centers = k)$cluster
    
    
    message("Now fitting VAE model")
    
    vae_model <- train_vae_model(train_dat = train_dat_sub[, useMarkers], 
                                       batch_size = batch_size, 
                                       cell_types = cell_types,
                                       dropout_rate = dropout_rate,
                                       useMarkers = useMarkers, 
                                       epochs = epochs, 
                                       latent_dim = latent_dim)
    
    message("Now predicting on training set")
    train_decoded <- decode_samples(new_samples = as.matrix(train_dat[, useMarkers]), 
                                    vae = vae_model$vae, batch_size = batch_size, 
                                    latent_dim = latent_dim)
    train_dat_norms <- train_dat_norms + train_decoded[[1]]
    train_dat_encoded <- train_dat_encoded + train_decoded$encoded %>%
      as.matrix()
  
    
    message("Now predicting on test set")
    
    test_decoded <- decode_samples(new_samples = as.matrix(test_dat[, useMarkers]), 
                                    vae = vae_model$vae, batch_size = batch_size, 
                                    latent_dim = latent_dim)
    test_dat_norms <- test_dat_norms + test_decoded[[1]]
    
    test_dat_encoded <- test_dat_encoded + test_decoded$encoded %>%
      as.matrix()
    
    message("All model fitting has completed")
  }
  return(list(train = train_dat_norms/length(kRange), test = test_dat_norms/length(kRange),
              train_encoded = train_dat_encoded/length(kRange), 
              test_encoded = test_dat_encoded/length(kRange)))
  
}


# Train a cell type classifier model using Caret
cellTypeClassifier <- function(train_x, test_x, model = "svmLinear"){
  library(doParallel)
  
  # cores <- 4
  # cl <- makePSOCKcluster(cores)
  # registerDoParallel(cl)
  # set up the grid
  fit.control <- trainControl(method = "cv", number = 3)
  
  # fit the model
  message('Fitting cell type model')
  set.seed(1994)  
  fit <- caret::train(cellTypes ~ ., data = train_x, method = model, 
               trControl = fit.control, trace = TRUE, preprocess = c('range'))
  print(fit)
  
  # predict on tes data
  message("Predicting on test")
  testCellTypes <- predict(fit, test_x)
  
  return(testCellTypes)
}

computeFeatures <- function(sce, featureType = "prop", cellTypeCol = "clusters",
                            sampleCol = "sample_id", logit = T, useMarkers, assay = "norm"){
  if(featureType == "prop"){
    features <- getProp(sce,
                         feature = cellTypeCol, imageID = sampleCol, logit = logit)
    colnames(features) <- paste0(colnames(features), "_logit")
    
  } else if(featureType == 'mean'){ # Marker mean per cell type
    # clusterCol <- rlang::sym(cellTypeCol)
    # sample_id <- rlang::sym(sampleCol)
    # get the colData
    coldat <- colData(sce) %>%
      as.data.frame()
    
    markerDf = 
      coldat |> 
      dplyr::select(sample_id, !!dplyr::sym(cellTypeCol)) |> 
      cbind(data.frame(t(assay(sce, assay)),check.names = F)) |>
      dplyr::select(sample_id, !!dplyr::sym(cellTypeCol),useMarkers)
    
    # Mean marker per cell type
    features = markerDf |> 
      group_by(sample_id, !!dplyr::sym(cellTypeCol)) |> 
      summarise_at(vars(-group_cols()), mean, na.rm = TRUE) |> 
      pivot_longer(-c(sample_id, !!dplyr::sym(cellTypeCol)), names_to = "markers") |> 
      pivot_wider(names_from = c(!!dplyr::sym(cellTypeCol), markers), values_from = value) |> 
      column_to_rownames("sample_id") %>% 
      replace(is.na(.), 0) 
  } 
  
  
  return(features)
}

# generate hierarchal tree using correlation
generateTree <- function(features, method = "ward"){
  # compute the correlation distance
  dists <- coop::pcor(features) %>%
    cor2dist() %>%
    as.dist()
  
  # generate the tree
  hcTree <- hclust(dists, method = method)
  
  # fix monotonicity
  # hcTree <- mhca::fixNonMonotHca(hcTree)
  hclustObj <- hcTree
  
  # Extract clusters from the dendrogram
  hcTree <- treekoR:::findChildren(ggtree(as.phylo(hcTree),
                                            ladderize = F, layout = 'dendrogram'))
  return(list(tree = hcTree, order = hclustObj$order))
}


# Generate hierarchal groups based on tree generated
generateGroups <- function(tree, nclust = 20, proportions, means, type = "all"){
  
  # Extract clusters and assign group numbers
  subGroups <- tree$data$clusters
  # print(subGroups)
  data <- tree$data
  data$group <- 1:length(subGroups)
  
  # Prepare variable groups for overlapping group Lasso using proportions
  var_groups<- list()
  n <- length(subGroups)
  # print(janitor::make_clean_names(colnames(proportions)))
  for(i in 1:length(subGroups)){
    var_groups[[i]] <- which(janitor::make_clean_names(colnames(proportions)) %in% 
                               janitor::make_clean_names(subGroups[[i]]))
   
  }
  
  
  if(type == "prop"){
    return(var_groups)
  }
  
  # Prepare variable groups for overlapping group Lasso using means
  last_max <- nclust
  off_set <- length(var_groups)
  for (i in 1:ncol(means)) {
    index <- off_set + i  # Calculate the new index in the list
    begin <- last_max + 1        # Start the next sequence right after the last max
    end <- begin      # Each sublist has exactly number of markers items
    
    # Add the new sequence to the resulting list at the new index
    var_groups[[index]] <- seq(begin, end)
    
    # Update last_max for the next iteration
    last_max <- end
  }
  
  # return the groupings
  return(var_groups)
}

visualiseModelTree <- function(fit, tree, type = "cluster", heatmap=TRUE,
                               training_data, nodes_to_remove = NULL,
                               title = "Feature tree"){
  
  coefs <- coef(fit) %>%
    as.matrix() %>%
    as.data.frame()
  coefs$feature <- rownames(coefs)
  
  # res_mp <- MP_gLasso(fit = fit, group = fit$grp.vec, lambda.type = "min", sort.type = "max")
  
  # head(res_mp$data) %>% print()
 
  p <- plotOverlapTree(model_fit = fit, training_data = training_data, 
                       tree_data = tree, nodes_to_remove = nodes_to_remove)
  if(heatmap == FALSE){
    p <- p + labs(title = "Proportions Tree") +
      # vexpand(.001, -1) + 
      theme(legend.title = element_text(color = "black", size = 16),
            plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
            legend.text = element_text(color = "black", size = 16))
    return(p)
  }
  
  # Extract the non-zero coefficients and non logit coefficients
  coefs$feature <- rownames(coefs)
  coefs  <- coefs %>%
    dplyr::filter(!feature %like% "logit")%>%
    dplyr::filter(!feature %like% "(Intercept)") %>%
    dplyr::filter(!feature %like% "Age") %>%
    dplyr::filter(!feature %like% "GenderM")
  
  # Splitting the full_string column
  splits <- strsplit(coefs$feature, "_", fixed = TRUE)
  
  # Creating new columns for cluster and marker
  if(type == "cluster"){
    coefs$Cell_Type <- sapply(splits, function(x) paste(x[1], x[2], sep = "_"))
    coefs$Marker <- sapply(splits, function(x) x[3])
  } else{
    coefs$Cell_Type <- gsub("_.*", "", coefs$feature)  # Remove everything after the underscore
    coefs$Marker <- sub(".*_", "", coefs$feature)  # Remove everything before and including the underscore
  }
  
  
  
  
  # print(coefs %>% head())
  
  # Pivot the data to wide format
  coefs_wide <- coefs %>%
    dplyr::select(Cell_Type, Marker, V1) %>%
    pivot_wider(names_from = Cell_Type, values_from = V1) %>%
    as.data.frame()
  
  rownames(coefs_wide) <- coefs_wide$Marker
  
  coefs_wide <- coefs_wide %>%
    dplyr::select(-Marker)
  
  # colnames(coefs_wide) <- janitor::make_clean_names(colnames(coefs_wide))
  colnames(coefs_wide) <- str_replace(colnames(coefs_wide), "logit", "")
  colnames(coefs_wide) <- str_replace(colnames(coefs_wide), "_", " ")
  
  
  val_range <- range(coefs_wide)
  # make sure it is centered at zero 
  # add a negative if the minimum is zero
  if(val_range[[1]] < 1e-16){
    val_range[[1]] <- -val_range[[2]]
  }
  
  print(val_range)
  colour_breaks <- c(val_range[1], 0, val_range[2])
  colours <- c("#4575B4", "white", "#D73027")
  
  
  coefs_wide <- coefs_wide[order(rownames(coefs_wide), decreasing = FALSE), ]
gheatmap(p, data = t(coefs_wide), offset = 0, hjust = 0.5, colnames_offset_y = -0.4, 
         colnames_offset_x = 0,
           color = 'black', font.size = 6, width = 2.5) +
    labs(title = title) +
    # vexpand(.001, -1) + 
    theme(legend.title = element_text(color = "black", size = 16),
          plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
          legend.text = element_text(color = "black", size = 16),
          legend.key.size = unit(2, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm')) + 
    scale_fill_gradientn(limits  = val_range, 
                         colours = colours[c(1, seq_along(colours), length(colours))],
                         values  = c(0, scales::rescale(colour_breaks, from = val_range), 1),
                         name = "Beta Means",
                         guide = guide_colorbar(order = 1)) 

  
}

getSigFeatures <- function(fit, type = "mean", mean = NULL, 
                prop = NULL, clinicaldata = NULL, outcome = NULL){
  
  # get the coefs from the model and filter 
  coefs_sub <- coef(fit) %>%
    as.matrix() %>%
    as.data.frame() %>%
    dplyr::mutate(feature = rownames(.)) %>%
    dplyr::filter(V1 != 0)
  
  if(type == "mean"){
    
    coefs_sub <- coefs_sub %>%
      dplyr::filter(!feature %like% '_logit') %>%
      dplyr::filter(!feature %like% 'Intercept') 
    
    sig_features <- mean %>%
      dplyr::select(coefs_sub$feature) %>%
      mutate(sample_id = as.integer(rownames(.))) %>%
      left_join(clinicaldata)
    
    stats <- sig_features %>% 
      group_by(!!dplyr::sym(outcome)) %>%
      dplyr::select(coefs_sub$feature) %>%
      as.data.frame() %>%
      melt() %>%
      as.data.frame() %>%
      mutate(feature = variable) %>%
      mutate(outcome = !!dplyr::sym(outcome)) %>%
      group_by(feature) %>%
      t_test(value ~ outcome) %>%
      adjust_pvalue(method = "fdr") %>%
      add_significance("p.adj") %>%
      add_xy_position(x = 1, dodge = 0.0) %>%
      # filter(p.adj < 0.05) %>%
      dplyr::mutate(p.adj = round(p.adj, 2))
  } else if(type == "prop"){
    coefs_sub <- coefs_sub %>%
      dplyr::filter(feature %like% '_logit') %>%
      dplyr::filter(!feature %like% 'Intercept') 
    
    sig_features <- prop %>%
      dplyr::select(coefs_sub$feature) %>%
      mutate(sample_id = as.integer(rownames(.))) %>%
      left_join(clinicaldata)
    
    stats <- sig_features %>% 
      group_by(!!dplyr::sym(outcome)) %>%
      dplyr::select(coefs_sub$feature) %>%
      as.data.frame() %>%
      melt() %>%
      as.data.frame() %>%
      mutate(feature = variable) %>%
      mutate(outcome = !!dplyr::sym(outcome)) %>%
      group_by(feature) %>%
      t_test(value ~ outcome) %>%
      adjust_pvalue(method = "fdr") %>%
      add_significance("p.adj") %>%
      add_xy_position(x = 1, dodge = 0.0) %>%
      # filter(p.adj < 0.05) %>%
      dplyr::mutate(p.adj = round(p.adj, 2))
  }
  
  colnames(sig_features) <- str_replace_all(colnames(sig_features), "_", " ")
  
  stats$feature <- str_replace_all(stats$feature, "_", " ")
  
  return(list(stats = stats, sig_features = sig_features ))

}

plotSigFeatures <- function(stats, sigFeatures, outcome, title = "", test = FALSE){
  
  
  p  <- sigFeatures %>% 
    dplyr::group_by(!!dplyr::sym(outcome)) %>%
    dplyr::select(stats$feature) %>%
    as.data.frame() %>%
    melt() %>%
    mutate(feature = variable) %>%
    ggboxplot(x = outcome, y = 'value',
              color = outcome, palette = 'jco',
              facet.by = 'feature', scale = 'free_y') +
    labs(title = title) +
    xlab("Gensini") +
    ylab("Mean Expression") +
    theme_minimal() + 
    theme_bw() +
    theme(axis.text.x = element_text(size = 16),
          plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
          axis.text.y = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.title.x = element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16),
          legend.position = 'None',
          strip.text.x = element_text(size = 16, colour = "black", face = 'bold')
    ) + stat_pvalue_manual(
      stats,  label = "p.adj", tip.length = 0,
      size = 6, color = 'red'
    ) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
  
  return(p)
}

plotAUC <- function(fit, x_test, y_test, title = ""){
  
  # get predictions from model
  preds <- predict(fit, X = x_test,
                          type = "response")
  preds_class <- predict(fit, X = x_test, type = "class")
  
  preds_df <- data.frame(Truth = y_test, Predicted = preds, Predicted_class = preds_class)
  
  # write.csv(data.frame(truth = condition, pred_test = pred_study_4),
  #           '../auc_data/fs_glasso_auc_study_4.csv')
  
  # get the TPR and FPR
  preds <- prediction(preds, y_test)
  preds_perf <- ROCR::performance(preds, "tpr", "fpr")
  
  # get the data for plotting
  plt_dat = data.frame(
    FPR = preds_perf@x.values[[1]],
    TPR = preds_perf@y.values[[1]]
  )
  
  # Compute the AUC
  auc_perf <- ROCR::performance(preds, measure = "auc")
  auc <- auc_perf@y.values[[1]]
  
  # generate ROC curve
  p.auc <- ggplot(plt_dat, aes(x = FPR, y = TPR)) +
    geom_line(colour = "blue") +
    labs(x = preds_perf@x.name, y = preds_perf@y.name) +
    geom_abline(slope = 1, intercept = 0) + theme_bw() +
    theme(
      plot.title = element_text(color="Black", size=16, face="bold", hjust = 0.5),
      plot.subtitle = element_text(color = "red", size = 16, hjust = 0.5),
      axis.title.x = element_text(color="Black", size=16),
      axis.title.y = element_text(color="Black", size=16),
      axis.text.y = element_text(size = 16),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 16),
      strip.text.x = element_text(size = 16),
      legend.title = element_text(size=16), #change legend title font size
      legend.text = element_text(size=16) #change legend text font size)
    )  + ggtitle(paste(title, round(auc, 2)))
  
  return(list(plot = p.auc, 
              preds = preds_df))
}

plotElbow <- function(aics){
  plotData <- data.frame(
    Alpha = seq(0.01, 0.1, 0.01),
    AIC = aics$aics[1:10]
  )
  # Sample dataframe
  df <- data.frame(value = c(10, 5, 20), row.names = c("C", "A", "B"))
  
  # View original dataframe
  print(df)
  
  # Order the dataframe by row names
  df_ordered <- df[order(rownames(df)), ]
  
  # View ordered dataframe
  print(df_ordered)
  ggpubr::ggline(
    plotData,
    x = "Alpha", y = 'AIC', group = 1, color = "steelblue"
  ) +
    geom_vline(xintercept = alpha*100, linetype = 1, color = "red") +
    labs(
      y = "AIC", x = "Alpha",
      title = "Optimal Alpha using the elbow method"
    ) +
    theme(
      plot.title = element_text(
        color = "Black", size = 14, face = "bold", hjust = 0.5
      ),
      axis.title.x = element_text(color = "Black", size = 14, face = "bold"),
      axis.title.y = element_text(color = "Black", size = 14, face = "bold")
    )
}

fitModel <- function(x_train, y_train, groups, alpha = NULL,
                     alpha_search = seq(0.01, 0.1, 0.01), 
                     penalty = 'cMCP', seed = '1994'){
  
  # if alpha is null, compute it
  message("Now computing alpha using the AIC method")
  if(is.null(alpha)){
    resAIC <- selectApha(x_train, y_train, alpha_search = alpha_search,
                         groups = groups, penalty = penalty, seed = seed)
    alpha <- alpha_search[computeElbow(resAIC$aics)]
    # alpha <- resAIC$best_alpha
  }
  
  
  
  message(paste("Optimal value for alpha is:", alpha))
  message("Now fitting final model")
  cvfit <- cv.grpregOverlap(x_train, y_train, groups, penalty = penalty,
                            family = 'binomial',  alpha = alpha, nfolds = 10,
                            seed = seed)
  # Plot cross-validation results
  par(mfrow=c(2,2))
  plot(cvfit, type = "all")
  
  # Fit overlapping group Lasso model
  fit <- grpregOverlap(x_train, y_train, groups,
                       family='binomial', alpha = alpha,
                       returnX.latent = T, returnOverlap = FALSE,
                       lambda = cvfit$lambda.min, penalty = penalty)
  message("Done")
  
  return(list(fit = fit))
}

plotHeatmap <- function(fit, type = "cluster", order = NULL,
                        markerOrder = NULL){
  
  coefs <- coef(fit) %>%
    as.matrix() %>%
    as.data.frame()
  coefs$feature <- rownames(coefs)
  
  # Extract the non-zero coefficients and non logit coefficients
  coefs$feature <- rownames(coefs)
  coefs  <- coefs %>%
    dplyr::filter(!feature %like% "logit")%>%
    dplyr::filter(!feature %like% "(Intercept)") %>%
    dplyr::filter(!feature %like% "Age") %>%
    dplyr::filter(!feature %like% "GenderM")
  
  # Splitting the full_string column
  splits <- strsplit(coefs$feature, "_", fixed = TRUE)
  
  # Creating new columns for cluster and marker
  if(type == "cluster"){
    coefs$Cell_Type <- sapply(splits, function(x) paste(x[1], x[2], sep = "_"))
    coefs$Marker <- sapply(splits, function(x) x[3])
  } else{
    coefs$Cell_Type <- gsub("_.*", "", coefs$feature)  # Remove everything after the underscore
    coefs$Marker <- sub(".*_", "", coefs$feature)  # Remove everything before and including the underscore
  }
  
  # print(coefs %>% head())
  
  # Pivot the data to wide format
  coefs_wide <- coefs %>%
    dplyr::select(Cell_Type, Marker, V1) %>%
    pivot_wider(names_from = Cell_Type, values_from = V1) %>%
    as.data.frame()
  
  rownames(coefs_wide) <- coefs_wide$Marker
  
  coefs_wide <- coefs_wide %>%
    dplyr::select(-Marker)
  
  # colnames(coefs_wide) <- janitor::make_clean_names(colnames(coefs_wide))
  colnames(coefs_wide) <- str_replace(colnames(coefs_wide), "logit", "")
  colnames(coefs_wide) <- str_replace(colnames(coefs_wide), "_", " ")
  
  
  val_range <- range(coefs_wide)
  print(val_range)
  colour_breaks <- c(val_range[1], 0, val_range[2])
  if(min(val_range[1]) < 0){
    colours <- c("#4575B4", "white", "#D73027")
  } else{
    colours <- c("white", "#4575B4", "#D73027")
  }
  
  if(!is.null(order)){
    # print(coefs_wide)
    coefs_wide <- coefs_wide[, order] %>%
      as.data.frame()
  }
  
  
  coefs_wide$Gene <- rownames(coefs_wide)
  coefs_wide_melted <- melt(coefs_wide)
  # print(coefs_wide_melted)
  colnames(coefs_wide_melted)[2:3] <- c("Cluster", "Value")
  coefs_wide_melted$Gene <- factor(coefs_wide_melted$Gene, 
                                   levels = unique(coefs_wide_melted$Gene[order(coefs_wide_melted$Gene, decreasing = TRUE)]))
  
  ggplot(coefs_wide_melted, aes(x = Cluster, y = Gene, fill = Value)) +
    geom_tile(color = "black",
              lwd = 0.05,
              linetype = 1) +
    ggtitle("Heatmap of celltype marker means") + 
    theme(legend.title = element_text(color = "black", size = 16),
          plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
          legend.text = element_text(color = "black", size = 16),
          axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                     hjust=1, size = 16, face = "bold"),
          axis.text.y = element_text(size = 16, face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) + 
    scale_fill_gradientn(limits  = val_range, 
                         colours = colours[c(1, seq_along(colours), length(colours))],
                         values  = c(0, scales::rescale(colour_breaks, from = val_range), 1),
                         name = "Beta Means") 
}


runPredictionModel <- function(sce_train, sce_test, numClusters, condition = "condition",
                               clinical_data_train, clinical_data_test){
  # cluster data into k clusters
  sce_train <- runFuseSOM(sce_train, numClusters = numClusters, assay = 'norm',
                          verbose = FALSE, clusterCol = 'clusters')
  
  # get training and testing data
  train_x <- reducedDim(sce_train, type = "VAE") %>%
    mutate(cellTypes = as.factor(sce_train$clusters))
  test_x <- reducedDim(sce_test, type = "VAE")
  
  # fit cell type classifier
  test_cluster <- cellTypeClassifier(train_x, test_x, model = 'lda')
  sce_test$clusters <- as.factor(test_cluster)
  
  # Compute training features
  data_logit_train <- computeFeatures(sce = sce_train, featureType = 'prop', 
                                      cellTypeCol = 'clusters', sampleCol = 'sample_id', 
                                      logit = T, useMarkers = useMarkers)
  row_names_train <- rownames(data_logit_train)
  condition_train <- factor(clinical_data_train[match(row_names_train,clinical_data_train$sample_id),
                                                condition])
  colnames(data_logit_train) <- as.factor(colnames(data_logit_train))
  markerMeanCellType_train <- computeFeatures(sce = sce_train, featureType = 'mean', 
                                              cellTypeCol = 'clusters', sampleCol = 'sample_id', 
                                              logit = T, useMarkers = useMarkers)
  colnames(markerMeanCellType_train) <- as.factor(colnames(markerMeanCellType_train))
  
  # compute testing features
  data_logit_test <- computeFeatures(sce = sce_test, featureType = 'prop', 
                                     cellTypeCol = 'clusters', sampleCol = 'sample_id', 
                                     logit = T, useMarkers = useMarkers)
  row_names_test <- rownames(data_logit_test)
  condition_test <- factor(clinical_data_test[match(row_names_test,clinical_data_test$sample_id),
                                              condition])
  colnames(data_logit_test) <- as.factor(colnames(data_logit_test))
  markerMeanCellType_test <- computeFeatures(sce = sce_test, featureType = 'mean', 
                                             cellTypeCol = 'clusters', sampleCol = 'sample_id', 
                                             logit = T, useMarkers = useMarkers)
  colnames(markerMeanCellType_test) <- as.factor(colnames(markerMeanCellType_test))
  
  
  # generate tree
  trees <- generateTree(features = data_logit_train, method = "ward")
  tree <- trees$tree
  order <- trees$order
  groups <- generateGroups(tree = tree, nclust = numClusters, 
                           proportions = data_logit_train, means = markerMeanCellType_train)
  
  # Combine the feature matrices
  X_train <- cbind(data_logit_train, markerMeanCellType_train) 
  
  # Scale the combined data
  scaleVals <- preProcess(X_train, method = c('scale'))
  X_train <- predict(scaleVals, X_train) %>%
    as.matrix()
  
  # Extract response variable
  y_train <- condition_train %>%
    as.numeric()
  # combine feature matrices for the test data
  X_test <- cbind(data_logit_test, markerMeanCellType_test) %>%
    as.matrix()
  X_test <- predict(scaleVals, X_test)
  
  # Fit the overlap group lasso model
  groups <- lapply(groups, sort)
  fit <- fitModel(x_train = X_train, y_train = y_train, groups = groups, penalty = "cMCP")
  
  test_auc <- plotAUC(fit = fit$fit, x_test = X_test, y_test = condition_test)
  return(test_auc)
}


runNormModel <- function(train, test, useMarkers, lambda = 0.1, numSamples = 2,
                         batch_size = 32, latent_dim = 16L){
  
  message("Now computing the reference samples")
  # get the reference samples
  res_ind <- ComputeReferenceSample(train, useMarkers, N = numSamples)
  
  # get the training and validation datasets
  train_dat <- train[which(train$sample_id %in% c(res_ind$topNSamples)), ]
  val_dat <- train[which(train$sample_id %in% c(res_ind$bottomNSamples)), ]
  
  message("Now fitting the MMD-VAE model")
  # fit the VAE model
  vae_model <- train_vae_model(train_dat = train_dat[, useMarkers], batch_size = batch_size, 
                               useMarkers = useMarkers, epochs = 100, lambda = lambda, l1_reg = 0,
                               val_dat = val_dat[, useMarkers], latent_dim = latent_dim)
  
  message("Now decoding the training and testing data")
  # Decode the training and testing datasets
  train_decoded <- decode_samples(new_samples = as.matrix(train[, useMarkers]), 
                                  vae = vae_model$vae)
  train_norm <- train_decoded[[1]]
  colnames(train_norm) <- useMarkers
  
  test_decoded <- decode_samples(new_samples = as.matrix(test[, useMarkers]), 
                                 vae = vae_model$vae)
  test_norm <- test_decoded[[1]]
  colnames(test_norm) <- useMarkers
  
  message("Now evaluating VAE model")
  train_vals <- vae_model$vae %>% evaluate(as.matrix(train[, useMarkers]), 
                                     as.matrix(train[, useMarkers]))
  
  test_vals <- vae_model$vae %>% evaluate(as.matrix(test[, useMarkers]), 
                                          as.matrix(test[, useMarkers]))
  
  message("Done")
  # create SCE objects for training and testing 
  sce_train <- SingleCellExperiment(assays = list(norm = t(train_norm[, useMarkers]),
                                                 raw = t(train[, useMarkers])
                                                 # norm_rescaled = t(df_norm_rescaled[, useMarkers])
  ),
  colData = train %>% dplyr::select(-useMarkers))
  reducedDims(sce_train) <- list(VAE=train_decoded$encoded %>%
                                  as.matrix() %>%
                                  as.data.frame())
  
  sce_test <- SingleCellExperiment(assays = list(norm = t(test_norm[, useMarkers]),
                                                   raw = t(test[, useMarkers])
                                                   # norm_rescaled = t(df_3_norm_rescaled[, useMarkers])
  ),
  colData = test %>% dplyr::select(-useMarkers))
  reducedDims(sce_test) <- list(VAE=test_decoded$encoded %>%
                                    as.matrix() %>%
                                    as.data.frame())
  
  return(list(sce_train = sce_train, sce_test = sce_test,
              train_vals = train_vals, test_vals = test_vals))
  
  
}