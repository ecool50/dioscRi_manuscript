sort_sublist <- function(x) {
  return(sort(x))
}

temp <- predict(fit, type = "vars")
coefs <- coef(fit) %>%
  as.matrix() %>%
  as.data.frame() %>%
  dplyr::mutate(feature = rownames(.)) %>%
  dplyr::filter(V1 != 0) %>%
  dplyr::filter(feature %like% '_logit') %>%
  dplyr::filter(!feature %like% 'Intercept') 

colnames(coefs)[[1]] <- 'coefficient'

grps <- list()
clusters_new <- list()
betas <- list()
for(i in 1:length(temp)){
  cur_feats <- temp[[i]]
  
  grp <- paste0('grp_',i)
  cluster <- colnames(X_train)[cur_feats]
  
  if(length(cur_feats) > 1){
    beta <- coefs[colnames(X_train)[cur_feats], 1] %>% abs() %>% mean
  } else{
    beta <- coefs[colnames(X_train)[cur_feats], 1]
  }
  
  
  grps[[i]] <- grp
  clusters_new[[i]] <- cluster
  betas[[i]] <- beta
}

res <- cbind(grps, clusters_new, betas) %>%
  as.data.frame() %>%
  dplyr::filter(clusters_new %like% "logit")
res$betas <- unlist(res$betas)


hc_dat <- tree$data
hc_dat$clusters_new <- lapply(hc_dat$clusters, function(x) gsub("_", " ", x))
hc_dat$clusters_new <- lapply(hc_dat$clusters_new, function(x) gsub(" logit", "", x))

# List of inner nodes to delete
nodes_to_delete <- c("Root", "Myeloid", "CD3", "CD4")

# Function to remove the specified nodes from each list
remove_nodes <- function(x, nodes) {
  return(setdiff(x, nodes))
}
# Apply the function to each sublist
hc_dat$clusters_new <- lapply(hc_dat$clusters_new, remove_nodes, nodes = nodes_to_delete)
hc_dat$clusters_new <- lapply(hc_dat$clusters_new, sort_sublist)


res$clusters_new <- lapply(res$clusters_new, function(x) gsub("_logit", "", x))
res$clusters_new <- lapply(res$clusters_new, function(x) gsub("_", " ", x))
res$clusters_new <- lapply(res$clusters_new, sort_sublist)
res <- res %>%
  dplyr::full_join(hc_dat, by = 'clusters_new') %>%
  dplyr::mutate(beta_internal = if_else(isTip == FALSE, betas, 0)) %>%
  dplyr::mutate(beta_leaf = if_else(isTip == TRUE, betas, 0)) %>%
  mutate_at(vars(betas, beta_leaf, beta_internal), ~replace_na(., 0))

# res$label <- str_replace(res$label, paste0('_','logit'), '')
# res$label <- str_replace(res$label, '_', " ")

# Create ggplot object for dendrogram visualization
p <- ggtree(res, layout = 'dendrogram', hang = 0) +
  geom_tiplab(as_ylab = T, geom = "text", size = 24, color = 'black') +
  geom_nodepoint(aes(subset = beta_internal != 0, color = beta_internal), size = 5) +
  geom_nodepoint(aes(color = beta_internal), size = 0) +
  scale_color_gradient2(mid = "white", high = muted("purple"),
                        midpoint=0, name = "Beta Internal") +
  new_scale_color() +
  geom_tippoint(aes(subset = beta_leaf != 0, color = beta_leaf), size = 5) + 
  scale_color_gradient2(low=muted("green"), high=muted("orange"), mid = "white",
                        midpoint=0, name = "Beta Leaf") +
  geom_tippoint(aes(color = beta_leaf), size = 0)
