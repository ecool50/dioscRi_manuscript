res_dat <- res_mp$data
res_dat$beta_name <- janitor::make_clean_names(res_dat$beta_name)
# res_dat$beta_name <- clean_strings(res_dat$beta_name)

# get the prop groups
res_dat_prop <- res_dat[res_dat$beta_name %like% 'prop', ]
res_dat_prop$group <- as.numeric(res_dat_prop$group)
coefs_prop <- coefs[coefs$feature %like% 'prop', ] %>%
  base::merge(data.frame(feature = colnames(data), 
                         beta = rep(0, ncol(data))), all = T) %>%
  replace(is.na(.), 0) %>%
  dplyr::select(-beta)

colnames(coefs_prop) <- c('feature', 'beta')

# get the logit groups
res_dat_logit <- res_dat[res_dat$beta_name %like% 'logit', ]
res_dat_logit$group <- as.numeric(res_dat_logit$group) - length(subs_prop)

## Highlight all the single groups for prop
h_prop_single <- res_dat_prop[table(res_dat_prop$groups) == 1, ]$beta_name
h_prop_single <- gsub("grp[0-9]+_", "", h_prop_single)
dend_prop <- hc_prop %>%
  as.dendrogram()
dend_prop <- color_labels(dend_prop, labels = h_prop_single,
                          col = rainbow(n = length(h_prop_single)))
plot(dend_prop)


## Highlight all the single groups for logit
h_logit_single <- res_dat_logit[table(res_dat_logit$groups) == 1, ]$beta_name
h_logit_single <- gsub("grp[0-9]+_", "", h_logit_single)
dend_logit <- hc_logit %>%
  as.dendrogram()
dend_logit <- color_labels(dend_logit, labels = h_logit_single,
                           col = rainbow(n = length(h_logit_single)))
plot(dend_logit)



res <- getHCGroups(hc)
var <- res$varGroup
group <- res$indexGroup
heights <- res$heights


df_info <- data.frame(height = heights, group = group, idx = var)

phylo <- as.phylo.hclust(hc)

plot(phylo, show.node.label = TRUE)
nodelabels()
tiplabels() # yellow
edgelabels() # green

subs <- partition_leaves(dend)

pathRoutes <- function(leaf) {
  which(sapply(subs, function(x) leaf %in% x))
}

paths <- lapply(subs[[1]], pathRoutes)


dend <- as.dendrogram(hc_prop)
# dend <- dend %>% set_labels(colnames(data)[dend %>% labels()])
labs <- dend %>% labels()
coefs_prop <- coefs_prop[match(labs, coefs_prop$feature),]

xy <- dend %>% 
  get_nodes_xy() %>%
  as.data.frame()
colnames(xy) <- c('x', 'y')

df_lab <- xy 
df_lab$group <- rep(1:nrow(df_lab))

df_lab <- base::merge(res_dat_prop, df_lab, all = T) %>%
  dplyr::select(group, beta, x, y) %>%
  replace(is.na(.), 0) %>%
  abs() %>%
  dplyr::group_by(group) %>%
  dplyr::summarise_at(vars(beta, x, y), list(mean))

is_leaf <- !is.na(dend %>% get_nodes_attr("leaf"))
df_lab_leaf <- xy[is_leaf, ]
coefs_prop <- cbind(coefs_prop, df_lab_leaf) %>%
  as.data.frame()
coefs_prop$beta <- round(coefs_prop$beta, 2)


is_internal_node <- is.na(dend %>% get_nodes_attr("leaf"))
is_internal_node[which.max(df_lab$y)] <- FALSE
df_lab <- df_lab[is_internal_node,]

plot(dend)
text(coefs_prop$x+.2, coefs_prop$y+.2, labels=format(coefs_prop$beta, digits=2), col="blue")
text(df_lab$x+.2, df_lab$y+.2, labels=format(df_lab$beta, digits=1), col="red")
