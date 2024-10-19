colourTreeCustom <- function(tree,
                       point_size = 3,
                       high = "#00c434",
                       low = "purple",
                       mid="ivory2"){
  
  tree$data$label <- ifelse(is.na(tree$data$label),
                            tree$data$node,
                            tree$data$label)
  
  tooltip <- paste("<b>Node</b>:", tree$data$label,
                   ", beta = ", round(tree$data$beta_final,3))
  
  tree_df <- tree$data
  
  # Add colour aes mapping to the first tree layer for the branch colours
  # tree$layers[[1]]$mapping <- modifyList(tree$layers[[1]]$mapping, aes(color=beta_color))
  # tree$layers[[2]]$mapping <- modifyList(tree$layers[[2]]$mapping, aes(color=beta_color))
  
  tree <- tree +
    geom_point_interactive(data = tree$data, 
                           aes(x,y, tooltip = tooltip, data_id = label),
                           size = point_size) +
    theme(legend.position = "right")
  
  return(tree)
}





temp <- colourTreeCustom(tree = hc_prop)

gf <- girafe(
  print(temp + theme(legend.position = "top")),
  width_svg=13,
  height_svg=9)

hover_css <- "fill:cyan; stroke:darkcyan; r:4pt;"
tooltip_css <- "border-style: solid; border-color: #c3c3c3; border-radius: 8px;
                  border-width: 0.5px; padding: 5px;box-shadow: 2px 2px 3px 0px #bbbbbb;
                  background-color: white; font: menu;"

# girafe_options(gf, opts_hover(css = hover_css), opts_tooltip(css = tooltip_css))
  # scale_color_manual(values=c("blue", "red","grey", "purple")) +
  
# geom_nodelab(aes(label=round(beta, 2)), hjust=-.35)
# geom_tiplab(aes(label=round(beta_leaf, 2)), hjust=1.5, vjust = 0.3)


