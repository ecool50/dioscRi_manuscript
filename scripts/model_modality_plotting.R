stats <- stats_mean$stats
means <- stats_mean$sig_features
means$`Gensini bin` <- as.factor(means$`Gensini bin`)
outcome <- "Gensini bin" 

p1 <- means %>% 
  dplyr::group_by(!!sym(outcome)) %>%
  dplyr::select(stats$feature) %>%
  as.data.frame() %>%
  melt() %>%
  mutate(feature = variable) %>%
  ggplot(aes(value, colour = `Gensini bin`)) +
  geom_density() +
  facet_wrap(~feature, scales = 'free_y') +
  labs(title = "") +
  xlab("Means") +
  ylab("Density") +
  theme_minimal() + 
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        strip.text.x = element_text(size = 16, colour = "black", face = 'bold')
  ) + 
  guides(color=guide_legend(title="Gensini"))

p2 <- means %>% 
  dplyr::group_by(!!sym(outcome)) %>%
  dplyr::select(stats$feature) %>%
  as.data.frame() %>%
  melt() %>%
  mutate(feature = variable) %>%
  ggplot(aes(value)) +
  geom_density() +
  facet_wrap(~feature, scales = 'free') +
  labs(title = "Modality plots means") +
  xlab("Means") +
  ylab("Density") +
  theme_minimal() + 
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        strip.text.x = element_text(size = 10, colour = "black", face = 'bold')
  )

stats_p <- stats_prop$stats
props <- stats_prop$sig_features
props$`Gensini bin` <- as.factor(props$`Gensini bin`)

p1.prop <- props %>% 
  dplyr::group_by(!!sym(outcome)) %>%
  dplyr::select(stats_p$feature) %>%
  as.data.frame() %>%
  melt() %>%
  mutate(feature = variable) %>%
  ggplot(aes(value, colour = `Gensini bin`)) +
  geom_density() +
  facet_wrap(~feature, scales = 'free') +
  labs(title = "Modality plots of logit proportions faceted by Gensini") +
  xlab("Logit Proportions") +
  ylab("Density") +
  theme_minimal() + 
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        strip.text.x = element_text(size = 10, colour = "black", face = 'bold')
  )

p2.prop <- props %>% 
  dplyr::group_by(!!sym(outcome)) %>%
  dplyr::select(stats_p$feature) %>%
  as.data.frame() %>%
  melt() %>%
  mutate(feature = variable) %>%
  ggplot(aes(value)) +
  geom_density() +
  facet_wrap(~feature, scales = 'free') +
  labs(title = "Modality plots of logit proportions") +
  xlab("Logit Proportions") +
  ylab("Density") +
  theme_minimal() + 
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        strip.text.x = element_text(size = 10, colour = "black", face = 'bold')
  )

pdf(file = "Documents/PhD/bioheart_analysis/Plots/MG_non_zero_coefs_modality_plots_full_model.pdf",   # The directory you want to save the file in
    width = 16, # The width of the plot in inches
    height = 10) # The height of the plot in inches

p1
# plotSigFeatures(stats = stats_prop$stats, sigFeatures = stats_prop$sig_features, outcome = "Gensini bin")
# p2
# cowplot::plot_grid(p1.prop,p2.prop)

dev.off()