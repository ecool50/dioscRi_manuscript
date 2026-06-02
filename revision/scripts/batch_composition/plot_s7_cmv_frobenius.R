# Supp Fig S7: CMV (multi-study) Frobenius distance vs study structure.
# Reads pre-computed Frobenius distances from
#   revision/results/batch_composition/cmv_frobenius_distances.csv
#   revision/results/batch_composition/cmv_frobenius_avgnorms.csv
# To regenerate from the 2.6 GB raw CMV CSV (not in the Zenodo archive),
# run revision/scripts/batch_composition/precompute_cmv_frobenius.R first
# (obtain DeepLearning_data.csv from https://github.com/hzc363/DeepLearningCyTOF
# and set CMV_DIR to its parent directory).
#
# Output: figures/supp_fig_cmv_frobenius_batch.pdf

suppressPackageStartupMessages({
  library(here); library(dplyr); library(reshape2)
  library(ggplot2); library(patchwork)
})

revision_root <- here::here("revision")
in_dir   <- file.path(revision_root, "results", "batch_composition")
dist_csv <- file.path(in_dir, "cmv_frobenius_distances.csv")
avg_csv  <- file.path(in_dir, "cmv_frobenius_avgnorms.csv")

if (!file.exists(dist_csv) || !file.exists(avg_csv)) {
  stop("Pre-computed CMV Frobenius CSVs not found. Run\n  Rscript ",
       "revision/scripts/batch_composition/precompute_cmv_frobenius.R\n",
       "first (requires the raw 2.6 GB DeepLearning_data.csv).")
}

dist_df <- read.csv(dist_csv, stringsAsFactors = FALSE)
avg_df  <- read.csv(avg_csv,  stringsAsFactors = FALSE)

# Panel A inputs: within-study vs between-study pairwise distances
pair_df <- dist_df %>%
  filter(Sample1 != Sample2) %>%
  mutate(type = ifelse(Study1 == Study2, "Within-study", "Between-study")) %>%
  dplyr::select(dist = Distance, type)

# Panel C input: distance matrix ordered by study (for the heatmap)
sample_order <- avg_df$sample_id[order(avg_df$study, avg_df$sample_id)]
hm_df <- dist_df %>%
  mutate(Sample1 = factor(Sample1, levels = sample_order),
         Sample2 = factor(Sample2, levels = sample_order))

# Panel D input: per-sample avg distance ranked from most to least representative
avg_df_ranked <- avg_df %>%
  arrange(avg_norm) %>%
  mutate(rank = row_number())

study_cols <- c('SDY112' = '#E69F00', 'SDY113' = '#56B4E9', 'SDY305' = '#009E73',
                'SDY311' = '#CC79A7', 'SDY315' = '#D55E00', 'SDY472' = '#0072B2',
                'SDY478' = '#F0E442', 'SDY515' = '#999999', 'SDY519' = '#E64B35')

p_box <- ggplot(pair_df, aes(x = type, y = dist, fill = type)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  scale_fill_manual(values = c('Within-study' = '#56B4E9',
                               'Between-study' = '#E69F00')) +
  labs(x = NULL, y = 'Frobenius distance', fill = NULL,
       title = 'Covariance distance: within vs between study') +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')

p_study_box <- ggplot(avg_df, aes(x = study, y = avg_norm, fill = study)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.4) +
  scale_fill_manual(values = study_cols) +
  labs(x = 'Study', y = 'Average Frobenius distance',
       title = 'Representativeness distribution per study') +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1))

p_heatmap <- ggplot(hm_df, aes(x = Sample1, y = Sample2, fill = Distance)) +
  geom_tile() + scale_fill_viridis_c() +
  labs(title = 'Frobenius distance matrix (ordered by study)',
       fill = 'Distance') +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank())

p_rank <- ggplot(avg_df_ranked, aes(x = rank, y = avg_norm, colour = study)) +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_colour_manual(values = study_cols) +
  labs(x = 'Sample rank (most to least representative)',
       y = 'Average Frobenius distance', colour = 'Study',
       title = 'Sample representativeness by study') +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

p_frob <- (p_box | p_study_box) / (p_heatmap | p_rank) +
  plot_layout(heights = c(0.8, 1)) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = 'bold', size = 14))

out <- file.path(revision_root, "figures/supp_fig_cmv_frobenius_batch.pdf")
ggsave(out, p_frob, width = 12, height = 10)
cat("Saved:", out, "\n")
