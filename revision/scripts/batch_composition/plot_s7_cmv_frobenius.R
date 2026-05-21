# Supp Fig S7: CMV (multi-study) Frobenius distance vs study structure.
# Same four panels as Supp Fig S6 but grouped by source study (SDY112..SDY515)
# instead of batch. PERMANOVA / permutation tests are computed elsewhere and
# not re-derived here.
#
# Output: figures/supp_fig_cmv_frobenius_batch.pdf

suppressPackageStartupMessages({
  library(here); library(data.table); library(dplyr); library(caret)
  library(dioscRi); library(reshape2); library(ggplot2); library(patchwork)
})

revision_root <- here::here("revision")
cmv_dir <- Sys.getenv("CMV_DIR",
  unset = "~/Documents/Academic/PhD/DeepLearning_CyTOF/DeepLearningCyTOF")

useMarkers <- c("TCRGD", "IGD", "HLADR", "CD94", "CD85J", "CD8", "CD56",
                "CD45RA", "CD4", "CD38", "CD33", "CD3", "CD28", "CD27",
                "CD25", "CD24", "CD20", "CD19", "CD161", "CD16", "CD14",
                "CD127", "CCR7")
split_and_extract <- function(string) strsplit(string, "\\.")[[1]][[2]]

cat("Reading raw CMV CSV (very large, ~2.8GB)...\n")
df <- fread(file.path(cmv_dir, 'DeepLearning_data.csv'), nThread = 7) %>%
  as.data.frame()
metaData <- read.csv(file.path(cmv_dir, 'cmv_metadata.csv'))
metaData <- metaData[, -1]
metaData$sample_id <- sapply(metaData$name, split_and_extract)
df$sample_id <- sapply(df$name, split_and_extract)
df <- df %>%
  left_join(metaData %>% dplyr::select(sample_id, study_accession, CMV_Ab) %>% distinct(),
            by = 'sample_id') %>%
  filter(study_accession != "SDY519")

df[, useMarkers] <- cyCombine::transform_asinh(df[, useMarkers],
                                               markers = useMarkers, derand = FALSE)
preProcValues <- preProcess(df[, useMarkers], method = c("range"))
df[, useMarkers] <- predict(preProcValues, df[, useMarkers])

samples <- unique(df$sample_id)
res <- computeReferenceSample(df, useMarkers, sampleCol = 'sample_id',
                              N = length(samples))
norm_mat <- res$Norms
rownames(norm_mat) <- samples; colnames(norm_mat) <- samples
sample_study <- df %>% dplyr::select(sample_id, study_accession) %>% distinct()

pair_df <- data.frame(dist = numeric(), type = character(), stringsAsFactors = FALSE)
for (i in 1:(length(samples) - 1)) {
  for (j in (i + 1):length(samples)) {
    d <- norm_mat[i, j]
    s_i <- sample_study$study_accession[sample_study$sample_id == samples[i]]
    s_j <- sample_study$study_accession[sample_study$sample_id == samples[j]]
    if (s_i == s_j) pair_df <- rbind(pair_df, data.frame(dist = d, type = 'Within-study'))
    else            pair_df <- rbind(pair_df, data.frame(dist = d, type = 'Between-study'))
  }
}

avg_df <- data.frame(
  sample_id = samples, avg_norm = res$avgNorms,
  study = sample_study$study_accession[match(samples, sample_study$sample_id)])

study_cols <- c('SDY112' = '#E69F00', 'SDY113' = '#56B4E9', 'SDY305' = '#009E73',
                'SDY311' = '#CC79A7', 'SDY315' = '#D55E00', 'SDY472' = '#0072B2',
                'SDY478' = '#F0E442', 'SDY515' = '#999999', 'SDY519' = '#E64B35')

p_box <- ggplot(pair_df, aes(x = type, y = dist, fill = type)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  scale_fill_manual(values = c('Within-study' = '#56B4E9', 'Between-study' = '#E69F00')) +
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

study_order <- order(sample_study$study_accession[match(samples, sample_study$sample_id)])
norm_ordered <- norm_mat[study_order, study_order]
hm_df <- melt(norm_ordered)
colnames(hm_df) <- c('Sample1', 'Sample2', 'Distance')
hm_df$Sample1 <- factor(hm_df$Sample1, levels = samples[study_order])
hm_df$Sample2 <- factor(hm_df$Sample2, levels = samples[study_order])

p_heatmap <- ggplot(hm_df, aes(x = Sample1, y = Sample2, fill = Distance)) +
  geom_tile() + scale_fill_viridis_c() +
  labs(title = 'Frobenius distance matrix (ordered by study)', fill = 'Distance') +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank())

avg_df <- avg_df %>% arrange(avg_norm)
avg_df$rank <- 1:nrow(avg_df)

p_rank <- ggplot(avg_df, aes(x = rank, y = avg_norm, colour = study)) +
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
