# Supp Fig S6: BioHEART-CT Frobenius distance vs batch structure.
# Reads the raw discovery cohort CSV, computes per-sample covariance Frobenius
# norms, and builds the four-panel figure (within/between batch, per-batch
# distribution, distance-matrix heatmap, ranked representativeness).
#
# Output: figures/supp_fig_bioheart_frobenius_batch.pdf

suppressPackageStartupMessages({
  library(here); library(data.table); library(dplyr); library(caret)
  library(dioscRi); library(reshape2); library(ggplot2); library(patchwork)
})

revision_root <- here::here("revision")
bioheart_root <- Sys.getenv("BIOHEART_ROOT",
                            unset = "~/Documents/Academic/PhD/bioheart_analysis")
raw_csv_candidates <- c(
  file.path(bioheart_root, "data/dioscRi_analysis_data/raw_files",
            "bioheart_ct_cytof_data_b4_mg.csv"),                 # Zenodo
  file.path(bioheart_root, "data/raw_data/bioheart_b4_clean.csv"),  # local cleaned
  file.path(bioheart_root, "data/raw_data/bioheart_ct_cytof_data_b4_mg.csv"))
raw_csv <- raw_csv_candidates[file.exists(raw_csv_candidates)][1]
if (is.na(raw_csv)) stop("Raw BioHEART CSV not found in any of:\n  ",
                         paste(raw_csv_candidates, collapse = "\n  "))

useMarkers <- c('HLA_DR', 'CD3', 'CD4', 'CD8a', 'CD25', 'CD127', 'FoxP3', 'CD27',
                'KLRG1', 'CD56', 'CD45RO', 'CD45RA', 'CD192_CCR2', 'CD194_CCR4',
                'CD196_CCR6', 'CD39', 'CD38', 'Ki67', 'CD183_CXCR3', 'CCR7',
                'CD19', 'CD20', 'IgD', 'CD14', 'CD304', 'CD141', 'CD1c_PE')

cat("Reading raw BioHEART CSV (large)...\n")
df <- fread(raw_csv, nThread = 7) %>% as.data.frame()
df[, useMarkers] <- cyCombine::transform_asinh(df[, useMarkers],
                                               markers = useMarkers, derand = FALSE)
preProcValues <- preProcess(df[, useMarkers], method = c("range"))
df[, useMarkers] <- predict(preProcValues, df[, useMarkers])
colnames(df)[colnames(df) %in% useMarkers] <-
  gsub("_", "-", colnames(df)[colnames(df) %in% useMarkers])
useMarkers_dash <- gsub("_", "-", useMarkers)

res <- computeReferenceSample(df, useMarkers_dash, sampleCol = 'sample_id',
                              N = length(unique(df$sample_id)))
norm_mat <- res$Norms
samples <- unique(df$sample_id)
rownames(norm_mat) <- samples; colnames(norm_mat) <- samples
sample_batch <- df %>% dplyr::select(sample_id, CyTOF.Batch) %>% distinct()

pair_df <- data.frame(dist = numeric(), type = character(), stringsAsFactors = FALSE)
for (i in 1:(length(samples) - 1)) {
  for (j in (i + 1):length(samples)) {
    d <- norm_mat[i, j]
    b_i <- sample_batch$CyTOF.Batch[sample_batch$sample_id == samples[i]]
    b_j <- sample_batch$CyTOF.Batch[sample_batch$sample_id == samples[j]]
    if (b_i == b_j) pair_df <- rbind(pair_df, data.frame(dist = d, type = 'Within-batch'))
    else            pair_df <- rbind(pair_df, data.frame(dist = d, type = 'Between-batch'))
  }
}

avg_df <- data.frame(
  sample_id = samples, avg_norm = res$avgNorms,
  batch = sample_batch$CyTOF.Batch[match(samples, sample_batch$sample_id)])

batch_cols <- c('4a' = '#E69F00', '4b' = '#56B4E9', '4c' = '#009E73',
                '4d' = '#CC79A7', '4e' = '#D55E00', '4f' = '#0072B2')

p_box <- ggplot(pair_df, aes(x = type, y = dist, fill = type)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  scale_fill_manual(values = c('Within-batch' = '#56B4E9', 'Between-batch' = '#E69F00')) +
  labs(x = NULL, y = 'Frobenius distance', fill = NULL,
       title = 'Covariance distance: within vs between batch') +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')

batch_order <- order(sample_batch$CyTOF.Batch)
norm_ordered <- norm_mat[batch_order, batch_order]
hm_df <- melt(norm_ordered)
colnames(hm_df) <- c('Sample1', 'Sample2', 'Distance')
hm_df$Sample1 <- factor(hm_df$Sample1, levels = samples[batch_order])
hm_df$Sample2 <- factor(hm_df$Sample2, levels = samples[batch_order])

p_heatmap <- ggplot(hm_df, aes(x = Sample1, y = Sample2, fill = Distance)) +
  geom_tile() + scale_fill_viridis_c() +
  labs(title = 'Frobenius distance matrix (ordered by batch)', fill = 'Distance') +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank())

avg_df <- avg_df %>% arrange(avg_norm)
avg_df$rank <- 1:nrow(avg_df)

p_rank <- ggplot(avg_df, aes(x = rank, y = avg_norm, colour = batch)) +
  geom_point(size = 2.5) + scale_colour_manual(values = batch_cols) +
  labs(x = 'Sample rank (most to least representative)',
       y = 'Average Frobenius distance', colour = 'Batch',
       title = 'Sample representativeness by batch') +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

p_batch_box <- ggplot(avg_df, aes(x = batch, y = avg_norm, fill = batch)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
  scale_fill_manual(values = batch_cols) +
  labs(x = 'Batch', y = 'Average Frobenius distance',
       title = 'Representativeness distribution per batch') +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')

p_frob <- (p_box | p_batch_box) / (p_heatmap | p_rank) +
  plot_layout(heights = c(0.8, 1)) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = 'bold', size = 14))

out <- file.path(revision_root, "figures/supp_fig_bioheart_frobenius_batch.pdf")
ggsave(out, p_frob, width = 12, height = 10)
cat("Saved:", out, "\n")
