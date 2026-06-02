# Precompute the CMV Frobenius distance matrix from the raw flat CSV
# (DeepLearning_data.csv, ~2.6 GB). The raw CSV is NOT in the Zenodo
# archive: obtain it from https://github.com/hzc363/DeepLearningCyTOF
# and set CMV_DIR to its parent directory.
#
# Outputs (small, committed to the repo so plot_s7_cmv_frobenius.R is
# reproducible from the Zenodo archive alone):
#   revision/results/batch_composition/cmv_frobenius_distances.csv
#   revision/results/batch_composition/cmv_frobenius_avgnorms.csv
#
# Re-run this script only if the raw CMV CSV changes.

suppressPackageStartupMessages({
  library(here); library(data.table); library(dplyr); library(caret)
  library(dioscRi); library(reshape2)
})

revision_root <- here::here("revision")
bioheart_root <- Sys.getenv("BIOHEART_ROOT",
  unset = "~/Documents/Academic/PhD/bioheart_analysis")
cmv_dir <- Sys.getenv("CMV_DIR",
  unset = "~/Documents/Academic/PhD/DeepLearning_CyTOF/DeepLearningCyTOF")
metadata_path <- file.path(bioheart_root, "data", "dioscRi_analysis_data",
                           "sce_objects", "CMV_Study_SDY519", "cmv_metadata.csv")

useMarkers <- c("TCRGD", "IGD", "HLADR", "CD94", "CD85J", "CD8", "CD56",
                "CD45RA", "CD4", "CD38", "CD33", "CD3", "CD28", "CD27",
                "CD25", "CD24", "CD20", "CD19", "CD161", "CD16", "CD14",
                "CD127", "CCR7")
split_and_extract <- function(string) strsplit(string, "\\.")[[1]][[2]]

cat("Reading raw CMV CSV (~2.6 GB)...\n")
df <- fread(file.path(cmv_dir, "DeepLearning_data.csv"), nThread = 7) %>%
  as.data.frame()
metaData <- read.csv(metadata_path)
metaData <- metaData[, -1]
metaData$sample_id <- sapply(metaData$name, split_and_extract)
df$sample_id <- sapply(df$name, split_and_extract)
df <- df %>%
  left_join(metaData %>%
              dplyr::select(sample_id, study_accession, CMV_Ab) %>%
              distinct(),
            by = "sample_id") %>%
  filter(study_accession != "SDY519")

df[, useMarkers] <- cyCombine::transform_asinh(df[, useMarkers],
                                                markers = useMarkers,
                                                derand = FALSE)
preProcValues <- preProcess(df[, useMarkers], method = c("range"))
df[, useMarkers] <- predict(preProcValues, df[, useMarkers])

cat("Computing pairwise Frobenius distances (full N x N)...\n")
samples <- unique(df$sample_id)
res <- computeReferenceSample(df, useMarkers, sampleCol = "sample_id",
                              N = length(samples))
norm_mat <- res$Norms
rownames(norm_mat) <- samples
colnames(norm_mat) <- samples
study_lookup <- df %>%
  dplyr::select(sample_id, study_accession) %>%
  distinct() %>%
  with(setNames(study_accession, sample_id))

out_dir <- file.path(revision_root, "results", "batch_composition")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dist_df <- melt(norm_mat) %>%
  setNames(c("Sample1", "Sample2", "Distance")) %>%
  mutate(Study1 = study_lookup[as.character(Sample1)],
         Study2 = study_lookup[as.character(Sample2)])
dist_path <- file.path(out_dir, "cmv_frobenius_distances.csv")
write.csv(dist_df, dist_path, row.names = FALSE)
cat("Wrote", dist_path, " (", nrow(dist_df), "pairs)\n")

avg_df <- data.frame(sample_id = samples,
                     avg_norm  = res$avgNorms,
                     study     = study_lookup[samples],
                     stringsAsFactors = FALSE)
avg_path <- file.path(out_dir, "cmv_frobenius_avgnorms.csv")
write.csv(avg_df, avg_path, row.names = FALSE)
cat("Wrote", avg_path, " (", nrow(avg_df), "samples)\n")
