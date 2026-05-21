# Batch Composition of Reference Samples [R2-2, R1-1]
# Run computeReferenceSample on BioHEART-CT discovery cohort
# and report which batches the selected reference samples come from.

rm(list = ls())
library(dioscRi)
library(tidyverse)
library(data.table)
library(caret)

bioheart_root <- Sys.getenv("BIOHEART_ROOT", unset = "~/Documents/Academic/PhD/bioheart_analysis")
revision_root <- here::here("revision")
base_dir <- bioheart_root  # alias: upstream inputs (raw data, sce objects, etc.)
useMarkers <- c('HLA_DR', 'CD3', 'CD4', 'CD8a', 'CD25', 'CD127', 'FoxP3', 'CD27',
                'KLRG1', 'CD56', 'CD45RO', 'CD45RA', 'CD192_CCR2', 'CD194_CCR4',
                'CD196_CCR6', 'CD39', 'CD38', 'Ki67', 'CD183_CXCR3', 'CCR7',
                'CD19', 'CD20', 'IgD', 'CD14', 'CD304', 'CD141', 'CD1c_PE')

# Load raw CSV (prefer the cleaned copy if present, fall back to the raw csv)
cat("Loading BioHEART-CT discovery cohort...\n")
csv_path <- file.path(bioheart_root, 'data/raw_data/bioheart_b4_clean.csv')
if (!file.exists(csv_path)) {
  csv_path <- file.path(bioheart_root, 'data/raw_data/bioheart_ct_cytof_data_b4_mg.csv')
}
df <- fread(csv_path, nThread = 7) %>% as.data.frame()

# Preprocess (same as pipeline: arcsinh + min-max)
df[, useMarkers] <- cyCombine::transform_asinh(df[, useMarkers],
                                                markers = useMarkers, derand = FALSE)
preProcValues <- preProcess(df[, useMarkers], method = c("range"))
df[, useMarkers] <- predict(preProcValues, df[, useMarkers])

# Run reference sample selection with K=2
cat("Computing reference samples (K=2)...\n")
res <- computeReferenceSample(df, useMarkers, sampleCol = 'sample_id', N = 2)

# Map samples to batches
sample_batch <- df %>%
  dplyr::select(sample_id, CyTOF.Batch) %>%
  distinct()

cat("\n=== Results ===\n")

cat("\nTop 2 reference samples (used for VAE training):\n")
top <- sample_batch %>% filter(sample_id %in% res$topNSamples)
print(top)

cat("\nBottom 2 samples (used for VAE validation):\n")
bottom <- sample_batch %>% filter(sample_id %in% res$bottomNSamples)
print(bottom)

cat("\nAll batches in discovery cohort:\n")
print(table(sample_batch$CyTOF.Batch))

cat("\nConclusion: Reference samples come from batches:",
    paste(top$CyTOF.Batch, collapse = " and "), "\n")

# ============================================================
# CMV
# ============================================================
cat("\n=== CMV ===\n")

split_and_extract <- function(string, n) {
  strsplit(string, "\\.")[[1]][[2]]
}

useMarkers_cmv <- c("TCRGD","IGD","HLADR","CD94","CD85J","CD8","CD56",
                     "CD45RA","CD4","CD38","CD33","CD3","CD28","CD27",
                     "CD25","CD24","CD20","CD19","CD161","CD16","CD14","CD127","CCR7")

cmv_dir <- file.path(
  Sys.getenv("DEEPLEARNING_CYTOF_ROOT",
             unset = "~/Documents/Academic/PhD/DeepLearning_CyTOF"),
  "DeepLearningCyTOF")

cat("Loading CMV data...\n")
dat_full <- fread(file.path(cmv_dir, 'DeepLearning_data.csv'), nThread = 7) %>%
  as.data.frame()
metaData_cmv <- read.csv(file.path(cmv_dir, 'cmv_metadata.csv'))
metaData_cmv <- metaData_cmv[, -1]
metaData_cmv$sample_id <- sapply(metaData_cmv$name, split_and_extract)

dat_full <- left_join(metaData_cmv, dat_full)

# Training set: exclude SDY519 and SDY515 (same as original pipeline)
train_cmv <- dat_full[-which(dat_full$study_accession %in% c('SDY519', 'SDY515')), ]

# Min-max scaling only (no arcsinh, matching original pipeline)
preProcValues_cmv <- preProcess(train_cmv[, useMarkers_cmv], method = c("range"))
train_cmv[, useMarkers_cmv] <- predict(preProcValues_cmv, train_cmv[, useMarkers_cmv])

cat("Computing reference samples (K=2)...\n")
res_cmv <- computeReferenceSample(train_cmv, useMarkers_cmv, sampleCol = 'sample_id', N = 2)

# Map to studies
sample_study <- train_cmv %>%
  dplyr::select(sample_id, study_accession) %>%
  distinct()

cat("\nTop 2 reference samples (used for VAE training):\n")
top_cmv <- sample_study %>% filter(sample_id %in% res_cmv$topNSamples)
print(top_cmv)

cat("\nBottom 2 samples (used for VAE validation):\n")
bottom_cmv <- sample_study %>% filter(sample_id %in% res_cmv$bottomNSamples)
print(bottom_cmv)

cat("\nAll studies in training data:\n")
print(table(sample_study$study_accession))

cat("\nConclusion: Reference samples come from studies:",
    paste(top_cmv$study_accession, collapse = " and "), "\n")

# ============================================================
# Combined batch / study composition table (CSV only)
# ============================================================
# BioHEART batch counts, annotated with which batches contain the K=2 refs
batch_counts_bh <- as.data.frame(table(sample_batch$CyTOF.Batch))
colnames(batch_counts_bh) <- c('Batch', 'Total samples')
batch_counts_bh$`Training ref`   <- if_else(batch_counts_bh$Batch %in% top$CyTOF.Batch,    'Yes', '')
batch_counts_bh$`Validation ref` <- if_else(batch_counts_bh$Batch %in% bottom$CyTOF.Batch, 'Yes', '')

# CMV study counts, annotated similarly
study_counts_cmv <- as.data.frame(table(sample_study$study_accession))
colnames(study_counts_cmv) <- c('Study', 'Total samples')
study_counts_cmv$`Training ref`   <- if_else(study_counts_cmv$Study %in% top_cmv$study_accession,    'Yes', '')
study_counts_cmv$`Validation ref` <- if_else(study_counts_cmv$Study %in% bottom_cmv$study_accession, 'Yes', '')

combined <- bind_rows(
  batch_counts_bh   %>% rename(Group = Batch) %>% mutate(Dataset = 'BioHEART-CT'),
  study_counts_cmv  %>% rename(Group = Study) %>% mutate(Dataset = 'CMV')
)
write.csv(combined,
          file.path(revision_root, 'results/batch_composition/reference_sample_batch_composition.csv'),
          row.names = FALSE)

cat("\nTables saved to figures/ and data/\n")
