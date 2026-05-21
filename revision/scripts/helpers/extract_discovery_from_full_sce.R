# Extract and summarize discovery/validation cohorts from full SCE
# Run on server: Rscript extract_discovery_from_full_sce.R
#
# Input:  sce_b4_b3_elijah.RData (contains sce_3 and sce_4)
# Output: sce_b4_full.rds (discovery, all cells)
#         sce_b3_full.rds (validation, all cells)

library(SingleCellExperiment)

cat("Loading full SCE...\n")
load("sce_b4_b3_elijah.RData")

# --- Discovery (sce_4) ---
cat("\n=== Discovery (sce_4) ===\n")
cat("Cells:", ncol(sce_4), "\n")
cat("Markers:", nrow(sce_4), "\n")
cat("Samples:", length(unique(sce_4$sample_id)), "\n")
cat("Assays:", assayNames(sce_4), "\n")
cat("colData:", names(colData(sce_4)), "\n")

tab_4 <- sort(table(sce_4$sample_id))
cat("Cells per sample - min:", min(tab_4), "max:", max(tab_4), "median:", median(tab_4), "\n")

cat("\nSaving sce_b4_full.rds...\n")
saveRDS(sce_4, "sce_b4_full.rds")

# --- Validation (sce_3) ---
cat("\n=== Validation (sce_3) ===\n")
cat("Cells:", ncol(sce_3), "\n")
cat("Markers:", nrow(sce_3), "\n")
cat("Samples:", length(unique(sce_3$sample_id)), "\n")
cat("Assays:", assayNames(sce_3), "\n")

tab_3 <- sort(table(sce_3$sample_id))
cat("Cells per sample - min:", min(tab_3), "max:", max(tab_3), "median:", median(tab_3), "\n")

cat("\nSaving sce_b3_full.rds...\n")
saveRDS(sce_3, "sce_b3_full.rds")

# --- Summary ---
cat("\n=== Summary ===\n")
write.csv(as.data.frame(tab_4), "sce_b4_cell_counts.csv", row.names = FALSE)
write.csv(as.data.frame(tab_3), "sce_b3_cell_counts.csv", row.names = FALSE)

cat("Done.\n")
