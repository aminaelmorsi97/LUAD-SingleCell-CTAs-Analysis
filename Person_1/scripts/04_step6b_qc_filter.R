# ============================================================================
# Step 6b - Apply QC filtering (thresholds from the original paper)
# Run from Terminal with:  Rscript step6b_qc_filter.R
#
# Thresholds match the original dataset publication:
# Kim N, Kim HK, Lee K, et al. "Single-cell RNA sequencing demonstrates the
# molecular and cellular reprogramming of metastatic lung adenocarcinoma."
# Nat Commun. 2020;11:2285. https://doi.org/10.1038/s41467-020-16164-1
#
# Quoted QC criteria (Methods section):
#   - mitochondrial genes: <= 20%
#   - UMI count (nCount_RNA): 100 to 150,000
#   - gene count (nFeature_RNA): 200 to 10,000
# ============================================================================

cat("=== Step 6b started:", format(Sys.time()), "===\n")

rm(list = ls())
gc(full = TRUE)

suppressPackageStartupMessages({
  library(Seurat)
  library(BPCells)
  library(dplyr)
})

setwd("~/Downloads/LUAD_Project")

# ----------------------------------------------------------------------
# 1. Load the on-disk Seurat object (lightweight - counts stay on disk)
# ----------------------------------------------------------------------
cat("Loading on-disk Seurat object...\n")
seurat_full <- readRDS("GSE131907_full_on_disk.rds")
meta <- seurat_full@meta.data
cat("  -> starting cells:", nrow(meta), "\n")

# ----------------------------------------------------------------------
# 2. Define thresholds (documented here for the QC report / manuscript)
# ----------------------------------------------------------------------
MIN_COUNT   <- 100
MAX_COUNT   <- 150000
MIN_FEATURE <- 200
MAX_FEATURE <- 10000
MAX_MT      <- 20   # percent

# ----------------------------------------------------------------------
# 3. Identify which cells pass/fail each criterion (for the QC report -
#    this lets us report exactly how many cells were removed and why,
#    which reviewers will want to see)
# ----------------------------------------------------------------------
meta <- meta %>%
  mutate(
    fail_count   = nCount_RNA < MIN_COUNT | nCount_RNA > MAX_COUNT,
    fail_feature = nFeature_RNA < MIN_FEATURE | nFeature_RNA > MAX_FEATURE,
    fail_mt      = percent_mt > MAX_MT,
    fail_any     = fail_count | fail_feature | fail_mt
  )

cat("\n--- Cells failing each individual criterion ---\n")
cat("  nCount_RNA out of [", MIN_COUNT, ",", MAX_COUNT, "]:", sum(meta$fail_count), "\n")
cat("  nFeature_RNA out of [", MIN_FEATURE, ",", MAX_FEATURE, "]:", sum(meta$fail_feature), "\n")
cat("  percent_mt >", MAX_MT, "%:", sum(meta$fail_mt), "\n")
cat("  TOTAL removed (any criterion failed):", sum(meta$fail_any), "\n")

cells_to_keep <- rownames(meta)[!meta$fail_any]
cat("\n  -> cells kept:", length(cells_to_keep),
    "(", round(100 * length(cells_to_keep) / nrow(meta), 1), "% of original )\n")

# ----------------------------------------------------------------------
# 4. Save a QC filtering report (deliverable for the team / manuscript)
# ----------------------------------------------------------------------
qc_report <- data.frame(
  criterion = c(
    paste0("nCount_RNA outside [", MIN_COUNT, ", ", MAX_COUNT, "]"),
    paste0("nFeature_RNA outside [", MIN_FEATURE, ", ", MAX_FEATURE, "]"),
    paste0("percent_mt > ", MAX_MT, "%"),
    "TOTAL removed (any criterion)",
    "Cells retained"
  ),
  n_cells = c(
    sum(meta$fail_count),
    sum(meta$fail_feature),
    sum(meta$fail_mt),
    sum(meta$fail_any),
    length(cells_to_keep)
  ),
  percent_of_original = round(100 * c(
    sum(meta$fail_count),
    sum(meta$fail_feature),
    sum(meta$fail_mt),
    sum(meta$fail_any),
    length(cells_to_keep)
  ) / nrow(meta), 2)
)
write.csv(qc_report, "GSE131907_QC_filtering_report.csv", row.names = FALSE)

# ----------------------------------------------------------------------
# 5. Subset the on-disk object to just the passing cells.
#    This slices the BPCells matrix on disk - does NOT load the full
#    208,506-cell matrix into RAM.
# ----------------------------------------------------------------------
cat("\nSubsetting to QC-passing cells...\n")
seurat_filtered <- subset(seurat_full, cells = cells_to_keep)
gc(full = TRUE)
cat("  -> filtered object:", ncol(seurat_filtered), "cells x",
    nrow(seurat_filtered), "genes\n")

# ----------------------------------------------------------------------
# 6. Write a fresh on-disk BPCells directory for the filtered data,
#    and save the Seurat object pointing at it.
# ----------------------------------------------------------------------
cat("\nWriting filtered on-disk BPCells directory...\n")
if (dir.exists("GSE131907_bpcells_filtered")) {
  unlink("GSE131907_bpcells_filtered", recursive = TRUE)
}

write_matrix_dir(
  mat = seurat_filtered[["RNA"]]$counts,
  dir = "GSE131907_bpcells_filtered",
  compress = TRUE
)
gc(full = TRUE)

saveRDS(seurat_filtered, "GSE131907_filtered_on_disk.rds")

cat("\n=== STEP 6b COMPLETE ===\n")
print(seurat_filtered)
cat("Saved:\n")
cat("  - GSE131907_QC_filtering_report.csv   (how many cells removed & why)\n")
cat("  - GSE131907_bpcells_filtered/          (on-disk counts, QC-passed cells)\n")
cat("  - GSE131907_filtered_on_disk.rds       (Seurat object, use THIS from now on)\n")
cat("=== Step 6b finished:", format(Sys.time()), "===\n")
