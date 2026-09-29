# ============================================================================
# Step 7 - Extract malignant cells from the FULL, QC-FILTERED dataset
#           (not the 50k sketch, and not the unfiltered raw object)
# Run from Terminal with:  Rscript step7_malignant_full.R
#
# WHY THIS EXISTS: an earlier version isolated malignant cells from the
# 50,000-cell sketch object only. That is fine for quick exploration, but
# NOT defensible for the actual CTA prevalence / DE analysis, per the
# team's own rule: "no single random/sketch subsample as the basis for
# DE or prevalence estimates." This script instead extracts malignant
# cells from GSE131907_filtered_on_disk.rds - the full dataset AFTER
# QC filtering (step6b), using thresholds from the original paper
# (Kim et al. 2020, Nat Commun).
# ============================================================================

cat("=== Step 7 started:", format(Sys.time()), "===\n")

rm(list = ls())
gc(full = TRUE)

suppressPackageStartupMessages({
  library(Seurat)
  library(BPCells)
  library(dplyr)
})

setwd("~/Downloads/LUAD_Project")

# ----------------------------------------------------------------------
# 1. Load ONLY the metadata table, not the counts matrix.
#    seurat_full here is on-disk (BPCells-backed), so @meta.data is a
#    lightweight data.frame - this does not pull the matrix into RAM.
# ----------------------------------------------------------------------
cat("Loading on-disk Seurat object (metadata only is used here)...\n")
seurat_full <- readRDS("GSE131907_filtered_on_disk.rds")
meta <- seurat_full@meta.data

# ----------------------------------------------------------------------
# 2. Confirm the annotation column and value exist before filtering
#    (adjust the column name below if your original annotation column
#    is not literally "Cell_subtype")
# ----------------------------------------------------------------------
stopifnot("Cell_subtype" %in% colnames(meta))

cat("Cell_subtype value counts:\n")
print(table(meta$Cell_subtype, useNA = "ifany"))

# IMPORTANT: the original paper (Kim et al. 2020, Nat Commun) does NOT use
# "Malignant cells" as the label for primary tumor (tLung) cells. Instead,
# primary tumor malignant cells are labeled as tumor cell states tS1/tS2/tS3
# (see paper Results: "upregulated signatures in tumor cell states 1, 2, 3").
# "Malignant cells" is used for metastatic-site tumor cells (mLN, PE, etc.).
# We must include ALL FOUR categories to get the complete malignant population
# across primary AND metastatic sites.
malignant_labels <- c("Malignant cells", "tS1", "tS2", "tS3")
malignant_cells <- rownames(meta)[meta$Cell_subtype %in% malignant_labels]
cat("\n-> malignant cells found in FULL dataset (Malignant cells + tS1/tS2/tS3):",
    length(malignant_cells), "\n")
cat("   Breakdown by subtype:\n")
print(table(meta$Cell_subtype[meta$Cell_subtype %in% malignant_labels]))
cat("   (expected ~31,000 total across early primary / advanced primary / mLN / brain mets,\n")
cat("    matching the original paper's reported counts)\n")

# ----------------------------------------------------------------------
# 3. Slice the on-disk BPCells matrix to just these cells.
#    This does NOT load the full matrix into RAM - it operates on the
#    on-disk representation and only materializes what's requested.
# ----------------------------------------------------------------------
cat("\nSubsetting Seurat object to malignant cells only...\n")
seurat_malignant_full <- subset(seurat_full, cells = malignant_cells)
gc(full = TRUE)

cat("  -> malignant object:", ncol(seurat_malignant_full), "cells x",
    nrow(seurat_malignant_full), "genes\n")

# ----------------------------------------------------------------------
# 4. Save both:
#    (a) a fresh on-disk BPCells directory just for malignant cells
#        (fast to re-load later without re-slicing)
#    (b) the Seurat object pointing at it
# ----------------------------------------------------------------------
cat("\nWriting malignant-only on-disk BPCells directory...\n")
if (dir.exists("GSE131907_bpcells_malignant")) {
  unlink("GSE131907_bpcells_malignant", recursive = TRUE)
}

write_matrix_dir(
  mat = seurat_malignant_full[["RNA"]]$counts,
  dir = "GSE131907_bpcells_malignant",
  compress = TRUE
)
gc(full = TRUE)

saveRDS(seurat_malignant_full, "GSE131907_malignant_full.rds")

cat("\n=== STEP 7 COMPLETE ===\n")
print(seurat_malignant_full)
cat("Saved:\n")
cat("  - GSE131907_bpcells_malignant/   (on-disk counts, full malignant cells)\n")
cat("  - GSE131907_malignant_full.rds   (Seurat object, use this - NOT the sketch-based one)\n")
cat("=== Step 7 finished:", format(Sys.time()), "===\n")
