# ============================================================================
# Step 5 - Sketch Generation (exploratory subsample, NOT used for final stats)
# Run this AFTER pipeline.R has finished successfully.
# Run from Terminal with:  Rscript step5_sketch.R
# ============================================================================

cat("=== Step 5 started:", format(Sys.time()), "===\n")

rm(list = ls())
gc(full = TRUE)
set.seed(123)
options(future.globals.maxSize = 8000 * 1024^2)

suppressPackageStartupMessages({
  library(Seurat)
  library(BPCells)
  library(Matrix)
})

setwd("~/Downloads/LUAD_Project")

# ----------------------------------------------------------------------
# 1. Load the on-disk Seurat object produced by pipeline.R
#    (This is lightweight to load - the counts stay on disk via BPCells,
#    only the metadata/pointers are read into RAM.)
# ----------------------------------------------------------------------
cat("Loading on-disk Seurat object...\n")
seurat_full <- readRDS("GSE131907_full_on_disk.rds")
cat("  -> loaded:", ncol(seurat_full), "cells x", nrow(seurat_full), "genes\n")

# ----------------------------------------------------------------------
# 2. Normalize + find variable features
#    These operate on the on-disk matrix without pulling it fully into RAM,
#    but still checkpoint with gc() between them just in case.
# ----------------------------------------------------------------------
cat("Normalizing data...\n")
seurat_full <- NormalizeData(seurat_full)
gc(full = TRUE)

cat("Finding variable features...\n")
seurat_full <- FindVariableFeatures(seurat_full)
gc(full = TRUE)

# ----------------------------------------------------------------------
# 3. Sketch - done once, in a single call, no manual chunking needed here
#    because SketchData works off the on-disk matrix directly and only
#    materializes the ~50,000-cell subsample in RAM, not the full 208,506.
# ----------------------------------------------------------------------
cat("Generating sketch (50,000 cells via leverage score)...\n")
seurat_sketch <- SketchData(
  object         = seurat_full,
  ncells         = 50000,
  method         = "LeverageScore",
  sketched.assay = "sketch"
)
gc(full = TRUE)

saveRDS(seurat_sketch, "GSE131907_sketch.rds")

cat("=== STEP 5 COMPLETE ===\n")
print(seurat_sketch)
cat("=== Step 5 finished:", format(Sys.time()), "===\n")
