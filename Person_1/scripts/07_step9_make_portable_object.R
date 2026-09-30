# ============================================================================
# Step 9 - Create a fully PORTABLE malignant-cell object (no BPCells needed)
# Run from Terminal:  Rscript step9_make_portable_object.R
#
# WHY THIS EXISTS:
#   The malignant-only object (GSE131907_malignant_full.rds) stores its counts
#   via BPCells - an on-disk pointer, not the actual numbers in the file. This
#   works fine on this machine, but teammates can't use it directly for
#   standard tools like pseudobulk aggregation or DESeq2/edgeR DE analysis,
#   which expect a normal in-memory matrix.
#
#   Since the malignant subset is small (~31,000 cells, not 208,506), we can
#   safely materialize the full counts matrix into memory and save a single,
#   completely standalone .rds file that works anywhere with zero setup:
#   just readRDS() and go.
# ============================================================================

cat("=== Step 9 started:", format(Sys.time()), "===\n")

rm(list = ls())
gc(full = TRUE)

suppressPackageStartupMessages({
  library(Seurat)
  library(BPCells)
  library(Matrix)
})

setwd("~/Downloads/LUAD_Project")

# ----------------------------------------------------------------------
# 1. Load the malignant-only on-disk object
# ----------------------------------------------------------------------
cat("Loading malignant on-disk Seurat object...\n")
seurat_malignant <- readRDS("GSE131907_malignant_full.rds")
cat("  -> loaded:", ncol(seurat_malignant), "cells x", nrow(seurat_malignant), "genes\n")

# ----------------------------------------------------------------------
# 2. Materialize the counts into a standard in-memory sparse matrix
#    (safe here because 31,136 cells is a manageable size, unlike the
#    full 208,506-cell dataset)
# ----------------------------------------------------------------------
cat("Materializing counts matrix into memory (this converts from BPCells\n")
cat("on-disk format to a standard dgCMatrix)...\n")

counts_bp <- seurat_malignant[["RNA"]]$counts
counts_in_memory <- as(counts_bp, "dgCMatrix")
gc(full = TRUE)
cat("  -> materialized matrix:", paste(dim(counts_in_memory), collapse = " x "), "\n")

# ----------------------------------------------------------------------
# 3. Build a brand-new, fully standalone Seurat object (no BPCells link)
# ----------------------------------------------------------------------
cat("Building portable Seurat object...\n")
meta <- seurat_malignant@meta.data

portable_object <- CreateSeuratObject(
  counts   = counts_in_memory,
  meta.data = meta,
  project  = "LUAD_malignant_portable"
)

# ----------------------------------------------------------------------
# 4. Save as one self-contained file
# ----------------------------------------------------------------------
saveRDS(portable_object, "GSE131907_malignant_portable.rds")

cat("\n=== STEP 9 COMPLETE ===\n")
print(portable_object)
cat("Saved: GSE131907_malignant_portable.rds\n")
cat("This single file is fully standalone - no BPCells folder needed.\n")
cat("Anyone on the team can just run: readRDS(\"GSE131907_malignant_portable.rds\")\n")
cat("=== Step 9 finished:", format(Sys.time()), "===\n")
