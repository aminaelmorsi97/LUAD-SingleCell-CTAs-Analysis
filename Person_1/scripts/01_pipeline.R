# ============================================================================
# GSE131907 - Data Engineering & QC Pipeline (Member 1)
# Run this from Terminal with:  Rscript pipeline.R
# (Do NOT run via RStudio's Run button - see notes below)
# ============================================================================

# ----------------------------------------------------------------------
# 0. Session log (helps debugging if it crashes again)
# ----------------------------------------------------------------------
cat("=== Pipeline started:", format(Sys.time()), "===\n")
cat("R version:", R.version.string, "\n")

# ----------------------------------------------------------------------
# 1. Clean Memory & Workspace First
# ----------------------------------------------------------------------
rm(list = ls())
gc(full = TRUE)
set.seed(123)
options(future.globals.maxSize = 8000 * 1024^2)

# ----------------------------------------------------------------------
# 2. Load Required Packages (install once manually, not inside the script)
# ----------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(BPCells)
  library(Matrix)
  library(dplyr)
})

# Log package versions -> helps catch Seurat/BPCells/Matrix mismatches
cat("Seurat:", as.character(packageVersion("Seurat")), "\n")
cat("BPCells:", as.character(packageVersion("BPCells")), "\n")
cat("Matrix:", as.character(packageVersion("Matrix")), "\n")

# ----------------------------------------------------------------------
# 3. Set Working Directory
# ----------------------------------------------------------------------
setwd("~/Downloads/LUAD_Project")

# ----------------------------------------------------------------------
# 4. Load Raw Metadata & Curate Table
# ----------------------------------------------------------------------
cat("Step 1: Loading and Curating Metadata...\n")
raw_metadata <- read.delim("GSE131907_cell_annotation.txt", header = TRUE, row.names = 1)
raw_metadata$cell_id <- rownames(raw_metadata)
write.csv(raw_metadata, "GSE131907_metadata_curated.csv", row.names = TRUE)
cat("  -> metadata rows:", nrow(raw_metadata), "\n")

# ----------------------------------------------------------------------
# 5. Convert Raw UMI Counts to On-Disk BPCells Matrix
#    THIS IS THE STEP THAT WAS CRASHING - most memory-hungry part.
# ----------------------------------------------------------------------
cat("Step 2: Converting counts matrix to BPCells on-disk storage...\n")

if (dir.exists("GSE131907_bpcells")) {
  unlink("GSE131907_bpcells", recursive = TRUE)
}

# --- OPTION A (preferred, if you have the original mtx/barcodes/features files
#     from GEO instead of a pre-built RDS) - this NEVER loads the full matrix
#     into RAM, it streams straight to disk. Uncomment if you have these files:
#
# counts_bpcells <- import_matrix_market(
#   mtx_path = "matrix.mtx",
#   outdir   = "GSE131907_bpcells"
# )

# --- OPTION B (what you currently have: a single .rds file, stored as a
#     data.frame - NOT a matrix. Converting the whole thing to a matrix in
#     one shot needs way too much RAM, so we do it column-chunk by
#     column-chunk instead, freeing memory after every chunk.) ---
cat("  Reading RDS file (this is the peak-memory moment)...\n")
raw_counts <- readRDS("GSE131907_raw_UMI_matrix.rds")
cat("  -> class:", class(raw_counts)[1], "\n")
cat("  -> dimensions:", paste(dim(raw_counts), collapse = " x "), "\n")

gene_names <- rownames(raw_counts)
cell_names <- colnames(raw_counts)
n_cols     <- ncol(raw_counts)

# Smaller chunk = safer on limited RAM, but more iterations (slower).
# 5000 cells per chunk is a reasonable starting point for 8-16GB machines.
chunk_size <- 5000
chunk_starts <- seq(1, n_cols, by = chunk_size)

cat("  Converting in", length(chunk_starts), "chunks of", chunk_size, "cells each...\n")

sparse_chunks <- vector("list", length(chunk_starts))

for (i in seq_along(chunk_starts)) {
  start_col <- chunk_starts[i]
  end_col   <- min(start_col + chunk_size - 1, n_cols)

  cat("    chunk", i, "of", length(chunk_starts),
      "(cols", start_col, "-", end_col, ") ...\n")

  # Convert only this slice to a plain matrix, then to sparse
  chunk_dense  <- as.matrix(raw_counts[, start_col:end_col, drop = FALSE])
  sparse_chunks[[i]] <- as(chunk_dense, "CsparseMatrix")

  rm(chunk_dense)
  gc(full = TRUE)
}

# Free the original (large) data.frame before assembling the final matrix
rm(raw_counts)
gc(full = TRUE)

cat("  Assembling chunks into one sparse matrix...\n")
raw_counts_sparse <- do.call(cbind, sparse_chunks)
rownames(raw_counts_sparse) <- gene_names
colnames(raw_counts_sparse) <- cell_names

rm(sparse_chunks)
gc(full = TRUE)
cat("  -> sparse matrix ready:", paste(dim(raw_counts_sparse), collapse = " x "), "\n")

cat("  Writing to BPCells on-disk directory...\n")
write_matrix_dir(
  mat = raw_counts_sparse,
  dir = "GSE131907_bpcells",
  compress = TRUE
)

# Free RAM immediately - critical before opening the on-disk handle
rm(raw_counts_sparse)
gc(full = TRUE)
cat("  -> BPCells write complete, RAM freed.\n")

counts_bpcells <- open_matrix_dir(dir = "GSE131907_bpcells")

# ----------------------------------------------------------------------
# 6. Create Full Seurat v5 Object (On-Disk Backed)
# ----------------------------------------------------------------------
cat("Step 3: Creating on-disk Seurat v5 object...\n")

common_cells <- intersect(colnames(counts_bpcells), rownames(raw_metadata))
raw_metadata <- raw_metadata[common_cells, ]
cat("  -> common cells:", length(common_cells), "\n")

seurat_full <- CreateSeuratObject(
  counts       = counts_bpcells[, common_cells],
  meta.data    = raw_metadata,
  project      = "LUAD_GSE131907",
  min.cells    = 3,
  min.features = 200
)
gc(full = TRUE)

# ----------------------------------------------------------------------
# 7. Calculate Mitochondrial QC Metrics & Export Summary
# ----------------------------------------------------------------------
cat("Step 4: Calculating QC metrics...\n")
seurat_full[["percent_mt"]] <- PercentageFeatureSet(seurat_full, pattern = "^MT-")

qc_summary <- seurat_full@meta.data %>%
  group_by(Sample) %>%
  summarise(
    cell_count      = n(),
    mean_nCount     = mean(nCount_RNA),
    median_nFeature = median(nFeature_RNA),
    mean_percent_mt = mean(percent_mt)
  )
write.csv(qc_summary, "GSE131907_sample_QC.csv", row.names = FALSE)

# Save the full on-disk object reference (lightweight - data stays on disk)
saveRDS(seurat_full, "GSE131907_full_on_disk.rds")
gc(full = TRUE)

# ----------------------------------------------------------------------
# 8. Final Inspection + sessionInfo (for reproducibility / Member 5's audit)
#    NOTE: Step 5 (the sketch/exploratory step) now lives in its own file,
#    step5_sketch.R, run separately AFTER this script finishes successfully.
#    See that file for details.
# ----------------------------------------------------------------------
cat("=== STEPS 1-4 COMPLETE (metadata, BPCells, Seurat object, QC) ===\n")
print(seurat_full)
capture.output(sessionInfo(), file = "sessionInfo.txt")
cat("=== Pipeline finished:", format(Sys.time()), "===\n")
cat("Next: run  Rscript step5_sketch.R  to generate the sketch assay.\n")
