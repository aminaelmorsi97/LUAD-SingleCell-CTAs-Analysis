# ============================================================================
# Step 13 - Doublet detection per sample (scDblFinder)
# Run from Terminal:  Rscript step13_doublet_detection.R
#
# WHY PER SAMPLE: doublets form during library prep of a single 10x run,
# so detection must happen within each sample independently, not on the
# pooled dataset (same principle as the per-sample CNV analysis).
#
# ONE-TIME SETUP (run once in R/RStudio):
#   if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#   BiocManager::install("scDblFinder")
# ============================================================================

cat("=== Step 13 started:", format(Sys.time()), "===\n")

rm(list = ls())
gc(full = TRUE)

suppressPackageStartupMessages({
  library(Seurat)
  library(BPCells)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(dplyr)
})

setwd("~/Downloads/LUAD_Project")

# ----------------------------------------------------------------------
# 1. Load the full QC-filtered object (need ALL cell types per sample,
#    not just malignant, since doublets can be malignant-immune pairs etc.)
# ----------------------------------------------------------------------
cat("Loading filtered on-disk Seurat object...\n")
seurat_full <- readRDS("GSE131907_filtered_on_disk.rds")
meta <- seurat_full@meta.data

all_samples <- sort(unique(meta$Sample))
cat("Total samples to process:", length(all_samples), "\n")

dir.create("doublet_results", showWarnings = FALSE)
results_list <- list()

# ----------------------------------------------------------------------
# 2. Loop per sample
# ----------------------------------------------------------------------
for (this_sample in all_samples) {

  out_file <- file.path("doublet_results", paste0(this_sample, "_doublets.csv"))
  if (file.exists(out_file)) {
    cat("Sample", this_sample, "already done, skipping.\n")
    results_list[[this_sample]] <- read.csv(out_file)
    next
  }

  cat("\n---", this_sample, "---", format(Sys.time()), "\n")
  cells_this_sample <- rownames(meta)[meta$Sample == this_sample]
  cat("  cells:", length(cells_this_sample), "\n")

  if (length(cells_this_sample) < 50) {
    cat("  -> too few cells, skipping doublet detection for this sample\n")
    next
  }

  tryCatch({
    counts_sub <- as(seurat_full[["RNA"]]$counts[, cells_this_sample], "dgCMatrix")

    sce <- SingleCellExperiment(assays = list(counts = counts_sub))
    sce <- scDblFinder(sce, verbose = FALSE)

    result_df <- data.frame(
      cell_id        = colnames(sce),
      Sample         = this_sample,
      doublet_score  = sce$scDblFinder.score,
      doublet_class  = sce$scDblFinder.class   # "singlet" or "doublet"
    )

    write.csv(result_df, out_file, row.names = FALSE)
    results_list[[this_sample]] <- result_df

    cat("  -> doublets found:", sum(result_df$doublet_class == "doublet"),
        "/", nrow(result_df), "\n")

    rm(counts_sub, sce)
    gc(full = TRUE)

  }, error = function(e) {
    cat("  !! ERROR on sample", this_sample, ":", conditionMessage(e), "\n")
  })
}

# ----------------------------------------------------------------------
# 3. Combine all results and merge into the master metadata
# ----------------------------------------------------------------------
cat("\nCombining all sample results...\n")
all_doublets <- bind_rows(results_list)
write.csv(all_doublets, "GSE131907_doublet_status.csv", row.names = FALSE)

cat("\n=== Overall doublet summary ===\n")
print(table(all_doublets$doublet_class))

cat("\n=== STEP 13 COMPLETE ===\n")
cat("Saved: GSE131907_doublet_status.csv (cell_id, Sample, doublet_score, doublet_class)\n")
cat("This can be merged into GSE131907_metadata_curated.csv by cell_id.\n")
cat("=== Step 13 finished:", format(Sys.time()), "===\n")
