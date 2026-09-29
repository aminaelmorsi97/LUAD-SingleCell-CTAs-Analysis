# ============================================================================
# Step 8a - Check metadata columns needed for CNV analysis
# Run from Terminal:  Rscript step8a_check_metadata.R
# This just prints information - no heavy computation, safe to run anywhere.
# ============================================================================

rm(list = ls())
suppressPackageStartupMessages({
  library(Seurat)
})

setwd("~/Downloads/LUAD_Project")

seurat_full <- readRDS("GSE131907_filtered_on_disk.rds")
meta <- seurat_full@meta.data

cat("=== Column names in metadata ===\n")
print(colnames(meta))

cat("\n=== Cell_subtype categories (need to identify immune cells as reference) ===\n")
print(table(meta$Cell_subtype))

# Try to guess a patient/sample identifier column
cat("\n=== Likely sample/patient columns (first few unique values) ===\n")
for (col in c("Sample", "Patient", "PatientID", "sample_id", "patient_id", "orig.ident")) {
  if (col %in% colnames(meta)) {
    cat("\n--", col, "--\n")
    print(head(unique(meta[[col]]), 15))
  }
}
