# ============================================================================
# Step 12 - Generate checksums for all key deliverable files
# Run from Terminal:  Rscript step12_checksums.R
# Fast - just hashes existing files, no heavy computation.
# ============================================================================

if (!requireNamespace("digest", quietly = TRUE)) install.packages("digest")
library(digest)

setwd("~/Downloads/LUAD_Project")

files_to_check <- c(
  "GSE131907_metadata_curated.csv",
  "GSE131907_sample_QC.csv",
  "GSE131907_QC_filtering_report.csv",
  "GSE131907_QC_percentiles.csv",
  "GSE131907_sample_to_patient_mapping.csv",
  "GSE131907_full_on_disk.rds",
  "GSE131907_filtered_on_disk.rds",
  "GSE131907_malignant_full.rds",
  "GSE131907_malignant_portable.rds",
  "GSE131907_sketch.rds"
)

# Only check files that actually exist (some may not be present on every
# machine, e.g. the full raw ones)
files_to_check <- files_to_check[file.exists(files_to_check)]

cat("Computing checksums (md5) for", length(files_to_check), "files...\n")

checksums <- data.frame(
  file = files_to_check,
  size_MB = round(file.info(files_to_check)$size / 1e6, 2),
  md5 = sapply(files_to_check, function(f) digest(f, algo = "md5", file = TRUE)),
  row.names = NULL
)

print(checksums)

write.csv(checksums, "CHECKSUMS.csv", row.names = FALSE)
cat("\nSaved: CHECKSUMS.csv\n")
cat("Anyone can verify their copy with: digest::digest(\"filename\", algo = \"md5\", file = TRUE)\n")
