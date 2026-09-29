# ============================================================================
# Step 11 - Build the Sample -> Patient_id mapping table (for Dr. Manwar)
# Run from Terminal:  Rscript step11_sample_patient_mapping.R
#
# SOURCE: GSE131907_Feature_Summary.xlsx (official GEO supplementary file)
# This file has a blank first row and the real column headers on row 2,
# so we skip row 1 when reading it.
# ============================================================================

cat("=== Step 11 started:", format(Sys.time()), "===\n")

if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
library(readxl)
library(dplyr)

setwd("~/Downloads/LUAD_Project")

# ----------------------------------------------------------------------
# 1. Read the file, skipping the blank first row so row 2 becomes headers
# ----------------------------------------------------------------------
raw <- read_excel("GSE131907_Feature_Summary.xlsx", sheet = 1, skip = 1)

cat("Column names found:\n")
print(colnames(raw))

cat("\nFirst few rows:\n")
print(head(raw))

# ----------------------------------------------------------------------
# 2. Build the Sample -> Patient mapping
#    (column names below assume "Patient" and "Sample" based on what we
#    saw - adjust if the printed column names above differ slightly)
# ----------------------------------------------------------------------
stopifnot("Patient id" %in% colnames(raw), "Samples" %in% colnames(raw))

mapping_df <- raw %>%
  select(Samples, `Patient id`) %>%
  rename(Sample = Samples, Patient_id = `Patient id`) %>%
  filter(!is.na(Sample), !is.na(Patient_id)) %>%
  distinct()

cat("\n=== Mapping table preview ===\n")
print(mapping_df)

# ----------------------------------------------------------------------
# 3. CRITICAL CHECK: verify each Sample maps to exactly ONE patient
#    (Dr. Manwar's explicit requirement - 1-to-1 mapping)
# ----------------------------------------------------------------------
dup_check <- mapping_df %>%
  count(Sample) %>%
  filter(n > 1)

if (nrow(dup_check) > 0) {
  cat("\n!! WARNING: the following samples map to MORE than one patient:\n")
  print(dup_check)
  cat("Do NOT send this file until this is resolved - check the source rows.\n")
} else {
  cat("\n-> Verified: every Sample maps to exactly one Patient_id. Good.\n")
}

# ----------------------------------------------------------------------
# 4. Compare against the samples we actually have malignant cells for
# ----------------------------------------------------------------------
meta <- read.csv("GSE131907_metadata_curated.csv", row.names = 1)
our_samples <- unique(meta$Sample)

missing_from_mapping <- setdiff(our_samples, mapping_df$Sample)
if (length(missing_from_mapping) > 0) {
  cat("\n!! NOTE: these Sample values from our data were NOT found in the\n")
  cat("   Feature Summary mapping table (check naming differences, e.g.\n")
  cat("   'LUNG_T06' vs 'T06', case sensitivity, etc.):\n")
  print(missing_from_mapping)
} else {
  cat("\n-> All samples in our dataset are covered by this mapping table.\n")
}

# ----------------------------------------------------------------------
# 5. Save
# ----------------------------------------------------------------------
write.csv(mapping_df, "GSE131907_sample_to_patient_mapping.csv", row.names = FALSE)

cat("\n=== STEP 11 COMPLETE ===\n")
cat("Saved: GSE131907_sample_to_patient_mapping.csv\n")
cat("=== Step 11 finished:", format(Sys.time()), "===\n")
