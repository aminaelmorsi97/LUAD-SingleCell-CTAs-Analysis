# ==============================================================================
# MASTER PIPELINE - Person 1: Data Engineering, QC & Identity
# GSE131907 CTA Project
#
# This script runs the full Person 1 pipeline end-to-end by sourcing each
# step script in order. Each step is also kept as its own file (in this same
# folder) so it can be re-run individually if something needs to be fixed or
# re-checked, without having to re-run the whole thing from scratch.
#
# HOW TO RUN:
#   From Terminal (recommended - some steps are memory/time intensive):
#     cd scripts/
#     Rscript 00_run_all.R
#
#   Or run each numbered script individually in the order shown below.
#
# NOT INCLUDED HERE: CNV validation scripts (CopyKAT/inferCNV) - those belong
# to Person 2's task/branch (malignancy/cnv), not the data engineering/QC
# scope of this branch (qc/gse131907).
# ==============================================================================

steps <- c(
  "01_pipeline.R",                       # Load raw data, build on-disk BPCells object, basic QC metrics
  "02_step5_sketch.R",                   # Generate 50k-cell sketch for fast exploration
  "03_step6a_qc_visualize.R",            # Plot QC metric distributions (before choosing thresholds)
  "04_step6b_qc_filter.R",               # Apply QC filtering (thresholds from Kim et al. 2020)
  "05_step8a_check_metadata.R",          # Inspect metadata columns/categories (diagnostic)
  "06_step7_malignant_full.R",           # Extract malignant cells (Malignant cells + tS1/tS2/tS3) from full filtered data
  "07_step9_make_portable_object.R",     # Build a fully standalone (non-BPCells) malignant object
  "08_step10_check_feature_summary.R",   # Inspect Feature_Summary.xlsx (diagnostic, for patient mapping)
  "09_step11_sample_patient_mapping.R",  # Build verified Sample -> Patient_id 1-to-1 mapping table
  "10_step13_doublet_detection.R",       # Per-sample doublet detection (scDblFinder)
  "11_step14_merge_doublets.R",          # Merge doublet status into the master metadata file
  "12_step12_checksums.R"                # Compute checksums for all deliverable files (run LAST, after any file changes)
)

cat("==============================================================\n")
cat("Person 1 master pipeline -", length(steps), "steps\n")
cat("==============================================================\n\n")

for (s in steps) {
  cat("\n>>> Running:", s, " -", format(Sys.time()), "\n")
  cat("--------------------------------------------------------------\n")
  source(s, echo = FALSE)
  cat("--------------------------------------------------------------\n")
  cat(">>> Finished:", s, "\n")
}

cat("\n==============================================================\n")
cat("ALL STEPS COMPLETE -", format(Sys.time()), "\n")
cat("==============================================================\n")
