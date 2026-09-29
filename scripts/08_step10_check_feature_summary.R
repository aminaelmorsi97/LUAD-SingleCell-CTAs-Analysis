# ============================================================================
# Step 10 - Check GSE131907_Feature_Summary.xlsx for patient/sample mapping
# Run from Terminal:  Rscript step10_check_feature_summary.R
# This is a light, fast check - safe to run anytime.
# ============================================================================

if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
library(readxl)

setwd("~/Downloads/LUAD_Project")

path <- "GSE131907_Feature_Summary.xlsx"

cat("=== Sheets in the file ===\n")
print(excel_sheets(path))

cat("\n=== Reading first sheet ===\n")
df <- read_excel(path, sheet = 1)

cat("\n=== Column names ===\n")
print(colnames(df))

cat("\n=== First 10 rows ===\n")
print(head(df, 10))

cat("\n=== Dimensions ===\n")
cat("Rows:", nrow(df), " Columns:", ncol(df), "\n")
