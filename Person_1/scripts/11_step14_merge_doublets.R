# ============================================================================
# Step 14 - Merge doublet status into the main metadata file
# Run from Terminal:  Rscript step14_merge_doublets.R
# Fast - just a table join, no heavy computation.
# ============================================================================

cat("=== Step 14 started:", format(Sys.time()), "===\n")

library(dplyr)

setwd("~/Downloads/LUAD_Project")

meta <- read.csv("GSE131907_metadata_curated.csv", row.names = 1)
doublets <- read.csv("GSE131907_doublet_status.csv")

cat("Metadata rows:", nrow(meta), "\n")
cat("Doublet results rows:", nrow(doublets), "\n")

meta$cell_id <- rownames(meta)

merged <- meta %>%
  left_join(doublets %>% select(cell_id, doublet_score, doublet_class),
            by = "cell_id")

cat("\nMerged rows:", nrow(merged), "\n")
cat("Cells with missing doublet info:", sum(is.na(merged$doublet_class)), "\n")

rownames(merged) <- merged$cell_id
write.csv(merged, "GSE131907_metadata_curated.csv", row.names = TRUE)

cat("\n=== STEP 14 COMPLETE ===\n")
cat("GSE131907_metadata_curated.csv now includes doublet_score and doublet_class.\n")
cat("=== Step 14 finished:", format(Sys.time()), "===\n")
