# ============================================================================
# Step 6a - QC Visualization (BEFORE choosing filtering thresholds)
# Run from Terminal with:  Rscript step6a_qc_visualize.R
# Produces plots + a distribution summary so thresholds can be chosen
# deliberately, not guessed.
# ============================================================================

cat("=== Step 6a started:", format(Sys.time()), "===\n")

rm(list = ls())
gc(full = TRUE)

suppressPackageStartupMessages({
  library(Seurat)
  library(BPCells)
  library(ggplot2)
  library(dplyr)
})

setwd("~/Downloads/LUAD_Project")

# ----------------------------------------------------------------------
# 1. Load the on-disk Seurat object (lightweight - counts stay on disk)
# ----------------------------------------------------------------------
cat("Loading on-disk Seurat object...\n")
seurat_full <- readRDS("GSE131907_full_on_disk.rds")
cat("  -> loaded:", ncol(seurat_full), "cells x", nrow(seurat_full), "genes\n")

meta <- seurat_full@meta.data

# ----------------------------------------------------------------------
# 2. Numeric summary (percentiles) - this is what actually helps pick
#    thresholds, more than eyeballing a plot alone.
# ----------------------------------------------------------------------
cat("\n--- nCount_RNA percentiles ---\n")
print(quantile(meta$nCount_RNA, probs = c(0.01, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99)))

cat("\n--- nFeature_RNA percentiles ---\n")
print(quantile(meta$nFeature_RNA, probs = c(0.01, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99)))

cat("\n--- percent_mt percentiles ---\n")
print(quantile(meta$percent_mt, probs = c(0.5, 0.75, 0.90, 0.95, 0.975, 0.99, 0.995, 0.999)))

# Save the summary to a file too, so it's easy to share with the team
qc_percentiles <- data.frame(
  metric = rep(c("nCount_RNA", "nFeature_RNA", "percent_mt"), each = 9),
  percentile = rep(c(0.01, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99), 3),
  value = c(
    quantile(meta$nCount_RNA, probs = c(0.01, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99)),
    quantile(meta$nFeature_RNA, probs = c(0.01, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99)),
    quantile(meta$percent_mt, probs = c(0.01, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99))
  )
)
write.csv(qc_percentiles, "GSE131907_QC_percentiles.csv", row.names = FALSE)

# ----------------------------------------------------------------------
# 3. Plots - saved to one PDF, one page per metric + a combined scatter
# ----------------------------------------------------------------------
cat("\nGenerating QC plots...\n")

pdf("GSE131907_QC_plots.pdf", width = 10, height = 7)

# Histograms
p1 <- ggplot(meta, aes(x = nCount_RNA)) +
  geom_histogram(bins = 100, fill = "steelblue") +
  scale_x_log10() +
  labs(title = "nCount_RNA distribution (log10 x-axis)", x = "nCount_RNA (UMIs/cell)")
print(p1)

p2 <- ggplot(meta, aes(x = nFeature_RNA)) +
  geom_histogram(bins = 100, fill = "darkorange") +
  labs(title = "nFeature_RNA distribution", x = "nFeature_RNA (genes/cell)")
print(p2)

p3 <- ggplot(meta, aes(x = percent_mt)) +
  geom_histogram(bins = 100, fill = "firebrick") +
  xlim(0, 50) +
  labs(title = "percent_mt distribution (capped at 50% for readability)", x = "% mitochondrial reads")
print(p3)

# Scatter: nCount vs nFeature, colored by percent_mt - classic QC diagnostic
p4 <- ggplot(meta, aes(x = nCount_RNA, y = nFeature_RNA, color = percent_mt)) +
  geom_point(alpha = 0.2, size = 0.3) +
  scale_x_log10() +
  scale_color_gradient(low = "grey80", high = "red", limits = c(0, 30)) +
  labs(title = "nCount vs nFeature, colored by % mito")
print(p4)

# Per-sample violin plots (if there aren't too many samples, this can be busy -
# check GSE131907_sample_QC.csv first if this looks too crowded)
if ("Sample" %in% colnames(meta) && length(unique(meta$Sample)) <= 40) {
  p5 <- VlnPlot(seurat_full, features = "nCount_RNA", group.by = "Sample", pt.size = 0) +
    theme(axis.text.x = element_text(angle = 90, size = 6)) + NoLegend()
  print(p5)

  p6 <- VlnPlot(seurat_full, features = "percent_mt", group.by = "Sample", pt.size = 0) +
    theme(axis.text.x = element_text(angle = 90, size = 6)) + NoLegend()
  print(p6)
}

dev.off()

cat("\n=== STEP 6a COMPLETE ===\n")
cat("Review these files before deciding on filtering thresholds:\n")
cat("  - GSE131907_QC_plots.pdf\n")
cat("  - GSE131907_QC_percentiles.csv\n")
cat("=== Step 6a finished:", format(Sys.time()), "===\n")
