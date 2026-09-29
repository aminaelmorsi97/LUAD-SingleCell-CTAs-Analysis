# Person 1 Pipeline - Data Engineering, QC & Identity (GSE131907 CTA Project)

Branch: `qc/gse131907` (per team convention - see project plan, section 7)

## How to run

**Option A - run everything in order:**
```bash
cd scripts/
Rscript 00_run_all.R
```

**Option B - run one step at a time** (recommended if you're re-running after a
crash, or only need to redo one part): run the numbered files in order,
`01_...` through `12_...`. Each file's header comment explains what it does
and why.

Run from a terminal (`Rscript file.R`), not by pasting into RStudio - several
steps process the full 208,506-cell dataset and are safer run this way (see
step file comments for details).

## Step order and purpose

| # | Script | Purpose | Key output |
|---|--------|---------|------------|
| 01 | `01_pipeline.R` | Load raw counts + metadata, build on-disk BPCells-backed Seurat v5 object (avoids loading the full matrix into RAM), compute per-cell QC metrics | `GSE131907_full_on_disk.rds`, `GSE131907_bpcells/` |
| 02 | `02_step5_sketch.R` | Generate a 50,000-cell leverage-score sketch for fast exploratory work (NOT used for final statistics - see note below) | `GSE131907_sketch.rds` |
| 03 | `03_step6a_qc_visualize.R` | Plot QC metric distributions before choosing filtering thresholds | `GSE131907_QC_plots.pdf`, `GSE131907_QC_percentiles.csv` |
| 04 | `04_step6b_qc_filter.R` | Apply QC filtering using thresholds from the original paper's Methods (percent.mt ≤ 20%, nCount 100-150,000, nFeature 200-10,000) - not arbitrary cutoffs | `GSE131907_filtered_on_disk.rds`, `GSE131907_QC_filtering_report.csv` |
| 05 | `05_step8a_check_metadata.R` | Diagnostic: print metadata columns and `Cell_subtype` categories before writing downstream logic | (console output) |
| 06 | `06_step7_malignant_full.R` | Extract malignant cells from the full, QC-filtered dataset. **Important**: the original paper labels primary-tumor malignant cells as `tS1`/`tS2`/`tS3` and metastatic-site malignant cells as `Malignant cells` - all four categories are included (31,136 cells total) | `GSE131907_malignant_full.rds`, `GSE131907_bpcells_malignant/` |
| 07 | `07_step9_make_portable_object.R` | Materialize the malignant-cell counts into a standard in-memory matrix and save one self-contained `.rds` file - no BPCells folder needed downstream | `GSE131907_malignant_portable.rds` |
| 08 | `08_step10_check_feature_summary.R` | Diagnostic: inspect `GSE131907_Feature_Summary.xlsx` for a Sample/Patient mapping | (console output) |
| 09 | `09_step11_sample_patient_mapping.R` | Build and verify the Sample -> Patient_id mapping (1-to-1, confirmed programmatically, not assumed) | `GSE131907_sample_to_patient_mapping.csv` |
| 10 | `10_step13_doublet_detection.R` | Per-sample doublet detection (scDblFinder) - doublets form within a single sample's library prep, so detection runs per sample, same principle as per-patient CNV analysis | `GSE131907_doublet_status.csv` |
| 11 | `11_step14_merge_doublets.R` | Merge doublet_score/doublet_class into the master metadata file (single source of truth) | Updated `GSE131907_metadata_curated.csv` |
| 12 | `12_step12_checksums.R` | Compute MD5 checksums for all deliverable files, so teammates can verify their copy matches exactly. Run this LAST, and again any time a file changes | `CHECKSUMS.csv` |

## Why the sketch (step 02) is not used for final statistics

The team's own quality rule: *"no single random/sketch subsample as basis
for DE or prevalence estimates."* The sketch is kept only for fast
exploratory work (e.g. quick UMAP checks); all final malignant-cell
extraction and downstream analysis uses the full, QC-filtered dataset
(step 06 onward).

## Not included in this branch

CNV validation scripts (CopyKAT/inferCNV) belong to the `malignancy/cnv`
branch (Person 2's task), not this data-engineering/QC scope.

## Requirements

```r
install.packages(c("Seurat", "Matrix", "dplyr", "ggplot2", "readxl", "digest"))
remotes::install_github("bnprks/BPCells/r")
BiocManager::install("scDblFinder")
```

## Reference

Kim N, Kim HK, Lee K, et al. Single-cell RNA sequencing demonstrates the
molecular and cellular reprogramming of metastatic lung adenocarcinoma.
*Nat Commun*. 2020;11:2285.
