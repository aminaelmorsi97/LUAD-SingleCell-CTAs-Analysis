# Person 1 Pipeline - Data Engineering, QC & Identity (GSE131907 CTA Project)

 person-1
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
=======
## 🧬 Project Overview

This repository contains the bioinformatic analysis pipeline for investigating
**Cancer-Testis Antigen (CTA)** expression dynamics during lung adenocarcinoma
(LUAD) metastatic progression, using single-cell RNA-seq data.

**Research question:** How does the prevalence and heterogeneity of
Cancer-Testis Antigens differ between patients, transcriptional states, and
disease sites in malignant cells of lung adenocarcinoma?

## 📊 Dataset Summary

- **Primary dataset (GSE131907):** 208,506 single cells, 58 samples, 44 patients
- **Malignant subset:** 31,136 cells across primary tumor and metastatic sites
  (`tLung`, `mLN`, `mBrain`, `PE`), labeled `Malignant cells` (metastatic
  sites) and `tS1`/`tS2`/`tS3` (primary tumor transcriptional states)
- **Quality control:** thresholds taken from the original publication's
  Methods (percent.mt ≤ 20%, nCount 100–150,000, nFeature 200–10,000)

## 👥 Team & Responsibilities

Work is divided by scientific function, on one shared pipeline, so that
discovery results can be applied independently to validation cohorts.

| # | Role | Branch | Key deliverable |
|---|------|--------|------------------|
| 1 | Data engineering, QC & identity | `person-1` | Portable Seurat object, Sample↔Patient_id mapping, QC report, doublet status |
| 2 | Malignant cell identification & CNV validation | `person-2` (`malignancy/cnv`) | Independent malignancy calls (CopyKAT/inferCNV), per-patient CNV sensitivity set |
| 3 | Transcriptional states & plasticity | `person-3` (`states/plasticity`) | Malignant state scores (EMT, Stemness, Plasticity, Proliferation, Lung differentiation) |
| 4 | CTA prevalence & statistics | `person-4` (`cta/pseudobulk`) | Patient-level CTA prevalence/DE, candidate ranking |
| 5 | External validation & reproducibility | `person-5` (`validation/gse123902`) | Independent replication in GSE123902 / GSE198291 |

## 📁 Repository Structure

```
person-1/   -> scripts/  (data engineering & QC pipeline, run in numbered order)
person-2/   -> CNV validation scripts and per-patient results
person-3/   -> malignant state scoring scripts
person-4/   -> CTA prevalence, pseudobulk, and DE analysis + final report
person-5/   -> external cohort validation scripts
```

Each branch is developed independently and merged via Pull Request after
review by another team member, per the team's quality rules (see below).


## ✅ Binding Quality Rules

- One approved pipeline; no member creates a separate parallel workflow
- Patient (not cell) is the unit of statistical inference
- Raw counts are used for pseudobulk; normalized data is for display only
- Discovery (original annotation) is kept separate from sensitivity analysis
  (inferCNV/CopyKAT)
- All thresholds, exclusions, package versions, and file names are documented
  without silent changes
- No result is called final before review by a second team member

## 🚀 How to Reproduce

See each branch's own README for step-by-step instructions. The data
engineering pipeline (`person-1/scripts/00_run_all.R`) is the required
starting point before any other branch's analysis can run.

## 📖 Reference

Kim N, Kim HK, Lee K, et al. Single-cell RNA sequencing demonstrates the
molecular and cellular reprogramming of metastatic lung adenocarcinoma.
*Nat Commun*. 2020;11:2285.

## 👥 Contributors

- Amina Elmorsi — Data Engineering & QC (Person 1)
- Nada Hossam - CNV validation scripts and per-patient results (Person 2)
- Abdallah Ali - malignant state scoring scripts (Person 3)
- Munawar Hraib - CTA prevalence, pseudobulk, and DE analysis (Person 4)
- Dilawer Chofan - external cohort validation scripts (Person 5)
main
