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
