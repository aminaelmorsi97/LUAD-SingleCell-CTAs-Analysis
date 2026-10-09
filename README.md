
## 🧬 Project Overview

This repository contains the bioinformatic analysis pipeline for investigating
**Cancer-Testis Antigen (CTA)** expression dynamics during lung adenocarcinoma
(LUAD) metastatic progression, using single-cell RNA-seq data.

**Research question:** How does the prevalence and heterogeneity of
Cancer-Testis Antigens differ between patients, transcriptional states, and
disease sites in malignant cells of lung adenocarcinoma?

## 📊 Dataset Summary

- **Primary discovery dataset (GSE131907):** 208,506 single cells, 58 samples, 44 patients
- **Malignant subset:** 31,136 cells across primary tumor and metastatic sites
  (`tLung`, `mLN`, `mBrain`, `PE`), labeled `Malignant cells` (metastatic
  sites) and `tS1`/`tS2`/`tS3` (primary tumor transcriptional states)
- **External validation cohorts:** 
  - `GSE123902` (analyzed by Amina Elmorsi)
  - `GSE198291` (analyzed by Abdallah Ali)
- **Quality control:** thresholds taken from the original publication's
  Methods (percent.mt ≤ 20%, nCount 100–150,000, nFeature 200–10,000)

## 📦 Data Access & Large File Downloads

Due to file size limitations on GitHub, large processing files and Seurat objects are hosted on Google Drive:

| File Name | Description | Download Link |
|-----------|-------------|---------------|
| `GSE131907_malignant_portable.rds` | Processed Seurat object of discovery malignant cells (~227 MB) | [Download from Google Drive](https://drive.google.com/file/d/1SHXKbHVjJeI5AjAIoCOSTctK2i3wAZ5g/view?usp=share_link) | 
> **Note:** Make sure to download these `.rds` files and place them in the working data directory before running downstream analysis scripts.

## 👥 Team & Responsibilities

Work is divided by scientific function on a unified pipeline, integrating discovery-stage single-cell analyses with independent validation across external cohorts.

| # | Team Member | Branch | Key Deliverable & Role |
|---|-------------|--------|------------------------|
| 1 | **Amina Elmorsi** | `person-1` (`data-qc/gse123902`) | Data engineering & QC pipeline, portable Seurat objects, plus independent validation on **GSE123902** |
| 2 | **Nada Hossam** | `person-2` (`malignancy/cnv`) | Independent malignancy calls (CopyKAT/inferCNV) & per-patient CNV sensitivity set |
| 3 | **Abdallah Ali** | `person-3` (`states/gse198291`) | Malignant state scoring (EMT, Stemness, Proliferation) & external validation on **GSE198291** |
| 4 | **Munawar Hraib** | `person-4` (`cta/synthesis`) | CTA prevalence/DE stats, candidate prioritization, pipeline integration, and manuscript finalization |

## 📁 Repository Structure
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

- **Amina Elmorsi** — Data Engineering, QC & External Validation (`GSE123902`)
- **Nada Hossam** — CNV Validation & Sensitivity Analysis (`inferCNV`)
- **Abdallah Ali** — Malignant State Scoring & External Validation (`GSE198291`)
- **Munawar Hraib** — CTA Analysis, Synthesis & Manuscript Finalization
