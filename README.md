# 🌐 Language / Jazyk
* [English](#english)
* [Čeština](README_cs.md)

---

# English
# Differentiation of *Penicillium* Species (Series *Viridicata*) Using MALDI-TOF MS & Machine Learning

This repository contains the bioinformatic workflow developed during an Erasmus+ research internship at the **Westerdijk Fungal Biodiversity Institute** (Utrecht, The Netherlands). The project focuses on automated species differentiation within the challenging *Penicillium* series *Viridicata* using MALDI-TOF mass spectrometry data and custom R-based statistical/machine learning pipelines.

---

## Project Overview

### The Problem
Traditional morphological identification of closely related *Penicillium* species is highly time-consuming and requires extreme expertise. While DNA sequencing provides accurate results, it is often too expensive and slow for routine high-throughput diagnostics. 

### The Solution
MALDI-TOF Mass Spectrometry offers a rapid, cost-effective alternative based on unique protein fingerprints. However, commercial software (such as Bruker's closed-source MBT Compass Explorer) has fixed parameters and often fails to cleanly differentiate visually similar species complexes, leading to messy, scattered dendrograms.

### This Project
To overcome these software limitations, this custom R pipeline was developed using open-source packages. By trimming raw spectra to optimal sectors, implementing strict peak alignment, and deploying machine learning models, the classification error rate (OOB) was drastically reduced from **40.68%** down to **8.7%**.

---

## Pipeline & Methodology

1. **Microbiology & Sample Prep:** 
   * Cultivation of strains from the CBS/DTO collection on MEA agar and Sabouraud Broth.
   * High-throughput protein extraction using a formic acid/acetonitrile protocol.
2. **Spectra Acquisition & Quality Control:**
   * Automated runs via `flexControl` (Bruker), generating 32 spectra per sample.
   * Manual filtering of flatlines and outliers in `flexAnalysis` (ensuring ≥ 20 high-quality spectra per strain).
3. **Bioinformatic Processing (The R Script):**
   * **Preprocessing:** Mass spectra trimming (2500–4700 Da), TIC normalization, and averaging.
   * **Feature Extraction:** Peak detection (MAD method, SNR = 3) and strict peak binning (tolerance = 0.0009).
   * **Statistical Analysis:** Kruskal-Wallis testing with Benjamini-Hochberg false-positive correction to detect highly specific species biomarkers.
   * **Machine Learning & Multivariate Stats:** Hierarchical clustering (Ward.D2), PCA, t-SNE, and Random Forest classification (500 trees).

---

## Data Protection & Anonymization (Important)

Due to intellectual property protection and internal regulations of the Westerdijk Institute, **raw mass spectrometry files are not included in this repository**.

To keep the pipeline fully functional and executable for demonstration purposes, a **mathematically modified and fully anonymized benchmark dataset** (`data_matrix.rds`) is provided, from which all hidden metadata has been stripped:
* It contains a subset of data (8 model groups labeled anonymously from `sample1` to `sample8`).
* All identifiers have been removed from the object attributes.
* The peak intensity values have been mathematically altered (random noise added) to prevent back-identification against reference databases, while maintaining data geometry for R computing.

Due to the absence of raw spectra, the initial preprocessing sections of `penicillium_maldi_analysis.R` are commented out. The script execution begins directly by loading the cleaned matrix (however, you need to load the R packages and set the working directory):

```R
data_matrix <- readRDS("data_matrix.rds")
```

---

## Repository Structure

* `R/` - Folder containing source code and processed data.
  * `data_matrix.rds` - Preprocessed data matrix ready for analysis.
  * `penicillium_maldi_analysis.R` - The complete R script containing data preprocessing, biomarker identification, PCA, t-SNE, and Random Forest modeling.
* `LICENSE` - Project license terms.
* `README.md` - Project documentation in English.
* `README_cs.md` - Project documentation in Czech.

---

## System Requirements & R Packages

### R Version
* **R (v4.4.3)** – This specific version was used to ensure full backward compatibility and stability with the required bioinformatic and mass spectrometry packages.

### Required Packages
To install all the necessary packages for preprocessing, machine learning, and visualization, run the following command in your R console:

```R
install.packages(c("readBrukerFlexData", "MALDIquant", "MALDIquantForeign", 
                   "MALDIrppa", "pheatmap", "Rtsne", "randomForest", 
                   "ggplot2", "ggrepel"))
```