# 🌐 Jazyk / Language
* [Čeština](#čeština)
* [English](#english)

---

# CZ: Diferenciace druhů rodu *Penicillium* (série *Viridicata*) za použití MALDI-TOF MS a strojového učení

Tento repozitář obsahuje bioinformatický skript vyvinutý během zahraniční výzkumné stáže Erasmus+ ve **Westerdijk Fungal Biodiversity Institute** (Utrecht, Nizozemsko). Projekt je zaměřen na automatizovanou druhovou diferenciaci uvnitř taxonomicky náročné série *Penicillium* *Viridicata* s využitím spektrálních dat z MALDI-TOF hmotnostní spektrometrie a pokročilých statistických modelů v jazyce R.

---

## 🔬 Přehled projektu

### Problém
Tradiční morfologická identifikace blízce příbuzných druhů rodu *Penicillium* je extrémně časově náročná a vyžaduje vysokou expertízu. Sekvenování DNA sice poskytuje přesné výsledky, ale pro rutinní diagnostiku s vysokým objemem vzorků je drahé a pomalé.

### Řešení
Hmotnostní spektrometrie MALDI-TOF nabízí rychlou, přesnou a nákladově efektivní alternativu založenou na jedinečných proteinových profilech ("otiscích prstů"). Komerční softwary (jako uzavřený MBT Compass Explorer od Brukeru) však mají fixně nastavené parametry a často selhávají při čistém rozlišení vizuálně podobných druhových komplexů, což vede k nepřehledným a roztříštěným dendrogramům.

### Tento projekt
K překonání limitů komerčního softwaru byl vytvořen tento vlastní bioinformatický kód v R za využití open-source knihoven. Ořezáním surových spekter na optimální sektory, implementací striktnějšího zarovnání píků a nasazením algoritmů strojového učení se podařilo drasticky snížit chybovost klasifikace (OOB error rate) ze **40,68 %** na pouhých **8,7 %**.

---

## 📊 Metodika a pracovní postup

1. **Mikrobiologie a příprava vzorků:** 
   * Kultivace 61 kmenů ze sbírky CBS/DTO na MEA agaru a Sabouraudově bujónu.
   * Vysokokapacitní proteinová extrakce za využití kyseliny mravenčí a acetonitrilu.
2. **Akvizice spekter a kontrola kvality:**
   * Automatické měření v programu `flexControl` (Bruker) – generování 32 spekter na jeden vzorek.
   * Manuální filtrace plochých spekter a odlehlých hodnot ve `flexAnalysis` (dosažení ≥ 20 kvalitních spekter na kmen).
3. **Bioinformatické zpracování (Skript v R):**
   * **Preprocessing:** Ořezání spektrálního rozsahu (2500–4700 Da), normalizace pomocí celkového iontového proudu (TIC) a průměrování spekter.
   * **Extrakce znaků:** Detekce píků (metoda MAD, SNR = 3) a striktní binování píků (tolerance = 0,0009).
   * **Statistická analýza:** Kruskal-Wallisův test s Benjamini-Hochbergovou korekcí falešné pozitivity pro identifikaci druhově specifických biomarkerů.
   * **Strojové učení a vícerozměrná statistika:** Hierarchické shlukování (Ward.D2), PCA, t-SNE a klasifikace pomocí Random Forest (500 stromů).

---

## 🛠️ Struktura repozitáře

* `penicillium_analysis.R` - Kompletní R skript obsahující preprocessing dat, identifikaci biomarkerů, PCA, t-SNE a Random Forest modelování.
* `README.md` - Dokumentace projektu.

*Poznámka: Surová data z hmotnostní spektrometrie (`~/R/MSP ready`) jsou uložena lokálně z důvodu institucionální ochrany dat a nejsou součástí tohoto repozitáře.*

---

## 📦 Systémové požadavky a R balíčky

### Verze R
* **R (v4.4.3)** – Tato specifická verze byla použita k zajištění plné zpětné kompatibility a stability s potřebnými bioinformatickými a hmotnostně-spektrometrickými balíčky.

### Potřebné knihovny
Pro instalaci všech balíčků nutných pro zpracování dat, strojové učení a vizualizaci spusťte v R konzoli následující příkaz:

```R
install.packages(c("readBrukerFlexData", "MALDIquant", "MALDIquantForeign", 
                   "MALDIrppa", "pheatmap", "Rtsne", "randomForest", 
                   "ggplot2", "ggrepel"))
```

<br>
<br>
<br>

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

<br>

# EN: Differentiation of *Penicillium* Species (Series *Viridicata*) Using MALDI-TOF MS & Machine Learning

This repository contains the bioinformatic workflow developed during an Erasmus+ research internship at the **Westerdijk Fungal Biodiversity Institute** (Utrecht, Netherlands). The project focuses on automated species differentiation within the challenging *Penicillium* series *Viridicata* using MALDI-TOF mass spectrometry data and custom R-based statistical/machine learning pipelines.

---

## 🔬 Project Overview

### The Problem
Traditional morphological identification of closely related *Penicillium* species is highly time-consuming and requires extreme expertise. While DNA sequencing provides accurate results, it is often too expensive and slow for routine high-throughput diagnostics. 

### The Solution
MALDI-TOF Mass Spectrometry offers a rapid, cost-effective alternative based on unique protein fingerprints. However, commercial software (such as Bruker's closed-source MBT Compass Explorer) has fixed parameters and often fails to cleanly differentiate visually similar species complexes, leading to messy, scattered dendrograms.

### This Project
To overcome these software limitations, this custom R pipeline was developed using open-source packages. By trimming raw spectra to optimal sectors, implementing strict peak alignment, and deploying machine learning models, the classification error rate (OOB) was drastically reduced from **40.68%** down to **8.7%**.

---

## 📊 Pipeline & Methodology

1. **Microbiology & Sample Prep:** 
   * Cultivation of 61 strains from the CBS/DTO collection on MEA agar and Sabouraud Broth.
   * High-throughput protein extraction using a formic acid/acetonitrile protocol.
2. **Spectra Acquisition & Quality Control:**
   * Automated runs via `flexControl` (Bruker), generating 32 spectra per sample.
   * Manual filtering of flatlines and outliers in `flexAnalysis` (ensuring ≥ 20 high-quality spectra per strain).
3. **Bioinformatic Processing (The R Script):**
   * **Preprocessing:** Mass spectral trimming (2500–4700 Da), TIC normalization, and averaging.
   * **Feature Extraction:** Peak detection (MAD method, SNR = 3) and strict peak binning (tolerance = 0.0009).
   * **Statistical Analysis:** Kruskal-Wallis testing with Benjamini-Hochberg false-positive correction to detect highly specific species biomarkers.
   * **Machine Learning & Multivariate Stats:** Hierarchical clustering (Ward.D2), PCA, t-SNE, and Random Forest classification (500 trees).

---

## 🛠️ Repository Structure

* `penicillium_analysis.R` - The complete R script containing data preprocessing, biomarker identification, PCA, t-SNE, and Random Forest modeling.
* `README.md` - Project documentation.

*Note: Raw mass spectrometry data files (`~/R/MSP ready`) are stored locally due to institutional data restrictions and are not included in this repository.*

---

## 📦 System Requirements & R Packages

### R Version
* **R (v4.4.3)** – This specific version was used to ensure full backwards compatibility and stability with the required bioinformatic and mass spectrometry packages.

### Required Packages
To install all the necessary packages for preprocessing, machine learning, and visualization, run the following command in your R console:

```R
install.packages(c("readBrukerFlexData", "MALDIquant", "MALDIquantForeign", 
                   "MALDIrppa", "pheatmap", "Rtsne", "randomForest", 
                   "ggplot2", "ggrepel"))
```
