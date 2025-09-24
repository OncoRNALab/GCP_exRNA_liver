---
# <center>GCP - Plasma cell‐Free transcriptome profiling in chronic liver disease</center> 
---

**This repository contains all the scripts and supplementary files required for reproducing the analysis and figures reported in [add paper citation] as part of the Grand Challenge Project: Overcoming the main current diagnostic challenges in hepatology practice. ⚠️ Given that our cohorts were originally defined using NAFLD/NASH criteria, we continue to use these terms throughout this manuscript. Readers should be aware that these correspond to MASLD and MASH under current nomenclature. ⚠️**

---

## Data Availability

- **Raw sequencing data:** [EGA](https://ega-archive.org), accession EGAD50000000775
- **Raw count tables:** [BioStudies Database](http://www.ebi.ac.uk/biostudies), accession S-BSST1614, and in `./Data/counts.zip`
- **Preprocessing pipelines:**
  - [ExtractUMI-picoV3-and-RNA-Exome-UMI-pipeline](https://github.com/OncoRNALab/ExtractUMI-picoV3-and-RNA-Exome-UMI-pipeline.git)
  - [circRNA](https://github.com/OncoRNALab/circRNA.git)
  - [totalexRNA (in development)](https://github.com/jasperverwilt/totalexRNA.git)

---

## Repository Structure

```text
.
├── CircRNA/           # Supplementary data for circRNA analysis
├── Data/              # Supplementary data for total RNA analysis (Diff. Abundance, GSEA)
│   ├── Suplementary_Table1.csv, Suplementary_Table2.csv   # Sample annotation for cohort 1 and 2
│   ├── CodeBook_SuplementaryTable1.xlsx, CodeBook_SuplementaryTable2.xlsx   # Clinical annotation descriptions
│   ├── counts.zip     # Cohort 1 and 2 raw count matrices
│   ├── OutliersC1.txt, OutliersC2.txt   # Outlier samples based on QC metrics
├── Deconv/            # Supplementary data for deconvolution analysis
├── Figures/           # Output figures from all analyses
├── *.Rmd, *.R         # Analysis scripts
├── README.md          # Project documentation
└── GCP_GitRepo.Rproj  # R project file
```

---

## Scripts

- **01_QC_plots.Rmd**: QC-metric analysis and box plots
- **02_Z_scores.Rmd**: Calculate z-scores for QC-metrics to identify outliers
- **circRNA_analysis.R**: circRNA differential abundance analysis
- **circRNAs_Volcano.Rmd**: circRNA volcano plots
- **03_Cohort1_Analysis_Cirrhosis.Rmd**: Differential abundance, volcano plots, and GSEA for cohort 1
- **04_Cohort2_Analysis_NashRisk.Rmd**: Differential abundance, volcano plots, and GSEA for cohort 2
- **05_EndoConc_NAFLD_C1.Rmd**: RNA concentration analysis for cohort 1
- **06_EndoConc_Fibrosis_C2.Rmd**: RNA concentration analysis for cohort 2
- **07_EndoConc_NashRisk_C2.Rmd**: RNA concentration analysis for cohort 2 (Nash risk)
- **CleanAnn.Rmd**: Data cleaning and annotation

---

## Figures

All output figures from the analyses are available in the `Figures/` directory:
- **Fig1.pdf** to **Fig7.pdf**: Main figures
- **ROC_cirrhosis_C1.png**, **ROC_NAFLD_C1.png**: ROC curves

---

## Supplementary Data

- **CircRNA/**: Data for circRNA analysis
- **Data/**: Data for total RNA analysis, sample annotation, clinical annotation, raw counts, and outlier lists
- **Deconv/**: Data for deconvolution analysis

---

## How to Reproduce the Analysis

1. Clone this repository
2. Download raw sequencing data and count tables from the provided links
3. Raw RNAseq FASTQ files can be processed using the above-cited preprocessing pipelines
4. Use the provided scripts and supplementary data to reproduce the analyses and figures
5. Refer to the pipelines for preprocessing steps

---

## Citation

Please cite the associated paper: [add paper citation]

---

## Contact

For questions or issues, please contact the repository maintainer.



