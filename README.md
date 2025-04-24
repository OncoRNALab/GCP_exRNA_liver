---
# <center>GCP - Plasma cell‐Free transcriptome profiling in chronic liver disease</center> 
---

**<div style="text-align: justify">This repository contains all the scripts and supplementary files required for reproducing the analysis and figures reported in [add paper citation] as part of the Grand Challenge Project: Overcoming the main current diagnostic challenges in hepatology practice.</div><br>**

<div style="text-align: justify">

The raw sequencing data has been deposited in the [EGA](https://ega-archive.org) portal under accession number EGAD50000000775. 

Raw count tables have been deposited in the [BioStudies Database](http://www.ebi.ac.uk/biostudies) under accession number S-BSST1614, and in this repository in the directory ./data/counts.zip 

The total RNA and circRNA preprocessing pipelines can be accessed via GitHub https://github.com/OncoRNALab/ExtractUMI-picoV3-and-RNA-Exome-UMI-pipeline.git  and https://github.com/OncoRNALab/circRNA.git  respectively.</div><br>

### Description

The repository is structured as follows:

```
.
├── CircRNA           # supplementary data needed to run circRNA analysis
├── DATA              # supplementary data needed to run total RNA analysis: Diff. Abundance, GSEA 
    ├── Suplementary_Table 1 and 2 # Sample annotation for cohort 1 and 2
    ├── CodeBook_ChemicalData      # description of clinical annotation in supl. files
    ├── Counts.zip    # Cohort1 and 2 raw count matrices.  
├── Deconv            # supplementary data needed to run Deconvolution analysis
├── Figures           # Output figures from all the analyses 

```
**Scripts**

- **QC_plots.Rmd** - script to do QC-metric analysis and reproduce QC-metrics box plots

- **Z_scores.Rmd** - script to calculate z scores for each QC-metric to identify outliers  

- **20210806_circRNA_analysis.R** - Script to perform the circRNA differential abundance analysis.

- **circRNAs_Volcano.Rmd** - Script to reproduce the circRNA volcano plots.

- **Cohort1_Analysis_Cirrhosis.Rmd** - Script to perform differential abundance analysis, volcano plots and GSEA for cohort 1 samples. 

- **Cohort2_Analysis_NashRisk.Rmd** - Script to perform differential abundance analysis, volcano plots and GSEA for cohort 2 samples.

- **EndoConc_NAFLD_C1.Rmd** - Script to perform RNA concentration analysis for cohort 1 samples 

- **EndoConc_Fibrosis_C2.Rmd** - Script to perform RNA concentration analysis for cohort 2 samples 



