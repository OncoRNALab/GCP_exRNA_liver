# Data Directory Documentation

This folder contains supplementary data files required to reproduce the figures and analyses in the main project. Each file is described below:

---

## File Descriptions

- **Suplementary_Table1.csv**: Sample annotation for cohort 1
- **Suplementary_Table2.csv**: Sample annotation for cohort 2
- **counts.zip**: Raw count matrices for cohorts 1 and 2, unzip before running analyses
- **CPMFilt300_NashRisk.csv**: CPM count matrix for NASH subcohort samples
- **dedup.csv**: Percentage of reads remaining after UMI deduplication (cohorts 1 and 2)
- **Ngenes.csv**: Number of protein-coding genes with counts above 10 (cohorts 1 and 2)
- **OutliersC1.txt**, **OutliersC2.txt**: List of sample IDs classified as outliers based on z-scores of 5 QC metrics
- **prctStar.csv**: Number and percentage of uniquely mapped reads from STAR alignment
- **Reads_Number.csv**: Number of reads at each step of the preprocessing workflow
- **SpikeRatio_C1.csv**, **SpikeRatio_C2.csv**: Endogenous and spike-in counts (cohorts 1 and 2)
- **Strand.csv**: Percentage of reads aligned to the correct strand
- **CodeBook_SuplementaryTable1.xlsx**: Contains the description of variables included in Suplementary_Table1.csv (cohort 1 sample annotation)
- **CodeBook_SuplementaryTable2.xlsx**: Contains the description of variables included in Suplementary_Table2.csv (cohort 2 sample annotation)

---

## Usage Notes

- These files are referenced by analysis scripts in the main repository.
- Ensure you have the correct file paths when running scripts.
- For more details on each file, refer to the main project README or script comments.

---

## Contact

For questions or issues regarding these data files, please contact the repository maintainer.


