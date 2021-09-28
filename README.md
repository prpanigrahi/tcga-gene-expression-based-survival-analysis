# tcga-gene-expression-based-survival-analysis
- R snippet to do survival analysis using Gene Expression cohort data from TCGA breast cancer project.
- Both R script and Rmd file provided. Rmd file is used to generate pdf report.

# About Data
- GeneExpression_log2Counts_TNBC_HER2pos_ERpos_Normal_normalized.txt
  - Gene expression data (61 genes, 1070 patients)
- sample_group_info.txt
  - Sample information data (For 1070 patients, their race, ethnicity, gender, ERpos/HER2pos/Normal/TNBC)
  - Samples are divided into 4 groups ERpos:805, HER2pos:37, Normal:112, TNBC:116
- TCGA-BRCA_clinical.csv
  - Clinical data (For 1070 patients, clinical information, data required for survival analysis

# About the script
- Script takes Gene Expression data of 61 key genes of 1070 patients grouped into 4 class. Using clinical information and survival data, the script plots survival plot of 4 subtypes.


