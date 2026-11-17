# Exploratory Data Analysis Report

## CropStress Intelligence

### 1. Introduction

Exploratory Data Analysis (EDA) was performed to understand the structure, experimental design, and distribution of the GSE288615 rice RNA-seq dataset.

The analysis focused on experimental group distribution, sampling conditions, raw count validation, and gene-level expression statistics.

### 2. Dataset Overview

| Attribute                   | Value                           |
| --------------------------- | ------------------------------- |
| Dataset                     | GSE288615                       |
| Organism                    | Oryza sativa                    |
| Cultivar                    | IR64                            |
| Original expression records | 55,986                          |
| Number of samples           | 72                              |
| Data type                   | Raw RNA-seq counts              |
| Tissue types                | Latest leaf, Second latest leaf |

The original count matrix contains 55,986 expression records across 72 samples.

### 3. Experimental Group Distribution

The dataset contains four treatment categories:

- Control (severe drought)
- Severe drought
- Control (rehydration)
- Rehydration after severe drought

Two tissue types were investigated:

- Latest leaf
- Second latest leaf

The sampling days include:

- Day 0
- Day 0 (2 hours)
- Day 3
- Day 6
- Day 9

The experimental group summary identified 24 treatment-day-tissue combinations. Each combination contains three samples.

This provides a total of 72 samples.

**Observation:** The experimental design includes separate severe drought and rehydration experiments. Their treatment groups and sampling times should be considered separately during downstream analysis.

### 4. Raw Count Matrix Validation

The raw count matrix was examined to verify its structure and numerical validity.

The following checks were performed:

- Confirmed that the count data were stored in a matrix or data-frame structure.
- Verified that count values were numeric.
- Checked for missing or infinite values.
- Confirmed that all count values were non-negative.
- Calculated overall count distribution statistics.
- Investigated zero and non-zero count values.

The validation checks completed successfully.

### 5. Gene-Level Expression Statistics

Gene-level statistics were calculated for the original count matrix.

The following measurements were generated:

| Statistic     | Description                             |
| ------------- | --------------------------------------- |
| Mean count    | Average expression count across samples |
| Median count  | Median expression count across samples  |
| Maximum count | Highest observed count                  |
| Zero samples  | Number of samples with zero counts      |

The results were saved in:

`results/tables/gene_count_statistics.csv`

These statistics provide an initial understanding of expression variability across the dataset.

### 6. Visualisations

The following visualisations were generated during the analysis.

**Experimental Group Distribution**

File: `results/figures/experimental_group_distribution.png`

This figure illustrates the distribution of samples across treatment groups and sampling days.

**Library Size Distribution**

File: `results/figures/library_size_distribution.png`

This figure supports the investigation of sequencing depth across samples.

**Detected Records by Sample**

File: `results/figures/detected_records_by_sample.png`

This figure illustrates the number of detected expression records across individual samples.

### 7. Important Observations

1. The dataset contains 72 samples and 55,986 original expression records.
2. The experimental design includes two tissue types and four treatment categories.
3. Each of the 24 observed treatment-day-tissue combinations contains three samples.
4. The raw count matrix passed the numerical validation checks.
5. Gene-level expression statistics were generated for further investigation.
6. Sample-level quality assessment identified an unusually low-count sample, `IR9AD2L1`, which requires further investigation.

### 8. Limitations

- The original count matrix does not contain an explicit gene identifier column.
- The EDA stage does not establish biological significance.
- No differential expression analysis was performed at this stage.
- The unusually low-count sample requires further investigation before final downstream analysis.
- The experimental groups must be handled according to their respective experimental designs.

### 9. Conclusion

The exploratory data analysis established an initial understanding of the GSE288615 dataset, including its experimental structure, sample distribution, and raw count characteristics.

The generated tables and figures provide a foundation for subsequent quality control, expression filtering, and differential expression analysis.

**Status:** Exploratory Data Analysis completed.
