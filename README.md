# 🌾 CropStress Intelligence

### An Agricultural Genomics Platform for Understanding Rice Responses to Environmental Stress

## 1. Project Overview

**CropStress Intelligence** is a bioinformatics project focused on analysing rice gene-expression data to investigate how plants respond to environmental stress, particularly severe drought and subsequent rehydration.

The project uses publicly available RNA-sequencing data to explore gene-expression patterns, assess sample quality, investigate gene annotations, and prepare the data for downstream differential expression analysis.

The long-term vision is to develop this work into an agricultural genomics platform that can support data-driven crop research and contribute to agricultural innovation, including applications relevant to Sri Lanka.

## 2. Project Objectives

- Analyse rice RNA-sequencing data associated with drought stress and rehydration.
- Extract and integrate experimental sample metadata.
- Assess RNA-seq sample quality and identify potential outliers.
- Investigate gene identifiers and available genome annotations.
- Perform exploratory data analysis.
- Apply preliminary gene-expression filtering.
- Prepare validated data for future differential expression analysis.

## 3. Dataset Description

This project uses the public dataset **GSE288615**, obtained from the NCBI Gene Expression Omnibus (GEO).

| Attribute                   | Description                        |
| --------------------------- | ---------------------------------- |
| Dataset ID                  | GSE288615                          |
| Organism                    | _Oryza sativa_                     |
| Cultivar                    | IR64                               |
| Data type                   | RNA-seq raw counts                 |
| Number of samples           | 72                                 |
| Original expression records | 55,986                             |
| Tissue types                | Latest leaf and second latest leaf |
| Experimental focus          | Severe drought and rehydration     |

The dataset includes control, severe drought, and rehydration-related experimental conditions.

## 4. Technologies and Tools

- **R** – Bioinformatics data processing and statistical analysis.
- **RStudio** – Development environment.
- **Git and GitHub** – Version control and project management.
- **NCBI GEO** – Public gene-expression dataset source.
- **MSU Rice Genome Annotation** – Gene annotation investigation.

## 5. Project Workflow

The project follows an organised R-based workflow.

| Stage | Analysis                                                                                       |
| ----- | ---------------------------------------------------------------------------------------------- |
| 01    | Extract dataset metadata                                                                       |
| 02    | Import raw RNA-seq count data                                                                  |
| 03    | Integrate sample metadata                                                                      |
| 04    | Perform sample quality control                                                                 |
| 05    | Investigate gene identifiers                                                                   |
| 06    | Investigate gene annotations                                                                   |
| 07    | Perform exploratory data analysis                                                              |
| 08    | Apply preliminary gene-expression filtering                                                    |
| 09–10 | Prepare and validate differential expression inputs                                            |
| 11–17 | Investigate original identifiers, annotation compatibility, and supplementary gene identifiers |

The individual analysis scripts are maintained in the `code/` directory.

## 6. Preliminary Expression Filtering Results

A preliminary expression filter was applied using the following criteria:

- Minimum raw count: 10
- Minimum number of samples meeting the threshold: 3

| Metric                      | Result |
| --------------------------- | -----: |
| Original expression records | 55,986 |
| Retained records            | 21,910 |
| Removed records             | 34,076 |
| Retention percentage        | 39.13% |

The resulting filtered expression matrix contains 21,910 records across 72 samples.

This is a preliminary filtering step and does not represent completed differential expression analysis.

## 7. Gene Annotation Investigation

The project investigates the compatibility of the original expression data with available rice genome annotation resources.

The MSU Rice Genome Annotation GFF3 file was examined.

A comparison identified that all 588 gene identifiers extracted from supplementary Table S2 were present in the investigated GFF3 annotation.

However, the original count matrix does not contain an explicit gene identifier column. Therefore, the correspondence between individual count-matrix rows and the GFF3 gene identifiers has not yet been established.

Gene identifiers have not been assigned to the complete count matrix based solely on row order.

## 8. Sample Quality Control

Sample quality assessment was performed using library-size and expression-related statistics.

One sample, `IR9AD2L1`, was identified with an unusually low library size and requires further investigation.

The sample has not been automatically removed from the dataset.

## 9. Repository Structure

```text
CropStress_Intelligence/
│
├── code/
│   ├── 00_run_all.R
│   ├── 01_extract_metadata.R
│   ├── 02_import_counts.R
│   ├── 03_integrate_metadata.R
│   ├── 04_sample_quality_control.R
│   ├── 05_gene_identifier_inspection.R
│   ├── 06_gene_annotation_inspection.R
│   ├── 07_exploratory_data_analysis.R
│   ├── 08_gene_expression_filtering.R
│   ├── 09_prepare_differential_expression.R
│   ├── 10_validate_differential_expression_inputs.R
│   └── ... additional investigation scripts
│
├── data/
│   ├── raw/
│   └── processed/
│
├── references/
│
├── reports/
│
├── results/
│   ├── figures/
│   └── tables/
│
├── CropStress_Intelligence.Rproj
├── .gitignore
└── README.md
```

## 10. Reports and Outputs

The repository contains:

- Dataset and methodology documentation.
- Initial sample quality control report.
- Quality control and filtering documentation.
- Exploratory data analysis figures.
- Gene-expression statistics and filtering summaries.
- Integrated sample metadata.
- Gene annotation investigation results.
- Supplementary gene identifier comparison results.

## 11. Current Limitations and Future Work

The following tasks remain for further development:

1. Establish reliable gene identifier mapping for the original count matrix.
2. Investigate the unusually low-count sample.
3. Complete differential expression analysis using an appropriate statistical method.
4. Perform biological interpretation of differentially expressed genes.
5. Investigate relevant biological pathways and drought-response mechanisms.

**Important:** Differential expression analysis and biological interpretation have not yet been completed.

## 12. Reproducibility

The project is organised into individual R scripts to support transparent and reproducible analysis.

To work with the project:

1. Clone or download this repository.
2. Open `CropStress_Intelligence.Rproj` in RStudio.
3. Ensure the required raw dataset and reference files are available in their respective directories.
4. Set the project root as the working directory.
5. Run the relevant scripts in the intended workflow order.

The `00_run_all.R` script executes the initial seven workflow stages. Later filtering and investigation scripts are maintained separately and should be executed as required.

## 13. Research Data Source

**NCBI Gene Expression Omnibus (GEO):** GSE288615

Dataset: Public rice RNA-seq expression data for drought and rehydration-related research.

---

_CropStress Intelligence – Exploring Agricultural Genomics Through Bioinformatics._
