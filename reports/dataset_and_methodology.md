# CropStress Intelligence
## Dataset and Methodology

### 1. Dataset Description

This project uses the public rice RNA-seq dataset GSE288615 obtained from NCBI Gene Expression Omnibus (GEO).

- Organism: Oryza sativa
- Cultivar: IR64
- Original expression records: 55,986
- Samples: 72
- Tissue types: Latest leaf and second latest leaf.

### 2. Data Processing

The analysis workflow includes:

1. Dataset import and validation.
2. Sample metadata extraction.
3. Metadata integration.
4. Sample quality assessment.
5. Gene identifier investigation.
6. Gene annotation inspection.
7. Exploratory data analysis.
8. Preliminary gene expression filtering.

### 3. Gene Expression Filtering

A preliminary filtering threshold was applied:

- Minimum raw count: 10
- Minimum number of samples: 3

Results:

| Metric | Value |
|---|---:|
| Original records | 55,986 |
| Retained records | 21,910 |
| Removed records | 34,076 |
| Retention | 39.13% |

### 4. Annotation Investigation

The MSU Rice Genome Annotation GFF3 file was investigated.

All 588 gene identifiers extracted from supplementary Table S2 were found in the GFF3 annotation.

However, the original count matrix does not contain explicit gene identifiers. The row-by-row correspondence between the count matrix and GFF3 annotation has not been established.

### 5. Current Limitations

- Gene identifier mapping requires further verification.
- Differential expression analysis has not yet been completed.
- Biological interpretation remains future work.

### 6. Reproducibility

R scripts are maintained in the `code/` directory. Processed data, figures and tables are stored in their respective project directories.