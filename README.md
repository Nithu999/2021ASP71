# CropStress Intelligence

### Agricultural Genomics Platform

## Project Overview

CropStress Intelligence is an agricultural genomics project focused on analysing gene-expression patterns associated with crop stress responses.

This project uses publicly available rice RNA-sequencing data to investigate gene-expression changes during drought stress and subsequent recovery.

The long-term vision is to develop a research-support platform that helps agricultural researchers explore crop stress responses using genomic data.

## Research Dataset

- **Source:** NCBI Gene Expression Omnibus (GEO)
- **Dataset:** GSE288615
- **Organism:** Rice (_Oryza sativa_)
- **Cultivar:** IR64
- **Technology:** RNA sequencing (RNA-seq)
- **Samples:** 72
- **Expression records:** 55,986
- **Reference genome:** IRGSP-1.0

## Technologies

- R
- RStudio
- Git and GitHub
- R Markdown
- ggplot2 (planned, if used in subsequent analyses)

## Current Progress

### Dataset Preparation

- Public GEO dataset selected and imported.
- GEO sample metadata extracted and integrated.
- Sample identifiers validated against the count matrix.
- Initial count-matrix validation completed.

### Quality Assessment

- Sample-level library sizes calculated.
- Number of detected expression records calculated.
- Library-size distribution visualised.
- Potential low-count sample identified: `IR9AD2L1`.
- Biological replicate comparison performed.
- Relative library size investigated.
- Count distribution investigated.
- Potential sample-level count anomaly documented.

### Research Documentation

- Initial dataset metadata documentation prepared.
- R Markdown report created.
- Quality-control tables and visualisations generated.

**Important:** The identified sample has not been excluded. Further investigation is required before making any exclusion decision.

## Project Structure

```text
CropStress_Intelligence/
│
├── code/
├── data/
│   ├── raw/
│   └── processed/
│
├── references/
├── reports/
│
└── results/
    ├── figures/
    └── tables/
```

## Project Status

Currently in the dataset validation and initial RNA-seq quality assessment stage.

The next stages will focus on:

1. Investigating potential sample-quality issues.
2. Validating gene identifiers.
3. Performing exploratory data analysis.
4. Applying appropriate normalization and filtering.
5. Preparing for differential gene-expression analysis.

## Long-Term Vision

Develop CropStress Intelligence into a research-support platform for agricultural genomics, with potential applications in crop stress analysis and agricultural decision support.
