# Dataset Metadata – GSE288615

## 1. Dataset Identification

- **GEO Accession:** GSE288615
- **Dataset Title:** Biomarker genes for model-based prediction of drought-stress perception levels in rice
- **Organism:** Oryza sativa (Rice)
- **Cultivar:** IR64
- **Data Type:** RNA-seq gene-expression raw counts
- **Number of Samples:** 72
- **Number of Expression Records:** 55,986
- **Data Source:** NCBI Gene Expression Omnibus (GEO)
- **Official Dataset URL:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE288615

## 2. Biological Objective

The original study investigated drought-stress perception in rice using time-series transcriptome and phenotypic analysis. It examined changes in gene expression during severe drought and after rehydration.

## 3. Experimental Overview

Rice cultivar IR64 was subjected to severe drought treatment for nine days, followed by irrigation for up to six days.

The study included control and drought-related experimental conditions, multiple collection time points, and two leaf positions.

## 4. Sequencing Information

- Sequencing technology: RNA sequencing.
- Instrument: Illumina NovaSeq 6000.
- Reference genome: IRGSP-1.0.
- Alignment software: STAR.
- Read-count quantification: featureCounts v2.0.0.

## 5. Initial Data Quality Observations

Initial inspection identified an unusually low-count sample:

- Sample ID: IR9AD2L1
- Total counts: 1,755
- Detected genes: 1,046

Official GEO metadata identifies this as a latest-leaf sample collected under severe drought conditions on Day 9, biological replicate 2.

The sample has been flagged for further quality-control investigation. No exclusion decision has been made.

## 6. Data Source and References

NCBI GEO accession GSE288615.

Sample record: GSM8771263.

The original study and associated metadata should be consulted when interpreting experimental conditions and biological results.
