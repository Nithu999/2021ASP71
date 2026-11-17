# Gene Annotation Investigation Report

## CropStress Intelligence

### 1. Introduction

Gene annotation is important in RNA-seq analysis because it connects gene-expression measurements with identifiable biological genes.

This investigation examined the availability and compatibility of gene identifiers in the GSE288615 rice RNA-seq dataset.

The investigation focused on the original count matrix, MSU Rice Genome Annotation GFF3 file, and supplementary gene identifiers.

### 2. Original Count Matrix

The original RNA-seq count matrix was examined to identify whether gene identifiers were included.

| Attribute                       |         Result |
| ------------------------------- | -------------: |
| Dataset                         |      GSE288615 |
| Organism                        |   Oryza sativa |
| Cultivar                        |           IR64 |
| Expression records              |         55,986 |
| Samples                         |             72 |
| Explicit gene identifier column | Not identified |

The original count matrix contains expression counts across 72 samples. However, an explicit gene identifier column was not identified during file inspection.

Therefore, the biological identity of each individual count-matrix row has not yet been established.

### 3. MSU Rice Genome Annotation Investigation

The reference annotation file examined was:

`references/osa1_r7.all_models.gff3.gz`

The GFF3 file was inspected to identify gene records and extract their identifiers.

| Attribute               | Result |
| ----------------------- | -----: |
| Total gene records      | 55,986 |
| Valid gene identifiers  | 55,986 |
| Unique gene identifiers | 55,986 |
| Missing identifiers     |      0 |
| Duplicate identifiers   |      0 |

The annotation file contains 55,986 unique gene identifiers.

### 4. Supplementary Gene Identifier Validation

The supplementary Excel workbook was also investigated:

`references/12870_2025_7175_MOESM2_ESM.xlsx`

The gene identifiers from supplementary Table S2 were compared with the identifiers extracted from the MSU7 GFF3 annotation.

| Attribute                      | Result |
| ------------------------------ | -----: |
| Supplementary gene identifiers |    588 |
| Matching identifiers           |    588 |
| Unmatched identifiers          |      0 |
| Matching percentage            |   100% |

All 588 supplementary gene identifiers were found in the investigated GFF3 annotation.

This demonstrates compatibility between the supplementary identifiers and the investigated annotation file for these 588 identifiers.

### 5. Count Matrix and Annotation Comparison

The original count matrix contains 55,986 expression records, while the investigated GFF3 annotation contains 55,986 gene records.

Although the dimensions match, this does not establish that the count-matrix rows correspond to the GFF3 gene identifiers in the same order.

The original count file does not provide an explicit gene identifier column that can be used to verify this correspondence.

Therefore, assigning all GFF3 identifiers to the count matrix based only on row order would not be scientifically justified.

### 6. Important Findings

1. The original count matrix contains 55,986 expression records across 72 samples.
2. An explicit gene identifier column was not identified in the original count matrix.
3. The MSU7 GFF3 annotation contains 55,986 unique gene identifiers.
4. All 588 supplementary Table S2 identifiers were found in the GFF3 annotation.
5. The matching record counts do not prove row-by-row correspondence between the original count matrix and the GFF3 annotation.

### 7. Limitations

The main limitation is the absence of explicit gene identifiers in the original count matrix.

Further investigation is required to establish a reliable mapping between expression records and biological gene identifiers.

The matching supplementary identifiers cannot be used to assign identifiers to all 55,986 count-matrix records.

### 8. Future Work

The following activities are recommended:

- Investigate the original dataset processing information and source files.
- Identify the exact reference annotation and gene identifier version used during read quantification.
- Locate a corresponding featureCounts output or other source containing the original gene identifiers.
- Establish and validate a reliable mapping before performing gene-level biological interpretation.

### 9. Conclusion

The investigation confirmed that the MSU7 GFF3 annotation contains 55,986 unique gene identifiers and that all 588 supplementary gene identifiers are compatible with this annotation.

However, the original count matrix does not contain an explicit gene identifier column. Consequently, a reliable mapping between all count-matrix records and the annotation has not yet been established.

This limitation must be resolved before making gene-specific biological interpretations.

**Status:** Gene annotation investigation completed; complete count-matrix gene mapping remains unresolved.
