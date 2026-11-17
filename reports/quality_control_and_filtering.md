# Quality Control and Gene Expression Filtering

## 1. Sample Quality Control

The dataset contains 72 RNA-seq samples.

Sample metadata was successfully integrated with the expression dataset.

The metadata includes:

- GEO accession
- Library name
- Cultivar
- Tissue
- Treatment
- Time point
- Replicate

One sample, `IR9AD2L1`, was identified with an unusually low library size of 1,755 reads compared with the dataset median of 13,707,163 reads.

The sample was flagged for further investigation but was not automatically removed.

## 2. Expression Data Assessment

The original expression matrix contains:

- 55,986 expression records
- 72 samples

The expression matrix was checked for:

- Numeric values
- Missing values
- Non-finite values
- Negative counts
- Integer count values
- Consistency between sample names and metadata

The validation checks confirmed that the expression data met the required input conditions.

## 3. Filtering Method

A preliminary expression filter was applied using:

- Count threshold: 10
- Minimum samples: 3

Records with counts of at least 10 in at least three samples were retained.

## 4. Filtering Results

| Metric               | Result |
| -------------------- | -----: |
| Original records     | 55,986 |
| Retained records     | 21,910 |
| Removed records      | 34,076 |
| Retention percentage | 39.13% |

The filtered matrix contains 21,910 records across 72 samples.

## 5. Outlier Investigation and Sensitivity Analysis

The sample `IR9AD2L1` was identified as an unusually low-count sample.

A sensitivity analysis was performed for the severe drought Day 9 comparison in latest leaf tissue.

Two analyses were compared:

1. Original analysis using all six samples (3 control and 3 treatment).
2. Sensitivity analysis excluding `IR9AD2L1` (3 control and 2 treatment).

### Sensitivity Analysis Results

| Metric                            | Original analysis | Without outlier |
| --------------------------------- | ----------------: | --------------: |
| Total samples                     |                 6 |               5 |
| Control samples                   |                 3 |               3 |
| Treatment samples                 |                 3 |               2 |
| Features tested                   |            18,700 |          18,700 |
| Significant features (FDR < 0.05) |                37 |           1,014 |

### Significant Feature Overlap

| Result                                | Number of features |
| ------------------------------------- | -----------------: |
| Significant in original analysis      |                 37 |
| Significant without outlier           |              1,014 |
| Significant in both analyses          |                 37 |
| Only significant in original analysis |                  0 |
| Only significant without outlier      |                977 |

All 37 significant features from the original analysis remained significant after excluding the low-count sample. The sensitivity analysis identified an additional 977 significant features.

### Interpretation

The substantial increase in significant features indicates that the Day 9 latest-leaf results are sensitive to the low-count sample.

However, excluding the sample reduces the treatment group to only two biological replicates. This limits statistical reliability.

Therefore, the original analysis including all six samples remains the primary analysis. The sensitivity analysis is retained as an exploratory investigation, and the original results have not been replaced.

## 6. Generated Outputs

The quality control and filtering workflow generated:

- `gene_expression_statistics.csv`
- `gene_expression_filtering_summary.csv`
- `gene_expression_filtering.png`
- `filtered_counts.rds`
- `preprocessing_sample_qc.csv`
- `preprocessing_filtering_summary.csv`
- `technical_feature_id_key.csv`
- `missing_padj_audit.csv`
- `outlier_sensitivity_summary.csv`
- `outlier_sensitivity_comparison.csv`

## 7. Conclusion

The quality control and preliminary filtering processes reduced the expression matrix from 55,986 to 21,910 records while retaining records that met the selected expression threshold.

The unusually low-count sample `IR9AD2L1` was investigated through sensitivity analysis rather than automatically removed. The results demonstrated that its inclusion substantially influences the severe drought Day 9 latest-leaf differential expression findings.

The original dataset and primary differential expression results were preserved for transparent reporting.

A limitation remains: the complete biological gene identifier mapping has not been validated. Therefore, the expression records currently use technical feature identifiers, and gene-specific biological interpretations cannot yet be made.
