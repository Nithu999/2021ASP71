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

One sample, IR9AD2L1, was identified with an unusually low library size and requires further investigation.

## 2. Expression Data Assessment

The original expression matrix contains:

- 55,986 expression records
- 72 samples

The expression matrix was checked for numeric values and non-negative counts.

## 3. Filtering Method

A preliminary expression filter was applied using:

- Count threshold: 10
- Minimum samples: 3

## 4. Filtering Results

| Metric               | Result |
| -------------------- | -----: |
| Original records     | 55,986 |
| Retained records     | 21,910 |
| Removed records      | 34,076 |
| Retention percentage | 39.13% |

The filtered matrix contains 21,910 records across 72 samples.

## 5. Generated Outputs

- `gene_expression_statistics.csv`
- `gene_expression_filtering_summary.csv`
- `gene_expression_filtering.png`
- `filtered_counts.rds`

## 6. Conclusion

The preliminary filtering process reduced the expression matrix while retaining records meeting the selected expression threshold.

Further differential expression analysis and biological interpretation are required.
