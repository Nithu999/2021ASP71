# Differential Expression Analysis

## 1. Objective

Differential expression analysis was performed to identify expression records showing differences between control and treatment conditions in rice under severe drought and rehydration.

## 2. Methodology

The analysis was conducted using DESeq2 in R.

The following settings were used:

- Statistical package: DESeq2
- Design formula: `~ condition`
- Size-factor estimation: `poscounts`
- Significance threshold: adjusted p-value (FDR) < 0.05
- Effect-size threshold for the stricter summary: absolute log2 fold change >= 1
- Differential expression comparisons: eight separate comparisons across two leaf tissues.

The comparisons included:

1. Severe drought Day 6
2. Severe drought Day 9
3. Rehydration Day 3
4. Rehydration Day 6

Each comparison initially included three control and three treatment samples.

## 3. Differential Expression Results

| Comparison           | Tissue             | Significant features (FDR < 0.05) | FDR < 0.05 and abs(log2FC) >= 1 |
| -------------------- | ------------------ | --------------------------------: | ------------------------------: |
| Severe drought Day 6 | Latest leaf        |                             2,183 |                             609 |
| Severe drought Day 9 | Latest leaf        |                                37 |                              37 |
| Rehydration Day 3    | Latest leaf        |                                 1 |                               1 |
| Rehydration Day 6    | Latest leaf        |                                60 |                              39 |
| Severe drought Day 6 | Second latest leaf |                               811 |                             447 |
| Severe drought Day 9 | Second latest leaf |                             1,099 |                             568 |
| Rehydration Day 3    | Second latest leaf |                                 7 |                               6 |
| Rehydration Day 6    | Second latest leaf |                               821 |                             299 |

## 4. Missing Adjusted P-values

An audit identified missing adjusted p-values in some comparisons.

A diagnostic analysis investigated independent filtering and Cook's distance handling.

The results indicated that independent filtering accounted for most missing adjusted p-values in the severe drought Day 9 latest-leaf and Day 6 second-latest-leaf comparisons.

Disabling filtering and Cook's cutoff was used only for diagnostic purposes. The original DESeq2 settings and results were preserved.

## 5. Outlier Sensitivity Analysis

The severe drought Day 9 latest-leaf comparison was investigated further because sample `IR9AD2L1` had an unusually low library size.

| Metric               | Original analysis | Without outlier |
| -------------------- | ----------------: | --------------: |
| Samples              |                 6 |               5 |
| Control samples      |                 3 |               3 |
| Treatment samples    |                 3 |               2 |
| Features tested      |            18,700 |          18,700 |
| Significant features |                37 |           1,014 |

All 37 significant features in the original analysis remained significant in the sensitivity analysis. An additional 977 features were significant only after excluding the sample.

This demonstrates that the results are sensitive to the low-count sample. However, the reduced treatment group contains only two replicates. Therefore, the sensitivity results are exploratory, and the original analysis remains the primary result.

## 6. Limitations

- The complete mapping between expression records and biological gene identifiers has not been validated.
- Results currently use technical feature identifiers rather than confirmed gene identifiers.
- Biological interpretation at the individual gene level is therefore not possible.
- Some adjusted p-values are missing because of DESeq2 filtering and outlier handling.
- The sensitivity analysis has only two treatment replicates after sample exclusion.

## 7. Conclusion

DESeq2 identified different numbers of significant expression features across drought and rehydration conditions and between leaf tissues.

The results provide a statistical basis for further investigation of rice responses to drought and rehydration. However, biological interpretation requires a validated mapping between expression records and gene identifiers.
