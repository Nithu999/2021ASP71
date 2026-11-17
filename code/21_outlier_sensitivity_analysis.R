# ============================================================
# CropStress Intelligence
# Stage 20.3: Outlier Sensitivity Analysis
# ============================================================
# Compare Day 9 severe drought results:
# 1. Original analysis: all six samples
# 2. Sensitivity analysis: excluding IR9AD2L1
#
# Original Stage 19 results will not be modified.
# ============================================================

cat("\nSTAGE 20.3: OUTLIER SENSITIVITY ANALYSIS\n")

library(DESeq2)

# STEP 1: Load data

counts <- readRDS(
  "data/processed/cleaned_filtered_counts.rds"
)

metadata <- read.csv(
  "results/tables/integrated_sample_metadata.csv",
  stringsAsFactors = FALSE
)

# STEP 2: Select Day 9 latest leaf samples

selected <- metadata[
  metadata$tissue == "latest leaf" &
    as.character(metadata$day) == "9" &
    metadata$treatment %in% c(
      "control (severe drought)",
      "severe drought"
    ),
  ,
  drop = FALSE
]

selected$condition <- ifelse(
  selected$treatment == "control (severe drought)",
  "Control",
  "Treatment"
)

selected$condition <- factor(
  selected$condition,
  levels = c("Control", "Treatment")
)

cat("\nOriginal sample groups:\n")
print(table(selected$condition))

# STEP 3: Confirm outlier

outlier_sample <- "IR9AD2L1"

if (!outlier_sample %in% selected$library_name) {
  stop("IR9AD2L1 was not found in the selected comparison.")
}

# STEP 4: Function to run DESeq2

run_deseq <- function(sample_metadata, counts_matrix) {

  sample_metadata <- as.data.frame(sample_metadata)

  rownames(sample_metadata) <- sample_metadata$library_name

  sample_metadata <- sample_metadata[
    colnames(counts_matrix),
    ,
    drop = FALSE
  ]

  dds <- DESeqDataSetFromMatrix(
    countData = round(counts_matrix),
    colData = sample_metadata,
    design = ~ condition
  )

  dds <- estimateSizeFactors(
    dds,
    type = "poscounts"
  )

  dds <- DESeq(
    dds,
    quiet = TRUE
  )

  results(
    dds,
    contrast = c("condition", "Treatment", "Control"),
    alpha = 0.05
  )
}

# STEP 5: Original analysis with all samples

original_counts <- counts[
  ,
  selected$library_name,
  drop = FALSE
]

original_keep <- rowSums(original_counts >= 10) >= 3

original_counts <- original_counts[
  original_keep,
  ,
  drop = FALSE
]

original_results <- run_deseq(
  selected,
  original_counts
)

# STEP 6: Sensitivity analysis without outlier

selected_without_outlier <- selected[
  selected$library_name != outlier_sample,
  ,
  drop = FALSE
]

cat("\nSample groups after excluding outlier:\n")
print(table(selected_without_outlier$condition))

sensitivity_counts <- counts[
  ,
  selected_without_outlier$library_name,
  drop = FALSE
]

sensitivity_keep <- rowSums(
  sensitivity_counts >= 10
) >= 3

sensitivity_counts <- sensitivity_counts[
  sensitivity_keep,
  ,
  drop = FALSE
]

sensitivity_results <- run_deseq(
  selected_without_outlier,
  sensitivity_counts
)

# STEP 7: Compare results

original_df <- as.data.frame(original_results)

original_df$feature_id <- rownames(original_df)

sensitivity_df <- as.data.frame(sensitivity_results)

sensitivity_df$feature_id <- rownames(sensitivity_df)

comparison <- merge(
  original_df,
  sensitivity_df,
  by = "feature_id",
  suffixes = c("_original", "_without_outlier")
)

# STEP 8: Calculate summary

summary_table <- data.frame(

  analysis = c(
    "Original: all samples",
    "Sensitivity: without IR9AD2L1"
  ),

  sample_count = c(
    nrow(selected),
    nrow(selected_without_outlier)
  ),

  control_count = c(
    sum(selected$condition == "Control"),
    sum(selected_without_outlier$condition == "Control")
  ),

  treatment_count = c(
    sum(selected$condition == "Treatment"),
    sum(selected_without_outlier$condition == "Treatment")
  ),

  features_tested = c(
    nrow(original_df),
    nrow(sensitivity_df)
  ),

  significant_features = c(
    sum(original_df$padj < 0.05, na.rm = TRUE),
    sum(sensitivity_df$padj < 0.05, na.rm = TRUE)
  )
)

# STEP 9: Save outputs

output_dir <- "results/tables/differential_expression/stage20"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  summary_table,
  file.path(output_dir, "outlier_sensitivity_summary.csv"),
  row.names = FALSE
)

write.csv(
  comparison,
  file.path(output_dir, "outlier_sensitivity_comparison.csv"),
  row.names = FALSE
)

cat("\nSENSITIVITY ANALYSIS SUMMARY:\n")
print(summary_table)

cat("\nStage 20.3 completed successfully.\n")
cat("Original Stage 19 results were preserved.\n")