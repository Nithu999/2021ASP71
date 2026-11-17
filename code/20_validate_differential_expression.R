# ============================================================
# CropStress Intelligence
# Stage 20: DESeq2 Diagnostic Analysis
# ============================================================
#
# Purpose:
# Investigate missing adjusted p-values caused by:
# 1. Independent filtering
# 2. Cook's distance outlier detection
#
# Original Stage 19 results will NOT be modified.
# ============================================================

cat("\nSTAGE 20: DESeq2 DIAGNOSTIC ANALYSIS\n")

library(DESeq2)

# ------------------------------------------------------------
# STEP 1: Load original data
# ------------------------------------------------------------

counts <- readRDS(
  "data/processed/cleaned_filtered_counts.rds"
)

metadata <- read.csv(
  "results/tables/integrated_sample_metadata.csv",
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# STEP 2: Define affected comparisons
# ------------------------------------------------------------

comparisons <- data.frame(

  experiment = c(
    "severe_drought",
    "severe_drought"
  ),

  day = c("9", "6"),

  tissue = c(
    "latest leaf",
    "second latest leaf"
  ),

  control_label = c(
    "control (severe drought)",
    "control (severe drought)"
  ),

  treatment_label = c(
    "severe drought",
    "severe drought"
  ),

  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# STEP 3: Create output directories
# ------------------------------------------------------------

dir.create(
  "results/tables/differential_expression/stage20",
  recursive = TRUE,
  showWarnings = FALSE
)

diagnostic_summary <- list()

# ------------------------------------------------------------
# STEP 4: Run diagnostic comparisons
# ------------------------------------------------------------

for (i in seq_len(nrow(comparisons))) {

  cmp <- comparisons[i, ]

  comparison_id <- paste(
    cmp$experiment,
    paste0("day_", cmp$day),
    gsub(" ", "_", cmp$tissue),
    sep = "__"
  )

  cat("\n----------------------------------------\n")
  cat("Comparison:", comparison_id, "\n")

  selected <- metadata[
    metadata$tissue == cmp$tissue &
      as.character(metadata$day) == cmp$day &
      metadata$treatment %in% c(
        cmp$control_label,
        cmp$treatment_label
      ),
    ,
    drop = FALSE
  ]

  selected$condition <- ifelse(
    selected$treatment == cmp$control_label,
    "Control",
    "Treatment"
  )

  selected$condition <- factor(
    selected$condition,
    levels = c("Control", "Treatment")
  )

  group_counts <- table(selected$condition)

  print(group_counts)

  if (group_counts["Control"] != 3 ||
      group_counts["Treatment"] != 3) {
    stop("Unexpected sample numbers: ", comparison_id)
  }

  selected_counts <- counts[
    ,
    selected$library_name,
    drop = FALSE
  ]

  # Same comparison-specific filtering as Stage 19
  keep <- rowSums(selected_counts >= 10) >= 3

  selected_counts <- selected_counts[
    keep,
    ,
    drop = FALSE
  ]

  rownames(selected) <- selected$library_name

  selected <- selected[
    colnames(selected_counts),
    ,
    drop = FALSE
  ]

  # ----------------------------------------------------------
  # STEP 5: Recreate DESeq2 dataset
  # ----------------------------------------------------------

  dds <- DESeqDataSetFromMatrix(
    countData = round(selected_counts),
    colData = selected,
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

  # ----------------------------------------------------------
  # STEP 6: Original default result behaviour
  # ----------------------------------------------------------

  original <- results(
    dds,
    contrast = c(
      "condition",
      "Treatment",
      "Control"
    ),
    alpha = 0.05
  )

  # ----------------------------------------------------------
  # STEP 7: Disable independent filtering
  # ----------------------------------------------------------

  no_independent_filtering <- results(
    dds,
    contrast = c(
      "condition",
      "Treatment",
      "Control"
    ),
    alpha = 0.05,
    independentFiltering = FALSE
  )

  # ----------------------------------------------------------
  # STEP 8: Disable Cook's cutoff
  # ----------------------------------------------------------

  no_cooks_cutoff <- results(
    dds,
    contrast = c(
      "condition",
      "Treatment",
      "Control"
    ),
    alpha = 0.05,
    cooksCutoff = FALSE
  )

  # ----------------------------------------------------------
  # STEP 9: Disable both options
  # ----------------------------------------------------------

  both_disabled <- results(
    dds,
    contrast = c(
      "condition",
      "Treatment",
      "Control"
    ),
    alpha = 0.05,
    independentFiltering = FALSE,
    cooksCutoff = FALSE
  )

  # ----------------------------------------------------------
  # STEP 10: Summarise missing values
  # ----------------------------------------------------------

  diagnostic_summary[[i]] <- data.frame(

    comparison = comparison_id,

    total_features = length(original$padj),

    original_missing_pvalue =
      sum(is.na(original$pvalue)),

    original_missing_padj =
      sum(is.na(original$padj)),

    no_filter_missing_padj =
      sum(is.na(no_independent_filtering$padj)),

    no_cooks_missing_pvalue =
      sum(is.na(no_cooks_cutoff$pvalue)),

    no_cooks_missing_padj =
      sum(is.na(no_cooks_cutoff$padj)),

    both_disabled_missing_pvalue =
      sum(is.na(both_disabled$pvalue)),

    both_disabled_missing_padj =
      sum(is.na(both_disabled$padj)),

    original_significant =
      sum(original$padj < 0.05, na.rm = TRUE),

    no_filter_significant =
      sum(
        no_independent_filtering$padj < 0.05,
        na.rm = TRUE
      ),

    stringsAsFactors = FALSE
  )

  # ----------------------------------------------------------
  # STEP 11: Save diagnostic result comparison
  # ----------------------------------------------------------

  comparison_df <- data.frame(

    feature_id = rownames(original),

    original_pvalue = original$pvalue,

    original_padj = original$padj,

    no_filter_padj =
      no_independent_filtering$padj,

    no_cooks_pvalue =
      no_cooks_cutoff$pvalue,

    no_cooks_padj =
      no_cooks_cutoff$padj,

    both_disabled_pvalue =
      both_disabled$pvalue,

    both_disabled_padj =
      both_disabled$padj

  )

  output_file <- file.path(
    "results/tables/differential_expression/stage20",
    paste0(
      comparison_id,
      "_diagnostic_comparison.csv"
    )
  )

  write.csv(
    comparison_df,
    output_file,
    row.names = FALSE
  )

  cat("Saved:", output_file, "\n")
}

# ------------------------------------------------------------
# STEP 12: Save diagnostic summary
# ------------------------------------------------------------

diagnostic_table <- do.call(
  rbind,
  diagnostic_summary
)

rownames(diagnostic_table) <- NULL

write.csv(
  diagnostic_table,
  "results/tables/differential_expression/stage20/DESeq2_diagnostic_summary.csv",
  row.names = FALSE
)

cat("\nDIAGNOSTIC SUMMARY:\n")

print(diagnostic_table)

cat(
  "\nStage 20 diagnostic analysis completed.\n",
  "Original Stage 19 results were not modified.\n"
)