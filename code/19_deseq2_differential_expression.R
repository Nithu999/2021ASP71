# ============================================================
# CropStress Intelligence
# Stage 19: Differential Expression Analysis using DESeq2
# ============================================================
# Run only after Stage 18 completes successfully.
#
# The analysis uses separate two-group comparisons at each
# selected time point and tissue. It does not combine the
# drought and rehydration experiments.
#
# IMPORTANT LIMITATION:
# The source count matrix has no explicit biological gene IDs.
# feature_XXXXX labels are technical row identifiers only.
# Results can be statistically tested at row level, but cannot
# be assigned biological gene names or interpreted as specific
# genes until the original row-to-gene mapping is recovered.

cat("\nSTAGE 19: DESeq2 DIFFERENTIAL EXPRESSION ANALYSIS\n")

if (!requireNamespace("DESeq2", quietly = TRUE)) {
  stop(
    "DESeq2 is not installed. In R, run:\n",
    'if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")\n',
    'BiocManager::install("DESeq2")\n',
    "Then rerun this script."
  )
}

library(DESeq2)

counts <- readRDS("data/processed/cleaned_filtered_counts.rds")
metadata <- read.csv(
  "results/tables/integrated_sample_metadata.csv",
  stringsAsFactors = FALSE
)

stopifnot(
  identical(colnames(counts), metadata$library_name),
  all(counts >= 0),
  all(counts == round(counts)),
  !is.null(rownames(counts)),
  !anyDuplicated(rownames(counts))
)

# Four time-specific comparisons in each of two tissues = eight analyses.
comparisons <- data.frame(
  experiment = c(
    "severe_drought", "severe_drought",
    "rehydration", "rehydration"
  ),
  day = c("6", "9", "3", "6"),
  control_label = c(
    "control (severe drought)",
    "control (severe drought)",
    "control (rehydration)",
    "control (rehydration)"
  ),
  treatment_label = c(
    "severe drought",
    "severe drought",
    "rehydration after severe drought",
    "rehydration after severe drought"
  ),
  stringsAsFactors = FALSE
)

tissues <- c("latest leaf", "second latest leaf")

dir.create("results/differential_expression", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures/differential_expression", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables/differential_expression", recursive = TRUE, showWarnings = FALSE)

summary_rows <- list()
analysis_index <- 1L

for (tissue_name in tissues) {
  for (i in seq_len(nrow(comparisons))) {

    cmp <- comparisons[i, ]
    comparison_id <- paste(
      cmp$experiment,
      paste0("day_", gsub("[^A-Za-z0-9]", "_", cmp$day)),
      gsub(" ", "_", tissue_name),
      sep = "__"
    )

    selected <- metadata[
      metadata$tissue == tissue_name &
        as.character(metadata$day) == cmp$day &
        metadata$treatment %in% c(cmp$control_label, cmp$treatment_label),
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

    cat("\n----------------------------------------\n")
    cat("Comparison:", comparison_id, "\n")
    cat("Tissue:", tissue_name, "| Day:", cmp$day, "\n")
    print(group_counts)

    if (!all(c("Control", "Treatment") %in% names(group_counts)) ||
        group_counts["Control"] != 3 ||
        group_counts["Treatment"] != 3) {
      stop(
        "Expected exactly 3 control and 3 treatment samples for ",
        comparison_id,
        ". Inspect metadata before continuing."
      )
    }

    selected_counts <- counts[
      ,
      selected$library_name,
      drop = FALSE
    ]

    # Comparison-specific pre-filter: count >=10 in at least 3 of 6 samples.
    keep <- rowSums(selected_counts >= 10) >= 3
    selected_counts <- selected_counts[keep, , drop = FALSE]

    if (nrow(selected_counts) < 2) {
      stop("Too few rows remain after filtering for ", comparison_id)
    }

    rownames(selected) <- selected$library_name
    selected <- selected[colnames(selected_counts), , drop = FALSE]

    dds <- DESeqDataSetFromMatrix(
      countData = round(selected_counts),
      colData = selected,
      design = ~ condition
    )

    # Positive-count size factors accommodate zero-heavy count matrices.
    dds <- estimateSizeFactors(dds, type = "poscounts")
    dds <- DESeq(dds, quiet = TRUE)

    res <- results(
      dds,
      contrast = c("condition", "Treatment", "Control"),
      alpha = 0.05
    )

    res_df <- as.data.frame(res)
    res_df$feature_id <- rownames(res_df)
    res_df <- res_df[
      order(res_df$padj, na.last = TRUE),
      c(
        "feature_id", "baseMean", "log2FoldChange",
        "lfcSE", "stat", "pvalue", "padj"
      )
    ]

    res_df$significant_FDR_0_05 <- (
      !is.na(res_df$padj) & res_df$padj < 0.05
    )
    res_df$significant_and_absLFC_ge_1 <- (
      !is.na(res_df$padj) &
        res_df$padj < 0.05 &
        abs(res_df$log2FoldChange) >= 1
    )

    result_file <- file.path(
      "results/tables/differential_expression",
      paste0(comparison_id, "_DESeq2_results.csv")
    )
    write.csv(res_df, result_file, row.names = FALSE)

    # Save normalized counts for traceability.
    norm_counts <- counts(dds, normalized = TRUE)
    norm_df <- data.frame(
      feature_id = rownames(norm_counts),
      norm_counts,
      check.names = FALSE
    )
    write.csv(
      norm_df,
      file.path(
        "results/differential_expression",
        paste0(comparison_id, "_normalized_counts.csv")
      ),
      row.names = FALSE
    )

    # MA plot
    png(
      file.path(
        "results/figures/differential_expression",
        paste0(comparison_id, "_MA.png")
      ),
      width = 1200, height = 900, res = 140
    )
    plotMA(res, main = comparison_id, ylim = c(-5, 5))
    dev.off()

    # Volcano plot using adjusted p-values
    volcano_p <- -log10(pmax(res_df$padj, .Machine$double.xmin))
    volcano_p[is.na(res_df$padj)] <- NA_real_

    png(
      file.path(
        "results/figures/differential_expression",
        paste0(comparison_id, "_volcano.png")
      ),
      width = 1200, height = 900, res = 140
    )
    plot(
      res_df$log2FoldChange,
      volcano_p,
      pch = 20,
      cex = 0.55,
      xlab = "Log2 fold change (Treatment vs Control)",
      ylab = "-Log10 adjusted p-value",
      main = comparison_id
    )
    abline(v = c(-1, 1), lty = 2)
    abline(h = -log10(0.05), lty = 2)
    dev.off()

    summary_rows[[analysis_index]] <- data.frame(
      comparison_id = comparison_id,
      experiment = cmp$experiment,
      day = cmp$day,
      tissue = tissue_name,
      control_n = unname(group_counts["Control"]),
      treatment_n = unname(group_counts["Treatment"]),
      rows_tested = nrow(res_df),
      significant_FDR_0_05 = sum(res_df$significant_FDR_0_05, na.rm = TRUE),
      significant_FDR_0_05_absLFC_ge_1 = sum(
        res_df$significant_and_absLFC_ge_1,
        na.rm = TRUE
      ),
      result_file = result_file,
      stringsAsFactors = FALSE
    )

    analysis_index <- analysis_index + 1L
    cat("Saved:", result_file, "\n")
  }
}

de_summary <- do.call(rbind, summary_rows)
write.csv(
  de_summary,
  "results/tables/differential_expression/DESeq2_analysis_summary.csv",
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  "results/differential_expression/sessionInfo.txt"
)

cat("\nDESeq2 analysis summary:\n")
print(de_summary)

cat(
  "\nStage 19 completed. These are row-level results using technical feature IDs.\n",
  "Do not claim gene-level biological findings until the original gene mapping is validated.\n"
)
