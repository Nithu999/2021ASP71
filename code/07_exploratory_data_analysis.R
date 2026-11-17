# ============================================================
# CropStress Intelligence
# Stage 7: Exploratory Data Analysis
# ============================================================


# STEP 38: Investigate experimental groups
# ============================================================

# STEP 38.1: Verify available data

cat("\nCount matrix dimensions:\n")
print(dim(counts))

cat("\nIntegrated metadata dimensions:\n")
print(dim(integrated_metadata))


# STEP 38.2: Inspect experimental factors

cat("\nTreatment groups:\n")
print(table(integrated_metadata$treatment))

cat("\nSampling days:\n")
print(table(integrated_metadata$day))

cat("\nTissue types:\n")
print(table(integrated_metadata$tissue))


# STEP 38.3: Investigate treatment and day combinations

cat("\nTreatment and day combinations:\n")

print(
  with(
    integrated_metadata,
    table(treatment, day)
  )
)


# STEP 38.4: Investigate treatment, day and tissue combinations

cat("\nExperimental group distribution:\n")

experimental_groups <- as.data.frame(
  with(
    integrated_metadata,
    table(treatment, day, tissue)
  )
)

experimental_groups <- experimental_groups[
  experimental_groups$Freq > 0,
]

print(experimental_groups)


# STEP 38.5: Save experimental group summary

write.csv(
  experimental_groups,
  file = "results/tables/experimental_group_summary.csv",
  row.names = FALSE
)

cat("\nExperimental group summary saved successfully.\n")


# STEP 39: Visualize experimental group distribution
# ============================================================

# STEP 39.1: Create treatment and day distribution plot

treatment_day_table <- table(
  integrated_metadata$treatment,
  integrated_metadata$day
)

png(
  filename = "results/figures/experimental_group_distribution.png",
  width = 1400,
  height = 900,
  res = 150
)

par(
  mar = c(10, 5, 4, 2)
)

barplot(
  height = treatment_day_table,
  beside = TRUE,
  col = c("skyblue", "orange", "lightgreen", "tomato"),
  legend.text = rownames(treatment_day_table),
  args.legend = list(
    x = "topright",
    cex = 0.7
  ),
  las = 2,
  ylab = "Number of samples",
  xlab = "Sampling day",
  main = "Distribution of Samples Across Treatment Groups"
)

dev.off()

cat("\nExperimental group distribution plot saved.\n")


# STEP 40: Investigate raw count distribution
# ============================================================

# STEP 40.1: Validate count matrix

cat("\nChecking count matrix:\n")

stopifnot(
  is.matrix(counts) || is.data.frame(counts)
)

stopifnot(
  is.numeric(as.matrix(counts))
)

stopifnot(
  all(is.finite(as.matrix(counts)))
)

stopifnot(
  all(as.matrix(counts) >= 0)
)

cat("\nCount matrix validation successful.\n")


# STEP 40.2: Calculate overall count summary

cat("\nOverall count summary:\n")

print(
  summary(
    as.numeric(as.matrix(counts))
  )
)

# STEP 40.3: Count zero and non-zero values

total_values <- length(as.matrix(counts))

zero_values <- sum(as.matrix(counts) == 0)

nonzero_values <- sum(as.matrix(counts) > 0)

cat("\nTotal count values:", total_values)

cat("\nZero count values:", zero_values)

cat("\nNon-zero count values:", nonzero_values)

cat(
  "\nPercentage of zero values:",
  round((zero_values / total_values) * 100, 2),
  "%\n"
)

# STEP 40.4: Calculate gene-level statistics

gene_statistics <- data.frame(
  
  gene_id = rownames(counts),
  
  mean_count = rowMeans(counts),
  
  median_count = apply(counts, 1, median),
  
  maximum_count = apply(counts, 1, max),
  
  zero_samples = rowSums(counts == 0)
  
)


cat("\nGene-level statistics summary:\n")

print(
  summary(
    gene_statistics[
      ,
      c(
        "mean_count",
        "median_count",
        "maximum_count",
        "zero_samples"
      )
    ]
  )
)


# STEP 40.5: Save gene-level statistics

write.csv(
  gene_statistics,
  file = "results/tables/gene_count_statistics.csv",
  row.names = FALSE
)

cat("\nGene-level statistics saved successfully.\n")


# ============================================================
# STAGE 7 COMPLETED
# ============================================================

cat("\n========================================\n")
cat("Stage 7: Exploratory Data Analysis\n")
cat("Completed successfully.\n")
cat("========================================\n")