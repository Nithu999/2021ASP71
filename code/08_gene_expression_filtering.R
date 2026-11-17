# ============================================================
# CropStress Intelligence
# Stage 8: Gene Expression Quality Assessment and Filtering
# ============================================================

# STEP 41: Validate input count matrix

cat("\nOriginal count matrix dimensions:\n")
print(dim(counts))

stopifnot(
  all(sapply(counts, is.numeric))
)

stopifnot(
  all(as.matrix(counts) >= 0)
)

# Convert to numeric matrix

count_matrix <- as.matrix(counts)

# STEP 42: Calculate expression statistics

expression_statistics <- data.frame(

  record_id = seq_len(nrow(count_matrix)),

  mean_count = rowMeans(count_matrix),

  median_count = apply(
    count_matrix,
    1,
    median
  ),

  maximum_count = apply(
    count_matrix,
    1,
    max
  ),

  zero_samples = rowSums(
    count_matrix == 0
  ),

  nonzero_samples = rowSums(
    count_matrix > 0
  ),

  samples_with_count_10_or_more = rowSums(
    count_matrix >= 10
  )

)

# STEP 43: Investigate expression distribution

cat("\nExpression statistics:\n")

print(summary(expression_statistics))

# STEP 44: Identify low-expression records

minimum_count <- 10

minimum_samples <- 3

keep_records <- rowSums(
  count_matrix >= minimum_count
) >= minimum_samples

# STEP 45: Generate filtered count matrix

filtered_counts <- count_matrix[
  keep_records,
  ,
  drop = FALSE
]

# STEP 46: Calculate filtering summary

original_records <- nrow(count_matrix)

retained_records <- nrow(filtered_counts)

removed_records <- original_records - retained_records

retention_percentage <- (
  retained_records / original_records
) * 100

filtering_summary <- data.frame(

  original_records = original_records,

  retained_records = retained_records,

  removed_records = removed_records,

  retention_percentage = round(
    retention_percentage,
    2
  ),

  minimum_count = minimum_count,

  minimum_samples = minimum_samples

)

cat("\nFiltering summary:\n")

print(filtering_summary)

# STEP 47: Save expression statistics

write.csv(
  expression_statistics,
  "results/tables/gene_expression_statistics.csv",
  row.names = FALSE
)

# STEP 48: Save filtering summary

write.csv(
  filtering_summary,
  "results/tables/gene_expression_filtering_summary.csv",
  row.names = FALSE
)

# STEP 49: Save filtered count matrix

saveRDS(
  filtered_counts,
  "data/processed/filtered_counts.rds"
)

# STEP 50: Visualize filtering results

png(
  "results/figures/gene_expression_filtering.png",
  width = 1200,
  height = 800,
  res = 120
)

barplot(
  c(
    Original = original_records,
    Retained = retained_records,
    Removed = removed_records
  ),
  main = "Gene Expression Filtering Results",
  ylab = "Number of Expression Records",
  col = c("steelblue", "forestgreen", "tomato")
)

dev.off()

# STEP 51: Final validation

stopifnot(
  ncol(filtered_counts) == ncol(count_matrix)
)

stopifnot(
  nrow(filtered_counts) <= nrow(count_matrix)
)

cat("\nFiltered matrix dimensions:\n")

print(dim(filtered_counts))

cat("\nFiltering completed successfully.\n")