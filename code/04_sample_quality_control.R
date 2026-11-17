# CropStress Intelligence
# Stage 4: Sample-level RNA-seq quality assessment
# Output: sample QC summaries, replicate comparisons and figures

# STEP 19: Investigate unusually low-count sample


target_sample <- "IR9AD2L1"


stopifnot(
  target_sample %in% colnames(counts)
)


library_sizes <- colSums(
  counts
)


detected_records <- colSums(
  counts > 0
)


sample_qc <- data.frame(

  library_name = names(library_sizes),

  library_size = as.numeric(library_sizes),

  detected_records = as.numeric(detected_records),

  stringsAsFactors = FALSE

)


cat("\nTarget sample quality information:\n")

print(
  sample_qc[
    sample_qc$library_name == target_sample,
  ]
)


cat("\nLibrary size summary:\n")

print(
  summary(sample_qc$library_size)
)


cat("\nDetected expression records summary:\n")

print(
  summary(sample_qc$detected_records)
)

# STEP 19.2: Compare sample with experimental group


target_metadata <- integrated_metadata[
  integrated_metadata$library_name == target_sample,
]


cat("\nTarget sample metadata:\n")

print(target_metadata)


same_group <- integrated_metadata[

  integrated_metadata$treatment == target_metadata$treatment &

    integrated_metadata$day == target_metadata$day &

    integrated_metadata$tissue == target_metadata$tissue,

]


cat("\nSamples from the same experimental group:\n")

print(
  same_group[
    ,
    c(
      "geo_accession",
      "library_name",
      "treatment",
      "day",
      "tissue",
      "replicate"
    )
  ]
)


group_qc <- merge(

  same_group,

  sample_qc,

  by = "library_name"

)


cat("\nQuality metrics for the experimental group:\n")

print(

  group_qc[
    ,
    c(
      "library_name",
      "replicate",
      "library_size",
      "detected_records"
    )
  ]

)

# STEP 19.3: Calculate relative library size


median_library_size <- median(
  sample_qc$library_size
)


target_library_size <- sample_qc$library_size[
  sample_qc$library_name == target_sample
]


relative_library_size <- (
  target_library_size / median_library_size
) * 100


cat("\nTarget sample library size:", target_library_size)

cat("\nMedian library size:", median_library_size)

cat(
  "\nTarget sample relative to median (%):",
  round(relative_library_size, 4),
  "%\n"
)


fold_difference <- median_library_size / target_library_size


cat(
  "\nMedian library size is",
  round(fold_difference, 2),
  "times larger than target sample.\n"
)

# STEP 19.4: Visualize biological replicate comparison


group_plot_data <- group_qc[
  ,
  c(
    "library_name",
    "replicate",
    "library_size"
  )
]


barplot(
  height = group_plot_data$library_size,

  names.arg = group_plot_data$library_name,

  las = 2,

  ylab = "Total raw counts",

  main = "Library Size Comparison: Severe Drought Day 9",

  cex.names = 0.8
)


png(
  filename = file.path(
    "results",
    "figures",
    "severe_drought_day9_replicates.png"
  ),

  width = 1200,

  height = 800,

  res = 120
)

barplot(
  height = group_plot_data$library_size,

  names.arg = group_plot_data$library_name,

  las = 2,

  ylab = "Total raw counts",

  main = "Library Size Comparison: Severe Drought Day 9",

  cex.names = 0.8
)

dev.off()

# STEP 19.5: Investigate target sample count distribution


target_counts <- counts[[target_sample]]


cat("\nTarget sample count distribution:\n")

print(summary(target_counts))


zero_counts <- sum(target_counts == 0)

cat("\nZero-count records:", zero_counts)


nonzero_counts <- sum(target_counts > 0)

cat("\nNon-zero records:", nonzero_counts)


top_counts <- sort(
  target_counts,
  decreasing = TRUE
)


cat("\nTop 10 expression counts:\n")

print(head(top_counts, 10))

# STEP 19.6: Compare count distributions across replicates


replicate_samples <- group_qc$library_name


replicate_distribution <- data.frame(

  library_name = replicate_samples,

  total_counts = NA,

  zero_records = NA,

  nonzero_records = NA,

  maximum_count = NA,

  mean_count = NA

)


for (i in seq_along(replicate_samples)) {

  sample_name <- replicate_samples[i]

  current_counts <- counts[[sample_name]]

  replicate_distribution$total_counts[i] <- sum(current_counts)

  replicate_distribution$zero_records[i] <- sum(
    current_counts == 0
  )

  replicate_distribution$nonzero_records[i] <- sum(
    current_counts > 0
  )

  replicate_distribution$maximum_count[i] <- max(
    current_counts
  )

  replicate_distribution$mean_count[i] <- mean(
    current_counts
  )

}


cat("\nBiological replicate count distribution:\n")

print(replicate_distribution)

# STEP 19.7: Save replicate count distribution


output_file <- file.path(
  "results",
  "tables",
  "replicate_count_distribution.csv"
)


write.csv(
  replicate_distribution,
  file = output_file,
  row.names = FALSE
)


cat("\nFile successfully saved:",
    file.exists(output_file))


cat("\nOutput location:", output_file)

