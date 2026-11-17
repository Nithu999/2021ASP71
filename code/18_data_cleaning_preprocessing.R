# ============================================================
# CropStress Intelligence
# Stage 18: Data Cleaning and Preprocessing
# ============================================================
# Run from the CropStress_Intelligence project root.
#
# Important:
# - This script does NOT remove the low-library sample IR9AD2L1.
# - It creates technical feature IDs only to track count-matrix rows.
# - These feature IDs are NOT biological gene identifiers.
# - No MSU7 identifiers are assigned by row order.

cat("\nSTAGE 18: DATA CLEANING AND PREPROCESSING\n")

# 1. Load input files
counts_file <- "data/raw/GSE288615_rawCounts_IR64_rehydration_All_samples.txt"
metadata_file <- "results/tables/integrated_sample_metadata.csv"

stopifnot(file.exists(counts_file), file.exists(metadata_file))

counts_df <- read.delim(
  counts_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

metadata <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# 2. Validate structure
cat("\nRaw count dimensions:\n")
print(dim(counts_df))

if (nrow(counts_df) != 55986L || ncol(counts_df) != 72L) {
  warning("The raw count dimensions differ from the previously observed 55,986 x 72.")
}

if (anyDuplicated(colnames(counts_df))) {
  stop("Duplicate sample names found in the count matrix.")
}

if (!all(vapply(counts_df, is.numeric, logical(1)))) {
  stop("At least one count column is not numeric. Inspect the raw input before continuing.")
}

count_matrix <- as.matrix(counts_df)
storage.mode(count_matrix) <- "numeric"

if (anyNA(count_matrix) || any(!is.finite(count_matrix))) {
  stop("Missing or non-finite count values detected.")
}
if (any(count_matrix < 0)) {
  stop("Negative count values detected.")
}
if (any(count_matrix != round(count_matrix))) {
  stop("Non-integer values detected. DESeq2 requires un-normalized integer counts.")
}

# 3. Validate sample metadata and sample order
required_metadata <- c(
  "geo_accession", "library_name", "cultivar", "tissue",
  "treatment", "day", "replicate"
)
missing_columns <- setdiff(required_metadata, colnames(metadata))
if (length(missing_columns) > 0) {
  stop("Missing metadata columns: ", paste(missing_columns, collapse = ", "))
}
if (anyDuplicated(metadata$library_name)) {
  stop("Duplicate library names found in metadata.")
}
if (!identical(colnames(count_matrix), metadata$library_name)) {
  stop("Count-matrix sample order does not exactly match metadata library_name order.")
}

# 4. Create traceable technical feature IDs.
# These are row labels only, not biological gene IDs.
feature_ids <- sprintf("feature_%05d", seq_len(nrow(count_matrix)))
rownames(count_matrix) <- feature_ids

feature_key <- data.frame(
  feature_id = feature_ids,
  original_count_matrix_row = seq_len(nrow(count_matrix)),
  stringsAsFactors = FALSE
)

# 5. Recalculate sample-level QC from the raw count matrix
sample_qc <- data.frame(
  library_name = colnames(count_matrix),
  library_size = colSums(count_matrix),
  detected_records = colSums(count_matrix > 0),
  zero_count_records = colSums(count_matrix == 0),
  stringsAsFactors = FALSE
)

sample_qc <- merge(
  metadata,
  sample_qc,
  by = "library_name",
  sort = FALSE
)

sample_qc$zero_count_percentage <- round(
  100 * sample_qc$zero_count_records / nrow(count_matrix),
  3
)

# Robust flag on log10 library size. A flag triggers review, not exclusion.
log_library_size <- log10(pmax(sample_qc$library_size, 1))
log_mad <- mad(log_library_size, constant = 1.4826)
if (is.finite(log_mad) && log_mad > 0) {
  sample_qc$robust_library_size_z <- (
    log_library_size - median(log_library_size)
  ) / log_mad
} else {
  sample_qc$robust_library_size_z <- NA_real_
}

sample_qc$low_library_flag <- (
  !is.na(sample_qc$robust_library_size_z) &
    sample_qc$robust_library_size_z < -3.5
)

sample_qc <- sample_qc[
  match(metadata$library_name, sample_qc$library_name),
]

# 6. Preliminary low-expression filtering
# Retain a row if raw count >=10 in at least 3 of the 72 samples.
keep <- rowSums(count_matrix >= 10) >= 3
cleaned_filtered_counts <- count_matrix[keep, , drop = FALSE]

filter_summary <- data.frame(
  original_records = nrow(count_matrix),
  retained_records = nrow(cleaned_filtered_counts),
  removed_records = nrow(count_matrix) - nrow(cleaned_filtered_counts),
  retention_percentage = round(100 * nrow(cleaned_filtered_counts) /
                                 nrow(count_matrix), 2),
  filtering_rule = "Raw count >= 10 in at least 3 of 72 samples",
  samples_retained = ncol(cleaned_filtered_counts),
  stringsAsFactors = FALSE
)

# 7. Save processed data and audit outputs
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

saveRDS(
  cleaned_filtered_counts,
  "data/processed/cleaned_filtered_counts.rds"
)

write.csv(
  sample_qc,
  "results/tables/preprocessing_sample_qc.csv",
  row.names = FALSE
)

write.csv(
  filter_summary,
  "results/tables/preprocessing_filtering_summary.csv",
  row.names = FALSE
)

write.csv(
  feature_key,
  "results/tables/technical_feature_id_key.csv",
  row.names = FALSE
)

# 8. Print key checks
cat("\nSample QC summary:\n")
print(summary(sample_qc$library_size))

cat("\nFlagged low-library samples (review only; not removed):\n")
print(sample_qc[
  sample_qc$low_library_flag,
  c("library_name", "treatment", "day", "tissue",
    "replicate", "library_size", "detected_records",
    "zero_count_percentage", "robust_library_size_z")
])

cat("\nFiltering summary:\n")
print(filter_summary)

cat("\nProcessed matrix dimensions:\n")
print(dim(cleaned_filtered_counts))

stopifnot(
  ncol(cleaned_filtered_counts) == 72L,
  all(cleaned_filtered_counts >= 0),
  all(cleaned_filtered_counts == round(cleaned_filtered_counts)),
  identical(colnames(cleaned_filtered_counts), metadata$library_name)
)

cat("\nStage 18 completed. No samples were removed.\n")
cat("Note: feature_XXXXX labels are technical row IDs, not gene IDs.\n")
