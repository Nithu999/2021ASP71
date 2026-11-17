# ============================================================
# Stage 10: Validate Differential Expression Inputs
# ============================================================

# Load metadata and filtered counts

metadata <- read.csv(
  "results/tables/integrated_sample_metadata.csv",
  stringsAsFactors = FALSE
)

filtered_counts <- readRDS(
  "data/processed/filtered_counts.rds"
)

# STEP 10.1: Check sample names

cat("\nSample name validation:\n")

print(identical(
  colnames(filtered_counts),
  metadata$library_name
))

# STEP 10.2: Check count matrix

cat("\nFiltered count matrix dimensions:\n")
print(dim(filtered_counts))

cat("\nAre counts integers?\n")
print(all(filtered_counts == round(filtered_counts)))

cat("\nAre counts non-negative?\n")
print(all(filtered_counts >= 0))

# STEP 10.3: Check gene identifiers

cat("\nFirst 10 row identifiers:\n")
print(head(rownames(filtered_counts), 10))

# STEP 10.4: Check sample groups

cat("\nTreatment and day combinations:\n")
print(with(metadata, table(treatment, day)))

# STEP 10.5: Check missing metadata

cat("\nMissing values:\n")
print(colSums(is.na(metadata)))

cat("\nInput validation completed.\n")