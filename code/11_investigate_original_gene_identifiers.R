# ============================================================
# Stage 11: Investigate Original Gene Identifiers
# ============================================================

# Original raw count file

raw_file <- paste0(
  "data/raw/",
  "GSE288615_rawCounts_IR64_rehydration_All_samples.txt"
)

# STEP 11.1: Inspect original file header

cat("\nOriginal file - first five lines:\n")

first_lines <- readLines(
  raw_file,
  n = 5
)

cat(first_lines, sep = "\n")

# STEP 11.2: Import a small preview

preview <- read.delim(
  raw_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  nrows = 5
)

cat("\nPreview dimensions:\n")
print(dim(preview))

# STEP 11.3: Inspect column names

cat("\nFirst 10 column names:\n")
print(head(names(preview), 10))

# STEP 11.4: Inspect first five records

cat("\nFirst five records:\n")
print(head(preview))

# STEP 11.5: Inspect first column

cat("\nFirst column structure:\n")
str(preview[[1]])

cat("\nFirst column values:\n")
print(head(preview[[1]], 10))

cat("\nInvestigation completed.\n")