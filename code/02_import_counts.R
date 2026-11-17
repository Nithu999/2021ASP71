# CropStress Intelligence
# Stage 2: Import and inspect the raw RNA-seq count matrix
# Input: data/raw/GSE288615_rawCounts_IR64_rehydration_All_samples.txt
# Output: counts object in the R environment

# STEP 12: Import and inspect raw RNA-seq count matrix


counts_file <- file.path(
  "data",
  "raw",
  "GSE288615_rawCounts_IR64_rehydration_All_samples.txt"
)

# Do not assign row.names yet.
# We first need to understand the file structure.

counts <- read.delim(
  counts_file,
  header = TRUE,
  check.names = FALSE
)


cat("\nCount matrix dimensions:\n")
print(dim(counts))


cat("\nFirst 10 column names:\n")
print(head(colnames(counts), 10))


cat("\nFirst six rows:\n")
print(head(counts))


cat("\nCount matrix structure:\n")
str(counts)


cat("\nAre all columns numeric?\n")

print(
  all(sapply(counts, is.numeric))
)


cat("\nExpression matrix dimensions:\n")
print(dim(counts))


cat("\nNumber of samples:\n")
print(ncol(counts))


cat("\nNumber of expression records:\n")
print(nrow(counts))


cat("\nSample identifiers:\n")
print(colnames(counts))


cat("\nFirst five expression records:\n")
print(head(counts, 5))


# STEP 13: Search GEO metadata for sample mapping information


mapping_lines <- grep(
  "Library name|Column name in rawCounts|IR0C1L1",
  soft_lines,
  value = TRUE,
  ignore.case = TRUE
)


cat("\nOfficial sample mapping information:\n")

print(head(mapping_lines, 30))


cat(
  "\nTotal matching entries:",
  length(mapping_lines),
  "\n"
)


