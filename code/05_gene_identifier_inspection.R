# CropStress Intelligence
# Stage 5: Preliminary inspection of gene identifiers
# Important: row numbers are not verified gene identifiers.

# STEP 20: Investigate gene identifiers

# STEP 20.1: Inspect count matrix structure

cat("\nCount matrix dimensions:\n")
print(dim(counts))

cat("\nColumn names:\n")
print(head(colnames(counts)))

cat("\nFirst six expression records:\n")
print(head(counts))

cat("\nData types:\n")
print(sapply(counts, class))

cat("\nRow names:\n")
print(head(rownames(counts), 10))

# STEP 20.2: Inspect original count file

count_file_connection <- file(
  counts_file,
  open = "r"
)

first_lines <- readLines(
  count_file_connection,
  n = 6
)

close(count_file_connection)

cat("\nFirst six lines of original count file:\n")

cat(
  first_lines,
  sep = "\n"
)

# STEP 20.3: Check row identifier uniqueness

cat("\nNumber of rows:", nrow(counts), "\n")

cat(
  "\nNumber of unique row names:",
  length(unique(rownames(counts))),
  "\n"
)

cat(
  "\nAre row names duplicated?",
  anyDuplicated(rownames(counts)) > 0,
  "\n"
)


